// src/cls_DepthResolutionSweep.cpp
#include "cls_DepthResolutionSweep.hpp"

#include <stdexcept>
#include <memory>
#include "cls_DetectorPanelArray.hpp"
#include "cls_Grid2dPillar.hpp"
#include "cls_Grid2dPillarParameters.hpp"
#include "cls_Grid2dVoxel.hpp"
#include "ns_mylogger.hpp"
#include "ns_mymacro.hpp"
#include "ns_geom_util.hpp"
#include "cls_Angle.hpp"
#include "ns_angle_util.hpp"
#include "ns_pathcalc.hpp"
#include "st_SignalNoiseStat.hpp"
#include "ns_iodir.hpp"

#include <cmath>
#include <limits>
#include <chrono>
#include "cls_PrefixSum2D.hpp"

// =======================================================================
// Anonymous namespace: Local helper functions
//
// These functions are internal to this translation unit and provide:
// - calc_ty_max_allowed: Computes elevation angle limit to ensure PL >= obj_size
// - mp_clamp_DL_nonnegative: OpenMP helper to clamp negative density lengths
// - st_add_DL: Single-thread DL addition (for use inside omp parallel regions)
// - st_clamp_DL_nonnegative: Single-thread DL clamping
// - st_recalc_peneflux_signal_from_DL: Single-thread flux/signal recomputation
// =======================================================================
namespace {

/// @brief Calculate ty_max_allowed such that PL >= obj_size
/// @param[in] elev_obj_from_det Vertical height from detector to object top [m]
/// @param[in] obj_size Object size (diameter) [m]
/// @param[in] unit Angle unit used by the detector panel
/// @return Maximum allowed ty value. Returns max double if no restriction needed.
/// @note If obj_size > elev_obj_from_det at all elevation angles, returns 0 (should skip this obj_size).
double calc_ty_max_allowed(double elev_obj_from_det, double obj_size, angle_util::AngleUnit unit)
{
  if (obj_size <= 0.0 || elev_obj_from_det <= 0.0) {
    return 0.0;  // invalid input, should skip
  }
  if (obj_size <= elev_obj_from_det) {
    return std::numeric_limits<double>::max();  // no restriction needed
  }

  // Geometry: For a vertical cylinder of height H and diameter D, the path length
  // PL at elevation angle theta satisfies: PL = D / cos(theta) when D < H.
  // When D >= H, there exists a maximum theta beyond which PL < D (ray exits top).
  // The critical angle: sin(theta_max) = H / D, giving PL_min = D at theta_max.
  const double sin_theta_max = elev_obj_from_det / obj_size;
  const double theta_max = std::asin(sin_theta_max);  // radians

  switch (unit) {
    case angle_util::AngleUnit::Tangent:
      return std::tan(theta_max);
    case angle_util::AngleUnit::Radian:
      return theta_max;
    case angle_util::AngleUnit::Degree:
      return theta_max * 180.0 / M_PI;
    default:
      return 0.0;  // unknown unit, should skip
  }
}

/// @brief Clamp negative DL values to zero in DetectorPanel
/// @param[in,out] panel DetectorPanel to process
/// @param[in] DL_min Minimum DL value (default 0.0)
/// @note Uses OpenMP for parallel processing
void mp_clamp_DL_nonnegative(DetectorPanel& panel, double DL_min = 0.0)
{
  const int nbinx = panel.get_nbinx();
  const int nbiny = panel.get_nbiny();

  #pragma omp parallel for collapse(2) schedule(static)
  for (int iy = 0; iy < nbiny; iy++) {
    for (int ix = 0; ix < nbinx; ix++) {
      DetectorElement& ele = panel.callDetectorElement(ix, iy);
      const double DL = ele.get_DL();
      if (DL < DL_min) {
        ele.set_DL(DL_min);
      }
    }
  }
}

// -----------------------------------------------------------------------
//  Single-thread helpers for use inside omp parallel for (Phase 2).
//
//  The existing mp_* functions contain internal #pragma omp parallel for.
//  To avoid nested OpenMP when the OUTER loop is parallelised, we use
//  these lightweight, OpenMP-free equivalents.
// -----------------------------------------------------------------------

/// @brief Add constant DL to all bins (single-thread version).
/// @note  Equivalent to pathcalc::naive::mp_add_DL without OpenMP.
void st_add_DL(DetectorPanel& panel, double delta_DL)
{
  const int nx = panel.get_nbinx();
  const int ny = panel.get_nbiny();
  for (int iy = 0; iy < ny; ++iy)
    for (int ix = 0; ix < nx; ++ix)
      panel.callDetectorElement(ix, iy).add_DL(delta_DL);
}

/// @brief Clamp negative DL to zero (single-thread version).
void st_clamp_DL_nonnegative(DetectorPanel& panel, double DL_min = 0.0)
{
  const int nx = panel.get_nbinx();
  const int ny = panel.get_nbiny();
  for (int iy = 0; iy < ny; ++iy)
    for (int ix = 0; ix < nx; ++ix) {
      DetectorElement& ele = panel.callDetectorElement(ix, iy);
      if (ele.get_DL() < DL_min) ele.set_DL(DL_min);
    }
}

/// @brief Recompute only peneflux and signal from updated DL.
///        omega, effective_area, exposure_time are assumed already set
///        (from Phase 1) and are NOT recomputed here.
///        Single-thread version for use inside omp parallel for.
/// @param[in,out] panel  Panel whose DL has been modified
/// @param[in]     ft     FluxTable for peneflux lookup
void st_recalc_peneflux_signal_from_DL(
    DetectorPanel& panel, const FluxTable& ft)
{
  const Grid2dXYZ g2flux = ft.get_g2_log_peneflux();
  const int nx = panel.get_nbinx();
  const int ny = panel.get_nbiny();
  for (int iy = 0; iy < ny; ++iy)
    for (int ix = 0; ix < nx; ++ix) {
      DetectorElement& ele = panel.callDetectorElement(ix, iy);
      // Only recompute peneflux (DL-dependent).
      // omega, effective_area, exposure_time are unchanged from Phase 1.
      const double peneflux = ele.calc_peneflux(g2flux);
      ele.set_peneflux(peneflux);
      // signal = peneflux * omega * area * time (omega/area/time already set)
      const double signal = ele.calc_signal();
      ele.set_signal(signal);
    }
}

} // anonymous namespace

// =======================================================
// static JSON utility function definitions
// =======================================================

double DepthResolutionSweep::assign_double(
  const nlohmann::json& j, const char* key)
{
  bool key_found = j.contains(key);
  if (!key_found){ THROW_ERROR2("Key not found in JSON:", key); }
  return j.at(key).get<double>();
}

std::vector<double> DepthResolutionSweep::assign_vec_double(
  const nlohmann::json& j, const char* key)
{
  bool key_found = j.contains(key);
  if (!key_found){ THROW_ERROR2("Key not found in JSON:", key); }
  return j.at(key).get<std::vector<double>>();
}

std::vector<std::array<double, 2>>
  DepthResolutionSweep::assign_vec_double2(
    const nlohmann::json& j, const char* key)
{
  std::vector<std::array<double, 2>> result;

  if (!j.contains(key))
    THROW_ERROR2("Key not found in JSON:", key);

  const auto& arr = j.at(key);
  if (!arr.is_array())
    THROW_ERROR2("Field must be an array of arrays: ", key);

  for (const auto& elem : arr) {
    if (!elem.is_array() || elem.size() != 2)
      THROW_ERROR("Each element must be an array of 2 numbers [x, y].");
    try {
      const double v0 = elem.at(0).get<double>();
      const double v1 = elem.at(1).get<double>();
      result.push_back({v0, v1});
    } catch (const nlohmann::json::type_error& e) {
      THROW_ERROR2("Non-numeric value found in ", key);
    }
  }
  return result;
}

Angle DepthResolutionSweep::assign_angle_rad(
  const nlohmann::json& j, const char* key)
{
  bool key_found = j.contains(key);
  if (!key_found){ THROW_ERROR2("Key not found in JSON:", key); }
  return Angle(j.at(key).get<double>(), Angle::Unit::Radian);
}

double DepthResolutionSweep::assign_double(
  const nlohmann::json* jptr, const char* key)
{
  if (!jptr || !jptr->contains(key)) { THROW_ERROR2("Key not found in JSON:", key); }
  return jptr->at(key).get<double>();
}

Angle DepthResolutionSweep::assign_angle_rad(
  const nlohmann::json* jptr, const char* key)
{
  if (!jptr || !jptr->contains(key)) { THROW_ERROR2("Key not found in JSON:", key); }
  return Angle(jptr->at(key).get<double>(), Angle::Unit::Radian);
}

// =======================================================
// class function definitions
// =======================================================

DepthResolutionSweep::Parameters
  DepthResolutionSweep::Parameters::from_json(
    const nlohmann::json& root, const std::string& section)
{
  Parameters prm;
  prm.section_name = section;

  if (!root.contains(section)) {
    LOG_WARN("JSON section '{}' not found. Using default parameters.", section);
    return prm;
  }

  const auto& js = root.at(section);

  // ---- Detector file loading ----
  // Priority: DEPTH_RESOLUTION_SWEEP.det_files > DETECTOR_PARAMETER_LISTS
  if (js.contains("det_files")) {
    prm.det_param_lists.assign_parameters(root, section);
    LOG_INFO("Using det_files from DEPTH_RESOLUTION_SWEEP section.");
  } else if (root.contains("DETECTOR_PARAMETER_LISTS")) {
    prm.det_param_lists.assign_parameters(root, "DETECTOR_PARAMETER_LISTS");
    LOG_INFO("Using det_files from DETECTOR_PARAMETER_LISTS section (fallback).");
  } else {
    THROW_ERROR("'det_files' not found in DEPTH_RESOLUTION_SWEEP or DETECTOR_PARAMETER_LISTS.");
  }
  if (prm.det_param_lists.get_n_det() <= 0) {
    THROW_ERROR("det_files must contain at least one detector file.");
  }
  LOG_INFO("Detected {} detector files for sweep.",
    prm.det_param_lists.get_n_det());

  // ---- Mode selection (cylinder vs DDL) ----
  if (js.contains("mode")) {
    const auto mode_s = js.at("mode").get<std::string>();
    const auto lower = [](std::string s){ for(auto& c:s) c = static_cast<char>(::tolower(c)); return s; };
    const std::string m = lower(mode_s);
    if (m == "cylinder") {
      prm.mode = Parameters::Mode::Cylinder;
      LOG_DEBUG("mode: cylinder");
    } else if (m == "ddl" || m == "delta density-length approximation"
            || m == "delta-density-length-approximation"
            || m == "baumkuchen" || m == "baumkuchen-approximation"
            || m == "baumkuchen approximation") {
      prm.mode = Parameters::Mode::DDL;
      LOG_DEBUG("mode: DDL (Delta Density-Length Approximation)");
    } else {
      THROW_ERROR2("Unknown mode. Use 'cylinder' or 'ddl'. Given:", mode_s);
    }
  } else {
    prm.mode = Parameters::Mode::Cylinder;
    LOG_DEBUG("mode: (default) cylinder");
  }

  // ---- Retrieve mode-specific subsections ----
  if (!js.contains("common")) THROW_ERROR("'common' section is required.");
  const auto& jcom = js.at("common");
  const nlohmann::json* jmode = nullptr;
  if (prm.mode == Parameters::Mode::Cylinder) {
    if (!js.contains("cylinder")) THROW_ERROR("'cylinder' section is required for cylinder mode.");
    jmode = &js.at("cylinder");
  } else {
    if (!js.contains("ddl")) THROW_ERROR("'ddl' section is required for DDL mode.");
    jmode = &js.at("ddl");
  }

  if (jcom.contains("output_ascii_prefix")) {
    prm.output_ascii_prefix = jcom.at("output_ascii_prefix").get<std::string>();
    LOG_DEBUG("output_ascii_prefix: {}", prm.output_ascii_prefix.string());
  }

  // ---- Common parameters: object center position ----
  LOG_DEBUG("Retrieve obj xy center (from 'common')");
  double x_cnt_obj = 0.0;
  double y_cnt_obj = 0.0;

  // x_cnt_obj
  x_cnt_obj = assign_double(jcom, "x_cnt_obj");

  // y_cnt_obj
  y_cnt_obj = assign_double(jcom, "y_cnt_obj");

  // set obj_center
  prm.obj_center = {x_cnt_obj, y_cnt_obj};
  LOG_DEBUG(" obj_center: x={} m, y={} m", prm.obj_center[0], prm.obj_center[1]);

  // ---- Common parameters: density and size limits ----
  prm.base_density = -1000.0; // default -1000 kg/m3
  prm.base_density = assign_double(jcom, "base_density");
  if( prm.base_density <= 0.0 )
    THROW_ERROR("base_density must be positive.");
  LOG_DEBUG("base_density: {} kg/m3", prm.base_density);

  // obj_size_upper_limit
  prm.obj_size_upper_limit = assign_double(jcom, "obj_size_upper_limit");
  if( prm.obj_size_upper_limit <= 0.0 )
    THROW_ERROR("obj_size_upper_limit must be positive.");
  LOG_DEBUG("obj_size_upper_limit: {} m", prm.obj_size_upper_limit);
  
  // obj_size_lower_limit
  prm.obj_size_lower_limit = assign_double(jcom, "obj_size_lower_limit");
  if( prm.obj_size_lower_limit <= 0.0 )
    THROW_ERROR("obj_size_lower_limit must be positive.");
  LOG_DEBUG("obj_size_lower_limit: {} m", prm.obj_size_lower_limit);

  // obj_size_upper_limit should be > obj_size_lower_limit
  if( prm.obj_size_upper_limit <= prm.obj_size_lower_limit )
    THROW_ERROR("obj_size_upper_limit must be greater than obj_size_lower_limit.");
  LOG_DEBUG(" obj_size_upper_limit > obj_size_lower_limit check passed.");

  // elev_center_step
  prm.elev_center_step = assign_double(jcom, "elev_center_step");
  if (prm.elev_center_step <= 0.0)
    THROW_ERROR("elev_center_step must be positive.");
  LOG_DEBUG("elev_center_step: {:.4f}", prm.elev_center_step);

  // angle_between_cut_factor
  prm.angle_between_cut_factor = assign_double(jcom, "angle_between_cut_factor");
  if (prm.angle_between_cut_factor <= 0.0)
    THROW_ERROR("angle_between_cut_factor must be positive.");
  LOG_DEBUG("angle_between_cut_factor: {}", prm.angle_between_cut_factor);

  // sweep_range_factor
  prm.sweep_range_factor = assign_double(jcom, "sweep_range_factor");
  if (prm.sweep_range_factor <= 0.0)
    THROW_ERROR("sweep_range_factor must be positive.");
  LOG_DEBUG("sweep_range_factor: {}", prm.sweep_range_factor);

  // ---- Common parameters: density variation and statistical thresholds ----
  prm.vec_delta_density = assign_vec_double(jcom, "vec_delta_density");
  if( prm.vec_delta_density.empty() ) 
    THROW_ERROR("vec_delta_density must contain at least one value.");
  const int n_delta_dens = static_cast<int>(prm.vec_delta_density.size());
  LOG_DEBUG("vec_delta_density size: {}", n_delta_dens);
  for (size_t i = 0; i < n_delta_dens; i++)
    LOG_DEBUG("vec_delta_density[{}]: {} kg/m3", i, prm.vec_delta_density[i]);

  // array of stat_alphas
  const auto vec_alpha_json = jcom.at("stat_alphas").get<std::vector<double>>();
  prm.vec_stat_alpha_value.clear();
  for( const auto alpha_val : vec_alpha_json ){
    if( alpha_val <= 0.0 || alpha_val >= 1.0 )
      THROW_ERROR2("Each stat_alpha value must be in (0,1).", alpha_val);
    prm.vec_stat_alpha_value.push_back( alpha_val );
    LOG_DEBUG("alpha value added: {}", alpha_val);
  }

  // both_side: statistical test sidedness (default: false = one-sided)
  if (jcom.contains("both_side")) {
    prm.both_side = jcom.at("both_side").get<bool>();
    LOG_DEBUG("both_side: {}", prm.both_side);
  } else {
    prm.both_side = false;
    LOG_DEBUG("both_side: (default) false (one-sided)");
  }

  // array of signal / noise amplifiers
  if (jcom.contains("signal_noise_amplifiers")) {
    const auto& arr = jcom.at("signal_noise_amplifiers"); // [[s,n], [s,n], ...]
    prm.vec_signal_noise_amplifiers.clear();
    for (const auto& pair_json : arr) {
      if (!pair_json.is_array() || pair_json.size() != 2)
        THROW_ERROR("Each element of signal_noise_amplifiers must be an array of 2 numbers [signal_amp, noise_amp].");

      const double signal_amp = pair_json.at(0).get<double>();
      const double noise_amp  = pair_json.at(1).get<double>();
      if (signal_amp <= 0.0)
        THROW_ERROR2("signal_amplifier must be positive.", signal_amp);
      if (noise_amp < 0.0)
        THROW_ERROR2("noise_amplifier must be non-negative.", noise_amp);

      prm.vec_signal_noise_amplifiers.push_back({signal_amp, noise_amp});
      LOG_DEBUG("signal_noise_amplifier added: signal={}, noise={}", signal_amp, noise_amp);
    }
  }

  // --- DDL mode specific ---
  if (prm.mode == Parameters::Mode::DDL) {
    prm.obj_size_step = assign_double(jmode, "obj_size_step");
    if (prm.obj_size_step <= 0.0)
      THROW_ERROR("obj_size_step must be positive.");
    LOG_DEBUG("obj_size_step: {}", prm.obj_size_step);

    prm.depth_step = assign_double(jmode, "depth_step");
    if (prm.depth_step <= 0.0)
      THROW_ERROR("depth_step must be positive (DDL).");
    LOG_DEBUG("depth_step (DDL): {}", prm.depth_step);

    // Additional consistency check
    if (prm.obj_size_step > (prm.obj_size_upper_limit - prm.obj_size_lower_limit))
      THROW_ERROR("For DDL, obj_size_step must be <= (upper - lower).");
  }

  // --- Cylinder mode specific ---
  if (prm.mode == Parameters::Mode::Cylinder) {
  }

  // --- Angle bin override (optional) ---
  prm.tf_override_angle_bin = false;
  if (jcom.contains("tf_override_angle_bin")) {
    prm.tf_override_angle_bin = jcom.at("tf_override_angle_bin").get<bool>();
  }
  LOG_DEBUG("tf_override_angle_bin: {}", prm.tf_override_angle_bin);

  if (prm.tf_override_angle_bin) {
    if (!jcom.contains("angle_bin_override")) {
      THROW_ERROR("'angle_bin_override' section is required when tf_override_angle_bin is true.");
    }
    const auto& jabo = jcom.at("angle_bin_override");

    prm.angle_bin_override.nbinx = jabo.at("nbinx").get<int>();
    prm.angle_bin_override.txmin = jabo.at("txmin").get<double>();
    prm.angle_bin_override.txmax = jabo.at("txmax").get<double>();
    prm.angle_bin_override.nbiny = jabo.at("nbiny").get<int>();
    prm.angle_bin_override.tymin = jabo.at("tymin").get<double>();
    prm.angle_bin_override.tymax = jabo.at("tymax").get<double>();

    if (prm.angle_bin_override.nbinx <= 0)
      THROW_ERROR("angle_bin_override.nbinx must be positive.");
    if (prm.angle_bin_override.nbiny <= 0)
      THROW_ERROR("angle_bin_override.nbiny must be positive.");
    if (prm.angle_bin_override.txmax <= prm.angle_bin_override.txmin)
      THROW_ERROR("angle_bin_override.txmax must be greater than txmin.");
    if (prm.angle_bin_override.tymax <= prm.angle_bin_override.tymin)
      THROW_ERROR("angle_bin_override.tymax must be greater than tymin.");

    LOG_INFO("Angle bin override enabled: nbinx={}, txmin={}, txmax={}, nbiny={}, tymin={}, tymax={}",
      prm.angle_bin_override.nbinx, prm.angle_bin_override.txmin, prm.angle_bin_override.txmax,
      prm.angle_bin_override.nbiny, prm.angle_bin_override.tymin, prm.angle_bin_override.tymax);
  }

  // --- Base distribution output (optional, default: false) ---
  prm.tf_out_det_PL_signal = false;
  if (jcom.contains("tf_out_det_PL_signal")) {
    prm.tf_out_det_PL_signal = jcom.at("tf_out_det_PL_signal").get<bool>();
  }
  LOG_DEBUG("tf_out_det_PL_signal: {}", prm.tf_out_det_PL_signal);

  return prm;
}

DepthResolutionSweep::DepthResolutionSweep(
  const DepthResolutionSweep::Parameters& params
, const FluxTable& ft_in
, const Grid2dPillar& g2pil_base
, const DetectorPanelArray& arrdet_in)
  : params_(params)
  , ft_(ft_in)
  , g2pil_base_(g2pil_base)
  , arrdet_base_(arrdet_in)
{
  LOG_INFO("DepthResolutionSweep initialized {} delta densities."
    , params_.vec_delta_density.size());
  LOG_INFO("DepthResolutionSweep received DetectorPanelArray with {} detectors.",
    arrdet_base_.get_n_det());
}

DepthResolutionSweep::~DepthResolutionSweep() = default;

/// @brief Calculate and set vec_vec_angbetween from arrdet_base_ and g2pil_base_ at params_.obj_center
void DepthResolutionSweep::set_vec_vec_angbetween( const double obj_size_upper_limit )
{
  const Eigen::Vector2d v2_pos_obj = params_.obj_center_vec();
  LOG_DEBUG("Calculating vec_vec_alpha based on obj_center ({}, {})",
      v2_pos_obj.x(), v2_pos_obj.y());
  
  LOG_DEBUG("Num of DetectorPanel is {}", arrdet_base_.get_n_det());
  vec_vec_angle_between_from_obj_center_.clear();
  
  const int n_det = arrdet_base_.get_n_det();
  vec_vec_angle_between_from_obj_center_.resize(n_det);

  for (int detid = 0; detid < n_det; ++detid) {
    LOG_DEBUG("Processing detid={}", detid);
    const DetectorPanel& panel = arrdet_base_.getDetectorPanel(detid);
    const Ray3d& ray3d = panel.get_ray3d();
    const double x_angle_min = panel.get_x_interval();
    const angle_util::AngleUnit angle_unit = panel.get_angle_unit();
    const std::string angle_unit_str = angle_util::to_string(angle_unit);
    LOG_DEBUG(" detid={} | angle_unit={}, x_angle_min={}", detid, angle_unit_str, x_angle_min);

    // For each detector, compute the maximum angular width (ang_between_max)
    // corresponding to the largest cylinder radius that fits within obj_size_upper_limit.
    // This uses the horizontal distance from detector to object center.
    const double radius_obj = 0.5 * obj_size_upper_limit;
    const Angle ang_between_max = geom_util::calc_horizontal_angle_from_cylinder_radius(
      ray3d, v2_pos_obj, radius_obj);
    LOG_DEBUG("  detid={} | ang_between_max = {:.6f} rad ({:.4f} deg)"
    , detid, ang_between_max.rad(), ang_between_max.deg());
    
    // Iterate through angular multiples of the detector's x-axis bin width.
    // Each ang_between corresponds to a cylinder radius that subtends that angle
    // as seen from the detector position.
    const int nbinx_half = panel.get_nbinx() / 2;
    LOG_DEBUG("  detid={} | nbinx_half = {}", detid, nbinx_half);

    int n_ang_between = 1;
    Angle ang_between(0.0, Angle::Unit::Radian);
    while(ang_between < ang_between_max &&  n_ang_between < nbinx_half) {
      const double x_angle = n_ang_between * x_angle_min;
      // Convert x_angle --> ang_between
      if      ( angle_unit == angle_util::AngleUnit::Tangent ) ang_between.setTangent(x_angle);
      else if ( angle_unit == angle_util::AngleUnit::Radian  ) ang_between.setRadian( x_angle);
      else if ( angle_unit == angle_util::AngleUnit::Degree  ) ang_between.setDegree( x_angle);
      else{ THROW_ERROR("Unknown angle unit."); }
      vec_vec_angle_between_from_obj_center_.at(detid).push_back(ang_between);
      LOG_DEBUG(" detid={}, | n_ang_between={}, x_angle={:.6f} (angle_unit={}), ang_between={:.6f} rad ({:.4f} deg)",
        detid, n_ang_between, x_angle, angle_unit_str, ang_between.rad(), ang_between.deg());
      n_ang_between++;
    } // while ang_between < ang_between_max, n_ang_between++
  } // for detid
}

void DepthResolutionSweep::run_mode_cylinder()
{
  LOG_INFO("....... ");
  LOG_INFO("  output path prefix: {}", params_.output_ascii_prefix.string());
  const Eigen::Vector2d v2_pos_obj = params_.obj_center_vec();
  LOG_INFO("  cyl center : ({}, {})", v2_pos_obj.x(), v2_pos_obj.y());

  // Set vec_vec_angle_between_from_obj_center_
  set_vec_vec_angbetween(params_.obj_size_upper_limit);

  // Make Grid2dPillar base as copy of g2pil_base_
  Grid2dPillar g2pil_base = this->g2pil_base_;

  // Set 1st uniform density of g2pil_base
  g2pil_base.set_uniform_density( params_.base_density );

  const int n_det = arrdet_base_.get_n_det();
  for (int detid = 0; detid < n_det; ++detid)
    LOG_INFO(" detid={}, size of vec_ang_between = {}"
      , detid, vec_vec_angle_between_from_obj_center_.at(detid).size());
  
  const int n_dens = static_cast<int>(params_.vec_delta_density.size());
  for( int i=0; i<n_dens; ++i )
    LOG_INFO("  delta dens : {}", params_.vec_delta_density.at(i));

  const int n_stat_alpha = static_cast<int>(params_.vec_stat_alpha_value.size());
  for( int i=0; i<n_stat_alpha; ++i )
    LOG_INFO("  statistical alpha : {}", params_.vec_stat_alpha_value.at(i));

  // Cylinder rotation angle is 0
  const Angle cyl_angle_rot(0.0, Angle::Unit::Radian);

  // Loop of detids
  for(int detid=0; detid<n_det; ++detid) {
    DetectorPanel panel = arrdet_base_.getDetectorPanelCopy(detid);
    panel.set_direction_to_v2_pos( v2_pos_obj );
    const Ray3d& ray3d = panel.get_ray3d();
    const double elev_det = ray3d.z();

    // vec_signal_noise_stat_result
    std::vector<SignalNoiseStatResult> vec_signal_stat_result;

    const auto& vec_ang_between = vec_vec_angle_between_from_obj_center_.at(detid);

    // ================================================================
    // Loop over angular sizes (spatial resolution candidates)
    // Each ang_between corresponds to a cylinder diameter (obj_size)
    // ================================================================
    for( const auto& ang_between : vec_ang_between ) {
      LOG_DEBUG(" detid={}, ang_between={:.6f} rad ({:.4f} deg)"
        , detid, ang_between.rad(), ang_between.deg());

      // Convert angular width to physical cylinder radius using the
      // detector-to-object geometry. r_cyl = distance * tan(ang_between/2) approximately.
      const double r_cyl = geom_util::calc_cylinder_radius_from_horizontal_angle_size(
        ray3d, v2_pos_obj, ang_between);
      LOG_DEBUG("  r_cyl = {} meters", r_cyl);
      const double obj_size = r_cyl * 2.0;

      // Skip objects smaller than the minimum detectable size
      if( obj_size < params_.obj_size_lower_limit ) {
        LOG_DEBUG("  obj_size {} m < obj_size_lower_limit {} m, skip this ang_between.",
          obj_size, params_.obj_size_lower_limit);
        continue;
      }

      // Set cylinder length vector
      Eigen::Vector2d v2_length;
      v2_length.x() = r_cyl; // diameter
      v2_length.y() = r_cyl; // diameter

      // Convert ang_between to detector's angle unit
      double ang_between_unit_value = 0.0;
      if( panel.get_angle_unit() == angle_util::AngleUnit::Tangent ) {
        ang_between_unit_value = std::tan( fabs(ang_between.rad()) );
        LOG_DEBUG("  Converted ang_between to tangent: {:.6f}", ang_between_unit_value);
      }
      else if( panel.get_angle_unit() == angle_util::AngleUnit::Radian ) {
        ang_between_unit_value = ang_between.rad();
        LOG_DEBUG("  Converted ang_between to radian: {:.6f}", ang_between_unit_value);
      }
      else if( panel.get_angle_unit() == angle_util::AngleUnit::Degree ) {
        ang_between_unit_value = ang_between.deg();
        LOG_DEBUG("  Converted ang_between to degree: {:.6f}", ang_between_unit_value);
      }
      else{
        THROW_ERROR("Unknown angle unit.");
      }

      // Get zmax of object
      const double elev_obj = g2pil_base.getMinimumZmaxCircle( 
        v2_pos_obj.x(), v2_pos_obj.y(), r_cyl );

      // Get the top of terrain
      const double elev_terrain = g2pil_base.getPillar(
        v2_pos_obj.x(), v2_pos_obj.y() ).get_zmax();

      // get the cylinder height
      const double elev_obj_from_det = elev_obj - elev_det;
      LOG_DEBUG(" zmax of g2pil_base at obj_center: {:.1f} meters", elev_obj);
      LOG_DEBUG(" height of cyl: {:.1f} meters", elev_obj_from_det);

      // get the depth from the terrain surface to the bottom of detector
      const double elev_terrain_from_det = elev_terrain - elev_det;
      LOG_DEBUG(" elev_terrain_from_det (terrain surface to bottom of detector): {:.1f} meters", elev_terrain_from_det);
      if (obj_size > elev_terrain_from_det) { continue; }

      // get the diff between terrain surface and obj zmax
      const double diff_elev = elev_terrain - elev_obj;

      // Position of cylinder top (at terrain surface)
      const Eigen::Vector3d v3_pos_obj_top{ v2_pos_obj.x(), v2_pos_obj.y(), elev_obj };

      const double tymin = panel.get_g2bg().get_ymin();
      const double tymax = panel.get_g2bg().get_ymax();

      // Define horizontal (tx) cut range: only bins within +-cut_factor * ang_between
      // contribute to the measurement window for this object size.
      const double cut_factor = params_.angle_between_cut_factor;
      const double tx_cut_min = -ang_between_unit_value * cut_factor;
      const double tx_cut_max =  ang_between_unit_value * cut_factor;

      // ---- BASELINE: Compute path lengths through terrain without anomaly ----
      DetectorPanel panel_cut_base = panel.cut(tx_cut_min, tx_cut_max, tymin, tymax);
      pathcalc::g2pil::mp_add_PLDL(panel_cut_base,g2pil_base);
      panel_cut_base.mp_calc_set_peneflux_signal_from_DL(ft_,false);
      const Eigen::Vector2d v2_pos_det = panel.get_ray3d().pos2d();

      // Sweep DetectorPanelArray over angular range
      const double sweep_range_factor = params_.sweep_range_factor;
      const double tx_lower = -ang_between_unit_value*sweep_range_factor;
      const double tx_upper =  ang_between_unit_value*sweep_range_factor;
      const double ty_lower = -ang_between_unit_value*sweep_range_factor;
      const double ty_upper =  ang_between_unit_value*sweep_range_factor;

      // Use params_.elev_center_step for ty_step
      const double ty_step = params_.elev_center_step;

      LOG_DEBUG("  Sweeping DetectorPanel detid={} | tx: [{:.6f}, {:.6f}], ty: [{:.6f}, {:.6f}], ty_step: {:.6f}",
        detid, tx_lower, tx_upper, ty_lower, ty_upper, ty_step);
      std::vector<SignalNoiseDepth> vec_snd_base =
        panel_cut_base.mp_get_signal_noise_sum_y_sweep(
          tx_lower, tx_upper, ty_lower, ty_upper, tymin, tymax, ty_step, v3_pos_obj_top, diff_elev);

      // ---- MODIFIED: Add cylindrical density anomaly and recompute ----
      for( const double& delta_dens : params_.vec_delta_density ) {
        LOG_DEBUG("  delta_dens = {} kg/m3", delta_dens);

        // Create cylinder object with density difference delta_dens
        const VerticalEllipticCylinderCapped
           vcyl( v2_pos_obj, elev_det, elev_obj_from_det
            , r_cyl, r_cyl, Angle(0.0, Angle::Unit::Radian), delta_dens );

        // Copy baseline panel for modification
        DetectorPanel panel_cut_mod(panel_cut_base);

        // pathcalc::vcyl::mp_add_DL computes the exact DL contribution from
        // a vertical elliptic cylinder using geometric intersection calculations.
        pathcalc::vcyl::mp_add_DL(panel_cut_mod,vcyl);

        // Calculate signal, noise
        panel_cut_mod.mp_calc_set_peneflux_signal_from_DL(ft_,false);

        LOG_DEBUG("  Sweeping DetectorPanel detid={} | tx: [{:.6f}, {:.6f}], ty: [{:.6f}, {:.6f}], ty_step: {:.6f}",
          detid, tx_lower, tx_upper, ty_lower, ty_upper, ty_step);
        std::vector<SignalNoiseDepth> vec_snd_modi =
          panel_cut_mod.mp_get_signal_noise_sum_y_sweep(
            tx_lower, tx_upper, ty_lower, ty_upper, tymin, tymax, ty_step, v3_pos_obj_top, diff_elev);

        // loop of signal/noise amplifier pairs
        for( const auto& [signal_amp, noise_amp] : params_.vec_signal_noise_amplifiers ){
          LOG_DEBUG("    Signal amplifier = {}", signal_amp);
          LOG_DEBUG("    Noise amplifier = {}", noise_amp);

          // Loop params_.vec_stat_alpha_value to perform and output significance tests
          for( const double stat_alpha_value : params_.vec_stat_alpha_value ){
            LOG_DEBUG("   Statistical alpha value = {}", stat_alpha_value);
              signal_noise_stat_util::eval_signal_significance(
                  vec_snd_base, vec_snd_modi
                , obj_size, delta_dens, stat_alpha_value, signal_amp
                , vec_signal_stat_result
                , params_.both_side); // true=two-sided, false(default)=one-sided

          } // for stat_alpha_value
        } // for signal/noise amplifier pairs
      } // for delta_dense
    } // for ang_between

    // output vec_signal_stat_result to file
    char cfnameout_stat[512];
    sprintf(cfnameout_stat, "%s_det%02d_signal_signifi.tmp"
      , params_.output_ascii_prefix.string().c_str(), detid);
    fs::path output_path_stat = iodir::make_pathout(cfnameout_stat);
    LOG_INFO(" Outputting signal significance results to {}"
      , output_path_stat.string());
    signal_noise_stat_util::out_signal_stat_result_vector(
      vec_signal_stat_result, output_path_stat);

    char cfnameout_stat_csv[512];
    sprintf(cfnameout_stat_csv, "%s_det%02d_signal_signifi-tmp.csv"
      , params_.output_ascii_prefix.string().c_str(), detid);
    fs::path output_path_stat_csv = iodir::make_pathout(cfnameout_stat_csv);
    LOG_INFO(" Outputting signal significance results to {}"
      , output_path_stat_csv.string());
    signal_noise_stat_util::out_signal_stat_result_vector_csv(
      vec_signal_stat_result, output_path_stat_csv);

  } // for detid  
}

// Execute sweep in DDL mode
void DepthResolutionSweep::run_mode_ddl()
{
  LOG_INFO("....... ");
  LOG_INFO("  output path prefix: {}", params_.output_ascii_prefix.string());
  const Eigen::Vector2d v2_pos_obj = params_.obj_center_vec();
  LOG_INFO("  obj center : ({}, {})", v2_pos_obj.x(), v2_pos_obj.y());

  // Set vec_vec_angle_between_from_obj_center_
  set_vec_vec_angbetween(params_.obj_size_upper_limit);

  // Make Grid2dPillar base as copy of g2pil_base_
  Grid2dPillar g2pil_base = this->g2pil_base_;

  // Set 1st uniform density of g2pil_base
  g2pil_base.set_uniform_density( params_.base_density );

  const int n_det = arrdet_base_.get_n_det();
  for (int detid = 0; detid < n_det; ++detid)
    LOG_INFO(" detid={}, size of vec_ang_between = {}", detid, vec_vec_angle_between_from_obj_center_.at(detid).size());
  
  const int n_dens = static_cast<int>(params_.vec_delta_density.size());
  for( int i=0; i<n_dens; ++i )
    LOG_INFO("  delta dens : {}", params_.vec_delta_density.at(i));

  const int n_stat_alpha = static_cast<int>(params_.vec_stat_alpha_value.size());
  for( int i=0; i<n_stat_alpha; ++i )
    LOG_INFO("  statistical alpha : {}", params_.vec_stat_alpha_value.at(i));

  // Cylinder rotation angle is 0
  const Angle cyl_angle_rot(0.0, Angle::Unit::Radian);

  // Loop of detids
  for(int detid=0; detid<n_det; ++detid) {
    DetectorPanel panel = arrdet_base_.getDetectorPanelCopy(detid);
    panel.set_direction_to_v2_pos( v2_pos_obj );
    const Ray3d& ray3d = panel.get_ray3d();
    const double elev_det = ray3d.z();

    // vec_signal_noise_stat_result
    std::vector<SignalNoiseStatResult> vec_signal_stat_result;

    const auto& vec_ang_between = vec_vec_angle_between_from_obj_center_.at(detid);

    // loop of alpha = obj_size
    for( const auto& ang_between : vec_ang_between ) {
      LOG_DEBUG_ND(" detid={}, ang_between={:.6f} rad ({:.4f} deg)"
        , detid, ang_between.rad(), ang_between.deg());

      // Calculate the radius of the tangent circle for each ang_between
      const double r_cyl = geom_util::calc_cylinder_radius_from_horizontal_angle_size(
        ray3d, v2_pos_obj, ang_between);
      LOG_DEBUG_ND("  r_cyl = {} meters", r_cyl);
      const double obj_size = r_cyl * 2.0;

      // if obj_size < obj_size_lower_limit, skip
      if( obj_size < params_.obj_size_lower_limit ) {
        LOG_DEBUG_ND("  obj_size {} m < obj_size_lower_limit {} m, skip this ang_between.",
          obj_size, params_.obj_size_lower_limit);
        continue;
      }

      // Convert ang_between to detector's angle unit
      double ang_between_unit_value = 0.0;
      if( panel.get_angle_unit() == angle_util::AngleUnit::Tangent ) {
        ang_between_unit_value = std::tan( fabs(ang_between.rad()) );
        LOG_DEBUG_ND("  Converted ang_between to tangent: {:.6f}", ang_between_unit_value);
      }
      else if( panel.get_angle_unit() == angle_util::AngleUnit::Radian ) {
        ang_between_unit_value = ang_between.rad();
        LOG_DEBUG_ND("  Converted ang_between to radian: {:.6f}", ang_between_unit_value);
      }
      else if( panel.get_angle_unit() == angle_util::AngleUnit::Degree ) {
        ang_between_unit_value = ang_between.deg();
        LOG_DEBUG_ND("  Converted ang_between to degree: {:.6f}", ang_between_unit_value);
      }
      else{
        THROW_ERROR("Unknown angle unit.");
      }
      
      // Get zmax of object
      const double elev_obj = g2pil_base.getMinimumZmaxCircle( 
        v2_pos_obj.x(), v2_pos_obj.y(), r_cyl );

      // Get the top of terrain
      const double elev_terrain = g2pil_base.getPillar(
        v2_pos_obj.x(), v2_pos_obj.y() ).get_zmax();

      // get the cylinder height
      const double elev_obj_from_det = elev_obj - elev_det;
      LOG_DEBUG(" zmax of g2pil_base at obj_center: {:.1f} meters", elev_obj);
      LOG_DEBUG(" height of cyl: {:.1f} meters", elev_obj_from_det);

      // get the depth from the terrain surface to the bottom of detector
      const double elev_terrain_from_det = elev_terrain - elev_det;
      LOG_DEBUG(" elev_terrain_from_det (terrain surface to bottom of detector): {:.1f} meters", elev_terrain_from_det);
      if (obj_size > elev_terrain_from_det) { continue; }

      // get the diff between terrain surface and obj zmax
      const double diff_elev = elev_terrain - elev_obj;

      // position of obj top
      const Eigen::Vector3d v3_pos_obj_top{ v2_pos_obj.x(), v2_pos_obj.y(), elev_obj };

      const double tymin = panel.get_g2bg().get_ymin();
      const double tymax_panel = panel.get_g2bg().get_ymax();
      // Limit tymax to avoid negative DL (where obj_size > PL)
      const double tymax_allowed = calc_ty_max_allowed(
        elev_obj_from_det, obj_size, panel.get_angle_unit());
      const double tymax = std::min(tymax_panel, tymax_allowed);

      // Check if valid ty range exists
      if (tymax <= tymin) {
        LOG_DEBUG_ND("  obj_size {} m exceeds PL at all elevations (elev_obj_from_det={} m), skip.",
          obj_size, elev_obj_from_det);
        continue;
      }
      if (tymax_allowed < tymax_panel) {
        LOG_DEBUG_ND("  tymax limited: panel={:.6f}, allowed={:.6f}, used={:.6f}",
          tymax_panel, tymax_allowed, tymax);
      }

      // panel cut range
      const double cut_factor = params_.angle_between_cut_factor;
      const double tx_cut_min = -ang_between_unit_value * cut_factor;
      const double tx_cut_max =  ang_between_unit_value * cut_factor;

      DetectorPanel panel_cut_base = panel.cut(tx_cut_min, tx_cut_max, tymin, tymax);
      pathcalc::g2pil::mp_add_PLDL(panel_cut_base,g2pil_base);
      panel_cut_base.mp_calc_set_peneflux_signal_from_DL(ft_,false);
      const Eigen::Vector2d v2_pos_det = panel.get_ray3d().pos2d();

      // Sweep DetectorPanelArray over angular range
      const double sweep_range_factor = params_.sweep_range_factor;
      const double tx_lower = -ang_between_unit_value*sweep_range_factor;
      const double tx_upper =  ang_between_unit_value*sweep_range_factor;
      const double ty_lower = -ang_between_unit_value*sweep_range_factor;
      const double ty_upper =  ang_between_unit_value*sweep_range_factor;

      // Use params_.elev_center_step for ty_step
      const double ty_step = params_.elev_center_step;

      LOG_DEBUG_ND("  Sweeping DetectorPanel detid={} | tx: [{:.6f}, {:.6f}], ty: [{:.6f}, {:.6f}], ty_step: {:.6f}",
        detid, tx_lower, tx_upper, ty_lower, ty_upper, ty_step);
      std::vector<SignalNoiseDepth> vec_snd_base =
        panel_cut_base.mp_get_signal_noise_sum_y_sweep(
          tx_lower, tx_upper, ty_lower, ty_upper, tymin, tymax, ty_step, v3_pos_obj_top, diff_elev);

      // ---- MODIFIED: DDL approximation (Baumkuchen model) ----
      for( const double& delta_dens : params_.vec_delta_density ) {
        LOG_DEBUG_ND("  delta_dens = {} kg/m3", delta_dens);
        // DDL approximation: delta_DL = delta_dens * obj_size (uniform across all bins)
        // This assumes the anomaly is a concentric ring ("Baumkuchen") structure
        // where every ray passes through the same effective thickness.
        const double delta_DL_obj = delta_dens * obj_size;

        // Copy baseline panel for modification
        DetectorPanel panel_cut_mod(panel_cut_base);

        // Add uniform delta_DL to all bins (approximation replaces exact geometry)
        pathcalc::naive::mp_add_DL(panel_cut_mod, delta_DL_obj);

        // Clamp negative DL to zero: can occur when obj_size > path_length at
        // high elevation angles where the ray exits the terrain before traversing
        // the full anomaly diameter.
        mp_clamp_DL_nonnegative(panel_cut_mod);

        // Calculate signal, noise
        panel_cut_mod.mp_calc_set_peneflux_signal_from_DL(ft_,false);

        LOG_DEBUG("  Sweeping DetectorPanel detid={} | tx: [{:.6f}, {:.6f}], ty: [{:.6f}, {:.6f}], ty_step: {:.6f}",
          detid, tx_lower, tx_upper, ty_lower, ty_upper, ty_step);
        std::vector<SignalNoiseDepth> vec_snd_modi =
          panel_cut_mod.mp_get_signal_noise_sum_y_sweep(
            tx_lower, tx_upper, ty_lower, ty_upper, tymin, tymax, ty_step, v3_pos_obj_top, diff_elev);

        // loop of signal/noise amplifier pairs
        for( const auto& [signal_amp, noise_amp] : params_.vec_signal_noise_amplifiers ){
          LOG_DEBUG_ND("    Signal amplifier = {}", signal_amp);
          LOG_DEBUG_ND("    Noise amplifier = {}", noise_amp);

          // Loop params_.vec_stat_alpha_value to perform and output significance tests
          for( const double stat_alpha_value : params_.vec_stat_alpha_value ){
            LOG_DEBUG_ND("   Statistical alpha value = {}", stat_alpha_value);
              signal_noise_stat_util::eval_signal_significance(
                  vec_snd_base, vec_snd_modi
                , obj_size, delta_dens, stat_alpha_value, signal_amp
                , vec_signal_stat_result
                , params_.both_side); // true=two-sided, false(default)=one-sided

          } // for stat_alpha_value
        } // for signal/noise amplifier pairs
      } // for delta_dense
    } // for ang_between

    // define sort condition
    std::vector<SignalNoiseStatSortCondition> vec_sort_condition = {
        {SignalNoiseStatSortKey::StatAlpha,   false}  // descending
      , {SignalNoiseStatSortKey::SignalAmp,   true}  // ascending
      , {SignalNoiseStatSortKey::DeltaDens,   false} // descending
      , {SignalNoiseStatSortKey::ObjSize,     true}  // ascending
      , {SignalNoiseStatSortKey::Depth_lower, true}  // ascending
    };

    // sort vec_signal_stat_result
    signal_noise_stat_util::sort_signal_stat_results(
      vec_signal_stat_result, vec_sort_condition);

    // output vec_signal_stat_result to file
    char cfnameout_stat[512];
    sprintf(cfnameout_stat, "%s_det%02d_signal_signifi.tmp"
      , params_.output_ascii_prefix.string().c_str(), detid);
    fs::path output_path_stat = iodir::make_pathout(cfnameout_stat);
    LOG_INFO(" Outputting signal significance results to {}"
      , output_path_stat.string());
    signal_noise_stat_util::out_signal_stat_result_vector(
      vec_signal_stat_result, output_path_stat);

    char cfnameout_stat_csv[512];
    sprintf(cfnameout_stat_csv, "%s_det%02d_signal_signifi-tmp.csv"
      , params_.output_ascii_prefix.string().c_str(), detid);
    fs::path output_path_stat_csv = iodir::make_pathout(cfnameout_stat_csv);
    LOG_INFO(" Outputting signal significance results to {}"
      , output_path_stat_csv.string());
    signal_noise_stat_util::out_signal_stat_result_vector_csv(
      vec_signal_stat_result, output_path_stat_csv);

  } // for detid
}

// =====================================================================
//
//   run_mode_ddl_v1() — Optimised DDL sweep
//
//   Key algorithmic differences from run_mode_ddl():
//
//   ┌──────────────────────────────────────────────────────────────────┐
//   │  OLD (run_mode_ddl):                                            │
//   │    for each ang_between:                                        │
//   │      panel.cut()                          ← new panel object    │
//   │      pathcalc::g2pil::mp_add_PLDL()       ← full ray-trace     │
//   │      mp_calc_set_peneflux_signal_from_DL() ← flux lookup        │
//   │      mp_get_signal_noise_sum_y_sweep()     ← O(steps × bins)    │
//   │      for each delta_dens:                                       │
//   │        copy panel + add DL + flux + y-sweep again               │
//   │                                                                 │
//   │  → Heavy functions called 160×/detector  (≈ 1 hour total)       │
//   ├──────────────────────────────────────────────────────────────────┤
//   │  NEW (run_mode_ddl_v1):                                         │
//   │    ONCE per detector:                                           │
//   │      pathcalc::g2pil::mp_add_PLDL()  on FULL panel  ← 1×       │
//   │      baseline flux + PrefixSum2D                     ← 1×       │
//   │    for each ang_between:                                        │
//   │      geometry only  (r_cyl, height, tymax)           ← cheap    │
//   │      y-sweep via prefix-sum query                    ← O(steps) │
//   │      for each delta_dens:                                       │
//   │        copy full panel + add DL + flux + PrefixSum2D            │
//   │        y-sweep via prefix-sum query                  ← O(steps) │
//   │                                                                 │
//   │  → mp_add_PLDL: 160× → 1×    (biggest saving)                  │
//   │  → panel.cut():  eliminated   (use full panel + index range)    │
//   │  → y-sweep sum:  O(1)/step    (prefix sum replaces 2D loop)     │
//   └──────────────────────────────────────────────────────────────────┘
//
// =====================================================================

namespace {

/// @brief Perform y-sweep using precomputed PrefixSum2D tables.
///
/// Replaces DetectorPanel::mp_get_signal_noise_sum_y_sweep() with O(1)
/// prefix-sum queries per step instead of iterating over all (ix, iy) bins.
///
/// The sweep iterates from high elevation (near surface, shallowest depth)
/// downward (deepest depth), computing signal/noise sums for rectangular
/// windows at each elevation step. Each step produces a SignalNoiseDepth
/// entry with elev_from_det_lower and elev_from_det_upper computed from the geometry:
///   depth_obj = elev_obj_from_det - tan(ty) * delta_horizontal
///   depth_out = depth_obj + diff_elev   (shifted to terrain-surface reference)
///   where elev_obj_from_det = elev_obj - elev_det
///
/// Steps where depth_obj < 0 are skipped.
/// This keeps the physical sweep limit in object-top reference while depth output
/// remains terrain-surface based (depth_out = depth_obj + diff_elev), so the
/// y=x cutoff is shifted by diff_elev when terrain is uneven.
///
/// @param[in] ps_signal   Precomputed signal prefix-sum table
/// @param[in] ps_noise    Precomputed noise  prefix-sum table
/// @param[in] g2bg        Bin group for coordinate-to-index conversion
/// @param[in] detid       Detector ID for result tagging
/// @param[in] tx_lower    Sweep window x left bound (angle unit)
/// @param[in] tx_upper    Sweep window x right bound (angle unit)
/// @param[in] ty_lower    Sweep window half-width, low side (angle unit)
/// @param[in] ty_upper    Sweep window half-width, high side (angle unit)
/// @param[in] ty_min      Panel lower y limit (angle unit)
/// @param[in] ty_max      Panel upper y limit, may be clamped by calc_ty_max_allowed
/// @param[in] ty_step     Elevation step between successive sweeps (angle unit)
/// @param[in] v3_pos_obj_top  3D position (x, y, z_surface) of anomaly top [m]
/// @param[in] v3_pos_det      3D position (x, y, z) of detector [m]
/// @param[in] angle_unit  Angle unit used for tx/ty values
/// @param[in] diff_elev   Elevation difference between terrain surface and object top [m].
///                        Output depth values are shifted by +diff_elev so that they represent
///                        depth from terrain surface (elev_terrain) rather than object top (elev_obj).
///                        Skip condition uses object-top depth: (elev_from_det_lower_obj < 0).
/// @return Vector of SignalNoiseDepth for each valid elevation step
///
/// @note Complexity: O(nsteps) where nsteps = floor((ty_max - ty_min - ty_len) / ty_step)
/// @note Units: all positions in meters, angles in the unit specified by angle_unit
/// @throws std::runtime_error If ty_step <= 0 or ty_max < ty_min
std::vector<SignalNoiseDepth> sweep_y_with_prefix_sum(
    const PrefixSum2D& ps_signal          // precomputed signal prefix sum
  , const PrefixSum2D& ps_noise           // precomputed noise  prefix sum
  , const Grid2dBinGroup& g2bg            // bin group (for coord→index)
  , const Detid detid                     // detector ID for result tagging
  , const double tx_lower                 // sweep window: x left  bound
  , const double tx_upper                 //                x right bound
  , const double ty_lower                 // sweep window half-width (low)
  , const double ty_upper                 // sweep window half-width (high)
  , const double ty_min                   // panel lower y limit
  , const double ty_max                   // panel upper y limit (may be clamped)
  , const double ty_step                  // elevation step between sweeps
  , const Eigen::Vector3d& v3_pos_obj_top // (x,y,z_surface) of anomaly top
  , const Eigen::Vector3d& v3_pos_det     // detector (x,y,z)
  , const angle_util::AngleUnit angle_unit
  , const double diff_elev)               // elev_terrain - elev_obj [m]
{
  std::vector<SignalNoiseDepth> vec_result;

  // --- Geometric quantities for depth calculation ---
  // depth = elev_obj_from_det - tan(ty) * delta_horizontal
  // where elev_obj_from_det = elev_obj - elev_det
  const Eigen::Vector2d v2_pos_obj = v3_pos_obj_top.head<2>();
  const double elev_obj = v3_pos_obj_top.z();
  const Eigen::Vector2d v2_pos_det = v3_pos_det.head<2>();
  const double elev_det = v3_pos_det.z();
  const Eigen::Vector2d v2_delta = v2_pos_obj - v2_pos_det;
  const double delta_horizontal = v2_delta.norm();
  const double elev_obj_from_det = elev_obj - elev_det;

  if (ty_step <= 0.0) THROW_ERROR("DepthResolutionSweep::run_mode_ddl_v1: ty_step must be positive");
  if (ty_max < ty_min) THROW_ERROR("DepthResolutionSweep::run_mode_ddl_v1: ty_max must be >= ty_min");

  // --- Number of sweep steps ---
  // Sweep starts at ty_max (highest elevation = shallowest depth) and moves down.
  const double ty_len = ty_upper - ty_lower;
  const double sweep_y_length = (ty_max - 0.5 * ty_len) - (ty_min + 0.5 * ty_len);
  int nsteps = static_cast<int>(std::floor(sweep_y_length / ty_step));
  if (nsteps <= 0) nsteps = 1;
  vec_result.reserve(nsteps);

  // --- Convert fixed tx bounds to bin indices (same for every step) ---
  const double eps = 1.0e-6;
  const double tx_eps = g2bg.get_x_interval() * eps;
  const double ty_eps = g2bg.get_y_interval() * eps;
  const int ix_lo = g2bg.get_ix(tx_lower + tx_eps);
  const int ix_hi = g2bg.get_ix(tx_upper - tx_eps);

  // --- Sweep from high elevation (ty_max) downward ---
  for (int i = 0; i < nsteps; ++i) {
    const double tycnt = ty_max - 0.5 * ty_len - static_cast<double>(i) * ty_step;
    const double ty_lower_slice = tycnt - 0.5 * ty_len;
    const double ty_upper_slice = tycnt + 0.5 * ty_len;

    // Depth calculation: elev_from_det_lower corresponds to ty_upper_slice (higher elevation = shallower)
    // elev_from_det_upper corresponds to ty_lower_slice (lower elevation = deeper)
    const double tan_y_low = angle_util::to_tangent(ty_lower_slice, angle_unit);
    const double tan_y_up  = angle_util::to_tangent(ty_upper_slice, angle_unit);
    const double elev_from_det_lower_obj = elev_obj_from_det - tan_y_up  * delta_horizontal;
    const double elev_from_det_upper_obj = elev_obj_from_det - tan_y_low * delta_horizontal;

    // Shift to terrain-based depth: add diff_elev (= elev_terrain - elev_obj)
    const double elev_from_det_lower = elev_from_det_lower_obj + diff_elev;
    const double elev_from_det_upper = elev_from_det_upper_obj + diff_elev;

    if (elev_from_det_lower_obj < 0.0) {
      LOG_WARN_ND("sweep_y_with_prefix_sum: elev_from_det_lower_obj < 0 (elev_from_det_lower_obj={}, elev_from_det_lower={}, elev_obj_from_det={}, diff_elev={}, tan_y_up={}, delta_horizontal={}). Object protrudes above object-top reference. Skipping."
        , elev_from_det_lower_obj, elev_from_det_lower, elev_obj_from_det, diff_elev, tan_y_up, delta_horizontal);
      continue;
    }

    // const double elev_from_det_diff = elev_from_det_upper - elev_from_det_lower;

    // Convert ty slice bounds to bin indices
    const int iy_lo = g2bg.get_iy(ty_lower_slice + ty_eps);
    const int iy_hi = g2bg.get_iy(ty_upper_slice - ty_eps);

    // ★ O(1) prefix-sum query  (replaces O(nbinx × nbiny_slice) loop)
    const double signal_sum = ps_signal.query(ix_lo, ix_hi, iy_lo, iy_hi);
    const double noise_sum  = ps_noise .query(ix_lo, ix_hi, iy_lo, iy_hi);

    // Package result (identical structure to original)
    SignalNoiseDepth snd;
    snd.detid = detid;
    snd.sn.signal = signal_sum;
    snd.sn.noise  = noise_sum;
    snd.sn.txty_range.set_xmin(tx_lower + tx_eps);
    snd.sn.txty_range.set_xmax(tx_upper - tx_eps);
    snd.sn.txty_range.set_ymin(ty_lower_slice + ty_eps);
    snd.sn.txty_range.set_ymax(ty_upper_slice - ty_eps);
    snd.elev_from_det_lower = elev_from_det_lower;
    snd.elev_from_det_upper = elev_from_det_upper;
    // snd.elev_from_det_lower = elev_obj_from_det;
    // snd.elev_from_det_upper = elev_obj_from_det + elev_from_det_diff;
    vec_result.push_back(snd);
  }

  return vec_result;
}

} // anonymous namespace


// =====================================================================
//  run_mode_ddl_v1 — main body
// =====================================================================
void DepthResolutionSweep::run_mode_ddl_v1()
{
  LOG_INFO("run_mode_ddl_v1 (optimised)");
  LOG_INFO("  output path prefix: {}", params_.output_ascii_prefix.string());

  const Eigen::Vector2d v2_pos_obj = params_.obj_center_vec();
  LOG_INFO("  obj center : ({}, {})", v2_pos_obj.x(), v2_pos_obj.y());

  // --- Prepare angle-between candidates (shared with v0) ---
  set_vec_vec_angbetween(params_.obj_size_upper_limit);

  // --- Copy and initialise the base density model ---
  Grid2dPillar g2pil_base = this->g2pil_base_;
  g2pil_base.set_uniform_density(params_.base_density);

  const int n_det = arrdet_base_.get_n_det();
  for (int detid = 0; detid < n_det; ++detid)
    LOG_INFO(" detid={}, n_ang_between={}", detid,
      vec_vec_angle_between_from_obj_center_.at(detid).size());

  // ================================================================
  //  OUTER LOOP : one detector at a time
  // ================================================================
  for (int detid = 0; detid < n_det; ++detid){ // detector loop start
    auto t_det_start = std::chrono::steady_clock::now();
    LOG_INFO("[v1] ===== detid={} START =====", detid);

    // --- Get a fresh FULL panel (no cut) and point it toward the object ---
    DetectorPanel panel_full = arrdet_base_.getDetectorPanelCopy(detid);
    panel_full.set_direction_to_v2_pos(v2_pos_obj);
    const Ray3d& ray3d = panel_full.get_ray3d();

    // ────────────────────────────────────────────────────────
    //  PHASE 1: ONE-TIME heavy computation per detector
    //
    //  ★ This is the key optimisation: mp_add_PLDL is called
    //    ONCE on the full panel instead of ~160 times on
    //    progressively wider cut panels.
    // ────────────────────────────────────────────────────────
    // mp_add_PLDL called ONCE on full panel (key optimization): computes path length
    // through terrain for all bins. This replaces ~160 calls in the original algorithm.
    LOG_INFO("[v1] detid={}: Computing PLDL on full panel ({}x{}) ...",
      detid, panel_full.get_nbinx(), panel_full.get_nbiny());
    pathcalc::g2pil::mp_add_PLDL(panel_full, g2pil_base);

    // Baseline flux / signal from the full-panel DL
    panel_full.mp_calc_set_peneflux_signal_from_DL(ft_, false);

    // Build prefix-sum tables from baseline signal/noise values.
    // These tables enable O(1) rectangular range sum queries during y-sweep.
    PrefixSum2D ps_signal_base, ps_noise_base;
    ps_signal_base.build_signal(panel_full);
    ps_noise_base .build_noise (panel_full);
    LOG_INFO("[v1] detid={}: Baseline prefix sums built.", detid);

    // Cache detector geometry for sweep helper
    const Detid det_id_tag  = panel_full.get_detid();
    const Eigen::Vector3d v3_pos_det = ray3d.pos();
    const double elev_det      = v3_pos_det.z();
    const angle_util::AngleUnit angle_unit = panel_full.get_angle_unit();
    const Grid2dBinGroup& g2bg_full = panel_full.get_g2bg();

    const auto& vec_ang_between =
        vec_vec_angle_between_from_obj_center_.at(detid);

    // ────────────────────────────────────────────────────────
    //  PHASE 2: Parallel sweep over angular sizes
    //
    //  ★ Key change from previous v1:
    //    - OpenMP is applied at the ang_between level (outer loop)
    //    - All inner operations use single-thread (st_*) helpers
    //      to avoid nested OpenMP overhead
    //    - delta_dens loop uses CUT panels (not full panels)
    //    - y_sweep uses PrefixSum2D for O(1)/step queries
    //
    //  Each ang_between iteration is fully independent:
    //    no shared writes → ideal for omp parallel for.
    // ────────────────────────────────────────────────────────
    const int n_ang = static_cast<int>(vec_ang_between.size());

    // Thread-local result storage: one vector per ang_between index.
    // Each thread writes only to its own per_ang_results[i_ang].
    std::vector<std::vector<SignalNoiseStatResult>> per_ang_results(n_ang);

    #pragma omp parallel for schedule(dynamic)
    for (int i_ang = 0; i_ang < n_ang; ++i_ang){ // omp parallel for ang_between start
      const auto& ang_between = vec_ang_between.at(i_ang);

      // --- Cylinder geometry (identical to v0) ---
      const double r_cyl =
          geom_util::calc_cylinder_radius_from_horizontal_angle_size(
            ray3d, v2_pos_obj, ang_between);
      const double obj_size = r_cyl * 2.0;

      // Skip too-small objects
      if (obj_size < params_.obj_size_lower_limit) continue;

      // --- Convert ang_between to detector's native angle unit ---
      double ang_between_unit_value = 0.0;
      if (angle_unit == angle_util::AngleUnit::Tangent)
        ang_between_unit_value = std::tan(std::fabs(ang_between.rad()));
      else if (angle_unit == angle_util::AngleUnit::Radian)
        ang_between_unit_value = ang_between.rad();
      else if (angle_unit == angle_util::AngleUnit::Degree)
        ang_between_unit_value = ang_between.deg();
      else
        THROW_ERROR("DepthResolutionSweep::run_mode_ddl_v1: Unknown angle unit.");

      // Get zmax of object
      const double elev_obj = g2pil_base.getMinimumZmaxCircle( 
        v2_pos_obj.x(), v2_pos_obj.y(), r_cyl );

      // Get the top of terrain
      const double elev_terrain = g2pil_base.getPillar(
        v2_pos_obj.x(), v2_pos_obj.y() ).get_zmax();

      // get the cylinder height
      const double elev_obj_from_det = elev_obj - elev_det;
      LOG_DEBUG(" zmax of g2pil_base at obj_center: {:.1f} meters", elev_obj);
      LOG_DEBUG(" height of cyl: {:.1f} meters", elev_obj_from_det);

      // get the depth from the terrain surface to the bottom of detector
      const double elev_terrain_from_det = elev_terrain - elev_det;
      LOG_DEBUG(" elev_terrain_from_det (terrain surface to bottom of detector): {:.1f} meters", elev_terrain_from_det);
      if (obj_size > elev_terrain_from_det) { continue; }

      // get the diff between terrain surface and obj zmax
      const double diff_elev = elev_terrain - elev_obj;
      LOG_DEBUG(" diff_elev (elev_terrain - elev_obj): {:.1f} meters", diff_elev);

      const Eigen::Vector3d v3_pos_obj_top{ v2_pos_obj.x(), v2_pos_obj.y(), elev_obj };

      // --- Elevation range ---
      // tymax may be clamped to avoid negative DL (occurs when obj_size exceeds
      // path length at high elevations).
      const double tymin = g2bg_full.get_ymin();
      const double tymax_panel = g2bg_full.get_ymax();
      const double tymax_allowed = calc_ty_max_allowed(
          elev_obj_from_det, obj_size, angle_unit);
      const double tymax = std::min(tymax_panel, tymax_allowed);

      if (tymax <= tymin) continue;  // obj_size > PL at all elevations

      // --- Sweep window (same formulas as v0) ---
      const double cut_factor = params_.angle_between_cut_factor;
      const double tx_cut_min = -ang_between_unit_value * cut_factor;
      const double tx_cut_max =  ang_between_unit_value * cut_factor;

      const double sf = params_.sweep_range_factor;
      const double tx_lower = -ang_between_unit_value * sf;
      const double tx_upper =  ang_between_unit_value * sf;
      const double ty_lower = -ang_between_unit_value * sf;
      const double ty_upper =  ang_between_unit_value * sf;
      const double ty_step  = params_.elev_center_step;

      // ★ Baseline y-sweep: O(1) per step via prefix sum on FULL panel.
      //   ps_signal_base / ps_noise_base are read-only → thread-safe.
      std::vector<SignalNoiseDepth> vec_snd_base =
          sweep_y_with_prefix_sum(
            ps_signal_base, ps_noise_base, g2bg_full,
            det_id_tag,
            tx_lower, tx_upper, ty_lower, ty_upper,
            tymin, tymax, ty_step,
            v3_pos_obj_top, v3_pos_det, angle_unit, diff_elev);

      // ──────────────────────────────────────────────────
      //  Cut panel for delta_dens loop.
      //
      //  ★ Key fix from failed v1: use CUT panel (small)
      //    instead of FULL panel (320×320) for the inner loop.
      //    panel_full is read-only here → cut() is thread-safe.
      //
      //  The cut panel inherits Phase 1's DL, omega, area,
      //  exposure_time from panel_full via deep copy.
      // ──────────────────────────────────────────────────
      DetectorPanel panel_cut_base = panel_full.cut(
          tx_cut_min, tx_cut_max, tymin, tymax);
      const Grid2dBinGroup& g2bg_cut = panel_cut_base.get_g2bg();

      // ──────────────────────────────────────────────────
      //  Inner loop: delta_dens variations (on CUT panel)
      //
      //  All operations use single-thread (st_*) helpers
      //  since we are inside an omp parallel for region.
      // ──────────────────────────────────────────────────
      for (const double delta_dens : params_.vec_delta_density) {
        // DDL approximation: delta_DL = delta_dens * obj_size (uniform across all bins)
        // that's the truth of the DDL approximation.
        const double delta_DL_obj = delta_dens * obj_size;

        // Copy CUT panel (small — much cheaper than copying FULL)
        DetectorPanel panel_cut_mod(panel_cut_base);

        // Add uniform DL shift (DDL approximation) — single-thread
        st_add_DL(panel_cut_mod, delta_DL_obj);

        // Clamp negative DL to zero — single-thread
        st_clamp_DL_nonnegative(panel_cut_mod);

        // Recompute peneflux + signal only (omega/area/time unchanged)
        // — single-thread, skips calc_solid_angle / calc_effective_area
        st_recalc_peneflux_signal_from_DL(panel_cut_mod, ft_);

        // Build modified prefix sums on CUT panel
        PrefixSum2D ps_signal_mod, ps_noise_mod;
        ps_signal_mod.build_signal(panel_cut_mod);
        ps_noise_mod.build_noise(panel_cut_mod);

        // ★ Modified y-sweep: O(1) per step via prefix sum on CUT panel.
        //   Uses g2bg_cut (cut panel's local coordinate system).
        std::vector<SignalNoiseDepth> vec_snd_modi =
            sweep_y_with_prefix_sum(
              ps_signal_mod, ps_noise_mod, g2bg_cut,
              det_id_tag,
              tx_lower, tx_upper, ty_lower, ty_upper,
              tymin, tymax, ty_step,
              v3_pos_obj_top, v3_pos_det, angle_unit, diff_elev);

        // ──────────────────────────────────────────
        //  Statistical significance tests
        //  (identical to v0 — lightweight)
        // ──────────────────────────────────────────
        for (const auto& [signal_amp, noise_amp] : params_.vec_signal_noise_amplifiers){
          for (const double stat_alpha : params_.vec_stat_alpha_value){
            signal_noise_stat_util::eval_signal_significance(
                vec_snd_base, vec_snd_modi,
                obj_size, delta_dens, stat_alpha, signal_amp,
                per_ang_results.at(i_ang),
                params_.both_side); // true=two-sided, false(default)=one-sided
          }
        }
      } // for delta_dens
    } // for i_ang (omp parallel for)

    // ──────────────────────────────────────────────────
    //  Merge thread-local results (sequential)
    // ──────────────────────────────────────────────────
    std::vector<SignalNoiseStatResult> vec_signal_stat_result;
    for (auto& v : per_ang_results) {
      vec_signal_stat_result.insert(
          vec_signal_stat_result.end(),
          std::make_move_iterator(v.begin()),
          std::make_move_iterator(v.end()));
    }

    // ────────────────────────────────────────────────────────
    //  Output: sort and write results (identical to v0)
    // ────────────────────────────────────────────────────────
    std::vector<SignalNoiseStatSortCondition> vec_sort_condition = {
        {SignalNoiseStatSortKey::StatAlpha,   false}
      , {SignalNoiseStatSortKey::SignalAmp,   true}
      , {SignalNoiseStatSortKey::DeltaDens,   false}
      , {SignalNoiseStatSortKey::ObjSize,     true}
      , {SignalNoiseStatSortKey::Depth_lower, true}
    };
    signal_noise_stat_util::sort_signal_stat_results(
        vec_signal_stat_result, vec_sort_condition);

    // ASCII output (.tmp)
    char cfnameout_stat[512];
    sprintf(cfnameout_stat, "%s_det%02d_signal_signifi.tmp",
      params_.output_ascii_prefix.string().c_str(), detid);
    fs::path output_path_stat = iodir::make_pathout(cfnameout_stat);
    LOG_INFO("[v1] Outputting results to {}", output_path_stat.string());
    signal_noise_stat_util::out_signal_stat_result_vector(
        vec_signal_stat_result, output_path_stat);

    // CSV output
    char cfnameout_csv[512];
    sprintf(cfnameout_csv, "%s_det%02d_signal_signifi-tmp.csv",
      params_.output_ascii_prefix.string().c_str(), detid);
    fs::path output_path_csv = iodir::make_pathout(cfnameout_csv);
    LOG_INFO("[v1] Outputting results to {}", output_path_csv.string());
    signal_noise_stat_util::out_signal_stat_result_vector_csv(
        vec_signal_stat_result, output_path_csv);

    auto t_det_end = std::chrono::steady_clock::now();
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        t_det_end - t_det_start).count();
    LOG_INFO("[v1] ===== detid={} DONE ({} ms) =====", detid, elapsed_ms);

  } // for detid

  LOG_INFO("[v1] run_mode_ddl_v1 complete.");
}


// Dispatch to run_mode_ddl() or run_mode_cylinder() depending on mode
void DepthResolutionSweep::run()
{
  LOG_INFO("dispatch by mode.");
  switch (params_.mode) {
    case Parameters::Mode::Cylinder:
      LOG_INFO("  mode: cylinder_v3");
      run_mode_cylinder();
      break;
    case Parameters::Mode::DDL:
      LOG_INFO("  mode: DDL (v1 optimised)");
      run_mode_ddl_v1();
      break;
    default:
      THROW_ERROR("Unknown mode in DepthResolutionSweep::run()");
  }
}
