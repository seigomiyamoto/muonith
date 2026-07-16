// cls_DetectorPanel.cpp
#include "cls_DetectorPanel.hpp"
#include "cls_DetectorPanelParameters.hpp"
#include "cls_EfficiencyModel.hpp"
#include "ns_myapp.hpp"
#include "ns_mymacro.hpp"
#include "ns_mylogger.hpp"
#include "cls_RectAngularBinGroup.hpp"
#include <memory>
#include "spdlog_pch.hpp"
#include "ns_mylogger.hpp"
#include "ns_stats_util.hpp"
#include "ns_geom_util.hpp"
#include <cmath>
#include <limits>
#include "ns_iodir.hpp"

//###################################################################################
//###################################################################################
// class DetectorPanel
//###################################################################################
//###################################################################################

// Inequality operator
// Name is not compared.
bool DetectorPanel::operator!=(const DetectorPanel& other) const
{
  #ifdef NODEBUG
    if (g2bg_ !=(other.g2bg_)) return true;
    if (detid_ != other.detid_) return true;
    if (ray3d != other.ray3d) return true;
    if (v3_det_length != other.v3_det_length) return true;
    if (n_unit != other.n_unit) return true;
    if (days != other.days) return true;
    if (n_element != other.n_element) return true;
    if (vec_vec_DetectorElement != other.vec_vec_DetectorElement) return true;
    if (rotation_type != other.rotation_type) return true;
    if (prm_bingrp != other.prm_bingrp) return true;
    if (angle_unit != other.angle_unit) return true;
  #else
    if (g2bg_ != other.g2bg_) { LOG_WARN("DetectorPanel: g2bg_ differs"); return true; }
    if (detid_ != other.detid_) { LOG_WARN("DetectorPanel: detid_ differs"); return true; }
    if (ray3d != other.ray3d) { LOG_WARN("DetectorPanel: ray3d differs"); return true; }
    if (v3_det_length != other.v3_det_length) { LOG_WARN("DetectorPanel: v3_det_length differs"); return true; }
    if (n_unit != other.n_unit) { LOG_WARN("DetectorPanel: n_unit differs"); return true; }
    if (days != other.days) { LOG_WARN("DetectorPanel: days differs"); return true; }
    if (n_element != other.n_element) { LOG_WARN("DetectorPanel: n_element differs"); return true; }
    if (vec_vec_DetectorElement != other.vec_vec_DetectorElement) { LOG_WARN("DetectorPanel: vec_vec_DetectorElement differs"); return true; }
    if (rotation_type != other.rotation_type) { LOG_WARN("DetectorPanel: rotation_type differs"); return true; }
    if (prm_bingrp != other.prm_bingrp) { LOG_WARN("DetectorPanel: prm_bingrp differs"); return true; }
    if (angle_unit != other.angle_unit) { LOG_WARN("DetectorPanel: angle_unit differs"); return true; }
  #endif
  return false;
};

// Move assignment operator
DetectorPanel& DetectorPanel::operator=(DetectorPanel&& other) noexcept {
  if (this != &other) {
    name = std::move(other.name);
    g2bg_ = std::move(other.g2bg_);
    pdic_ = other.pdic_;
    detid_ = other.detid_;
    ray3d = std::move(other.ray3d);
    v3_det_length = other.v3_det_length;
    n_unit = other.n_unit;
    days = other.days;
    n_element = other.n_element;
    vec_vec_DetectorElement = std::move(other.vec_vec_DetectorElement);
    rotation_type = other.rotation_type;
    prm_bingrp = std::move(other.prm_bingrp);
    angle_unit = other.angle_unit;
    other.pdic_ = nullptr;
  }
  return *this;
}

// Move constructor
DetectorPanel::DetectorPanel(DetectorPanel&& other) noexcept
  : name(std::move(other.name))
  , g2bg_(std::move(other.g2bg_))
  , pdic_(other.pdic_)
  , detid_(other.detid_)
  , ray3d(std::move(other.ray3d))
  , v3_det_length(other.v3_det_length)
  , n_unit(other.n_unit)
  , days(other.days)
  , n_element(other.n_element)
  , vec_vec_DetectorElement(std::move(other.vec_vec_DetectorElement))
  , rotation_type(other.rotation_type)
  , prm_bingrp(std::move(other.prm_bingrp))
  , angle_unit(other.angle_unit)
{
  other.pdic_ = nullptr;
}

// Detector panel constructor (OpenMP version)
DetectorPanel::DetectorPanel(
    const DetectorPanel::Parameters &detprm
  , const Detid detid_in
  , const Uqid uqid_min_tmp, const Uqid uqid_max_tmp
  , std::vector<UqidInfo>& vec_UqidInfo
  , DetectorIndexContainer &dic ) : DetectorPanel()
{
  
  // assignment of detector name
  name = detprm.name;

  // exposure time
  days = detprm.days;

  // shared angle unit
  angle_unit = detprm.angle_unit;

  // read bin list

  LOG_DEBUG("set detid_ = {}", detid_in);
  set_detid(detid_in);

  g2bg_.set_detid(detid_in);

  LOG_DEBUG("set pointer of DetectorIndexContainer dic to DetectorPanel and Grid2dBinGroup");
  this->set_dic_ptr(&dic);
  g2bg_.set_dic_ptr(&dic);
  
  const double tx_interval = ( detprm.txmax - detprm.txmin )/((double) detprm.nbinx );
  const double ty_interval = ( detprm.tymax - detprm.tymin )/((double) detprm.nbiny );
  
  // assignment of x and y axis
  const std::string name_x_axis = "x_axis_" + name;
  const std::string name_y_axis = "y_axis_" + name;
  
  // setting parameter of Grid1d x_axis and y_axis
  const Grid1d x_axis_tmp(name_x_axis, detprm.nbinx, detprm.txmin, detprm.txmax, tx_interval);
  const Grid1d y_axis_tmp(name_y_axis, detprm.nbiny, detprm.tymin, detprm.tymax, ty_interval);
  
  set_x_axis(x_axis_tmp);
  get_x_axis().out_info(spdlog::level::debug);

  set_y_axis(y_axis_tmp);
  get_y_axis().out_info(spdlog::level::debug);

  // detector rotation
  const Angle  roll( detprm.roll );
  const Angle pitch( detprm.pitch);
  const Angle   yaw( detprm.yaw  );

  // rotation matrix should be defined by rotation_type
  // yaw is clockwise, so change the sign
  const Eigen::Matrix3d rotationMatrix =
    angle_util::make_rotation_3d_matrix_ZYX(
      roll, pitch, -yaw, detprm.rotation_type);

PRINTF(detprm.roll.deg());
roll.out();
// PRINT_MAT(rollAngle.matrix());

PRINTF(detprm.pitch.deg());
pitch.out();
// PRINT_MAT(pitchAngle.matrix());

PRINTF(detprm.yaw.deg());
yaw.out();
// PRINT_MAT(yawAngle.matrix());

PRINT_MAT(rotationMatrix);

  // initialzation of v3_direction, yaw=0 at north
  const Eigen::Vector3d v3_dir0(0,1,0);

  LOG_INFO("det_id={:02d}, v3_dir0 = ({:7.4f}, {:7.4f}, {:7.4f})"
    ,detid_in,v3_dir0(0),v3_dir0(1),v3_dir0(2));

// PRINT_MAT(quat.matrix());
  // rotation of v3_direction
  // const Eigen::Vector3d v3_dir1 = quat * v3_dir0;
  const Eigen::Vector3d v3_dir1 = rotationMatrix * v3_dir0;
  // v3_direction = quat.matrix() * v3_direction;
  // PRINT_MAT(v3_dir1);
  LOG_INFO("det_id={:02d}, v3_dir1 = ({:7.4f}, {:7.4f}, {:7.4f})"
    ,detid_in,v3_dir1(0),v3_dir1(1),v3_dir1(2));
  this->ray3d.set_dir(v3_dir1);

  // allocate memory for std::vector<std::vector<DetectorElement>>
  vec_vec_memory_allocate();

  // reserve memory for DetectorElement::vec_tf_in_PL
  reserve_vec_tf_in_PL( detprm.n_reserve_vec_tf_in_PL );

  // later this is used for detector element loop again.
  const Eigen::Vector3d v3_pos_ele(detprm.v3_position);

  // detector position assignment
  // v3_position = v3_pos_ele;
  set_position(v3_pos_ele);

  // assignment of detector size
  set_v3_det_length(detprm.v3_det_length);

  LOG_INFO("det_id={:02d}, v3_pos = ({:7.4f}, {:7.4f}, {:7.4f})"
    ,detid_in,v3_pos_ele(0),v3_pos_ele(1),v3_pos_ele(2));

  // assignment of n_unit
  n_unit = detprm.n_unit;

  // day to sec factor
  constexpr double factor_day2sec = 24.*3600.;

  // number of DetectorElements made in this function
  int n_element_tmp = 0;
  
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  // const int n_element_predicted = nbinx * nbiny;

  const double x_interval = get_x_interval();
  const double y_interval = get_y_interval();

  const int n_element = nbinx * nbiny;

  // memory allocation for vec_UqidInfo
  vec_UqidInfo.resize(n_element);

  // y-axis loop
  #pragma omp parallel for collapse(2)
  for(int iy=0;iy<nbiny;iy++){
    // x-axis loop
    for(int ix=0;ix<nbinx;ix++){
      // calculation of tymin, tymax
      double tymin = detprm.tymin + ((double)iy)*y_interval;
      double tymax = tymin + y_interval;
      double tycnt = 0.5*(tymin + tymax);
      // At this point, uqid and id_in_this_detector are determined
      // as functions of uqid_min_tmp, ix, iy, and nbinx.
      Inthis id_in_this_detector = 0 + ix + iy * nbinx; 
      Uqid uqid = uqid_min_tmp + id_in_this_detector;

      DetectorElement& ele = callDetectorElement(ix,iy);

      // initilization of DetectorElement
      ele.reset();
      
      // setting of id_detector and id_in_this_detector
      ele.set_detid(detid_in);
      ele.set_id_in_this_detector(id_in_this_detector);

      // assign detector_index 
      ele.set_uqid(uqid);
      
      // calculation of txmin, txmax
      double txmin = detprm.txmin + ((double)ix)*x_interval;
      double txmax = txmin + x_interval;
      double txcnt = 0.5*(txmin + txmax);
      
      // assignment of txmin, .... , tymax 
      // callDetectorElement(ix,iy).set_txmin_txmax_tymin_tymax(txmin,txmax,tymin,tymax);
      ele.set_txmin_txmax_tymin_tymax(txmin,txmax,tymin,tymax,detprm.angle_unit);
      
      // define v3_dir_ele for later use
      Eigen::Vector3d v3_dir_ele(0,1,0);

      // calc v3_direction of DetectorElement
      if ( detprm.angle_unit == DetectorElement::AngleUnit::Tangent ){
        v3_dir_ele = angle_util::tx_ty_vertical_to_v3_tangent(txcnt, tycnt);
      }
      else if ( detprm.angle_unit == DetectorElement::AngleUnit::Degree ){
        const Angle txcnt_angle = Angle(txcnt, Angle::Unit::Degree);
        const Angle tycnt_angle = Angle(tycnt, Angle::Unit::Degree);
        v3_dir_ele = angle_util::tx_ty_vertical_to_v3(txcnt_angle, tycnt_angle);
      }
      else if ( detprm.angle_unit == DetectorElement::AngleUnit::Radian ){
        const Angle txcnt_angle = Angle(txcnt, Angle::Unit::Radian);
        const Angle tycnt_angle = Angle(tycnt, Angle::Unit::Radian);
        v3_dir_ele = angle_util::tx_ty_vertical_to_v3(txcnt_angle, tycnt_angle);
      }
      else{
        THROW_ERROR("angle_unit is not correct");
      }

      // rotation of unit vector by detector
      // v3_dir_ele = quat * v3_dir_ele;
      Eigen::Vector3d v3_dir_ele_rot = rotationMatrix * v3_dir_ele;
      ele.set_v3_direction(v3_dir_ele_rot);

      // position assignment of element
      ele.set_v3_position(v3_pos_ele);

      // assignment of exposure time
      ele.set_exposure_time_sec( days*factor_day2sec );

      // assignment of solid angle
      const double omega = ele.calc_solid_angle();
      ele.set_solid_angle(omega);

      // assignment of effective_area
      const double effective_area = ele.calc_effective_area(ray3d.dir(),v3_det_length,get_n_unit());
      ele.set_effective_area_m2(effective_area);

      //
      // Add uqid_info to vec_UqidInfo
      //
      UqidInfo uinfo;
      uinfo.uqid = uqid;
      uinfo.detid = detid_in;
      uinfo.ixiy = Ixiy{ix,iy};
      uinfo.inthis = id_in_this_detector;
      uinfo.uqig = UqigNotAssigned; // Not yet assigned to any group
      uinfo.uqig_avail = UqigAvailNotAssigned; // Not yet assigned to any available group
      uinfo.is_avail = false;  // At this stage, availability is not yet determined
      vec_UqidInfo.at(id_in_this_detector) = uinfo;

    } // xloop end
  } // yloop end
  set_n_element(n_element);

  // Compare size of vec_UqidInfo with n_element
  assert(vec_UqidInfo.size() == n_element);

  const Uqid uqid_max_of_vec = vec_UqidInfo.back().uqid;

  // uqid consistency check. if not consistent, throw error
  if( uqid_max_of_vec != uqid_max_tmp ){
    LOG_ERROR("uqid_max_of_vec = {}, uqid_max_tmp = {}", uqid_max_of_vec, uqid_max_tmp);
    THROW_ERROR_NAME("uqid_max_of_vec != uqid_max_tmp");
  }

  // when detprm.filepath_bin_list is not default value or null,
  // read Grid2dBinGroup data file
  if( detprm.tf_read_bin_list==false ){
    hoge("g2bg_.read_bin_list2 was not called");
  }else{
    hoge("g2bg_.read_bin_list2 was called");
    g2bg_.read_bin_list2(detprm.filepath_bin_list);
    const int n_merged_list = pdic_->get_n_group(detid_in);
    if( n_merged_list==0 ){
      THROW_ERROR_NAME("g2bg_.read_bin_list2 n_merged_list==0");
    }
    hoge(n_merged_list);
  }
}

const DetectorElement& DetectorPanel::getDetectorElement( const int ix, const int iy ) const
{
  //macro, if out of range, throw error
  check_ix_inside(ix);
  check_iy_inside(iy);
  return vec_vec_DetectorElement.at(iy).at(ix);
}

const DetectorElement& DetectorPanel::getDetectorElement( const double tx, const double ty ) const
{
  const int ix = get_ix(tx);
  const int iy = get_iy(ty);
  //macro, if out of range, throw error
  check_ix_inside(ix);
  check_iy_inside(iy);
  return getDetectorElement(ix,iy);
}

// call DetectorElement
DetectorElement& DetectorPanel::callDetectorElement(const int ix, const int iy)
{
  //macro, if out of range, throw error
  check_ix_inside(ix);
  check_iy_inside(iy);
  return vec_vec_DetectorElement.at(iy).at(ix);
}

// call DetectorElement
DetectorElement& DetectorPanel::callDetectorElement(const double tx, const double ty)
{
  const int ix = get_ix(tx);
  const int iy = get_iy(ty);
  //macro, if out of range, throw error
  check_ix_inside(ix);
  check_iy_inside(iy);
  return callDetectorElement(ix,iy);
}

// allocate memory for std::vector<std::vector<DetectorElement>> vec_vec_DetectorElement;
void DetectorPanel::vec_vec_memory_allocate()
{
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  vec_vec_DetectorElement.resize(nbiny); // allocate y direction
  for(int iy=0;iy<get_nbiny();iy++){
    vec_vec_DetectorElement.at(iy).resize(nbinx); // allocate x direction
  }
  g2bg_.init_vec_vec();
}

// reserve memory for DetectorElement::vec_tf_in_PL
void DetectorPanel::reserve_vec_tf_in_PL(const int n_reserve_vec_tf_in_length)
{
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  #pragma omp parallel for collapse(2) schedule(static)
  for(int iy=0;iy<get_nbiny();iy++){
    for(int ix=0;ix<get_nbinx();ix++){
      DetectorElement& ele = callDetectorElement(ix,iy);
      ele.reserve_vec_tf_in_PL(n_reserve_vec_tf_in_length);
    }
  }
}

// apply calc_approx_volue to DetectorElement(ix,iy)
double DetectorPanel::calc_approx_volume( const int ix, const int iy ) const
{
  const DetectorElement& ele = getDetectorElement(ix,iy);
  return ele.calc_approx_volume();
}

// Calculate calc_approx_volume for all {(ix,iy)} belonging to igroup and return the sum
double DetectorPanel::calc_approx_volume_sum_grouped( const Igroup igroup ) const
{
  const std::vector<Ixiy> vec_ixiy = g2bg_.get_vec_ixiy(igroup);
  double volume_sum = 0.0;
  for(const auto [ix,iy] : vec_ixiy){
    volume_sum += calc_approx_volume(ix,iy);
  }
  return volume_sum;
}

// Compute nmuon_grouped for a given density.
double DetectorPanel::calc_nmuon_grouped(
  const Igroup igroup, const double dens, const FluxTable &ft ) const
{
  // First, get {(ix, iy)} belonging to igroup.
  const std::vector<Ixiy> vec_ixiy = g2bg_.get_vec_ixiy(igroup);

  // Loop over vec_ixiy.
  double nmuon_grouped = 0.0;
  for(const auto [ix,iy] : vec_ixiy){
    const DetectorElement& ele = getDetectorElement(ix,iy);
    const double costhz = ele.get_ray3d().vz();
    const double PL = ele.get_PL();
    const double DL = PL * dens;
    const double peneflux = ft.get_peneflux(costhz,DL);
    const double SOT = ele.get_SOT();
    const double nmuon = peneflux * SOT;
    nmuon_grouped += nmuon;
  }
  return nmuon_grouped;
}

// Compute nmuon_grouped for a given density (parallel).
double DetectorPanel::mp_calc_nmuon_grouped(
  const Igroup igroup, const double dens, const FluxTable &ft ) const
{
  // Get {(ix, iy)} belonging to igroup.
  const std::vector<Ixiy> vec_ixiy = g2bg_.get_vec_ixiy(igroup);

  double nmuon_grouped = 0.0;
  #pragma omp parallel for reduction(+: nmuon_grouped)
  for (size_t i = 0; i < vec_ixiy.size(); i++) {
    const auto [ix, iy] = vec_ixiy[i];
    const DetectorElement& ele = getDetectorElement(ix, iy);
    const double costhz = ele.get_ray3d().vz();
    const double PL = ele.get_PL();
    const double DL = PL * dens;
    const double peneflux = ft.get_peneflux(costhz, DL);
    const double SOT = ele.get_SOT();
    const double nmuon = peneflux * SOT;
    nmuon_grouped += nmuon;
  }
  return nmuon_grouped;
}

// Return nmuon_grouped results over a density range.
std::vector<std::array<double,2>> DetectorPanel::calc_vec_dens_nmuon_grouped(
  const Igroup igroup
, const double dens_min, const double dens_max, const double dens_step
, const FluxTable &ft ) const
{
  const int num_dens = (int) ((dens_max-dens_min)/dens_step) + 1;

  // Storage for 1D grid-search results over density.
  std::vector<std::array<double,2>> vec_dens_nmuon_grouped(num_dens);
  
  // Run 1D grid search over density and store results in vec_dens_nmuon_grouped.
  #pragma omp parallel for
  for(int idens=0;idens<num_dens;idens++){
    const double dens = dens_min + ((double)idens)*dens_step;
    const double nmuon_grouped = calc_nmuon_grouped(igroup,dens,ft);
    vec_dens_nmuon_grouped.at(idens) = {dens,nmuon_grouped};
  }
  return vec_dens_nmuon_grouped;
}

// Perform adaptive grid search to minimize |calc_nmuon_grouped - target_value|.
std::array<double, 2> DetectorPanel::calc_dens_nmuon_grouped(
  const Igroup igroup, 
  const double dens_min, const double dens_max, 
  const std::vector<double> &vec_dens_step, const double range_factor,
  const FluxTable &ft, 
  const double target_value)
{
  std::string msg;

  // Sort vec_dens_step in descending order (copy because input is const).
  std::vector<double> sorted_dens_step = vec_dens_step;
  std::sort(sorted_dens_step.begin(), sorted_dens_step.end(), std::greater<double>());

  double best_dens = dens_min;
  double best_error = std::numeric_limits<double>::max();
  double current_min = dens_min;
  double current_max = dens_max;

  // Search for each step size in sorted_dens_step.
  for (size_t i = 0; i < sorted_dens_step.size(); ++i) {
    double current_step = sorted_dens_step[i];
    int num_steps = static_cast<int>((current_max - current_min) / current_step) + 1;
    for (int j = 0; j < num_steps; ++j) {
      double dens = current_min + j * current_step;
      double nmuon_grouped = mp_calc_nmuon_grouped(igroup, dens, ft);
      double error = fabs(nmuon_grouped - target_value);
      if (error < best_error) {
        best_error = error;
        best_dens = dens;
      }
    }
    
    // Update the search range if a next step size exists.
    if (i + 1 < sorted_dens_step.size()) {
      // Center the new search range around current best_dens.
      double new_range_min = best_dens - range_factor * sorted_dens_step[i];
      double new_range_max = best_dens + range_factor * sorted_dens_step[i];

      // Check if the new lower bound is smaller than current_min.
      if (new_range_min < current_min) {
        if (current_min == dens_min) {
          new_range_min = dens_min;
        } else {
          LOG_WARN("new current_min ({}) is smaller than the old current_min ({})"
            , new_range_min, current_min);
        }
      }
      // Check if the new upper bound is larger than current_max.
      if (new_range_max > current_max) {
        if (current_max == dens_max) {
          new_range_max = dens_max;
        } else {
          LOG_WARN("new current_max ({}) is greater than the old current_max ({})"
            , new_range_max, current_max);
        }
      }
      
      // Update the search range (respect dens_min/dens_max).
      current_min = std::max(dens_min, new_range_min);
      current_max = std::min(dens_max, new_range_max);
    }
  }

  return {best_dens, best_error};
}

// Perform adaptive grid search to minimize |calc_nmuon_grouped - target_value|.
std::array<double, 2> DetectorPanel::calc_dens_nmuon_grouped_adaptive_search(
  const Igroup igroup, 
  const double dens_min, const double dens_max, 
  const double initial_step, const double min_step,
  const double refinement_factor, 
  const int num_refinements,
  const FluxTable &ft, 
  const double target_value)
{
  double best_dens = dens_min;
  double best_error = std::numeric_limits<double>::max();
  double current_step = initial_step;
  double current_min = dens_min;
  double current_max = dens_max;

  for (int iter = 0; iter < num_refinements; ++iter) {
    int num_steps = static_cast<int>((current_max - current_min) / current_step) + 1;
    for (int i = 0; i < num_steps; ++i) {
      double dens = current_min + i * current_step;
      double nmuon_grouped = mp_calc_nmuon_grouped(igroup, dens, ft);
      double error = fabs(nmuon_grouped - target_value);
      if (error < best_error) {
        best_error = error;
        best_dens = dens;
      }
    }
    
    // Compute the refined step size.
    double new_step = current_step / refinement_factor;
    // If new_step is below min_step, use min_step and finish.
    if (new_step < min_step) {
      current_step = min_step;
      int num_steps = static_cast<int>((current_max - current_min) / current_step) + 1;
      for (int i = 0; i < num_steps; ++i) {
        double dens = current_min + i * current_step;
        double nmuon_grouped = mp_calc_nmuon_grouped(igroup, dens, ft);
        double error = fabs(nmuon_grouped - target_value);
        if (error < best_error) {
          best_error = error;
          best_dens = dens;
        }
      }
      break;
    } else {
      current_step = new_step;
      // Narrow the search range around the best point.
      current_min = std::max(dens_min, best_dens - current_step);
      current_max = std::min(dens_max, best_dens + current_step);
    }
  }

  return {best_dens, best_error};
}

// Efficiency-uncertainty variance of the grouped count (count^2), panel-local.
// Mirrors DetectorPanelArray::get_vecxf_var_eff_all for a single group.
double DetectorPanel::calc_var_eff_group(
  const Igroup igroup, const bool tf_independent) const
{
  // Accumulate sigma_eff_i * b_i over the merge-group elements.
  double sum_contrib = 0.0; // for common mode (square of sum)
  double sum_sq      = 0.0; // for independent mode (sum of squares)
  const auto vec_ixiy = g2bg_.get_vec_ixiy(igroup);
  for (const auto& [ix, iy] : vec_ixiy) {
    const DetectorElement& ele = getDetectorElement(ix, iy);
    // efficiency sigma from the (lower, upper) band half-width
    const double sigma_eff = 0.5 * (ele.get_eff_upp() - ele.get_eff_low());
    // base count b_i: efficiency-free signal (get_signal() may be eff-applied)
    const double base_count = ele.calc_signal();
    const double contrib = sigma_eff * base_count;
    sum_contrib += contrib;
    sum_sq      += contrib * contrib;
  }
  return tf_independent ? sum_sq : sum_contrib * sum_contrib;
}

// Compute grouped-bin densities and set vec_dens_group_center/upper/lower.
void DetectorPanel::calc_set_proj_dens_grouped(
  const bool tf_signal_poisson
, const Igroup igroup, const double dens_min, const double dens_max
, const double dens_step, const double sigma
, const FluxTable &ft
, const bool tf_eff, const bool tf_eff_independent )
{
  // get det_id
  const int det_id = get_detid();

  // get is_avail
  const bool tf_avail = g2bg_.is_avail_group(igroup);
  if( ! tf_avail ){
    // if not available, set not_found_dens to vec_dens_group_center, vec_dens_group_upper, vec_dens_group_lower
    g2bg_.set_dens_group_lower(igroup,not_found_dens);
    g2bg_.set_dens_group_center(igroup,not_found_dens);
    g2bg_.set_dens_group_upper(igroup,not_found_dens);
    g2bg_.set_delta_nmuon_group_lower(igroup,not_found_delta_nmuon);
    g2bg_.set_delta_nmuon_group_center(igroup,not_found_delta_nmuon);
    g2bg_.set_delta_nmuon_group_upper(igroup,not_found_delta_nmuon);
    return;
  }

  // Get nmuon_group_center.
  // signal: Poisson switch applies. noise: always det floor + Poisson(poi bucket).
  const double sig_term_center
   = tf_signal_poisson ? get_signal_group_poisson(igroup) : get_signal_group(igroup);
  const double noi_term_center
   = get_noise_det_group(igroup) + get_noise_poi_group_poisson(igroup);
  double nmuon_group_center = sig_term_center + noi_term_center;

  if(nmuon_group_center<=0.0){
    // if nmuon_group_center is less than 0, set not_found_dens to vec_dens_group_center, vec_dens_group_upper, vec_dens_group_lower
    g2bg_.set_dens_group_lower(igroup,not_found_dens);
    g2bg_.set_dens_group_center(igroup,not_found_dens);
    g2bg_.set_dens_group_upper(igroup,not_found_dens);
    g2bg_.set_delta_nmuon_group_lower(igroup,not_found_delta_nmuon);
    g2bg_.set_delta_nmuon_group_center(igroup,not_found_delta_nmuon);
    g2bg_.set_delta_nmuon_group_upper(igroup,not_found_delta_nmuon);
    return;
  }

  // Define nmuon_group_upper and nmuon_group_lower.
  // be careful, upper/lower means density error
  // Add efficiency variance in quadrature to the Poisson term when tf_eff is on:
  // sigma_total = sqrt(N_g + var_eff_g), matching the C_N diagonal form (id-6bgzem).
  const double var_eff
    = tf_eff ? calc_var_eff_group(igroup, tf_eff_independent) : 0.0;
  const double nmuon_group_sigma = sqrt(nmuon_group_center + var_eff);
  double nmuon_group_upper
    = nmuon_group_center - sigma * nmuon_group_sigma;
  if (nmuon_group_upper < 0.0) nmuon_group_upper = 0.0;

  double nmuon_group_lower = nmuon_group_center + sigma * nmuon_group_sigma;

  // Find densities closest to signal_group and to nmuon_group_upper/lower.
  const int num_dens = (int) ((dens_max-dens_min)/dens_step) + 1;

  // Initialize delta_nmuon_min.
  double center_delta_nmuon_min = std::numeric_limits<double>::max();
  double upper_delta_nmuon_min = std::numeric_limits<double>::max();
  double lower_delta_nmuon_min = std::numeric_limits<double>::max();

  // Initialize density values.
  double dens_center = not_found_dens;
  double dens_upper = not_found_dens;
  double dens_lower = not_found_dens;

  // 1D grid-search results over density.
  std::vector<std::array<double,2>> vec_dens_nmuon_grouped
   = calc_vec_dens_nmuon_grouped(igroup,dens_min,dens_max,dens_step,ft);

  // Compute dens_center/upper/lower.
  for(int idens=0;idens<num_dens;idens++){
    const double dens =vec_dens_nmuon_grouped.at(idens).at(0);
    const double nmuon_grouped = vec_dens_nmuon_grouped.at(idens).at(1);
    const double delta_nmuon_lower  = fabs(nmuon_grouped - nmuon_group_lower);
    const double delta_nmuon_center = fabs(nmuon_grouped - nmuon_group_center);
    const double delta_nmuon_upper  = fabs(nmuon_grouped - nmuon_group_upper);

    // for debug
    if(delta_nmuon_lower==0.0 ){
      LOG_WARN_ND(
        "delta_nmuon_lower==0.0 : det_id={:02d}, igroup={:5d}, dens={:6.0f}, nmuon_grouped={:E}, nmuon_group_lower={:E}, nmuon_group_center={:E}, nmuon_group_upper={:E}"
        ,det_id,igroup,dens,nmuon_grouped,nmuon_group_lower,nmuon_group_center,nmuon_group_upper);
    }
    if(delta_nmuon_center==0.0 ){
      LOG_WARN_ND(
        "delta_nmuon_center==0.0 : det_id={:02d}, igroup={:5d}, dens={:6.0f}, nmuon_grouped={:E}, nmuon_group_lower={:E}, nmuon_group_center={:E}, nmuon_group_upper={:E}"
        ,det_id,igroup,dens,nmuon_grouped,nmuon_group_lower,nmuon_group_center,nmuon_group_upper);
    }
    if(delta_nmuon_upper==0.0 ){
      LOG_WARN_ND(
        "delta_nmuon_upper==0.0 : det_id={:02d}, igroup={:5d}, dens={:6.0f}, nmuon_grouped={:E}, nmuon_group_lower={:E}, nmuon_group_center={:E}, nmuon_group_upper={:E}"
        ,det_id,igroup,dens,nmuon_grouped,nmuon_group_lower,nmuon_group_center,nmuon_group_upper);
    }

    // update dens_lower
    if( delta_nmuon_lower < lower_delta_nmuon_min ){
      lower_delta_nmuon_min = delta_nmuon_lower;
      dens_lower = dens;
    }
    // update dens_center
    if( delta_nmuon_center < center_delta_nmuon_min ){
      center_delta_nmuon_min = delta_nmuon_center;
      dens_center = dens;
    }
    // update dens_upper
    if( delta_nmuon_upper < upper_delta_nmuon_min ){
      upper_delta_nmuon_min = delta_nmuon_upper;
      dens_upper = dens;
    }
  }

  // Set dens_center/upper/lower into vec_dens_group_center/upper/lower.
  g2bg_.set_dens_group_lower(igroup,dens_lower);
  g2bg_.set_dens_group_center(igroup,dens_center);
  g2bg_.set_dens_group_upper(igroup,dens_upper);

  // Set delta_nmuon_* into vec_delta_nmuon_group_*.
  g2bg_.set_delta_nmuon_group_lower(igroup,lower_delta_nmuon_min);
  g2bg_.set_delta_nmuon_group_center(igroup,center_delta_nmuon_min);
  g2bg_.set_delta_nmuon_group_upper(igroup,upper_delta_nmuon_min);

  LOG_DEBUG_ND(
    "calc_set_proj_dens_grouped : igroup={:5d}, dens_lower={:6.0f}, dens_center={:6.0f}, dens_upper={:6.0f}"
    ,igroup,dens_lower,dens_center,dens_upper);
}

// Compute projected group densities via adaptive grid search and store results.
void DetectorPanel::calc_set_proj_dens_grouped(
  const bool tf_signal_poisson
, const Igroup igroup, const double sigma
, const double dens_min, const double dens_max
, const std::vector<double> &vec_dens_step, const double range_factor
, const FluxTable &ft
, const bool tf_eff, const bool tf_eff_independent )
{
  const int det_id = get_detid();

  // If the group is unavailable, set not_found values.
  if (!g2bg_.is_avail_group(igroup)) {
    g2bg_.set_dens_group_lower(igroup, not_found_dens);
    g2bg_.set_dens_group_center(igroup, not_found_dens);
    g2bg_.set_dens_group_upper(igroup, not_found_dens);
    g2bg_.set_delta_nmuon_group_lower(igroup, not_found_delta_nmuon);
    g2bg_.set_delta_nmuon_group_center(igroup, not_found_delta_nmuon);
    g2bg_.set_delta_nmuon_group_upper(igroup, not_found_delta_nmuon);
    return;
  }

  // Get nmuon_group_center.
  // signal: Poisson switch applies. noise: always det floor + Poisson(poi bucket).
  const double sig_term_center
   = tf_signal_poisson ? g2bg_.get_signal_group_poisson(igroup) : g2bg_.get_signal_group(igroup);
  const double noi_term_center
   = g2bg_.get_noise_det_group(igroup) + g2bg_.get_noise_poi_group_poisson(igroup);
  double nmuon_group_center = sig_term_center + noi_term_center;
  {
    // noerror baseline: signal mean + det floor + poi mean (no Poisson fluctuation)
    const double nmuon_group_center_noerror = g2bg_.get_signal_group(igroup)
      + g2bg_.get_noise_det_group(igroup) + g2bg_.get_noise_poi_group(igroup);
    const double diff = nmuon_group_center - nmuon_group_center_noerror;
    LOG_TRACE_ND("calc_set_proj_dens_grouped : det_id={:02d}, igroup={:4d}, nmuon_poisson={:5.0f}, nmuon_noerror={:5.0f}, diff={:5.0f}"
      ,det_id,igroup,nmuon_group_center,nmuon_group_center_noerror,diff);
  }

  if (nmuon_group_center <= 0.0) {
    g2bg_.set_dens_group_lower(igroup, not_found_dens);
    g2bg_.set_dens_group_center(igroup, not_found_dens);
    g2bg_.set_dens_group_upper(igroup, not_found_dens);
    g2bg_.set_delta_nmuon_group_lower(igroup, not_found_delta_nmuon);
    g2bg_.set_delta_nmuon_group_center(igroup, not_found_delta_nmuon);
    g2bg_.set_delta_nmuon_group_upper(igroup, not_found_delta_nmuon);
    return;
  }

  // ! for debug
  const auto [txmin, txmax, tymin, tymax] = g2bg_.get_xmin_xmax_ymin_ymax(igroup);

  // Define upper/lower target values.
  // Add efficiency variance in quadrature to the Poisson term when tf_eff is on.
  const double var_eff
    = tf_eff ? calc_var_eff_group(igroup, tf_eff_independent) : 0.0;
  const double nmuon_group_sigma = sqrt(nmuon_group_center + var_eff);
  double nmuon_group_upper = nmuon_group_center - sigma * nmuon_group_sigma;
  if (nmuon_group_upper < 0.0) nmuon_group_upper = 0.0;
  double nmuon_group_lower = nmuon_group_center + sigma * nmuon_group_sigma;

  // Run adaptive grid search for each target value.
  std::array<double, 2> res_lower
   = calc_dens_nmuon_grouped(
      igroup, dens_min, dens_max,
      vec_dens_step, range_factor,
      ft, nmuon_group_lower);

  std::array<double, 2> res_center
    = calc_dens_nmuon_grouped(
        igroup, dens_min, dens_max,
        vec_dens_step, range_factor,
        ft, nmuon_group_center);

  std::array<double, 2> res_upper
    = calc_dens_nmuon_grouped(
        igroup, dens_min, dens_max,
        vec_dens_step, range_factor,
        ft, nmuon_group_upper);

  double dens_lower  = res_lower[0];
  double dens_center = res_center[0];
  double dens_upper  = res_upper[0];

  double lower_delta_nmuon_min  = res_lower[1];
  double center_delta_nmuon_min = res_center[1];
  double upper_delta_nmuon_min  = res_upper[1];


  // Store results in member variables.
  g2bg_.set_dens_group_lower(igroup, dens_lower);
  g2bg_.set_dens_group_center(igroup, dens_center);
  g2bg_.set_dens_group_upper(igroup, dens_upper);
  g2bg_.set_delta_nmuon_group_lower(igroup, lower_delta_nmuon_min);
  g2bg_.set_delta_nmuon_group_center(igroup, center_delta_nmuon_min);
  g2bg_.set_delta_nmuon_group_upper(igroup, upper_delta_nmuon_min);

  // LOG_DEBUG_ND(
  //   "calc_set_proj_dens_grouped : det_id={:02d} igroup={:5d}, dens_lower={:6.0f}, dens_center={:6.0f}, dens_upper={:6.0f}",
  //   det_id, igroup, dens_lower, dens_center, dens_upper);
}

void DetectorPanel::calc_set_proj_dens_grouped_adaptive_search(
  const Igroup igroup
, const double sigma , const double dens_min, const double dens_max
, const double dens_step_init, const double dens_step_min
, const double dens_step_refinement_factor, const int dens_step_num_refinements
, const FluxTable &ft
, const bool tf_eff, const bool tf_eff_independent )
{
  const int det_id = get_detid();

  // If the group is unavailable, set not_found values.
  if (!g2bg_.is_avail_group(igroup)) {
    g2bg_.set_dens_group_lower(igroup, not_found_dens);
    g2bg_.set_dens_group_center(igroup, not_found_dens);
    g2bg_.set_dens_group_upper(igroup, not_found_dens);
    g2bg_.set_delta_nmuon_group_lower(igroup, not_found_delta_nmuon);
    g2bg_.set_delta_nmuon_group_center(igroup, not_found_delta_nmuon);
    g2bg_.set_delta_nmuon_group_upper(igroup, not_found_delta_nmuon);
    return;
  }

  // Get nmuon_group_center.
  const double nmuon_group_center = g2bg_.get_signal_group(igroup);
  if (nmuon_group_center <= 0.0) {
    g2bg_.set_dens_group_lower(igroup, not_found_dens);
    g2bg_.set_dens_group_center(igroup, not_found_dens);
    g2bg_.set_dens_group_upper(igroup, not_found_dens);
    g2bg_.set_delta_nmuon_group_lower(igroup, not_found_delta_nmuon);
    g2bg_.set_delta_nmuon_group_center(igroup, not_found_delta_nmuon);
    g2bg_.set_delta_nmuon_group_upper(igroup, not_found_delta_nmuon);
    return;
  }

  // Define upper/lower target values.
  // Add efficiency variance in quadrature to the Poisson term when tf_eff is on.
  const double var_eff
    = tf_eff ? calc_var_eff_group(igroup, tf_eff_independent) : 0.0;
  const double nmuon_group_sigma = sqrt(nmuon_group_center + var_eff);
  double nmuon_group_upper = nmuon_group_center - sigma * nmuon_group_sigma;
  if (nmuon_group_upper < 0.0) nmuon_group_upper = 0.0;
  double nmuon_group_lower = nmuon_group_center + sigma * nmuon_group_sigma;

  // Parameters for the adaptive grid search.
  double initial_step = dens_step_init; // Initial step size.
  // Example: set min_step to one tenth of dens_step (adjust as needed).
  double min_step = dens_step_min; // Minimum step size.
  double refinement_factor = dens_step_refinement_factor; // Step shrink factor per iteration.
  int num_refinements = dens_step_num_refinements; // Number of refinements.

  // Run adaptive grid search for each target value.
  std::array<double, 2> res_lower
    = calc_dens_nmuon_grouped_adaptive_search(
      igroup, dens_min, dens_max
    , initial_step, min_step,refinement_factor
    , num_refinements, ft, nmuon_group_lower);

  std::array<double, 2> res_center
    = calc_dens_nmuon_grouped_adaptive_search(
      igroup, dens_min, dens_max
    , initial_step, min_step,refinement_factor
    , num_refinements, ft, nmuon_group_center);
  
  std::array<double, 2> res_upper
    = calc_dens_nmuon_grouped_adaptive_search(
      igroup, dens_min, dens_max
    , initial_step, min_step,refinement_factor
    , num_refinements, ft, nmuon_group_upper);

  double dens_lower  = res_lower[0];
  double dens_center = res_center[0];
  double dens_upper  = res_upper[0];

  double lower_delta_nmuon_min  = res_lower[1];
  double center_delta_nmuon_min = res_center[1];
  double upper_delta_nmuon_min  = res_upper[1];

  // Store results in member variables.
  g2bg_.set_dens_group_lower(igroup, dens_lower);
  g2bg_.set_dens_group_center(igroup, dens_center);
  g2bg_.set_dens_group_upper(igroup, dens_upper);
  g2bg_.set_delta_nmuon_group_lower(igroup, lower_delta_nmuon_min);
  g2bg_.set_delta_nmuon_group_center(igroup, center_delta_nmuon_min);
  g2bg_.set_delta_nmuon_group_upper(igroup, upper_delta_nmuon_min);

  // LOG_DEBUG(
  //   "calc_set_proj_dens_grouped : det_id={:02d} igroup={:5d}, dens_lower={:6.0f}, dens_center={:6.0f}, dens_upper={:6.0f}"
  //   , det_id, igroup, dens_lower, dens_center, dens_upper);
}

// Compute grouped-bin densities and set vec_dens_group_lower/center/upper.
void DetectorPanel::mp_calc_set_proj_dens_grouped(
  const bool tf_signal_poisson
, const double dens_min, const double dens_max, const double dens_step
, const double sigma, const FluxTable &ft
, const bool tf_eff, const bool tf_eff_independent )
{
  // get det_id
  const int det_id = get_detid();

  // get number of groups
  const int num_group = get_dic().get_n_group(det_id);
  LOG_INFO("num_group={}",num_group);

  // rezise vec_dens_group_center, vec_dens_group_upper, vec_dens_group_lower
  g2bg_.resize_vec_dens_group_lower(num_group,not_found_dens);
  g2bg_.resize_vec_dens_group_center(num_group,not_found_dens);
  g2bg_.resize_vec_dens_group_upper(num_group,not_found_dens);

  // rezise vec_delta_nmuon_group_center, vec_delta_nmuon_group_upper, vec_delta_nmuon_group_lower
  g2bg_.resize_vec_delta_nmuon_group_lower(num_group,not_found_delta_nmuon);
  g2bg_.resize_vec_delta_nmuon_group_center(num_group,not_found_delta_nmuon);
  g2bg_.resize_vec_delta_nmuon_group_upper(num_group,not_found_delta_nmuon);

  #pragma omp parallel for
  for(Igroup igroup=0;igroup<num_group;igroup++){
    calc_set_proj_dens_grouped(tf_signal_poisson
      ,igroup,dens_min,dens_max,dens_step,sigma,ft,tf_eff,tf_eff_independent);
  }
}

// Compute grouped-bin densities and set vec_dens_group_lower/center/upper.
void DetectorPanel::mp_calc_set_proj_dens_grouped(
  const bool tf_signal_poisson
, const double dens_min, const double dens_max
, const std::vector<double> &vec_dens_step, const double range_factor
, const double sigma, const FluxTable &ft
, const bool tf_eff, const bool tf_eff_independent )
{
  const int det_id = get_detid();

  // Get the number of groups.
  const int num_group = get_dic().get_n_group(det_id);
  LOG_INFO("detid={}, num_group={}", det_id, num_group);
  if( num_group <= 0 ) {
    LOG_ERROR("No groups available for det_id={}", det_id);
    THROW_ERROR("num_group <= 0");
  }

  // Resize result vectors.
  g2bg_.resize_vec_dens_group_lower(num_group, not_found_dens);
  g2bg_.resize_vec_dens_group_center(num_group, not_found_dens);
  g2bg_.resize_vec_dens_group_upper(num_group, not_found_dens);

  g2bg_.resize_vec_delta_nmuon_group_lower(num_group, not_found_delta_nmuon);
  g2bg_.resize_vec_delta_nmuon_group_center(num_group, not_found_delta_nmuon);
  g2bg_.resize_vec_delta_nmuon_group_upper(num_group, not_found_delta_nmuon);

  #pragma omp parallel for
  for (Igroup igroup = 0; igroup < num_group; igroup++) {
    fprintf(stderr
      , "Processing calc_set_proj_dens_grouped det_id=%2d, igroup=%5d/%5d\r"
      , det_id, igroup, num_group);
    calc_set_proj_dens_grouped(
        tf_signal_poisson, igroup, sigma, dens_min, dens_max
      , vec_dens_step, range_factor, ft, tf_eff, tf_eff_independent);
  }
}

// Compute grouped-bin densities and set vec_dens_group_lower/center/upper.
void DetectorPanel::mp_calc_set_proj_dens_grouped_adaptive_search(
  const double dens_min, const double dens_max
, const double dens_step_init, const double dens_step_min
, const double dens_step_refinement_factor, const int dens_step_num_refinements
, const double sigma, const FluxTable &ft
, const bool tf_eff, const bool tf_eff_independent )
{
  // get det_id
  const int det_id = get_detid();

  // get number of groups
  const int num_group = get_dic().get_n_group(det_id);
  LOG_INFO("num_group={}",num_group);

  // rezise vec_dens_group_center, vec_dens_group_upper, vec_dens_group_lower
  g2bg_.resize_vec_dens_group_lower(num_group,not_found_dens);
  g2bg_.resize_vec_dens_group_center(num_group,not_found_dens);
  g2bg_.resize_vec_dens_group_upper(num_group,not_found_dens);

  // rezise vec_delta_nmuon_group_center, vec_delta_nmuon_group_upper, vec_delta_nmuon_group_lower
  g2bg_.resize_vec_delta_nmuon_group_lower(num_group,not_found_delta_nmuon);
  g2bg_.resize_vec_delta_nmuon_group_center(num_group,not_found_delta_nmuon);
  g2bg_.resize_vec_delta_nmuon_group_upper(num_group,not_found_delta_nmuon);

  #pragma omp parallel for
  for(Igroup igroup=0;igroup<num_group;igroup++){
    calc_set_proj_dens_grouped_adaptive_search(
      igroup, sigma, dens_min, dens_max
    , dens_step_init, dens_step_min
    , dens_step_refinement_factor, dens_step_num_refinements
    , ft, tf_eff, tf_eff_independent);
  }
}

// calc_signal_group without efficiency
double DetectorPanel::mp_calc_signal_noeff_group(const Igroup igroup) const
{
  const Detid detid_tmp = get_detid();
  const int num_group = get_dic().get_n_group(detid_tmp);
  if (igroup < 0 || igroup >= num_group) {
    THROW_ERROR_NAME("Invalid group index: " + std::to_string(igroup));
  }

  const auto& vec_ixiy = g2bg_.get_vec_ixiy(igroup);
  std::vector<const DetectorElement*> vec_elem;
  vec_elem.reserve(vec_ixiy.size());
  for (int i = 0; i < (int)vec_ixiy.size(); ++i) {
    auto [ix, iy] = vec_ixiy[i];
    const DetectorElement& ele = getDetectorElement(ix, iy);
    vec_elem.push_back(&ele);
  }

  double signal_group = 0.0;
  #pragma omp parallel for schedule(static) reduction(+:signal_group)
  for (int i = 0; i < (int)vec_elem.size(); ++i) {
    signal_group += vec_elem.at(i)->calc_signal();
  }
  return signal_group;
}


// calc_signal_group with efficiency
double DetectorPanel::mp_calc_signal_eff_group(const Igroup igroup) const
{
  const int num_group = get_dic().get_n_group(detid_);
  if (igroup < 0 || igroup >= num_group) {
    THROW_ERROR_NAME("Invalid group index: " + std::to_string(igroup));
  }

  const auto& vec_ixiy = g2bg_.get_vec_ixiy(igroup);
  std::vector<const DetectorElement*> vec_elem;
  vec_elem.reserve(vec_ixiy.size());
  for (int i = 0; i < (int)vec_ixiy.size(); ++i) {
    auto [ix, iy] = vec_ixiy[i];
    const auto& ele = getDetectorElement(ix, iy);
    vec_elem.push_back(&ele);
  }

  double signal_group = 0.0;
  #pragma omp parallel for schedule(static) reduction(+:signal_group)
  for (int i = 0; i < (int)vec_elem.size(); ++i) {
    signal_group += vec_elem.at(i)->calc_signal_eff();
  }
  return signal_group;
}


// simple unweighted mean of element efficiency triplets over one group
std::array<double,3> DetectorPanel::calc_eff_mean_group(const Igroup igroup) const
{
  // No get_n_group() range check here: countOne() rescans all pairs per call,
  // which made this O(n_group * n_element) (id-c8hy15). An out-of-range igroup
  // still throws below because get_vec_ixiy() returns an empty vector for it.
  if (igroup < 0) {
    THROW_ERROR_NAME("Invalid group index: " + std::to_string(igroup));
  }
  const auto& vec_ixiy = g2bg_.get_vec_ixiy(igroup);
  if (vec_ixiy.empty()) {
    THROW_ERROR_NAME("Group has no elements: igroup=" + std::to_string(igroup));
  }

  double sum_low = 0.0, sum_cnt = 0.0, sum_upp = 0.0;
  for (int i = 0; i < (int)vec_ixiy.size(); ++i) {
    auto [ix, iy] = vec_ixiy[i];
    const DetectorElement& ele = getDetectorElement(ix, iy);
    sum_low += ele.get_eff_low();
    sum_cnt += ele.get_eff_cnt();
    sum_upp += ele.get_eff_upp();
  }
  const double n_ele = static_cast<double>(vec_ixiy.size());
  return { sum_low/n_ele, sum_cnt/n_ele, sum_upp/n_ele };
}



// output header
void DetectorPanel::out_header( FILE *fout ) const
{
  fprintf(fout,"# %d %s ",get_detid(), name.c_str());
  fprintf(fout,"%d %7.4lf %7.4lf ", g2bg_.get_nbinx(), g2bg_.get_xmin(), g2bg_.get_xmax());
  fprintf(fout,"%d %7.4lf %7.4lf ", g2bg_.get_nbiny(), g2bg_.get_ymin(), g2bg_.get_ymax());

  fprintf(fout,"%7.4lf %7.4lf %7.4lf "
  , ray3d.vx()
  , ray3d.vy()
  , ray3d.vz() );
  fprintf(fout,"%.2lf %.2lf %.2lf "
  , ray3d.x()
  , ray3d.y()
  , ray3d.z() );
  fprintf(fout,"%.2lf %.2lf %.2lf "
  , v3_det_length.x()
  , v3_det_length.y()
  , v3_det_length.z() );
  fprintf(fout,"%.4lf\n",days);
}

// out_all
void DetectorPanel::out_all( FILE *fout ) const
{
  // output header info
  out_header(fout);
  
  LOG_DEBUG("vec_vec_DetectorElement.size()={}"
  ,vec_vec_DetectorElement.size());

  //output vecvec_DetectorElement
  int ix,iy;
  for(iy=0;iy<get_nbiny();iy++){
    for(ix=0;ix<get_nbinx();ix++){
      // const DetectorElement& ele = getDetectorElement(ix,iy);
      getDetectorElement(ix,iy).out(fout);
    }
  }
}

void DetectorPanel::out_all( const std::filesystem::path& pathout ) const
{
  LOG_INFO("outputing {} ...",pathout.string());
  if( pathout.empty() ) THROW_ERROR_NAME("fname is not assigned");
  FILE *fout = fopen(pathout.c_str(),"wt");
  if(fout == NULL) THROW_ERROR_NAME("Cannot open file : " + pathout.string());
  out_all(fout);
  myapp::close(fout,pathout);
}

// output txcnt, tycnt, signal to FILE *fout
void DetectorPanel::out_axay_signal( FILE *fout ) const
{
  int ix,iy;
  for(iy=0;iy<get_nbiny();iy++){
    for(ix=0;ix<get_nbinx();ix++){
      const DetectorElement& ele = getDetectorElement(ix,iy);
      const double txcnt = ele.get_tx();
      const double tycnt = ele.get_ty();
      const double signal = ele.get_signal();
      fprintf(fout,"%7.4lf %7.4lf %E\n",txcnt,tycnt,signal);
    }
  }
}

// output txcnt, tycnt, signal to std::filesystem::path
void DetectorPanel::out_axay_signal( const std::filesystem::path& pathout ) const
{
  LOG_INFO("outputing {} ...",pathout.string());
  if( pathout.empty() ) THROW_ERROR_NAME("fname is not assigned");
  FILE *fout = fopen(pathout.c_str(),"wt");
  if(fout == NULL) THROW_ERROR_NAME("Cannot open file : " + pathout.string());
  out_axay_signal(fout);
  myapp::close(fout,pathout);
}

// Output (tx, ty, arbitrary data) for all DetectorElements.
void DetectorPanel::out_txtyData(
    const std::string& suffix
  , std::function<double(const DetectorElement&)> data_extractor
  , const std::string& prefix ) const
{
  char fname[512];
  std::snprintf( fname, sizeof(fname)
    , "%s%s_txty%s_det%02d.tmp"
    , prefix.c_str(), get_name().c_str(), suffix.c_str()
    , static_cast<int>(get_detid())
  );

  fs::path pathout = iodir::make_pathout(fname);
  FILE* fout = std::fopen(pathout.c_str(), "wt");
  if (!fout) {
    THROW_ERROR_NAME("Cannot open file : " + pathout.string());
  }

  for (int iy = 0; iy < get_nbiny(); ++iy) {
    for (int ix = 0; ix < get_nbinx(); ++ix) {
      const DetectorElement& ele = getDetectorElement(ix, iy);
      const double tx = ele.get_tx();
      const double ty = ele.get_ty();
      const double data = data_extractor(ele);
      std::fprintf(fout, "%9.6lf %9.6lf %E\n", tx, ty, data);
    }
  }

  myapp::close(fout, pathout);
}

void DetectorPanel::out_txtySignal() const {
  out_txtyData("_signal",
    [](const DetectorElement& ele) { return ele.get_signal(); });
}

void DetectorPanel::out_txtySignal(const std::string& prefix) const {
  out_txtyData("_signal",
    [](const DetectorElement& ele) { return ele.get_signal(); },
    prefix);
}

void DetectorPanel::out_txtyPL() const {
  out_txtyData("_PL",
    [](const DetectorElement& ele) { return ele.get_PL(); });
}

void DetectorPanel::out_txtyPL(const std::string& prefix) const {
  out_txtyData("_PL",
    [](const DetectorElement& ele) { return ele.get_PL(); },
    prefix);
}

void DetectorPanel::out_txtyDL() const {
  out_txtyData("_DL",
    [](const DetectorElement& ele) { return ele.get_DL(); });
}

void DetectorPanel::out_txtyDL(const std::string& prefix) const {
  out_txtyData("_DL",
    [](const DetectorElement& ele) { return ele.get_DL(); },
    prefix);
}

// read efficinecy table from file
void DetectorPanel::read_efficiency_table( const std::filesystem::path &path_in )
{
  LOG_DEBUG("Reading efficiency table from {} ...", path_in.string());
  std::ifstream ifs(path_in);
  if (!ifs.is_open()) {
    THROW_ERROR_NAME("Failed to open file: " + path_in.string());
  }
  std::string line;
  size_t line_no = 0;
  size_t warn_count = 0;  
  size_t count_set = 0;

  while (std::getline(ifs, line)) {
    ++line_no;
    if (line.empty()) continue;

    std::istringstream iss(line);
    EffBin bin;
    if (!(iss >> bin.xlow >> bin.xup >> bin.ylow >> bin.yup
              >> bin.eff_low >> bin.eff_cnt >> bin.eff_upp)) {
      THROW_ERROR_NAME("Format error in line " + std::to_string(line_no));
    }
    const double xcnt = (bin.xlow + bin.xup) * 0.5;
    const double ycnt = (bin.ylow + bin.yup) * 0.5;

    // if xcnt, ycnt is out of range, skip
    if (xcnt < g2bg_.get_xmin() || xcnt > g2bg_.get_xmax() ||
        ycnt < g2bg_.get_ymin() || ycnt > g2bg_.get_ymax()) {
      if (warn_count < 5) {
        LOG_WARN("xcnt or ycnt out of range in line {}: xcnt={}, ycnt={}", line_no, xcnt, ycnt);
        ++warn_count;
      }
      continue;
    }

    DetectorElement& ele = callDetectorElement(xcnt, ycnt);
    ele.set_eff_low(bin.eff_low);
    ele.set_eff_cnt(bin.eff_cnt);
    ele.set_eff_upp(bin.eff_upp);
    ++count_set;
  }
  LOG_INFO("Set efficiency for {} elements from {}", count_set, path_in.string());
}

// Evaluate the efficiency model at every bin center of this panel's own grid
// (g2bg_), so no external table and no grid mismatch is possible.
// Uses OpenMP: each (ix, iy) writes only to its own DetectorElement (no race).
void DetectorPanel::mp_assign_efficiency_model( const EfficiencyModel &model )
{
  if (!model.get_tf_loaded()) {
    THROW_ERROR_NAME("EfficiencyModel is not loaded");
  }
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  const Grid1d &x_axis = get_x_axis();
  const Grid1d &y_axis = get_y_axis();

  #pragma omp parallel for collapse(2) schedule(static)
  for (int iy = 0; iy < nbiny; iy++) {
    for (int ix = 0; ix < nbinx; ix++) {
      const double tx = x_axis.get_center_value(ix);
      const double ty = y_axis.get_center_value(iy);
      const EfficiencyModel::Band band = model.eval(tx, ty);
      DetectorElement &ele = callDetectorElement(ix, iy);
      ele.set_eff_low(band.eff_low);
      ele.set_eff_cnt(band.eff_cnt);
      ele.set_eff_upp(band.eff_upp);
    }
  }
  LOG_INFO("Set efficiency for {} elements from the efficiency model", nbinx * nbiny);
}

// Get summed signal over angle-bin regions defined by BinList.
Eigen::VectorXf DetectorPanel::get_vecxf_nmuoned() const
{
  const int num_groups = get_dic().get_n_group(detid_);
  Eigen::VectorXf vecxf_merged_signal(num_groups);

  // bin group loop
  for(int imerged=0;imerged<num_groups;imerged++){
    // initilize signal_group
    double signal_group = 0.0;

    // pair loop
    const std::vector<Ixiy> vec_ixiy = g2bg_.get_vec_ixiy(imerged);
    for(const auto [ix,iy] : vec_ixiy){
      const DetectorElement& ele = getDetectorElement(ix,iy);
      double signal = ele.get_signal();
      fprintf(stderr,"ix=%4d, iy=%4d, signal=%E\r",ix,iy,signal);
      signal_group += signal;
    }

    // store signal_group to vecxf_merged_signal
    vecxf_merged_signal(imerged) =  signal_group;
  }
  return vecxf_merged_signal;
}

// set direction from Eigen::Vector2d v2_pos
void DetectorPanel::set_direction_to_v2_pos( const Eigen::Vector2d v2_pos )
{
  const Eigen::Vector3d v3_dir_src = this->get_ray3d().dir();
  const Eigen::Vector3d v3_dir_dst = geom_util::calc_v3_dir(
                                      this->get_ray3d().pos2d(), v2_pos );
  this->set_direction( v3_dir_dst );

  const Eigen::Matrix3d rotmat
    = geom_util::calc_rotmat_from_v3_to_v3_z_dominant( v3_dir_src, v3_dir_dst );

  // parallel for loop of ix, iy
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  #pragma omp parallel for collapse(2) schedule(static)
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      DetectorElement& ele = callDetectorElement(ix,iy);
      const Eigen::Vector3d v3_dir_src_ele = ele.get_ray3d().dir();
      const Eigen::Vector3d v3_dir_dst_ele = rotmat * v3_dir_src_ele;
      ele.set_v3_direction( v3_dir_dst_ele );
    }
  }
}


// set the exposure time (sec) for all DetectorElement
void DetectorPanel::set_exposure_time_sec( const double time_sec_in )
{
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  #pragma omp parallel for collapse(2) schedule(static)
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      DetectorElement& ele = callDetectorElement(ix,iy);
      ele.set_exposure_time_sec(time_sec_in);
    }
  }
}

// copy only DetectortElement::vec_tf_in_PL
void DetectorPanel::mp_set_vec_tf_in_PL( const DetectorPanel& panel_in )
{
  const int nbinx = this->get_nbinx();
  const int nbiny = this->get_nbiny();
  #pragma omp parallel for collapse(2) schedule(static)
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      DetectorElement& ele = this->callDetectorElement(ix,iy);
      const DetectorElement& ele_in = panel_in.getDetectorElement(ix,iy);
      ele.set_vec_tf_in_PL(ele_in.get_vec_tf_in_PL());
    }
  }
}

void DetectorPanel::mp_calc_set_peneflux_signal_from_DL(
  const FluxTable &ft, const bool tf_apply_eff
, const bool tf_apply_eff_cnt )
{
  // flux table
  const Grid2dXYZ g2flux = ft.get_g2_log_peneflux();

  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  const int n_threads = omp_get_max_threads();
  #pragma omp parallel for collapse(2) schedule(static)
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      if(ix%n_threads*10==0) fprintf(stderr
          ,"DetectorPanel::mp_calc_set_peneflux_signal_from_DL ix=%d, iy=%d....\r"
          ,ix,iy);
      // call DetectorElement pointer
      DetectorElement& ele = callDetectorElement(ix,iy);

      // assign peneflux to *ele
      const double peneflux = ele.calc_peneflux(g2flux);
      ele.set_peneflux(peneflux);

      // calc and assign solid angle to *ele
      const double omega = ele.calc_solid_angle();
      ele.set_solid_angle(omega);
      
      // calc & assign effective_area to *ele
      const double effective_area = ele.calc_effective_area(
        this->get_v3_dir(),this->get_v3_det_length(),this->get_n_unit()
      );
      ele.set_effective_area_m2(effective_area);

      double signal = 0.0;
      if(tf_apply_eff_cnt){
        // central deterministic efficiency: signal = b_i * eff_cnt (not the dice path)
        signal = ele.calc_signal() * ele.get_eff_cnt();
      }else if(tf_apply_eff){
        signal = ele.calc_signal_eff();
      }else{
        signal = ele.calc_signal();
      }

      ele.set_signal(signal);
    }
  }
}

void DetectorPanel::mp_set_noise(
  const double flux_proport_ratio_floor
, const double flux_proport_ratio_poisson
, const double SOT_proport_noise_ratio_floor
, const double SOT_proport_noise_ratio_poisson
, const double user_defined_noise_flux_ratio
, const std::vector<std::filesystem::path> vec_path_user_defined_noise_flux
, const double DL_thres)
{
  LOG_INFO("detid_={:02d} flux_floor={:E}, flux_poisson={:E}, SOT_floor={:E}, SOT_poisson={:E}, user_defined_noise_flux_ratio={:E}"
    ,get_detid(),flux_proport_ratio_floor,flux_proport_ratio_poisson
    ,SOT_proport_noise_ratio_floor,SOT_proport_noise_ratio_poisson,user_defined_noise_flux_ratio);
  if(  flux_proport_ratio_floor <= 0.0
    && flux_proport_ratio_poisson <= 0.0
    && SOT_proport_noise_ratio_floor <= 0.0
    && SOT_proport_noise_ratio_poisson <= 0.0
    && user_defined_noise_flux_ratio <= 0.0){
    LOG_WARN("all noise ratios are not positive: flux_floor={}, flux_poisson={}, SOT_floor={}, SOT_poisson={}, user_defined={}"
      ,flux_proport_ratio_floor,flux_proport_ratio_poisson
      ,SOT_proport_noise_ratio_floor,SOT_proport_noise_ratio_poisson,user_defined_noise_flux_ratio);
    return;
  }

  LOG_INFO("calc total signal (DL>DLthes={}) in this DetectorPanel id={:02d} : {}"
    ,DL_thres,get_detid(),this->get_name());
  const double total_signal_DL_thres = this->mp_calc_total_signal_DLthres(DL_thres);

  LOG_INFO("calc total SOT (DL>DLthes={}) in this DetectorPanel id={:02d} : {}"
    ,DL_thres,get_detid(),this->get_name());
  const double total_SOT_DL_thres = this->mp_calc_total_SOT_DLthres(DL_thres);

  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();

  #pragma omp parallel for collapse(2) schedule(static)
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      DetectorElement& ele = callDetectorElement(ix,iy);

      double DL = ele.get_DL();
      if( DL < DL_thres ) continue;

      double signal = ele.get_signal();
      double SOT = ele.get_SOT();

      // flux-proportional parts (scale with signal): floor vs poisson bucket
      double flux_floor_noise = 0.0;
      if( flux_proport_ratio_floor > 0.0 ){
        flux_floor_noise = signal * flux_proport_ratio_floor;
      }
      double flux_poisson_noise = 0.0;
      if( flux_proport_ratio_poisson > 0.0 ){
        flux_poisson_noise = signal * flux_proport_ratio_poisson;
      }

      // SOT-proportional (angle-independent) parts: floor vs poisson bucket
      double SOT_floor_noise = 0.0;
      if( SOT_proport_noise_ratio_floor > 0.0 ){
        SOT_floor_noise
          = SOT/total_SOT_DL_thres * SOT_proport_noise_ratio_floor * total_signal_DL_thres;
      }
      double SOT_poisson_noise = 0.0;
      if( SOT_proport_noise_ratio_poisson > 0.0 ){
        SOT_poisson_noise
          = SOT/total_SOT_DL_thres * SOT_proport_noise_ratio_poisson * total_signal_DL_thres;
      }

      // user-defined noise (reserved; not implemented yet 2025-06-03)
      double user_defined_noise = 0.0;

      // det = deterministic floor sum (no fluctuation); poi = Poisson-bucket sum
      const double noise_det = flux_floor_noise + SOT_floor_noise;
      const double noise_poi = flux_poisson_noise + SOT_poisson_noise + user_defined_noise;
      ele.set_noise_det(noise_det);
      ele.set_noise_poi(noise_poi);
    }
  }
}

// initilize path length
void DetectorPanel::mp_initPL( const double value_in )
{
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  #pragma omp parallel for collapse(2) schedule(static)
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      DetectorElement& ele = callDetectorElement(ix,iy);
      ele.set_PL(value_in);
    }
  }
}

// initilize DL
void DetectorPanel::mp_initDL( const double value_in )
{
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  #pragma omp parallel for collapse(2) schedule(static)
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      DetectorElement& ele = callDetectorElement(ix,iy);
      ele.set_DL(value_in);
    }
  }
}

// initilize PL and PL for all DetectorElement
void DetectorPanel::mp_init_PLDL( const double PL_in, const double DL_in)
{
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  #pragma omp parallel for collapse(2) schedule(static)
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      DetectorElement& ele = callDetectorElement(ix,iy);
      ele.set_PL(PL_in);
      ele.set_DL(DL_in);
    }
  }
}

// initilize peneflux
void DetectorPanel::mp_init_peneflux( const double value_in )
{
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  #pragma omp parallel for collapse(2) schedule(static)
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      DetectorElement& ele = callDetectorElement(ix,iy);
      ele.set_peneflux(value_in);
    }
  }
}

// initilize signal
void DetectorPanel::mp_init_signal( const double value_in )
{
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  #pragma omp parallel for collapse(2) schedule(static)
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      DetectorElement& ele = callDetectorElement(ix,iy);
      ele.set_signal(value_in);
    }
  }
}

// initilize len, DL, pene, sig, noi for all DetectorElement
void DetectorPanel::mp_init_PL_DL_pene_sig_noi(
    const double PL_in, const double DL_in
  , const double peneflux_in
  , const double sig_in, const double noi_in )
{
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  #pragma omp parallel for collapse(2) schedule(static)
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      DetectorElement& ele = callDetectorElement(ix,iy);
      ele.set_PL(PL_in);
      ele.set_DL(DL_in);
      ele.set_peneflux(peneflux_in);
      ele.set_signal(sig_in);
      ele.set_noise(noi_in);
    }
  }
}


// initilize noise
void DetectorPanel::mp_init_noise( const double value_in )
{
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      DetectorElement& ele = callDetectorElement(ix,iy);
      ele.set_noise(value_in);
    }
  }
}

// @brief get set of uqid from igroup
std::set<Uqid> DetectorPanel::get_set_uqid_by_igroup( const Igroup igroup ) const
{
  std::set<Uqid> set_uqid;
  const std::vector<Ixiy> vec_ixiy = g2bg_.get_vec_ixiy(igroup);
  for(const auto [ix,iy] : vec_ixiy){
    const DetectorElement& ele = getDetectorElement(ix,iy);
    const Uqid uqid = ele.get_uqid();
    set_uqid.insert(uqid);
  }
  return set_uqid;
}

// @brief get vector of uqid from igroup
std::vector<Uqid> DetectorPanel::get_vec_uqid_by_igroup( const Igroup igroup ) const
{
  const std::vector<Ixiy> vec_ixiy = g2bg_.get_vec_ixiy(igroup);
  std::vector<Uqid> vec_uqid;
  vec_uqid.resize(vec_ixiy.size());
  #pragma omp parallel for
  for(int i=0;i<vec_ixiy.size();i++){
    const auto [ix,iy] = vec_ixiy.at(i);
    const DetectorElement& ele = getDetectorElement(ix,iy);
    const Uqid uqid = ele.get_uqid();
    vec_uqid.at(i) = uqid;
  }
  return vec_uqid;
}

// return std::vector<std::set<unique_index>>, the index of vector is igroup
// with transforming ixiy into uqid
std::vector<std::set<Uqid>> DetectorPanel::get_vec_set_uqid() const
{
  std::vector<std::set<Uqid>> vec_set_uqid;
  std::vector<Uqid> vec_igroup = get_dic().get_vec_igroup(detid_);
  const size_t n_igroup = vec_igroup.size();
  LOG_DEBUG("name={}, n_igroup={}",get_name(),n_igroup);

  for( const auto &igroup : vec_igroup ){
    std::set<Uqid> set_uqid = get_set_uqid_by_igroup(igroup);
    vec_set_uqid.push_back(set_uqid);
    const int n_uqid = set_uqid.size();
    // mylogger::g_logger->trace("igroup={}, n_uqid={}",igroup,n_uqid);
  }
  return vec_set_uqid;
};

// get Eigen::VectorXf of PL, index is id_in_this_detector
Eigen::VectorXf DetectorPanel::mp_get_vecxf_PL() const
{
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  Eigen::VectorXf vecxf_PL(nbinx*nbiny);
  #pragma omp parallel for collapse(2) schedule(static)
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      const DetectorElement& ele = getDetectorElement(ix,iy);
      const double PL = ele.get_PL();
      const int inthis = ele.get_id_in_this_detector();
      vecxf_PL(inthis) =  PL;
    }
  }
  return vecxf_PL;
}

// get Eigen::VectorXf of DL, index is id_in_this_detector
Eigen::VectorXf DetectorPanel::mp_get_vecxf_DL() const
{
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  Eigen::VectorXf vecxf_DL(nbinx*nbiny);
  #pragma omp parallel for collapse(2) schedule(static)
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      const DetectorElement& ele = getDetectorElement(ix,iy);
      const double DL = ele.get_DL();
      const int inthis = ele.get_id_in_this_detector();
      vecxf_DL(inthis) = DL;
    }
  }
  return vecxf_DL;
}

// @brief get the signal values of all DetectorElement as Eigen::VectorXf in order of id_in_this_detector
Eigen::VectorXf DetectorPanel::get_vecxf_signal() const
{
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  Eigen::VectorXf vecxf_signal(nbinx*nbiny);
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      const DetectorElement& ele = getDetectorElement(ix,iy);
      const double signal = ele.get_signal();
      const int inthis = ele.get_id_in_this_detector();
      vecxf_signal(inthis) =  signal;
    }
  }
  return vecxf_signal;
}

// @brief get the noise values of all DetectorElement as Eigen::VectorXf in order of id_in_this_detector
Eigen::VectorXf DetectorPanel::get_vecxf_noise() const
{
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  Eigen::VectorXf vecxf_noise(nbinx*nbiny);
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      const DetectorElement& ele = getDetectorElement(ix,iy);
      const double noise = ele.get_noise();
      const int inthis = ele.get_id_in_this_detector();
      vecxf_noise(inthis) =  noise;
    }
  }
  return vecxf_noise;
}

// get the diagonal SpMatf of efficiency with index of (id_inthis,id_inthis)
SpMatf DetectorPanel::get_spmat_eff_cnt() const
{
  using Trip = Eigen::Triplet<float>;
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  const int n_elem = nbinx * nbiny;

  // Prepare Triplet buffer.
  std::vector<Trip> triplets;
  triplets.reserve(n_elem);

  for (int iy = 0; iy < nbiny; ++iy) {
    for (int ix = 0; ix < nbinx; ++ix) {
      const DetectorElement& ele = getDetectorElement(ix, iy);
      const float  eff           = static_cast<float>(ele.get_eff_cnt());
      const Inthis  id_inthis    = ele.get_id_in_this_detector();
      triplets.emplace_back(id_inthis, id_inthis, eff);
    }
  }

  // Convert to sparse matrix.
  SpMatf spmat_eff(n_elem, n_elem);
  spmat_eff.setFromTriplets(triplets.begin(), triplets.end());

  return spmat_eff;
}

// get the diagonal SpMatf of eff(random sampling value) with index of (id_inthis,id_inthis)
SpMatf DetectorPanel::get_spmat_eff_sample() const
{
  using Trip = Eigen::Triplet<float>;
  const int nbinx    = get_nbinx();
  const int nbiny    = get_nbiny();
  const int n_elem   = nbinx * nbiny;

  // Prepare Triplet buffer.
  std::vector<Trip> triplets;
  triplets.reserve(n_elem);

  for (int iy = 0; iy < nbiny; ++iy) {
    for (int ix = 0; ix < nbinx; ++ix) {
      const DetectorElement& ele = getDetectorElement(ix, iy);
      const float           eff_sample = static_cast<float>(ele.calc_eff_sample());
      const Inthis           id_inthis = ele.get_id_in_this_detector();
      triplets.emplace_back(id_inthis, id_inthis, eff_sample);
    }
  }

  // Convert to sparse matrix.
  SpMatf spmat_eff(n_elem, n_elem);
  spmat_eff.setFromTriplets(triplets.begin(), triplets.end());
  return spmat_eff;
}

// @brief get the total signal of all DetectorElement as double
double DetectorPanel::mp_calc_total_signal_DLthres(const double DL_thres) const
{
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  double total_signal = 0.0;
  #pragma omp parallel for collapse(2) reduction(+:total_signal) schedule(static)
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      const DetectorElement& ele = getDetectorElement(ix,iy);
      double DL = ele.get_DL();
      if( DL < DL_thres ) continue;
      total_signal += ele.get_signal();
    }
  }
  return total_signal;
}

// @brief get the total SOT of all DetectorElement as double
double DetectorPanel::mp_calc_total_SOT_DLthres(const double DL_thres) const
{
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  double total_SOT = 0.0;
  #pragma omp parallel for collapse(2) reduction(+:total_SOT) schedule(static)
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      const DetectorElement& ele = getDetectorElement(ix,iy);
      double DL = ele.get_DL();
      if( DL < DL_thres ) continue;
      total_SOT += ele.get_SOT();
    }
  }
  return total_SOT;
}

// get_subtract of *this and other DetectorPanel
DetectorPanel DetectorPanel::get_subtract(const DetectorPanel &other) const
{
  if( this->get_nbinx() != other.get_nbinx() || this->get_nbiny() != other.get_nbiny() ){
    THROW_ERROR_NAME("DetectorPanel::get_subtract: nbinx or nbiny is not same");
  }
  DetectorPanel result(*this);
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();

  // subtract signal and noise of each DetectorElement, and vec_vec_signal, vec_vec_noise_poi
  #pragma omp parallel for collapse(2) schedule(static)
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      DetectorElement& ele = result.callDetectorElement(ix,iy);
      const DetectorElement& ele_other = other.getDetectorElement(ix,iy);
      ele.set_signal(ele.get_signal() - ele_other.get_signal());
      // subtract det/poi noise channels separately to preserve the floor/poisson split
      ele.set_noise_det(ele.get_noise_det() - ele_other.get_noise_det());
      ele.set_noise_poi(ele.get_noise_poi() - ele_other.get_noise_poi());
    }
  }

  // subtract of vec_signal_group, vec_noise_poi_group,  in base class Grid2dBinGroup
  auto [subtracted_vec_vec_signal, subtracted_vec_vec_noise_poi]
    = g2bg_.get_vec_vec_signal_noise_subtracted(other.get_g2bg());

  auto [subtracted_vec_signal_group, subtracted_vec_noise_poi_group]
    = g2bg_.get_vec_signal_noise_group_subtracted(other.get_g2bg());

  result.g2bg_.set_vec_vec_signal(subtracted_vec_vec_signal);
  result.g2bg_.set_vec_vec_noise_poi(subtracted_vec_vec_noise_poi);
  result.g2bg_.set_vec_signal_group(subtracted_vec_signal_group);
  result.g2bg_.set_vec_noise_poi_group(subtracted_vec_noise_poi_group);

  // subtract the deterministic floor (det) noise channel: per-bin then per-group
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      const double noise_det_sub
        = g2bg_.get_noise_det(ix,iy) - other.get_g2bg().get_noise_det(ix,iy);
      result.g2bg_.set_noise_det(ix,iy,noise_det_sub);
    }
  }
  const int n_group_det = result.g2bg_.get_n_group();
  result.g2bg_.resize_vec_noise_det_group(n_group_det, 0.0);
  for(int ig=0; ig<n_group_det; ig++){
    result.g2bg_.set_noise_det_group(ig, result.g2bg_.calc_noise_det_group(ig));
  }

  return result;
}

// get the sum of signal in the specified range
SignalNoiseSum 
  DetectorPanel::mp_get_signal_noise_sum_range_xy(
    const double tx_lower,const double tx_upper
  , const double ty_lower,const double ty_upper ) const
{
  if( tx_upper < tx_lower ) THROW_ERROR3("tx_upper must be >= tx_lower", tx_upper, tx_lower);
  if( ty_upper < ty_lower ) THROW_ERROR3("ty_upper must be >= ty_lower", ty_upper, ty_lower);

  SignalNoiseSum result;

  const double eps = 1.0e-6;
  const double tx_eps = this->get_x_interval() * eps;
  const double ty_eps = this->get_y_interval() * eps;

  const double tx_lower_eps = tx_lower + tx_eps;
  const double tx_upper_eps = tx_upper - tx_eps;
  const double ty_lower_eps = ty_lower + ty_eps;
  const double ty_upper_eps = ty_upper - ty_eps;

  result.txty_range.set_xmin(tx_lower_eps);
  result.txty_range.set_xmax(tx_upper_eps);
  result.txty_range.set_ymin(ty_lower_eps);
  result.txty_range.set_ymax(ty_upper_eps);

  // Get the range of bins in x and y directions
  const int ix_lower = get_ix(tx_lower_eps);
  const int ix_upper = get_ix(tx_upper_eps);
  const int iy_lower = get_iy(ty_lower_eps);
  const int iy_upper = get_iy(ty_upper_eps);

  // Loop over the bins and accumulate the signal & noise
  // - collapse(2) flattens the two nested loops
  // - reduction on both signal_sum and noise_sum to avoid data races
  double signal_sum = 0.0;
  double noise_sum  = 0.0;
  #pragma omp parallel for collapse(2) reduction(+:signal_sum, noise_sum) schedule(static)
  for (int iy = iy_lower; iy <= iy_upper; ++iy) {
    for (int ix = ix_lower; ix <= ix_upper; ++ix) {
      const DetectorElement& ele = getDetectorElement(ix, iy);
      signal_sum += ele.get_signal();
      noise_sum  += ele.get_noise();
    }
  }
  result.signal = signal_sum;
  result.noise  = noise_sum;

  return result;
}

// get the sum of signal in the specified range (y-sweep)
std::vector<SignalNoiseDepth> 
  DetectorPanel::mp_get_signal_noise_sum_y_sweep(
      const double tx_lower, const double tx_upper
    , const double ty_lower, const double ty_upper
    , const double ty_min, const double ty_max, const double ty_step
    , const Eigen::Vector3d& v3_pos_obj_top
    , const double diff_elev ) const
{
  std::vector<SignalNoiseDepth> vec_result;

  const Detid detid = this->get_detid();
  const Eigen::Vector3d v3_pos_det = this->get_ray3d().pos();
  const Eigen::Vector2d v2_pos_det = v3_pos_det.head<2>();
  const Eigen::Vector2d v2_pos_obj = v3_pos_obj_top.head<2>();
  const Eigen::Vector2d v2_delta = v2_pos_obj - v2_pos_det;
  const double delta_horizontal = v2_delta.norm();
  const double elev_obj_from_det = v3_pos_obj_top.z() - v3_pos_det.z();
  
  // calc y_tan_cnt, y_tan_low, y_tan_up
  const angle_util::AngleUnit angle_unit = this->get_angle_unit();

  // Safety checks.
  if (ty_step <= 0.0) THROW_ERROR("ty_step must be positive");
  if (ty_max < ty_min) THROW_ERROR("ty_max must be greater than or equal to ty_min");

  // Number of sweep steps (from upper to lower). Use 1 step for zero width.
  const double ty_len = ty_upper - ty_lower;
  const double sweep_y_length = (ty_max-0.5*ty_len) - (ty_min+0.5*ty_len);
  int nsteps = static_cast<int>(std::floor(sweep_y_length / ty_step));
  if (nsteps <= 0) nsteps = 1;

  vec_result.reserve(nsteps);

  // Sweep y from larger (ty_max) to smaller (ty_min).
  for (int i = 0; i < nsteps; ++i) {
    const double tycnt = ty_max - 0.5 * ty_len - static_cast<double>(i) * ty_step;
    const double ty_lower_slice = tycnt - 0.5 * ty_len;
    const double ty_upper_slice = tycnt + 0.5 * ty_len;
    const double tan_y_low = angle_util::to_tangent(ty_lower_slice, angle_unit);
    const double tan_y_up  = angle_util::to_tangent(ty_upper_slice, angle_unit);
    const double elev_obj_from_det_lower = elev_obj_from_det - tan_y_up  * delta_horizontal;
    const double elev_obj_from_det_upper = elev_obj_from_det - tan_y_low * delta_horizontal;

    // Shift to terrain-based depth: add diff_elev (= elev_terrain - elev_obj)
    const double elev_from_det_lower = elev_obj_from_det_lower + diff_elev;
    const double elev_from_det_upper = elev_obj_from_det_upper + diff_elev;

    if( elev_obj_from_det_lower < 0.0 ) continue; // skip: object protrudes above object-top reference

    const SignalNoiseSum sn = mp_get_signal_noise_sum_range_xy(
      tx_lower, tx_upper, ty_lower_slice, ty_upper_slice);

    SignalNoiseDepth snd;
    snd.detid = detid;
    snd.sn = sn;
    snd.elev_from_det_lower = elev_from_det_lower;
    snd.elev_from_det_upper = elev_from_det_upper;

    vec_result.push_back(snd);
  }

  return vec_result;
}

DetectorPanel DetectorPanel::cut(
    const double x_lower, const double x_upper
  , const double y_lower, const double y_upper
  , const double x_eps, const double y_eps) const
{
  // Get a Grid2dBinGroup already cut for axes and signal/noise arrays.
  Grid2dBinGroup g2bg_cut = g2bg_.cut(x_lower, x_upper, y_lower, y_upper, x_eps, y_eps);

  const Grid1d::RangeIndices range_x = get_x_axis().calc_range_indices(x_lower, x_upper, x_eps);
  const Grid1d::RangeIndices range_y = get_y_axis().calc_range_indices(y_lower, y_upper, y_eps);

  const int nbinx_new = g2bg_cut.get_nbinx();
  const int nbiny_new = g2bg_cut.get_nbiny();
  const int expected_nbinx = (range_x.end >= range_x.start) ? (range_x.end - range_x.start + 1) : 0;
  const int expected_nbiny = (range_y.end >= range_y.start) ? (range_y.end - range_y.start + 1) : 0;

  if (expected_nbinx != nbinx_new || expected_nbiny != nbiny_new) {
    THROW_ERROR("DetectorPanel::cut : inconsistent bin counts between source and result axes");
  }

  DetectorPanel result(*this);
  result.set_name(this->get_name() + "_cut");

  g2bg_cut.set_name(g2bg_.get_name() + "_cut");
  g2bg_cut.set_dic_ptr(pdic_);
  result.g2bg_ = g2bg_cut;
  result.g2bg_.set_detid(detid_);
  result.g2bg_.set_dic_ptr(pdic_);

  result.set_n_element(nbinx_new * nbiny_new);

  std::vector<std::vector<DetectorElement>> new_elements(nbiny_new);
  int next_inthis = 0;

  for (int iy_new = 0; iy_new < nbiny_new; ++iy_new) {
    new_elements.at(iy_new).resize(nbinx_new);
    const int iy_org = range_y.start + iy_new;
    for (int ix_new = 0; ix_new < nbinx_new; ++ix_new) {
      const int ix_org = range_x.start + ix_new;

      DetectorElement ele = vec_vec_DetectorElement.at(iy_org).at(ix_org);
      ele.set_id_in_this_detector(next_inthis++);
      ele.set_uqid(UqidNotAssigned);
      new_elements.at(iy_new).at(ix_new) = ele;
    }
  }

  result.vec_vec_DetectorElement = std::move(new_elements);

  return result;
}

void DetectorPanel::copy_signal_noise_to_g2bg(
    const double signal_init, const double noise_init
  , const bool is_avail_init
  , const double PL_thres
  , const double DL_thres
  , const double signal_under_thres
  , const double noise_under_thres
  , const bool is_avail_under_thres )
{
  g2bg_.init_vec_vec(signal_init,noise_init,is_avail_init);

  g2bg_.set_name( "g2bg_" + this->get_name() );
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();

  // Initialize g2bg_.map_igroup_ixiy and map_ixiy_igroup.
  // detector element loop
  #pragma omp parallel for collapse(2) schedule(static)
  for (int iy = 0; iy < nbiny; ++iy) {
    for (int ix = 0; ix < nbinx; ++ix) {
      bool is_avail = is_avail_init;
      
      // get the pointer of DetectorElement
      const DetectorElement& ele = getDetectorElement(ix,iy);
      const double PL = ele.get_PL();
      const double DL = ele.get_DL();
      double signal = ele.get_signal();
      double noise_poi = ele.get_noise_poi();  // Poisson-bucket side
      double noise_det = ele.get_noise_det();  // deterministic floor side

      // if PL is smaller than PL_thres,
      // or DL is smaller than DL_thres,
      if( PL < PL_thres || DL < DL_thres ){
        // set is_avail to under_thres values
        is_avail = is_avail_under_thres;

        // if signal_under_thres is negative, do not change signal
        // if signal_under_thres is positive or zero, set signal to signal_under_thres
        if(signal_under_thres>=0.0 ) signal = signal_under_thres;

        // if noise_under_thres is negative, do not change noise
        // if noise_under_thres is positive or zero, override the total noise:
        // route it to the Poisson bucket (poi) and zero the deterministic floor (det)
        if(noise_under_thres>=0.0 ){
          noise_poi = noise_under_thres;
          noise_det = 0.0;
        }
      }

      // Set signal/noise(poi)/noise_det/is_available to g2bg_ vec_vec buffers.
      g2bg_.set_signal(ix,iy,signal);
      g2bg_.set_noise(ix,iy,noise_poi);
      g2bg_.set_noise_det(ix,iy,noise_det);
      g2bg_.set_is_avail(ix,iy,is_avail);
    }
  }  // detector element loop end
}

int DetectorPanel::assign_naive_igroup( const int igroup_start )
{
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  // detector element loop
  int igroup_last = 0;
  for (int iy = 0; iy < nbiny; ++iy) {
    for (int ix = 0; ix < nbinx; ++ix) {
      DetectorElement& ele = callDetectorElement(ix, iy);
      // Determine igroup uniquely here.
      const Igroup igroup = igroup_start + ix + nbinx * iy;
      g2bg_.insert_map(igroup, {ix, iy});

      // update igroup_last
      igroup_last = igroup;
    }
  }  // detector element loop end
  // Return the last igroup.
  return igroup_last;
}

int DetectorPanel::assign_1st_igroup(
    const int igroup_start
  , const int nx_div_init, const int ny_div_init )
{
  const int nbinx = get_nbinx(); // Minimum bin count in x (equals DetectorElement count).
  const int nbiny = get_nbiny(); // Minimum bin count in y.
  // Compute minimum bins per bin group.
  const int nx_per_group = std::ceil(static_cast<double>(nbinx) / nx_div_init);
  const int ny_per_group = std::ceil(static_cast<double>(nbiny) / ny_div_init);

  // Initialize g2bg_.map_igroup_ixiy and map_ixiy_igroup.
  // detector element loop
  int igroup_last = igroup_start;
  for (int iy = 0; iy < nbiny; ++iy) {
    // Compute y-direction group partition.
    const int iy_div = iy / ny_per_group;

    for (int ix = 0; ix < nbinx; ++ix) {
      // Compute x-direction group partition.
      const int ix_div = ix / nx_per_group;
      
      // Determine igroup uniquely here.
      const Igroup igroup = igroup_start + nx_div_init * iy_div + ix_div;
      
      // Register igroup in the detector-specific bimap.
      g2bg_.insert_map(igroup, {ix, iy});
      igroup_last = igroup;
    }
  }  // detector element loop end

  // Return the last igroup.
  return igroup_last;
}

// Perform bin grouping as defined by the ASCII file.
// Parallelized example of grouping_by_bin_list.
void DetectorPanel::grouping_by_bin_list(const std::filesystem::path &path_in)
{
  // 1) Allocate memory.
  g2bg_.init_vec_vec();

  // 2) Set name based on filename.
  const std::string basename = path_in.stem().string();
  g2bg_.set_name("g2bg_" + basename);

  // 3) Load RectAngularBinGroup and get bin info.
  const RectAngularBinGroup rect_bin_group(path_in);
  rect_bin_group.check_no_overlap();
  rect_bin_group.check_no_void();

  // 4) Prepare a vector to store results temporarily.
  //    Store nbinx * nbiny entries of ((ix, iy), igroup).
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  std::vector<std::pair<Ixiy, Igroup>> vec_tmpResult;
  vec_tmpResult.resize(nbinx * nbiny);

  // 5) Compute (ix, iy) -> igroup in parallel and store in vec_tmpResult.
  //    (Do not insert into the map yet.)
  #pragma omp parallel for collapse(2) schedule(static)
  for (int iy = 0; iy < nbiny; ++iy) {
    for (int ix = 0; ix < nbinx; ++ix) {
      const int idx = iy * nbinx + ix;
      const DetectorElement& ele = getDetectorElement(ix, iy);
      const double tx = ele.get_tx();
      const double ty = ele.get_ty();
      const Igroup igroup = rect_bin_group.get_index(tx, ty);

      vec_tmpResult.at(idx) = std::make_pair(std::array<int, 2>{ix, iy}, igroup);
    }
  }

  // 6) Insert into the map on a single thread.
  //    (Avoids locking overhead.)
  for (const auto &entry : vec_tmpResult) {
    const auto ixiy   = entry.first;
    const int igroup = entry.second;
    g2bg_.insert_map(igroup, ixiy);
  }
}

// save vec_vec_DetectorElement to std::ofstream &ofs
void DetectorPanel::save_vec_vec_DetectorElement( std::ofstream &ofs ) const
{
  const int nbiny = get_nbiny();
  const int nbinx = get_nbinx();
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      const DetectorElement& ele = getDetectorElement(ix,iy);
      ele.save(ofs); // Binary save happens only here.
    }
  }
}

// load vec_vec_DetectorElement from std::ifstream &ifs
void DetectorPanel::load_vec_vec_DetectorElement( std::ifstream &ifs )
{
  // Must call g2bg_.load first.
  vec_vec_memory_allocate();
  const int nbiny = get_nbiny();
  const int nbinx = get_nbinx();

  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      DetectorElement& ele = callDetectorElement(ix,iy);
      ele.load(ifs);
    }
  }
}

// save all variables of DetectorPanel to std::ofstream &ofs
void DetectorPanel::save( std::ofstream &ofs ) const
{
  g2bg_.save(ofs); // save Grid2dBinGroup variables
  io_binary::write_binary(ofs,detid_); // save detid_
  io_binary::write_string(ofs,name); // save name
  ray3d.save(ofs); // save ray3d
  io_binary::write_vec3d(ofs,v3_det_length); // save v3_det_length
  io_binary::write_binary(ofs,n_unit); // save n_unit
  io_binary::write_binary(ofs,days); // save days
  io_binary::write_binary(ofs, static_cast<int>(angle_unit)); // save angle_unit
  io_binary::write_binary(ofs,n_element); // save n_element
  save_vec_vec_DetectorElement(ofs); // save vec_vec_DetectorElement
  prm_bingrp.save(ofs); // save prm_bingrp
}

// load all variables of DetectorPanel from std::ifstream &ifs
void DetectorPanel::load( std::ifstream &ifs )
{
  g2bg_.load(ifs); // load Grid2dBinGroup variables
  detid_ = io_binary::read_binary<int>(ifs); // load detid_
  name = io_binary::read_string(ifs); // load name
  ray3d.load(ifs); // load ray3d
  v3_det_length = io_binary::read_vec3d(ifs); // load v3_det_length
  n_unit = io_binary::read_binary<double>(ifs); // load n_unit
  days = io_binary::read_binary<double>(ifs); // load days
  angle_unit = static_cast<DetectorElement::AngleUnit>(io_binary::read_binary<int>(ifs)); // load angle_unit
  n_element = io_binary::read_binary<int>(ifs); // load n_element
  load_vec_vec_DetectorElement(ifs); // load vec_vec_DetectorElement
  prm_bingrp.load(ifs); // load prm_bingrp
  if(ifs.fail()) THROW_ERROR("ifs.fail() is true");
}

// save all variables of DetectorPanel to fs::path
void DetectorPanel::save( const fs::path& pathout ) const
{
  std::ofstream ofs = io_binary::open_ofstream(pathout);
  // Write architecture info first.
  io_binary::ArchitectureInfo currentInfo = io_binary::get_current_architecture_info();
  io_binary::write_architecture_info(ofs, currentInfo);
  // Then write the data.
  save(ofs); ofs.close();
};

// load all variables of DetectorPanel from fs::path
void DetectorPanel::load( const fs::path &path_in )
{
  std::ifstream ifs = io_binary::open_ifstream(path_in);
  // 1) Read architecture info stored in the file.
  io_binary::ArchitectureInfo fileInfo = io_binary::read_architecture_info(ifs);
  // 2) Check compatibility with current environment (throw if mismatch).
  io_binary::check_architecture_compatibility_or_throw(fileInfo);
  // 3) Read the data.
  load(ifs); ifs.close();
};
