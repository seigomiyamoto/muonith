// src/cls_FluxTable.cpp

#include "cls_FluxTable.hpp"
#include "ns_myapp.hpp"
#include "ns_param_constants.hpp"

void FluxTable::assign_parameters(
  const nlohmann::json& js, const std::string& section_name)
{
  if (!js.contains(section_name)) {
    LOG_ERROR("Section '{}' not found in JSON config.", section_name);
    THROW_ERROR_NAME("JSON section not found: " + section_name);
  }

  const auto& sec = js.at(section_name);
  std::string key;

  // path_log_peneflux
  key = TOSTRING(pathin_log_peneflux);
  if (sec.contains(key))
    pathin_log_peneflux = sec.at(key).get<std::string>();

  // tf_xcnt_peneflux
  key = TOSTRING(tf_xcnt_peneflux);
  if (sec.contains(key))
    tf_xcnt_peneflux = sec.at(key).get<bool>();

  // tf_ycnt_peneflux
  key = TOSTRING(tf_ycnt_peneflux);
  if (sec.contains(key))
    tf_ycnt_peneflux = sec.at(key).get<bool>();

  // pathin_dFdR_R_costhz
  key = TOSTRING(pathin_dFdR_R_costhz);
  if (sec.contains(key))
    pathin_dFdR_R_costhz = sec.at(key).get<std::string>();

  // tf_xcnt_dFdR
  key = TOSTRING(tf_xcnt_dFdR);
  if (sec.contains(key))
    tf_xcnt_dFdR = sec.at(key).get<bool>();

  // tf_ycnt_dFdR
  key = TOSTRING(tf_ycnt_dFdR);
  if (sec.contains(key))
    tf_ycnt_dFdR = sec.at(key).get<bool>();

  // tf_check_dFdR_divergence
  key = TOSTRING(tf_check_dFdR_divergence);
  if (sec.contains(key))
    tf_check_dFdR_divergence = sec.at(key).get<bool>();

  // dFdR_divergence_threshold
  key = TOSTRING(dFdR_divergence_threshold);
  if (sec.contains(key))
    dFdR_divergence_threshold = sec.at(key).get<double>();

  // tf_dFdR_divergence_fatal
  key = TOSTRING(tf_dFdR_divergence_fatal);
  if (sec.contains(key))
    tf_dFdR_divergence_fatal = sec.at(key).get<bool>();
}

void FluxTable::load_tables()
{
  if (pathin_log_peneflux != "none") {
    g2_log_penef_R_costhz = Grid2dXYZ(pathin_log_peneflux, tf_xcnt_peneflux, tf_ycnt_peneflux);
    g2_log_penef_R_costhz->set_name("log_peneflux");
    // DL clamp bounds for peneflux reads: explicit PARAM_CONSTANTS values win;
    // otherwise follow this table's own y-axis (table-derived defaults).
    const double DL_lo = param_constants::DL_min_is_set()
      ? param_constants::DL_min() : g2_log_penef_R_costhz->get_ymin();
    const double DL_hi = param_constants::DL_max_is_set()
      ? param_constants::DL_max() : g2_log_penef_R_costhz->get_ymax();
    g2_log_penef_R_costhz->set_DL_clamp_bounds(DL_lo, DL_hi);
    LOG_INFO("FluxTable::load_tables: DL clamp bounds [{:.4E}, {:.4E}] kg/m^2 (min: {}, max: {})",
      DL_lo, DL_hi,
      param_constants::DL_min_is_set() ? "json" : "table",
      param_constants::DL_max_is_set() ? "json" : "table");
  } else {
    LOG_WARN("pathin_log_peneflux is 'none'. Skipping load.");
  }

  if (pathin_dFdR_R_costhz != "none") {
    g2_dFdR_R_costhz = Grid2dXYZ(pathin_dFdR_R_costhz, tf_xcnt_dFdR, tf_ycnt_dFdR);
    g2_dFdR_R_costhz->set_name("dFdR");
    // Guard against loading a divergent (numerical-noise) dF/dR table.
    if (tf_check_dFdR_divergence && g2_dFdR_R_costhz.has_value()) {
      check_dFdR_divergence(
        *g2_dFdR_R_costhz, dFdR_divergence_threshold,
        tf_dFdR_divergence_fatal, "FluxTable::load_tables");
    }
  } else {
    LOG_WARN("pathin_dFdR_R_costhz is 'none'. Skipping load.");
  }

  LOG_INFO("Loaded Grid2dXYZ tables.");
}

//----------------------------------------
// Divergence (numerical-noise) guard shared by the generation and load paths.
// A healthy dF/dR table is finite and strictly negative for every (costhz, R)
// bin because F(R) decreases monotonically with range R. Positive or
// sign-flipping values are the numerical noise that once slipped through to
// the 3D reconstruction (issue 90gy73). Non-finite values are always fatal;
// the divergence signal aborts when fatal == true and only warns otherwise.
void check_dFdR_divergence(
  const Grid2dXYZ& g2_dFdR_R_costhz
, const double threshold
, const bool fatal
, const std::string& context )
{
  // vv[iy][ix]: outer iy = range R (nbiny), inner ix = costhz (nbinx)
  const std::vector<std::vector<double>> vv = g2_dFdR_R_costhz.get_vec_vec_z();
  const std::size_t nbiny = vv.size();
  if (nbiny == 0 || vv.front().empty()) {
    LOG_WARN("{}: dF/dR table is empty; skipping divergence check.", context);
    return;
  }
  const std::size_t nbinx = vv.front().size();

  // Small tolerance for boundary rounding; a healthy table has zero positives.
  const double pos_tol = 0.01;

  std::size_t total = 0;          // points visited (finite + non-finite)
  std::size_t pos = 0;            // positive dF/dR points (invariant: expected 0)
  std::size_t nonfinite = 0;      // NaN / Inf points
  double worst_flip_ratio = 0.0;  // worst per-costhz sign-flip ratio
  std::size_t bad_slices = 0;     // costhz slices above the threshold

  // Aggregate per costhz slice: fix ix, scan range R = iy.
  for (std::size_t ix = 0; ix < nbinx; ++ix) {
    std::size_t flips = 0;  // sign flips along R in this slice
    std::size_t pts = 0;    // points in this slice
    int prev_sign = 0;      // previous non-zero sign
    for (std::size_t iy = 0; iy < nbiny; ++iy) {
      if (ix >= vv[iy].size()) continue;  // guard against a ragged row
      const double z = vv[iy][ix];
      ++total;
      ++pts;
      if (!std::isfinite(z)) { ++nonfinite; continue; }
      if (z > 0.0) ++pos;
      const int sign = (z > 0.0) ? 1 : ((z < 0.0) ? -1 : 0);
      if (sign != 0 && prev_sign != 0 && sign != prev_sign) ++flips;
      if (sign != 0) prev_sign = sign;
    }
    if (pts > 1) {
      const double ratio =
        static_cast<double>(flips) / static_cast<double>(pts - 1);
      if (ratio > worst_flip_ratio) worst_flip_ratio = ratio;
      if (ratio > threshold) ++bad_slices;
    }
  }

  const double pos_ratio =
    (total > 0) ? static_cast<double>(pos) / static_cast<double>(total) : 0.0;

  // Non-finite values are never acceptable, regardless of the fatal flag.
  if (nonfinite > 0) {
    THROW_ERROR(
      "{}: dF/dR table has {} non-finite value(s) out of {} points (NaN/Inf). Table is invalid.",
      context, nonfinite, total);
  }

  // Divergence signal: positive values or many sign flips per costhz slice.
  const bool divergent = (pos_ratio > pos_tol) || (worst_flip_ratio > threshold);
  if (divergent) {
    if (fatal) {
      THROW_ERROR(
        "{}: dF/dR table looks divergent. pos_ratio={:.3e} (expected ~0, all-negative), "
        "worst_sign_flip_ratio={:.3e}, threshold={:.3e}, bad_slices={}/{}. "
        "Known cause: dF/dE axis flags not matching the grid (e.g. tf_ycnt_dFdE=false "
        "for a bin-centered log10KE axis) — set tf_xcnt_dFdE/tf_ycnt_dFdE in the "
        "FLUX_BUILDER config to match and rebuild the table "
        "(docs/user-guide/flux-tables/build-and-install.md).",
        context, pos_ratio, worst_flip_ratio, threshold, bad_slices, nbinx);
    }
    LOG_WARN(
      "{}: dF/dR table looks divergent but continuing (fatal flag off). "
      "pos_ratio={:.3e}, worst_sign_flip_ratio={:.3e}, threshold={:.3e}, bad_slices={}/{}. "
      "Known cause: dF/dE axis flags not matching the grid (tf_xcnt_dFdE/tf_ycnt_dFdE); "
      "see docs/user-guide/flux-tables/build-and-install.md.",
      context, pos_ratio, worst_flip_ratio, threshold, bad_slices, nbinx);
    return;
  }

  LOG_INFO(
    "{}: dF/dR divergence check passed. pos_ratio={:.3e}, worst_sign_flip_ratio={:.3e}, slices={}",
    context, pos_ratio, worst_flip_ratio, nbinx);
}

bool FluxTable::operator!=(const FluxTable& other) const
{
#ifdef NODEBUG
  if (g2_log_penef_R_costhz != other.g2_log_penef_R_costhz) return true;
  if (g2_dFdR_R_costhz != other.g2_dFdR_R_costhz) return true;
  if (pathin_log_peneflux != other.pathin_log_peneflux) return true;
  if (tf_xcnt_peneflux != other.tf_xcnt_peneflux) return true;
  if (tf_ycnt_peneflux != other.tf_ycnt_peneflux) return true;
  if (pathin_dFdR_R_costhz != other.pathin_dFdR_R_costhz) return true;
  if (tf_xcnt_dFdR != other.tf_xcnt_dFdR) return true;
  if (tf_ycnt_dFdR != other.tf_ycnt_dFdR) return true;
#else
  if (g2_log_penef_R_costhz != other.g2_log_penef_R_costhz) {
    LOG_WARN("FluxTable: g2_log_penef_R_costhz differs");
    return true;
  }
  if (g2_dFdR_R_costhz != other.g2_dFdR_R_costhz) {
    LOG_WARN("FluxTable: g2_dFdR_R_costhz differs");
    return true;
  }
  if (pathin_log_peneflux != other.pathin_log_peneflux) {
    LOG_WARN("FluxTable: pathin_log_peneflux differs");
    return true;
  }
  if (tf_xcnt_peneflux != other.tf_xcnt_peneflux) {
    LOG_WARN("FluxTable: tf_xcnt_peneflux differs");
    return true;
  }
  if (tf_ycnt_peneflux != other.tf_ycnt_peneflux) {
    LOG_WARN("FluxTable: tf_ycnt_peneflux differs");
    return true;
  }
  if (pathin_dFdR_R_costhz != other.pathin_dFdR_R_costhz) {
    LOG_WARN("FluxTable: pathin_dFdR_R_costhz differs");
    return true;
  }
  if (tf_xcnt_dFdR != other.tf_xcnt_dFdR) {
    LOG_WARN("FluxTable: tf_xcnt_dFdR differs");
    return true;
  }
  if (tf_ycnt_dFdR != other.tf_ycnt_dFdR) {
    LOG_WARN("FluxTable: tf_ycnt_dFdR differs");
    return true;
  }
#endif
  return false;
}

const Grid2dXYZ& FluxTable::get_g2_log_peneflux() const
{
  if (!g2_log_penef_R_costhz.has_value()) {
    LOG_ERROR("g2_log_penef_R_costhz is not loaded");
    LOG_ERROR("please do FluxTable::load_tables() before using this function");
    THROW_ERROR("peneflux table is not loaded");
  }
  return *g2_log_penef_R_costhz;
};

const Grid2dXYZ& FluxTable::get_g2_dFdR_R_costhz() const
{
  if (!g2_dFdR_R_costhz.has_value()) {
    LOG_ERROR("g2_dFdR_R_costhz is not loaded");
    LOG_ERROR("please do FluxTable::load_tables() before using this function");
    THROW_ERROR("dFdR table is not loaded");
  }
  return *g2_dFdR_R_costhz;
};



double FluxTable::get_peneflux(const double costhz, const double DL) const {
  if (!g2_log_penef_R_costhz.has_value()) {
    LOG_ERROR("g2_log_penef_R_costhz is not loaded");
    LOG_ERROR("please insert FluxTable::load_tables() before using this function");
    THROW_ERROR("log_peneflux table is not loaded");
  }
  return g2_log_penef_R_costhz->get_bilinear_interpolated_peneflux(costhz, DL);
}


double FluxTable::calc_dFdR( const DetectorElement& ele, const double DL_prior ) const
{
  if (!g2_dFdR_R_costhz.has_value()) {
    LOG_ERROR("g2_dFdR_R_costhz is not available");
    THROW_ERROR("dFdR table is not available in this FluxTable instance");
  }

  // Extract cos(theta_z) from detector element's direction
  const double costhz = ele.get_ray3d().vz();

  // Get Range bounds from g2_dFdR_R_costhz table
  constexpr double eps = 1.0e-8;
  const double Rmin = this->g2_dFdR_R_costhz->get_ymin() * ( 1.0 + eps );
  const double Rmax = this->g2_dFdR_R_costhz->get_ymax() * ( 1.0 - eps );

  // Clamp DL to valid Range bounds
  double DL = std::clamp(DL_prior, Rmin, Rmax);

  // Interpolate dF/dR from cos(theta_z) and DL using bilinear interpolation
  const double dFdR = g2_dFdR_R_costhz->get_bilinear_interpolated_z_value(costhz,DL);

  return dFdR;
}
