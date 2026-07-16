
// src/cls_FluxTableBuilder.cpp
#include "cls_FluxTableBuilder.hpp"
#include "ns_myapp.hpp"
#include "ns_mylogger.hpp"

namespace {

void read_bool_with_legacy_key(
  const json& section,
  const std::string& current_key,
  const std::string& legacy_key,
  bool& value)
{
  if (section.contains(current_key)) {
    value = section.at(current_key).get<bool>();
    if (section.contains(legacy_key)) {
      if (section.at(legacy_key).is_boolean()
          && section.at(legacy_key).get<bool>() != value) {
        LOG_WARN(
          "FluxTableBuilder: deprecated JSON key '{}' is ignored because '{}' is set to a different value",
          legacy_key, current_key);
      } else {
        LOG_WARN(
          "FluxTableBuilder: deprecated JSON key '{}' is ignored because '{}' is set",
          legacy_key, current_key);
      }
    }
    return;
  }

  if (section.contains(legacy_key)) {
    value = section.at(legacy_key).get<bool>();
    LOG_WARN(
      "FluxTableBuilder: deprecated JSON key '{}' is used; please rename it to '{}'",
      legacy_key, current_key);
  }
}

} // namespace

//----------------------------------------
// Constructor: initialize from JSON
FluxTableBuilder::FluxTableBuilder(
  const json& js, const std::string& section_name )
{
  // Set default instance name
  name = "FluxTableBuilder";

  const auto& section = js.at(section_name);
  std::string key;

  // Load pathin_logR if present
  key = TOSTRING(pathin_logR);
  if (section.contains(key)) {
    pathin_logR = section.at(key).get<std::string>();
  }

  // Load tf_xcnt_logR flag if present
  key = TOSTRING(tf_xcnt_logR);
  if (section.contains(key)) {
    tf_xcnt_logR = section.at(key).get<bool>();
  }

  // Load pathin_log_dFdE if present
  key = TOSTRING(pathin_log_dFdE);
  if (section.contains(key)) {
    pathin_log_dFdE = section.at(key).get<std::string>();
  }

  read_bool_with_legacy_key(
    section, TOSTRING(tf_xcnt_dFdE), "tf_xcnt_log_dFdE", tf_xcnt_dFdE);
  read_bool_with_legacy_key(
    section, TOSTRING(tf_ycnt_dFdE), "tf_ycnt_log_dFdE", tf_ycnt_dFdE);

  // Build interpolation and spectrum grids from input files
  g1_logR_logE = Grid1dXZ(pathin_logR, tf_xcnt_logR);
  g2_log_dFdE  = Grid2dXYZ(pathin_log_dFdE, tf_xcnt_dFdE, tf_ycnt_dFdE);

  // Name the grids for logging purposes
  g1_logR_logE.set_name("logR");
  g2_log_dFdE.set_name("log_dFdE");

  // Load output file paths if present
  key = TOSTRING(pathout_log_penef_R_costhz);
  if (section.contains(key)) {
    pathout_log_penef_R_costhz = section.at(key).get<std::string>();
  }
  key = TOSTRING(pathout_dFdR_R_costhz);
  if (section.contains(key)) {
    pathout_dFdR_R_costhz = section.at(key).get<std::string>();
  }

  // Load integration parameters: dlogE, Rmin, Rmax, dR
  key = TOSTRING(dlogE);
  if (section.contains(key)) {
    dlogE = section.at(key).get<double>();
  }
  key = TOSTRING(Rmin);
  if (section.contains(key)) {
    Rmin = section.at(key).get<double>();
  }
  key = TOSTRING(Rmax);
  if (section.contains(key)) {
    Rmax = section.at(key).get<double>();
  }
  key = TOSTRING(dR);
  if (section.contains(key)) {
    dR = section.at(key).get<double>();
  }

  // Rmin default: when JSON omits Rmin, derive it from the dF/dE table's
  // lowest energy via the range table (R at Emin), then snap it up onto the
  // Rmax/dR lattice so that (Rmax - Rmin) / dR stays an exact integer
  // (Grid1d rejects non-integer bin counts).
  if (!section.contains(TOSTRING(Rmin))) {
    double log10_Emin = g2_log_dFdE.get_ymin();
    const auto arr = g1_logR_logE.get_zmin_x_zmin_zmax_x_zmax();
    if (log10_Emin < arr[1]) {
      LOG_WARN("FluxTableBuilder: dFdE log10(Emin)={:.4f} is below the range-table "
        "domain min {:.4f}; using the domain edge instead.", log10_Emin, arr[1]);
      log10_Emin = arr[1];
    }
    const double Rmin_derived
      = std::pow(10.0, g1_logR_logE.get_linear_interpolated_z(log10_Emin));
    const int nbin_tmp = static_cast<int>(std::floor((Rmax - Rmin_derived) / dR));
    Rmin = Rmax - static_cast<double>(nbin_tmp) * dR;
    LOG_INFO("FluxTableBuilder: Rmin not given in JSON; derived R(Emin)={:.4E} kg/m^2, "
      "snapped up to Rmin={:.4E} on the Rmax/dR lattice.", Rmin_derived, Rmin);
  }

  // Load file-output flags for axis centers
  key = TOSTRING(tf_costhz_cnt);
  if (section.contains(key)) {
    tf_costhz_cnt = section.at(key).get<bool>();
  }
  key = TOSTRING(tf_range_cnt);
  if (section.contains(key)) {
    tf_range_cnt = section.at(key).get<bool>();
  }

  // Load dF/dR divergence-check knobs (all optional)
  key = TOSTRING(tf_check_dFdR_divergence);
  if (section.contains(key)) {
    tf_check_dFdR_divergence = section.at(key).get<bool>();
  }
  key = TOSTRING(dFdR_divergence_threshold);
  if (section.contains(key)) {
    dFdR_divergence_threshold = section.at(key).get<double>();
  }
  key = TOSTRING(tf_dFdR_divergence_fatal);
  if (section.contains(key)) {
    tf_dFdR_divergence_fatal = section.at(key).get<bool>();
  }
}

//----------------------------------------
// Build cutoff-energy table for each range bin
void FluxTableBuilder::build_log10Ecut_table(const Grid1d& g1_range_linear)
{
  // Extract the grid boundaries: zmin, x_zmin, zmax, x_zmax
  const auto arr = g1_logR_logE.get_zmin_x_zmin_zmax_x_zmax();
  const double logRmin = arr[0];
  const double logEmin = arr[1];
  const double logRmax = arr[2];
  const double logEmax = arr[3];

  // Number of range bins
  const int N = g1_range_linear.get_nbin();

  // Initialize cutoff arrays
  vec_log10Ecut.assign(N, -std::numeric_limits<double>::infinity());
  vec_valid_cut.assign(N, false);

  // Sequential loop over each range bin to avoid UB from exceptions in parallel region
  for (int i = 0; i < N; ++i) {
    // Get the lower range value for this bin and compute its log10
    double R = g1_range_linear.get_lower_value(i);
    double log10R = std::log10(R);

    // inverse interpolation to find the corresponding log10Ecut
    auto vecEcut = g1_logR_logE.get_linear_interpolated_x(log10R);
    if (vecEcut.empty()) {
      THROW_ERROR(
        "FluxTableBuilder::build_log10Ecut_table: No Ecut found for log10R={} at bin {}", log10R, i);
    }

    // The first element is the cutoff energy
    double log10Ecut = vecEcut.front();

    // Skip if cutoff outside valid energy bounds
    if (log10Ecut < logEmin || log10Ecut >= logEmax) {
      THROW_ERROR(
        "FluxTableBuilder::build_log10Ecut_table: Cutoff {} out of bounds [{},{}] at bin {}",
        log10Ecut, logEmin, logEmax, i);
    }

    // Store valid cutoff
    vec_log10Ecut[i] = log10Ecut;
    vec_valid_cut[i] = true;
  }

  LOG_INFO("Completed build_log10Ecut_table.");
}

//----------------------------------------
// Build log10-scale penetrating flux table
Grid2dXYZ FluxTableBuilder::build_log_peneflux()
{
  LOG_INFO("Starting build_log_peneflux.");

  // Create x-axis (costhz) and y-axis (range) for the output table
  Grid1d axis_costhz = make_costhz_axis(g2_log_dFdE);  // x-axis: costhz
  Grid1d axis_range  = make_lin_range_axis();          // y-axis: range
  axis_costhz.set_name("costhz");
  axis_range.set_name("range");

  // Initialize output table: z = log10(penetrating flux)
  Grid2dXYZ g2_pene(axis_costhz, axis_range);
  g2_pene.set_name("log10_peneflux");

  // Get number of bins
  int Ncx = g2_pene.get_nbinx();  // number of costhz bins
  int Nry = g2_pene.get_nbiny();  // number of range bins

  // Determine valid energy integration range
  double logE_interval = g2_log_dFdE.get_y_interval();
  double logEmin = g2_log_dFdE.get_ymin();
  double logEmax = g2_log_dFdE.get_ymax() - logE_interval;

  std::atomic<int> counter{0};
  constexpr double EPS = 1e-12;

  // Loop over costhz bins
  #pragma omp parallel for schedule(dynamic)
  for (int ix = 0; ix < Ncx; ++ix) {
    double costhz = g2_pene.get_xcnt(ix);

    // Get dFdE profile as a function of logE for this costhz
    Grid1dXZ g1_log_dFdE_logE = g2_log_dFdE.get_Grid1dXZ_interp(costhz);
    if (g1_log_dFdE_logE.get_nbin() == 0) {
      THROW_ERROR("No data for costhz={}", costhz);
    }

    // Loop over range bins
    for (int iy = 0; iy < Nry; ++iy) {
      // Skip if cutoff is not valid for this range
      if (!vec_valid_cut[iy]) {
        g2_pene.set_z(ix, iy, -std::numeric_limits<double>::infinity());
        continue;
      }

      // Get cutoff energy for this range bin
      double logEcut = vec_log10Ecut[iy];

      // Integration bounds [lo, hi] in log10(E)
      // Ensure logEcut does not fall below the table minimum to avoid invalid integration
      double lo = std::max(logEcut, logEmin);

      // Ensure upper bound does not exceed table maximum to avoid referencing out-of-bounds data
      // Subtract EPS to prevent accessing data outside the dFdE table range
      double hi = std::min(logEmax, g1_log_dFdE_logE.get_max()) - EPS;

      // Skip if integration range is invalid
      if (lo >= hi) {
        g2_pene.set_z(ix, iy, -std::numeric_limits<double>::infinity());
        continue;
      }

      // Number of steps in trapezoidal integration
      int nb = static_cast<int>(std::floor((hi - lo) / dlogE));
      if (nb <= 0) {
        g2_pene.set_z(ix, iy, -std::numeric_limits<double>::infinity());
        continue;
      }

      // Initialize integration
      double integral = 0.0;
      double prev_logE = lo;
      double prev_E = std::pow(10.0, prev_logE);
      double prev_val = std::pow(10.0, g1_log_dFdE_logE.get_linear_interpolated_z(prev_logE));

      // Trapezoidal integration loop, F = \integral dFdE(E) dE (from E_cut to E_max）
      for (int k = 1; k <= nb; ++k) {
        // log10(E) for current step
        double logE = lo + k * dlogE;
        
        // Stop if we exceed upper bound
        if (logE >= hi) break;

        // Convert logE to E in linear scale
        double E = std::pow(10.0, logE);

        // Interpolate log10(dFdE) at logE and convert to linear dFdE
        double val = std::pow(10.0, g1_log_dFdE_logE.get_linear_interpolated_z(logE));

        // Apply trapezoidal rule:
        // Area = 0.5 * (prev + current) * (E - prev_E)        
        integral += 0.5 * (prev_val + val) * (E - prev_E);

        // Update previous point for next iteration
        prev_logE = logE;
        prev_E = E;
        prev_val = val;
      }

      // Store result in log10 scale or mark invalid
      double log_flux = (integral > 0)
        ? std::log10(integral)
        : -std::numeric_limits<double>::infinity();
      g2_pene.set_z(ix, iy, log_flux);
    }

    // Optional progress print to stderr
    ++counter;
    fprintf(stderr, "Progress: %5d/%5d\r", counter.load(), Ncx);
  }

  LOG_INFO("Completed build_log_peneflux.");
  return g2_pene;
}

//----------------------------------------
// Build linear-scale dF/dR from log10 penetrating flux
Grid2dXYZ FluxTableBuilder::build_dFdR_from_peneflux(
  const Grid2dXYZ& g2_log_peneflux_R_costhz )
{
  LOG_INFO("Starting build_dFdR_from_peneflux.");

  // Convert log10 flux to linear scale
  Grid2dXYZ g2_linear = g2_log_peneflux_R_costhz.make_pow10_z();
  g2_linear.set_name("g2_F_R_costhz");

  // Output axis info for debugging
  g2_linear.get_x_axis().out_info(spdlog::level::info);
  g2_linear.get_y_axis().out_info(spdlog::level::info);

  // Compute derivative along the Y (range) axis
  Grid2dXYZ g2_dFdR = g2_linear.make_z_dzdy();
  g2_dFdR.set_name("dFdR_R_costhz");

  // Guard against a divergent (numerical-noise) dF/dR table before returning it.
  if (tf_check_dFdR_divergence) {
    check_dFdR_divergence(
      g2_dFdR, dFdR_divergence_threshold,
      tf_dFdR_divergence_fatal, "FluxTableBuilder::build_dFdR_from_peneflux");
  }

  LOG_INFO("Completed build_dFdR_from_peneflux.");
  return g2_dFdR;
}

//----------------------------------------
// Output built tables to files
void FluxTableBuilder::out_built_tables(
  const Grid2dXYZ& g2_log_penef_R_costhz
, const Grid2dXYZ& g2_dFdR_R_costhz ) const
{
  LOG_INFO("Creating output directories.");

  // Ensure output directories exist
  std::filesystem::create_directories(
    pathout_log_penef_R_costhz.parent_path());
  std::filesystem::create_directories(
    pathout_dFdR_R_costhz.parent_path());

  LOG_INFO("Writing flux tables to files (g2zbin format).");

  // Resolve g2zbin output paths (replace extension)
  auto path_pene = pathout_log_penef_R_costhz;
  path_pene.replace_extension(".g2zbin");
  auto path_dFdR = pathout_dFdR_R_costhz;
  path_dFdR.replace_extension(".g2zbin");

  LOG_INFO("Peneflux g2zbin: {}", path_pene.string());
  LOG_INFO("dFdR g2zbin: {}", path_dFdR.string());

  // Write g2zbin binary files
  g2_log_penef_R_costhz.save_g2zbin(path_pene);
  g2_dFdR_R_costhz.save_g2zbin(path_dFdR);

  // Prepare ASCII output paths
  auto ascii_pene = pathout_log_penef_R_costhz;
  ascii_pene.replace_extension(".txt");
  auto ascii_dFdR = pathout_dFdR_R_costhz;
  ascii_dFdR.replace_extension(".txt");

  LOG_INFO("Peneflux ASCII: {}", ascii_pene.string());
  LOG_INFO("dFdR ASCII: {}", ascii_dFdR.string());

  // Write ASCII tables
  g2_log_penef_R_costhz.out_xy(
    ascii_pene, tf_costhz_cnt, tf_range_cnt);
  g2_dFdR_R_costhz.out_xy(
    ascii_dFdR, tf_costhz_cnt, tf_range_cnt);

  LOG_INFO("All tables written.");
}

//----------------------------------------


Grid1d FluxTableBuilder::make_costhz_axis(const Grid2dXYZ& g2_logdFdE)
{
  Grid1d g1x = g2_logdFdE.get_x_axis();
  g1x.set_name("costhz");
  return g1x;
}

Grid1d FluxTableBuilder::make_lin_range_axis()
{
  const int nbin = static_cast<int>((Rmax - Rmin) / dR);
  return Grid1d("range", nbin, Rmin, Rmax, dR);
}
