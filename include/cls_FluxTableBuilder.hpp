/// @file cls_FluxTableBuilder.hpp
/// @brief Builder class for constructing flux tables from energy and range data
/// @details
/// This class constructs integrated tables of penetrating flux and dF/dR from differential flux data.
///
/// **Workflow:**
/// 1. Load input data (range-energy interpolation table and differential flux spectrum)
/// 2. Build cutoff-energy table for each range bin (build_log10Ecut_table)
/// 3. Build log10-scale penetrating flux table (build_log_peneflux)
/// 4. Build linear-scale dF/dR table (build_dFdR_from_peneflux)
/// 5. Output all tables to binary and ASCII files (out_built_tables)
///
/// **Units:**
/// - Energy: MeV (logarithmic scale: log10(E/MeV))
/// - Range: kg/m²
/// - Flux: particles/(cm²·s·sr) (log10 scale for penetrating flux)
/// - dF/dR: particles/(cm²·s·sr)/(kg/m²)
///
/// **Thread-safety:** Not thread-safe. Methods that use OpenMP are marked explicitly.
///
/// **Typical usage:**
/// @code
/// FluxTableBuilder builder(json_config, "flux_section");
/// Grid1d range_axis = builder.make_lin_range_axis();
/// builder.build_log10Ecut_table(range_axis);
/// Grid2dXYZ peneflux = builder.build_log_peneflux();
/// Grid2dXYZ dFdR = builder.build_dFdR_from_peneflux(peneflux);
/// builder.out_built_tables(peneflux, dFdR);
/// @endcode
#pragma once

#include <string>
#include <optional>
#include <filesystem>
#include "cls_Grid1dXZ.hpp"
#include "cls_Grid2dXYZ.hpp"
#include "cls_Grid1d.hpp"
#include "cls_FluxTable.hpp"
#include "ns_mylogger.hpp"
namespace fs = std::filesystem;

/// @class FluxTableBuilder
/// @brief Builds integrated tables of peneflux and dF/dR from dFdE and range data
/// @ingroup detectorClasses
/// @details
/// **Responsibility:**
/// This class integrates differential energy flux (dF/dE) over energy to produce:
/// - Penetrating flux as a function of range and angle (costhz)
/// - Rate of flux change with respect to range (dF/dR)
///
/// **Invariants:**
/// - g1_logR_logE must contain valid interpolation data (loaded from file or set manually)
/// - g2_log_dFdE must contain valid spectrum data (loaded from file or set manually)
/// - dlogE > 0 (integration step size in log10(E) scale)
/// - Rmin < Rmax and dR > 0 (range axis parameters)
/// - vec_log10Ecut and vec_valid_cut have the same size after build_log10Ecut_table
///
/// **Key use cases:**
/// 1. Atmospheric cosmic ray flux computation at various depths
/// 2. Radiation shielding analysis with energy-dependent cutoffs
/// 3. Pre-computation of lookup tables for detector simulations
///
/// **Minimal usage example:**
/// @code
/// // Load from JSON configuration
/// FluxTableBuilder builder(json_obj, "flux_builder");
///
/// // Build range axis
/// Grid1d range_axis = builder.make_lin_range_axis();
///
/// // Execute build pipeline
/// builder.build_log10Ecut_table(range_axis);
/// auto peneflux = builder.build_log_peneflux();
/// auto dFdR = builder.build_dFdR_from_peneflux(peneflux);
///
/// // Write output files
/// builder.out_built_tables(peneflux, dFdR);
/// @endcode
class FluxTableBuilder {
private:
  /// @brief Instance name (optional)
  std::string name;

  /// @brief Path to input log10(Range) vs log10(Ecut) table
  fs::path pathin_logR{"none"};
  /// @brief Use center values when loading pathin_logR (true = use x_center_value, false= use x_lower_value)
  bool tf_xcnt_logR{false};

  /// @brief Path to input log10(dF/dE) table
  fs::path pathin_log_dFdE{"none"};

  /// @brief Use x-axis center when loading pathin_log_dFdE, false= use x_lower_value, x=costhz
  bool tf_xcnt_dFdE{true};

  /// @brief Use y-axis center when loading pathin_log_dFdE, false= use y_lower_value, y=log10E
  bool tf_ycnt_dFdE{true};

  /// @brief Interpolation data for log10(Range) → log10(Ecut) (y = log10R, x = log10E)
  Grid1dXZ g1_logR_logE;

  /// @brief Energy spectrum data (x = costhz, y = log10E, z = log10(dF/dE))
  Grid2dXYZ g2_log_dFdE;

  /// @brief Cutoff energies (log10Ecut) for each range bin
  std::vector<double> vec_log10Ecut;

  /// @brief Validity flags for vec_log10Ecut
  std::vector<bool>   vec_valid_cut;

  /// @brief Output file path for peneflux binary (log scale)
  fs::path pathout_log_penef_R_costhz{"none"};
  
  /// @brief Output file path for dF/dR binary
  fs::path pathout_dFdR_R_costhz{"none"};

  /// @brief Energy integration step size (log10 scale)
  double dlogE;

  /// @brief Minimum range for integration (kg/m²)
  double Rmin = 5000.0;

  /// @brief Maximum range for integration (kg/m²)
  double Rmax = 4999000.0;

  /// @brief Range step size for integration (kg/m²)
  double dR   = 1000.0;

  /// @brief Output costhz axis as center values when writing files
  bool tf_costhz_cnt = true;

  /// @brief Output range axis as center values when writing files
  bool tf_range_cnt  = false;

  /// @brief Enable the dF/dR divergence (numerical-noise) check after building
  bool tf_check_dFdR_divergence{true};
  /// @brief Sign-flip ratio threshold per costhz slice for the divergence check
  /// @note Default 0.05; see @ref check_dFdR_divergence() for the threshold calibration.
  double dFdR_divergence_threshold{0.05};
  /// @brief If true, abort on a divergent dF/dR table; if false, only warn
  bool tf_dFdR_divergence_fatal{true};

public:
  //============================================================================
  /// @name Constructors
  ///@{

  /// @brief Default constructor
  FluxTableBuilder() = default;

  /// @brief Copy constructor
  FluxTableBuilder(const FluxTableBuilder &other) = default;

  /// @brief Move constructor
  FluxTableBuilder(FluxTableBuilder &&other) noexcept = default;

  /// @brief Destructor
  ~FluxTableBuilder() = default;

  /// @brief Initialize parameters from JSON
  /// @param[in] js JSON object containing configuration
  /// @param[in] section_name Name of the section in js
  /// @throws std::exception If required JSON keys are missing or invalid
  /// @note Expected JSON keys: pathin_logR, tf_xcnt_logR, pathin_log_dFdE, tf_xcnt_dFdE,
  ///       tf_ycnt_dFdE, pathout_log_penef_R_costhz, pathout_dFdR_R_costhz,
  ///       dlogE, Rmin, Rmax, dR, tf_costhz_cnt, tf_range_cnt,
  ///       tf_check_dFdR_divergence, dFdR_divergence_threshold,
  ///       tf_dFdR_divergence_fatal (all optional with defaults)
  ///       Deprecated aliases tf_xcnt_log_dFdE and tf_ycnt_log_dFdE are still accepted.
  FluxTableBuilder(const nlohmann::json& js, const std::string& section_name);

  ///@}
  //============================================================================
  /// @name getter functions
  ///@{

  /// @brief Get instance name
  std::string get_name() const { return name; }

  /// @brief Get reference to log10(Range) vs log10(Ecut) data
  const Grid1dXZ& get_g1_logR() const { return g1_logR_logE; }

  /// @brief Get reference to log10(dF/dE) data
  const Grid2dXYZ& get_g2_log_dFdE() const { return g2_log_dFdE; }

  /// @brief Get energy integration step size
  double get_dlogE() const { return dlogE; }

  /// @brief Get minimum range
  double get_Rmin()   const { return Rmin; }

  /// @brief Get maximum range
  double get_Rmax()   const { return Rmax; }

  /// @brief Get range step size
  double get_dR()     const { return dR; }

  ///@}
  //============================================================================
  /// @name Setters
  ///@{

  /// @brief Set instance name
  void set_name(const std::string& name_in) { name = name_in; }

  /// @brief Set log10(Range) vs log10(Ecut) data
  void set_g1_logR(const Grid1dXZ& g1_in)   { g1_logR_logE = g1_in; }

  /// @brief Set log10(dF/dE) data
  void set_g2_log_dFdE(const Grid2dXYZ& g2_in) { g2_log_dFdE = g2_in; }

  /// @brief Set energy integration step size (log10 scale)
  void set_dlogE(double d) { dlogE = d; }

  /// @brief Set minimum range
  void set_Rmin(double v)  { Rmin = v; }

  /// @brief Set maximum range
  void set_Rmax(double v)  { Rmax = v; }

  /// @brief Set range step size
  void set_dR(double v)    { dR   = v; }

  ///@}
  //============================================================================
  /// @name Core Build Methods
  ///@{

  /// @brief Build the cutoff-energy table for each range bin
  /// @param[in] g1_range_linear Linear range axis (Grid1d) in units of kg/m²
  /// @throws std::runtime_error If interpolation fails or cutoff is out of energy bounds
  /// @note This method populates vec_log10Ecut and vec_valid_cut using inverse interpolation
  ///       from g1_logR_logE. Cutoff energies outside [logEmin, logEmax) are rejected.
  ///       Thread-safety: Sequential execution (removed OpenMP to allow safe exception throwing)
  void build_log10Ecut_table(const Grid1d& g1_range_linear);

  /// @brief Build log10-scale penetrating flux table
  /// @details For each costhz, extracts the dFdE profile and integrates
  ///          above the cutoff energy using the trapezoidal rule.
  ///          The integration computes: F = ∫[E_cut to E_max] (dF/dE) dE
  /// @return Grid2dXYZ of log10(penetrating flux) with x=costhz, y=range, z=log10(flux)
  /// @note Uses OpenMP parallel for loop over costhz bins (schedule=dynamic).
  ///       Invalid bins (no data or invalid range) are set to -infinity.
  ///       Integration step size is controlled by dlogE (log10 scale).
  ///       Thread-safety: Safe for concurrent reads of g2_log_dFdE and vec_log10Ecut.
  ///       Units: flux in particles/(cm²·s·sr)
  /// @note Complexity: O(Ncosthz * Nrange * Nsteps) where Nsteps ≈ (logEmax - logEcut)/dlogE
  Grid2dXYZ build_log_peneflux();

  /// @brief Build linear-scale dF/dR table from penetrating flux
  /// @param[in] g2_log_peneflux_R_costhz Input log10 penetrating flux table
  /// @return Grid2dXYZ of linear-scale dF/dR with x=costhz, y=range, z=dF/dR
  /// @note This method first converts log10(flux) to linear flux (10^z),
  ///       then computes the derivative ∂F/∂R along the range (y) axis.
  ///       Units: dF/dR in particles/(cm²·s·sr)/(kg/m²)
  /// @note Complexity: O(Ncosthz * Nrange)
  Grid2dXYZ build_dFdR_from_peneflux(const Grid2dXYZ& g2_log_peneflux_R_costhz);

  /// @brief Output built tables in binary and ASCII formats
  /// @param[in] g2_log_penef_R_costhz Log10-scale penetrating flux
  /// @param[in] g2_dFdR_R_costhz      Linear-scale dF/dR
  /// @note Creates output directories if they don't exist.
  ///       Writes binary files to pathout_log_penef_R_costhz and pathout_dFdR_R_costhz.
  ///       Writes ASCII files with .txt extension.
  ///       Uses tf_costhz_cnt and tf_range_cnt flags to determine axis output format.
  void out_built_tables(
    const Grid2dXYZ& g2_log_penef_R_costhz
  , const Grid2dXYZ& g2_dFdR_R_costhz ) const;

  ///@}
  //============================================================================
  /// @name Utility Methods (static)
  ///@{

  /// @brief Create costhz axis for peneflux
  /// @param[in] g2_logdFdE Input dFdE table
  /// @return Grid1d representing costhz axis (extracted from g2_logdFdE's x-axis)
  /// @note This is a static utility that extracts and copies the x-axis from the input grid.
  static Grid1d make_costhz_axis(const Grid2dXYZ& g2_logdFdE);

  /// @brief Create linear range axis from Rmin to Rmax with step dR
  /// @return Grid1d representing range axis in kg/m²
  /// @note Number of bins computed as: nbin = (Rmax - Rmin) / dR
  ///       Uses Rmin, Rmax, and dR member variables.
  Grid1d make_lin_range_axis();

  ///@} ----------------------------------------------------------------
};
