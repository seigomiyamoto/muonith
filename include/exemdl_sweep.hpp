/// @file exemdl_sweep.hpp
/// @brief Parameter sweep orchestration for 3D density reconstruction.
/// @details
/// Provides `exemdl::sweep` namespace for sweeping sigma_rho, corr_length,
/// sigma_rho_diag, and uniform_prior_density parameters.
///
/// Sweep nesting (outer = expensive, inner = cheap):
/// - Outer loop: uniform_prior_density (re-runs compute_prior + build_observation_matrix + invert_density)
/// - Inner loop: sigma_rho x corr_length x sigma_rho_diag (re-runs invert_density only)
///
/// @note Thread safety: Not thread-safe (uses global iodir output directory).
/// @note Units: densities in kg/m^3, correlation length in meters.
#pragma once

#include <array>
#include <string>
#include <vector>
#include <filesystem>
#include <nlohmann/json.hpp>

namespace fs = std::filesystem;

namespace exemdl::sweep {

  /// @brief A single parameter combination for one sweep iteration.
  struct SweepPoint {
    double sigma_rho;              ///< Density prior std dev (kg/m^3).
    double corr_length;            ///< Spatial correlation length (meters).
    double sigma_rho_diag;         ///< Diagonal covariance std dev (kg/m^3).
    /// Density parameters: [prior_density, shell_upper, shell_lower, shell_lateral] (kg/m^3).
    std::array<double, 4> density_quad;
    int index;                     ///< 0-based sweep index.
    std::string label;             ///< Human-readable label, e.g., "sr300_cl70_srd300_upd2000".

    /// @brief Get uniform prior density (kg/m^3).
    double prior_density() const { return density_quad[0]; }
    /// @brief Get upper shell density (kg/m^3).
    double shell_density_upper() const { return density_quad[1]; }
    /// @brief Get lower shell density (kg/m^3).
    double shell_density_lower() const { return density_quad[2]; }
    /// @brief Get lateral shell density (kg/m^3).
    double shell_density_lateral() const { return density_quad[3]; }
  };

  /// @brief Configuration for parameter sweep, parsed from JSON5.
  struct SweepConfig {
    bool tf_exec = false;                        ///< Enable sweep execution.
    int base_index = 0;                          ///< Index in NAGAINV_PARAMETERS to use as base.
    std::vector<double> vec_sigma_rho;           ///< Values to sweep (kg/m^3).
    std::vector<double> vec_corr_length;         ///< Values to sweep (meters).
    /// Values to sweep: each element is [prior_density, shell_upper, shell_lower, shell_lateral] (kg/m^3).
    /// Empty = use base.
    std::vector<std::array<double, 4>> vec_density_quad;
    bool link_diag_to_sigma = true;              ///< If true, sigma_rho_diag = sigma_rho.
    std::vector<double> vec_sigma_rho_diag;      ///< Values to sweep when link_diag_to_sigma=false.
    std::string module8_mode = "none";           ///< "none"|"all"|"last"|"selected".
    std::vector<int> module8_indices;            ///< Indices for "selected" mode.

    /// @brief Parse sweep configuration from JSON root.
    /// @param[in] root JSON root object.
    /// @param[in] section Section name (default: "NAGAINV_PARAM_SWEEP").
    /// @return Parsed configuration. If section is absent, tf_exec=false.
    static SweepConfig from_json(const nlohmann::json& root,
        const std::string& section = "NAGAINV_PARAM_SWEEP");

    /// @brief Generate all sweep points as Cartesian product.
    /// @return Ordered vector of sweep points, grouped by uniform_prior_density.
    std::vector<SweepPoint> generate_points() const;

    /// @brief Check if analyze_errors should run for a given sweep index.
    /// @param[in] sweep_index 0-based index of current sweep point.
    /// @param[in] total_count Total number of sweep points.
    /// @return True if analyze_errors should execute.
    bool should_run_module8(int sweep_index, int total_count) const;

    /// @brief Total number of parameter combinations.
    int total_combinations() const;
  };

  /// @brief Result of a single sweep iteration.
  struct SweepResult {
    SweepPoint point;              ///< Parameter values used.
    bool module8_executed = false;  ///< Whether analyze_errors ran for this point.
  };

  /// @brief Execute the full parameter sweep.
  /// @param[in] json_path Path to JSON5 configuration file.
  /// @param[in] seed_given Whether seed is provided via CLI.
  /// @param[in] seed_value Seed value.
  /// @param[in] resume_from Path to checkpoint directory to resume from.
  ///            Empty path means start from scratch (the init, load_parameters, and build_geometry stages).
  /// @param[in] end_stage Pipeline stop stage. Valid: {3,4,5,6,7,8}.
  ///            CLI(--end-stage) > JSON5(end_stage) > default(8).
  ///            Negative value means "resolve from JSON5 or default".
  /// @return Vector of sweep results.
  /// @throws std::runtime_error on configuration or execution errors.
  /// @note Modifies global iodir output directory during execution, restores on completion.
  std::vector<SweepResult> run_sweep(
      const fs::path& json_path,
      bool seed_given = false,
      unsigned seed_value = 0,
      const fs::path& resume_from = "",
      int end_stage = 8);

  /// @brief Resolve an end_stage token to its integer stage number.
  /// @details Accepts either a plain integer string (e.g. "7") or a feature-name
  ///          token: build_geometry(3), trace_path_lengths(4), compute_prior(5),
  ///          build_observation_matrix(6), invert_density(7), analyze_errors(8).
  /// @param[in] token Integer string or feature-name token.
  /// @return Stage number. Range validity ({3..8}) is checked later by run_sweep.
  /// @throws std::runtime_error If token is neither a plain integer nor a known feature name.
  int resolve_end_stage(const std::string& token);

  /// @brief Write CSV summary of sweep results.
  /// @param[in] results Sweep results to summarize.
  /// @param[in] output_path Output CSV file path.
  void write_sweep_summary_csv(
      const std::vector<SweepResult>& results,
      const fs::path& output_path);

} // namespace exemdl::sweep
