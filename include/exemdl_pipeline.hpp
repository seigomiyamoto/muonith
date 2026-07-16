/// @file exemdl_pipeline.hpp
/// @brief Pipeline state management for muon tomography modules.
/// @details
/// Provides `exemdl::pipeline` namespace with State structure and functions
/// for checkpoint save/load, enabling:
/// - Module execution from intermediate checkpoints
/// - Iterative estimation (reconstruct -> prior -> reconstruct)
/// - Nested parameter sweeps (uniform_prior_density -> sigma_rho x corr_length)
///
/// @note Coordinate system: z-up, right-handed. All densities in kg/m^3.
/// @note Thread safety: State is not thread-safe for concurrent writes.
#pragma once

#include <array>
#include <filesystem>
#include <optional>
#include <fstream>
#include <nlohmann/json.hpp>
#include <Eigen/Sparse>

#include "exemdl_load_parameters.hpp"
#include "exemdl_init_json_logger_seed.hpp"
#include "exemdl_build_geometry.hpp"
#include "exemdl_build_detector.hpp"
#include "exemdl_run_inversion.hpp"
#include "cls_Grid2dPillar.hpp"
#include "cls_Grid3dVoxel.hpp"
#include "cls_DetectorPanelArray.hpp"

namespace fs = std::filesystem;

/// @namespace exemdl::pipeline
/// @brief Pipeline state management for modular muon tomography execution.
namespace exemdl::pipeline {

  /// @brief Sparse matrix type alias (float precision).
  using SpMatf = Eigen::SparseMatrix<float>;

  /// @brief Pipeline state version for binary compatibility check.
  /// @note v5: vec_signal_group_poisson / vec_noise_poi_group_poisson added to
  ///           Grid2dBinGroup::save/load (older checkpoints are rejected).
  /// @note v6: Grid1d::save/load became pass-through (canonical min/max stored
  ///           as-is; tf_half_shift removed).  v5 checkpoints stored pre-shift
  ///           coordinates and are rejected.
  /// @note v7: DetectorElement::save/load gained eff_low/eff_cnt/eff_upp and
  ///           Grid2dBinGroup::save/load gained vec_eff_low/cnt/upp_group
  ///           (older checkpoints are rejected).
  static constexpr uint32_t PIPELINE_VERSION = 7;

  //======================================================================
  /// @brief Geometry state from the init, load_parameters, and build_geometry stages.
  struct GeoState {
    Grid2dPillar g2pil_naive;          ///< 2D terrain grid (naive).
    Grid3dVoxel  g3vox_input;          ///< 3D voxel grid (input density).
    Grid3dVoxel  g3vox_merged_input;   ///< 3D voxel grid (merged input).
    Grid2dPillar g2pil_shell_upper;    ///< Upper shell terrain (above g3vox).
    Grid2dPillar g2pil_shell_lower;    ///< Lower shell terrain (below g3vox).
    Grid2dPillar g2pil_shell_lateral;  ///< Lateral shell terrain (AABB-outside columns).
    bool has_shell_upper = false;      ///< True if upper shell is valid.
    bool has_shell_lower = false;      ///< True if lower shell is valid.
    bool has_shell_lateral = false;    ///< True if lateral shell is valid.
  };

  //======================================================================
  /// @brief Detector state from the trace_path_lengths and compute_prior stages (optional).
  struct DetState {
    DetectorPanelArray arrdet_g2pil_naive; ///< DetectorPanelArray from g2pil_naive.
    DetectorPanelArray arrdet_g3vox_input; ///< DetectorPanelArray from g3vox_input.
    std::vector<SpMatf> vec_spmat_PL;      ///< Path-length sparse matrices per detector.
    exemdl::build_prior::ShellPL shell_pl; ///< Shell path-length vectors (density-independent).
    exemdl::build_prior::PriorInfoAll prior_info_all; ///< Prior info (lower/center/upper).
    std::array<double, 3> avr_dens_lower_center_upper{}; ///< Volume-weighted average density.
    /// Density quadruplet: [prior, shell_upper, shell_lower, shell_lateral] (kg/m^3).
    std::array<double, 4> density_quad{};
  };

  //======================================================================
  /// @brief Matrix state from the build_observation_matrix stage (optional).
  struct MatState {
    Eigen::MatrixXf mat_dNdD_grouped_lower;  ///< dN/dD matrix (lower bound).
    Eigen::MatrixXf mat_dNdD_grouped_center; ///< dN/dD matrix (center).
    Eigen::MatrixXf mat_dNdD_grouped_upper;  ///< dN/dD matrix (upper bound).
  };

  //======================================================================
  /// @brief Pipeline state containing results from all modules.
  /// @details
  /// Holds geometry data (init/load_parameters/build_geometry), detector data (trace_path_lengths/compute_prior), and matrix data (build_observation_matrix).
  /// Optional fields (det, mat) are populated only after the corresponding module runs.
  struct State {
    // init / load_parameters results
    exemdl::load_parameters::AppParameters app_params; ///< Application parameters.
    nlohmann::json js;                                  ///< Original JSON configuration.
    fs::path path_json;                                 ///< Path to JSON file.

    // build_geometry results
    GeoState geom;                                      ///< Geometry state.

    // trace_path_lengths / compute_prior results (optional)
    std::optional<DetState> det;                        ///< Detector state (if trace_path_lengths/compute_prior ran).

    // build_observation_matrix results (optional)
    std::optional<MatState> mat;                        ///< Matrix state (if build_observation_matrix ran).

    int completed_module = 0;                           ///< Last completed module number.
  };

  //======================================================================
  /// @name Pipeline functions
  /// @{

  /// @brief Save pipeline state to directory.
  /// @param[in] state State to save.
  /// @param[in] dir Output directory path.
  /// @throws std::runtime_error if files cannot be created.
  /// @note Writes `app_params.json` (human-readable copy of the full JSON5
  ///       configuration) alongside geometry and detector binaries.
  ///       Grid1d::save() stores canonical (already shifted) coordinates
  ///       as-is, so no shift information is needed to read them back.
  void save(const State& state, const fs::path& dir);

  /// @brief Load pipeline state from directory.
  /// @param[in] dir Input directory path.
  /// @return Loaded State.
  /// @throws std::runtime_error if files cannot be read or version mismatch.
  /// @note Checkpoint binaries hold canonical coordinates and are loaded
  ///       pass-through; a save/load roundtrip is bit-identical.
  State load(const fs::path& dir);

  /// @brief Initialize pipeline state by running init, load_parameters, and build_geometry.
  /// @param[in] json_path Path to JSON configuration file.
  /// @param[in] seed_given Whether seed is explicitly provided via CLI (default: false).
  /// @param[in] seed_value Seed value (only used if seed_given is true, default: 0).
  /// @return State with completed_module = 3.
  /// @throws std::runtime_error if JSON loading, validation, or geometry build fails.
  /// @note Runs: init, then load_parameters, then build_geometry.
  /// @note Thread safety: Not thread-safe (configures logger, OpenMP, Eigen threads).
  State init(const fs::path& json_path, bool seed_given = false, unsigned seed_value = 0);

  /// @brief Trace path lengths (density-independent, expensive raytrace).
  /// @details Builds arrdet_g2pil_naive, arrdet_g3vox_input, vec_spmat_PL, and shell_pl.
  ///          This is the expensive part that only needs to run once per geometry.
  /// @param[in,out] state Pipeline state (must have completed_module >= 3).
  /// @throws std::runtime_error if state.completed_module < 3.
  /// @note After execution, state.det is partially populated (no prior_info_all yet).
  /// @note Uses OpenMP internally for parallel path length computation.
  void run_trace_path_lengths(State& state);

  /// @brief Compute prior (density-dependent, cheap).
  /// @details Rebuilds prior_info_all from saved PL data + density_quad.
  ///          No ray tracing — uses linear combination of saved path lengths.
  /// @param[in,out] state Pipeline state (must have det with shell_pl and vec_spmat_PL).
  /// @param[in] density_quad Density quadruplet [prior, upper, lower, lateral] (kg/m^3).
  /// @throws std::runtime_error if state.det is empty.
  /// @note run_trace_path_lengths must be called first.
  /// @note After execution, state.det->prior_info_all is populated and state.completed_module = 5.
  void run_compute_prior(State& state, const std::array<double, 4>& density_quad);

  /// @brief Build observation matrix using prior_info_all's DL vectors.
  /// @param[in,out] state Pipeline state (must have completed_module >= 5).
  /// @throws std::runtime_error if state.completed_module < 5 or state.det is empty.
  /// @note After execution, state.mat is populated and state.completed_module = 6.
  /// @note Uses lower/center/upper DL vectors from prior_info_all respectively.
  void run_build_observation_matrix(State& state);

  /// @brief Build observation matrix with specified DL vectors.
  /// @param[in,out] state Pipeline state (must have completed_module >= 5).
  /// @param[in] vec_vecxf_DL Density-length vectors for each detector.
  /// @throws std::runtime_error if state.completed_module < 5 or state.det is empty.
  /// @note After execution, state.mat is populated and state.completed_module = 6.
  /// @note The same DL is used for lower/center/upper matrices (for iterative estimation).
  void run_build_observation_matrix(State& state, const std::vector<Eigen::VectorXf>& vec_vecxf_DL);

  /// @brief Calculate density-length vectors from reconstruction result.
  /// @param[in] state Pipeline state (must have state.det populated).
  /// @param[in] vecxf_density Reconstructed density vector (uqiv_avail indexed).
  /// @return Density-length vectors for each detector.
  /// @throws std::runtime_error if state.det is empty.
  /// @note Formula: DL = PL * density (sparse matrix-vector product).
  /// @note Used in iterative estimation to update prior DL from reconstruction result.
  std::vector<Eigen::VectorXf> calc_vec_vecxf_DL(
    const State& state, const Eigen::VectorXf& vecxf_density);

  /// @brief Invert density (MAP) with specified regularization parameters.
  /// @param[in] state Pipeline state (must have completed_module >= 6).
  /// @param[in] sigma_rho Density uncertainty parameter (kg/m^3).
  /// @param[in] corr_length Correlation length parameter (meters).
  /// @return Inversion results for lower/center/upper priors.
  /// @throws std::runtime_error if state.completed_module < 6, state.det or state.mat is empty.
  /// @note This function does NOT modify state (const reference).
  /// @note Can be called multiple times with different sigma_rho/corr_length values.
  exemdl::run_inversion::InversionResultsAll run_invert_density(
    const State& state, double sigma_rho, double corr_length);

  /// @brief Invert density (MAP) with sigma_rho_diag override.
  /// @param[in] state Pipeline state (must have completed_module >= 6).
  /// @param[in] sigma_rho Density uncertainty parameter (kg/m^3).
  /// @param[in] corr_length Correlation length parameter (meters).
  /// @param[in] sigma_rho_diag Diagonal covariance std dev (kg/m^3).
  ///            If negative, falls back to sigma_rho (auto behavior).
  /// @return Inversion results for lower/center/upper priors.
  /// @throws std::runtime_error if state.completed_module < 6, state.det or state.mat is empty.
  /// @note This function does NOT modify state (const reference).
  exemdl::run_inversion::InversionResultsAll run_invert_density(
    const State& state, double sigma_rho, double corr_length, double sigma_rho_diag);

  /// @brief Export invert_density inputs and results for external comparison tools.
  /// @details
  /// Writes self-describing binaries into the "rec" subdirectory of the current
  /// default output directory. Every .bin uses the same layout as mat/
  /// (uint64 rows/cols header, raw float32 column-major; see
  /// io_binary::write_matxf_stream), so external readers can reuse the mat/
  /// loader unchanged. Written files:
  /// - vec_nmuon_obs.bin: observed muon counts actually used by invert_density.
  ///   The Poisson values are precomputed reads stored in the detector panels,
  ///   so this reproduces the identical vector passed to the inversion.
  /// - vec_nmuon_prior_{lower,center,upper}.bin: prior muon count vectors.
  /// - vec_rho_prior_{lower,center,upper}.bin: prior density per voxel [kg/m^3].
  /// - vec_rho_true.bin: input (ground-truth) density per voxel [kg/m^3].
  /// - vec_rho_rec_center.bin / vec_dens_err_center.bin: reconstructed density
  ///   and its 1-sigma uncertainty per voxel (lower/upper variants are written
  ///   only when independently reconstructed; see InversionResultsAll flags).
  /// - mat_voxel_xyz.bin: voxel center coordinates, n_voxel x 3 [m].
  /// - manifest.json: file list, row/column order contract (identical to
  ///   mat/manifest.json), solver settings (sigma_rho, corr_length,
  ///   sigma_rho_diag, nmuon thresholds, anisotropy, efficiency-variance flag),
  ///   and chi2 of the center reconstruction.
  /// @param[in] state Pipeline state (must have det, mat, and geometry populated).
  /// @param[in] inv_result Inversion results from invert_density.
  /// @param[in] sigma_rho Density uncertainty parameter actually used (kg/m^3).
  /// @param[in] corr_length Correlation length actually used (meters).
  /// @param[in] sigma_rho_diag Diagonal covariance std dev actually used (kg/m^3).
  ///            Negative means invert_density fell back to sigma_rho.
  /// @throws std::runtime_error If preconditions fail, vector sizes disagree
  ///         with the matrix dimensions, or an output file cannot be opened.
  /// @note This function does NOT modify state and adds no new dependencies.
  void save_recon_io(
    const State& state,
    const exemdl::run_inversion::InversionResultsAll& inv_result,
    double sigma_rho, double corr_length, double sigma_rho_diag);

  /// @brief Analyze errors (NagaInvLooper + error calculation).
  /// @param[in] state Pipeline state (must have det populated).
  /// @param[in] inv_result Inversion results from the density inversion stage.
  /// @param[in] output_prefix Prefix for output filenames (e.g., "001").
  ///            Empty string for default filenames (backward compatible).
  /// @throws std::runtime_error if inv_result.center.pNagainv is null or state.det is empty.
  /// @note Appends prefix to error cross-section filenames to avoid collision during sweep.
  void run_analyze_errors(
    const State& state,
    const exemdl::run_inversion::InversionResultsAll& inv_result,
    const std::string& output_prefix = "");

  /// @}

} // namespace exemdl::pipeline
