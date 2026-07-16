/// @file exemdl_run_inversion.hpp
/// @brief Linear inversion execution module
/// @details
/// This module provides functions for executing muon radiography linear inversion
/// with prior error analysis. It builds NagaInv instances, runs density reconstruction
/// using BLAS/LAPACK-accelerated solvers, and outputs cross-sectional visualizations.
///
/// Workflow:
/// 1. Create NagaInvManager with inversion parameters
/// 2. Build NagaInv instances for lower/center/upper prior models
/// 3. Execute density reconstruction via BLAS/LAPACK-based linear solver
/// 4. Write reconstruction results and differences to voxel grids
/// 5. Output cross-sections for visualization and analysis
///
/// @note Thread-safety: Functions in this module are not thread-safe due to
/// internal use of NagaInvManager state and file I/O operations.
///
/// @note Units: Density values are in g/cm³. Spatial coordinates (x, y, z) are in meters.
#pragma once

#include <array>
#include <Eigen/Dense>
#include "cls_Grid3dVoxel.hpp"
#include "cls_NagaInvManager.hpp"
#include "cls_NagaInv.hpp"
#include "exemdl_load_parameters.hpp"
#include "exemdl_build_detector.hpp"
#include "exemdl_calc_matrix.hpp"

/// @namespace exemdl::run_inversion
/// @brief Namespace for inversion execution related functions
/// @ingroup ExecModule
namespace exemdl::run_inversion {

  /// @brief Arguments for building linear inversion
  /// @details
  /// This structure aggregates all necessary inputs for running the inversion process.
  /// It contains observed muon data, prior models, sensitivity matrices, and configuration
  /// parameters.
  ///
  /// Typical usage:
  /// @code
  /// BuildArgs args{
  ///   .index_run = 0,
  ///   .prior_info_all = prior_info,
  ///   .g3vox_input = input_grid,
  ///   .res_mat = matrix_results,
  ///   .vecxf_nmuon_obs = observed_muons,
  ///   .tf_use_dens_input = true,  // for simulation
  ///   .params = app_params
  /// };
  /// auto results = build_all(args);
  /// @endcode
  struct BuildArgs {
    const int                             index_run; ///< Run index for multi-run scenarios
    const exemdl::build_prior::PriorInfoAll& prior_info_all; ///< Prior information for lower/center/upper bounds
    const Grid3dVoxel&                    g3vox_input; ///< Input model voxel grid (ground truth for simulations)
    const exemdl::calc_matrix::BuildResult&  res_mat; ///< Sensitivity matrix results for lower, center, upper priors
    const Eigen::VectorXf&                vecxf_nmuon_obs; ///< Observed muon count vector
    const bool tf_use_dens_input = true; ///< Whether to use input density. True for simulations with known input, false for real observation data with unknown input.
    const exemdl::load_parameters::AppParameters& params;   ///< Common application parameters
    const std::array<double,4> density_quad{}; ///< Shell density [prior, upper, lower, lateral]. Default {0,0,0,0} means no shell classification.
    Eigen::VectorXf vecxf_var_eff{}; ///< Per-row efficiency variance added to the C_N diagonal. Empty (default) disables it.
  };

  /// @brief Container for linear inversion execution results
  /// @details
  /// Holds all outputs from a single inversion run, including the reconstruction
  /// result, voxel grids with various density fields, and the NagaInv instance
  /// used for the reconstruction.
  ///
  /// @note Memory ownership: The pNagainv unique pointer owns the NagaInv instance.
  /// Voxel grids (g3vox_*) are value types and own their data.
  struct InversionResults {
    NagaInv::ReconstResult reconst_res; ///< Reconstruction result with density vectors
    std::unique_ptr<NagaInv> pNagainv; ///< NagaInv instance (owned). May be null if cloned from center.
    Grid3dVoxel::CrossSectionZParameters prm_zcross; ///< Cross-section output parameters
    Grid3dVoxel g3vox_rec;         ///< Voxel grid with reconstructed density (g/cm³)
    Grid3dVoxel g3vox_delta_prior; ///< Voxel grid with difference: reconstructed - prior density (g/cm³)
    Grid3dVoxel g3vox_diff_real;   ///< Voxel grid with difference: reconstructed - input density (g/cm³)
    Grid3dVoxel g3vox_dens_err;    ///< Voxel grid with reconstruction density uncertainty (diagonal of covariance matrix, g/cm³)
    double chi2_muon = 0.0;        ///< Muon chi-square of the reconstructed model: (N_obs - N(rho'))^T C_N^-1 (N_obs - N(rho')) [dimensionless]
    double ndf_muon  = 0.0;        ///< Degrees of freedom = number of observation bins N_obs (Nishiyama et al. 2017). chi2/ndf is derived (not stored)
    double p_eff_muon = 0.0;       ///< Effective number of parameters p_eff = trace(R) (reference correction; 0 unless tf_calc_chi2ndf)
  };

  /// @brief Build inversion results for a single prior model
  /// @details
  /// Executes the complete inversion workflow for a single prior model:
  /// 1. Creates NagaInvManager and NagaInv instance
  /// 2. Runs BLAS/LAPACK-based density reconstruction
  /// 3. Writes cross-sections for input, prior, reconstructed, and difference fields
  /// 4. Returns InversionResults with all outputs
  ///
  /// This function performs file I/O to output cross-sections and logs progress
  /// via spdlog.
  ///
  /// @param args Build arguments containing observed data and configuration
  /// @param name_NagaInv Name identifier for the NagaInv instance (e.g., "lower", "center", "upper")
  /// @param vecxf_nmuon_prior_in Prior muon count vector
  /// @param mat_grouped_dNdD_prior_in Grouped prior sensitivity matrix (dN/dD). Each row corresponds to a detector group.
  /// @param g3vox_prior Prior voxel grid information
  /// @return InversionResults containing reconstruction outputs and voxel grids
  ///
  /// @note Thread-safety: Not thread-safe due to file I/O and NagaInvManager state.
  /// @note Performance: Uses BLAS/LAPACK for linear algebra. Performance scales with voxel count and detector count.
  InversionResults
    build_inv_res(
        const BuildArgs& args
      , const std::string& name_NagaInv
      , const Eigen::VectorXf& vecxf_nmuon_prior_in
      , const Eigen::MatrixXf& mat_grouped_dNdD_prior_in
      , const Grid3dVoxel& g3vox_prior);

  /// @brief Container for all inversion reconstruction results (lower/center/upper priors)
  /// @details
  /// Aggregates reconstruction results for all three prior models (lower, center, upper).
  /// When tf_prior_error is false, only the center reconstruction is performed,
  /// and lower/upper slots contain copies of the center result for downstream compatibility.
  ///
  /// @note Validity flags: has_lower and has_upper indicate whether independent
  /// reconstructions were performed for lower/upper priors. If false, those slots
  /// contain mirrored center results.
  struct InversionResultsAll {
    InversionResults lower;   ///< Reconstruction result for lower bound prior model
    InversionResults center;  ///< Reconstruction result for center prior model
    InversionResults upper;   ///< Reconstruction result for upper bound prior model
    bool has_lower = false;   ///< Whether lower bound reconstruction is valid (independent run, not mirrored)
    bool has_upper = false;   ///< Whether upper bound reconstruction is valid (independent run, not mirrored)
  };

  /// @brief Run inversion process for all prior models
  /// @details
  /// Executes the complete inversion workflow:
  /// - If tf_prior_error is false: only center prior is reconstructed, and lower/upper
  ///   slots are populated with mirrored center results.
  /// - If tf_prior_error is true: all three priors (lower, center, upper) are independently
  ///   reconstructed.
  ///
  /// Each reconstruction calls build_inv_res internally, which performs BLAS/LAPACK-based
  /// density reconstruction and outputs cross-sections.
  ///
  /// @param args Build arguments containing observed data, priors, and configuration
  /// @return InversionResultsAll containing results for lower, center, and upper priors
  ///
  /// @note Thread-safety: Not thread-safe due to file I/O and NagaInvManager state.
  /// @note Performance: When tf_prior_error is true, runtime is approximately 3x longer
  /// due to three independent reconstructions.
  InversionResultsAll build_all( const BuildArgs& args );

}; // namespace exemdl::run_inversion
