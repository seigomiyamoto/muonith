/// @file exemdl_build_detector.hpp
/// @brief Detector array construction module for muon tomography.
/// @details Provides functions for building detector panel arrays (DetectorPanelArray)
///          with signal/noise characteristics and detection efficiency applied.
///          This module handles multi-stage construction:
///          1. Build arrdet_g2pil_naive from 2D terrain grid (Grid2dPillar)
///          2. Build arrdet_g3vox_input from 3D voxel grid (Grid3dVoxel) with sparse matrix generation
///          3. Compute volume-weighted average density for prior estimation
///          4. Build prior info (lower/center/upper) for inversion
///
/// @note Coordinate system: z-up, right-handed. All densities in kg/m^3.
/// @note Thread safety: Functions use OpenMP internally (mp_* calls). Not reentrant.
#pragma once

#include <filesystem>
#include <fstream>
#include <vector>
#include <utility>
#include <Eigen/Sparse>
#include "cls_DetectorPanelArray.hpp"
#include "cls_Grid2dPillar.hpp"
#include "cls_Grid3dVoxel.hpp"
#include "cls_FluxTable.hpp"
#include "exemdl_load_parameters.hpp"
#include "exemdl_build_prior.hpp"

/// @namespace exemdl::build_detector
/// @brief Namespace for building detector panel arrays from geometry data.
/// @ingroup ExecModule
namespace exemdl::build_detector {

  /// @brief Sparse matrix type alias (float precision).
  using SpMatf = Eigen::SparseMatrix<float>;

  /// @brief Input parameters for the detector build pipeline (the trace_path_lengths stage).
  struct BuildArgs {
    const exemdl::load_parameters::AppParameters& app_params; ///< Application parameter set.
    const Grid2dPillar& g2pil_naive; ///< Grid2dPillar derived from terrain data.
    Grid3dVoxel&        g3vox_input; ///< Grid3dVoxel with input density (merged or not).
    const Grid2dPillar& g2pil_shell_upper; ///< Upper shell Grid2dPillar (terrain above g3vox).
    const Grid2dPillar& g2pil_shell_lower; ///< Lower shell Grid2dPillar (terrain below g3vox).
    const Grid2dPillar& g2pil_shell_lateral; ///< Lateral shell Grid2dPillar (terrain surrounding g3vox).
    const bool has_shell_upper; ///< True if upper shell has valid cuboids.
    const bool has_shell_lower; ///< True if lower shell has valid cuboids.
    const bool has_shell_lateral; ///< True if lateral shell has valid cuboids.
    const nlohmann::json& js_proj_dens; ///< JSON object for "PROJ_DENS_EVAL_GROUPED".
  };

  /// @brief Build a DetectorPanelArray linked with Grid2dPillar (2D terrain).
  /// @param[in] app_params Application parameters.
  /// @param[in] g2pil_naive Grid2dPillar based on terrain data.
  /// @param[in] arrdet_template Empty DetectorPanelArray template for construction.
  /// @param[in] js_proj_dens JSON object for projected density evaluation.
  /// @param[in] ft_prior FluxTable for prior flux data.
  /// @return Constructed DetectorPanelArray with path/density lengths computed.
  /// @note Uses OpenMP for parallel path length calculation.
  DetectorPanelArray build_arrdet_g2pil(
      const exemdl::load_parameters::AppParameters& app_params
    , const Grid2dPillar& g2pil_naive
    , const DetectorPanelArray& arrdet_template
    , const nlohmann::json& js_proj_dens
    , const FluxTable&      ft_prior
  );

  /// @brief Build a DetectorPanelArray linked with Grid3dVoxel, with sparse matrix generation.
  /// @param[in] app_params Application parameters.
  /// @param[in,out] g3vox_input Merged Grid3dVoxel with input density.
  /// @param[in] g2pil_shell_upper Pre-built upper Grid2dPillar shell (terrain above g3vox).
  /// @param[in] g2pil_shell_lower Pre-built lower Grid2dPillar shell (terrain below g3vox).
  /// @param[in] has_shell_upper True if upper shell has valid cuboids.
  /// @param[in] has_shell_lower True if lower shell has valid cuboids.
  /// @param[in] g2pil_shell_lateral Pre-built lateral Grid2dPillar shell (terrain surrounding g3vox).
  /// @param[in] has_shell_lateral True if lateral shell has valid cuboids.
  /// @param[in] opt_arrdet_g2pil_to_copy_bimap Optional arrdet_g2pil to copy bimap to return value.
  /// @param[in] arrdet_template Empty DetectorPanelArray template.
  /// @return Tuple of (constructed DetectorPanelArray, sparse path-length matrices per detector,
  ///         flat per-element PL through n_hit_det-filtered voxels).
  /// @note Uses OpenMP for parallel computation.
  std::tuple<DetectorPanelArray, std::vector<SpMatf>, Eigen::VectorXf>
    build_arrdet_g3vox(
        const exemdl::load_parameters::AppParameters& app_params
      , Grid3dVoxel& g3vox_input
      , const Grid2dPillar& g2pil_shell_upper
      , const Grid2dPillar& g2pil_shell_lower
      , const Grid2dPillar& g2pil_shell_lateral
      , const bool has_shell_upper
      , const bool has_shell_lower
      , const bool has_shell_lateral
      , const std::optional<DetectorPanelArray>& opt_arrdet_g2pil_to_copy_bimap
      , const DetectorPanelArray& arrdet_template
    );

  /// @brief Compute volume-weighted average density (lower, center, upper bounds).
  /// @param[in,out] arrdet_g3vox_input DetectorPanelArray generated from g3vox_input.
  /// @param[in] js_proj_dens JSON object for "PROJ_DENS_EVAL_GROUPED".
  /// @param[in] ft_prior FluxTable for prior flux data.
  /// @param[in] tf_prior_error If true, compute lower/center/upper; else center only.
  /// @param[in] tf_eff If true, add efficiency-uncertainty variance in quadrature to the
  ///            Poisson term of the projected-density error band (gated by tf_eff_cn_diag).
  /// @param[in] tf_eff_independent If true, treat per-element efficiency uncertainties as
  ///            independent (sum of squares); else common mode (square of sum).
  /// @note noise term is always det floor + Poisson(poi bucket), driven by NOISE_PARAMETERS ratios.
  /// @return Array of [lower, center, upper] volume-weighted average densities (kg/m^3).
  /// @note Uses OpenMP for parallel density computation.
  std::array<double, 3>
    calc_volume_weighted_average_density_lower_center_upper(
      DetectorPanelArray& arrdet_g3vox_input
    , const nlohmann::json& js_proj_dens
    , const FluxTable& ft_prior
    , bool tf_prior_error
    , bool tf_eff = false
    , bool tf_eff_independent = false
  );

  // Type aliases from exemdl::build_prior (for backward compatibility).
  using PriorInfo    = exemdl::build_prior::PriorInfo;
  using PriorInfoAll = exemdl::build_prior::PriorInfoAll;
  using ShellPL      = exemdl::build_prior::ShellPL;

  /// @brief Execution result of the trace_path_lengths stage (detector build pipeline).
  struct BuildResult {
    DetectorPanelArray  arrdet_g2pil_naive; ///< DetectorPanelArray generated from g2pil_naive.
    DetectorPanelArray  arrdet_g3vox_input; ///< DetectorPanelArray generated from g3vox_input.
    std::vector<SpMatf> vec_spmat_PL;       ///< Path-length sparse matrices per detector.
    PriorInfoAll prior_info_all;            ///< Prior Grid3dVoxel and DetectorPanelArray triplet.
    ShellPL shell_pl;                       ///< Shell path-length vectors (density-independent).
    std::array<double, 3> avr_dens_lower_center_upper; ///< Volume-weighted average density [lower, center, upper] (kg/m^3).
  };

  /// @brief Execute full detector build pipeline with sparse matrix generation.
  /// @param[in] args Input parameters and geometry information.
  /// @return BuildResult containing all construction outputs.
  /// @details Orchestrates: build_arrdet_g2pil -> build_arrdet_g3vox ->
  ///          calc_volume_weighted_average_density -> build_prior_info_all.
  BuildResult build_all( const BuildArgs& args );

} // namespace exemdl::build_detector
