/// @file exemdl_calc_matrix.hpp
/// @brief Path length calculation and observation matrix generation for muon tomography
/// @details
/// This module provides functions for computing muon path lengths through voxels and
/// constructing the observation matrix (dN/dD) used in muon tomography density inversion.
///
/// ## Typical Workflow
/// 1. Build detector arrays and voxel grids using exemdl::build_detector
/// 2. Prepare prior density-length (DL) vectors for each detector
/// 3. Call build_all() to generate observation matrices for lower/center/upper bounds
/// 4. Use the resulting matrices in inversion algorithms
///
/// ## Key Concepts
/// - **NPL Matrix (dN/dD):** Observation matrix relating muon counts to density
/// - **Prior DL Vectors:** Density-length products computed from prior density model
/// - **Sparse PL Matrices:** Path length matrices stored in sparse format for efficiency
///
/// ## Thread Safety
/// - build_matrix() and build_all() are NOT thread-safe; they modify shared state
/// - Each call should be made from a single thread
///
/// @see exemdl::build_detector for detector and voxel grid construction
/// @see calc_dNdD for underlying matrix computation routines
#pragma once

#include <vector>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include "exemdl_load_parameters.hpp"
#include "cls_DetectorPanelArray.hpp"
#include "cls_Grid2dPillar.hpp"
#include "cls_Grid3dVoxel.hpp"
#include "cls_FluxTable.hpp"
#include "ns_pathcalc.hpp"
#include "cls_MatrixBuildParameters.hpp"
#include "exemdl_build_detector.hpp"

/// @namespace exemdl::calc_matrix
/// @brief Namespace for path calculation and observation matrix generation
/// @ingroup ExecModule
namespace exemdl::calc_matrix {

  /// @brief Sparse matrix type alias (single precision)
  using SpMatf = Eigen::SparseMatrix<float>;

  /// @brief Input arguments for observation matrix construction
  /// @details Aggregates references to detector build results, voxel grid, and parameters.
  /// All references must remain valid throughout the build_matrix() or build_all() call.
  struct BuildArgs {
    /// @brief Results from exemdl::build_detector module
    exemdl::build_detector::BuildResult& res_built_det;
    /// @brief Merged voxel grid covering the imaging volume
    Grid3dVoxel& g3vox_input;
    /// @brief Application parameters including path and detector settings
    const exemdl::load_parameters::AppParameters& app_params;
  };

  /// @brief Build a single observation matrix (dN/dD) from prior DL vectors
  /// @param[in,out] args Build arguments containing detector results and settings
  /// @param[in] vec_vecxf_DL_prior_in Prior density-length vectors for each detector
  /// @return Grouped observation matrix (dN/dD) with rows = groups, cols = voxels
  /// @note Thread-safe: No. Modifies internal state of args.res_built_det.
  /// @note Complexity: O(n_detectors * n_voxels * n_rays) approximately
  Eigen::MatrixXf build_matrix(BuildArgs& args,
    const std::vector<Eigen::VectorXf>& vec_vecxf_DL_prior_in);

  /// @brief Output results from observation matrix construction
  /// @details Contains observation matrices computed for three prior density scenarios
  /// (lower bound, center, upper bound) to support uncertainty quantification.
  struct BuildResult {
    /// @brief Observation matrix computed with lower-bound prior density
    Eigen::MatrixXf mat_dNdD_grouped_lower;
    /// @brief Observation matrix computed with center prior density
    Eigen::MatrixXf mat_dNdD_grouped_center;
    /// @brief Observation matrix computed with upper-bound prior density
    Eigen::MatrixXf mat_dNdD_grouped_upper;
  };

  /// @brief Build observation matrices for all prior scenarios (lower/center/upper)
  /// @param[in,out] args Build arguments containing detector results and settings
  /// @return BuildResult containing three observation matrices
  /// @note Thread-safe: No. Calls build_matrix() three times sequentially.
  /// @note This is a convenience function that wraps three build_matrix() calls.
  BuildResult build_all(BuildArgs& args);

} // namespace exemdl::calc_matrix
