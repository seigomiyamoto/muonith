/// @file st_NagaInvLooper.hpp
/// @brief Naga inversion detector disable mode looper
/// @details
/// Provides systematic studies by running muon tomography reconstruction (NagaInv)
/// with individual detectors disabled. This allows quantifying the contribution
/// and sensitivity of each detector to the overall reconstruction quality.
///
/// **Workflow:**
/// 1. Clone the reference NagaInv instance
/// 2. For each detector, disable its observations and rerun reconstruction
/// 3. Compare results to quantify detector importance
/// 4. Output cross-section images and error metrics
///
/// **Thread-safety:** Not thread-safe. Operates sequentially on detector sets.
#pragma once

#include <string>
#include <vector>
#include <memory>
#include <unordered_map>
#include <Eigen/Dense>

#include "ns_myapp.hpp"
#include "cls_Grid3dVoxel.hpp"
#include "cls_NagaInv.hpp"
#include "cls_NagaInvParameters.hpp"
#include "cls_DetectorIndexContainer.hpp"

/// @brief Looper for running NagaInv reconstruction with detectors systematically disabled
/// @details
/// This structure performs sensitivity analysis by disabling individual detectors
/// and re-running reconstruction. All members are const references to avoid deep copies
/// during setup; internal methods clone data as needed for each reconstruction run.
///
/// **Typical usage:**
/// @code
/// NagaInvLooper looper{
///   .name = "sensitivity_study",
///   .nagainv_org = base_reconstructor,
///   .prm_nagainv = params,
///   .g3vox_org = voxel_grid,
///   .prm_zcross = cross_section_params,
///   .dic = detector_index_container
/// };
/// auto results = looper.exec_disable_mode_all(prior_density, all_detector_result);
/// looper.output_disable_mode_all_as_cross_section(results, all_detector_result);
/// @endcode
struct NagaInvLooper{
  /// @brief Instance name used as prefix for output files
  std::string name;

  /// @brief Reference NagaInv instance (will be cloned for each disabled-detector run)
  const NagaInv& nagainv_org;

  /// @brief Parameters for NagaInv reconstruction (regularization, iteration limits, etc.)
  const NagaInvParameters& prm_nagainv;

  /// @brief Reference voxel grid (will be cloned for output operations)
  const Grid3dVoxel &g3vox_org;

  /// @brief Parameters for cross-section output (z-slice position, resolution, etc.)
  const Grid3dVoxel::CrossSectionZParameters &prm_zcross;

  /// @brief Detector index container mapping Detid to matrix rows
  const DetectorIndexContainer& dic;

  /// @brief Run reconstruction with one detector disabled
  /// @param[in] detid_disabled Detector ID to disable
  /// @param[in] vecxf_dens_prior Prior density (initial guess for reconstruction)
  /// @param[in] input_vecxf_rec_all Optional reference reconstruction using all detectors
  /// @param[in] precomputed_mat_cov_dens_inv Optional precomputed inverse of mat_cov_dens.
  ///            When provided, avoids redundant O(n^3) inversion of the density covariance matrix.
  /// @return Reconstruction result with disabled detector
  /// @details
  /// Clones nagainv_org, removes observations from detid_disabled by zeroing corresponding
  /// rows in mat_dNdD and mat_cov_muon, then runs reconstruction. The reference reconstruction
  /// (if provided) is used to compute vecxf_diff_from_real in the result.
  /// @throws std::runtime_error If reconstruction fails (logged as ERROR)
  NagaInv::ReconstResult exec_disable_mode(const Detid &detid_disabled
    , const Eigen::VectorXf &vecxf_dens_prior
    , const std::optional<Eigen::VectorXf>& input_vecxf_rec_all
    , const std::optional<Eigen::MatrixXf>& precomputed_mat_cov_dens_inv
      = std::nullopt);

  /// @brief Run reconstruction with each detector disabled sequentially
  /// @param[in] vecxf_dens_prior Prior density (initial guess for reconstruction)
  /// @param[in] input_vecxf_rec_all Optional reference reconstruction using all detectors
  /// @return Vector of reconstruction results, one per detector in dic
  /// @details
  /// Iterates over all detectors in dic.get_set_detid(), calling exec_disable_mode
  /// for each. Skips detectors that fail reconstruction (logs ERROR and continues).
  /// Results are ordered by iteration order over the set (not guaranteed to match Detid order).
  /// @note Returns empty vector if dic contains no detectors (logs WARN)
  std::vector<NagaInv::ReconstResult> exec_disable_mode_all(
    const Eigen::VectorXf &vecxf_dens_prior
  , const std::optional<Eigen::VectorXf>& input_vecxf_rec_all);

  /// @brief Output reconstruction results as cross-section images
  /// @param[in] vec_reconst_res Reconstruction results from exec_disable_mode_all
  /// @param[in] input_vecxf_rec_all Optional reference reconstruction using all detectors
  /// @details
  /// For each result in vec_reconst_res (indexed 0, 1, 2, ...), writes three cross-section files:
  /// - `{name}_disable_{idx:02}`: Reconstructed density
  /// - `{name}_disable_{idx:02}_delta_prior`: Difference from prior
  /// - `{name}_disable_{idx:02}_diff_from_all`: Difference from reference (if provided)
  ///
  /// The index corresponds to iteration order from exec_disable_mode_all, not necessarily Detid order.
  /// @note Skips diff_from_all output if input_vecxf_rec_all is not provided (logs WARN per detector)
  void output_disable_mode_all_as_cross_section(
    const std::vector<NagaInv::ReconstResult> &vec_reconst_res
  , const std::optional<Eigen::VectorXf>& input_vecxf_rec_all = std::nullopt) const;

  /// @brief Calculate maximum deviation across disabled-detector reconstructions
  /// @param[in] n_vecxf_dens Number of density elements (voxels)
  /// @param[in] vec_reconst_disable_mode_res Reconstruction results from exec_disable_mode_all
  /// @return Vector (size n_vecxf_dens) where each element is the maximum absolute deviation
  ///         of that voxel's density across all disabled-detector reconstructions
  /// @details
  /// For each voxel i, computes max_j |reconst_res[j].vecxf_delta_dens_prior[i]|,
  /// representing the worst-case density error caused by disabling any single detector.
  /// @note This method is declared but not implemented in st_NagaInvLooper.cpp.
  ///       Actual implementation exists in ns_mix_error namespace.
  Eigen::VectorXf calc_max_disabled_error(
      const int n_vecxf_dens
    , const std::vector<NagaInv::ReconstResult> &vec_reconst_disable_mode_res);

  /// @brief Calculate combined error from prior, statistics, and disabled-detector sensitivity
  /// @param[in] vecxf_prior Prior density error distribution
  /// @param[in] vecxf_mat_cov_dens_dash Statistical error from muon statistics
  /// @param[in] vecxf_disabled_error Systematic error from calc_max_disabled_error
  /// @return Vector (same size as inputs) of combined errors via root-sum-square:
  ///         sqrt(prior^2 + stat^2 + disabled^2) element-wise
  /// @details
  /// Combines three independent error sources using quadrature addition.
  /// All input vectors must have the same size (number of voxels).
  /// @note This method is declared but not implemented in st_NagaInvLooper.cpp.
  ///       Actual implementation exists in ns_mix_error namespace.
  Eigen::VectorXf calc_mixed_error(const Eigen::VectorXf &vecxf_prior
    , const Eigen::VectorXf &vecxf_mat_cov_dens_dash
    , const Eigen::VectorXf &vecxf_disabled_error);

};
