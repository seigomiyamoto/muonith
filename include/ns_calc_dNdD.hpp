/// @file ns_calc_dNdD.hpp
/// @brief Functions for dN/dD matrix construction from path length matrices
/// @details Provides functions to create muon number per density matrices from path length matrices,
/// with support for both dense and sparse representations and angle bin grouping.
/// Extracted from pathcalc::matrix namespace as an independent top-level namespace.
#pragma once

#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "cls_DetectorPanelArray.hpp"
#include "cls_Grid3dVoxel.hpp"
#include "cls_FluxTable.hpp"

namespace pathcalc {
  class Parameters; // forward declaration
  class MatrixBuildParameters; // forward declaration
};

//===================================================
/// @namespace calc_dNdD
/// @brief Functions for path length matrix generation
/// @details Provides functions to create muon number per density matrices from path length matrices,
/// with support for both dense and sparse representations and angle bin grouping.
//===================================================
namespace calc_dNdD {
  using SpMatf = Eigen::SparseMatrix<float>;

  /// @brief Create mat_dNdD matrix from mat_PL and dFdR_R_costhz in a detector unit
  /// @details Creates muon number per density matrix from path length matrix.
  /// @param panel Reference to DetectorPanel
  /// @param mat_PL Path length matrix
  /// @param vecxf_init_dens Initial density vector
  /// @param ft FluxTable instance
  /// @return Muon number per density matrix
  /// @note Uses OpenMP
  Eigen::MatrixXf mp_make_mat_dNdD(
    const DetectorPanel& panel, const Eigen::MatrixXf &mat_PL
  , const Eigen::VectorXf &vecxf_init_dens, const FluxTable &ft );

  /// @brief Create spmat_dNdD matrix from spmat_PL and dFdR_R_costhz in a detector unit
  /// @details Creates muon number per density matrix from path length matrix (sparse version).
  /// @param panel Reference to DetectorPanel
  /// @param spmat_PL Sparse path length matrix
  /// @param vecxf_init_dens Initial density vector
  /// @param ft FluxTable instance
  /// @return Muon number per density matrix (sparse)
  /// @note Uses OpenMP
  SpMatf mp_make_spmat_dNdD(
    const DetectorPanel& panel, const SpMatf &spmat_PL
  , const Eigen::VectorXf &vecxf_init_dens, const FluxTable &ft );


  /// @brief Convert mat_dNdD to number per density matrix grouped by angle bins
  /// @param panel Pointer to DetectorPanel containing merging information
  /// @param mat_dNdD Eigen::MatrixXf of number path length matrix
  /// @return Eigen::MatrixXf storing grouped number per density matrix
  Eigen::MatrixXf make_grouped_mat_dNdD(
      const DetectorPanel& panel, const Eigen::MatrixXf &mat_dNdD );

  /// @brief Convert mat_dNdD to number per density matrix grouped by angle bins
  /// @param panel Pointer to DetectorPanel containing merging information
  /// @param mat_dNdD Eigen::MatrixXf of number path length matrix
  /// @return Eigen::MatrixXf storing grouped number per density matrix
  /// @note Uses OpenMP, but not significantly faster
  Eigen::MatrixXf mp_make_grouped_mat_dNdD(
      const DetectorPanel& panel, const Eigen::MatrixXf &mat_dNdD );

  /// @brief Convert spmat_dNdD to number per density matrix (dense matrix) grouped by angle bins
  /// @param panel Pointer to DetectorPanel containing merging information
  /// @param spmat_dNdD Sparse number path length matrix
  /// @return Eigen::MatrixXf storing grouped number per density matrix
  /// @note Uses OpenMP. Uses InnerIterator to visit only non-zero elements
  ///       instead of converting sparse columns to dense vectors.
  Eigen::MatrixXf mp_make_grouped_mat_dNdD(
      const DetectorPanel& panel, const SpMatf &spmat_dNdD );

  /// @brief Legacy version of mp_make_grouped_mat_dNdD (sparse overload)
  /// @details Converts each sparse column to dense VectorXf before grouping.
  ///          Retained for regression testing.
  /// @note Uses OpenMP
  Eigen::MatrixXf mp_make_grouped_mat_dNdD_legacy(
      const DetectorPanel& panel, const SpMatf &spmat_dNdD );

  /// @brief Calculate mat_PL for g3vox for all detectors and calculate expected signal value from initial density
  /// @details Workflow for each detector:
  /// 1. Execute g3vox::mp_make_mat_PL(det)
  /// 2. Execute panel.mp_calc_set_peneflux_signal_from_DL(ft) to calculate signal amount
  /// 3. Execute panel.copy_signal_noise_to_g2bg() to initialize BinGroup
  /// 4. Execute panel.auto_divide_by_signal_noise_group_all(panel.get_prm_bingrp()) to optimize BinGroup
  /// 5. Execute panel.calc_set_vec_signal_noise_group() to calculate signal/noise_poi_group and set to vec_signal/noise_poi_group
  /// 6. Store calculated mat_PL in vector
  /// 7. Repeat above for all detectors
  /// @param arrdet Reference to DetectorPanelArray
  /// @param g3vox Reference to Grid3dVoxel
  /// @param ft FluxTable instance
  /// @param prm pathcalc::Parameters
  /// @param tf_apply_eff Whether to apply efficiency
  /// @return Vector of mat_PL for each detector
  /// @note Uses OpenMP
  std::vector<Eigen::MatrixXf> calc_mat_PL_and_set_signal(
      DetectorPanelArray &arrdet, Grid3dVoxel &g3vox
    , const FluxTable &ft, const pathcalc::Parameters &prm
    , const bool tf_apply_eff );

  /// @brief Perform make_grouped_mat_dNdD for all detectors and create final matrix
  /// @param arrdet Reference to DetectorPanelArray
  /// @param g3vox Reference to Grid3dVoxel
  /// @param vec_spmat_PL Vector of mat_PL for each detector
  /// @param vec_vecxf_DL_prior Reference to density length information at prior density
  /// @param ft Instance of FluxTable
  /// @return Eigen::MatrixXf storing the final number per density matrix grouped by angle bins for all detectors (alldet)
  /// @note Uses OpenMP
  Eigen::MatrixXf
    make_grouped_alldet_mat_dNdD(
      DetectorPanelArray &arrdet, Grid3dVoxel &g3vox
    , const std::vector<SpMatf> &vec_spmat_PL
    , const std::vector<Eigen::VectorXf>& vec_vecxf_DL_prior
    , const FluxTable &ft );

  /// @brief Perform make_grouped_spmat_dNdD for all detectors and create final matrix
  /// @param arrdet Reference to DetectorPanelArray
  /// @param g3vox Reference to Grid3dVoxel
  /// @param vec_spmat_PL Vector of mat_PL for each detector
  /// @param vec_vecxf_DL_prior Reference to density length information at prior density
  /// @param ft Instance of FluxTable
  /// @param tf_apply_eff Whether to apply efficiency
  /// @return Eigen::MatrixXf storing the final number per density matrix grouped by angle bins for all detectors (alldet)
  /// @note Uses OpenMP
  Eigen::MatrixXf
    make_grouped_alldet_mat_dNdD_sprs(
      DetectorPanelArray &arrdet, Grid3dVoxel &g3vox
    , const std::vector<SpMatf> &vec_spmat_PL
    , const std::vector<Eigen::VectorXf>& vec_vecxf_DL_prior
    , const FluxTable &ft
    , const bool tf_apply_eff
    , const bool tf_apply_eff_cnt = false );

  /// @brief Factory function to create numbered path length matrix containing information for all detectors and voxels
  /// @param prm_mat Reference to pathcalc::MatrixBuildParameters
  /// @return Eigen::MatrixXf storing numbered path length matrix
  /// @note Reference values (arrdet, g3vox, vec_spmat_PL) stored in prm_mat will be modified.
  /// @details \image html images/chart_pathcalc_matrix_create_large_dNdD_matrix.dio.png
  Eigen::MatrixXf create_grouped_mat_dNdD_alldet( const pathcalc::MatrixBuildParameters& prm_mat );

  /// @brief Factory function to create numbered path length matrix containing information for all detectors and voxels
  /// @param prm_mat Reference to pathcalc::MatrixBuildParameters
  /// @return Eigen::MatrixXf storing numbered path length matrix
  /// @note Sparse matrix version.
  /// Reference values (arrdet, g3vox, vec_spmat_PL) stored in prm_mat will be modified.
  /// @details \image html images/chart_pathcalc_matrix_create_large_dNdD_matrix.dio.png
  Eigen::MatrixXf create_grouped_mat_dNdD_alldet_sprs( const pathcalc::MatrixBuildParameters& prm_mat );

};
