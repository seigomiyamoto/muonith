/// @file ns_detector_indexing.hpp
/// @brief Detector indexing utilities for muon tomography
/// @details
/// This module provides helper functions for detector element indexing and ID management
/// in muon tomography applications. It supports:
/// - Column extraction from matrices based on index maps
/// - Row and block selection from matrices
/// - Detector ID to unique index mapping
/// - Matrix subsetting by disabling specific detector IDs
///
/// Thread-safety: Functions using OpenMP are documented individually.
/// Coordinate system: Not applicable (operates on abstract detector indices).
/// Units: Dimensionless indices.
#pragma once

#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>
#include <chrono>
#include <thread>

#define _USE_MATH_DEFINES // for M_PI, must be defined before including cmath
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Sparse>

// spdlog logger
#include "spdlog_pch.hpp"
#include "ns_mylogger.hpp"
#include "ns_mymacro.hpp"
#include "ns_type_definitions.hpp"
#include "st_NmuonVectors.hpp"
#include "cls_Grid3d.hpp"

/// @brief Sparse matrix type definition
using SpMatf = Eigen::SparseMatrix<float>;

/// @brief for type definitions
using namespace index_type_definitions;

/// @brief Namespace for detector indexing related functions
/// @details
/// Provides utilities for managing detector element indices and extracting subsets
/// of matrices corresponding to available or disabled detector elements. Functions
/// support both dense and sparse matrix representations using Eigen.
///
/// @ingroup detectorClasses
namespace detector_indexing {

  /// @brief Extract columns specified by map_uqiv_old_new_avail values and create a new matrix
  /// @details Sparse matrix version
  /// @param[in] mat_org Original sparse matrix
  /// @param[in] map_uqiv_old_new_avail Map from old column indices to new column indices
  /// @return New sparse matrix with columns reordered/extracted according to map
  /// @note Thread-safety: Yes (read-only operation)
  SpMatf get_unavaible_cols_removed_matrix(
    const SpMatf &mat_org, const std::map<Grid3d::Uqiv,Grid3d::Uqiv> &map_uqiv_old_new_avail);
    
  /// @brief Vector version of get_unavaible_cols_removed_matrix
  /// @param[in] vec_mat Vector of dense matrices
  /// @param[in] map_uqiv_old_new_avail Map from old column indices to new column indices
  /// @return Vector of new matrices with columns extracted according to map
  /// @note Thread-safety: Yes (read-only operation)
  std::vector<Eigen::MatrixXf> get_unavaible_cols_removed_matrix(
    const std::vector<Eigen::MatrixXf> &vec_mat,
    const std::map<Grid3d::Uqiv,Grid3d::Uqiv> &map_uqiv_old_new_avail);
  
  /// @brief Vector<SpMatf> version of get_unavaible_cols_removed_matrix
  /// @param[in] vec_mat_old Vector of sparse matrices
  /// @param[in] map_uqiv_old_new_avail Map from old column indices to new column indices
  /// @return Vector of new sparse matrices with columns extracted according to map
  /// @note Thread-safety: Yes (read-only operation)
  std::vector<SpMatf> get_unavaible_cols_removed_matrix(
    const std::vector<SpMatf> &vec_mat_old, const std::map<Grid3d::Uqiv,Grid3d::Uqiv> &map_uqiv_old_new_avail);

  /// @brief Extract columns specified by map_uqiv_old_new_avail values and create a new matrix
  /// @param[in] mat_org Original dense matrix
  /// @param[in] map_uqiv_old_new_avail Map from old column indices to new column indices
  /// @return New matrix with columns extracted according to map
  /// @note Thread-safety: Yes (read-only operation)
  Eigen::MatrixXf get_unavaible_cols_removed_matrix(
    const Eigen::MatrixXf &mat_org, const std::map<Grid3d::Uqiv,Grid3d::Uqiv> &map_uqiv_old_new_avail);

  /// @brief select rows from mat based on indices
  /// @param[in] mat Source matrix
  /// @param[in] keep_indices Row indices to select
  /// @return New matrix containing selected rows
  /// @throws std::runtime_error If any index is out of bounds
  /// @note Uses OpenMP parallel for loop. Data race-free (each thread writes to distinct output row).
  /// @note Complexity: O(rows * cols) where rows = keep_indices.size()
  Eigen::MatrixXf selectRows(
    const Eigen::MatrixXf& mat, const std::vector<int>& keep_indices);

  /// @brief select square block from mat based on indices
  /// @param[in] mat Source matrix (must be square or at least contain the indexed block)
  /// @param[in] keep_indices Row and column indices to select
  /// @return New square matrix containing selected rows and columns
  /// @note Uses OpenMP parallel for with collapse(2). Data race-free (each thread writes to distinct output element).
  /// @note Complexity: O(n^2) where n = keep_indices.size()
  /// @note No bounds checking performed; caller must ensure indices are valid
  Eigen::MatrixXf selectSquareBlock(
    const Eigen::MatrixXf& mat, const std::vector<int>& keep_indices );

  /// @brief get the unique, sorted Detid set from sorted_Detid_UqigAvail
  /// @param[in] sorted_Detid_UqigAvail Sorted map of detector ID to unique available indices
  /// @return Set of unique detector IDs
  /// @note Thread-safety: Yes (read-only operation)
  /// @note Complexity: O(n) where n = sorted_Detid_UqigAvail.size()
  std::set<Detid> get_set_Detid(
    const SortedDetidUqigSet& sorted_Detid_UqigAvail);

  /// @brief get the vector of unique available indices for a given detector ID
  /// @param[in] detid_in Target detector ID
  /// @param[in] sorted_Detid_UqigAvail Sorted map of detector ID to unique available indices
  /// @return std::vector<UqigAvail> Vector of UqigAvail corresponding to the target detector ID
  /// @note Thread-safety: Yes (read-only operation)
  /// @note Complexity: O(n) where n = sorted_Detid_UqigAvail.size() in worst case; early exit when detid > detid_in
  std::vector<UqigAvail> get_vec_uqig_avail(
      const Detid detid_in
    , const SortedDetidUqigSet& sorted_Detid_UqigAvail);

/// @brief get the vector of newly assigned UqigAvail_new by disabling by detid_in
/// @details Creates a mapping from new consecutive indices to original UqigAvail indices,
/// excluding all elements corresponding to detid_in.
/// @param[in] sorted_Detid_UqigAvail Sorted map of detector ID to unique available indices
/// @param[in] detid_in Target detector ID to disable
/// @param[in] uqig_disabled_avail_start Start index for new consecutive numbering (default: 0)
/// @return MapUqigAvailDisabled Map of newly assigned indices after disabling specified detector ID
/// @note Thread-safety: Yes (read-only operation)
/// @note Complexity: O(n) where n = sorted_Detid_UqigAvail.size()
  MapUqigAvailDisabled
    get_map_uqig_disabled_avail(
      const SortedDetidUqigSet& sorted_Detid_UqigAvail
    , const Detid detid_in
    , const UqigAvailDisabled uqig_disabled_avail_start=0 );

  /// @brief make and get the disabled mat_dNdD by detid
  /// @details Extracts rows from mat_dNdD_org corresponding to all detector elements except those
  /// belonging to detid_in. Useful for computing sensitivities with one detector disabled.
  /// @param[in] sorted_Detid_UqigAvail Sorted map of detector ID to unique available indices
  /// @param[in] mat_dNdD_org Original sensitivity matrix (rows = detector elements, cols = voxels)
  /// @param[in] detid_in Target detector ID to disable
  /// @param[in] uqig_disabled_avail_start Start index for assignment (default: 0)
  /// @return Matrix obtained by extracting rows from mat_dNdD based on the uqig of specified detector ID
  /// @note Must return a copy, not a reference
  /// @note Thread-safety: Uses OpenMP internally via selectRows
  Eigen::MatrixXf get_disable_mat_dNdD(
      const SortedDetidUqigSet& sorted_Detid_UqigAvail
    , const Eigen::MatrixXf& mat_dNdD_org
    , const Detid detid_in
    , const UqigAvailDisabled uqig_disabled_avail_start=0 );

  /// @brief make and get the disable_detid_vec_nmuon
  /// @details Extracts muon count vectors (observed and prior) corresponding to all detector elements
  /// except those belonging to detid_in.
  /// @param[in] sorted_Detid_UqigAvail Sorted map of detector ID to unique available indices
  /// @param[in] vecxf_nmuon_obs Observed muon counts per detector element
  /// @param[in] vecxf_nmuon_prior Prior muon counts per detector element
  /// @param[in] detid_in Target detector ID to disable
  /// @param[in] uqig_disabled_avail_start Start value for disabled uqig_avail (default: 0)
  /// @return NmuonVectors Structure containing vec_obs and vec_prior (subsetted vectors)
  /// @note Thread-safety: Yes (read-only operation)
  NmuonVectors get_disable_vec_nmuons(
      const SortedDetidUqigSet& sorted_Detid_UqigAvail
    , const Eigen::VectorXf &vecxf_nmuon_obs
    , const Eigen::VectorXf &vecxf_nmuon_prior
    , const Detid detid_in
    , const UqigAvailDisabled uqig_disabled_avail_start=0 );

  /// @brief make and get the disabled covariance matrix by detid
  /// @details Extracts square block from mat_cov_muon_org corresponding to all detector elements
  /// except those belonging to detid_in. Useful for computing uncertainties with one detector disabled.
  /// @param[in] sorted_Detid_UqigAvail Sorted map of detector ID to unique available indices
  /// @param[in] mat_cov_muon_org Original muon covariance matrix (square, size = total detector elements)
  /// @param[in] detid_in Target detector ID to disable
  /// @param[in] uqig_disabled_avail_start Start index for assignment (default: 0)
  /// @return Eigen::MatrixXf Covariance matrix excluding rows/columns corresponding to detid_in
  /// @note Must absolutely return a copy, not a reference
  /// @note Thread-safety: Uses OpenMP internally via selectSquareBlock
  Eigen::MatrixXf get_disable_mat_cov_muon(
      const SortedDetidUqigSet& sorted_Detid_UqigAvail
    , const Eigen::MatrixXf& mat_cov_muon_org
    , const Detid detid_in
    , const UqigAvailDisabled uqig_disabled_avail_start=0 );

};
