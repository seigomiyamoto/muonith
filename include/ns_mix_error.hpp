/// @file ns_mix_error.hpp
/// @brief Mixed error calculation utilities
/// @details
/// This module provides utilities for computing combined errors and uncertainties
/// from multiple independent error sources in muon tomography reconstruction.
///
/// The mixed error combines three error components using root of sum of squares (RSS):
/// - Prior error: Uncertainty in the prior density distribution
/// - Statistical covariance error: Error derived from measurement statistics
/// - Disabled mode error: Systematic error from reconstruction in disabled mode
///
/// ## Workflow
/// 1. Calculate maximum disabled error using calc_max_disabled_error()
/// 2. Combine all error sources using calc_mixed_error()
///
/// ## Thread safety
/// All functions are thread-safe for read-only access to input parameters.
/// No internal state is maintained.
///
/// @note All error vectors must have the same size when passed to calc_mixed_error().
#pragma once
#include <vector>

#include "ns_type_definitions.hpp"
#include "cls_NagaInv.hpp"

#include <Eigen/Dense>

/// @brief for type definitions
using namespace index_type_definitions;

/// @brief Namespace for mixed error calculation functions
/// @details
/// Provides utilities to compute combined errors from multiple independent sources
/// in muon tomography density reconstruction. Uses root of sum of squares (RSS) method
/// to combine prior, statistical, and systematic error components.
namespace mix_error {

  //=================================================================
  /// @name mix_error functions 
  /// @details Functions to calculate mixed error from three error sources:
  /// @{
  
  /// @brief Calculates the maximum disabled error from reconstruction results.
  /// @details Computes the element-wise maximum of absolute differences from prior density
  /// across all reconstruction results in disable mode.
  /// @param n_vecxf_dens Number of elements in the density vector
  /// @param vec_reconst_disable_mode_res Reconstruction results in disable mode
  /// @return Vector of maximum absolute errors for each density element
  /// @throws std::runtime_error If all elements of the result are zero, indicating an error condition
  Eigen::VectorXf calc_max_disabled_error( const int n_vecxf_dens
    , const std::vector<NagaInv::ReconstResult> &vec_reconst_disable_mode_res);

  /// @brief Calculates combined error from multiple error sources.
  /// @details Computes the mixed error by taking the root of sum of squares (RSS)
  /// of three error components: prior error, statistical covariance error, and disabled mode error.
  /// @param vecxf_prior_error Error distribution from prior density
  /// @param vecxf_mat_cov_dens_dash Density error derived from statistical covariance
  /// @param vecxf_disabled_error Density error from reconstruction results in disable mode
  /// @return Vector of combined errors computed as sqrt(prior^2 + cov^2 + disabled^2)
  /// @throws std::runtime_error If input vectors have incompatible sizes
  Eigen::VectorXf calc_mixed_error(
      const Eigen::VectorXf& vecxf_prior_error
    , const Eigen::VectorXf& vecxf_mat_cov_dens_dash
    , const Eigen::VectorXf& vecxf_disabled_error);


  ///@} ------------------------------------------------------------------


};
