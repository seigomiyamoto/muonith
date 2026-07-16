/// @file st_NmuonVectors.hpp
/// @brief Observed and prior muon count vectors
/// @details
/// This file defines a simple data structure for holding observed and prior
/// muon count vectors as Eigen::VectorXf. It is used to pass both vectors
/// together throughout the inversion workflow.
///
/// **Thread Safety:** Yes. This is a plain data structure with no internal state.
/// Users are responsible for synchronizing concurrent access to shared instances.
///
/// **Memory Layout:** Standard Eigen::VectorXf storage (column-major by default).
#pragma once
#include <Eigen/Dense>

/// @brief Container for observed and prior muon count vectors
/// @details
/// This struct holds two muon count vectors:
/// - `vec_obs`: Observed muon counts from detector measurements
/// - `vec_prior`: Prior (expected) muon counts from forward model or initial guess
///
/// Both vectors should have the same size, corresponding to the number of
/// detector bins or measurement points.
///
/// **Usage Example:**
/// @code
/// NmuonVectors nmuon;
/// nmuon.vec_obs = Eigen::VectorXf::Zero(100);
/// nmuon.vec_prior = Eigen::VectorXf::Ones(100) * 50.0f;
/// @endcode
///
/// @note Thread-safety: This is a plain data struct. Concurrent access to
/// the same instance requires external synchronization.
struct NmuonVectors {
  Eigen::VectorXf vec_obs;   ///< Observed muon counts from detector measurements
  Eigen::VectorXf vec_prior; ///< Prior (expected) muon counts from forward model
};
