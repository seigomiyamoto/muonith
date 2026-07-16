// src/ns_mix_error.cpp
#include "ns_mix_error.hpp"
#include "ns_mylogger.hpp"
#include "ns_type_definitions.hpp"
#include "ns_mymacro.hpp"

Eigen::VectorXf mix_error::calc_max_disabled_error(
    const int n_vecxf_dens
  , const std::vector<NagaInv::ReconstResult> &vec_reconst_disable_mode_res)
{
  constexpr float kEps = 1e-6f;
  Eigen::VectorXf vecxf_max_disabled_error = Eigen::VectorXf::Zero(n_vecxf_dens);
  for (const auto &reconst_res : vec_reconst_disable_mode_res) {
    // If vecxf_delta_dens_prior is (near) zero, warn and skip
    if (reconst_res.vecxf_delta_dens_prior.cwiseAbs().maxCoeff() <= kEps) {
      LOG_WARN("vecxf_delta_dens_prior is (near) zero, skipping.");
      continue;
    }
    vecxf_max_disabled_error = vecxf_max_disabled_error.cwiseMax(
      reconst_res.vecxf_diff_from_real.cwiseAbs());
  }
  // If all elements of vecxf_max_disabled_error are zero, throw error
  if (vecxf_max_disabled_error.isZero(0)) {
    LOG_ERROR("all elements of vecxf_max_disabled_error are zero, something is wrong.");
    THROW_ERROR("mix_error::calc_max_disabled_error: all elements of vecxf_max_disabled_error are zero, something is wrong.");
  }

  return vecxf_max_disabled_error;
}

Eigen::VectorXf mix_error::calc_mixed_error(
    const Eigen::VectorXf& vecxf_prior_error
  , const Eigen::VectorXf& vecxf_mat_cov_dens_dash
  , const Eigen::VectorXf& vecxf_disabled_error)
{
  // Validate input vector sizes
  const int n_prior = vecxf_prior_error.size();
  const int n_cov = vecxf_mat_cov_dens_dash.size();
  const int n_disabled = vecxf_disabled_error.size();

  if (n_prior != n_cov || n_prior != n_disabled) {
    THROW_ERROR("mix_error::calc_mixed_error: Input vector size mismatch. prior={}, cov={}, disabled={}", n_prior, n_cov, n_disabled);
  }

  Eigen::VectorXf vecxf_mixed_error = Eigen::VectorXf::Zero(n_prior);

  // Add squared prior error
  vecxf_mixed_error.array() += vecxf_prior_error.array().square();

  // Add squared statistical covariance error
  vecxf_mixed_error.array() += vecxf_mat_cov_dens_dash.array().square();

  // Add squared disabled mode error
  vecxf_mixed_error.array() += vecxf_disabled_error.array().square();

  // Take element-wise square root
  vecxf_mixed_error.array() = vecxf_mixed_error.array().sqrt();

  return vecxf_mixed_error;
}

