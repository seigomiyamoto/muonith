// src/exemdl_calc_rec_error.cpp
#include "exemdl_calc_rec_error.hpp"
#include "ns_mymacro.hpp"
#include "ns_mix_error.hpp"
#include "spdlog_pch.hpp"

void exemdl::calc_rec_error::check_sizes_eq(const char* what, Eigen::Index a, Eigen::Index b) {
  if (a != b) THROW_ERROR("calc_rec_error::check_sizes_eq: {} size mismatch: {} vs {}", what, a, b);
}

// Core implementation shared by both overloads
static void write_disable_and_errors_impl(
    const NagaInvLooper& looper
  , const std::vector<NagaInv::ReconstResult>& vec_reconst_res_disabled
  , const exemdl::run_inversion::InversionResultsAll& invResAll
  , const exemdl::build_prior::PriorInfoAll& prior_info_all
  , const Grid3dVoxel& g3vox_merged_input
  , const std::string& output_prefix )
{
  // --- Visualization output (disable test) ---
  const Eigen::VectorXf& vecxf_dens_rec = invResAll.center.reconst_res.vecxf_dens_rec;
  LOG_INFO("NagaInv looper, output_disable_mode_all_as_cross_section");
  looper.output_disable_mode_all_as_cross_section(vec_reconst_res_disabled, vecxf_dens_rec);

  // --- Error vector calculation ---
  const Eigen::VectorXf vecxf_stat =
    invResAll.center.reconst_res.vecxf_diag_sqrt_cov_dens;
  const Eigen::VectorXf vecxf_prior_low =
    prior_info_all.lower.g3vox.get_vecxf_density();
  const Eigen::VectorXf vecxf_prior_upp =
    prior_info_all.upper.g3vox.get_vecxf_density();

  exemdl::calc_rec_error::check_sizes_eq("prior low/up", vecxf_prior_low.size(), vecxf_prior_upp.size());
  exemdl::calc_rec_error::check_sizes_eq("stat vs rec",  vecxf_stat.size(),       vecxf_dens_rec.size());

  Eigen::VectorXf vecxf_prior_err =
    (vecxf_prior_upp - vecxf_prior_low).array().abs() * 0.5f;
  exemdl::calc_rec_error::check_sizes_eq("prior err vs rec", vecxf_prior_err.size(), vecxf_dens_rec.size());

  const int n = static_cast<int>(vecxf_prior_err.size());
  const Eigen::VectorXf vecxf_disabled_err =
    mix_error::calc_max_disabled_error(n, vec_reconst_res_disabled);
  exemdl::calc_rec_error::check_sizes_eq("disabled vs rec", vecxf_disabled_err.size(), vecxf_dens_rec.size());

  const Eigen::VectorXf vecxf_mixed_err =
    mix_error::calc_mixed_error(vecxf_prior_err, vecxf_stat, vecxf_disabled_err);

  // --- Cross-section output ---
  auto& prm_zcross = invResAll.center.prm_zcross;
  auto prefixed = [&](const std::string& base) -> std::string {
    return output_prefix.empty() ? base : fmt::format("{}_{}", base, output_prefix);
  };
  LOG_INFO("output stat/prior/disabled/mixed error as cross section (prefix='{}')", output_prefix);
  g3vox_merged_input.write_density_to_cross_section(prefixed("g3vox_stat_error"),     vecxf_stat,         prm_zcross);
  g3vox_merged_input.write_density_to_cross_section(prefixed("g3vox_prior_error"),    vecxf_prior_err,    prm_zcross);
  g3vox_merged_input.write_density_to_cross_section(prefixed("g3vox_disabled_error"), vecxf_disabled_err, prm_zcross);
  g3vox_merged_input.write_density_to_cross_section(prefixed("g3vox_mixed_error"),    vecxf_mixed_err,    prm_zcross);
}

void exemdl::calc_rec_error::write_disable_and_errors(
    const NagaInvLooper& looper
  , const std::vector<NagaInv::ReconstResult>& vec_reconst_res_disabled
  , const exemdl::run_inversion::InversionResultsAll& invResAll
  , const exemdl::build_detector::BuildResult& detRes
  , const Grid3dVoxel& g3vox_merged_input
  , const std::string& output_prefix )
{
  write_disable_and_errors_impl(
    looper, vec_reconst_res_disabled, invResAll,
    detRes.prior_info_all, g3vox_merged_input, output_prefix);
}

void exemdl::calc_rec_error::write_disable_and_errors(
    const NagaInvLooper& looper
  , const std::vector<NagaInv::ReconstResult>& vec_reconst_res_disabled
  , const exemdl::run_inversion::InversionResultsAll& invResAll
  , const exemdl::build_prior::PriorInfoAll& prior_info_all
  , const Grid3dVoxel& g3vox_merged_input
  , const std::string& output_prefix )
{
  write_disable_and_errors_impl(
    looper, vec_reconst_res_disabled, invResAll,
    prior_info_all, g3vox_merged_input, output_prefix);
}
