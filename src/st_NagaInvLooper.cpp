/// @file st_NagaInvLooper.cpp
/// @brief Implementation of NagaInvLooper disable-mode reconstruction
#include "st_NagaInvLooper.hpp"
#include "ns_eigen_blas.hpp"

NagaInv::ReconstResult NagaInvLooper::exec_disable_mode(
  const Detid &detid_disabled, const Eigen::VectorXf &vecxf_dens_prior
  , const std::optional<Eigen::VectorXf>& input_vecxf_rec_all
  , const std::optional<Eigen::MatrixXf>& precomputed_mat_cov_dens_inv)
{
  LOG_INFO("detid_disabled={}", detid_disabled);
  // Extract matrices from nagainv_org and disable rows corresponding to detid_disabled
  Eigen::MatrixXf mat_dNdD =
    dic.get_disabled_mat_dNdD(detid_disabled, nagainv_org.get_mat_dNdD_ref());

  Eigen::MatrixXf mat_cov_muon
    = dic.get_disabled_mat_cov_muon(detid_disabled, nagainv_org.get_mat_cov_muon_ref());

  const auto [vec_nmuon_obs, vec_nmuon_prior]
     = dic.get_disabled_vec_nmuons(detid_disabled
        , nagainv_org.get_vecxf_nmuon_obs_ref(), nagainv_org.get_vecxf_nmuon_prior_ref());

  // Clone NagaInv with disabled detector data
  NagaInv nagainv = nagainv_org.clone_with_modified_data(
    &mat_dNdD, &mat_cov_muon, &vec_nmuon_obs, &vec_nmuon_prior);

  const bool tf_save_tmp_data_bin = false;
  const bool tf_verbose = false;
  // Run reconstruction with disabled detector
  NagaInv::ReconstResult reconst_res
    = nagainv.mp_reconst_density_float(
        prm_nagainv
      , vecxf_dens_prior
      , input_vecxf_rec_all
      , tf_save_tmp_data_bin
      , tf_verbose
      , precomputed_mat_cov_dens_inv);

  if (!reconst_res.is_valid) {
    LOG_ERROR("reconstruction failed.");
  }

  return reconst_res;
}

std::vector<NagaInv::ReconstResult> NagaInvLooper::exec_disable_mode_all(
    const Eigen::VectorXf &vecxf_dens_prior
  , const std::optional<Eigen::VectorXf>& input_vecxf_rec_all)
{
  std::vector<NagaInv::ReconstResult> vec_reconst_res;

  // Iterate over all detectors, disabling each one sequentially
  std::set<Detid> set_detid = dic.get_set_detid();
  if (set_detid.empty()) {
    LOG_WARN("detector set is empty, nothing to do.");
    return vec_reconst_res;
  }
  LOG_INFO("detector set size={}", set_detid.size());

  // Precompute mat_cov_dens inverse once (shared across all disabled-detector runs)
  LOG_INFO("Precomputing mat_cov_dens_inv (reused for all {} iterations)...", set_detid.size());
  const Eigen::MatrixXf mat_cov_dens_inv
    = eigen_blas::getInverseMatrixFloat(nagainv_org.get_mat_cov_dens_ref());

  vec_reconst_res.reserve(set_detid.size());

  for (const Detid &detid : set_detid) {
    LOG_INFO("disabling detid={} and reconstructing...", detid);
    NagaInv::ReconstResult reconst_res
      = exec_disable_mode(detid, vecxf_dens_prior, input_vecxf_rec_all,
          mat_cov_dens_inv);

    if (!reconst_res.is_valid) {
      LOG_ERROR("reconstruction failed for detid={}", detid);
      continue;
    }

    LOG_INFO("detid={} reconstruction completed.", detid);
    vec_reconst_res.push_back(std::move(reconst_res));
  }
  LOG_INFO("for '{}': finished.", name);

  return vec_reconst_res;
}

void NagaInvLooper::output_disable_mode_all_as_cross_section(
    const std::vector<NagaInv::ReconstResult> &vec_reconst_res
  , const std::optional<Eigen::VectorXf>& input_vecxf_rec_all ) const
{
  for(size_t idx=0; idx<vec_reconst_res.size(); idx++) {
    const NagaInv::ReconstResult &reconst_res = vec_reconst_res.at(idx);
    LOG_INFO("outputing reconst_res.vecxf_dens_rec for idx={}", idx);
    std::string detid_str = fmt::format("{:02}", idx);
    std::string prefix = name + "_disable_" + detid_str;
    g3vox_org.write_density_to_cross_section(
      prefix , reconst_res.vecxf_dens_rec, prm_zcross);

    LOG_INFO("outputing reconst_res.vecxf_delta_dens_prior for idx={}", idx);
    g3vox_org.write_density_to_cross_section(
      prefix + "_delta_prior", reconst_res.vecxf_delta_dens_prior, prm_zcross);

    LOG_INFO("outputing diff density from all_detector reconst result for idx={}", idx);
    if (input_vecxf_rec_all.has_value()) {
      g3vox_org.write_density_to_cross_section(
        prefix + "_diff_from_all", reconst_res.vecxf_diff_from_real, prm_zcross);
    } else {
      LOG_WARN("input_vecxf_rec_all is not provided, skipping diff output for idx={}", idx);
    }
  }
}