/// @file exemdl_calc_matrix.cpp
/// @brief Implementation of observation matrix construction for muon tomography

#include "exemdl_calc_matrix.hpp"
#include "cls_MatrixBuildParameters.hpp"
#include <spdlog/spdlog.h>
#include "cls_PathLengthParameters.hpp"
#include "ns_calc_dNdD.hpp"

Eigen::MatrixXf exemdl::calc_matrix::build_matrix(
  BuildArgs& args, const std::vector<Eigen::VectorXf>& vec_vecxf_DL_prior_in)
{
  // Central efficiency (deterministic eff_cnt) for A, gated by the SAME flag as the C_N diagonal
  // (tf_eff_cn_diag), so the matrix A and the forward count share the same central efficiency.
  const bool tf_eff_cn_diag =
      !args.app_params.vec_prm_nagainv.empty()
   && args.app_params.vec_prm_nagainv.front().get_tf_eff_cn_diag();

  // Build observation matrix (dN/dD) from path lengths and prior density-length vectors
  pathcalc::MatrixBuildParameters prm_mat(
      "prm_pathmatrix"
    , args.app_params.prm_path.tf_load_bin_obs_mat_dNdD
    , args.app_params.prm_path.tf_save_bin_obs_mat_dNdD
    , args.app_params.prm_path.path_bin_obs_mat_dNdD
    , args.res_built_det.arrdet_g3vox_input
    , args.g3vox_input
    , args.res_built_det.vec_spmat_PL
    , vec_vecxf_DL_prior_in
    , args.app_params.ft_prior
    , args.app_params.prm_det.tf_apply_eff
    , tf_eff_cn_diag
  );

  LOG_INFO("Creating grouped observation matrix (mat_dNdD_grouped)");
  Eigen::MatrixXf mat_dNdD_grouped
    = calc_dNdD::create_grouped_mat_dNdD_alldet_sprs(prm_mat);

  return mat_dNdD_grouped;
}

exemdl::calc_matrix::BuildResult exemdl::calc_matrix::build_all(BuildArgs& args)
{
  BuildResult res;

  LOG_INFO("Building observation matrix for lower-bound prior");
  res.mat_dNdD_grouped_lower
    = build_matrix(args, args.res_built_det.prior_info_all.lower.vec_vecxf_DL);

  LOG_INFO("Building observation matrix for center prior");
  res.mat_dNdD_grouped_center
    = build_matrix(args, args.res_built_det.prior_info_all.center.vec_vecxf_DL);

  LOG_INFO("Building observation matrix for upper-bound prior");
  res.mat_dNdD_grouped_upper
    = build_matrix(args, args.res_built_det.prior_info_all.upper.vec_vecxf_DL);

  return res;
}