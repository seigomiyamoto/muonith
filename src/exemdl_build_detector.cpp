/// @file exemdl_build_detector.cpp
/// @brief Implementation of detector array construction module.
/// @details Implements functions for building DetectorPanelArray objects
///          from Grid2dPillar (2D terrain) and Grid3dVoxel (3D voxel) data.
///          Handles path/density length computation, signal/noise calculation,
///          and prior estimation with lower/center/upper density bounds.
///
/// @note All densities in kg/m^3. Coordinate system: z-up, right-handed.
/// @note Uses OpenMP for parallel computation via mp_* method calls.
#include "exemdl_build_detector.hpp"
#include "cls_FluxTable.hpp"
#include "cls_Grid2dPillar.hpp"
#include "ns_mylogger.hpp"
#include "ns_io_binary.hpp"
#include "ns_pathcalc.hpp"
#include <spdlog/spdlog.h>
#include <Eigen/Sparse>
#include "ns_iodir.hpp"

namespace fs = std::filesystem;

DetectorPanelArray exemdl::build_detector::build_arrdet_g2pil(
    const exemdl::load_parameters::AppParameters& app_params
  , const Grid2dPillar& g2pil_naive
  , const DetectorPanelArray& arrdet_template
  , const nlohmann::json& js_proj_dens
  , const FluxTable&      ft_prior )
{
  LOG_INFO("Building arrdet_g2pil_naive via DetectorPanelArray::create");
  const bool tf_apply_eff = app_params.prm_det.tf_apply_eff;
  DetectorPanelArray::PillarBuildParams prm(
      app_params.prm_path.tf_load_arrdet_g2pil
    , app_params.prm_path.tf_save_arrdet_g2pil
    , app_params.prm_bingroup.tf_run_1st_grouping
    , app_params.prm_bingroup.tf_run_auto_grouping
    , app_params.prm_path.path_arrdet_g2pil_bin
    , arrdet_template
    , app_params.prm_bingroup
    , g2pil_naive
    , app_params.ft_real
    , tf_apply_eff
    , app_params.prm_noise
  );
  DetectorPanelArray ret_arrdet_g2pil_naive = DetectorPanelArray::create(prm);
  ret_arrdet_g2pil_naive.set_name("arrdet_g2pil_naive");
  ret_arrdet_g2pil_naive.display_status();
  LOG_INFO("mp_calc_set_proj_density_all for {}, PL_thres={}, DL_thres={}"
    , ret_arrdet_g2pil_naive.get_name()
    , app_params.prm_bingroup.PL_thres
    , app_params.prm_bingroup.DL_thres
  );
  // Compute projected average density
  ret_arrdet_g2pil_naive.mp_calc_set_proj_density_all(
    app_params.prm_bingroup.PL_thres, app_params.prm_bingroup.DL_thres);

  // Also compute grouped density here
  LOG_INFO("Running mp_calc_set_proj_dens_grouped_all on arrdet_g2pil_naive (inside build_arrdet_g2pil)");
  ret_arrdet_g2pil_naive.mp_calc_set_proj_dens_grouped_all(
    js_proj_dens, ft_prior);
  return ret_arrdet_g2pil_naive;
}

std::tuple<DetectorPanelArray, std::vector<SpMatf>, Eigen::VectorXf>
  exemdl::build_detector::build_arrdet_g3vox(
      const exemdl::load_parameters::AppParameters& app_params
    , Grid3dVoxel& g3vox_input
    , const Grid2dPillar& g2pil_shell_upper
    , const Grid2dPillar& g2pil_shell_lower
    , const Grid2dPillar& g2pil_shell_lateral
    , const bool has_shell_upper
    , const bool has_shell_lower
    , const bool has_shell_lateral
    , const std::optional<DetectorPanelArray>& opt_arrdet_g2pil_to_copy_bimap
    , const DetectorPanelArray& arrdet_template)
{
  g3vox_input.set_name("g3vox_input");

  LOG_INFO("Building arrdet_g3vox_input via DetectorPanelArray::create_sprs");
  const bool tf_apply_eff = app_params.prm_det.tf_apply_eff;
  // Central efficiency (deterministic eff_cnt) gated by the SAME flag as the C_N diagonal
  // (tf_eff_cn_diag), so the forward count and the error band always turn on/off together.
  const bool tf_apply_eff_cnt =
      !app_params.vec_prm_nagainv.empty()
   && app_params.vec_prm_nagainv.front().get_tf_eff_cn_diag();
  constexpr unsigned uqivStart = 0;
  DetectorPanelArray::VoxelBuildParams prm(
      app_params.prm_path.tf_load_arrdet_g3vox
    , app_params.prm_path.tf_save_arrdet_g3vox
    , app_params.prm_bingroup.tf_run_1st_grouping
    , app_params.prm_bingroup.tf_run_auto_grouping
    , app_params.prm_path.path_arrdet_g3vox_bin
    , app_params.prm_path.path_vec_spmat_PL_bin
    , arrdet_template
    , app_params.prm_bingroup
    , g3vox_input
    , app_params.ft_real
    , app_params.prm_path
    , app_params.prm_g3vox
    , g2pil_shell_upper
    , g2pil_shell_lower
    , g2pil_shell_lateral
    , has_shell_upper
    , has_shell_lower
    , has_shell_lateral
    , uqivStart
    , tf_apply_eff
    , app_params.prm_noise
    , tf_apply_eff_cnt
  );
  // create_sprs returns 3-tuple: (arrdet, vec_spmat_PL, vecxf_non_rec_vox_PL)
  std::tuple<DetectorPanelArray, std::vector<SpMatf>, Eigen::VectorXf>
    retTuple = DetectorPanelArray::create_sprs(prm, opt_arrdet_g2pil_to_copy_bimap);
  auto& arrdet_g3vox_input = std::get<0>(retTuple);
  arrdet_g3vox_input.set_name("arrdet_g3vox_input");
  arrdet_g3vox_input.display_status();
  SLEEP_MSEC(500);

  // Compute signal/noise for each group
  LOG_INFO("Computing signal/noise_poi_group for arrdet_g3vox_input and setting vec_signal/noise_poi_group");
  arrdet_g3vox_input.mp_calc_vec_signal_noise_group_all();

  LOG_INFO("mp_calc_set_proj_density_all for {}, PL_thres={}, DL_thres={}"
    , arrdet_g3vox_input.get_name()
    , app_params.prm_bingroup.PL_thres
    , app_params.prm_bingroup.DL_thres
  );
  arrdet_g3vox_input.mp_calc_set_proj_density_all(
    app_params.prm_bingroup.PL_thres, app_params.prm_bingroup.DL_thres);
  if (app_params.prm_det.tf_out_txty_ascii) {
    arrdet_g3vox_input.out_txtyPL();
    arrdet_g3vox_input.out_txtyDens();
    arrdet_g3vox_input.out_txtySignal();
  } else {
    LOG_INFO("Skipped txty text output (DETECTOR_PARAMETER_LISTS.tf_out_txty_ascii=false)");
  }
  return retTuple;
}

std::array<double, 3>
  exemdl::build_detector::calc_volume_weighted_average_density_lower_center_upper(
    DetectorPanelArray& arrdet_g3vox_input
  , const nlohmann::json& js_proj_dens
  , const FluxTable& ft_prior
  , const bool tf_prior_error
  , const bool tf_eff
  , const bool tf_eff_independent)
{
  // Check whether to execute
  if(!js_proj_dens["tf_exec"].get<bool>()){
    LOG_INFO("js_proj_dens.tf_exec == false; skipping volume-weighted average density computation");
    return {
        Grid2dBinGroup::dens_group_lower_init
      , Grid2dBinGroup::dens_group_center_init
      , Grid2dBinGroup::dens_group_upper_init
    };
  }

  LOG_INFO("Computing projected average density for all DetectorPanels");
  auto start = time_now;
  // Pass efficiency-uncertainty flags so the projected-density error band adds
  // var_eff in quadrature to the Poisson term (gated by tf_eff_cn_diag).
  arrdet_g3vox_input.mp_calc_set_proj_dens_grouped_all(
    js_proj_dens, ft_prior, tf_eff, tf_eff_independent);
  auto end = time_now;
  std::string msg = "mp_calc_set_proj_dens_grouped_all done.";
  myapp::cast_time_msec(spdlog::level::debug,msg,start,end);

  LOG_INFO("Computing volume-weighted average density from all DetectorPanel projected densities");
  start = time_now;
  
  std::array<double, 3> arr_dens = {
    Grid2dBinGroup::dens_group_lower_init
  , Grid2dBinGroup::dens_group_center_init
  , Grid2dBinGroup::dens_group_upper_init
  };

  if (!tf_prior_error) {
    const double dens_center = arrdet_g3vox_input.calc_volume_weighted_average_density(
      "vol_wei_avr_dens_detail");
    arr_dens = {dens_center, dens_center, dens_center};
    msg = "calc_volume_weighted_average_density (center only) done.";
  } else {
    arr_dens = arrdet_g3vox_input.calc_volume_weighted_average_density_lower_center_upper();
    msg = "calc_volume_weighted_average_density_lower_center_upper done.";
  }

  end = time_now;
  myapp::cast_time_msec(spdlog::level::debug,msg,start,end);

  if (!tf_prior_error) {
    LOG_INFO("Computed volume-weighted average density (center only): {:.4E} kg/m^3", arr_dens[1]);
  } else {
    LOG_INFO("Computed volume-weighted average density [lower, center, upper]: {:.4E}, {:.4E}, {:.4E} kg/m^3"
      , arr_dens[0], arr_dens[1], arr_dens[2]
    );
  }

  const fs::path pathout = iodir::make_pathout("prior_error_volume_weighted_average_density.tmp");
  FILE *fout = myapp::get_fout(pathout);
  LOG_INFO("Writing result to text file: {}", pathout.string());
  fprintf(fout, "# Volume-weighted average density: lower, center, upper [kg/m^3]\n");
  fprintf(fout, "%.4E  %.4E  %.4E\n", arr_dens[0], arr_dens[1], arr_dens[2]);
  myapp::close(fout,pathout);

  return arr_dens;
}

exemdl::build_detector::BuildResult
  exemdl::build_detector::build_all( const BuildArgs& args )
{
  // Return value
  BuildResult ret;

  // Build base DetectorPanelArray
  LOG_INFO("Building base DetectorPanelArray");
  DetectorPanelArray retArrdet(args.app_params.prm_det);
  retArrdet.set_name("arrdet_template");

  // Build Grid2dPillar linked arrdet
  ret.arrdet_g2pil_naive = build_arrdet_g2pil(
      args.app_params
    , args.g2pil_naive
    , retArrdet
    , args.js_proj_dens
    , args.app_params.ft_prior
  );

  // Build Grid3dVoxel linked arrdet
  std::tuple<DetectorPanelArray, std::vector<SpMatf>, Eigen::VectorXf>
    tuple_res = build_arrdet_g3vox(
        args.app_params
      , args.g3vox_input
      , args.g2pil_shell_upper
      , args.g2pil_shell_lower
      , args.g2pil_shell_lateral
      , args.has_shell_upper
      , args.has_shell_lower
      , args.has_shell_lateral
      , ret.arrdet_g2pil_naive
      , retArrdet
    );

  ret.arrdet_g3vox_input        = std::move(std::get<0>(tuple_res));
  ret.vec_spmat_PL              = std::move(std::get<1>(tuple_res));
  // stash filtered-out voxel PL; final destination is ret.shell_pl.vecxf_non_rec_vox_PL
  Eigen::VectorXf vecxf_non_rec_vox_PL = std::move(std::get<2>(tuple_res));

  LOG_INFO("Computing volume-weighted average density (with bounds); volume obtained from arrdet_g2pil in create_sprs");
  // Efficiency-uncertainty flags: same gate as the C_N diagonal (tf_eff_cn_diag),
  // so the projected-density error band turns on together with the forward count.
  const bool tf_eff_band =
      !args.app_params.vec_prm_nagainv.empty()
   && args.app_params.vec_prm_nagainv.front().get_tf_eff_cn_diag();
  const bool tf_eff_band_independent =
      tf_eff_band
   && args.app_params.vec_prm_nagainv.front().get_tf_eff_cn_diag_independent();
  ret.avr_dens_lower_center_upper = calc_volume_weighted_average_density_lower_center_upper(
      ret.arrdet_g3vox_input
    , args.js_proj_dens
    , args.app_params.ft_prior
    , args.app_params.tf_run_inversion_prior_error
    , tf_eff_band
    , tf_eff_band_independent
  );
  LOG_INFO("avr_dens_lower_center_upper = {:.4E}, {:.4E}, {:.4E} kg/m^3"
    , ret.avr_dens_lower_center_upper[0]
    , ret.avr_dens_lower_center_upper[1]
    , ret.avr_dens_lower_center_upper[2]
  );
  if (args.app_params.prm_det.tf_out_g2bg_ascii) {
    LOG_INFO("Outputting arrdet_g3vox_input g2bg after projected density calculation");
    ret.arrdet_g3vox_input.out_g2bg_all();
  } else {
    LOG_INFO("Skipped g2bg text output (DETECTOR_PARAMETER_LISTS.tf_out_g2bg_ascii=false)");
  }

  double uniform_prior_density = Grid2dBinGroup::dens_group_center_init;
  // If avr_dens_lower_center_upper is initial value, set NagaInvParameters::uniform_prior_density
  if (ret.avr_dens_lower_center_upper == std::array<double, 3>{
      Grid2dBinGroup::dens_group_lower_init
    , Grid2dBinGroup::dens_group_center_init
    , Grid2dBinGroup::dens_group_upper_init })
  {
    LOG_INFO("avr_dens_lower_center_upper is initial value; setting NagaInvParameters::uniform_prior_density");
    if (args.app_params.vec_prm_nagainv.empty()) {
      LOG_WARN("vec_prm_nagainv is EMPTY; fallback to Grid2dBinGroup::dens_group_center_init={}",
               Grid2dBinGroup::dens_group_center_init);
      uniform_prior_density = Grid2dBinGroup::dens_group_center_init;
    } else {
      const NagaInvParameters& prm_nagainv0 = args.app_params.vec_prm_nagainv.front();
      uniform_prior_density = prm_nagainv0.get_uniform_prior_density();
    }
    LOG_INFO("Setting uniform prior density: {}", uniform_prior_density);
  }else{
    uniform_prior_density = ret.avr_dens_lower_center_upper[1];
    LOG_INFO("Using center value of avr_dens_lower_center_upper as uniform prior density: {}"
      , uniform_prior_density);
  }

  LOG_INFO("Building arrdet_g3vox_prior");

  // Resolve per-shell densities (fallback to uniform_prior_density)
  double shell_density_upper   = uniform_prior_density;
  double shell_density_lower   = uniform_prior_density;
  double shell_density_lateral = uniform_prior_density;
  if (!args.app_params.vec_prm_nagainv.empty()) {
    const auto& prm0 = args.app_params.vec_prm_nagainv.front();
    shell_density_upper   = prm0.get_shell_density_upper(uniform_prior_density);
    shell_density_lower   = prm0.get_shell_density_lower(uniform_prior_density);
    shell_density_lateral = prm0.get_shell_density_lateral(uniform_prior_density);
  }
  LOG_INFO("Shell densities: upper={}, lower={}, lateral={}", shell_density_upper, shell_density_lower, shell_density_lateral);

  ret.prior_info_all = exemdl::build_prior::build_prior_info_all(
      args.app_params
    , args.g2pil_shell_upper
    , args.g2pil_shell_lower
    , args.g2pil_shell_lateral
    , args.has_shell_upper
    , args.has_shell_lower
    , args.has_shell_lateral
    , args.g3vox_input
    , ret.avr_dens_lower_center_upper
    , shell_density_upper
    , shell_density_lower
    , shell_density_lateral
    , ret.arrdet_g3vox_input
    , ret.arrdet_g2pil_naive
  );

  return ret;
}
