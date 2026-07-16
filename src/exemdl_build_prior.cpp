/// @file exemdl_build_prior.cpp
/// @brief Implementation of prior estimation module.
/// @details Implements functions for building prior density information
///          (lower/center/upper bounds) used in the inversion stage.
///          Extracted from exemdl_build_detector.cpp.
///
/// @note All densities in kg/m^3. Coordinate system: z-up, right-handed.
/// @note Uses OpenMP for parallel computation via mp_* method calls.
#include "exemdl_build_prior.hpp"
#include "cls_FluxTable.hpp"
#include "cls_Grid2dPillar.hpp"
#include "ns_mylogger.hpp"
#include "ns_io_binary.hpp"
#include "ns_pathcalc.hpp"
#include <spdlog/spdlog.h>
#include <Eigen/Sparse>
#include "ns_iodir.hpp"

exemdl::build_prior::PriorInfo
  exemdl::build_prior::build_arrdet_g3vox_prior(
      const std::string& suffix
    , const exemdl::load_parameters::AppParameters& app_params
    , const Grid2dPillar& g2pil_shell_upper
    , const Grid2dPillar& g2pil_shell_lower
    , const Grid2dPillar& g2pil_shell_lateral
    , const bool has_shell_upper
    , const bool has_shell_lower
    , const bool has_shell_lateral
    , const Grid3dVoxel&  g3vox_input
    , const double voxel_density
    , const double shell_density_upper
    , const double shell_density_lower
    , const double shell_density_lateral
    , const DetectorPanelArray& arrdet_g3vox_input
    , const DetectorPanelArray& arrdet_g2pil_naive)
{
  LOG_INFO("Preparing g2pil_shell_upper_copy (density={})", shell_density_upper);
  Grid2dPillar g2pil_shell_upper_copy(g2pil_shell_upper);
  g2pil_shell_upper_copy.set_name("g2pil_shell_upper_copy");
  g2pil_shell_upper_copy.set_uniform_density(shell_density_upper);

  Grid2dPillar g2pil_shell_lower_copy(g2pil_shell_lower);
  g2pil_shell_lower_copy.set_name("g2pil_shell_lower_copy");
  g2pil_shell_lower_copy.set_uniform_density(shell_density_lower);

  Grid2dPillar g2pil_shell_lateral_copy(g2pil_shell_lateral);
  g2pil_shell_lateral_copy.set_name("g2pil_shell_lateral_copy");
  g2pil_shell_lateral_copy.set_uniform_density(shell_density_lateral);

  LOG_INFO("Preparing g3vox_prior (density={})", voxel_density);
  Grid3dVoxel g3vox_prior(g3vox_input);
  g3vox_prior.set_name("g3vox_prior"+suffix);
  g3vox_prior.set_uniform_density_ixiyiz(voxel_density);

  LOG_INFO("Initializing arrdet_g3vox_prior");
  DetectorPanelArray arrdet_g3vox_prior(arrdet_g3vox_input);
  arrdet_g3vox_prior.set_name("arrdet_g3vox_prior"+suffix);
  arrdet_g3vox_prior.mp_init_PL_DL_pene_sig_noi_all();
  arrdet_g3vox_prior.display_status();

  // Check if detector positions are within the grid range
  arrdet_g3vox_prior.check_detector_positions_in_grid_xy(
    g3vox_prior.get_xmin(), g3vox_prior.get_xmax(),
    g3vox_prior.get_ymin(), g3vox_prior.get_ymax());

  LOG_INFO("Computing path/density lengths from g2pil_shell for arrdet_g3vox_prior");
  const double BL_max = app_params.prm_path.BL_max;
  constexpr double eps = 1.0e-6;
  if (app_params.prm_path.tf_add_shell) {
    LOG_INFO("---------------------------------------------------");
    if (has_shell_upper) {
      LOG_INFO("Adding path/density lengths from g2pil_shell_upper_copy");
      pathcalc::g2pil::mp_add_PLDL(arrdet_g3vox_prior, g2pil_shell_upper_copy, BL_max, eps);
    }
    if (has_shell_lower) {
      LOG_INFO("Adding path/density lengths from g2pil_shell_lower_copy");
      pathcalc::g2pil::mp_add_PLDL(arrdet_g3vox_prior, g2pil_shell_lower_copy, BL_max, eps);
    }
    if (has_shell_lateral) {
      LOG_INFO("Adding path/density lengths from g2pil_shell_lateral_copy");
      pathcalc::g2pil::mp_add_PLDL(arrdet_g3vox_prior, g2pil_shell_lateral_copy, BL_max, eps);
    }
  } else {
    LOG_INFO("tf_add_shell==false; skipping g2pil_shell path/density lengths");
  }

  // Add g3vox PL/DL to arrdet_g3vox
  LOG_INFO("---------------------------------------------------");
  LOG_INFO("instance_name={}, adding path/density lengths from g3vox"
    , arrdet_g3vox_prior.get_name());

  pathcalc::g3vox::mp_add_PLDL(arrdet_g3vox_prior, g3vox_prior, BL_max, eps);
  arrdet_g3vox_prior.set_tf_calc_PL_DL(true);
  LOG_INFO("---------------------------------------------------");

  // Compute and set signal for arrdet_g3vox
  LOG_INFO("Computing and setting signal for instance_name={}", arrdet_g3vox_prior.get_name() );
  arrdet_g3vox_prior.mp_calc_set_peneflux_signal_from_DL(
    app_params.ft_prior, app_params.prm_det.tf_apply_eff
  );

  LOG_INFO("Setting BinGroup map for arrdet_g3vox...");
  arrdet_g3vox_prior.set_dic_bimap_all_from(arrdet_g2pil_naive);
  arrdet_g3vox_prior.display_status();
  SLEEP_MSEC(500);

  LOG_INFO("Copying detector element signal/noise to Grid2dBinGroup vec_vec_signal/noise");
  arrdet_g3vox_prior.mp_copy_signal_noise_to_g2bg_all(app_params.prm_bingroup);
  arrdet_g3vox_prior.display_status();

  LOG_INFO("Computing signal/noise_poi_group for each group in arrdet_g3vox and setting vec_signal/noise_poi_group");
  arrdet_g3vox_prior.mp_calc_vec_signal_noise_group_all();
  SLEEP_MSEC(500);

  // regrouping via set_dic_bimap_all_from reset the group efficiencies to -1
  LOG_INFO("Set group efficiency (simple mean of element efficiencies) of {}", arrdet_g3vox_prior.get_name());
  arrdet_g3vox_prior.calc_set_eff_group();

  LOG_INFO("mp_calc_set_proj_density_all for {}, PL_thres={}, DL_thres={}"
    , arrdet_g3vox_prior.get_name()
    , app_params.prm_bingroup.PL_thres
    , app_params.prm_bingroup.DL_thres
  );
  arrdet_g3vox_prior.mp_calc_set_proj_density_all(
    app_params.prm_bingroup.PL_thres, app_params.prm_bingroup.DL_thres);
  if (app_params.prm_det.tf_out_txty_ascii) {
    arrdet_g3vox_prior.out_txtyPL();
    arrdet_g3vox_prior.out_txtyDens();
    arrdet_g3vox_prior.out_txtySignal();
    SLEEP_MSEC(500);
  } else {
    LOG_INFO("Skipped txty text output (DETECTOR_PARAMETER_LISTS.tf_out_txty_ascii=false)");
  }

  if (app_params.prm_det.tf_out_g2bg_ascii) {
    LOG_INFO("Outputting g2bg for arrdet_g3vox");
    arrdet_g3vox_prior.out_g2bg_all();
    SLEEP_MSEC(500);
  } else {
    LOG_INFO("Skipped g2bg text output (DETECTOR_PARAMETER_LISTS.tf_out_g2bg_ascii=false)");
  }

  PriorInfo ret;
  ret.vec_vecxf_DL = arrdet_g3vox_prior.get_vec_vecxf_DL();
  ret.vecxf_nmuon = arrdet_g3vox_prior.get_vecxf_nmuon_all();
  ret.g3vox = std::move(g3vox_prior);

  return ret;
}

exemdl::build_prior::PriorInfoAll
  exemdl::build_prior::build_prior_info_all(
      const exemdl::load_parameters::AppParameters& app_params
    , const Grid2dPillar& g2pil_shell_upper
    , const Grid2dPillar& g2pil_shell_lower
    , const Grid2dPillar& g2pil_shell_lateral
    , const bool has_shell_upper
    , const bool has_shell_lower
    , const bool has_shell_lateral
    , const Grid3dVoxel&  g3vox_input
    , const std::array<double,3>& avr_dens_lower_center_upper
    , const double shell_density_upper
    , const double shell_density_lower
    , const double shell_density_lateral
    , const DetectorPanelArray& arrdet_g3vox_input
    , const DetectorPanelArray& arrdet_g2pil_naive)
{
  LOG_INFO("Building PriorInfo prior_lower");
  PriorInfo prior_lower = build_arrdet_g3vox_prior(
      "_lower"
    , app_params
    , g2pil_shell_upper
    , g2pil_shell_lower
    , g2pil_shell_lateral
    , has_shell_upper
    , has_shell_lower
    , has_shell_lateral
    , g3vox_input
    , avr_dens_lower_center_upper[0]
    , shell_density_upper
    , shell_density_lower
    , shell_density_lateral
    , arrdet_g3vox_input
    , arrdet_g2pil_naive
  );
  prior_lower.g3vox.set_name("g3vox_prior_lower");

  LOG_INFO("Building PriorInfo prior_center");
  PriorInfo prior_center = build_arrdet_g3vox_prior(
      "_center"
    , app_params
    , g2pil_shell_upper
    , g2pil_shell_lower
    , g2pil_shell_lateral
    , has_shell_upper
    , has_shell_lower
    , has_shell_lateral
    , g3vox_input
    , avr_dens_lower_center_upper[1]
    , shell_density_upper
    , shell_density_lower
    , shell_density_lateral
    , arrdet_g3vox_input
    , arrdet_g2pil_naive
  );
  prior_center.g3vox.set_name("g3vox_prior_center");

  LOG_INFO("Building PriorInfo prior_upper");
  PriorInfo prior_upper = build_arrdet_g3vox_prior(
      "_upper"
    , app_params
    , g2pil_shell_upper
    , g2pil_shell_lower
    , g2pil_shell_lateral
    , has_shell_upper
    , has_shell_lower
    , has_shell_lateral
    , g3vox_input
    , avr_dens_lower_center_upper[2]
    , shell_density_upper
    , shell_density_lower
    , shell_density_lateral
    , arrdet_g3vox_input
    , arrdet_g2pil_naive
  );
  prior_upper.g3vox.set_name("g3vox_prior_upper");

  PriorInfoAll prior_info_all;
  prior_info_all.lower = std::move(prior_lower);
  prior_info_all.center = std::move(prior_center);
  prior_info_all.upper = std::move(prior_upper);

  return prior_info_all;
}

//----------------------------------------------------------------------
// PriorInfo save/load
//----------------------------------------------------------------------

void exemdl::build_prior::PriorInfo::save(std::ofstream& ofs) const
{
  io_binary::write_vec_vecxf(ofs, vec_vecxf_DL);
  io_binary::write_vecxf_stream(ofs, vecxf_nmuon);
  g3vox.save(ofs);
}

void exemdl::build_prior::PriorInfo::load(std::ifstream& ifs)
{
  vec_vecxf_DL = io_binary::read_vec_vecxf(ifs);
  vecxf_nmuon = io_binary::read_vecxf_stream(ifs);
  g3vox.load(ifs);
}

//----------------------------------------------------------------------
// PriorInfoAll save/load
//----------------------------------------------------------------------

void exemdl::build_prior::PriorInfoAll::save(std::ofstream& ofs) const
{
  lower.save(ofs);
  center.save(ofs);
  upper.save(ofs);
}

void exemdl::build_prior::PriorInfoAll::load(std::ifstream& ifs)
{
  lower.load(ifs);
  center.load(ifs);
  upper.load(ifs);
}

//----------------------------------------------------------------------
// ShellPL save/load
//----------------------------------------------------------------------

void exemdl::build_prior::ShellPL::save(const std::filesystem::path& dir) const
{
  // Persist the four density-independent PL vectors to disk so compute_prior
  // (rebuild_prior_from_shell_PL) can rerun with different densities
  // without redoing ray tracing.
  std::filesystem::create_directories(dir);
  io_binary::out_vecxf_bin(dir / "PL_shell_upper.bin", vecxf_upper_PL);
  io_binary::out_vecxf_bin(dir / "PL_shell_lower.bin", vecxf_lower_PL);
  io_binary::out_vecxf_bin(dir / "PL_shell_lateral.bin", vecxf_lateral_PL);
  // vecxf_non_rec_vox_PL covers the PL contribution of mountain-interior voxels
  // dropped by the n_hit_det filter (see ShellPL doc-comment). Without it,
  // rebuild_prior_from_shell_PL under-estimates prior DL.
  io_binary::out_vecxf_bin(dir / "PL_non_rec_voxel.bin", vecxf_non_rec_vox_PL);
  LOG_INFO("ShellPL::save: saved 4 binary files to {}", dir.string());
}

exemdl::build_prior::ShellPL
  exemdl::build_prior::ShellPL::load(const std::filesystem::path& dir)
{
  ShellPL ret;
  ret.vecxf_upper_PL   = io_binary::read_vecxf_bin(dir / "PL_shell_upper.bin");
  ret.vecxf_lower_PL   = io_binary::read_vecxf_bin(dir / "PL_shell_lower.bin");
  ret.vecxf_lateral_PL = io_binary::read_vecxf_bin(dir / "PL_shell_lateral.bin");
  // Older checkpoints (written before vecxf_non_rec_vox_PL was introduced) lack
  // PL_non_rec_voxel.bin. Treat that as an empty contribution rather than
  // failing the load, so old runs can still be inspected. New runs will
  // produce the file from build_all and rebuild_prior_from_shell_PL adds
  // it only when size() > 0.
  const std::filesystem::path nrv_path = dir / "PL_non_rec_voxel.bin";
  if (std::filesystem::exists(nrv_path)) {
    ret.vecxf_non_rec_vox_PL = io_binary::read_vecxf_bin(nrv_path);
    LOG_INFO("ShellPL::load: loaded 4 binary files from {}", dir.string());
  } else {
    LOG_WARN("ShellPL::load: PL_non_rec_voxel.bin not found in {}; "
             "vecxf_non_rec_vox_PL left empty (legacy checkpoint)", dir.string());
  }
  return ret;
}

//----------------------------------------------------------------------
// build_shell_PL
//----------------------------------------------------------------------

exemdl::build_prior::ShellPL
  exemdl::build_prior::build_shell_PL(
      const exemdl::load_parameters::AppParameters& app_params
    , const DetectorPanelArray& arrdet_g3vox_input
    , const Grid2dPillar& g2pil_shell_upper
    , const Grid2dPillar& g2pil_shell_lower
    , const Grid2dPillar& g2pil_shell_lateral
    , bool has_shell_upper
    , bool has_shell_lower
    , bool has_shell_lateral)
{
  const double BL_max = app_params.prm_path.BL_max;
  constexpr double eps = 1.0e-6;

  ShellPL ret;

  if (app_params.prm_path.tf_add_shell) {
    if (has_shell_upper) {
      LOG_INFO("build_shell_PL: computing PL for shell_upper");
      ret.vecxf_upper_PL = pathcalc::g2pil::mp_calc_PL(
        arrdet_g3vox_input, g2pil_shell_upper, BL_max, eps);
    }
    if (has_shell_lower) {
      LOG_INFO("build_shell_PL: computing PL for shell_lower");
      ret.vecxf_lower_PL = pathcalc::g2pil::mp_calc_PL(
        arrdet_g3vox_input, g2pil_shell_lower, BL_max, eps);
    }
    if (has_shell_lateral) {
      LOG_INFO("build_shell_PL: computing PL for shell_lateral");
      ret.vecxf_lateral_PL = pathcalc::g2pil::mp_calc_PL(
        arrdet_g3vox_input, g2pil_shell_lateral, BL_max, eps);
    }
  } else {
    LOG_INFO("build_shell_PL: tf_add_shell==false; all shell PL set to empty");
  }

  return ret;
}

//----------------------------------------------------------------------
// rebuild_prior_from_shell_PL
//----------------------------------------------------------------------
exemdl::build_prior::PriorInfo
  exemdl::build_prior::rebuild_prior_from_shell_PL(
      const std::string& suffix
    , const exemdl::load_parameters::AppParameters& app_params
    , const ShellPL& shell_pl
    , const std::vector<SpMatf>& vec_spmat_PL
    , const Grid3dVoxel& g3vox_input
    , double voxel_density
    , double shell_density_upper
    , double shell_density_lower
    , double shell_density_lateral
    , const DetectorPanelArray& arrdet_g3vox_input
    , const DetectorPanelArray& arrdet_g2pil_naive)
{
  LOG_INFO("rebuild_prior_from_shell_PL: suffix={}", suffix);

  // 1. Create g3vox_prior from g3vox_input with uniform density
  Grid3dVoxel g3vox_prior(g3vox_input);
  g3vox_prior.set_name("g3vox_prior" + suffix);
  g3vox_prior.set_uniform_density_ixiyiz(voxel_density);

  // 2. Initialize arrdet_g3vox_prior
  DetectorPanelArray arrdet_g3vox_prior(arrdet_g3vox_input);
  arrdet_g3vox_prior.set_name("arrdet_g3vox_prior" + suffix);
  arrdet_g3vox_prior.mp_init_PL_DL_pene_sig_noi_all();
  arrdet_g3vox_prior.display_status();

  arrdet_g3vox_prior.check_detector_positions_in_grid_xy(
    g3vox_prior.get_xmin(), g3vox_prior.get_xmax(),
    g3vox_prior.get_ymin(), g3vox_prior.get_ymax());

  // 3. Reconstruct PL/DL from shell PL + sparse matrices (no ray tracing)
  //    Phase A: Eigen vector ops per panel, packed into single VectorXf
  //    Phase B: Set PL/DL via Uqid with OpenMP
  LOG_INFO("rebuild_prior_from_shell_PL: reconstructing PL/DL from cached shell PL");
  const float f_shell_upper_dens = static_cast<float>(shell_density_upper);
  const float f_shell_lower_dens = static_cast<float>(shell_density_lower);
  const float f_shell_lateral_dens = static_cast<float>(shell_density_lateral);
  const float f_vox_dens = static_cast<float>(voxel_density);

  // Phase A: Compute PL/DL into single contiguous VectorXf (one alloc each)
  const int n_det = arrdet_g3vox_prior.get_n_det();
  int total_ele = 0;
  std::vector<int> offsets(n_det);
  for (Detid detid = 0; detid < n_det; ++detid) {
    offsets[detid] = total_ele;
    total_ele += arrdet_g3vox_prior.getDetectorPanel(detid).get_n_element();
  }

  Eigen::VectorXf vecxf_pl_all = Eigen::VectorXf::Zero(total_ele);
  Eigen::VectorXf vecxf_dl_all = Eigen::VectorXf::Zero(total_ele);
  for (Detid detid = 0; detid < n_det; ++detid) {
    const int n_ele = arrdet_g3vox_prior.getDetectorPanel(detid).get_n_element();
    const int off = offsets[detid];
    // .segment(off, n_ele): writable view of [off, off+n_ele); writes hit vecxf_*_all in-place.
    auto seg_pl = vecxf_pl_all.segment(off, n_ele);
    auto seg_dl = vecxf_dl_all.segment(off, n_ele);

    if (shell_pl.vecxf_upper_PL.size() > 0) {
      // .segment here is a read-only slice over the same [off, off+n_ele) window.
      seg_pl += shell_pl.vecxf_upper_PL.segment(off, n_ele);
      seg_dl += shell_pl.vecxf_upper_PL.segment(off, n_ele) * f_shell_upper_dens;
    }
    if (shell_pl.vecxf_lower_PL.size() > 0) {
      seg_pl += shell_pl.vecxf_lower_PL.segment(off, n_ele);                          // read slice
      seg_dl += shell_pl.vecxf_lower_PL.segment(off, n_ele) * f_shell_lower_dens;     // read slice
    }
    if (shell_pl.vecxf_lateral_PL.size() > 0) {
      seg_pl += shell_pl.vecxf_lateral_PL.segment(off, n_ele);                        // read slice
      seg_dl += shell_pl.vecxf_lateral_PL.segment(off, n_ele) * f_shell_lateral_dens; // read slice
    }

    if (detid < static_cast<int>(vec_spmat_PL.size())) {
      const SpMatf& spmatf_pl = vec_spmat_PL[detid];
      // Step 1: build [1, 1, ..., 1] of length = voxel count (= spmatf_pl.cols()).
      Eigen::VectorXf vecxf_ones_voxel = Eigen::VectorXf::Ones(spmatf_pl.cols());
      // Step 2: spmatf_pl * vecxf_ones_voxel sums each row -> per-element total PL.
      Eigen::VectorXf vecxf_pl_vox = spmatf_pl * vecxf_ones_voxel;
      seg_pl += vecxf_pl_vox;                 // PL: add geometric path length
      seg_dl += vecxf_pl_vox * f_vox_dens;    // DL: PL weighted by voxel density
    }

    // Add PL of voxels dropped by n_hit_det filter (missing from spmatf_pl above).
    // Same voxel density (mountain interior is uniform).
    if (shell_pl.vecxf_non_rec_vox_PL.size() > 0) {
      seg_pl += shell_pl.vecxf_non_rec_vox_PL.segment(off, n_ele);                  // read slice
      seg_dl += shell_pl.vecxf_non_rec_vox_PL.segment(off, n_ele) * f_vox_dens;     // read slice
    }
  }

  // Phase B: Set PL/DL via Uqid (OpenMP-friendly flat iteration)
  const auto& uqid_mgr = arrdet_g3vox_prior.get_dic().getUqidMgr();
  const auto vec_uqid = uqid_mgr.get_vecUqid_all();
  const int n_uqid = static_cast<int>(vec_uqid.size());
  #pragma omp parallel for schedule(static)
  for (int i = 0; i < n_uqid; ++i) {
    const Uqid uqid = vec_uqid[i];
    const auto info = uqid_mgr.getInfo(uqid);
    const Detid detid = info.detid;
    const int nbinx = arrdet_g3vox_prior.getDetectorPanel(detid).get_nbinx();
    const int flat = offsets[detid] + info.ixiy[1] * nbinx + info.ixiy[0];

    DetectorElement& ele = arrdet_g3vox_prior.callDetectorElement(uqid);
    ele.add_PL(static_cast<double>(vecxf_pl_all(flat)));
    ele.add_DL(static_cast<double>(vecxf_dl_all(flat)));
  }
  arrdet_g3vox_prior.set_tf_calc_PL_DL(true);
  LOG_INFO("---------------------------------------------------");

  // 4. Post-processing (same as build_arrdet_g3vox_prior)
  LOG_INFO("Computing and setting signal for instance_name={}", arrdet_g3vox_prior.get_name());
  arrdet_g3vox_prior.mp_calc_set_peneflux_signal_from_DL(
    app_params.ft_prior, app_params.prm_det.tf_apply_eff);

  LOG_INFO("Setting BinGroup map for arrdet_g3vox...");
  arrdet_g3vox_prior.set_dic_bimap_all_from(arrdet_g2pil_naive);
  arrdet_g3vox_prior.display_status();
  SLEEP_MSEC(500);

  LOG_INFO("Copying detector element signal/noise to Grid2dBinGroup vec_vec_signal/noise");
  arrdet_g3vox_prior.mp_copy_signal_noise_to_g2bg_all(app_params.prm_bingroup);
  arrdet_g3vox_prior.display_status();

  LOG_INFO("Computing signal/noise_poi_group for each group in arrdet_g3vox and setting vec_signal/noise_poi_group");
  arrdet_g3vox_prior.mp_calc_vec_signal_noise_group_all();
  SLEEP_MSEC(500);

  // regrouping via set_dic_bimap_all_from reset the group efficiencies to -1
  LOG_INFO("Set group efficiency (simple mean of element efficiencies) of {}", arrdet_g3vox_prior.get_name());
  arrdet_g3vox_prior.calc_set_eff_group();

  LOG_INFO("mp_calc_set_proj_density_all for {}, PL_thres={}, DL_thres={}"
    , arrdet_g3vox_prior.get_name()
    , app_params.prm_bingroup.PL_thres
    , app_params.prm_bingroup.DL_thres);
  arrdet_g3vox_prior.mp_calc_set_proj_density_all(
    app_params.prm_bingroup.PL_thres, app_params.prm_bingroup.DL_thres);
  if (app_params.prm_det.tf_out_txty_ascii) {
    arrdet_g3vox_prior.out_txtyPL();
    arrdet_g3vox_prior.out_txtyDens();
    arrdet_g3vox_prior.out_txtySignal();
    SLEEP_MSEC(500);
  } else {
    LOG_INFO("Skipped txty text output (DETECTOR_PARAMETER_LISTS.tf_out_txty_ascii=false)");
  }

  if (app_params.prm_det.tf_out_g2bg_ascii) {
    LOG_INFO("Outputting g2bg for arrdet_g3vox");
    arrdet_g3vox_prior.out_g2bg_all();
    SLEEP_MSEC(500);
  } else {
    LOG_INFO("Skipped g2bg text output (DETECTOR_PARAMETER_LISTS.tf_out_g2bg_ascii=false)");
  }

  // PriorInfo keeps only DL/nmuon/g3vox, so the array itself would be lost here.
  // Persist it as a DetectorPanelArray binary, the only form plot_det_arrdet.py reads.
  if (app_params.prm_det.tf_save_arrdet_prior) {
    const fs::path path_arrdet_prior =
      iodir::make_pathout("det/" + arrdet_g3vox_prior.get_name() + ".bin");
    LOG_INFO("Saving arrdet_g3vox_prior to {}", path_arrdet_prior.string());
    arrdet_g3vox_prior.save(path_arrdet_prior);
  } else {
    LOG_INFO("Skipped arrdet_g3vox_prior binary "
             "(DETECTOR_PARAMETER_LISTS.tf_save_arrdet_prior=false)");
  }

  PriorInfo ret;
  ret.vec_vecxf_DL = arrdet_g3vox_prior.get_vec_vecxf_DL();
  ret.vecxf_nmuon = arrdet_g3vox_prior.get_vecxf_nmuon_all();
  ret.g3vox = std::move(g3vox_prior);

  return ret;
}
