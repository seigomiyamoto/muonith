/// @file cls_DetectorPanelArray.cpp
/// @brief Implementation of DetectorPanelArray class
/// @details
/// Implements the DetectorPanelArray container which manages multiple detector panels,
/// provides unified indexing, and handles signal/noise calculations across all panels.
///
/// Key implementation details:
/// - Uses OpenMP for parallel operations (build_vec_panel, grouping, signal calculations)
/// - Maintains DetectorIndexContainer (dic_) for unified element access
/// - Supports both cuboid and voxel geometry backends
/// - Binary serialization includes architecture compatibility checks
///
/// @note Thread Safety: OpenMP-parallelized methods (mp_*) are internally parallel but
///       not safe for concurrent external calls. Sequential external usage required.
#include <cassert>
#include <iostream>
#include <fstream>
#include <filesystem>
namespace fs = std::filesystem;

#include "cls_DetectorPanelArray.hpp"

#include "ns_myapp.hpp"
#include "ns_mymacro.hpp"
#include "ns_mylogger.hpp"
#include "ns_pathcalc.hpp"
#include "spdlog_pch.hpp"
#include "ns_mylogger.hpp"
#include "ns_eigen_blas.hpp"
#include "ns_detector_indexing.hpp"
#include "ns_geom_util.hpp"

#include "ns_iodir.hpp"

#include <nlohmann/json.hpp>
using json = nlohmann::json;
//###########################################################################################
//###########################################################################################
// class DetectorPanelArray
//###########################################################################################
//###########################################################################################

bool DetectorPanelArray::operator!=(const DetectorPanelArray &other) const
{
  #ifdef NODEBUG
    if (n_all_element != other.n_all_element) return true;
    if (vec_parameter_file_path != other.vec_parameter_file_path) return true;
    if (vec_panel != other.vec_panel) return true;
    if (dic_ != other.dic_) return true;
    if (flags_ != other.flags_) return true;
    if (set_uqiv != other.set_uqiv) return true;
  #else
    if (n_all_element != other.n_all_element) { LOG_WARN("DetectorPanelArray: n_all_element differs"); return true; }
    if (vec_parameter_file_path != other.vec_parameter_file_path) { LOG_WARN("DetectorPanelArray: vec_parameter_file_path differs"); return true; }
    if (vec_panel != other.vec_panel) { LOG_WARN("DetectorPanelArray: vec_panel differs"); return true; }
    if (dic_ != other.dic_) { LOG_WARN("DetectorPanelArray: dic_ differs"); return true; }
    if (flags_ != other.flags_) { LOG_WARN("DetectorPanelArray: flags_ differs"); return true; }
    if (set_uqiv != other.set_uqiv) { LOG_WARN("DetectorPanelArray: set_uqiv differs"); return true; }
  #endif
  return false;
}

// copy constructor
DetectorPanelArray::DetectorPanelArray(const DetectorPanelArray &org)
  : name(org.name)
  , n_all_element(org.n_all_element)
  , vec_parameter_file_path(org.vec_parameter_file_path)
  , vec_panel(org.vec_panel)
  , dic_(org.dic_)
  , flags_(org.flags_)
  , set_uqiv(org.set_uqiv)
{
  // Update each DetectorPanel to point to the new dic_ after copying
  for (DetectorPanel& panel : vec_panel) {
    panel.set_dic_ptr(&dic_);
    panel.call_g2bg().set_dic_ptr(&dic_);
  }
}

// move assignment operator
DetectorPanelArray& DetectorPanelArray::operator=(DetectorPanelArray&& other) noexcept
{
  if (this != &other) {
    name                    = std::move(other.name);
    n_all_element           = other.n_all_element;
    vec_parameter_file_path = std::move(other.vec_parameter_file_path);
    vec_panel               = std::move(other.vec_panel);
    dic_                    = std::move(other.dic_);
    flags_                  = std::move(other.flags_);
    set_uqiv                = std::move(other.set_uqiv);
    // Update pdic_ in vec_panel to point to the new dic_ after move
    for (DetectorPanel& panel : vec_panel) panel.set_dic_ptr(&dic_);
  }
  return *this;
}

// move constructor
DetectorPanelArray::DetectorPanelArray(DetectorPanelArray&& other) noexcept
  : name(std::move(other.name))
  , n_all_element(other.n_all_element)
  , vec_parameter_file_path(std::move(other.vec_parameter_file_path))
  , vec_panel(std::move(other.vec_panel))
  , dic_(std::move(other.dic_))
  , flags_(std::move(other.flags_))
  , set_uqiv(std::move(other.set_uqiv))
{
  // Update pdic_ in vec_panel to point to the new dic_ after move
  for (DetectorPanel& panel : vec_panel) panel.set_dic_ptr(&dic_);
}

// Factory function that calls DetectorPanelArray constructor based on conditions
DetectorPanelArray DetectorPanelArray::create(
  const DetectorPanelArray::PillarBuildParams &prm_arrdet_cub )
{
  std::chrono::system_clock::time_point start, end;
  start = time_now;
  const bool tf_file_found = fs::exists(prm_arrdet_cub.path_arrdet_bin);

  if( prm_arrdet_cub.tf_load_arrdet_g2pil == true && tf_file_found == true ){
    // When tf_load_arrdet_g2pil == true and file exists,
    // construct data from DetectorPanelArray binary file.
    LOG_INFO("tf_load_arrdet_g2pil = true.");
    return DetectorPanelArray(prm_arrdet_cub.path_arrdet_bin); // File-loading constructor
  }

  // When tf_load_arrdet_g2pil == true but file not found, continue with false path
  if (prm_arrdet_cub.tf_load_arrdet_g2pil == true && tf_file_found == false) {
    LOG_WARN(
      "tf_load_arrdet_g2pil = true but file {} not found. Proceeding with calculation."
      , prm_arrdet_cub.path_arrdet_bin.string());
  }

  // The following handles cases when file does not exist or tf_load_arrdet_g2pil == false
  LOG_INFO("tf_load_arrdet_g2pil = false or file not found.");

  // Copy the built arrdet_built to arrdet_g2pil.
  LOG_INFO("Copy arrdet_built to arrdet_g2pil.");
  DetectorPanelArray arrdet_g2pil(prm_arrdet_cub.arrdet_built);
  LOG_DEBUG("DetectorPanelArray arrdet_g2pil(prm_arrdet_cub.arrdet_built");

  // Set Grid2dBinGroup::Parameters to arrdet_g2pil.
  LOG_INFO("Set Grid2dBinGroup::Parameters to arrdet_g2pil.");
  arrdet_g2pil.set_parameters_all(prm_arrdet_cub.prm_bingroup);
  LOG_DEBUG("(prm_arrdet_cub.prm_bingroup)");
    
  // Check if detector positions are within the grid range
  arrdet_g2pil.check_detector_positions_in_grid_xy(
    prm_arrdet_cub.g2pil.get_xmin(), prm_arrdet_cub.g2pil.get_xmax(),
    prm_arrdet_cub.g2pil.get_ymin(), prm_arrdet_cub.g2pil.get_ymax());

  // Step-by-step timing within create(PillarBuildParams)
  std::chrono::system_clock::time_point t_step, t_step_end;

  // Calculate path length and density length of arrdet_g2pil from g2pil.
  LOG_INFO("Calculate path length and density length of arrdet_g2pil from g2pil.");
  // LOG_DEBUG("(arrdet_g2pil,prm_arrdet_cub.g2pil)");
  // pathcalc::g2pil::mp_add_PLDL(arrdet_g2pil,prm_arrdet_cub.g2pil);

  LOG_INFO("adding PLDL with vec_tf_in_PL");
  const int ixiy_dist_thres = 1;
  pathcalc::g2pil::mp_add_PLDL_with_vec_tf_in_PL(
    arrdet_g2pil, prm_arrdet_cub.g2pil, ixiy_dist_thres);

  // Set signal to arrdet_g2pil.
  t_step = time_now;
  LOG_INFO("Set signal to arrdet_g2pil.");
  arrdet_g2pil.mp_calc_set_peneflux_signal_from_DL(
    prm_arrdet_cub.ft, prm_arrdet_cub.tf_apply_eff);
  t_step_end = time_now;
  myapp::cast_time_msec(spdlog::level::info,"[g2pil] mp_calc_set_peneflux_signal_from_DL",t_step,t_step_end);

  // give the noise distribution from prm_noise to arrdet_g2pil
  t_step = time_now;
  LOG_INFO("arrdet_g2pil.mp_set_noise_all");
  const double DL_thres = prm_arrdet_cub.prm_bingroup.DL_thres;
  arrdet_g2pil.mp_set_noise_all(prm_arrdet_cub.prm_noise, DL_thres);
  t_step_end = time_now;
  myapp::cast_time_msec(spdlog::level::info,"[g2pil] mp_set_noise_all",t_step,t_step_end);

  // set detector element signal to Grid2dBinGroup vec_vec_signal;
  t_step = time_now;
  LOG_INFO("set detector element signal/noise/is_avail to Grid2dBinGroup vec_vec_signal/noise/is_avail");
  arrdet_g2pil.mp_copy_signal_noise_to_g2bg_all(prm_arrdet_cub.prm_bingroup);
  t_step_end = time_now;
  myapp::cast_time_msec(spdlog::level::info,"[g2pil] mp_copy_signal_noise_to_g2bg_all",t_step,t_step_end);

  // From vec_vec_signal, automatically set appropriate bin divisions and group minimum angle elements.
  if( prm_arrdet_cub.tf_run_1st_grouping == false ){
    LOG_INFO("(g2pil) tf_run_1st_grouping == false. So auto grouping is not executed.");
    SLEEP_MSEC(500);
  } // tf_1st_grouping==true below
  else{
    // assign 1st igroup for all detectors
    t_step = time_now;
    LOG_INFO("assign 1st igroup for all detectors");
    arrdet_g2pil.mp_assign_1st_igroup_all(
      prm_arrdet_cub.prm_bingroup.tf_run_1st_grouping, prm_arrdet_cub.prm_bingroup.igroup_start
    , prm_arrdet_cub.prm_bingroup.nx_div_init, prm_arrdet_cub.prm_bingroup.ny_div_init
    , prm_arrdet_cub.prm_bingroup.signal_noise_group_trig );
    t_step_end = time_now;
    myapp::cast_time_msec(spdlog::level::info,"[g2pil] mp_assign_1st_igroup_all",t_step,t_step_end);
  }

  if( prm_arrdet_cub.tf_run_auto_grouping == false ){
    LOG_INFO("(g2pil) tf_run_auto_grouping == false. So auto grouping is not executed.");

    if( prm_arrdet_cub.prm_bingroup.n_detector_grouping_manual > 0 ){
      // Manual grouping from user-supplied bin-list files (tx/ty rectangles).
      t_step = time_now;
      arrdet_g2pil.mp_grouping_by_bin_list_all(prm_arrdet_cub.prm_bingroup);
      t_step_end = time_now;
      myapp::cast_time_msec(spdlog::level::info,"[g2pil] mp_grouping_by_bin_list_all",t_step,t_step_end);
    }

    // Call build_index_container().
    if ( !(prm_arrdet_cub.tf_load_arrdet_g2pil && tf_file_found) ) {
      arrdet_g2pil.build_index_container(); // Build dic_
    }

    if( prm_arrdet_cub.prm_bingroup.n_detector_grouping_manual > 0 ){
      LOG_INFO("Allocate memory for vec_signal/noise_poi_group in each Group of arrdet_g2pil.");
      arrdet_g2pil.mp_allocate_vec_value_group_all();
      LOG_INFO("Copy each signal_group value of the generated BinGroup to vec_signal_group.");
      arrdet_g2pil.mp_calc_vec_signal_noise_group_all();
      LOG_INFO("Set group efficiency (simple mean of element efficiencies) of arrdet_g2pil.");
      arrdet_g2pil.calc_set_eff_group();
    }

    SLEEP_MSEC(500);
  }else{
    // tf_run_1st_grouping && tf_run_auto_grouping == true
    t_step = time_now;
    LOG_INFO("(g2pil) tf_run_1st_grouping && tf_run_auto_grouping == true.");
    LOG_INFO("(g2pil) From vec_vec_signal, automatically set appropriate bin divisions and group minimum angle elements.");
    arrdet_g2pil.mp_auto_grouping_by_signal_noise_group_alldet(prm_arrdet_cub.prm_bingroup);
    t_step_end = time_now;
    myapp::cast_time_msec(spdlog::level::info,"[g2pil] mp_auto_grouping_by_signal_noise_group_alldet",t_step,t_step_end);

    // Call build_index_container().
    t_step = time_now;
    if ( !(prm_arrdet_cub.tf_load_arrdet_g2pil && tf_file_found) ) {
      arrdet_g2pil.build_index_container(); // Build dic_
    }
    t_step_end = time_now;
    myapp::cast_time_msec(spdlog::level::info,"[g2pil] build_index_container",t_step,t_step_end);

    t_step = time_now;
    LOG_INFO("Allocate memory for vec_signal/noise_poi_group in each Group of arrdet_g2pil.");
    arrdet_g2pil.mp_allocate_vec_value_group_all();
    t_step_end = time_now;
    myapp::cast_time_msec(spdlog::level::info,"[g2pil] mp_allocate_vec_value_group_all",t_step,t_step_end);

    t_step = time_now;
    LOG_INFO("Copy each signal_group value of the generated BinGroup to vec_signal_group.");
    arrdet_g2pil.mp_calc_vec_signal_noise_group_all();
    t_step_end = time_now;
    myapp::cast_time_msec(spdlog::level::info,"[g2pil] mp_calc_vec_signal_noise_group_all",t_step,t_step_end);

    t_step = time_now;
    LOG_INFO("Set group efficiency (simple mean of element efficiencies) of arrdet_g2pil.");
    arrdet_g2pil.calc_set_eff_group();
    t_step_end = time_now;
    myapp::cast_time_msec(spdlog::level::info,"[g2pil] calc_set_eff_group",t_step,t_step_end);
  }

  // Check that no group is smaller than ixlen_min x iylen_min (covers all grouping routes).
  arrdet_g2pil.mp_check_group_ixiylen_min_all(
    prm_arrdet_cub.prm_bingroup.ixlen_min, prm_arrdet_cub.prm_bingroup.iylen_min );

  // save arrdet_g2pil to binary file if prm_arrdet_cub.tf_save_arrdet_g2pil==true
  if( prm_arrdet_cub.tf_save_arrdet_g2pil == true ){
    LOG_INFO("save arrdet_g2pil to binary file.");
    arrdet_g2pil.save(prm_arrdet_cub.path_arrdet_bin);
    LOG_DEBUG("(prm_arrdet_cub.path_arrdet_bin)");
  }

  end = time_now;
  std::string msg = "DetectorPanelArray::create return arrdet_g2pil";
  myapp::cast_time_msec(spdlog::level::info, msg, start, end);

  arrdet_g2pil.set_name("arrdet_g2pil_create");
  arrdet_g2pil.display_status();

  return arrdet_g2pil;
}

/// @brief Factory function that calls DetectorPanelArray constructor based on conditions
std::tuple<DetectorPanelArray, std::vector<SpMatf>>
  DetectorPanelArray::create(
    const DetectorPanelArray::VoxelBuildParams &prm_arrdet_vox )
{
  const bool tf_file_found = fs::exists(prm_arrdet_vox.path_arrdet_bin);
  // When tf_load_arrdet_g3vox == true, construct data from binary file if file exists
  if (prm_arrdet_vox.tf_load_arrdet_g3vox == true && tf_file_found == true) {
    LOG_INFO("tf_load_arrdet_g3vox = true. ");
    LOG_INFO(
      "So load DetectorPanelArray object from binary file {}",
      prm_arrdet_vox.path_arrdet_bin.string()
    );

    DetectorPanelArray arrdet_g3vox(prm_arrdet_vox.path_arrdet_bin); // File-loading constructor

    LOG_INFO(
      "also read vec_spmat_PL from binary file {}",
      prm_arrdet_vox.path_vec_spmat_PL_bin.string()
    );
    std::vector<SpMatf> vec_spmat_PL
      = io_binary::read_vec_spmatf(prm_arrdet_vox.path_vec_spmat_PL_bin);
    LOG_INFO("vec_spmat_PL.size()={}", vec_spmat_PL.size());

    LOG_INFO("Finally, return arrdet_g3vox, vec_spmat_PL");
    return std::make_tuple(arrdet_g3vox, vec_spmat_PL);
  }

  // When tf_load_arrdet_g3vox == true but file not found, continue with false path
  if (prm_arrdet_vox.tf_load_arrdet_g3vox == true && tf_file_found == false) {
    LOG_WARN("tf_load_arrdet_g3vox = true.\n");
    LOG_WARN("but file {} not found. So program will calculate arrdet_g3vox and vec_spmat_PL."
    , prm_arrdet_vox.path_arrdet_bin.string());
  }
  //
  // Process when tf_load_arrdet_g3vox == false or file not found
  //

  LOG_INFO("tf_load_arrdet_g3vox = false or file not found. calculate arrdet_g3vox and vec_spmat_PL.");

  std::chrono::system_clock::time_point t_create_start, t_create_end;
  t_create_start = time_now;

  LOG_INFO("====================================================================================");
  LOG_INFO("Create arrdet_g3vox from arrdet_built.");
  DetectorPanelArray arrdet_g3vox(prm_arrdet_vox.arrdet_built);

  // When tf_add_shell==true, add path length and density length of g2pil_shell (upper and lower)
  if(  prm_arrdet_vox.prm_pathcalc.tf_add_shell == true ){
    LOG_INFO("----------------------------------------");
    if( prm_arrdet_vox.has_shell_upper ){
      LOG_INFO("Add path length and density length of g2pil_shell_upper");
      pathcalc::g2pil::mp_add_PLDL(
        arrdet_g3vox,
        prm_arrdet_vox.g2pil_shell_upper
      );
    }
    if( prm_arrdet_vox.has_shell_lower ){
      LOG_INFO("Add path length and density length of g2pil_shell_lower");
      pathcalc::g2pil::mp_add_PLDL(
        arrdet_g3vox,
        prm_arrdet_vox.g2pil_shell_lower
      );
    }
    if( prm_arrdet_vox.has_shell_lateral ){
      LOG_INFO("Add path length and density length of g2pil_shell_lateral");
      pathcalc::g2pil::mp_add_PLDL(
        arrdet_g3vox,
        prm_arrdet_vox.g2pil_shell_lateral
      );
    }
  }

  LOG_INFO(" Calculate peneflux and signal of arrdet_g3vox.");
  arrdet_g3vox.mp_calc_set_peneflux_signal_from_DL(
    prm_arrdet_vox.ft, prm_arrdet_vox.tf_apply_eff, prm_arrdet_vox.tf_apply_eff_cnt
  );

  LOG_INFO("----------------------------------------");
  LOG_INFO("Calculate vec_mat");

  // Display prm.reference_matPL_sparse, prm.epsilon_matPL_sparse
  LOG_DEBUG("prm.reference_matPL_sparse={:.1E}, prm.epsilon_matPL_sparse={:.1E}"
    , prm_arrdet_vox.prm_pathcalc.reference_matPL_sparse
    , prm_arrdet_vox.prm_pathcalc.epsilon_matPL_sparse);

  // Set signal to arrdet_g3vox and create observation matrix for each detector.
  LOG_INFO("Set signal to arrdet_g3vox and create observation matrix for each detector");
  std::vector<SpMatf> vec_spmat_PL_allvox =
    pathcalc::g3vox::mp_make_mat_PL(
      arrdet_g3vox,
      prm_arrdet_vox.g3vox,
      prm_arrdet_vox.prm_pathcalc
    );

  // give the noise distribution from prm_noise to arrdet_g2pil
  LOG_INFO("arrdet_g3vox.mp_set_noise_all");
  const double DL_thres = prm_arrdet_vox.prm_bingroup.DL_thres;
  arrdet_g3vox.mp_set_noise_all(prm_arrdet_vox.prm_noise, DL_thres);

  LOG_INFO("Execute get_g2bg().calc_vec_signal_group() for all detectors in arrdet_g3vox,");
  LOG_INFO("Calculate signal_group for each Group and set to vec_signal_group.");
  arrdet_g3vox.mp_calc_vec_signal_noise_group_all();

  const int n_vec = vec_spmat_PL_allvox.size();
  LOG_DEBUG(" matrix size check: vec_spmat_PL_allvox.size()={}", vec_spmat_PL_allvox.size());
  for (int i = 0; i < n_vec; i++) {
    const int n_row = vec_spmat_PL_allvox.at(i).rows();
    const int n_col = vec_spmat_PL_allvox.at(i).cols();
    if (n_row == 0) THROW_ERROR("n_row == 0");
    if (n_col == 0) THROW_ERROR("n_col == 0");
    LOG_DEBUG("vec_spmat_PL_allvox.at({}).rows()={}, cols={}", i, n_row, n_col);
  }

  LOG_INFO("Assign new unique_index to g3vox with n_hit_det filter.");
  const bool tf_exist = true;
  const std::map<Grid3d::Uqiv,Grid3d::Uqiv> map_uqiv_old_new_avail
    = prm_arrdet_vox.g3vox.re_assign_uqiv_by_nhit_det(
        tf_exist,
        prm_arrdet_vox.prm_g3vox,
        prm_arrdet_vox.uqiv_start
      );
  LOG_INFO(" Remove columns for voxels that did not meet filter conditions for n_hit_det.");
  std::vector<SpMatf> vec_spmat_PL
    = detector_indexing::get_unavaible_cols_removed_matrix(vec_spmat_PL_allvox, map_uqiv_old_new_avail);

  LOG_DEBUG("size of vec_spmat_PL_allvox BEFORE:");
  eigen_blas::disp_vec_spmatf(vec_spmat_PL_allvox);
  LOG_DEBUG("size of vec_spmat_PL AFTER:");
  eigen_blas::disp_vec_spmatf(vec_spmat_PL);
  SLEEP_MSEC(500);

  LOG_INFO(" Release memory of vec_spmat_PL_allvox.");
  myapp::clearVectorMemory(vec_spmat_PL_allvox);

  LOG_INFO("Set Grid2dBinGroup::Parameters to arrdet_g3vox_input.");
  arrdet_g3vox.set_parameters_all(prm_arrdet_vox.prm_bingroup);
  LOG_DEBUG("(prm_arrdet_vox.prm_bingroup)");

  LOG_INFO("set all detector element signal to Grid2dBinGroup vec_vec_signal");
  arrdet_g3vox.mp_copy_signal_noise_to_g2bg_all(prm_arrdet_vox.prm_bingroup);
  LOG_DEBUG("(prm_arrdet_vox.prm_bingroup)");


  // From vec_vec_signal, automatically set appropriate bin divisions and group minimum angle elements.
  if( prm_arrdet_vox.tf_run_1st_grouping == false ){
    LOG_INFO("(g3vox) tf_run_1st_grouping == false. So auto grouping is not executed.");
    SLEEP_MSEC(500);
  } // tf_1st_grouping==true below
  else{
    // assign 1st igroup for all detectors
    LOG_INFO("assign 1st igroup for all detectors");
    arrdet_g3vox.mp_assign_1st_igroup_all(
      prm_arrdet_vox.prm_bingroup.tf_run_1st_grouping, prm_arrdet_vox.prm_bingroup.igroup_start
    , prm_arrdet_vox.prm_bingroup.nx_div_init, prm_arrdet_vox.prm_bingroup.ny_div_init
    , prm_arrdet_vox.prm_bingroup.signal_noise_group_trig );
  }
  if( prm_arrdet_vox.tf_run_auto_grouping == false ){
    LOG_INFO("(g3vox) tf_run_auto_grouping == false. So auto grouping is not executed.");
    if( prm_arrdet_vox.prm_bingroup.n_detector_grouping_manual > 0 ){
      // Manual grouping from user-supplied bin-list files (tx/ty rectangles).
      arrdet_g3vox.mp_grouping_by_bin_list_all(prm_arrdet_vox.prm_bingroup);
    }
    SLEEP_MSEC(500);
  }else{
    // tf_run_1st_grouping && tf_run_auto_grouping == true
    LOG_INFO("(g3vox) tf_run_1st_grouping && tf_run_auto_grouping == true.");
    LOG_INFO("(g3vox) From vec_vec_signal, automatically set appropriate bin divisions and group minimum angle elements.");
    arrdet_g3vox.mp_auto_grouping_by_signal_noise_group_alldet(prm_arrdet_vox.prm_bingroup);

    LOG_INFO("Allocate memory for vec_signal/noise_poi_group in each Group of arrdet_g3vox.");
    arrdet_g3vox.mp_allocate_vec_value_group_all();

    LOG_INFO("Set group efficiency (simple mean of element efficiencies) of arrdet_g3vox.");
    arrdet_g3vox.calc_set_eff_group();
  }

  if (prm_arrdet_vox.tf_save_arrdet_g3vox == true) {
    LOG_INFO("save arrdet_g3vox to binary file.");
    arrdet_g3vox.save(prm_arrdet_vox.path_arrdet_bin);
    LOG_DEBUG("(prm_arrdet_vox.path_arrdet_bin)");

    LOG_INFO("save vec_spmat_PL to binary file.");
    io_binary::write_vec_spmatf(prm_arrdet_vox.path_vec_spmat_PL_bin, vec_spmat_PL);
    LOG_DEBUG("(prm_arrdet_vox.path_vec_spmat_PL_bin, vec_spmat_PL)");
  }

  LOG_INFO("get the set of new uqiv from map_uqiv_old_new_avail");
  arrdet_g3vox.set_set_uqiv( tuple_int::get_set_value(map_uqiv_old_new_avail) );


  // Call build_index_container().
  if (!(prm_arrdet_vox.tf_load_arrdet_g3vox && tf_file_found)) {
    arrdet_g3vox.build_index_container(); // Build dic_
  }

  if( prm_arrdet_vox.tf_run_1st_grouping == true ){
    LOG_INFO("Calculate vec_signal/noise_poi_group in each Group of arrdet_g3vox.");
    arrdet_g3vox.mp_calc_vec_signal_noise_group_all();
  }

  arrdet_g3vox.set_name("arrdet_g3vox_create");
  arrdet_g3vox.display_status();

  t_create_end = time_now;
  myapp::cast_time_msec(spdlog::level::info,"DetectorPanelArray::create(VoxelBuildParams) finished",t_create_start,t_create_end);

  LOG_INFO("return arrdet_g3vox, vec_spmat_PL");
  return std::make_tuple(arrdet_g3vox, vec_spmat_PL);
}

/// @brief Factory function that calls DetectorPanelArray constructor based on conditions
std::tuple<DetectorPanelArray, std::vector<SpMatf>, Eigen::VectorXf>
  DetectorPanelArray::create_sprs(
    const DetectorPanelArray::VoxelBuildParams &prm_arrdet_vox
  , const std::optional<DetectorPanelArray>& opt_arrdet_built )
{
  const bool tf_file_found = fs::exists(prm_arrdet_vox.path_arrdet_bin);
  // When tf_load_arrdet_g3vox == true, construct data from binary file if file exists
  if (prm_arrdet_vox.tf_load_arrdet_g3vox == true && tf_file_found == true) {
    LOG_INFO("tf_load_arrdet_g3vox = true. ");
    LOG_INFO(
      "So load DetectorPanelArray object from binary file {}",
      prm_arrdet_vox.path_arrdet_bin.string()
    );

    DetectorPanelArray arrdet_g3vox(prm_arrdet_vox.path_arrdet_bin); // File-loading constructor

    LOG_INFO(
      "also read vec_spmat_PL from binary file {}",
      prm_arrdet_vox.path_vec_spmat_PL_bin.string()
    );
    std::vector<SpMatf> vec_spmat_PL
      = io_binary::read_vec_spmatf(prm_arrdet_vox.path_vec_spmat_PL_bin);
    LOG_INFO("vec_spmat_PL.size()={}", vec_spmat_PL.size());

    LOG_INFO("Finally, return arrdet_g3vox, vec_spmat_PL");
    // cache-load path: vecxf_non_rec_vox_PL is empty here; obtain via ShellPL::load
    return std::make_tuple(arrdet_g3vox, vec_spmat_PL, Eigen::VectorXf());
  }

  // When tf_load_arrdet_g3vox == true but file not found, continue with false path
  if (prm_arrdet_vox.tf_load_arrdet_g3vox == true && tf_file_found == false) {
    LOG_WARN("tf_load_arrdet_g3vox = true.\n");
    LOG_WARN("but file {} not found. So program will calculate arrdet_g3vox and vec_spmat_PL."
    , prm_arrdet_vox.path_arrdet_bin.string());
  }
  //
  // Process when tf_load_arrdet_g3vox == false or file not found
  //

  LOG_INFO("tf_load_arrdet_g3vox = false or file not found. calculate arrdet_g3vox and vec_spmat_PL.");

  std::chrono::system_clock::time_point t_create_start, t_create_end;
  t_create_start = time_now;

  LOG_INFO("====================================================================================");
  LOG_INFO("Create arrdet_g3vox from arrdet_built.");
  DetectorPanelArray arrdet_g3vox(prm_arrdet_vox.arrdet_built);

  // When tf_add_shell==true, add path length and density length of g2pil_shell (upper and lower)
  if(  prm_arrdet_vox.prm_pathcalc.tf_add_shell == true ){
    LOG_INFO("----------------------------------------");
    if( prm_arrdet_vox.has_shell_upper ){
      LOG_INFO("Add path length and density length of g2pil_shell_upper");
      pathcalc::g2pil::mp_add_PLDL(
        arrdet_g3vox,
        prm_arrdet_vox.g2pil_shell_upper
      );
    }
    if( prm_arrdet_vox.has_shell_lower ){
      LOG_INFO("Add path length and density length of g2pil_shell_lower");
      pathcalc::g2pil::mp_add_PLDL(
        arrdet_g3vox,
        prm_arrdet_vox.g2pil_shell_lower
      );
    }
    if( prm_arrdet_vox.has_shell_lateral ){
      LOG_INFO("Add path length and density length of g2pil_shell_lateral");
      pathcalc::g2pil::mp_add_PLDL(
        arrdet_g3vox,
        prm_arrdet_vox.g2pil_shell_lateral
      );
    }
  }

  LOG_INFO("----------------------------------------");
  LOG_INFO("Calculate vec_mat");

  // Display prm.reference_matPL_sparse, prm.epsilon_matPL_sparse
  LOG_DEBUG("prm.reference_matPL_sparse={:.1E}, prm.epsilon_matPL_sparse={:.1E}"
    , prm_arrdet_vox.prm_pathcalc.reference_matPL_sparse
    , prm_arrdet_vox.prm_pathcalc.epsilon_matPL_sparse);

  LOG_DEBUG("Set PL and DL to arrdet_g3vox and create observation matrix for each detector");
  std::vector<SpMatf> vec_spmat_PL_allvox =
    pathcalc::g3vox::mp_make_vec_spmat_PL(
      arrdet_g3vox,
      prm_arrdet_vox.g3vox,
      prm_arrdet_vox.prm_pathcalc
    );

  const int n_vec = vec_spmat_PL_allvox.size();
  LOG_DEBUG(" matrix size check: vec_spmat_PL_allvox.size()={}", vec_spmat_PL_allvox.size());
  for (int i = 0; i < n_vec; i++) {
    const int n_row = vec_spmat_PL_allvox.at(i).rows();
    const int n_col = vec_spmat_PL_allvox.at(i).cols();
    if (n_row == 0) THROW_ERROR("n_row == 0");
    if (n_col == 0) THROW_ERROR("n_col == 0");
    LOG_DEBUG("vec_spmat_PL_allvox.at({}).rows()={}, cols={}", i, n_row, n_col);
  }

  LOG_INFO("Assign new unique_index to g3vox with n_hit_det filter.");
  const bool tf_exist = true;
  const std::map<Grid3d::Uqiv,Grid3d::Uqiv> map_uqiv_old_new_avail
    = prm_arrdet_vox.g3vox.re_assign_uqiv_by_nhit_det(
        tf_exist,
        prm_arrdet_vox.prm_g3vox,
        prm_arrdet_vox.uqiv_start
      );
  LOG_INFO(" Remove columns for voxels that did not meet filter conditions for n_hit_det.");
  std::vector<SpMatf> vec_spmat_PL
    = detector_indexing::get_unavaible_cols_removed_matrix(vec_spmat_PL_allvox, map_uqiv_old_new_avail);

  LOG_DEBUG("size of vec_spmat_PL_allvox BEFORE:");
  eigen_blas::disp_vec_spmatf(vec_spmat_PL_allvox);
  LOG_DEBUG("size of vec_spmat_PL AFTER:");
  eigen_blas::disp_vec_spmatf(vec_spmat_PL);
  SLEEP_MSEC(500);

  // non_rec_vox = voxels filtered out by n_hit_det (not in rec target).
  // capture their per-element PL contribution before allvox is released.
  const int n_det_non_rec_vox = static_cast<int>(vec_spmat_PL_allvox.size());
  // build per-detector offsets into the flat output vector
  int total_n_ele = 0;
  std::vector<int> offsets_non_rec_vox(n_det_non_rec_vox);
  for (int detid = 0; detid < n_det_non_rec_vox; ++detid) {
    offsets_non_rec_vox[detid] = total_n_ele;
    total_n_ele += static_cast<int>(vec_spmat_PL_allvox[detid].rows());
  }
  // flat output: one slot per detector element across all detectors
  Eigen::VectorXf vecxf_non_rec_vox_PL = Eigen::VectorXf::Zero(total_n_ele);
  for (int detid = 0; detid < n_det_non_rec_vox; ++detid) {
    const SpMatf& spmatf_all  = vec_spmat_PL_allvox[detid]; // before n_hit_det filter
    const SpMatf& spmatf_filt = vec_spmat_PL[detid];        // after  n_hit_det filter
    const int n_ele = static_cast<int>(spmatf_all.rows());

    // Step 1: build [1, 1, ..., 1] of length = voxel count for each spmatf.
    Eigen::VectorXf vecxf_ones_all  = Eigen::VectorXf::Ones(spmatf_all.cols());
    Eigen::VectorXf vecxf_ones_filt = Eigen::VectorXf::Ones(spmatf_filt.cols());

    // Step 2: spmatf * vecxf_ones sums each row -> per-element total PL.
    Eigen::VectorXf vecxf_pl_all  = spmatf_all  * vecxf_ones_all;   // PL over all voxels
    Eigen::VectorXf vecxf_pl_filt = spmatf_filt * vecxf_ones_filt;  // PL over rec-target voxels only

    // (all-voxel PL) - (rec-target PL) = PL through the n_hit_det-dropped voxels.
    // .segment(start, len): writable view of [start, start+len); LHS assigns in-place.
    vecxf_non_rec_vox_PL.segment(offsets_non_rec_vox[detid], n_ele)
      = vecxf_pl_all - vecxf_pl_filt;
  }
  LOG_INFO("vecxf_non_rec_vox_PL computed (size={}, sum={:.3E})",
           total_n_ele, vecxf_non_rec_vox_PL.sum());

  LOG_INFO(" Release memory of vec_spmat_PL_allvox.");
  myapp::clearVectorMemory(vec_spmat_PL_allvox);

  LOG_INFO(" Calculate peneflux and signal of arrdet_g3vox.");
  arrdet_g3vox.mp_calc_set_peneflux_signal_from_DL(
    prm_arrdet_vox.ft, prm_arrdet_vox.tf_apply_eff, prm_arrdet_vox.tf_apply_eff_cnt
  );

  // give the noise distribution from prm_noise to arrdet_g2pil
  LOG_INFO("arrdet_g3vox.mp_set_noise_all");
  const double DL_thres = prm_arrdet_vox.prm_bingroup.DL_thres;
  arrdet_g3vox.mp_set_noise_all(prm_arrdet_vox.prm_noise, DL_thres);

  LOG_INFO("Set Grid2dBinGroup::Parameters to arrdet_g3vox_input.");
  arrdet_g3vox.set_parameters_all(prm_arrdet_vox.prm_bingroup);
  LOG_DEBUG("(prm_arrdet_vox.prm_bingroup)");

  LOG_INFO("set all detector element signal to Grid2dBinGroup vec_vec_signal");
  arrdet_g3vox.mp_copy_signal_noise_to_g2bg_all(prm_arrdet_vox.prm_bingroup);
  LOG_DEBUG("(prm_arrdet_vox.prm_bingroup)");

  if( opt_arrdet_built.has_value() ){
    // When opt_arrdet_built is valid, set dic to arrdet_g3vox.
    LOG_INFO("set dic to arrdet_g3vox");
    arrdet_g3vox.set_dic_bimap_all_from(opt_arrdet_built.value());
    LOG_DEBUG("(opt_arrdet_built.value())");
  } else {
    // From vec_vec_signal, automatically set appropriate bin divisions and group minimum angle elements.
    if( prm_arrdet_vox.tf_run_1st_grouping == false ){
      LOG_INFO("(g3vox) tf_run_1st_grouping == false. So auto grouping is not executed.");
      SLEEP_MSEC(500);
    } // tf_1st_grouping==true below
    else{
      // assign 1st igroup for all detectors
      LOG_INFO("assign 1st igroup for all detectors");
      arrdet_g3vox.mp_assign_1st_igroup_all( 
        prm_arrdet_vox.prm_bingroup.tf_run_1st_grouping, prm_arrdet_vox.prm_bingroup.igroup_start
      , prm_arrdet_vox.prm_bingroup.nx_div_init, prm_arrdet_vox.prm_bingroup.ny_div_init
      , prm_arrdet_vox.prm_bingroup.signal_noise_group_trig );
    }
    if( prm_arrdet_vox.tf_run_auto_grouping == false ){
      LOG_INFO("(g3vox) tf_run_auto_grouping == false. So auto grouping is not executed.");
      if( prm_arrdet_vox.prm_bingroup.n_detector_grouping_manual > 0 ){
        // Manual grouping from user-supplied bin-list files (tx/ty rectangles).
        arrdet_g3vox.mp_grouping_by_bin_list_all(prm_arrdet_vox.prm_bingroup);
      }
      SLEEP_MSEC(500);
    }else{
      // tf_run_1st_grouping && tf_run_auto_grouping == true
      LOG_INFO("(g3vox) tf_run_1st_grouping && tf_run_auto_grouping == true.");
      LOG_INFO("(g3vox) From vec_vec_signal, automatically set appropriate bin divisions and group minimum angle elements.");
      arrdet_g3vox.mp_auto_grouping_by_signal_noise_group_alldet(prm_arrdet_vox.prm_bingroup);

    }
    // Check that no group is smaller than ixlen_min x iylen_min (covers all grouping routes).
    arrdet_g3vox.mp_check_group_ixiylen_min_all(
      prm_arrdet_vox.prm_bingroup.ixlen_min, prm_arrdet_vox.prm_bingroup.iylen_min );
    LOG_INFO(" Build dic_ of arrdet_g3vox.");
    arrdet_g3vox.build_index_container(); // Build dic_

  }
  LOG_INFO("Allocate memory for vec_signal/noise_poi_group in each Group of arrdet_g3vox.");
  arrdet_g3vox.mp_allocate_vec_value_group_all();

  // group efficiency needs completed auto grouping; skip when grouping was not run
  if( arrdet_g3vox.get(FlgProg::tf_auto_grouping) ){
    LOG_INFO("Set group efficiency (simple mean of element efficiencies) of arrdet_g3vox.");
    arrdet_g3vox.calc_set_eff_group();
  }else{
    LOG_WARN("calc_set_eff_group is skipped because auto grouping is not done.");
  }

  LOG_INFO("Calculate vec_tf_in_PL");
  if( opt_arrdet_built.has_value() ){
    // When opt_arrdet_built is valid, set vec_tf_in_PL to arrdet_g3vox.
    LOG_INFO("set vec_tf_in_PL to arrdet_g3vox");
    arrdet_g3vox.mp_copy_vec_tf_in_PL_all(opt_arrdet_built.value());
    LOG_DEBUG("(opt_arrdet_built.value())");
  }else{
    LOG_WARN("vec_tf_in_PL is not set because opt_arrdet_built is not set.");
  }

  if( prm_arrdet_vox.tf_run_1st_grouping == true ){
    LOG_INFO("Execute get_g2bg().calc_vec_signal_group() for all detectors in arrdet_g3vox,");
    LOG_INFO("Calculate signal_group for each Group and set to vec_signal_group.");
    arrdet_g3vox.mp_calc_vec_signal_noise_group_all();
  }
  arrdet_g3vox.display_status();

  if (prm_arrdet_vox.tf_save_arrdet_g3vox == true) {
    LOG_INFO("save arrdet_g3vox to binary file.");
    arrdet_g3vox.save(prm_arrdet_vox.path_arrdet_bin);
    LOG_DEBUG("(prm_arrdet_vox.path_arrdet_bin)");

    LOG_INFO("save vec_spmat_PL to binary file.");
    io_binary::write_vec_spmatf(prm_arrdet_vox.path_vec_spmat_PL_bin, vec_spmat_PL);
    LOG_DEBUG("(prm_arrdet_vox.path_vec_spmat_PL_bin, vec_spmat_PL)");
  }

  LOG_INFO("get the set of new uqiv from map_uqiv_old_new_avail");
  arrdet_g3vox.set_set_uqiv( tuple_int::get_set_value(map_uqiv_old_new_avail) );

  arrdet_g3vox.set_name("arrdet_g3vox");
  arrdet_g3vox.display_status();

  t_create_end = time_now;
  myapp::cast_time_msec(spdlog::level::info,"DetectorPanelArray::create_sprs(VoxelBuildParams) finished",t_create_start,t_create_end);

  LOG_INFO("return arrdet_g3vox, vec_spmat_PL, vecxf_non_rec_vox_PL");
  return std::make_tuple(arrdet_g3vox, vec_spmat_PL, vecxf_non_rec_vox_PL);
}

void DetectorPanelArray::out_vec_panel( ) const
{
  fs::path fullpath;
  for(const auto& panel : vec_panel){
    fs::path pathout = iodir::make_pathout( name + "_" + panel.get_name() + ".tmp" );
    LOG_DEBUG("pathout={}",pathout.string());
    panel.out_all(pathout);
  }
}

int DetectorPanelArray::build_vec_panel(
  const Uqid uqid_start, const Detid detid_start )
{
  if(get(FlgProg::tf_built_panels)){
    LOG_WARN("already built vec_panel.");
    SLEEP_MSEC(500);
  }

  // call mutable ref of UqidManager
  UqidManager &uqid_mgr = dic_.callUqidMgr();

  LOG_INFO("uqid_start={}, detid_start={}",uqid_start,detid_start);

  // memory allocation
  const int n_det = vec_parameter_file_path.size(); // number of detector panels
  LOG_INFO("n_det={}",n_det);

  // reserve memory for vec_panel
  vec_panel.reserve(n_det);

  // set the detid_max_ of dic_
  dic_.set_detid_max(detid_start + n_det - 1);

  // count number of all elements
  size_t n_all_element_sum = 0;
  // for(const auto& rc_fname : vec_parameter_file_path){

  // Load all json files first
  std::vector<DetectorPanel::Parameters> vec_prm;
  vec_prm.reserve(n_det);
  for (auto& path : vec_parameter_file_path) {
    vec_prm.push_back( DetectorPanel::Parameters{ myapp::load_json(path) } );
  }

  // Apply angle_bin_override if enabled
  if (tf_override_angle_bin) {
    LOG_INFO("Applying angle_bin_override: nbinx={}, nbiny={}", nbinx_override, nbiny_override);
    for (auto& prm : vec_prm) {
      prm.nbinx = nbinx_override;
      prm.txmin = txmin_override;
      prm.txmax = txmax_override;
      prm.nbiny = nbiny_override;
      prm.tymin = tymin_override;
      prm.tymax = tymax_override;
    }
  }

  // Create storage vector for the number of DetectorElements in each detector
  std::vector<int> vec_n_element(n_det, 0);

  for(int i=0;i<n_det;i++){
    const Detid detid = i + detid_start;
    const DetectorPanel::Parameters &prm = vec_prm.at(i);
    std::string det_name = prm.name;
    LOG_DEBUG("det_name={}",det_name);  

    const int n_element = prm.nbinx * prm.nbiny;
    vec_n_element.at(i) = n_element; // store number of elements for each detector
    n_all_element_sum += static_cast<size_t>(n_element);
  }
  this->set_n_all_element(n_all_element_sum);

  //
  // make accumlative sum of vector of vec_n_element
  // 
  std::vector<int> vec_n_element_accum(n_det, 0);
  for(int i=1;i<n_det;i++){
    vec_n_element_accum.at(i) = vec_n_element_accum.at(i-1) + vec_n_element.at(i-1);
  }

  //
  // allocate memory for unique index map container
  //
  LOG_DEBUG("n_det={}, n_all_element_sum={}",n_det,n_all_element_sum);
  // Set Uqid min/max once (via DetectorIndexContainer)
  const Uqid uqid_min = uqid_start;
  const Uqid uqid_max = uqid_min + static_cast<Uqid>(n_all_element_sum) - 1;
  uqid_mgr.set_uqid_min(uqid_min);
  uqid_mgr.set_uqid_max(uqid_max);
  
  // This std::vector<std::vector<UqidInfo>> vec_vec_UqidInfo must be thread-safe.
  std::vector<std::vector<UqidInfo>> vec_vec_UqidInfo(n_det);

  // resize vec_panel
  vec_panel.resize(n_det);

  // loop of detid
  #pragma omp parallel for schedule(static)
  for(int i=0;i<n_det;i++){
    const Detid detid = i + detid_start;
    const fs::path json_path = vec_parameter_file_path.at(i);
    const DetectorPanel::Parameters &prm = vec_prm.at(i);
    std::string det_name = prm.name;
    LOG_DEBUG("det_name={}",det_name);  

    // Set uqid_min, uqid_max for each detector here.
    const Uqid uqid_min_tmp = uqid_start + vec_n_element_accum.at(i);
    const Uqid uqid_max_tmp = uqid_min_tmp + static_cast<Uqid>(vec_n_element.at(i)) - 1;

    // build DetectorPanel with assigning unique index and making the uomap
    // after this function, detid become +1,
    // and unique_index +the number of detector element
    LOG_DEBUG("making panel with uqid_max_tmp={}, detid={}",uqid_max_tmp,detid);

    // Declare std::vector<UqidInfo>& vec_UqidInfo
    std::vector<UqidInfo> vec_UqidInfo_tmp;
    vec_UqidInfo_tmp.reserve(vec_n_element.at(i)); // reserve memory for UqidInfo

    // using OpenMP to parallelize the loop
    // dic_ is simply passed as a pointer.
    DetectorPanel panel(prm, detid, uqid_min_tmp, uqid_max_tmp, vec_UqidInfo_tmp, dic_);

    // Assign to thread-safe vec_vec_UqidInfo
    vec_vec_UqidInfo.at(i) = std::move(vec_UqidInfo_tmp);
    
    LOG_DEBUG("panel.get_n_element()={}",panel.get_n_element()); 
    LOG_DEBUG("panel.get_name()={}",panel.get_name());
    LOG_DEBUG("panel.get_detid()={}",panel.get_detid());

    // path_eff_model (JSON5 efficiency model) takes precedence over the
    // legacy path_eff_table (7-column text table) during migration.
    if (prm.eff_model.get_tf_loaded()) {
      panel.mp_assign_efficiency_model(prm.eff_model);
      LOG_DEBUG("panel.mp_assign_efficiency_model({}) done", prm.path_eff_model.string());
    } else {
      const std::string eff_table_path = prm.path_eff_table;
      panel.read_efficiency_table(eff_table_path);
      LOG_DEBUG("panel.read_efficiency_table({}) done",eff_table_path);
    }
    
    vec_panel.at(i) = std::move(panel); // move panel to vec_panel

  } // end of detid loop
  LOG_INFO("n_all_element_sum={}",n_all_element_sum);

  //
  // Set vec_vec_UqidInfo to UqidManager here.
  // ** This cannot be parallelized. **
  //
  LOG_DEBUG("Initialize bimap of uqid_mgr");
  uqid_mgr.initialize(uqid_min);

  LOG_DEBUG("Allocate memory for uqid_mgr");
  uqid_mgr.reserveAdditional(n_all_element_sum);

  LOG_DEBUG("Finally, set vec_vec_UqidInfo to UqidManager.");
  for(size_t i=0;i<n_det;i++){
    const Detid detid = i + detid_start;
    std::vector<UqidInfo> &vec_UqidInfo = vec_vec_UqidInfo.at(i);
    for(size_t j=0;j<vec_UqidInfo.size();j++){
      UqidInfo &uqid_info = vec_UqidInfo.at(j);
      // insert to UqidManager
      uqid_mgr.insertInit(uqid_info);
    }
  }
  // set uqid_min, uqid_max to UqidManager
  uqid_mgr.set_uqid_min(uqid_min);
  uqid_mgr.set_uqid_max(uqid_max);

  // set flags
  set(FlgProg::tf_built_panels,true);
  dic_.set_built_uqid(true);

  return n_all_element_sum;
}

// build Detector Index Container
// Return value: maximum values of uqig, uqig_avail
std::tuple<Uqig,UqigAvail> DetectorPanelArray::build_index_container(
  const Uqig uqig_min, const UqigAvail uqig_avail_min )
{
  LOG_INFO(
    "uqig_min={}, uqig_avail_min={}"
    ,uqig_min, uqig_avail_min);

  // check if uqid_maps is built, and auto grouping all, if not throw error
  // uqid_maps.check_build_uqid();
  dic_.check_built_uqid();
  check_auto_grouping();

  // call instance of mutable GroupManager and UqidManager
  GroupManager &grp_mgr = dic_.callGrpMgr();
  UqidManager &uqid_mgr = dic_.callUqidMgr();

  // if tf_built_uqig is true, LOG_WARN
  if(dic_.is_built_uqig()) LOG_WARN("==true");
  
  // if tf_build_uqigAvail is true, LOG_WARN
  if(dic_.is_built_uqigAvail()) LOG_WARN("==true");

  // initialize GroupManager
  grp_mgr.initialize();

  // unique index of igroup for all detectors
  Uqig uqig = uqig_min;
  UqigAvail uqig_avail = uqig_min;
  grp_mgr.set_uqig_min(uqig);
  grp_mgr.set_uqig_avail_min(uqig_avail);

  // get number of detectors
  const int n_det = get_n_det();

  // for each detector
  for(Detid detid=0;detid<n_det;detid++){
    DetectorPanel& panel = callDetectorPanel(detid);
    Grid2dBinGroup &g2bg = panel.call_g2bg();
    // const std::vector<Igroup> vec_igroup = g2bg.get_vec_igroup();
    const std::set<Igroup> set_igroup = g2bg.get_set_igroup();
    const int n_group = set_igroup.size();
    g2bg.allocate_vec_value_group();
    // for each igroup, get set of uqid, assign to uqig
    // and insert to uqid_maps.map_uqig_uqids
    LOG_DEBUG("build_index_container detid={}, n_group={}",detid,n_group);


    // loop for each igroup
    for(const int &igroup : set_igroup){
      // check if igroup is available
      const bool is_avail = g2bg.is_avail_group(igroup);

      const DetIgroup det_igroup = {detid,igroup};
      const std::vector<Uqid> vec_uqid = panel.get_vec_uqid_by_igroup(igroup);
      // const int n_uqid = set_uqid.size();
      // mylogger::g_logger->trace("detid={}, igroup={}, n_uqid={}",detid,igroup,n_uqid);

      // make temporary instance of GroupInfo and set members
      GroupInfo ginfo;
      ginfo.uqig = uqig;
      ginfo.detid = detid;
      ginfo.igroup = igroup;
      ginfo.is_avail = is_avail;

      // If avail is true, assign uqig_avail, otherwise set NotAssigned
      UqigAvail current_uqigAvail = is_avail
        ? uqig_avail
        : UqigAvailNotAssigned;
      ginfo.uqig_avail = current_uqigAvail;

      // #pragma omp parallel for
      for(int idx = 0; idx < vec_uqid.size(); ++idx){
        const Uqid uqid = vec_uqid.at(idx);
        UqidInfo &uinfo = uqid_mgr.callInfo(uqid);
        if(uinfo.detid != detid){
          LOG_ERROR(
            "build_index_container : uinfo.detid={} != detid={} for uqid={}"
            ,uinfo.detid, detid, uqid
          );
          THROW_ERROR_NAME("uqid_info.detid != detid");
        }
        uinfo.igroup       = igroup;
        uinfo.uqig         = uqig;
        uinfo.uqig_avail   = current_uqigAvail;
        uinfo.is_avail     = is_avail;
      }
      // insert to map_uqig_uqids, map_uqig_detigroup
      // uqid_maps.insert_to_map_uqig_detigroup(uqig,det_igroup);
      grp_mgr.insert_to_Bimaps(ginfo);
      
      // If is_avail==true, increment uqig_avail by one group
      if( is_avail ) uqig_avail++;
      
      // increment uqig
      uqig++;
    }
  }
  LOG_INFO("n_uqig={}",uqig-uqig_min);
  // return ther number of unique igroup

  // uqid_maps.set_uqig_max(uqig-1);
  // uqid_maps.set_uqig_avail_max(uqig_avail-1);
  grp_mgr.set_uqig_max(uqig-1);
  grp_mgr.set_uqig_avail_max(uqig_avail-1);

  // set tf_build_uqig
  dic_.set_built_uqig(true);

  // set tf_build_uqigAvail
  dic_.set_built_uqigAvail(true);

  // build uqigAvailIndexer
  LOG_INFO("build uqigAvailIndexer");
  dic_.build_uqigAvailIndexer();

  // return maximum uqig and uqig_avail
  LOG_INFO(
    "uqig_max={}, uqig_avail_max={}"
    ,uqig-1, uqig_avail-1
  );
  return std::make_tuple(uqig-1,uqig_avail-1);
}

// Allocate vec_value_group for all DetectorPanels
void DetectorPanelArray::mp_allocate_vec_value_group_all( )
{
  // check if uqid maps is built, if not throw error
  check_auto_grouping();

  // allocate vec_value_group
  const int n_det = get_n_det();
  #pragma omp parallel for
  for(Detid detid=0;detid<n_det;detid++){
    DetectorPanel& panel = callDetectorPanel(detid);
    panel.call_g2bg().allocate_vec_value_group();
  }
  set(FlgProg::tf_calc_signal_group,false);
  set(FlgProg::tf_calc_noise_group,false);
  set(FlgProg::tf_calc_proj_dens_group,false);
}

// 2022-12-05 16:26:29
// return arrdetay = this - arrdet_sub
void DetectorPanelArray::add_sub_path_and_DL(
  const DetectorPanelArray &arrdet_in, const double sign )
{ 
  // check if uqid maps is built, if not throw error
  dic_.check_built_uqid();

  LOG_INFO("return arrdetay = this - arrdet_sub");
  if( this->get_n_all_element() != arrdet_in.get_n_all_element() )
    THROW_ERROR_NAME("this->get_n_all_element() != arrdet_in.get_n_all_element()");

  const Uqid uqid_min = get_uqid_min();
  const Uqid uqid_max = get_uqid_max();
  for(Uqid uqid=uqid_min; uqid<=uqid_max; uqid++ ){
    DetectorElement& ele_this = this->callDetectorElement(uqid);
    const DetectorElement& ele_in = arrdet_in.getDetectorElement(uqid);

    // calc rocklength difference and assign the value
    ele_this.set_PL(
      ele_this.get_PL() + sign * ele_in.get_PL()
    );

    // calc density length difference and assign the value
    ele_this.set_DL(
      ele_this.get_DL() + sign * ele_in.get_DL()
    );
  }
}

/// @brief get the exposure time (sec) of Detid=detid_in
double DetectorPanelArray::get_exposure_time_sec_panel( const Detid detid_in ) const
{
  const DetectorPanel& panel = getDetectorPanel(detid_in);
  const double time_sec = panel.get_exposure_time_sec();
  return time_sec;
}


Eigen::VectorXf DetectorPanelArray::get_vecxf_nmuon_all() const
{
  // ——————————————————————————————————
  // 1) Prerequisite check
  //——————————————————————————————————
  dic_.check_built_uqid();  // Check if uqid range is properly built

  //——————————————————————————————————
  // 2) Get group information and row indexer
  //——————————————————————————————————
  const auto& gm = dic_.getGrpMgr();
  const auto& ai = dic_.getAvailIndexer();

  // Get all GroupInfo at once (each element contains uqig_avail, detid, igroup, members)
  std::vector<GroupInfo> vec_group_info = gm.get_vecGroupInfo_all();

  // Output vector length = number of indexer rows
  const int n_rows = ai.size();
  if (n_rows <= 0) {
    THROW_ERROR_NAME("n_rows <= 0, n_rows={}", n_rows);
  }

  Eigen::VectorXf vec(n_rows);
  vec.setZero(); // Initialize

  //——————————————————————————————————
  // 3) Scan each group information and fill signal+noise
  //——————————————————————————————————
  // Scan all group information and copy corresponding rows
  const std::vector<Index> vec_index = ai.get_vecIndex();
  const size_t n_index = vec_index.size();

  if (n_index != n_rows) {
    LOG_ERROR("n_index != n_rows, n_index={}, n_rows={}", n_index, n_rows);
    THROW_ERROR_NAME("n_index != n_rows, n_index=" + std::to_string(n_index)
      + ", n_rows=" + std::to_string(n_rows));
  }

  #pragma omp parallel for
  for (size_t i = 0; i < n_index; ++i) {
    const Index row = vec_index.at(i);
    const GroupInfo& ginfo = dic_.getGroupInfo_by_Index( row );

    // Extract the corresponding detector panel
    const DetectorPanel& panel = getDetectorPanel(ginfo.detid);

    // Get signal_group, noise_poi_group
    double sig = panel.get_signal_group(ginfo.igroup);
    double noi = panel.get_noise_poi_group( ginfo.igroup );

    vec(row) = sig + noi;
  }

  return vec;
}

/// @brief Returns vector with log10(signal_group + noise_poi_group) arranged in uqig_avail order
/// @param nmuon_group_thres  Threshold below which log_value_under_thres is used
/// @param log_value_under_thres  Log value to fill when below threshold
/// @throws std::out_of_range When uqid or group is not built
Eigen::VectorXf DetectorPanelArray::get_vecxf_log_nmuon_all(
  const double nmuon_group_thres
, const double log_value_under_thres) const
{
  // Check uqid→group construction
  dic_.check_built_uqid();
  check_auto_grouping();

  // Get each manager/indexer reference
  const auto& gm = dic_.getGrpMgr();
  const auto& ai = dic_.getAvailIndexer();

  // Output vector length = number of indexer rows
  const int n_rows = ai.size();
  Eigen::VectorXf vec(n_rows);
  vec.setZero();

  // Process all groups in order
  for(const auto& g : gm.get_vecGroupInfo_all()){
    // Get row index
    const int row = ai.getIndex(g.uqig_avail);

    // Panel & group number
    const DetectorPanel& panel = getDetectorPanel(g.detid);
    const double sig = panel.get_signal_group(g.igroup);
    const double noi = panel.get_noise_poi_group( g.igroup );
    const double total = sig + noi;

    // Check threshold and set log value
    vec(row) = (total < nmuon_group_thres)
             ? log_value_under_thres
             : std::log10(total);
  }

  return vec;
}

/// @brief Returns per-row (signal,noise) counts.
/// @param tf_signal_poisson whether to apply Poisson error to signal counts
/// @note noise term is always det floor + Poisson(poi bucket); no noise switch.
/// @note Mirrors get_vecxf_nmuon_all()'s row-indexer + OpenMP structure.
/// @throws std::out_of_range When uqid or group is not built
Eigen::VectorXf DetectorPanelArray::get_vecxf_nmuon_poisson_all(
  const bool tf_signal_poisson) const
{
  // Prerequisite check (same as get_vecxf_nmuon_all)
  dic_.check_built_uqid();

  const auto& ai = dic_.getAvailIndexer();

  // Output vector length = number of indexer rows
  const int n_rows = ai.size();
  if (n_rows <= 0) {
    THROW_ERROR_NAME("n_rows <= 0, n_rows={}", n_rows);
  }

  Eigen::VectorXf vec(n_rows);
  vec.setZero();

  const std::vector<Index> vec_index = ai.get_vecIndex();
  const size_t n_index = vec_index.size();
  if (n_index != n_rows) {
    LOG_ERROR("n_index != n_rows, n_index={}, n_rows={}", n_index, n_rows);
    THROW_ERROR_NAME("n_index != n_rows, n_index=" + std::to_string(n_index)
      + ", n_rows=" + std::to_string(n_rows));
  }

  // Per-component Poisson selection. The Poisson values are precomputed reads,
  // so calling them inside the parallel loop is data-race free.
  #pragma omp parallel for
  for (size_t i = 0; i < n_index; ++i) {
    const Index row = vec_index.at(i);
    const GroupInfo& ginfo = dic_.getGroupInfo_by_Index( row );
    const DetectorPanel& panel = getDetectorPanel(ginfo.detid);
    const double sig = tf_signal_poisson
      ? panel.get_signal_group_poisson(ginfo.igroup)
      : panel.get_signal_group(ginfo.igroup);
    // noise: always det floor + Poisson(poi bucket)
    const double noi = panel.get_noise_det_group(ginfo.igroup)
      + panel.get_noise_poi_group_poisson(ginfo.igroup);
    vec(row) = sig + noi;
  }

  return vec;
}

// Build per-row efficiency variance for the C_N diagonal (efficiency uncertainty).
// Uses the same AvailIndexer row order as get_vecxf_nmuon_poisson_all, so the
// resulting vector is aligned row-for-row with the C_N diagonal.
Eigen::VectorXf DetectorPanelArray::get_vecxf_var_eff_all(
  const bool tf_independent) const
{
  // Prerequisite check (same as get_vecxf_nmuon_poisson_all)
  dic_.check_built_uqid();

  const auto& ai = dic_.getAvailIndexer();

  const int n_rows = ai.size();
  if (n_rows <= 0) {
    THROW_ERROR_NAME("n_rows <= 0, n_rows={}", n_rows);
  }

  Eigen::VectorXf vec(n_rows);
  vec.setZero();

  const std::vector<Index> vec_index = ai.get_vecIndex();
  const size_t n_index = vec_index.size();
  if (n_index != n_rows) {
    LOG_ERROR("n_index != n_rows, n_index={}, n_rows={}", n_index, n_rows);
    THROW_ERROR_NAME("n_index != n_rows, n_index=" + std::to_string(n_index)
      + ", n_rows=" + std::to_string(n_rows));
  }

  // Per-row efficiency variance. Same indexer as get_vecxf_nmuon_poisson_all.
  #pragma omp parallel for
  for (size_t i = 0; i < n_index; ++i) {
    const Index row = vec_index.at(i);
    const GroupInfo& ginfo = dic_.getGroupInfo_by_Index( row );
    const DetectorPanel& panel = getDetectorPanel(ginfo.detid);

    // Accumulate sigma_eff_i * b_i over the merge group elements.
    double sum_contrib = 0.0; // for common mode (square of sum)
    double sum_sq      = 0.0; // for independent mode (sum of squares)
    const auto vec_ixiy = panel.get_g2bg().get_vec_ixiy(ginfo.igroup);
    for (const auto& [ix, iy] : vec_ixiy) {
      const DetectorElement& ele = panel.getDetectorElement(ix, iy);
      // efficiency sigma from the (lower, upper) band half-width
      const double sigma_eff = 0.5 * (ele.get_eff_upp() - ele.get_eff_low());
      // base count b_i: efficiency-free signal (get_signal() may be eff-applied)
      const double base_count = ele.calc_signal();
      const double contrib = sigma_eff * base_count;
      sum_contrib += contrib;
      sum_sq      += contrib * contrib;
    }

    vec(row) = static_cast<float>(
      tf_independent ? sum_sq : sum_contrib * sum_contrib);
  }

  return vec;
}

/// @brief Returns vector with signal arranged in uqid order for all detectors
Eigen::VectorXf
DetectorPanelArray::get_vecxf_signal_all() const
{
  if (!dic_.is_uqid_range_valid())
    THROW_ERROR_NAME("dic_.is_uqid_range_valid() == false, uqid range is not set");

  const Uqid uqid_min = dic_.get_uqid_min();
  const Uqid uqid_max = dic_.get_uqid_max();
  const int  n       = static_cast<int>(uqid_max - uqid_min + 1);

  Eigen::VectorXf vec(n);
  vec.setZero();

  for (Uqid u = uqid_min; u <= uqid_max; ++u) {
    // Get uqid → detid, ix, iy
    const DetIxiy detixiy = dic_.get_detidixiy(u);

    // Get pointer to the element and read signal
    const DetectorElement& ele = getDetectorElement(detixiy);

    // Assign considering offset
    vec[static_cast<int>(u - uqid_min)] = ele.get_signal();
  }

  return vec;
}

/// @brief Returns vector with signal for each bin group arranged in uqig_avail order for all detectors
/// @throws std::out_of_range When uqid or group is not built
Eigen::VectorXf DetectorPanelArray::get_vecxf_nmuoned_all() const
{
  // Check uqid range and group construction
  dic_.check_built_uqid();
  check_auto_grouping();

  // Manager/indexer reference
  const auto& gm = dic_.getGrpMgr();
  const auto& ai = dic_.getAvailIndexer();

  // Output length = number of unique uqig_avail
  const int n_rows = ai.size();
  Eigen::VectorXf vec(n_rows);
  vec.setZero();

  // Scan all groups and fill signal
  for (const auto& g : gm.get_vecGroupInfo_all()) {
    int row = ai.getIndex(g.uqig_avail);
    const DetectorPanel& panel = getDetectorPanel(g.detid);
    double sig = panel.get_g2bg().calc_signal_group(g.igroup);
    vec(row) = sig;
  }

  return vec;
}

/// @brief Returns 3-element array with eff_low/cnt/upp_group collected in uqig_avail order for all detectors
/// @throws std::out_of_range When dic_ is not built or grouping is incomplete
std::array<Eigen::VectorXf,3>
DetectorPanelArray::get_vecxf_eff_low_cnt_upp_group_all() const
{
  // Initialization check
  dic_.check_built_uqid();
  check_auto_grouping();

  // Get manager/indexer
  const auto& gm = dic_.getGrpMgr();
  const auto& ai = dic_.getAvailIndexer();

  // Output vector length = number of unique uqig_avail
  const int n = ai.size();
  Eigen::VectorXf vec_low(n), vec_cnt(n), vec_upp(n);
  vec_low.setZero();
  vec_cnt.setZero();
  vec_upp.setZero();

  // Get and set signal for each group
  for (const auto& g : gm.get_vecGroupInfo_all()) {
    int row = ai.getIndex(g.uqig_avail);
    const DetectorPanel& panel = getDetectorPanel(g.detid);
    vec_low(row) = panel.get_g2bg().get_eff_low_group(g.igroup);
    vec_cnt(row) = panel.get_g2bg().get_eff_cnt_group(g.igroup);
    vec_upp(row) = panel.get_g2bg().get_eff_upp_group(g.igroup);
  }

  return { vec_low, vec_cnt, vec_upp };
}

std::vector<Eigen::VectorXf> DetectorPanelArray::get_vec_vecxf_DL() const
{
  LOG_INFO("Return std::vector of vecxf_DL for all detectors.");
  const int n_det = get_n_det();
  std::vector<Eigen::VectorXf> vec_vecxf_DL;
  for(Detid detid=0;detid<n_det;detid++){
    const DetectorPanel& panel = getDetectorPanel(detid);
    vec_vecxf_DL.push_back(panel.mp_get_vecxf_DL());
  }
  return vec_vecxf_DL;
}

// @brief Returns std::vector of vecxf_DL_prior for all detectors assuming uniform prior density.
std::vector<Eigen::VectorXf>
  DetectorPanelArray::get_vec_vecxf_PL_times_dens( const double dens ) const
{
  LOG_INFO("Return std::vector of vecxf_PL_times_dens for all detectors.");
  const int n_det = get_n_det();
  std::vector<Eigen::VectorXf> vec_vecxf_PL_times_dens;
  for(Detid detid=0;detid<n_det;detid++){
    const DetectorPanel& panel = getDetectorPanel(detid);
    vec_vecxf_PL_times_dens.push_back(panel.mp_get_vecxf_PL()*dens);
  }
  return vec_vecxf_PL_times_dens;
}

/// @brief Return unique_index (uqid) from det_id, ix, iy
Uqid DetectorPanelArray::get_uqid(
  const Detid det_id, const int ix, const int iy ) const
{
  const DetIxiy detixiy = { det_id, ix, iy };
  return dic_.getUqidMgr().getUqidByDetIxiy(detixiy);
};

/// @brief Get det_id from uqid
Detid DetectorPanelArray::get_detid(const Uqid uqid_in) const
{
  return dic_.getUqidMgr().getInfo(uqid_in).detid;
};


/// @brief Convert mat(irow=uqig) to mat(irow=uqig_avail)
/// @param mat Input matrix (row = uqig)
/// @return Output matrix (row = uqig_avail)
/// @throws std::out_of_range When uqid/dic is not built or grouping is incomplete
Eigen::MatrixXf
DetectorPanelArray::get_matrix_uqig_avail_from_matrix_uqig(
  const Eigen::MatrixXf& mat) const
{
  // Check uqid range and group construction
  dic_.check_built_uqid();
  dic_.check_built_uqig();
  dic_.check_built_uqigAvail();
  check_auto_grouping();

  // Reference each manager/indexer
  const GroupManager& grp_mgr = dic_.getGrpMgr();
  const UqigAvailIndexer& indexer = dic_.getAvailIndexer();

  // Output row count = number of unique uqig_avail, column count same as input
  const int n_avail = indexer.size();
  const int n_cols  = mat.cols();
  // Eigen::MatrixXf mat_ret = Eigen::MatrixXf::Zero(n_avail, n_cols);

  // Scan all group information and copy corresponding rows
  const std::vector<Index> vec_index = indexer.get_vecIndex();
  const size_t n_index = vec_index.size();
  Eigen::MatrixXf mat_ret = Eigen::MatrixXf::Zero(n_index, n_cols);

  #pragma omp parallel for
  for (size_t i = 0; i < n_index; ++i) {
    const GroupInfo& ginfo = dic_.getGroupInfo_by_Index( vec_index.at(i) );
    const Index index_row = indexer.getIndex( ginfo.uqig_avail );
    mat_ret.row( index_row ) = mat.row( ginfo.uqig );
  }

  return mat_ret;
}

// @brief Calculate volume-weighted average density for all groups
// @details Get volume and density for each group and return volume-weighted average
// @return Volume-weighted average density
// @throws std::out_of_range When dic_ is not built or grouping is incomplete
double DetectorPanelArray::calc_volume_weighted_average_density(
  const fs::path& prefix_out )
{
  // Check uqid/group construction
  dic_.check_built_uqid();
  check_auto_grouping();

  // Get group list
  const auto& gm = dic_.getGrpMgr();

  double sum_vol  = 0.0;
  double sum_mass = 0.0;

  std::string fname = prefix_out.string() + "_cnt.tmp";
  fs::path pathout = iodir::make_pathout(fname);
  LOG_DEBUG("Output detailed average density and volume to {}", pathout.string());
  FILE *fout = myapp::get_fout(pathout);
  fprintf(fout,"# detid igroup vol dens mass\n");

  // Scan each group
  for (const auto& g : gm.get_vecGroupInfo_all()) {
    // Corresponding panel and group number
    DetectorPanel& panel = callDetectorPanel(g.detid);
    Grid2dBinGroup & g2bg = panel.call_g2bg();

    // Calculate and set volume
    double vol = panel.calc_approx_volume_sum_grouped(g.igroup);
    int ngrp = get_dic().get_n_group(g.detid);
    g2bg.resize_vec_volume_group(ngrp, g2bg.volume_group_init);
    g2bg.set_volume_group(g.igroup, vol);

    // Get cell center density
    double dens = panel.get_g2bg().get_dens_group_center(g.igroup);
    double mass = vol * dens;

    fprintf(fout,"%2d %4d %6.4E %6.4E %6.4E\n",g.detid,g.igroup,vol,dens,mass);

    sum_vol  += vol;
    sum_mass += mass;
  }
  myapp::close(fout,pathout);

  if(sum_vol == 0.0){
    LOG_WARN("sum_vol == 0.0, returning initial density value.");  
    return Grid2dBinGroup::dens_group_center_init;
  }

  set(FlgProg::tf_calc_volume_group, true);

  if( !std::isfinite(sum_vol) )
    THROW_ERROR("sum_vol is non-finite.");
  
  const double dens_ave = sum_mass / sum_vol;

  if (!std::isfinite(dens_ave))
    THROW_ERROR("dens_ave is non-finite.");

  return dens_ave;
}

/// @brief Calculate volume-weighted average density.
/// @details Get volume and density for each group and return lower, center, and upper bounds of volume-weighted average.
/// @return 3-element array (lower, center, upper density bounds)
std::array<double, 3>
  DetectorPanelArray::calc_volume_weighted_average_density_lower_center_upper(
    const fs::path& prefix_out )
{
  // Check uqid/group construction
  dic_.check_built_uqid();
  check_auto_grouping();

  // Get group list
  const auto& gm = dic_.getGrpMgr();

  double sum_vol  = 0.0;
  double sum_mass_low = 0.0;
  double sum_mass_cnt = 0.0;
  double sum_mass_upp = 0.0;

  const fs::path pathout_low = iodir::make_pathout(prefix_out.string() + "_low.tmp");
  const fs::path pathout_cnt = iodir::make_pathout(prefix_out.string() + "_cnt.tmp");
  const fs::path pathout_upp = iodir::make_pathout(prefix_out.string() + "_upp.tmp");

  FILE *fout_low = myapp::get_fout(pathout_low);
  FILE *fout_cnt = myapp::get_fout(pathout_cnt);
  FILE *fout_upp = myapp::get_fout(pathout_upp);

  // Scan each group
  for (const auto& ginfo : gm.get_vecGroupInfo_all()) {
    if(!ginfo.is_avail) continue; // Skip if group is not available

    // Corresponding panel and group number
    DetectorPanel& panel = callDetectorPanel(ginfo.detid);
    Grid2dBinGroup & g2bg = panel.call_g2bg();

    // Calculate and set volume
    double vol = panel.calc_approx_volume_sum_grouped(ginfo.igroup);

    // check vol is not NaN
    if( std::isnan(vol) )
      THROW_ERROR("DetectorPanelArray::calc_volume_weighted_average_density_lower_center_upper: vol is NaN for detid={}, igroup={}",ginfo.detid,ginfo.igroup);

    int ngrp = get_dic().get_n_group(ginfo.detid);
    g2bg.resize_vec_volume_group(ngrp, g2bg.volume_group_init);
    g2bg.set_volume_group(ginfo.igroup, vol);

    // Get cell center density
    const double dens_low = panel.get_g2bg().get_dens_group_lower(ginfo.igroup);
    const double dens_cnt = panel.get_g2bg().get_dens_group_center(ginfo.igroup);
    const double dens_upp = panel.get_g2bg().get_dens_group_upper(ginfo.igroup);

    // check for NaN sources
    if (!std::isfinite(dens_low)) 
      THROW_ERROR("dens_low is non-finite. detid={}, igroup={}, dens_low={}"
      , ginfo.detid, ginfo.igroup, dens_low);
        
    if (!std::isfinite(dens_cnt))
      THROW_ERROR("dens_cnt is non-finite. detid={}, igroup={}, dens_cnt={}"
      , ginfo.detid, ginfo.igroup, dens_cnt);
    if (!std::isfinite(dens_upp)) 
      THROW_ERROR("dens_upp is non-finite. detid={}, igroup={}, dens_upp={}"
      , ginfo.detid, ginfo.igroup, dens_upp);

    const double mass_low = vol * dens_low;
    const double mass_cnt = vol * dens_cnt;
    const double mass_upp = vol * dens_upp;

    fprintf(fout_low,"%2d %4d %6.4E %6.4E %6.4E\n",ginfo.detid,ginfo.igroup,vol,dens_low,mass_low);
    fprintf(fout_cnt,"%2d %4d %6.4E %6.4E %6.4E\n",ginfo.detid,ginfo.igroup,vol,dens_cnt,mass_cnt);
    fprintf(fout_upp,"%2d %4d %6.4E %6.4E %6.4E\n",ginfo.detid,ginfo.igroup,vol,dens_upp,mass_upp);

    sum_vol  += vol;
    sum_mass_low += mass_low;
    sum_mass_cnt += mass_cnt;
    sum_mass_upp += mass_upp;
  }
  myapp::close(fout_low,pathout_low);
  myapp::close(fout_cnt,pathout_cnt);
  myapp::close(fout_upp,pathout_upp);

  // sum_vol shouldn't be NaN here
  if (!std::isfinite(sum_vol)) THROW_ERROR("sum_vol is non-finite.");

  // when sum_vol == 0.0, return initial density values
  if (sum_vol == 0.0) {
    LOG_WARN("sum_vol == 0.0, return initial density values.");
    return Grid2dBinGroup::arr_dens_group_init;
  }

  const double dens_ave_low = sum_mass_low / sum_vol;
  const double dens_ave_cnt = sum_mass_cnt / sum_vol;
  const double dens_ave_upp = sum_mass_upp / sum_vol;

  return { dens_ave_low, dens_ave_cnt, dens_ave_upp };
}

/// @brief get volume sum of detid
double DetectorPanelArray::get_volume_sum(const Detid detid) const
{
  // check if uqid &uqig maps is built, if not throw error
  dic_.check_built_uqid();
  check_auto_grouping();

  double sum_vol = 0.0;
  // get volume
  const DetectorPanel& panel = getDetectorPanel(detid);
  const std::vector<Igroup> vec_igroup = panel.get_dic().get_vec_igroup(detid);
  for(const int &igroup : vec_igroup){
    double vol = panel.get_g2bg().get_volume_group(igroup);    
    // check vol is not NaN
    sum_vol += vol;
  }
  if( std::isnan(sum_vol) )
    THROW_ERROR("get_volume_sum: sum_vol is NaN for detid={}",detid);
  return sum_vol;
}

// @brief disp volume sum of all DetectorPanel
// @return volume sum of all DetectorPanel
double DetectorPanelArray::disp_volume_sum_all(
  spdlog::level::level_enum spdlog_level ) const
{
  const int n_det = get_n_det();
  double sum_vol = 0.0;
  for(Detid detid=0;detid<n_det;detid++){
    const double vol = get_volume_sum(detid);
    mylogger::g_logger->log(spdlog_level,"disp_volume_sum_all | detid={}, vol={:E}",detid,vol);
    sum_vol += vol;
  }
  mylogger::g_logger->log(spdlog_level,"disp_volume_sum_all | sum_vol={:E}",sum_vol);
  return sum_vol;
}

// set group efficiency triplets as the simple mean of element efficiencies
void DetectorPanelArray::calc_set_eff_group()
{
  // grouping must be complete (group buffers are allocated per grouping result)
  check_auto_grouping();

  if( get(FlgProg::tf_set_eff_group) ){
    LOG_WARN("tf_set_eff_group is already set");
    SLEEP_MSEC(500);
  }

  const int n_det = get_n_det();
  for(Detid detid=0;detid<n_det;detid++){
    DetectorPanel& panel = callDetectorPanel(detid);
    Grid2dBinGroup& g2bg = panel.call_g2bg();
    const int n_group = g2bg.get_n_group();
    for(Igroup igroup=0;igroup<n_group;igroup++){
      // simple mean of the element triplets: (eff_low, eff_cnt, eff_upp)
      const std::array<double,3> arr_eff = panel.calc_eff_mean_group(igroup);
      g2bg.set_eff_low_group(igroup, arr_eff[0]);
      g2bg.set_eff_cnt_group(igroup, arr_eff[1]);
      g2bg.set_eff_upp_group(igroup, arr_eff[2]);
    }
  }

  set(FlgProg::tf_set_eff_group, true);
}

/// @brief set PL and DL for all DetectorElement in all DetectorPanel
void DetectorPanelArray::init_PLDL_all( const double PL_init, const double DL_init)
{
  LOG_INFO("");
  LOG_INFO("Set PL and DL to all DetectorElements");
  const int n_det = get_n_det();
  for(Detid detid=0;detid<n_det;detid++){
    DetectorPanel& panel = callDetectorPanel(detid);
    LOG_DEBUG("detid={}",detid);
    panel.mp_init_PLDL(PL_init, DL_init);
  }
}

/// @brief display the status of this instance
void DetectorPanelArray::display_status(const spdlog::level::level_enum level) const
{
  using namespace mymacro; // for char tf()
  mylogger::g_logger->log(level,"============================================");
  mylogger::g_logger->log(level,"DetectorPanelArray:: name = {}", get_name());
  mylogger::g_logger->log(level,"============================================");

  // Verify that std::bitset size and flag_names size match
  const size_t flag_names_size = FlgProgNames.size();
  static_assert(flag_names_size == static_cast<size_t>(FlgProg::COUNT),
    "Mismatch between FLAG_LIST and flag_names");

  // Output each flag of flag_ to log
  for (size_t i = 0; i < flag_names_size; ++i) {
    mylogger::g_logger->log(level, "{}: {}", FlgProgNames[i], tf( flags_.test(i) ) );
  }

  // Output various states to log
  mylogger::g_logger->log(level,"tf_built_uqid: {}", tf(dic_.is_built_uqid()));
  mylogger::g_logger->log(level,"tf_built_uqig: {}", tf(dic_.is_built_uqig()));
  mylogger::g_logger->log(level,"tf_built_indexer: {}", tf(dic_.is_built_indexer()));
  mylogger::g_logger->log(level,"n_det: {}", get_n_det());
  mylogger::g_logger->log(level,"n_all_element: {}", get_n_all_element());
  mylogger::g_logger->log(level,"uqid_min: {}", get_uqid_min());
  mylogger::g_logger->log(level,"uqid_max: {}", get_uqid_max());
  mylogger::g_logger->log(level,"uqig_min: {}", dic_.get_uqig_min());
  mylogger::g_logger->log(level,"uqig_max: {}", dic_.get_uqig_max());
  mylogger::g_logger->log(level,"uqig_avail_min: {}", dic_.get_uqig_avail_min());
  mylogger::g_logger->log(level,"uqig_avail_max: {}", dic_.get_uqig_avail_max());
  mylogger::g_logger->log(level,"---------------------------------------------");
}
/// @brief display the status of this instance
void DetectorPanelArray::display_status(FILE *fout) const
{
  using namespace mymacro;
  fprintf(fout,"============================================\n");
  fprintf(fout,"DetectorPanelArray:: name = %s\n", get_name().c_str());
  fprintf(fout,"============================================\n");

  // Verify that std::bitset size and flag_names size match
  const size_t flag_names_size = FlgProgNames.size();
  static_assert(flag_names_size == static_cast<size_t>(FlgProg::COUNT),
    "Mismatch between FLAG_LIST and flag_names");

  // Output each flag of flag_ to log
  for (size_t i = 0; i < flag_names_size; ++i) {
    std::fprintf(fout, "%s: %s\n",
                  FlgProgNames[i],
                  flags_.test(i) ? "true" : "false");
  }

  fprintf(fout,"tf_built_uqid: %s\n", tf(dic_.is_built_uqid()));
  fprintf(fout,"tf_built_uqig: %s\n", tf(dic_.is_built_uqig()));
  fprintf(fout,"tf_built_indexer: %s\n", tf(dic_.is_built_indexer()));
  fprintf(fout,"n_det: %d\n", get_n_det());
  fprintf(fout,"n_all_element: %d\n", get_n_all_element());
  fprintf(fout,"uqid_min: %d\n", get_uqid_min());
  fprintf(fout,"uqid_max: %d\n", get_uqid_max());
  fprintf(fout,"uqig_min: %d\n", dic_.get_uqig_min());
  fprintf(fout,"uqig_max: %d\n", dic_.get_uqig_max());
  fprintf(fout,"uqig_avail_min: %d\n", dic_.get_uqig_avail_min());
  fprintf(fout,"uqig_avail_max: %d\n", dic_.get_uqig_avail_max());
  fprintf(fout,"---------------------------------------------\n");
}

void DetectorPanelArray::check_detector_positions_in_grid_xy(
  double xmin, double xmax, double ymin, double ymax,
  double margin_factor) const
{
  // Ensure xmin < xmax and ymin < ymax
  if (xmin > xmax) std::swap(xmin, xmax);
  if (ymin > ymax) std::swap(ymin, ymax);

  // Calculate grid dimensions and expanded range
  const double x_range = xmax - xmin;
  const double y_range = ymax - ymin;
  const double x_margin = x_range * margin_factor;
  const double y_margin = y_range * margin_factor;

  const double x_check_min = xmin - x_margin;
  const double x_check_max = xmax + x_margin;
  const double y_check_min = ymin - y_margin;
  const double y_check_max = ymax + y_margin;

  LOG_INFO("Checking detector positions against grid range");
  LOG_INFO("Grid range: x=[{:.1f}, {:.1f}], y=[{:.1f}, {:.1f}]", xmin, xmax, ymin, ymax);
  LOG_INFO("Check range (margin_factor={:.1f}): x=[{:.1f}, {:.1f}], y=[{:.1f}, {:.1f}]",
           margin_factor, x_check_min, x_check_max, y_check_min, y_check_max);

  const int n_det = get_n_det();
  std::vector<std::string> errors;

  for (Detid detid = 0; detid < n_det; ++detid) {
    const DetectorPanel& panel = getDetectorPanel(detid);
    const Eigen::Vector3d pos = panel.get_v3_pos();
    const double det_x = pos.x();
    const double det_y = pos.y();
    const double det_z = pos.z();

    LOG_DEBUG("Detector {}: name={}, position=({:.1f}, {:.1f}, {:.1f})",
              detid, panel.get_name(), det_x, det_y, det_z);

    bool x_ok = (det_x >= x_check_min && det_x <= x_check_max);
    bool y_ok = (det_y >= y_check_min && det_y <= y_check_max);

    if (!x_ok || !y_ok) {
      std::string msg = fmt::format(
        "Detector {} '{}' at ({:.1f}, {:.1f}, {:.1f}) is outside grid range. "
        "Grid: x=[{:.1f}, {:.1f}], y=[{:.1f}, {:.1f}]. "
        "Distance from grid center: dx={:.1f}, dy={:.1f}",
        detid, panel.get_name(), det_x, det_y, det_z,
        xmin, xmax, ymin, ymax,
        det_x - (xmin + xmax) / 2.0,
        det_y - (ymin + ymax) / 2.0);
      errors.push_back(msg);
      LOG_ERROR("{}", msg);
    }
  }

  if (!errors.empty()) {
    std::string full_msg = fmt::format(
      "DetectorPanelArray::check_detector_positions_in_grid_xy: {} of {} detectors are outside the grid range. "
      "This likely indicates a coordinate system mismatch between detector parameters and grid parameters.",
      errors.size(), n_det);
    THROW_ERROR("{}", full_msg);
  }

  LOG_INFO("All {} detectors are within the expected grid range", n_det);
}

/// @brief set the exposure time for all DetectorElement
void DetectorPanelArray::set_exposure_time_sec( const double time_sec_in )
{
  LOG_INFO("");
  LOG_INFO("Set exposure_time_sec to all DetectorElements");
  const int n_det = get_n_det();
  for(Detid detid=0;detid<n_det;detid++){
    DetectorPanel& panel = callDetectorPanel(detid);
    LOG_DEBUG("detid={}",detid);
    panel.set_exposure_time_sec(time_sec_in);
  }
}

// Assign penetrating_muonflux and signal muon count to each DetectorElement
// in vec_vec_DetectorElement
// openmp version
void DetectorPanelArray::mp_calc_set_peneflux_signal_from_DL(
  const FluxTable &ft, const bool tf_apply_eff
, const bool tf_apply_eff_cnt )
{
  LOG_INFO("");
  LOG_INFO("Assign penetrating_muonflux and signal muon count to each DetectorElement in vec_vec_DetectorElement");
  if( !get(FlgProg::tf_calc_PL_DL) ){
    LOG_ERROR("PL and DL are not calculated yet, call mp_calc_PL_DL_all() first.");
    THROW_ERROR_NAME("PL and DL are not calculated yet.");
  }
  const int n_det = get_n_det();
  for(Detid detid=0;detid<n_det;detid++){
    DetectorPanel& panel = callDetectorPanel(detid);
    LOG_DEBUG("detid={}",detid);
    panel.mp_calc_set_peneflux_signal_from_DL(ft,tf_apply_eff,tf_apply_eff_cnt);
  }
  set(FlgProg::tf_calc_signal, true);
}

// @brief Assign noise to each DetectorElement in vec_vec_DetectorElement
// @details openmp version
void DetectorPanelArray::mp_set_noise_all( const NoiseParameters &ndist, const double DL_thres)
{
  if(ndist.get_tf_exec() == false){
    LOG_INFO("");
    LOG_INFO("ndist.get_tf_exec() == false, skip setting noise");
    set(FlgProg::tf_calc_noise, true);
    return;
  }

  if( get(FlgProg::tf_calc_noise) ){
    LOG_WARN("tf_calc_noise is true");
    SLEEP_MSEC(500);
  }

  for(Detid detid=0;detid<get_n_det();detid++){
    DetectorPanel& panel = callDetectorPanel(detid);
    panel.mp_set_noise( ndist.get_flux_proport_ratio_floor()
      , ndist.get_flux_proport_ratio_poisson()
      , ndist.get_SOT_proport_noise_ratio_floor()
      , ndist.get_SOT_proport_noise_ratio_poisson()
      , ndist.get_user_defined_noise_flux_ratio()
      , ndist.get_vec_path_user_defined_noise_flux()
      , DL_thres );
  }
  set(FlgProg::tf_calc_noise, true);
  display_status(spdlog::level::info);
  LOG_INFO("done");
}

// Execute DetectorPanel::copy_signal_noise_to_g2bg for all DetectorPanels
void DetectorPanelArray::mp_copy_signal_noise_to_g2bg_all(
  const double signal_init, const double noise_init, const bool is_avail_init
, const double PL_thres, const double DL_thres, const double signal_under_thres, const double noise_under_thres
, const bool is_avail_under_thres )
{
  LOG_INFO("");
  LOG_INFO("Execute DetectorPanel::copy_signal_noise_to_g2bg for all DetectorPanels");
  const int n_det = get_n_det();
  // #pragma omp parallel for
  for(Detid detid=0;detid<n_det;detid++){
    LOG_DEBUG("detid={}",detid);
    DetectorPanel& panel = callDetectorPanel(detid);
    panel.copy_signal_noise_to_g2bg(
      signal_init, noise_init, is_avail_init
      , PL_thres, DL_thres
      , signal_under_thres, noise_under_thres, is_avail_under_thres);
  }
}

void DetectorPanelArray::mp_assign_1st_igroup_all(
  const bool tf_run_1st_grouping
, const int igroup_start
, const int nx_div_init, const int ny_div_init
, const double signal_noise_thres )
{
  const int n_det = get_n_det();
  const int n_group_1st = nx_div_init * ny_div_init;
  if( tf_run_1st_grouping ){
    LOG_INFO("Execute DetectorPanel::assign_1st_igroup for all detectors");
    #pragma omp parallel for
    for(Detid detid=0;detid<n_det;detid++){
      DetectorPanel& panel = callDetectorPanel(detid);
      LOG_DEBUG("assign_1st_igroup detid={}",detid);
      panel.assign_1st_igroup(igroup_start,nx_div_init,ny_div_init);
      // check min statistics per group
      bool tf_all_larger = 
        panel.get_g2bg().is_all_signal_noise_group_larger_than_thres(signal_noise_thres);
      if( ! tf_all_larger ){
        LOG_ERROR(
          "detid={} has group(s) with signal/noise < {} after assign_1st_igroup",
          detid, signal_noise_thres);
        THROW_ERROR_NAME("Some group has insufficient statistics after assign_1st_igroup");
      }
    }
  }else{
    // if tf_run_auto_grouping==false
    LOG_INFO("Execute DetectorPanel::assign_naive_igroup for all detectors");
    #pragma omp parallel for
    for(Detid detid=0;detid<n_det;detid++){
      DetectorPanel& panel = callDetectorPanel(detid);
      panel.assign_naive_igroup(igroup_start);
    }
  }
}

// @brief Copy DetectorPanel::vec_tf_in_PL for all DetectorPanels
void DetectorPanelArray::mp_copy_vec_tf_in_PL_all( const DetectorPanelArray &arrdet_in )
{
  for(int det_id=0;det_id<get_n_det();det_id++){
    const DetectorPanel& panel_in = arrdet_in.getDetectorPanel(det_id);
    DetectorPanel& panel = callDetectorPanel(det_id);
    panel.mp_set_vec_tf_in_PL( panel_in ); 
  }
}

// @brief Execute DetectorPanel::grouping_by_bin_list for all DetectorPanels
// @param prm_bingrp
void DetectorPanelArray::mp_grouping_by_bin_list_all(
  const Grid2dBinGroup::Parameters &prm_bingrp )
{
  if( get(FlgProg::tf_auto_grouping) ){
    LOG_WARN(
      "tf_auto_grouping_done is true, do nothing.");
    SLEEP_MSEC(500);
    return;
  }

  LOG_INFO("");
  LOG_INFO("Execute DetectorPanel::grouping_by_bin_list for all DetectorPanels");
  const int n_det = get_n_det();

  // NOTE: dic_ must NOT be re-initialized here. The uqid index built at panel
  // construction is required by the subsequent build_index_container() call.
  if( static_cast<int>(prm_bingrp.vec_tf_read_bin_group_list.size()) != n_det )
    THROW_ERROR(
      "DetectorPanelArray::mp_grouping_by_bin_list_all: vec_tf_read_bin_group_list.size()={} != n_det={}"
      , prm_bingrp.vec_tf_read_bin_group_list.size(), n_det);
  if( static_cast<int>(prm_bingrp.vec_file_path_bin_group_list.size()) != n_det )
    THROW_ERROR(
      "DetectorPanelArray::mp_grouping_by_bin_list_all: vec_file_path_bin_group_list.size()={} != n_det={}"
      , prm_bingrp.vec_file_path_bin_group_list.size(), n_det);
  int n_manual = 0; // number of detectors with a manual bin-group list
  for(Detid detid=0;detid<n_det;detid++)
    if( prm_bingrp.vec_tf_read_bin_group_list.at(detid) ) n_manual++;
  if( n_manual != prm_bingrp.n_detector_grouping_manual )
    LOG_WARN("n_detector_grouping_manual={} != number of true flags={}"
      , prm_bingrp.n_detector_grouping_manual, n_manual);

  #pragma omp parallel for
  for(Detid detid=0;detid<n_det;detid++){
    if( ! prm_bingrp.vec_tf_read_bin_group_list.at(detid) ){
      LOG_DEBUG("detid={} skip",detid);
      continue;
    }
    DetectorPanel& panel = callDetectorPanel(detid);
    panel.grouping_by_bin_list(prm_bingrp.vec_file_path_bin_group_list.at(detid));
  }
  set(FlgProg::tf_auto_grouping, true);
}



// Execute DetectorPanel::auto_divide_by_signal_noise_group_all for all DetectorPanels.
void DetectorPanelArray::mp_auto_grouping_by_signal_noise_group_alldet(
  const double signal_group_trg, const int nloop_limit,
  const int ixlen_min, const int iylen_min, const bool tf_prefer_split_x )
{
  LOG_INFO("");
  LOG_INFO("Execute DetectorPanel::auto_divide_by_signal_noise_group_all for all DetectorPanels");

  // if tf_auto_grouping_done is true, display warning.
  if ( get(FlgProg::tf_auto_grouping) ){
    LOG_WARN("tf_auto_grouping_done is already true");
    SLEEP_MSEC(500);
  }

  const int n_det = get_n_det();

  #pragma omp parallel for
  for(Detid detid=0;detid<n_det;detid++){
    DetectorPanel& panel = callDetectorPanel(detid);
    Grid2dBinGroup & g2bg = panel.call_g2bg();
    LOG_DEBUG(
      "thread_id={}, Grind2d::auto_divide_by_signal_noise_group_all detid={}"
      ,omp_get_thread_num(),detid);
    g2bg.auto_divide_by_signal_noise_group_all(
      signal_group_trg,nloop_limit,ixlen_min,iylen_min,tf_prefer_split_x);

    // char cfilename[256];
    // sprintf(cfilename, "bimap_final_det%02d.tmp", detid);
    // fs::path pathout = iodir::make_pathout(cfilename);
    // g2bg.write_bimap_header_to_ascii(pathout);
  }
  // set tf_auto_grouping_done to true
  set(FlgProg::tf_auto_grouping, true);

}

// Check that every group of every DetectorPanel satisfies ixlen >= ixlen_min
// and iylen >= iylen_min, regardless of which grouping route built the groups.
void DetectorPanelArray::mp_check_group_ixiylen_min_all(
  const int ixlen_min, const int iylen_min ) const
{
  LOG_INFO("");
  LOG_INFO("Check that no group is smaller than ixlen_min={} x iylen_min={} for all DetectorPanels"
    , ixlen_min, iylen_min);
  const int n_det = get_n_det();
  int n_det_violate = 0;
  #pragma omp parallel for reduction(+:n_det_violate)
  for(Detid detid=0;detid<n_det;detid++){
    const DetectorPanel& panel = getDetectorPanel(detid);
    // record violations only; throwing out of an OpenMP loop is not allowed
    if( ! panel.get_g2bg().is_all_group_ixiylen_larger_than_min(ixlen_min,iylen_min) )
      n_det_violate++;
  }
  if( n_det_violate > 0 )
    THROW_ERROR(
      "DetectorPanelArray::mp_check_group_ixiylen_min_all: {} detector(s) have group(s) "
      "smaller than ixlen_min={} x iylen_min={}. See LOG_ERROR lines above."
      , n_det_violate, ixlen_min, iylen_min);
}

/// @brief Execute mp_calc_set_proj_dens_grouped for all DetectorPanels
void DetectorPanelArray::mp_calc_set_proj_dens_grouped_all(
  const bool tf_signal_poisson
, const double dens_min, const double dens_max
, const double dens_step, const double sigma
, const FluxTable &ft
, const bool tf_eff, const bool tf_eff_independent )
{
  assert(dens_min >= 0.0);
  assert(dens_min < dens_max);
  assert(dens_step > 0.0);

  LOG_INFO("...");
  const int n_det = get_n_det();

  #pragma omp parallel for
  for(Detid detid=0;detid<n_det;detid++){
    DetectorPanel& panel = callDetectorPanel(detid);
    panel.mp_calc_set_proj_dens_grouped(
      tf_signal_poisson,dens_min,dens_max,dens_step,sigma,ft,tf_eff,tf_eff_independent);
    // vec_vec_ret.at(detid).resize(vec_ret.size());
  }
  LOG_INFO("done");
}

void DetectorPanelArray::mp_calc_set_proj_dens_grouped_all(
  const nlohmann::json &js, const FluxTable &ft
, const bool tf_eff, const bool tf_eff_independent )
{
  const bool tf_exec = js["tf_exec"];
  if( !tf_exec ){
    LOG_INFO("tf_exec=false");
    return;
  }
  const bool tf_signal_poisson = js.value("tf_signal_poisson", false);
  const double dens_min = js["dens_min"];
  const double dens_max = js["dens_max"];
  const std::vector<double> vec_dens_step = js["dens_steps"].get<std::vector<double>>();
  const double range_factor = js["range_factor"];
  const double sigma = js["sigma"];
  LOG_DEBUG("tf_exec={}, tf_signal_poisson={}, dens_min={}, dens_max={}, vec_dens_step.size()={}, range_factor={}, sigma={}, tf_eff={}, tf_eff_independent={}",
    tf_exec, tf_signal_poisson, dens_min, dens_max, vec_dens_step.size(), range_factor, sigma, tf_eff, tf_eff_independent
  );
  return mp_calc_set_proj_dens_grouped_all(
      tf_signal_poisson, dens_min, dens_max
    , vec_dens_step, range_factor, sigma, ft, tf_eff, tf_eff_independent );
}

/// @brief Execute mp_calc_set_proj_dens_grouped for all DetectorPanels
/// @note adaptive search version
void DetectorPanelArray::mp_calc_set_proj_dens_grouped_all(
  const bool tf_signal_poisson
, const double dens_min, const double dens_max
, const std::vector<double> &vec_dens_step, const double range_factor
, const double sigma, const FluxTable &ft
, const bool tf_eff, const bool tf_eff_independent )
{
  assert(dens_min >= 0.0);
  assert(dens_min < dens_max);
  if( vec_dens_step.size() > 1 ){
    LOG_INFO("vec_dens_step.size() > 1, adaptive search version");
    // check vec_dens_step is all positive
    for(const double dens_step : vec_dens_step){
      if( dens_step <= 0.0 ) THROW_ERROR("dens_step <= 0.0");
    }
    // check range_factor is positive
    if( range_factor <= 1.0 ) THROW_ERROR("range_factor <= 1.0");
  }else{
      LOG_INFO("vec_dens_step.size() == 1, uniform step search version");
      // check dens_step is positive
      const double dens_step = vec_dens_step.at(0);
      if( dens_step <= 0.0 ) THROW_ERROR("dens_step <= 0.0");
  }

  LOG_INFO("...");
  const int n_det = get_n_det();

  #pragma omp parallel for
  for(Detid detid=0;detid<n_det;detid++){
    DetectorPanel& panel = callDetectorPanel(detid);
    if(vec_dens_step.size() > 1){
      panel.mp_calc_set_proj_dens_grouped(
        tf_signal_poisson,dens_min,dens_max
        , vec_dens_step,range_factor,sigma,ft,tf_eff,tf_eff_independent);
    }else{
      panel.mp_calc_set_proj_dens_grouped(
        tf_signal_poisson,dens_min,dens_max
        , vec_dens_step.at(0),sigma,ft,tf_eff,tf_eff_independent);
    }
  }
  set(FlgProg::tf_calc_proj_dens_group, true);
  LOG_INFO("done");
}


/// @brief Calculate and set proj_density for all DetectorElements
/// @param PL_thres  Projection threshold
/// @param DL_thres  Distance threshold
/// @throws std::out_of_range When dic_ is not built
void DetectorPanelArray::mp_calc_set_proj_density_all(
  const double PL_thres , const double DL_thres)
{
  // Check uqid range
  dic_.check_built_uqid();

  // Get all uqid list
  const std::vector<Uqid> vec_Uqid = dic_.getUqidMgr().get_vecUqid_all();
  const size_t n = vec_Uqid.size();

  #pragma omp parallel for schedule(static)
  for (size_t i = 0; i < n; ++i) {
    // Get uqid→DetectorElement* via call wrapper
    DetectorElement& ele = callDetectorElement(vec_Uqid[i]);
    // Calculate and set proj_density
    ele.calc_set_proj_density(PL_thres, DL_thres);
  }
}

std::map<Detid, FILE*>
  DetectorPanelArray::get_outfile_map_and_out_header(
    const std::string& prefix, const std::string& suffix) const
{
  // Check if UqidManager is built
  dic_.check_built_uqid();

  const int n_det = get_n_det();
  if (n_det <= 0) {
    THROW_ERROR_NAME("No detectors found in DetectorPanelArray");
  }

  // Generate file map
  std::map<int, FILE*> file_map;
  char fname[512];

  // Generate file_map for each detector and output header
  for (Detid detid = 0; detid < n_det; ++detid) {
    // Generate file name
    std::snprintf(
      fname, sizeof(fname),
      "%s_%s_det%02d%s.tmp",
      prefix.c_str(), name.c_str(), detid, suffix.c_str());

    auto file_path  = fs::path(fname);
    auto full_path  = iodir::make_pathout(file_path);

    // Open file
    FILE* fout = std::fopen(full_path.string().c_str(), "wt");
    if (!fout) {
      THROW_ERROR_NAME("Cannot open file: " + full_path.string());
    }
    file_map[detid] = fout;

    // Output header
    const DetectorPanel& panel = getDetectorPanel(detid);
    panel.get_g2bg().out_g2_header(fout);
  }

  return file_map;
}

void DetectorPanelArray::close_g2bg_outfile_map(std::map<Detid, FILE*>& file_map) const
{
  for (auto& [detid, fout] : file_map) {
    if (fout != nullptr) {
      fclose(fout);
      fout = nullptr; // Avoid dangling pointer by setting to nullptr
    }
  }
  file_map.clear(); // Clear the map
}

/// @brief Output signal_noise_group etc. to file for each det_id for all uqig_avail
void DetectorPanelArray::out_g2bg_all() const
{
  // Prerequisite check
  dic_.check_built_uqid();
  check_auto_grouping();
  dic_.check_built_uqig();
  dic_.check_built_uqigAvail();

  // Generate file map and output header
  auto file_map = get_outfile_map_and_out_header("g2bg", "");

  // Get each Manager/Indexer
  const auto& gm = dic_.getGrpMgr();
  const auto& ai = dic_.getAvailIndexer();

  // Scan all groups
  for (const auto& g : gm.get_vecGroupInfo_all()) {
    // Index index = ai.getIndex(g.uqig_avail);
    Detid detid = g.detid;
    Igroup igroup = g.igroup;

    // Output file
    FILE* fout = file_map.at(detid);
    if (!fout)
      THROW_ERROR_NAME("File pointer for detid " + std::to_string(detid) + " is nullptr");

    // Get panel
    const DetectorPanel& panel = getDetectorPanel(detid);
    const Grid2dBinGroup& g2bg = panel.get_g2bg();

    // Coordinate boundaries
    auto [xmin, xmax, ymin, ymax] = g2bg.get_xmin_xmax_ymin_ymax(igroup);

    // Signal noise
    double sg  = g2bg.calc_signal_group(igroup);
    double ng  = g2bg.calc_noise_poi_group(igroup);
    double sng = sg + ng;
    bool   ok  = g2bg.is_avail_group(igroup);

    // Density lower/center/upper
    double dgl = g2bg.get_dens_group_lower(igroup);
    double dgc = g2bg.get_dens_group_center(igroup);
    double dgu = g2bg.get_dens_group_upper(igroup);

    // Δnmuon group
    double dnl = g2bg.get_delta_nmuon_group_lower(igroup);
    double dnc = g2bg.get_delta_nmuon_group_center(igroup);
    double dnu = g2bg.get_delta_nmuon_group_upper(igroup);

    // Volume and efficiency
    double vol = g2bg.get_volume_group(igroup);
    double el  = g2bg.get_eff_low_group(igroup);
    double ec  = g2bg.get_eff_cnt_group(igroup);
    double eu  = g2bg.get_eff_upp_group(igroup);

    // Output (following original format)
    fprintf(fout,
      "%7d %7.4lf %7.4lf %7.4lf %7.4lf "
      "%E %E %E "
      "%6.0lf %6.0lf %6.0lf "
      "%E %E %E "
      "%E %7.4lf %7.4lf %7.4lf %d\n"
      , g.uqig_avail, xmin, xmax, ymin, ymax
      , sg, ng,  sng
      , dgl, dgc, dgu
      , dnl, dnc, dnu
      , vol, el,  ec,  eu, static_cast<int>(ok)
    );
  }

  // Close all files
  close_g2bg_outfile_map(file_map);
}

/// @brief Call get_g2bg().out2 for all DetectorPanels and output vecxf
/// @param vecxf         Vector holding values to output (row = uqig_avail)
/// @param vecxf_name    Vector name also used for file name prefix
/// @throws std::out_of_range When dic_ is not built or grouping is incomplete
void DetectorPanelArray::out_g2bg_all(
  const Eigen::VectorXf& vecxf, const std::string& vecxf_name) const
{
  // Check uqid & group construction
  dic_.check_built_uqid();
  check_auto_grouping();

  auto file_map = get_outfile_map_and_out_header(vecxf_name);

  const auto& gm = dic_.getGrpMgr();
  const auto& ai = dic_.getAvailIndexer();

  for (const auto& gi : gm.get_vecGroupInfo_all()) {
    const Detid detid     = gi.detid;
    const Igroup igroup     = gi.igroup;
    const int uqig_avail = gi.uqig_avail;

    FILE* fout = file_map.at(detid);
    if (!fout)
      THROW_ERROR_NAME("File pointer for detid: " + std::to_string(detid) + " is null");

    const DetectorPanel& panel = getDetectorPanel(detid);
    const Grid2dBinGroup& g2bg = panel.get_g2bg();
    auto [xmin, xmax, ymin, ymax] = g2bg.get_xmin_xmax_ymin_ymax(igroup);
    const bool is_avail = g2bg.is_avail_group(igroup);

    // Retrieve vecxf row using avail_indexer
    const int row = ai.getIndex(uqig_avail);
    const double value = vecxf(row);

    // Output one line
    fprintf(fout,
      "%7d %7.4lf %7.4lf %7.4lf %7.4lf %E %d\n",
      uqig_avail,
      xmin, xmax, ymin, ymax,
      value,
      (int)is_avail
    );
  }

  close_g2bg_outfile_map(file_map);
}

void DetectorPanelArray::out_txtyData( const std::string& suffix
  , std::function<double(const DetectorElement&)> data_extractor) const
{
  // Check uqid range and grouping construction
  dic_.check_built_uqid();
  check_auto_grouping();

  // Get detid list
  std::set<Detid> set_detid = dic_.get_set_detid();
  if( set_detid.size() == 0 ){
    LOG_ERROR(
      "No detid found in dic_, please check you ran DetectorPanelArray::build_index_container()");
    THROW_ERROR_NAME("No detid found in dic_");
  }
  // Open file for each detid
  std::map<int, FILE*> file_map;
  char fname[512];
  for (Detid detid : set_detid) {
    std::snprintf(
      fname, sizeof(fname),
      "%s_txty%s_det%02d.tmp",
      get_name().c_str(), suffix.c_str(), detid);

    const fs::path pathout = iodir::make_pathout(fname);
    FILE* fout = std::fopen(pathout.c_str(), "wt");
    if (!fout) THROW_ERROR_NAME("Cannot open file: " + pathout.string());
    file_map[detid] = fout;
  }

  // Scan all uqid
  const auto uqids = dic_.getUqidMgr().get_vecUqid_all();
  for (Uqid u : uqids) {
    // uqid → detid, ix, iy
    DetIxiy di = dic_.get_detidixiy(u);

    // Get element
    const DetectorElement& ele = getDetectorElement(di);
    double tx    = ele.get_tx();
    double ty    = ele.get_ty();
    double data  = data_extractor(ele);

    // Output
    FILE* fout = file_map.at(di[0]);
    fprintf(fout, "%7.4lf %7.4lf %E\n", tx, ty, data);
  }

  // Close all files
  for (auto& [detid, fout] : file_map) {
    if (std::fclose(fout) == EOF)
      THROW_ERROR_NAME("fclose failed for detid: " + std::to_string(detid));
  }
}

void DetectorPanelArray::mp_calc_vec_signal_noise_group_all()
{
  if( get(FlgProg::tf_calc_signal_group) ){
    LOG_WARN("tf_calc_signal_group is already true, do nothing.");
    SLEEP_MSEC(500);
  }
  if( get(FlgProg::tf_calc_noise_group) ){
    LOG_WARN("tf_calc_noise_group is already true, do nothing.");
    SLEEP_MSEC(500);
  }

  LOG_INFO("Copy signal/noise_poi_group to vec_signal/noise_poi_group for all DetectorPanels");
  const int n_det = get_n_det();
  // #pragma omp parallel for // ! calc_set_vec_signal_noise_group uses openmp.
  for(Detid detid=0;detid<n_det;detid++){
    DetectorPanel& panel = callDetectorPanel(detid);
    LOG_DEBUG(
      "det_name={}, detid={}"
      ,panel.get_name(),detid);

    panel.call_g2bg().calc_set_vec_signal_noise_group();
  }
  set(FlgProg::tf_calc_signal_group, true);
  set(FlgProg::tf_calc_noise_group, true);
}

void DetectorPanelArray::mp_set_signal_noise_group_poisson_all()
{
  LOG_INFO("Assign different Poisson errors to signal/noise_poi_group for all DetectorPanels");
  const int n_det = get_n_det();
  #ifdef NODEBUG
    #pragma omp parallel for
  #else
    // do nothing
  #endif
  for(Detid detid=0;detid<n_det;detid++){
    DetectorPanel& panel = callDetectorPanel(detid);
    LOG_DEBUG(
      "det_name={}, detid={}"
      ,panel.get_name(),detid);

    panel.call_g2bg().set_diff_poisson_again();
  }
}

DetectorPanelArray DetectorPanelArray::get_subtract(
  const DetectorPanelArray &other) const
{
  // Check if dic_ of self and other are properly built
  dic_.check_built_uqid();
  other.dic_.check_built_uqid();

  // Copy this and create result
  DetectorPanelArray result(*this);
  result.name = this->get_name() + "_subtract_" + other.get_name();

  // Also copy DetectorIndexContainer
  result.set_dic(other.get_dic());

  // Clear old panels and fill with new difference panels
  result.vec_panel.clear();
  for (const auto &panel : vec_panel) {
    Detid detid = panel.get_detid();
    const DetectorPanel& otherPanel = other.getDetectorPanel(detid);
    // Execute subtract for each panel
    DetectorPanel sub = panel.get_subtract(otherPanel);
    result.vec_panel.push_back(std::move(sub));
  }

  return result;
}

void DetectorPanelArray::allocate_vec_value_group_all()
{
  const int n_det = get_n_det();
  LOG_INFO("Execute Grid2dBinGroup::allocate_vec_value_group for all DetectorPanels.");
  #pragma omp parallel for
  for(Detid detid=0; detid<n_det; detid++){
    DetectorPanel& panel = callDetectorPanel(detid);
    panel.call_g2bg().allocate_vec_value_group();
  }
}

void DetectorPanelArray::set_dic_bimap_all_from( const DetectorPanelArray& src )
{
  // english Copy the entire index container from src to this.
  set_dic(src.get_dic());

  // Copy bimap_ of each DetectorPanel from src.
  mp_copy_bimap_all(src);

  // Allocate memory size for vec_xxxx_group
  allocate_vec_value_group_all();

  // group buffers were re-initialized, so group efficiency must be recomputed
  set(FlgProg::tf_set_eff_group, false);
}

void DetectorPanelArray::copy_g2bg(
  const Detid det_id, const DetectorPanel& panel_in )
{
  DetectorPanel& panel = callDetectorPanel(det_id);
  panel.call_g2bg().copy_g2bg( panel_in.get_g2bg() );
}

void DetectorPanelArray::mp_copy_g2bg_all( const DetectorPanelArray &arrdet_in )
{
  LOG_INFO("Copy g2bg_ for all DetectorPanels.");
  const int n_det = get_n_det();
  #pragma omp parallel for
  for(Detid detid=0;detid<n_det;detid++){
    DetectorPanel& panel = callDetectorPanel(detid);
    const DetectorPanel& panel_in = arrdet_in.getDetectorPanel(detid);
    panel.copy_g2bg(panel_in);
  }
}

void DetectorPanelArray::mp_copy_bimap_all( const DetectorPanelArray &arrdet_src )
{
  LOG_INFO("Copy bimap_ for all DetectorPanels.");
  const int n_det = get_n_det();
  #pragma omp parallel for
  for(Detid detid=0;detid<n_det;detid++){
    DetectorPanel& panel = callDetectorPanel(detid);
    const DetectorPanel& panel_src = arrdet_src.getDetectorPanel(detid);
    panel.copy_bimap_from(panel_src);
  }
  const bool tf_auto_grouping_src = arrdet_src.get(FlgProg::tf_auto_grouping);
  LOG_DEBUG("tf_auto_grouping_src={}", tf_auto_grouping_src ? "true" : "false");
  set(FlgProg::tf_auto_grouping, tf_auto_grouping_src);
}

void DetectorPanelArray::mp_init_PL_DL_pene_sig_noi_all()
{
  LOG_INFO("Execute DetectorPanel::init_len_DL_pene_sig_noi for all DetectorPanels.");
  const int n_det = get_n_det();
  #pragma omp parallel for
  for(Detid detid=0; detid<n_det; detid++){
    callDetectorPanel(detid).mp_init_PL_DL_pene_sig_noi();
  }
  set(FlgProg::tf_calc_PL_DL, false);
  set(FlgProg::tf_calc_signal, false);
  set(FlgProg::tf_calc_noise, false);
}

// write vec_parameter_file_path to std::ofstream &ofs with binary mode
void DetectorPanelArray::write_vec_parameter_file_path( std::ofstream &ofs ) const
{
  LOG_INFO("write vec_parameter_file_path to std::ofstream &ofs with binary mode");
  const std::size_t n_path = vec_parameter_file_path.size();
  io_binary::write_binary(ofs,n_path);
  for(const auto& path : vec_parameter_file_path ) io_binary::write_path(ofs,path);
  if( ofs.fail() ) THROW_ERROR_NAME("ofs.fail()");
}

// read vec_parameter_file_path from std::ifstream &ifs with binary mode
void DetectorPanelArray::read_vec_parameter_file_path( std::ifstream &ifs )
{
  LOG_INFO("read vec_parameter_file_path from std::ifstream &ifs with binary mode");
  const std::size_t n_path = io_binary::read_binary<std::size_t>(ifs);
  vec_parameter_file_path.resize(n_path);
  for(auto& path : vec_parameter_file_path ) path = io_binary::read_path(ifs);
  if( ifs.fail() ) THROW_ERROR_NAME("ifs.fail()");
}

// write vec_panel to std::ofstream &ofs with binary mode
void DetectorPanelArray::write_vec_panel( std::ofstream &ofs ) const
{
  LOG_INFO("write vec_panel to std::ofstream &ofs with binary mode");
  const std::size_t n_det = vec_panel.size();
  io_binary::write_binary(ofs,n_det);
  for(const auto& panel : vec_panel ) panel.save(ofs);
  if( ofs.fail() ) THROW_ERROR_NAME("ofs.fail()");
}

// read vec_panel from std::ifstream &ifs with binary mode
void DetectorPanelArray::read_vec_panel( std::ifstream &ifs )
{
  LOG_INFO("read vec_panel from std::ifstream &ifs with binary mode");
  const std::size_t n_det = io_binary::read_binary<std::size_t>(ifs);
  vec_panel.resize(n_det);
  for(auto& panel : vec_panel ) panel.load(ifs);
  if( ifs.fail() ) THROW_ERROR_NAME("ifs.fail()");
}

// write all private members to std::ofstream &ofs with binary mode
void DetectorPanelArray::save(std::ofstream &ofs) const
{
  LOG_INFO("- write all private members to binary file");

  io_binary::write_string(ofs, name);
  io_binary::write_binary(ofs, n_all_element);
  write_vec_parameter_file_path(ofs);
  write_vec_panel(ofs);
  dic_.save(ofs);

  // Convert bitset to unsigned long and write
  unsigned long bits = flags_.to_ulong();
  io_binary::write_binary(ofs, bits);

  io_binary::write_set_int(ofs, set_uqiv);

  if (ofs.fail()) THROW_ERROR_NAME("ofs.fail() in DetectorPanelArray::save");
}

// read all private members from std::ifstream &ifs with binary mode
void DetectorPanelArray::load(std::ifstream &ifs)
{
  LOG_INFO("- read all private members from binary file");

  name = io_binary::read_string(ifs);
  n_all_element = io_binary::read_binary<int>(ifs);
  read_vec_parameter_file_path(ifs);
  read_vec_panel(ifs);
  dic_.load(ifs);

  // Restore bit pattern to flag group
  unsigned long bits = io_binary::read_binary<unsigned long>(ifs);
  flags_ = std::bitset<static_cast<size_t>(FlgProg::COUNT)>(bits);

  set_uqiv = io_binary::read_set_int(ifs);

  if (ifs.fail()) THROW_ERROR_NAME("ifs.fail() in DetectorPanelArray::load");
}


//##################################################################################
// @struct DetectorPanelArray::PillarBuildParams
// @brief Structure collecting parameters needed when calculating path length, signal count, etc. of Grid2dPillar in DetectorPanelArray
//##################################################################################

/// @brief Inequality operator
bool DetectorPanelArray::PillarBuildParams::operator!=(const PillarBuildParams& other) const {
  #ifdef NODEBUG
  if (tf_load_arrdet_g2pil != other.tf_load_arrdet_g2pil) return true;
  if (tf_save_arrdet_g2pil != other.tf_save_arrdet_g2pil) return true;
  if (tf_run_1st_grouping != other.tf_run_1st_grouping) return true;
  if (tf_run_auto_grouping != other.tf_run_auto_grouping) return true;
  if (path_arrdet_bin != other.path_arrdet_bin) return true;
  if (arrdet_built != other.arrdet_built) return true;
  if (prm_bingroup != other.prm_bingroup) return true;
  if (g2pil != other.g2pil) return true;
  if (ft != other.ft) return true;
  if (tf_apply_eff != other.tf_apply_eff) return true;
  if (prm_noise != other.prm_noise) return true;
  #else
  if (tf_load_arrdet_g2pil != other.tf_load_arrdet_g2pil) { LOG_WARN("PillarBuildParams: tf_load_arrdet_g2pil differs"); return true; }
  if (tf_save_arrdet_g2pil != other.tf_save_arrdet_g2pil) { LOG_WARN("PillarBuildParams: tf_save_arrdet_g2pil differs"); return true; }
  if (tf_run_1st_grouping != other.tf_run_1st_grouping) { LOG_WARN("PillarBuildParams: tf_run_1st_grouping differs"); return true; }
  if (tf_run_auto_grouping != other.tf_run_auto_grouping) { LOG_WARN("PillarBuildParams: tf_run_auto_grouping differs"); return true; }
  if (path_arrdet_bin != other.path_arrdet_bin) { LOG_WARN("PillarBuildParams: path_arrdet_bin differs"); return true; }
  if (arrdet_built != other.arrdet_built) { LOG_WARN("PillarBuildParams: arrdet_built differs"); return true; }
  if (prm_bingroup != other.prm_bingroup) { LOG_WARN("PillarBuildParams: prm_bingroup differs"); return true; }
  if (g2pil != other.g2pil) { LOG_WARN("PillarBuildParams: g2pil differs"); return true; }
  if (ft != other.ft) { LOG_WARN("PillarBuildParams: ft differs"); return true; }
  if (tf_apply_eff != other.tf_apply_eff) { LOG_WARN("PillarBuildParams: tf_apply_eff differs"); return true; }
  if (prm_noise != other.prm_noise) { LOG_WARN("PillarBuildParams: prm_noise differs"); return true; }
  #endif
  return false;
}

//##################################################################################
// @struct DetectorPanelArray::VoxelBuildParams
// @brief Structure collecting parameters needed when calculating path length, signal count, etc. of Grid3dVoxel in DetectorPanelArray
//##################################################################################

// @brief Inequality operator
bool DetectorPanelArray::VoxelBuildParams::operator!=(const VoxelBuildParams& other) const {
  #ifdef NODEBUG
    if (tf_load_arrdet_g3vox != other.tf_load_arrdet_g3vox) return true;
    if (tf_save_arrdet_g3vox != other.tf_save_arrdet_g3vox) return true;
    if (tf_run_1st_grouping != other.tf_run_1st_grouping) return true;
    if (tf_run_auto_grouping != other.tf_run_auto_grouping) return true;
    if (path_arrdet_bin != other.path_arrdet_bin) return true;
    if (path_vec_spmat_PL_bin != other.path_vec_spmat_PL_bin) return true;
    if (arrdet_built != other.arrdet_built) return true;
    if (prm_bingroup != other.prm_bingroup) return true;
    if (g3vox != other.g3vox) return true;
    if (ft != other.ft) return true;
    if (prm_pathcalc != other.prm_pathcalc) return true;
    if (prm_g3vox != other.prm_g3vox) return true;
    if (g2pil_shell_upper != other.g2pil_shell_upper) return true;
    if (g2pil_shell_lower != other.g2pil_shell_lower) return true;
    if (g2pil_shell_lateral != other.g2pil_shell_lateral) return true;
    if (has_shell_upper != other.has_shell_upper) return true;
    if (has_shell_lower != other.has_shell_lower) return true;
    if (has_shell_lateral != other.has_shell_lateral) return true;
    if (uqiv_start != other.uqiv_start) return true;
    if (tf_apply_eff != other.tf_apply_eff) return true;
    if (tf_apply_eff_cnt != other.tf_apply_eff_cnt) return true;
    if (prm_noise != other.prm_noise) return true;
  #else
    if (tf_load_arrdet_g3vox != other.tf_load_arrdet_g3vox) { LOG_WARN("VoxelBuildParams: tf_load_arrdet_g3vox differs"); return true; }
    if (tf_save_arrdet_g3vox != other.tf_save_arrdet_g3vox) { LOG_WARN("VoxelBuildParams: tf_save_arrdet_g3vox differs"); return true; }
    if (tf_run_1st_grouping != other.tf_run_1st_grouping) { LOG_WARN("VoxelBuildParams: tf_run_1st_grouping differs"); return true; }
    if (tf_run_auto_grouping != other.tf_run_auto_grouping) { LOG_WARN("VoxelBuildParams: tf_run_auto_grouping differs"); return true; }
    if (path_arrdet_bin != other.path_arrdet_bin) { LOG_WARN("VoxelBuildParams: path_arrdet_bin differs"); return true; }
    if (path_vec_spmat_PL_bin != other.path_vec_spmat_PL_bin) { LOG_WARN("VoxelBuildParams: path_vec_spmat_PL_bin differs"); return true; }
    if (arrdet_built != other.arrdet_built) { LOG_WARN("VoxelBuildParams: arrdet_built differs"); return true; }
    if (prm_bingroup != other.prm_bingroup) { LOG_WARN("VoxelBuildParams: prm_bingroup differs"); return true; }
    if (g3vox != other.g3vox) { LOG_WARN("VoxelBuildParams: g3vox differs"); return true; }
    if (ft != other.ft) { LOG_WARN("VoxelBuildParams: ft differs"); return true; }
    if (prm_pathcalc != other.prm_pathcalc) { LOG_WARN("VoxelBuildParams: prm_pathcalc differs"); return true; }
    if (prm_g3vox != other.prm_g3vox) { LOG_WARN("VoxelBuildParams: prm_g3vox differs"); return true; }
    if (g2pil_shell_upper != other.g2pil_shell_upper) { LOG_WARN("VoxelBuildParams: g2pil_shell_upper differs"); return true; }
    if (g2pil_shell_lower != other.g2pil_shell_lower) { LOG_WARN("VoxelBuildParams: g2pil_shell_lower differs"); return true; }
    if (g2pil_shell_lateral != other.g2pil_shell_lateral) { LOG_WARN("VoxelBuildParams: g2pil_shell_lateral differs"); return true; }
    if (has_shell_upper != other.has_shell_upper) { LOG_WARN("VoxelBuildParams: has_shell_upper differs"); return true; }
    if (has_shell_lower != other.has_shell_lower) { LOG_WARN("VoxelBuildParams: has_shell_lower differs"); return true; }
    if (has_shell_lateral != other.has_shell_lateral) { LOG_WARN("VoxelBuildParams: has_shell_lateral differs"); return true; }
    if (uqiv_start != other.uqiv_start) { LOG_WARN("VoxelBuildParams: uqiv_start differs"); return true; }
    if (tf_apply_eff != other.tf_apply_eff) { LOG_WARN("VoxelBuildParams: tf_apply_eff differs"); return true; }
    if (tf_apply_eff_cnt != other.tf_apply_eff_cnt) { LOG_WARN("VoxelBuildParams: tf_apply_eff_cnt differs"); return true; }
    if (prm_noise != other.prm_noise) { LOG_WARN("VoxelBuildParams: prm_noise differs"); return true; }
  #endif
  return false;
}
