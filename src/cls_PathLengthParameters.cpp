/// @file cls_PathLengthParameters.cpp
/// @brief Implementation of pathcalc::Parameters class
/// @details
/// Provides JSON parameter loading and comparison operators for path length
/// calculation configuration.
///
/// @see cls_PathLengthParameters.hpp

#include "cls_PathLengthParameters.hpp"

#include "ns_mymacro.hpp"
#include "ns_mylogger.hpp"
#include "ns_param_util.hpp"
#include "ns_param_constants.hpp"
#include <cassert>

bool pathcalc::Parameters::operator!=(const Parameters& other) const
{
#ifdef NODEBUG
  if (PL_min != other.PL_min) return true;
  if (PL_max != other.PL_max) return true;
  if (PL_pit != other.PL_pit) return true;
  if (BL_max != other.BL_max) return true;
  if (tf_load_arrdet_g2pil != other.tf_load_arrdet_g2pil) return true;
  if (tf_save_arrdet_g2pil != other.tf_save_arrdet_g2pil) return true;
  if (path_arrdet_g2pil_bin != other.path_arrdet_g2pil_bin) return true;
  if (tf_load_arrdet_g3vox != other.tf_load_arrdet_g3vox) return true;
  if (tf_save_arrdet_g3vox != other.tf_save_arrdet_g3vox) return true;
  if (tf_add_shell != other.tf_add_shell) return true;
  if (path_arrdet_g3vox_bin != other.path_arrdet_g3vox_bin) return true;
  if (path_vec_spmat_PL_bin != other.path_vec_spmat_PL_bin) return true;
  if (tf_add_PLDL != other.tf_add_PLDL) return true;
  if (tf_incr_nhit_det != other.tf_incr_nhit_det) return true;
  if (tf_incr_nhit_ele != other.tf_incr_nhit_ele) return true;
  if (reference_matPL_sparse != other.reference_matPL_sparse) return true;
  if (epsilon_matPL_sparse != other.epsilon_matPL_sparse) return true;
  if (tf_load_bin_obs_mat_dNdD != other.tf_load_bin_obs_mat_dNdD) return true;
  if (tf_save_bin_obs_mat_dNdD != other.tf_save_bin_obs_mat_dNdD) return true;
  if (path_bin_obs_mat_dNdD != other.path_bin_obs_mat_dNdD) return true;
#else
  if (PL_min != other.PL_min) { LOG_WARN("PL_min differs"); return true; }
  if (PL_max != other.PL_max) { LOG_WARN("PL_max differs"); return true; }
  if (PL_pit != other.PL_pit) { LOG_WARN("PL_pit differs"); return true; }
  if (BL_max != other.BL_max) { LOG_WARN("BL_max differs"); return true; }
  if (tf_load_arrdet_g2pil != other.tf_load_arrdet_g2pil) { LOG_WARN("tf_load_arrdet_g2pil differs"); return true; }
  if (tf_save_arrdet_g2pil != other.tf_save_arrdet_g2pil) { LOG_WARN("tf_save_arrdet_g2pil differs"); return true; }
  if (path_arrdet_g2pil_bin != other.path_arrdet_g2pil_bin) { LOG_WARN("path_arrdet_g2pil_bin differs"); return true; }
  if (tf_load_arrdet_g3vox != other.tf_load_arrdet_g3vox) { LOG_WARN("tf_load_arrdet_g3vox differs"); return true; }
  if (tf_save_arrdet_g3vox != other.tf_save_arrdet_g3vox) { LOG_WARN("tf_save_arrdet_g3vox differs"); return true; }
  if (tf_add_shell != other.tf_add_shell) { LOG_WARN("tf_add_shell differs"); return true; }
  if (path_arrdet_g3vox_bin != other.path_arrdet_g3vox_bin) { LOG_WARN("path_arrdet_g3vox_bin differs"); return true; }
  if (path_vec_spmat_PL_bin != other.path_vec_spmat_PL_bin) { LOG_WARN("path_vec_spmat_PL_bin differs"); return true; }
  if (tf_add_PLDL != other.tf_add_PLDL) { LOG_WARN("tf_add_PLDL differs"); return true; }
  if (tf_incr_nhit_det != other.tf_incr_nhit_det) { LOG_WARN("tf_incr_nhit_det differs"); return true; }
  if (tf_incr_nhit_ele != other.tf_incr_nhit_ele) { LOG_WARN("tf_incr_nhit_ele differs"); return true; }
  if (reference_matPL_sparse != other.reference_matPL_sparse) { LOG_WARN("reference_matPL_sparse differs"); return true; }
  if (epsilon_matPL_sparse != other.epsilon_matPL_sparse) { LOG_WARN("epsilon_matPL_sparse differs"); return true; }
  if (tf_load_bin_obs_mat_dNdD != other.tf_load_bin_obs_mat_dNdD) { LOG_WARN("tf_load_bin_obs_mat_dNdD differs"); return true; }
  if (tf_save_bin_obs_mat_dNdD != other.tf_save_bin_obs_mat_dNdD) { LOG_WARN("tf_save_bin_obs_mat_dNdD differs"); return true; }
  if (path_bin_obs_mat_dNdD != other.path_bin_obs_mat_dNdD) { LOG_WARN("path_bin_obs_mat_dNdD differs"); return true; }
#endif
  return false;
}

void pathcalc::Parameters::assign_parameters(
  const nlohmann::json& js, const std::string &section_name)
{
  LOG_INFO("section_name={}", section_name);

  // Instance name
  param_util::read_json_value_required(js, section_name, TOSTRING(name), name);

  // Path length range parameters [meters]
  param_util::read_json_value_required(js, section_name, TOSTRING(PL_min), PL_min);
  param_util::read_json_value_required(js, section_name, TOSTRING(PL_max), PL_max);
  param_util::read_json_value_required(js, section_name, TOSTRING(PL_pit), PL_pit);

  assert(PL_min < PL_max);
  assert(PL_pit > 0.0);

  // Beam length warning for ray tracing (optional, default from param_constants)
  const auto& section = js.at(section_name);
  BL_max = section.value("BL_max", param_constants::BL_max_default());
  LOG_INFO("BL_max={} (<=0 disables warning)", BL_max);

  // g2pil I/O parameters
  param_util::read_json_value_required(js, section_name, TOSTRING(tf_load_arrdet_g2pil), tf_load_arrdet_g2pil);
  param_util::read_json_value_required(js, section_name, TOSTRING(tf_save_arrdet_g2pil), tf_save_arrdet_g2pil);
  param_util::read_json_path_required(js, section_name, TOSTRING(path_arrdet_g2pil_bin), path_arrdet_g2pil_bin);

  // g3vox I/O parameters
  param_util::read_json_value_required(js, section_name, TOSTRING(tf_load_arrdet_g3vox), tf_load_arrdet_g3vox);
  param_util::read_json_value_required(js, section_name, TOSTRING(tf_save_arrdet_g3vox), tf_save_arrdet_g3vox);
  param_util::read_json_value_required(js, section_name, TOSTRING(tf_add_shell), tf_add_shell);
  param_util::read_json_path_required(js, section_name, TOSTRING(path_arrdet_g3vox_bin), path_arrdet_g3vox_bin);
  param_util::read_json_path_required(js, section_name, TOSTRING(path_vec_spmat_PL_bin), path_vec_spmat_PL_bin);

  // g3vox calculation control parameters
  param_util::read_json_value_required(js, section_name, TOSTRING(tf_add_PLDL), tf_add_PLDL);
  param_util::read_json_value_required(js, section_name, TOSTRING(tf_incr_nhit_det), tf_incr_nhit_det);
  param_util::read_json_value_required(js, section_name, TOSTRING(tf_incr_nhit_ele), tf_incr_nhit_ele);
  param_util::read_json_value_required(js, section_name, TOSTRING(reference_matPL_sparse), reference_matPL_sparse);
  param_util::read_json_value_required(js, section_name, TOSTRING(epsilon_matPL_sparse), epsilon_matPL_sparse);

  // Observation matrix I/O parameters
  param_util::read_json_value_required(js, section_name, TOSTRING(tf_load_bin_obs_mat_dNdD), tf_load_bin_obs_mat_dNdD);
  param_util::read_json_value_required(js, section_name, TOSTRING(tf_save_bin_obs_mat_dNdD), tf_save_bin_obs_mat_dNdD);
  param_util::read_json_path_required(js, section_name, TOSTRING(path_bin_obs_mat_dNdD), path_bin_obs_mat_dNdD);
}
