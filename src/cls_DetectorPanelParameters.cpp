// src/cls_DetectorPanelParameters.cpp

#include "cls_DetectorPanelParameters.hpp"
#include "cls_DetectorPanel.hpp"
#include "ns_angle_util.hpp"
#include "ns_myapp.hpp"
#include "ns_mylogger.hpp"
#include "ns_mymacro.hpp"
#include "ns_param_util.hpp"
#include "spdlog_pch.hpp"

// ##############################################
// ##############################################
//  class DetectorPanel::Parameters
// ##############################################
// ##############################################

// assign_parameters
void DetectorPanel::Parameters::assign_parameters(
    const nlohmann::json &js,
    const std::string &section_name
) {
  LOG_INFO("...section_name={}", section_name);

  // Read basic parameters using param_util
  param_util::read_json_value(js, section_name, TOSTRING(name), name);
  param_util::read_json_value(
      js,
      section_name,
      TOSTRING(tf_read_bin_list),
      tf_read_bin_list
  );

  if (tf_read_bin_list == true) {
    param_util::read_json_value(
        js,
        section_name,
        TOSTRING(filepath_bin_list),
        filepath_bin_list
    );
  }

  // angle_unit
  std::string angle_unit_str;
  param_util::read_json_value(js, section_name, TOSTRING(angle_unit), angle_unit_str);
  try {
    angle_unit = angle_util::angle_unit_from_string(angle_unit_str);
  } catch (const std::exception &e) {
    THROW_ERROR(
        "DetectorPanel::Parameters::assign_parameters: invalid angle_unit string: {} "
        "({})",
        angle_unit_str,
        e.what()
    );
  }

  // Read grid parameters
  param_util::read_json_value(js, section_name, TOSTRING(nbinx), nbinx);
  param_util::read_json_value(js, section_name, TOSTRING(txmin), txmin);
  param_util::read_json_value(js, section_name, TOSTRING(txmax), txmax);
  param_util::read_json_value(js, section_name, TOSTRING(nbiny), nbiny);
  param_util::read_json_value(js, section_name, TOSTRING(tymin), tymin);

  if (angle_unit == DetectorElement::AngleUnit::Degree) {
    assert(-90 <= tymin && tymin <= 90);
  } else if (angle_unit == DetectorElement::AngleUnit::Radian) {
    assert(-M_PI / 2 <= tymin && tymin <= M_PI / 2);
  }

  param_util::read_json_value(js, section_name, TOSTRING(tymax), tymax);

  // Read detector dimensions
  double length_hori = 0.0, length_vert = 0.0, length_dept = 0.0;
  param_util::read_json_value(js, section_name, TOSTRING(length_hori), length_hori);
  param_util::read_json_value(js, section_name, TOSTRING(length_vert), length_vert);
  param_util::read_json_value(js, section_name, TOSTRING(length_dept), length_dept);
  v3_det_length = Eigen::Vector3d(length_hori, length_vert, length_dept);

  param_util::read_json_value(js, section_name, TOSTRING(n_unit), n_unit);

  // rotation type
  std::string rotation_type_str;
  param_util::read_json_value(
      js,
      section_name,
      TOSTRING(rotation_type),
      rotation_type_str
  );
  if (rotation_type_str == "LOCAL" || rotation_type_str == "local" ||
      rotation_type_str == "Local") {
    rotation_type = angle_util::Rotation3dType::LOCAL;
  } else if (rotation_type_str == "GLOBAL" || rotation_type_str == "global" ||
             rotation_type_str == "Global") {
    rotation_type = angle_util::Rotation3dType::GLOBAL;
  } else {
    THROW_ERROR(
        "DetectorPanel::Parameters::assign_parameters: rotation_type must be LOCAL or "
        "GLOBAL, got: {}",
        rotation_type_str
    );
  }
  LOG_DEBUG("rotation_type_str={}", rotation_type_str);

  // Read rotation angles
  double yaw_deg = 0.0, roll_deg = 0.0, pitch_deg = 0.0;
  param_util::read_json_value(js, section_name, TOSTRING(yaw_deg), yaw_deg);
  param_util::read_json_value(js, section_name, TOSTRING(pitch_deg), pitch_deg);
  param_util::read_json_value(js, section_name, TOSTRING(roll_deg), roll_deg);
  yaw.setDegree(yaw_deg);
  pitch.setDegree(pitch_deg);
  roll.setDegree(roll_deg);

  // Read position
  double x = 0.0, y_val = 0.0, z = 0.0;
  param_util::read_json_value(js, section_name, TOSTRING(x), x);
  param_util::read_json_value(js, section_name, TOSTRING(y), y_val);
  param_util::read_json_value(js, section_name, TOSTRING(z), z);
  v3_position = Eigen::Vector3d(x, y_val, z);

  // Read remaining parameters
  param_util::read_json_value(js, section_name, TOSTRING(days), days);
  param_util::read_json_value(
      js,
      section_name,
      TOSTRING(n_reserve_vec_tf_in_PL),
      n_reserve_vec_tf_in_PL
  );
  param_util::read_json_value(js, section_name, TOSTRING(path_eff_table), path_eff_table);
  param_util::read_json_value(js, section_name, TOSTRING(path_eff_model), path_eff_model);

  // path_eff_model takes precedence over the legacy path_eff_table (see
  // DetectorPanelArray::build_vec_panel). The model coefficients are loaded
  // here so the panel-build stage only evaluates them.
  if (path_eff_model != "none") {
    eff_model.load(path_eff_model);
  }

  assert(nbinx > 0);
  assert(nbiny > 0);
  assert(!(txmin == 0.0 && txmax == 0.0));
  assert(!(tymin == 0.0 && tymax == 0.0));
  assert(txmin < txmax);
  assert(tymin < tymax);
  assert(length_hori > 0.0);
  assert(length_vert > 0.0);
  assert(length_dept >= 0.0);
  assert(n_unit > 0.0);
  assert(n_reserve_vec_tf_in_PL >= 0);
  LOG_INFO("assign_parameters : cheking all parameters assigned correctly...OK");
}

// ##################################################################################
// ##################################################################################
//  DetectorPanelParameterLists
//  Parameter lists for std::vector<DetectorPanel>
// ##################################################################################
// ##################################################################################

// setting parameters from RunCard
void DetectorPanel::ParameterLists::assign_parameters(
    const nlohmann::json &js,
    const std::string &section_name
) {
  LOG_INFO("...section_name={}", section_name);

  // Read basic parameters
  param_util::read_json_value(js, section_name, TOSTRING(name), name);

  // Retrieve vec_parameter_file_path directly from "det_files" array
  vec_parameter_file_path.clear();
  if (js.contains(section_name) && js.at(section_name).contains(TOSTRING(det_files))) {
    const auto &det_files = js.at(section_name).at(TOSTRING(det_files));
    if (det_files.is_array()) {
      for (const auto &file : det_files)
        vec_parameter_file_path.push_back(file.get<std::string>());
    }
  }

  // Check that all files in vec_parameter_file_path exist
  for (const auto &path : vec_parameter_file_path)
    myapp::filecheck(path);

  // tf_apply_eff
  if (js.contains(section_name) && js.at(section_name).contains("tf_apply_eff"))
    tf_apply_eff = js.at(section_name).at("tf_apply_eff").get<bool>();
  else
    tf_apply_eff = false;

  // Reject the pre-rename keys instead of silently ignoring them
  if (js.contains(section_name) && js.at(section_name).contains("tf_out_txty"))
    THROW_ERROR(
        "DetectorPanel::ParameterLists::assign_parameters: {} contains the obsolete key "
        "tf_out_txty. Rename it to tf_out_txty_ascii.",
        section_name);
  if (js.contains(section_name) && js.at(section_name).contains("tf_out_g2bg"))
    THROW_ERROR(
        "DetectorPanel::ParameterLists::assign_parameters: {} contains the obsolete key "
        "tf_out_g2bg. Rename it to tf_out_g2bg_ascii.",
        section_name);

  // tf_out_txty_ascii / tf_out_g2bg_ascii: optional, default true = keep writing the ASCII dumps
  if (js.contains(section_name) && js.at(section_name).contains("tf_out_txty_ascii"))
    tf_out_txty_ascii = js.at(section_name).at("tf_out_txty_ascii").get<bool>();
  else
    tf_out_txty_ascii = true;

  if (js.contains(section_name) && js.at(section_name).contains("tf_out_g2bg_ascii"))
    tf_out_g2bg_ascii = js.at(section_name).at("tf_out_g2bg_ascii").get<bool>();
  else
    tf_out_g2bg_ascii = true;

  // tf_save_arrdet_prior: optional, default false = do not write the ~160 MB prior binary
  if (js.contains(section_name) && js.at(section_name).contains("tf_save_arrdet_prior"))
    tf_save_arrdet_prior = js.at(section_name).at("tf_save_arrdet_prior").get<bool>();
  else
    tf_save_arrdet_prior = false;
}

void DetectorPanel::ParameterLists::set_angle_bin_override(
    int nbinx, double txmin, double txmax,
    int nbiny, double tymin, double tymax)
{
  tf_override_angle_bin = true;
  nbinx_override = nbinx;
  txmin_override = txmin;
  txmax_override = txmax;
  nbiny_override = nbiny;
  tymin_override = tymin;
  tymax_override = tymax;
  LOG_INFO("angle_bin_override set: nbinx={}, tx=[{},{}], nbiny={}, ty=[{},{}]",
           nbinx, txmin, txmax, nbiny, tymin, tymax);
}
