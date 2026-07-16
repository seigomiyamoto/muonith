// cls_Grid2dPillarParameters.cpp
#include "cls_Grid2dPillar.hpp"
#include "cls_Grid2dPillarParameters.hpp"
#include "ns_myapp.hpp"
#include "ns_mymacro.hpp"
#include "ns_mylogger.hpp"
#include "ns_param_util.hpp"
#include "spdlog_pch.hpp"

//####################################
// parameters for Grid2dPillar
//####################################

void Grid2dPillar::Parameters::assign_parameters(
  const nlohmann::json& js, const std::string &section_name)
{
  LOG_INFO("...section name={}", section_name);

  // Read parameters from JSON
  param_util::read_json_value(js, section_name, TOSTRING(name), name);
  param_util::read_json_value(js, section_name, TOSTRING(initial_uniform_density), initial_uniform_density);
  param_util::read_json_path(js, section_name, TOSTRING(path_dem), path_dem);
  param_util::read_json_value(js, section_name, TOSTRING(tolerance_ratio), tolerance_ratio);
  param_util::read_json_value(js, section_name, TOSTRING(zmin), zmin);
  param_util::read_json_value(js, section_name, TOSTRING(tf_shift_x), tf_shift_x);
  param_util::read_json_value(js, section_name, TOSTRING(tf_shift_y), tf_shift_y);

  // validate tolerance_ratio > 0
  if (tolerance_ratio <= 0.0) {
    THROW_ERROR("Grid2dPillar::Parameters::assign_parameters: tolerance_ratio must be positive, got {}", tolerance_ratio);
  }

  // validate initial_uniform_density > 0
  if (initial_uniform_density <= 0.0) {
    THROW_ERROR("Grid2dPillar::Parameters::assign_parameters: initial_uniform_density must be positive, got {}", initial_uniform_density);
  }

  // validate path_dem is not empty
  if (path_dem.string().empty()) {
    THROW_ERROR("Grid2dPillar::Parameters::assign_parameters: path_dem must not be empty");
  }

  LOG_INFO("reading additional density structure parameters...");

  // read vertical cylinder parameters
  if(js.contains(section_name) && js.at(section_name).contains(section_name_cylinder)){
    const auto& cyl_array = js.at(section_name).at(section_name_cylinder);
    const auto before_count = vec_vertical_cylinder_parameters.size();
    for (const auto &cyl_js : cyl_array) {
      Grid2dPillar::VerticalCylinderParameters prm;
      prm.assign_parameters(cyl_js, section_name);
      vec_vertical_cylinder_parameters.push_back(prm);
    }
    if (!cyl_array.empty() && vec_vertical_cylinder_parameters.size() == before_count) {
      LOG_WARN("'{}' lists {} entries, but none were parsed.",
        section_name_cylinder, cyl_array.size());
    }
  }

  // read vertical dike parameters
  if(js.contains(section_name) && js.at(section_name).contains(section_name_dike)){
    const auto& dike_array = js.at(section_name).at(section_name_dike);
    const auto before_count = vec_vertical_dike_parameters.size();
    for (const auto &dike_js : dike_array) {
      Grid2dPillar::VerticalDikeParameters prm;
      prm.assign_parameters(dike_js, section_name);
      vec_vertical_dike_parameters.push_back(prm);
    }
    if (!dike_array.empty() && vec_vertical_dike_parameters.size() == before_count) {
      LOG_WARN("'{}' lists {} entries, but none were parsed.",
        section_name_dike, dike_array.size());
    }
  }

  // read vertical checkerboard parameters
  if(js.contains(section_name) && js.at(section_name).contains(section_name_checkerboard)){
    const auto& checker_array = js.at(section_name).at(section_name_checkerboard);
    const auto before_count = vec_vertical_checkerboard_parameters.size();
    for (const auto &checker_js : checker_array) {
      Grid2dPillar::VerticalCheckerBoardParameters prm;
      prm.assign_parameters(checker_js, section_name);
      vec_vertical_checkerboard_parameters.push_back(prm);
    }
    if (!checker_array.empty() && vec_vertical_checkerboard_parameters.size() == before_count) {
      LOG_WARN("'{}' lists {} entries, but none were parsed.",
        section_name_checkerboard, checker_array.size());
    }
  }

  LOG_INFO("reading additional density structure parameters...done");
}

//####################################
// Grid2dPillar::VerticalCylinderParameters
//####################################

void Grid2dPillar::VerticalCylinderParameters::assign_parameters(
  const nlohmann::json& js, const std::string &section_name)
{
  LOG_INFO("...section_name={}", section_name);

  std::string entry_name = "unnamed_vertical_cylinder";
  const std::string class_name = "VerticalCylinderParameters";

  // Read name first to use it in subsequent warnings
  if(js.contains(TOSTRING(name))) {
    entry_name = js.at(TOSTRING(name)).get<std::string>();
    name = entry_name;
  } else {
    LOG_WARN("{}: key '{}' missing in section '{}' (name='{}').",
      class_name, TOSTRING(name), section_name, entry_name);
    name = entry_name;
  }

  // Read parameters with warnings
  param_util::read_json_value_warn(js, TOSTRING(tf_exec), tf_exec, class_name, entry_name);
  param_util::read_json_value_warn(js, TOSTRING(delta_density), delta_density, class_name, entry_name);

  // Read vector components
  double x_pos_cnt = v2_pos_cnt.x();
  double y_pos_cnt = v2_pos_cnt.y();
  double x_length = v2_length.x();
  double y_length = v2_length.y();

  param_util::read_json_value_warn(js, TOSTRING(x_pos_cnt), x_pos_cnt, class_name, entry_name);
  param_util::read_json_value_warn(js, TOSTRING(y_pos_cnt), y_pos_cnt, class_name, entry_name);
  param_util::read_json_value_warn(js, TOSTRING(x_length), x_length, class_name, entry_name);
  param_util::read_json_value_warn(js, TOSTRING(y_length), y_length, class_name, entry_name);

  v2_pos_cnt.x() = x_pos_cnt;
  v2_pos_cnt.y() = y_pos_cnt;
  v2_length.x() = x_length;
  v2_length.y() = y_length;

  // Read rotation angle
  double rotation_angle_deg = 0.0;
  param_util::read_json_value_warn(js, TOSTRING(rotation_angle_deg), rotation_angle_deg, class_name, entry_name);

  // Set rotation angle: negate to convert from input convention (clockwise-positive)
  // to internal convention (counterclockwise-positive, right-hand rule around z-up)
  angle_rot.setDegree(-rotation_angle_deg);
}

//####################################
// Grid2dPillar::VerticalDikeParameters
//####################################

void Grid2dPillar::VerticalDikeParameters::assign_parameters(
  const nlohmann::json& js, const std::string &section_name)
{
  LOG_INFO("... section_name={}", section_name);

  std::string entry_name = "unnamed_vertical_dike";
  const std::string class_name = "VerticalDikeParameters";

  // Read name first to use it in subsequent warnings
  if(js.contains(TOSTRING(name))) {
    entry_name = js.at(TOSTRING(name)).get<std::string>();
    name = entry_name;
  } else {
    LOG_WARN("{}: key '{}' missing in section '{}' (name='{}').",
      class_name, TOSTRING(name), section_name, entry_name);
    name = entry_name;
  }

  // Read parameters with warnings
  param_util::read_json_value_warn(js, TOSTRING(tf_exec), tf_exec, class_name, entry_name);
  param_util::read_json_value_warn(js, TOSTRING(delta_density), delta_density, class_name, entry_name);

  // Read position and length parameters
  double x_pos_cnt = 0.0;
  double y_pos_cnt = 0.0;
  double x_length = 0.0;
  double y_length = 0.0;

  param_util::read_json_value_warn(js, TOSTRING(x_pos_cnt), x_pos_cnt, class_name, entry_name);
  param_util::read_json_value_warn(js, TOSTRING(y_pos_cnt), y_pos_cnt, class_name, entry_name);
  param_util::read_json_value_warn(js, TOSTRING(x_length), x_length, class_name, entry_name);
  param_util::read_json_value_warn(js, TOSTRING(y_length), y_length, class_name, entry_name);

  // calc xmin, xmax
  const double xmin = x_pos_cnt - 0.5 * x_length;
  const double xmax = x_pos_cnt + 0.5 * x_length;
  aabb2d.set_xmin(xmin);
  aabb2d.set_xmax(xmax);

  // calc ymin, ymax
  const double ymin = y_pos_cnt - 0.5 * y_length;
  const double ymax = y_pos_cnt + 0.5 * y_length;
  aabb2d.set_ymin(ymin);
  aabb2d.set_ymax(ymax);

  // Read rotation angle
  double rotation_angle_deg = 0.0;
  param_util::read_json_value_warn(js, TOSTRING(rotation_angle_deg), rotation_angle_deg, class_name, entry_name);

  // Set rotation angle: negate to convert from input convention (clockwise-positive)
  // to internal convention (counterclockwise-positive, right-hand rule around z-up)
  angle_rot.setDegree(-rotation_angle_deg);
}

//####################################
// Grid2dPillar::VerticalCheckerBoardParameters
//####################################

void Grid2dPillar::VerticalCheckerBoardParameters::assign_parameters(const nlohmann::json& js, const std::string &section_name)
{
  LOG_INFO("...section_name={}",
    section_name);

  std::string entry_name = "unnamed_vertical_checkerboard";
  const std::string class_name = "VerticalCheckerBoardParameters";

  // Read name first to use it in subsequent warnings
  if(js.contains(TOSTRING(name))) {
    entry_name = js.at(TOSTRING(name)).get<std::string>();
    name = entry_name;
  } else {
    LOG_WARN("{}: key '{}' missing in section '{}' (name='{}').",
      class_name, TOSTRING(name), section_name, entry_name);
    name = entry_name;
  }

  // Read parameters with warnings
  param_util::read_json_value_warn(js, TOSTRING(tf_exec), tf_exec, class_name, entry_name);
  param_util::read_json_value_warn(js, TOSTRING(delta_density_offset), delta_density_offset, class_name, entry_name);
  param_util::read_json_value_warn(js, TOSTRING(delta_density), delta_density, class_name, entry_name);

  // Read AABB bounds
  double xmin = 0.0;
  double xmax = 0.0;
  double ymin = 0.0;
  double ymax = 0.0;

  param_util::read_json_value_warn(js, TOSTRING(xmin), xmin, class_name, entry_name);
  param_util::read_json_value_warn(js, TOSTRING(xmax), xmax, class_name, entry_name);
  param_util::read_json_value_warn(js, TOSTRING(ymin), ymin, class_name, entry_name);
  param_util::read_json_value_warn(js, TOSTRING(ymax), ymax, class_name, entry_name);

  aabb2d.set_xmin(xmin);
  aabb2d.set_xmax(xmax);
  aabb2d.set_ymin(ymin);
  aabb2d.set_ymax(ymax);

  // Read center position and length
  double xcnt = 0.0;
  double ycnt = 0.0;
  double xlen = 0.0;
  double ylen = 0.0;

  param_util::read_json_value_warn(js, TOSTRING(xcnt), xcnt, class_name, entry_name);
  param_util::read_json_value_warn(js, TOSTRING(ycnt), ycnt, class_name, entry_name);
  param_util::read_json_value_warn(js, TOSTRING(xlen), xlen, class_name, entry_name);
  param_util::read_json_value_warn(js, TOSTRING(ylen), ylen, class_name, entry_name);

  v2_pos_cnt.x() = xcnt;
  v2_pos_cnt.y() = ycnt;
  v2_length.x() = xlen;
  v2_length.y() = ylen;
}
