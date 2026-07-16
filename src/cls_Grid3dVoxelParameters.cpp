/// @file cls_Grid3dVoxelParameters.cpp
/// @brief Implementation of Grid3dVoxel parameter classes

#include "cls_Grid3dVoxel.hpp"
#include "ns_myapp.hpp"
#include "ns_mymacro.hpp"
#include "ns_mylogger.hpp"
#include "cls_Grid3dVoxelParameters.hpp"
#include "ns_param_util.hpp"
#include <nlohmann/json.hpp>
#include "spdlog_pch.hpp"

//#######################################
// parameters for Grid3dVoxel::Parameters
//#######################################
bool Grid3dVoxel::Parameters::operator!=(const Parameters& other) const
{
  #ifdef NODEBUG
  if( xmin != other.xmin ) return true;
  if( xmax != other.xmax ) return true;
  if( x_pitch != other.x_pitch ) return true;
  if( ymin != other.ymin ) return true;
  if( ymax != other.ymax ) return true;
  if( y_pitch != other.y_pitch ) return true;
  if( zmin != other.zmin ) return true;
  if( zmax != other.zmax ) return true;
  if( z_pitch != other.z_pitch ) return true;
  if( n_hit_det_min != other.n_hit_det_min ) return true;
  if( n_hit_det_max != other.n_hit_det_max ) return true;
  if( n_hit_ele_min != other.n_hit_ele_min ) return true;
  if( n_hit_ele_max != other.n_hit_ele_max ) return true;
  if( use_grid2d_of_g2pil != other.use_grid2d_of_g2pil ) return true;
  if( tf_build_g3vox != other.tf_build_g3vox ) return true;
  if( vec_checkerboard_prm != other.vec_checkerboard_prm ) return true;
  if( vec_ellipsoid_prm != other.vec_ellipsoid_prm ) return true;
  if( vec_cylinder_prm != other.vec_cylinder_prm ) return true;
  if( vec_cuboid_prm != other.vec_cuboid_prm ) return true;
  if( prm_merge != other.prm_merge ) return true;
  if( prm_reconst_voxels != other.prm_reconst_voxels ) return true;
  #else
  if( xmin != other.xmin ) { LOG_WARN("xmin differs"); return true; }
  if( xmax != other.xmax ) { LOG_WARN("xmax differs"); return true; }
  if( x_pitch != other.x_pitch ) { LOG_WARN("x_pitch differs"); return true; }
  if( ymin != other.ymin ) { LOG_WARN("ymin differs"); return true; }
  if( ymax != other.ymax ) { LOG_WARN("ymax differs"); return true; }
  if( y_pitch != other.y_pitch ) { LOG_WARN("y_pitch differs"); return true; }
  if( zmin != other.zmin ) { LOG_WARN("zmin differs"); return true; }
  if( zmax != other.zmax ) { LOG_WARN("zmax differs"); return true; }
  if( z_pitch != other.z_pitch ) { LOG_WARN("z_pitch differs"); return true; }
  if( n_hit_det_min != other.n_hit_det_min ) { LOG_WARN("n_hit_det_min differs"); return true; }
  if( n_hit_det_max != other.n_hit_det_max ) { LOG_WARN("n_hit_det_max differs"); return true; }
  if( n_hit_ele_min != other.n_hit_ele_min ) { LOG_WARN("n_hit_ele_min differs"); return true; }
  if( n_hit_ele_max != other.n_hit_ele_max ) { LOG_WARN("n_hit_ele_max differs"); return true; }
  if( use_grid2d_of_g2pil != other.use_grid2d_of_g2pil ) { LOG_WARN("use_grid2d_of_g2pil differs"); return true; }
  if( tf_build_g3vox != other.tf_build_g3vox ) { LOG_WARN("tf_build_g3vox differs"); return true; }
  if( vec_checkerboard_prm != other.vec_checkerboard_prm ) { LOG_WARN("vec_checkerboard_prm differs"); return true; }
  if( vec_ellipsoid_prm != other.vec_ellipsoid_prm ) { LOG_WARN("vec_ellipsoid_prm differs"); return true; }
  if( vec_cylinder_prm != other.vec_cylinder_prm ) { LOG_WARN("vec_cylinder_prm differs"); return true; }
  if( vec_cuboid_prm != other.vec_cuboid_prm ) { LOG_WARN("vec_cuboid_prm differs"); return true; }
  if( prm_merge != other.prm_merge ) { LOG_WARN("prm_merge differs"); return true; }
  if( prm_reconst_voxels != other.prm_reconst_voxels ) { LOG_WARN("prm_reconst_voxels differs"); return true; }
  #endif
  return false;
}

void Grid3dVoxel::Parameters::assign_parameters(const nlohmann::json& js, const std::string &section_name)
{
  LOG_INFO("...section_name = {}", section_name);

  // Read name
  param_util::read_json_value(js, section_name, TOSTRING(name), name);

  // use_grid2d_of_g2pil (optional)
  bool prev_use_grid2d = use_grid2d_of_g2pil;
  param_util::read_json_value(js, section_name, TOSTRING(use_grid2d_of_g2pil), use_grid2d_of_g2pil);
  if(prev_use_grid2d == use_grid2d_of_g2pil && !js.at(section_name).contains(TOSTRING(use_grid2d_of_g2pil))){
    LOG_INFO("{} not found. use default value = {}",
      TOSTRING(use_grid2d_of_g2pil), use_grid2d_of_g2pil);
  }
  if( use_grid2d_of_g2pil ){
    LOG_INFO("use_grid2d_of_g2pil=true, skip reading xmin/xmax/x_pitch/ymin/ymax/y_pitch");
  }

  // tf_build_g3vox (optional)
  bool prev_tf_build = tf_build_g3vox;
  param_util::read_json_value(js, section_name, TOSTRING(tf_build_g3vox), tf_build_g3vox);
  if(prev_tf_build == tf_build_g3vox && !js.at(section_name).contains(TOSTRING(tf_build_g3vox))){
    LOG_INFO("{} not found. use default value = {}",
      TOSTRING(tf_build_g3vox), tf_build_g3vox);
  }
  if( tf_build_g3vox == false ){
    LOG_INFO("tf_build_g3vox=false, skip reading other Grid3dVoxel parameters.");
    return;
  }

  // Read grid parameters conditionally
  if( !use_grid2d_of_g2pil ){
    param_util::read_json_value(js, section_name, TOSTRING(xmin), xmin);
    param_util::read_json_value(js, section_name, TOSTRING(xmax), xmax);
    param_util::read_json_value(js, section_name, TOSTRING(x_pitch), x_pitch);
    if(x_pitch <= 0) THROW_ERROR("Grid3dVoxel::Parameters::assign_parameters: x_pitch must be positive. x_pitch={}", x_pitch);
    param_util::read_json_value(js, section_name, TOSTRING(ymin), ymin);
    param_util::read_json_value(js, section_name, TOSTRING(ymax), ymax);
    param_util::read_json_value(js, section_name, TOSTRING(y_pitch), y_pitch);
    if(y_pitch <= 0) THROW_ERROR("Grid3dVoxel::Parameters::assign_parameters: y_pitch must be positive. y_pitch={}", y_pitch);
  }

  // Read z parameters (always)
  param_util::read_json_value(js, section_name, TOSTRING(zmin), zmin);
  param_util::read_json_value(js, section_name, TOSTRING(zmax), zmax);
  param_util::read_json_value(js, section_name, TOSTRING(z_pitch), z_pitch);
  if(z_pitch <= 0) THROW_ERROR("Grid3dVoxel::Parameters::assign_parameters: z_pitch must be positive. z_pitch={}", z_pitch);

  // Read hit count parameters
  param_util::read_json_value(js, section_name, TOSTRING(n_hit_det_min), n_hit_det_min);
  if(n_hit_det_min < 0) THROW_ERROR("Grid3dVoxel::Parameters::assign_parameters: n_hit_det_min must be non-negative. n_hit_det_min={}", n_hit_det_min);
  param_util::read_json_value(js, section_name, TOSTRING(n_hit_det_max), n_hit_det_max);
  if(n_hit_det_max < n_hit_det_min) THROW_ERROR("Grid3dVoxel::Parameters::assign_parameters: n_hit_det_max must be >= n_hit_det_min. n_hit_det_min={}, n_hit_det_max={}", n_hit_det_min, n_hit_det_max);
  param_util::read_json_value(js, section_name, TOSTRING(n_hit_ele_min), n_hit_ele_min);
  if(n_hit_ele_min < 0) THROW_ERROR("Grid3dVoxel::Parameters::assign_parameters: n_hit_ele_min must be non-negative. n_hit_ele_min={}", n_hit_ele_min);
  param_util::read_json_value(js, section_name, TOSTRING(n_hit_ele_max), n_hit_ele_max);
  if(n_hit_ele_max < n_hit_ele_min) THROW_ERROR("Grid3dVoxel::Parameters::assign_parameters: n_hit_ele_max must be >= n_hit_ele_min. n_hit_ele_min={}, n_hit_ele_max={}", n_hit_ele_min, n_hit_ele_max);

  // Read tf_end_after_merged
  param_util::read_json_value(js, section_name, TOSTRING(tf_end_after_merged), tf_end_after_merged);

  // check all parameters assigned correctly
  if( !use_grid2d_of_g2pil ){
    if(xmin >= xmax) THROW_ERROR("Grid3dVoxel::Parameters::assign_parameters: xmin must be < xmax. xmin={}, xmax={}", xmin, xmax);
    if(ymin >= ymax) THROW_ERROR("Grid3dVoxel::Parameters::assign_parameters: ymin must be < ymax. ymin={}, ymax={}", ymin, ymax);
  }
  if(zmin >= zmax) THROW_ERROR("Grid3dVoxel::Parameters::assign_parameters: zmin must be < zmax. zmin={}, zmax={}", zmin, zmax);
  if(n_hit_det_min >= n_hit_det_max) THROW_ERROR("Grid3dVoxel::Parameters::assign_parameters: n_hit_det_min must be < n_hit_det_max. n_hit_det_min={}, n_hit_det_max={}", n_hit_det_min, n_hit_det_max);
  if(n_hit_ele_min >= n_hit_ele_max) THROW_ERROR("Grid3dVoxel::Parameters::assign_parameters: n_hit_ele_min must be < n_hit_ele_max. n_hit_ele_min={}, n_hit_ele_max={}", n_hit_ele_min, n_hit_ele_max);
  LOG_INFO("assign_parameters : checking all parameters assigned correctly...OK");

  // read checkerboard parameters
  if( js.at(section_name).contains(section_name_checkerboard) ){
    auto checker_array = js.at(section_name).at(section_name_checkerboard);
    for (const auto &chk : checker_array) {
      CheckerBoard3dParameters prm(chk);
      vec_checkerboard_prm.push_back(prm);
    }
  }
  LOG_INFO("read 3d_checkerboard parameters... num of checkerboard = {}"
    , vec_checkerboard_prm.size());

  // read ellipsoid parameters
  if( js.at(section_name).contains(section_name_ellipsoid) ){
    auto ellipsoid_array = js.at(section_name).at(section_name_ellipsoid);
    for (const auto &ell : ellipsoid_array) {
      EllipsoidParameters prm(ell);
      vec_ellipsoid_prm.push_back(prm);
    }
  }
  LOG_INFO("read ellipsoid parameters... num of ellipsoid = {}"
    , vec_ellipsoid_prm.size());

  // read cylinder parameters
  if( js.at(section_name).contains(section_name_cylinder) ){
    auto cylinder_array = js.at(section_name).at(section_name_cylinder);
    for (const auto &cyl : cylinder_array) {
      CylinderParameters prm(cyl);
      vec_cylinder_prm.push_back(prm);
    }
  }
  LOG_INFO("read cylinder parameters...num of cylinder = {}"
    , vec_cylinder_prm.size());

  // read cuboid parameters
  if( js.at(section_name).contains(section_name_cuboid) ){
    auto cuboid_array = js.at(section_name).at(section_name_cuboid);
    for (const auto &cub : cuboid_array) {
      CuboidParameters prm(cub);
      vec_cuboid_prm.push_back(prm);
    }
  }
  LOG_INFO("read cuboid parameters...num of cuboid = {}"
    , vec_cuboid_prm.size());

  // read merge parameters
  if( js.at(section_name).contains(section_name_merge) ){
    auto merge_obj = js.at(section_name).at(section_name_merge);
    prm_merge.assign_parameters(merge_obj);
  }
  LOG_INFO("read merge parameters...OK");

  // read reconst_voxels parameters
  if( js.at(section_name).contains(section_name_reconst_voxels) ){
    auto rv_obj = js.at(section_name).at(section_name_reconst_voxels);
    prm_reconst_voxels.assign_parameters(rv_obj);
    // Fallback: if center not specified, use merge_center
    if( !rv_obj.contains("x_aabb_cnt") ) prm_reconst_voxels.x_aabb_cnt = prm_merge.x_merge_center;
    if( !rv_obj.contains("y_aabb_cnt") ) prm_reconst_voxels.y_aabb_cnt = prm_merge.y_merge_center;
    if( !rv_obj.contains("x_cyl_cnt") )  prm_reconst_voxels.x_cyl_cnt  = prm_merge.x_merge_center;
    if( !rv_obj.contains("y_cyl_cnt") )  prm_reconst_voxels.y_cyl_cnt  = prm_merge.y_merge_center;
  }
  LOG_INFO("read reconst_voxels parameters...OK");
}

//####################################
// CheckerBoard3dParameters
//####################################
bool Grid3dVoxel::CheckerBoard3dParameters::operator!=(const CheckerBoard3dParameters& other) const
{
  #ifdef NODEBUG
  if( tf_exec != other.tf_exec ) return true;
  if( delta_density_offset != other.delta_density_offset ) return true;
  if( delta_density != other.delta_density ) return true;
  if( aabb3d != other.aabb3d ) return true;
  if( v3_pos_cnt != other.v3_pos_cnt ) return true;
  if( v3_len_interval_mult != other.v3_len_interval_mult ) return true;
  if( tf_snap_to_grid != other.tf_snap_to_grid ) return true;
  if( v3_len_cells != other.v3_len_cells ) return true;
  if( region_type != other.region_type ) return true;
  if( radius_x_meters != other.radius_x_meters ) return true;
  if( radius_y_meters != other.radius_y_meters ) return true;
  #else
  if( tf_exec != other.tf_exec ){ LOG_WARN("CheckerBoard3dParameters : tf_exec differs"); return true; }
  if( delta_density_offset != other.delta_density_offset ){ LOG_WARN("CheckerBoard3dParameters : delta_density_offset differs"); return true; }
  if( delta_density != other.delta_density ){ LOG_WARN("CheckerBoard3dParameters : delta_density differs"); return true; }
  if( aabb3d != other.aabb3d ){ LOG_WARN("CheckerBoard3dParameters : aabb3d differs"); return true; }
  if( v3_pos_cnt != other.v3_pos_cnt ){ LOG_WARN("CheckerBoard3dParameters : v3_pos_cnt differs"); return true; }
  if( v3_len_interval_mult != other.v3_len_interval_mult ){ LOG_WARN("CheckerBoard3dParameters : v3_len_interval_mult differs"); return true; }
  if( tf_snap_to_grid != other.tf_snap_to_grid ){ LOG_WARN("CheckerBoard3dParameters : tf_snap_to_grid differs"); return true; }
  if( v3_len_cells != other.v3_len_cells ){ LOG_WARN("CheckerBoard3dParameters : v3_len_cells differs"); return true; }
  if( region_type != other.region_type ){ LOG_WARN("CheckerBoard3dParameters : region_type differs"); return true; }
  if( radius_x_meters != other.radius_x_meters ){ LOG_WARN("CheckerBoard3dParameters : radius_x_meters differs"); return true; }
  if( radius_y_meters != other.radius_y_meters ){ LOG_WARN("CheckerBoard3dParameters : radius_y_meters differs"); return true; }
  #endif
  return false;
}

void Grid3dVoxel::CheckerBoard3dParameters::assign_parameters(const nlohmann::json& js)
{
  LOG_INFO("");

  // Read basic parameters
  param_util::read_json_direct(js, TOSTRING(tf_exec), tf_exec);
  param_util::read_json_direct(js, TOSTRING(name), name);
  param_util::read_json_direct(js, TOSTRING(delta_density_offset), delta_density_offset);
  param_util::read_json_direct(js, TOSTRING(delta_density), delta_density);

  // --- Reject deprecated direct AABB specification ---
  if (js.contains("xmin") || js.contains("xmax") ||
      js.contains("ymin") || js.contains("ymax") ||
      js.contains("zmin") || js.contains("zmax")) {
    THROW_ERROR("CheckerBoard3dParameters::assign_parameters: "
      "xmin/xmax/ymin/ymax/zmin/zmax are deprecated. "
      "Use xlen_cells/ylen_cells/zlen_cells instead "
      "(e.g. xlen_cells=6, ylen_cells=6, zlen_cells=4).");
  }

  // --- AABB specification: cells-based (required) ---
  if (!js.contains("xlen_cells") || !js.contains("ylen_cells") || !js.contains("zlen_cells")) {
    THROW_ERROR("CheckerBoard3dParameters::assign_parameters: "
      "xlen_cells, ylen_cells, and zlen_cells are required.");
  }
  {
    int xc = js.at("xlen_cells").get<int>();
    int yc = js.at("ylen_cells").get<int>();
    int zc = js.at("zlen_cells").get<int>();
    if (xc <= 0 || yc <= 0 || zc <= 0) {
      THROW_ERROR("CheckerBoard3dParameters::assign_parameters: "
        "xlen_cells/ylen_cells/zlen_cells must be positive. values=({},{},{})", xc, yc, zc);
    }
    v3_len_cells = Eigen::Vector3i(xc, yc, zc);
    // AABB is left at default (0,0,0)-(0,0,0); computed later in add_density_structure_all
  }

  // Read center position using Vector3d helper
  param_util::read_vector3d_direct(js, TOSTRING(xcnt), TOSTRING(ycnt), TOSTRING(zcnt), v3_pos_cnt);

  // Reject deprecated parameters (xlen/ylen/zlen)
  if (js.contains("xlen") || js.contains("ylen") || js.contains("zlen")) {
    THROW_ERROR("CheckerBoard3dParameters::assign_parameters: xlen/ylen/zlen are deprecated. Use xlen_interval_mult/ylen_interval_mult/zlen_interval_mult instead.");
  }

  // Read interval multipliers (required parameters)
  if (!js.contains("xlen_interval_mult")) {
    THROW_ERROR("CheckerBoard3dParameters::assign_parameters: xlen_interval_mult is required.");
  }
  if (!js.contains("ylen_interval_mult")) {
    THROW_ERROR("CheckerBoard3dParameters::assign_parameters: ylen_interval_mult is required.");
  }
  if (!js.contains("zlen_interval_mult")) {
    THROW_ERROR("CheckerBoard3dParameters::assign_parameters: zlen_interval_mult is required.");
  }

  int xlen_mult = js.at("xlen_interval_mult").get<int>();
  int ylen_mult = js.at("ylen_interval_mult").get<int>();
  int zlen_mult = js.at("zlen_interval_mult").get<int>();
  v3_len_interval_mult = Eigen::Vector3i(xlen_mult, ylen_mult, zlen_mult);

  // Validate interval multipliers (must be positive)
  if (v3_len_interval_mult.x() <= 0) {
    THROW_ERROR("CheckerBoard3dParameters::assign_parameters: xlen_interval_mult must be positive. value={}", v3_len_interval_mult.x());
  }
  if (v3_len_interval_mult.y() <= 0) {
    THROW_ERROR("CheckerBoard3dParameters::assign_parameters: ylen_interval_mult must be positive. value={}", v3_len_interval_mult.y());
  }
  if (v3_len_interval_mult.z() <= 0) {
    THROW_ERROR("CheckerBoard3dParameters::assign_parameters: zlen_interval_mult must be positive. value={}", v3_len_interval_mult.z());
  }

  // Read snap-to-grid flag (default: true)
  param_util::read_json_direct(js, TOSTRING(tf_snap_to_grid), tf_snap_to_grid);

  // Read region type (default: "aabb")
  param_util::read_json_direct(js, TOSTRING(region_type), region_type);
  if (region_type != "aabb" && region_type != "cylinder") {
    THROW_ERROR("CheckerBoard3dParameters::assign_parameters: "
      "region_type must be \"aabb\" or \"cylinder\". region_type={}", region_type);
  }

  // Read elliptic cylinder radii (required when region_type == "cylinder")
  if (region_type == "cylinder") {
    if (!js.contains("radius_x_meters") || !js.contains("radius_y_meters")) {
      THROW_ERROR("CheckerBoard3dParameters::assign_parameters: "
        "radius_x_meters and radius_y_meters are required when region_type is \"cylinder\".");
    }
    radius_x_meters = js.at("radius_x_meters").get<double>();
    radius_y_meters = js.at("radius_y_meters").get<double>();
    if (radius_x_meters <= 0.0 || radius_y_meters <= 0.0) {
      THROW_ERROR("CheckerBoard3dParameters::assign_parameters: "
        "radius_x_meters/radius_y_meters must be positive. values=({},{})",
        radius_x_meters, radius_y_meters);
    }
  }
}

//####################################
// EllipsoidParameters
//####################################
bool Grid3dVoxel::EllipsoidParameters::operator!=(const EllipsoidParameters& other) const
{
  #ifdef NODEBUG
  if( tf_exec != other.tf_exec ) return true;
  if( delta_density != other.delta_density ) return true;
  if( v3_pos_cnt != other.v3_pos_cnt ) return true;
  if( v3_length != other.v3_length ) return true;
  if( theta_x != other.theta_x ) return true;
  if( theta_y != other.theta_y ) return true;
  if( theta_z != other.theta_z ) return true;
  if( rotation_type != other.rotation_type ) return true;
  #else
  if( tf_exec != other.tf_exec ){ LOG_WARN("EllipsoidParameters : tf_exec differs"); return true; }
  if( delta_density != other.delta_density ){ LOG_WARN("EllipsoidParameters : delta_density differs"); return true; }
  if( v3_pos_cnt != other.v3_pos_cnt ){ LOG_WARN("EllipsoidParameters : v3_pos_cnt differs"); return true; }
  if( v3_length != other.v3_length ){ LOG_WARN("EllipsoidParameters : v3_length differs"); return true; }
  if( theta_x != other.theta_x ){ LOG_WARN("EllipsoidParameters : theta_x differs"); return true; }
  if( theta_y != other.theta_y ){ LOG_WARN("EllipsoidParameters : theta_y differs"); return true; }
  if( theta_z != other.theta_z ){ LOG_WARN("EllipsoidParameters : theta_z differs"); return true; }
  if( rotation_type != other.rotation_type ){ LOG_WARN("EllipsoidParameters : rotation_type differs"); return true; }
  #endif
  return false;
}

void Grid3dVoxel::EllipsoidParameters::assign_parameters(const nlohmann::json& js)
{
  LOG_INFO("");

  // Read basic parameters
  param_util::read_json_direct(js, TOSTRING(tf_exec), tf_exec);
  param_util::read_json_direct(js, TOSTRING(name), name);
  param_util::read_json_direct(js, TOSTRING(delta_density), delta_density);

  // Read center position and length using Vector3d helper
  param_util::read_vector3d_direct(js, TOSTRING(xcnt), TOSTRING(ycnt), TOSTRING(zcnt), v3_pos_cnt);
  param_util::read_vector3d_direct(js, TOSTRING(xlen), TOSTRING(ylen), TOSTRING(zlen), v3_length);

  // Read rotation angles
  double theta_x_deg = 0.0, theta_y_deg = 0.0, theta_z_deg = 0.0;
  param_util::read_json_direct(js, TOSTRING(theta_x_deg), theta_x_deg);
  param_util::read_json_direct(js, TOSTRING(theta_y_deg), theta_y_deg);
  param_util::read_json_direct(js, TOSTRING(theta_z_deg), theta_z_deg);
  theta_x.setDegree(theta_x_deg);
  theta_y.setDegree(theta_y_deg);
  theta_z.setDegree(theta_z_deg);

  // Read rotation type
  std::string rotation_type_str;
  param_util::read_json_direct(js, TOSTRING(rotation_type), rotation_type_str);
  if(rotation_type_str == "LOCAL" || rotation_type_str == "local" || rotation_type_str == "Local"){
    rotation_type = angle_util::Rotation3dType::LOCAL;
  }
  else if(rotation_type_str == "GLOBAL" || rotation_type_str == "global" || rotation_type_str == "Global"){
    rotation_type = angle_util::Rotation3dType::GLOBAL;
  }
  else{
    THROW_ERROR("EllipsoidParameters::assign_parameters: Invalid rotation_type. Must be LOCAL or GLOBAL. rotation_type={}", rotation_type_str);
  }

}

//####################################
// CylinderParameters
//####################################
bool Grid3dVoxel::CylinderParameters::operator!=(const CylinderParameters& other) const
{
  #ifdef NODEBUG
  if( tf_exec != other.tf_exec ) return true;
  if( delta_density != other.delta_density ) return true;
  if( v3_pos_cnt != other.v3_pos_cnt ) return true;
  if( v3_length != other.v3_length ) return true;
  if( theta_x != other.theta_x ) return true;
  if( theta_y != other.theta_y ) return true;
  if( theta_z != other.theta_z ) return true;
  if( rotation_type != other.rotation_type ) return true;
  #else
  if( tf_exec != other.tf_exec ){ LOG_WARN("CylinderParameters : tf_exec differs"); return true; }
  if( delta_density != other.delta_density ){ LOG_WARN("CylinderParameters : delta_density differs"); return true; }
  if( v3_pos_cnt != other.v3_pos_cnt ){ LOG_WARN("CylinderParameters : v3_pos_cnt differs"); return true; }
  if( v3_length != other.v3_length ){ LOG_WARN("CylinderParameters : v3_length differs"); return true; }
  if( theta_x != other.theta_x ){ LOG_WARN("CylinderParameters : theta_x differs"); return true; }
  if( theta_y != other.theta_y ){ LOG_WARN("CylinderParameters : theta_y differs"); return true; }
  if( theta_z != other.theta_z ){ LOG_WARN("CylinderParameters : theta_z differs"); return true; }
  if( rotation_type != other.rotation_type ){ LOG_WARN("CylinderParameters : rotation_type differs"); return true; }
  #endif
  return false;
}

void Grid3dVoxel::CylinderParameters::assign_parameters(const nlohmann::json& js)
{
  LOG_INFO("");

  // Read basic parameters
  param_util::read_json_direct(js, TOSTRING(tf_exec), tf_exec);
  param_util::read_json_direct(js, TOSTRING(name), name);
  param_util::read_json_direct(js, TOSTRING(delta_density), delta_density);

  // Read center position and length using Vector3d helper
  param_util::read_vector3d_direct(js, TOSTRING(xcnt), TOSTRING(ycnt), TOSTRING(zcnt), v3_pos_cnt);
  param_util::read_vector3d_direct(js, TOSTRING(xlen), TOSTRING(ylen), TOSTRING(zlen), v3_length);
  if(v3_length.x() <= 0) THROW_ERROR("CylinderParameters::assign_parameters: xlen must be positive. xlen={}", v3_length.x());
  if(v3_length.y() <= 0) THROW_ERROR("CylinderParameters::assign_parameters: ylen must be positive. ylen={}", v3_length.y());
  if(v3_length.z() <= 0) THROW_ERROR("CylinderParameters::assign_parameters: zlen must be positive. zlen={}", v3_length.z());

  // Read rotation angles
  double theta_x_deg = 0.0, theta_y_deg = 0.0, theta_z_deg = 0.0;
  param_util::read_json_direct(js, TOSTRING(theta_x_deg), theta_x_deg);
  param_util::read_json_direct(js, TOSTRING(theta_y_deg), theta_y_deg);
  param_util::read_json_direct(js, TOSTRING(theta_z_deg), theta_z_deg);
  theta_x.setDegree(theta_x_deg);
  theta_y.setDegree(theta_y_deg);
  theta_z.setDegree(theta_z_deg);

  // Read rotation type
  std::string rotation_type_str;
  param_util::read_json_direct(js, TOSTRING(rotation_type), rotation_type_str);
  if(rotation_type_str == "LOCAL" || rotation_type_str == "local" || rotation_type_str == "Local"){
    rotation_type = angle_util::Rotation3dType::LOCAL;
  }
  else if(rotation_type_str == "GLOBAL" || rotation_type_str == "global" || rotation_type_str == "Global"){
    rotation_type = angle_util::Rotation3dType::GLOBAL;
  }
  else{
    THROW_ERROR("CylinderParameters::assign_parameters: Invalid rotation_type. Must be LOCAL or GLOBAL. rotation_type={}", rotation_type_str);
  }

}

//####################################
// CuboidParameters
//####################################
bool Grid3dVoxel::CuboidParameters::operator!=(const CuboidParameters& other) const
{
  #ifdef NODEBUG
  if( tf_exec != other.tf_exec ) return true;
  if( delta_density != other.delta_density ) return true;
  if( v3_pos_cnt != other.v3_pos_cnt ) return true;
  if( v3_length != other.v3_length ) return true;
  if( theta_x != other.theta_x ) return true;
  if( theta_y != other.theta_y ) return true;
  if( theta_z != other.theta_z ) return true;
  if( rotation_type != other.rotation_type ) return true;
  #else
  if( tf_exec != other.tf_exec ){ LOG_WARN("CuboidParameters : tf_exec differs"); return true; }
  if( delta_density != other.delta_density ){ LOG_WARN("CuboidParameters : delta_density differs"); return true; }
  if( v3_pos_cnt != other.v3_pos_cnt ){ LOG_WARN("CuboidParameters : v3_pos_cnt differs"); return true; }
  if( v3_length != other.v3_length ){ LOG_WARN("CuboidParameters : v3_length differs"); return true; }
  if( theta_x != other.theta_x ){ LOG_WARN("CuboidParameters : theta_x differs"); return true; }
  if( theta_y != other.theta_y ){ LOG_WARN("CuboidParameters : theta_y differs"); return true; }
  if( theta_z != other.theta_z ){ LOG_WARN("CuboidParameters : theta_z differs"); return true; }
  if( rotation_type != other.rotation_type ){ LOG_WARN("CuboidParameters : rotation_type differs"); return true; }
  #endif
  return false;
}

void Grid3dVoxel::CuboidParameters::assign_parameters(const nlohmann::json& js)
{
  LOG_INFO("");

  // Read basic parameters
  param_util::read_json_direct(js, TOSTRING(tf_exec), tf_exec);
  param_util::read_json_direct(js, TOSTRING(name), name);
  param_util::read_json_direct(js, TOSTRING(delta_density), delta_density);

  // Read center position and length using Vector3d helper
  param_util::read_vector3d_direct(js, TOSTRING(xcnt), TOSTRING(ycnt), TOSTRING(zcnt), v3_pos_cnt);
  param_util::read_vector3d_direct(js, TOSTRING(xlen), TOSTRING(ylen), TOSTRING(zlen), v3_length);
  if(v3_length.x() <= 0) THROW_ERROR("CuboidParameters::assign_parameters: xlen must be positive. xlen={}", v3_length.x());
  if(v3_length.y() <= 0) THROW_ERROR("CuboidParameters::assign_parameters: ylen must be positive. ylen={}", v3_length.y());
  if(v3_length.z() <= 0) THROW_ERROR("CuboidParameters::assign_parameters: zlen must be positive. zlen={}", v3_length.z());

  // Read rotation angles
  double theta_x_deg = 0.0, theta_y_deg = 0.0, theta_z_deg = 0.0;
  param_util::read_json_direct(js, TOSTRING(theta_x_deg), theta_x_deg);
  param_util::read_json_direct(js, TOSTRING(theta_y_deg), theta_y_deg);
  param_util::read_json_direct(js, TOSTRING(theta_z_deg), theta_z_deg);
  theta_x.setDegree(theta_x_deg);
  theta_y.setDegree(theta_y_deg);
  theta_z.setDegree(theta_z_deg);

  // Read rotation type
  std::string rotation_type_str;
  param_util::read_json_direct(js, TOSTRING(rotation_type), rotation_type_str);
  if(rotation_type_str == "LOCAL" || rotation_type_str == "local" || rotation_type_str == "Local"){
    rotation_type = angle_util::Rotation3dType::LOCAL;
  }
  else if(rotation_type_str == "GLOBAL" || rotation_type_str == "global" || rotation_type_str == "Global"){
    rotation_type = angle_util::Rotation3dType::GLOBAL;
  }
  else{
    THROW_ERROR("CuboidParameters::assign_parameters: Invalid rotation_type. Must be LOCAL or GLOBAL. rotation_type={}", rotation_type_str);
  }

}

//####################################
// MergeParameters
//####################################
bool Grid3dVoxel::MergeParameters::operator!=(const MergeParameters& other) const
{
  #ifdef NODEBUG
  if( tf_exec != other.tf_exec ) return true;
  if( x_merge_center != other.x_merge_center ) return true;
  if( y_merge_center != other.y_merge_center ) return true;
  if( z_merge_center != other.z_merge_center ) return true;
  if( x_merge_factor != other.x_merge_factor ) return true;
  if( y_merge_factor != other.y_merge_factor ) return true;
  if( z_merge_factor != other.z_merge_factor ) return true;
  #else
  if( tf_exec != other.tf_exec ){ LOG_WARN("MergeParameters : tf_exec differs"); return true; }
  if( x_merge_center != other.x_merge_center ){ LOG_WARN("MergeParameters : x_merge_center differs"); return true; }
  if( y_merge_center != other.y_merge_center ){ LOG_WARN("MergeParameters : y_merge_center differs"); return true; }
  if( z_merge_center != other.z_merge_center ){ LOG_WARN("MergeParameters : z_merge_center differs"); return true; }
  if( x_merge_factor != other.x_merge_factor ){ LOG_WARN("MergeParameters : x_merge_factor differs"); return true; }
  if( y_merge_factor != other.y_merge_factor ){ LOG_WARN("MergeParameters : y_merge_factor differs"); return true; }
  if( z_merge_factor != other.z_merge_factor ){ LOG_WARN("MergeParameters : z_merge_factor differs"); return true; }
  #endif
  return false;
}

void Grid3dVoxel::MergeParameters::assign_parameters(const nlohmann::json& js)
{
  LOG_INFO("");

  // Read basic parameters
  param_util::read_json_direct(js, TOSTRING(tf_exec), tf_exec);
  param_util::read_json_direct(js, TOSTRING(name), name);

  // Read merge center coordinates
  param_util::read_json_direct(js, TOSTRING(x_merge_center), x_merge_center);
  param_util::read_json_direct(js, TOSTRING(y_merge_center), y_merge_center);
  param_util::read_json_direct(js, TOSTRING(z_merge_center), z_merge_center);

  // Read merge factors
  param_util::read_json_direct(js, TOSTRING(x_merge_factor), x_merge_factor);
  if(x_merge_factor < 2) THROW_ERROR("MergeParameters::assign_parameters: x_merge_factor must be >= 2. x_merge_factor={}", x_merge_factor);
  param_util::read_json_direct(js, TOSTRING(y_merge_factor), y_merge_factor);
  if(y_merge_factor < 2) THROW_ERROR("MergeParameters::assign_parameters: y_merge_factor must be >= 2. y_merge_factor={}", y_merge_factor);
  param_util::read_json_direct(js, TOSTRING(z_merge_factor), z_merge_factor);
  if(z_merge_factor < 2) THROW_ERROR("MergeParameters::assign_parameters: z_merge_factor must be >= 2. z_merge_factor={}", z_merge_factor);
}

//####################################
// ReconstVoxelsParameters
//####################################
bool Grid3dVoxel::ReconstVoxelsParameters::operator!=(const ReconstVoxelsParameters& other) const
{
  #ifdef NODEBUG
  if( tf_aabb != other.tf_aabb ) return true;
  if( x_aabb_cnt != other.x_aabb_cnt ) return true;
  if( y_aabb_cnt != other.y_aabb_cnt ) return true;
  if( x_aabb_meters != other.x_aabb_meters ) return true;
  if( y_aabb_meters != other.y_aabb_meters ) return true;
  if( aabb_zmin_mode != other.aabb_zmin_mode ) return true;
  if( aabb_zmin_value != other.aabb_zmin_value ) return true;
  if( aabb_zmax != other.aabb_zmax ) return true;
  if( tf_cylinder != other.tf_cylinder ) return true;
  if( x_cyl_cnt != other.x_cyl_cnt ) return true;
  if( y_cyl_cnt != other.y_cyl_cnt ) return true;
  if( cylinder_radius_x_meters != other.cylinder_radius_x_meters ) return true;
  if( cylinder_radius_y_meters != other.cylinder_radius_y_meters ) return true;
  #else
  if( tf_aabb != other.tf_aabb ){ LOG_WARN("ReconstVoxelsParameters : tf_aabb differs"); return true; }
  if( x_aabb_cnt != other.x_aabb_cnt ){ LOG_WARN("ReconstVoxelsParameters : x_aabb_cnt differs"); return true; }
  if( y_aabb_cnt != other.y_aabb_cnt ){ LOG_WARN("ReconstVoxelsParameters : y_aabb_cnt differs"); return true; }
  if( x_aabb_meters != other.x_aabb_meters ){ LOG_WARN("ReconstVoxelsParameters : x_aabb_meters differs"); return true; }
  if( y_aabb_meters != other.y_aabb_meters ){ LOG_WARN("ReconstVoxelsParameters : y_aabb_meters differs"); return true; }
  if( aabb_zmin_mode != other.aabb_zmin_mode ){ LOG_WARN("ReconstVoxelsParameters : aabb_zmin_mode differs"); return true; }
  if( aabb_zmin_value != other.aabb_zmin_value ){ LOG_WARN("ReconstVoxelsParameters : aabb_zmin_value differs"); return true; }
  if( aabb_zmax != other.aabb_zmax ){ LOG_WARN("ReconstVoxelsParameters : aabb_zmax differs"); return true; }
  if( tf_cylinder != other.tf_cylinder ){ LOG_WARN("ReconstVoxelsParameters : tf_cylinder differs"); return true; }
  if( x_cyl_cnt != other.x_cyl_cnt ){ LOG_WARN("ReconstVoxelsParameters : x_cyl_cnt differs"); return true; }
  if( y_cyl_cnt != other.y_cyl_cnt ){ LOG_WARN("ReconstVoxelsParameters : y_cyl_cnt differs"); return true; }
  if( cylinder_radius_x_meters != other.cylinder_radius_x_meters ){ LOG_WARN("ReconstVoxelsParameters : cylinder_radius_x_meters differs"); return true; }
  if( cylinder_radius_y_meters != other.cylinder_radius_y_meters ){ LOG_WARN("ReconstVoxelsParameters : cylinder_radius_y_meters differs"); return true; }
  #endif
  return false;
}

void Grid3dVoxel::ReconstVoxelsParameters::assign_parameters(const nlohmann::json& js)
{
  LOG_INFO("");

  // Read AABB mask parameters
  param_util::read_json_direct(js, TOSTRING(tf_aabb), tf_aabb);
  param_util::read_json_direct(js, TOSTRING(x_aabb_cnt), x_aabb_cnt);
  param_util::read_json_direct(js, TOSTRING(y_aabb_cnt), y_aabb_cnt);
  param_util::read_json_direct(js, TOSTRING(x_aabb_meters), x_aabb_meters);
  param_util::read_json_direct(js, TOSTRING(y_aabb_meters), y_aabb_meters);
  param_util::read_json_direct(js, TOSTRING(aabb_zmin_mode), aabb_zmin_mode);
  param_util::read_json_direct(js, TOSTRING(aabb_zmin_value), aabb_zmin_value);
  param_util::read_json_direct(js, TOSTRING(aabb_zmax), aabb_zmax);
  if(tf_aabb && x_aabb_meters <= 0.0) THROW_ERROR("ReconstVoxelsParameters::assign_parameters: x_aabb_meters must be > 0 when tf_aabb is true. x_aabb_meters={}", x_aabb_meters);
  if(tf_aabb && y_aabb_meters <= 0.0) THROW_ERROR("ReconstVoxelsParameters::assign_parameters: y_aabb_meters must be > 0 when tf_aabb is true. y_aabb_meters={}", y_aabb_meters);
  if(aabb_zmin_mode != "g3vox_zmin" && aabb_zmin_mode != "manual") THROW_ERROR("ReconstVoxelsParameters::assign_parameters: aabb_zmin_mode must be \"g3vox_zmin\" or \"manual\". aabb_zmin_mode={}", aabb_zmin_mode);
  if(aabb_zmin_mode == "manual" && aabb_zmin_value >= aabb_zmax) THROW_ERROR("ReconstVoxelsParameters::assign_parameters: aabb_zmin_value must be < aabb_zmax. aabb_zmin_value={}, aabb_zmax={}", aabb_zmin_value, aabb_zmax);

  // Read cylinder mask parameters
  param_util::read_json_direct(js, TOSTRING(tf_cylinder), tf_cylinder);
  param_util::read_json_direct(js, TOSTRING(x_cyl_cnt), x_cyl_cnt);
  param_util::read_json_direct(js, TOSTRING(y_cyl_cnt), y_cyl_cnt);
  param_util::read_json_direct(js, TOSTRING(cylinder_radius_x_meters), cylinder_radius_x_meters);
  param_util::read_json_direct(js, TOSTRING(cylinder_radius_y_meters), cylinder_radius_y_meters);
  if(tf_cylinder && cylinder_radius_x_meters <= 0.0) THROW_ERROR("ReconstVoxelsParameters::assign_parameters: cylinder_radius_x_meters must be > 0 when tf_cylinder is true. cylinder_radius_x_meters={}", cylinder_radius_x_meters);
  if(tf_cylinder && cylinder_radius_y_meters <= 0.0) THROW_ERROR("ReconstVoxelsParameters::assign_parameters: cylinder_radius_y_meters must be > 0 when tf_cylinder is true. cylinder_radius_y_meters={}", cylinder_radius_y_meters);
}
