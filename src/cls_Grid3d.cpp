// cls_Grid3d.cpp
#include "cls_Grid3d.hpp"

#include <limits>

#include "ns_myapp.hpp"
#include "ns_mymacro.hpp"
#include "ns_mylogger.hpp"
#include "spdlog_pch.hpp"

//##############################################
//##############################################
// class Grid3d
//##############################################
//##############################################

Grid3d::Grid3d( const Grid1d &x_axis_in, const Grid1d &y_axis_in, const Grid1d &z_axis_in )
  : name("g3_" + x_axis_in.get_name()
         + "_" + y_axis_in.get_name()
         + "_" + z_axis_in.get_name() )
  , x_axis(x_axis_in)
  , y_axis(y_axis_in)
  , z_axis(z_axis_in)
{
  // cache grid AABB so ray tracing skips per-ray AABB construction
  rebuild_cached_aabb3d();
}

/// @brief Copy assignment operator.
Grid3d& Grid3d::operator=(const Grid3d& other)
{
  if(this == &other) return *this;
  name = other.name;
  x_axis = other.x_axis;
  y_axis = other.y_axis;
  z_axis = other.z_axis;
  // axes changed — keep cached AABB consistent
  rebuild_cached_aabb3d();
  return *this;
}

/// @brief Inequality operator.
/// @remark Name is intentionally ignored.
bool Grid3d::operator!=(const Grid3d& other) const
{
  // Name is intentionally ignored.
  // if (name != other.name) return true;
  if (x_axis != other.x_axis) return true;
  if (y_axis != other.y_axis) return true;
  if (z_axis != other.z_axis) return true;
  return false;
}

// Rebuild cached AABB3d from current axis min/max values.
// Called after any axis change so ray tracing always uses up-to-date bounds.
void Grid3d::rebuild_cached_aabb3d()
{
  cached_aabb3d_ = AABB3d(
    Eigen::Vector3d(x_axis.get_min(), y_axis.get_min(), z_axis.get_min()),
    Eigen::Vector3d(x_axis.get_max(), y_axis.get_max(), z_axis.get_max()));
}

// Build DDA traversal parameters from a valid ray-AABB intersection result.
// Extracts the common ~80 lines shared by get_hit_voxel_index() variants
// and traverse_ray_with_pathlength(), eliminating code duplication.
Grid3d::VoxelTraversalParameters3d Grid3d::build_traversal_params(
    const Ray3d &ray3d, double eps, double tmin, double tmax) const
{
  VoxelTraversalParameters3d prm;

  const auto [vx_sign, vy_sign, vz_sign] = ray3d.get_sign();
  prm.vx_sign = vx_sign;
  prm.vy_sign = vy_sign;
  prm.vz_sign = vz_sign;

  prm.step_x = get_x_interval() * vx_sign;
  prm.step_y = get_y_interval() * vy_sign;
  prm.step_z = get_z_interval() * vz_sign;

  const double t_eps =
      std::max({get_x_interval(), get_y_interval(), get_z_interval()}) * eps;

  const Eigen::Vector3d v3_pos_min = ray3d.pos() + (tmin + t_eps) * ray3d.dir();
  const Eigen::Vector3d v3_pos_max = ray3d.pos() + (tmax - t_eps) * ray3d.dir();

  if (tmin <= 0 && tmax >= 0) {
    prm.x_start = ray3d.x();
    prm.y_start = ray3d.y();
    prm.z_start = ray3d.z();
  }
  if (tmin > 0 && tmax > 0) {
    prm.x_start = v3_pos_min.x();
    prm.y_start = v3_pos_min.y();
    prm.z_start = v3_pos_min.z();
  }

  int ix_start = get_ix(prm.x_start);
  int iy_start = get_iy(prm.y_start);
  int iz_start = get_iz(prm.z_start);

  check_ix_inside(ix_start);
  check_iy_inside(iy_start);
  check_iz_inside(iz_start);

  int ix_end = get_ix(v3_pos_max.x());
  int iy_end = get_iy(v3_pos_max.y());
  int iz_end = get_iz(v3_pos_max.z());

  if (ix_end == Grid1d::OUT_OF_RANGE_UPPER) ix_end = get_nbinx() - 1;
  if (ix_end == Grid1d::OUT_OF_RANGE_LOWER) ix_end = 0;
  if (iy_end == Grid1d::OUT_OF_RANGE_UPPER) iy_end = get_nbiny() - 1;
  if (iy_end == Grid1d::OUT_OF_RANGE_LOWER) iy_end = 0;
  if (iz_end == Grid1d::OUT_OF_RANGE_UPPER) iz_end = get_nbinz() - 1;
  if (iz_end == Grid1d::OUT_OF_RANGE_LOWER) iz_end = 0;

  double x_offset = 0.0;
  if (vx_sign > 0) x_offset = get_xup(prm.x_start) - prm.x_start;
  if (vx_sign < 0) x_offset = prm.x_start - get_xlow(prm.x_start);

  double y_offset = 0.0;
  if (vy_sign > 0) y_offset = get_yup(prm.y_start) - prm.y_start;
  if (vy_sign < 0) y_offset = prm.y_start - get_ylow(prm.y_start);

  double z_offset = 0.0;
  if (vz_sign > 0) z_offset = get_zup(prm.z_start) - prm.z_start;
  if (vz_sign < 0) z_offset = prm.z_start - get_zlow(prm.z_start);

  const Eigen::Vector3d v3_dir = ray3d.dir();
  const double dir_norm = v3_dir.norm();
  constexpr double kDirEps = 1.0e-12;
  const double inf = std::numeric_limits<double>::infinity();

  const double path_factor_dx =
      (std::abs(v3_dir.x()) > kDirEps) ? (dir_norm / std::abs(v3_dir.x())) : inf;
  const double path_factor_dy =
      (std::abs(v3_dir.y()) > kDirEps) ? (dir_norm / std::abs(v3_dir.y())) : inf;
  const double path_factor_dz =
      (std::abs(v3_dir.z()) > kDirEps) ? (dir_norm / std::abs(v3_dir.z())) : inf;

  prm.t_delta_x = (vx_sign == 0.0) ? inf : prm.step_x * path_factor_dx;
  prm.t_delta_y = (vy_sign == 0.0) ? inf : prm.step_y * path_factor_dy;
  prm.t_delta_z = (vz_sign == 0.0) ? inf : prm.step_z * path_factor_dz;

  prm.t_max_x0 = (vx_sign == 0.0) ? inf : x_offset * path_factor_dx * vx_sign;
  prm.t_max_y0 = (vy_sign == 0.0) ? inf : y_offset * path_factor_dy * vy_sign;
  prm.t_max_z0 = (vz_sign == 0.0) ? inf : z_offset * path_factor_dz * vz_sign;

  prm.manhattan_distance =
        std::abs(ix_end - ix_start)
      + std::abs(iy_end - iy_start)
      + std::abs(iz_end - iz_start);

  return prm;
}

double Grid3d::get_xcnt( const double value ) const
{
  const int index = get_ix(value);
  check_ix_inside(index);
  return get_xcnt(index);
}

double Grid3d::get_ycnt( const double value ) const
{
  const int index = get_iy(value);
  check_iy_inside(index);
  return get_ycnt(index);
}

double Grid3d::get_zcnt( const double value ) const
{
  const int index = get_iz(value);
  check_iz_inside(index);
  return get_zcnt(index);
}

double Grid3d::get_xup( const double value ) const
{
  const int index = get_ix(value);
  check_ix_inside(index);
  return x_axis.get_upper_value(index);
}

double Grid3d::get_yup( const double value ) const
{
  const int index = get_iy(value);
  check_iy_inside(index);
  return y_axis.get_upper_value(index);
}

double Grid3d::get_zup( const double value ) const
{
  const int index = get_iz(value);
  check_iz_inside(index);
  return z_axis.get_upper_value(index);
}

double Grid3d::get_xlow( const double value ) const
{
  const int index = get_ix(value);
  check_ix_inside(index);
  return x_axis.get_lower_value(index);
}

double Grid3d::get_ylow( const double value ) const
{
  const int index = get_iy(value);
  check_iy_inside(index);
  return y_axis.get_lower_value(index);
}

double Grid3d::get_zlow( const double value ) const
{
  const int index = get_iz(value);
  check_iz_inside(index);
  return z_axis.get_lower_value(index);
}

void Grid3d::check_ix_inside( const int ix ) const
{
  if(is_ix_inside(ix)==true) return;
  THROW_ERROR("Grid3d::check_ix_inside: ix out of range. ix={}", ix);
}

void Grid3d::check_iy_inside( const int iy ) const
{
  if(is_iy_inside(iy)==true) return;
  THROW_ERROR("Grid3d::check_iy_inside: iy out of range. iy={}", iy);
}

void Grid3d::check_iz_inside( const int iz ) const
{
  if(is_iz_inside(iz)==true) return;
  THROW_ERROR("Grid3d::check_iz_inside: iz out of range. iz={}", iz);
}

void Grid3d::check_ix_inside( const int ix, const std::string &filename, const std::string &function_name, const int line ) const
{
  if(is_ix_inside(ix)==false){
    const int nbinx = get_nbinx();
    LOG_ERROR("ix out of range. file={}, function={}, line={}, ix={}, valid=[0, {}]"
    ,filename,function_name,line,ix,nbinx-1);
    THROW_ERROR("Grid3d::check_ix_inside: ix out of range. ix={}", ix);
  }
}

void Grid3d::check_iy_inside( const int iy, const std::string &filename, const std::string &function_name, const int line ) const
{
  if(is_iy_inside(iy)==false){
    const int nbiny = get_nbiny();
    LOG_ERROR("iy out of range. file={}, function={}, line={}, iy={}, valid=[0, {}]"
    ,filename,function_name,line,iy,nbiny-1);
    THROW_ERROR("Grid3d::check_iy_inside: iy out of range. iy={}", iy);
  }
}

void Grid3d::check_iz_inside( const int iz, const std::string &filename, const std::string &function_name, const int line ) const
{
  if(is_iz_inside(iz)==false){
    const int nbinz = get_nbinz();
    LOG_ERROR("iz out of range. file={}, function={}, line={}, iz={}, valid=[0, {}]"
    ,filename,function_name,line,iz,nbinz-1);
    THROW_ERROR("Grid3d::check_iz_inside: iz out of range. iz={}", iz);
  }
}

/// @brief Convert coordinates to (ix, iy, iz).
Grid3d::Ixiyiz Grid3d::get_ixiyiz(
  const double x_in, const double y_in, const double z_in ) const
{
  const int ix = get_ix(x_in);
  const int iy = get_iy(y_in);
  const int iz = get_iz(z_in);
  check_ix_inside(ix);
  check_iy_inside(iy);
  check_iz_inside(iz);
  return {ix,iy,iz};
}


Eigen::Vector3d Grid3d::get_v3_AABB_min( const int ix, const int iy, const int iz ) const
{
  const Eigen::Vector3d v3_AABB_min(
    get_xlow(ix)
  , get_ylow(iy)
  , get_zlow(iz));
  return v3_AABB_min;
}

Eigen::Vector3d Grid3d::get_v3_AABB_max( const int ix, const int iy, const int iz ) const
{
  const Eigen::Vector3d v3_AABB_max(
    get_xup(ix)
  , get_yup(iy)
  , get_zup(iz));
  return v3_AABB_max;
}

AABB3d Grid3d::get_AABB3d( const int ix, const int iy, const int iz ) const
{
  const Eigen::Vector3d v3_AABB_min = get_v3_AABB_min(ix,iy,iz);
  const Eigen::Vector3d v3_AABB_max = get_v3_AABB_max(ix,iy,iz);
  const AABB3d aabb3d(v3_AABB_min,v3_AABB_max);
  return aabb3d;
}

// [Legacy] Original implementation — constructs a local AABB on every call.
// Retained for regression testing against the optimized version.
std::vector<Grid3d::Ixiyiz>
  Grid3d::get_hit_voxel_index_legacy( const Ray3d &ray3d, const double eps ) const
{
  std::vector<Ixiyiz> vec_ixiyiz;
  VoxelTraversalParameters3d prm;

  const Eigen::Vector3d v3_AABB_min(get_xmin(),get_ymin(),get_zmin());
  const Eigen::Vector3d v3_AABB_max(get_xmax(),get_ymax(),get_zmax());
  const AABB3d aabb3d(v3_AABB_min,v3_AABB_max);

  // Algorithm overview:
  // 1) Intersect the ray with the grid AABB to get entry/exit parameters.
  // 2) Choose a start point inside the AABB (origin or entry point).
  // 3) Compute start/end voxel indices and step directions.
  // 4) Compute DDA increments (t_delta) and first boundary distances (t_max).
  // 5) Traverse voxels until the ray exits the grid.

  auto [tf_hit_volume,tmin,tmax] = ray3d.is_intersect(aabb3d);

  if( tf_hit_volume == false ) return vec_ixiyiz;
  if( tmin < 0 && tmax < 0 ) return vec_ixiyiz;

  const auto [vx_sign, vy_sign, vz_sign] = ray3d.get_sign();

  prm.vx_sign = vx_sign;
  prm.vy_sign = vy_sign;
  prm.vz_sign = vz_sign;

  const double step_x = get_x_interval() * vx_sign;
  const double step_y = get_y_interval() * vy_sign;
  const double step_z = get_z_interval() * vz_sign;

  prm.step_x = step_x;
  prm.step_y = step_y;
  prm.step_z = step_z;

  const double t_eps = std::max({get_x_interval(), get_y_interval(), get_z_interval()}) * eps;

  const Eigen::Vector3d v3_pos_min = ray3d.pos() + ( tmin + t_eps ) * ray3d.dir();
  const Eigen::Vector3d v3_pos_max = ray3d.pos() + ( tmax - t_eps ) * ray3d.dir();
  double x_start, y_start, z_start;
  if( tmin <= 0 && tmax >= 0 ){
    x_start = ray3d.x();
    y_start = ray3d.y();
    z_start = ray3d.z();
  }
  if( tmin > 0 && tmax > 0 ){
    x_start = v3_pos_min.x();
    y_start = v3_pos_min.y();
    z_start = v3_pos_min.z();
  }

  prm.x_start = x_start;
  prm.y_start = y_start;
  prm.z_start = z_start;

  int ix_start = get_ix(x_start);
  int iy_start = get_iy(y_start);
  int iz_start = get_iz(z_start);

  check_ix_inside(ix_start);
  check_iy_inside(iy_start);
  check_iz_inside(iz_start);

  double x_end = v3_pos_max.x();
  double y_end = v3_pos_max.y();
  double z_end = v3_pos_max.z();

  int ix_end = get_ix(x_end);
  int iy_end = get_iy(y_end);
  int iz_end = get_iz(z_end);

  // End point might be on the boundary (OUT_OF_RANGE_UPPER) due to floating-point precision
  // Clamp to valid range instead of throwing error
  if( ix_end == Grid1d::OUT_OF_RANGE_UPPER ) ix_end = get_nbinx() - 1;
  if( ix_end == Grid1d::OUT_OF_RANGE_LOWER ) ix_end = 0;
  if( iy_end == Grid1d::OUT_OF_RANGE_UPPER ) iy_end = get_nbiny() - 1;
  if( iy_end == Grid1d::OUT_OF_RANGE_LOWER ) iy_end = 0;
  if( iz_end == Grid1d::OUT_OF_RANGE_UPPER ) iz_end = get_nbinz() - 1;
  if( iz_end == Grid1d::OUT_OF_RANGE_LOWER ) iz_end = 0;

  double x_offset = 0.0;
  if( vx_sign > 0 ) x_offset = get_xup(x_start) - x_start;
  if( vx_sign < 0 ) x_offset = x_start - get_xlow(x_start);

  double y_offset = 0.0;
  if( vy_sign > 0 ) y_offset = get_yup(y_start) - y_start;
  if( vy_sign < 0 ) y_offset = y_start - get_ylow(y_start);

  double z_offset = 0.0;
  if( vz_sign > 0 ) z_offset = get_zup(z_start) - z_start;
  if( vz_sign < 0 ) z_offset = z_start - get_zlow(z_start);

  const Eigen::Vector3d v3_dir = ray3d.dir();
  const double dir_norm = v3_dir.norm();
  constexpr double kDirEps = 1.0e-12;
  const double inf = std::numeric_limits<double>::infinity();

  const double path_factor_dx =
    (std::abs(v3_dir.x()) > kDirEps) ? (dir_norm / std::abs(v3_dir.x())) : inf;
  const double path_factor_dy =
    (std::abs(v3_dir.y()) > kDirEps) ? (dir_norm / std::abs(v3_dir.y())) : inf;
  const double path_factor_dz =
    (std::abs(v3_dir.z()) > kDirEps) ? (dir_norm / std::abs(v3_dir.z())) : inf;

  const double t_delta_x = (vx_sign == 0.0) ? inf : step_x * path_factor_dx;
  const double t_delta_y = (vy_sign == 0.0) ? inf : step_y * path_factor_dy;
  const double t_delta_z = (vz_sign == 0.0) ? inf : step_z * path_factor_dz;

  prm.t_delta_x = t_delta_x;
  prm.t_delta_y = t_delta_y;
  prm.t_delta_z = t_delta_z;

  const double t_max_x0 = (vx_sign == 0.0) ? inf : x_offset * path_factor_dx * vx_sign;
  const double t_max_y0 = (vy_sign == 0.0) ? inf : y_offset * path_factor_dy * vy_sign;
  const double t_max_z0 = (vz_sign == 0.0) ? inf : z_offset * path_factor_dz * vz_sign;

  prm.t_max_x0 = t_max_x0;
  prm.t_max_y0 = t_max_y0;
  prm.t_max_z0 = t_max_z0;

  const int manhattan_distance =
      abs( ix_end - ix_start )
    + abs( iy_end - iy_start )
    + abs( iz_end - iz_start );

  prm.manhattan_distance = manhattan_distance;

  return get_hit_index_by_ray_tracing_algorithm(prm);
}

// [Legacy] Original BL_max overload — performs double ray-AABB intersection.
// First intersection for BL_max check, second inside get_hit_voxel_index_legacy().
// Retained for regression testing.
std::vector<Grid3d::Ixiyiz>
  Grid3d::get_hit_voxel_index_legacy( const Ray3d &ray3d, const double BL_max, const double eps ) const
{
  if( BL_max <= 0.0 ){
    return get_hit_voxel_index_legacy(ray3d, eps);
  }

  const Eigen::Vector3d v3_AABB_min(get_xmin(),get_ymin(),get_zmin());
  const Eigen::Vector3d v3_AABB_max(get_xmax(),get_ymax(),get_zmax());
  const AABB3d aabb3d(v3_AABB_min,v3_AABB_max);
  auto [tf_hit_volume,tmin,tmax] = ray3d.is_intersect(aabb3d);

  if( tf_hit_volume && tmax > BL_max ){
    LOG_WARN_ND("Ray beam length ({:.1f}m) exceeds BL_max ({:.1f}m). Ray3d={}", tmax, BL_max, ray3d.to_string());
  }

  return get_hit_voxel_index_legacy(ray3d, eps);
}

/// @brief Helper function for 3D DDA traversal.
std::vector<Grid3d::Ixiyiz>
  Grid3d::get_hit_index_by_ray_tracing_algorithm(
    const Grid3d::VoxelTraversalParameters3d &prm) const
{
  std::vector<Ixiyiz> vec_ixiyiz;
  vec_ixiyiz.reserve(prm.manhattan_distance + 1);

  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  const int nbinz = get_nbinz();

  int ix = get_ix(prm.x_start);
  int iy = get_iy(prm.y_start);
  int iz = get_iz(prm.z_start);

  double t_max_x = prm.t_max_x0;
  double t_max_y = prm.t_max_y0;
  double t_max_z = prm.t_max_z0;

  vec_ixiyiz.push_back({ix,iy,iz});

  for(int t=0; t<=prm.manhattan_distance; t++){
    if( fabs(t_max_x) < fabs(t_max_y) ){
      if( fabs(t_max_x) < fabs(t_max_z) ){
        t_max_x += prm.t_delta_x;
        ix += (int)prm.vx_sign;
        if( ix < 0 || ix > nbinx-1 ) return vec_ixiyiz;
      } else {
        t_max_z += prm.t_delta_z;
        iz += (int)prm.vz_sign;
        if( iz < 0 || iz > nbinz-1 ) return vec_ixiyiz;
      }
    } else {
      if( fabs(t_max_y) < fabs(t_max_z) ){
        t_max_y += prm.t_delta_y;
        iy += (int)prm.vy_sign;
        if( iy < 0 || iy > nbiny-1 ) return vec_ixiyiz;
      } else {
        t_max_z += prm.t_delta_z;
        iz += (int)prm.vz_sign;
        if( iz < 0 || iz > nbinz-1 ) return vec_ixiyiz;
      }
    }
    vec_ixiyiz.push_back({ix,iy,iz});
  }
  return vec_ixiyiz;
}

// Optimized get_hit_voxel_index — uses cached_aabb3d_ instead of
// constructing a temporary AABB from 6 getters on every ray (opt #4/#5).
std::vector<Grid3d::Ixiyiz>
  Grid3d::get_hit_voxel_index( const Ray3d &ray3d, const double eps ) const
{
  std::vector<Ixiyiz> vec_ixiyiz;

  // use pre-built AABB — avoids 6 getter calls + 2 Vector3d + 1 AABB3d per ray
  auto [tf_hit_volume,tmin,tmax] = ray3d.is_intersect(cached_aabb3d_);

  if( tf_hit_volume == false ) return vec_ixiyiz;
  if( tmin < 0 && tmax < 0 ) return vec_ixiyiz;

  // delegate DDA setup to shared helper — eliminates duplicated ~80-line block
  const auto prm = build_traversal_params(ray3d, eps, tmin, tmax);

  return get_hit_index_by_ray_tracing_algorithm(prm);
}

// Optimized BL_max overload — performs a single ray-AABB intersection using
// cached_aabb3d_, then reuses the result for DDA setup. The legacy version
// performed the intersection twice: once for BL_max check, once inside the
// delegated get_hit_voxel_index(ray3d, eps) call (opt #4).
std::vector<Grid3d::Ixiyiz>
  Grid3d::get_hit_voxel_index( const Ray3d &ray3d, const double BL_max, const double eps ) const
{
  if( BL_max <= 0.0 ){
    return get_hit_voxel_index(ray3d, eps);
  }

  std::vector<Ixiyiz> vec_ixiyiz;

  // single intersection test using cached AABB — replaces two separate tests
  auto [tf_hit_volume,tmin,tmax] = ray3d.is_intersect(cached_aabb3d_);

  // BL_max warning check — done with the same intersection result
  if( tf_hit_volume && tmax > BL_max ){
    LOG_WARN_ND("Ray beam length ({:.1f}m) exceeds BL_max ({:.1f}m). Ray3d={}", tmax, BL_max, ray3d.to_string());
  }

  if( tf_hit_volume == false ) return vec_ixiyiz;
  if( tmin < 0 && tmax < 0 ) return vec_ixiyiz;

  // delegate DDA setup to shared helper — eliminates duplicated ~80-line block
  const auto prm = build_traversal_params(ray3d, eps, tmin, tmax);

  return get_hit_index_by_ray_tracing_algorithm(prm);
}

//-----------------------------------------------------------------
void Grid3d::out_header( FILE *fout )
{
  if (!fout) {
    throw std::runtime_error("Grid3d::out_header: fout is null.");
  }
  fprintf(fout,
    "# Grid3d name=%s nbinx=%d xmin=%E xmax=%E nbiny=%d ymin=%E ymax=%E nbinz=%d zmin=%E zmax=%E\n",
    name.c_str(),
    get_nbinx(), get_xmin(), get_xmax(),
    get_nbiny(), get_ymin(), get_ymax(),
    get_nbinz(), get_zmin(), get_zmax());
}

//-----------------------------------------------------------------
void Grid3d::out_grid_2d(const std::filesystem::path& pathout ) const
{
  FILE *fout = myapp::get_fout(pathout);
  for(int iy=0;iy<get_nbiny();iy++){
    for(int ix=0;ix<get_nbinx();ix++){
      double xmin = get_xlow(ix);
      double xmax = get_xup(ix);
      double ymin = get_ylow(iy);
      double ymax = get_yup(iy);
      fprintf(fout,"%d %d %E %E %E %E\n"
      ,ix,iy,xmin,xmax,ymin,ymax);
    }
  }
  myapp::close(fout,pathout);
}

void Grid3d::out_grid_3d(const std::filesystem::path& pathout) const
{
  FILE *fout = myapp::get_fout(pathout);
  int ix,iy,iz;
  double xmin,xmax,ymin,ymax,zmin,zmax;
  for(iz=0;iz<get_nbinz();iz++){
    for(iy=0;iy<get_nbiny();iy++){
      for(ix=0;ix<get_nbinx();ix++){
        xmin = get_xlow(ix);
        xmax = get_xup(ix);
        ymin = get_ylow(iy);
        ymax = get_yup(iy);
        zmin = get_zlow(iz);
        zmax = get_zup(iz);
        fprintf(fout,"%d %d %d %E %E %E %E %E %E\n"
        ,ix,iy,iz,xmin,xmax,ymin,ymax,zmin,zmax);
      }
    }
  }
  myapp::close(fout,pathout);
}

void Grid3d::display_voxel_index( FILE *fout,const double &x, const double &y, const double &z) const {
  const int ix = get_ix(x);
  const int iy = get_iy(y);
  const int iz = get_iz(z);
  const double x_lower = get_xlow(ix);
  const double y_lower = get_ylow(iy);
  const double z_lower = get_zlow(iz);
  const double x_upper = get_xup(ix);
  const double y_upper = get_yup(iy);
  const double z_upper = get_zup(iz);
  fprintf(fout,"%s : ix=%d, iy=%d, iz=%d | x:%.1lf-%.1lf y:%.1lf-%.1lf z:%.1lf-%.1lf \n"
  ,__FUNCTION__,ix,iy,iz,x_lower,x_upper,y_lower,y_upper,z_lower,z_upper);
}

void Grid3d::display_voxel_coordinate( FILE *fout, const int &ix, const int &iy, const int &iz) const {
  const double x_lower = get_xlow(ix);
  const double y_lower = get_ylow(iy);
  const double z_lower = get_zlow(iz);
  const double x_upper = get_xup(ix);
  const double y_upper = get_yup(iy);
  const double z_upper = get_zup(iz);
  fprintf(fout,"%s : ix=%d, iy=%d, iz=%d | x:%.1lf-%.1lf y:%.1lf-%.1lf z:%.1lf-%.1lf \n"
  ,__FUNCTION__,ix,iy,iz,x_lower,x_upper,y_lower,y_upper,z_lower,z_upper);
}

/// @brief calc x,y,z distance between two voxel
double Grid3d::calc_dist_dxyz(
  const Ixiyiz &ixiyiz0, const Ixiyiz &ixiyiz1 ) const
{  
  const auto [ix0,iy0,iz0] = ixiyiz0;
  const auto [ix1,iy1,iz1] = ixiyiz1;
  const int delta_ix = abs(ix1-ix0);
  const int delta_iy = abs(iy1-iy0);
  const int delta_iz = abs(iz1-iz0);
  const double delta_x = get_x_interval()*(double)delta_ix;
  const double delta_y = get_y_interval()*(double)delta_iy;
  const double delta_z = get_z_interval()*(double)delta_iz;
  return sqrt( delta_x*delta_x + delta_y*delta_y + delta_z*delta_z );
}

/// @brief calc x,y distance between two voxel
double Grid3d::calc_dist_dxy(
  const Ixiyiz &ixiyiz0, const Ixiyiz &ixiyiz1 ) const
{  
  const auto [ix0,iy0,iz0] = ixiyiz0;
  const auto [ix1,iy1,iz1] = ixiyiz1;
  const int delta_ix = abs(ix1-ix0);
  const int delta_iy = abs(iy1-iy0);
  const double delta_x = get_x_interval()*(double)delta_ix;
  const double delta_y = get_y_interval()*(double)delta_iy;
  return sqrt( delta_x*delta_x + delta_y*delta_y);
}

  /// @brief calc z distance between two voxel
double Grid3d::calc_dist_dz(
  const Ixiyiz &ixiyiz0, const Ixiyiz &ixiyiz1 ) const
{  
  const auto [ix0,iy0,iz0] = ixiyiz0;
  const auto [ix1,iy1,iz1] = ixiyiz1;
  const int delta_iz = abs(iz1-iz0);
  const double delta_z = get_z_interval()*(double)delta_iz;
  return delta_z;
}

//===============================================
// binary file I/O
//===============================================

/// @brief Save Grid3d to a binary stream.
void Grid3d::save( std::ofstream &ofs ) const
{
  io_binary::write_string(ofs,name);
  x_axis.save(ofs);
  y_axis.save(ofs);
  z_axis.save(ofs);
}

/// @brief Load Grid3d from a binary stream.
void Grid3d::load( std::ifstream &ifs, double tolerance_ratio )
{
  name = io_binary::read_string(ifs);
  x_axis.load(ifs, tolerance_ratio);
  y_axis.load(ifs, tolerance_ratio);
  z_axis.load(ifs, tolerance_ratio);
  rebuild_cached_aabb3d(); // keep cache consistent after loading axes
}
