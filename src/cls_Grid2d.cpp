// cls_Grid2d.cpp
#include "cls_Grid2d.hpp"
#include "ns_myapp.hpp"
#include "ns_mymacro.hpp"
#include "ns_mylogger.hpp"
#include "spdlog_pch.hpp"
#include "ns_mylogger.hpp"

#include "ns_type_definitions.hpp"
using namespace index_type_definitions;

#include <cmath>  // std::lround for robust bin-count rounding

//##############################################
//##############################################
// class Grid2d
//##############################################
//##############################################

/// @brief Inequality operator
bool Grid2d::operator!=(const Grid2d& other) const
{
  #ifdef NODEBUG
    if (x_axis != other.x_axis) return true;
    if (y_axis != other.y_axis) return true;
  #else
    if (x_axis != other.x_axis) { LOG_WARN("x_axis != other.x_axis"); return true; }
    if (y_axis != other.y_axis) { LOG_WARN("y_axis != other.y_axis"); return true; }
  #endif
  return false;
}

std::array<double,4> Grid2d::get_xy_lower_upper( const int ix, const int iy ) const
{
  double xmin = get_xlow(ix);
  double xmax = get_xup(ix);
  double ymin = get_ylow(iy);
  double ymax = get_yup(iy);
  return {xmin,xmax,ymin,ymax};
}

// output out_g2_header
void Grid2d::out_g2_header( FILE *fout ) const
{
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  const double xmin0 = get_xmin();
  const double xmax0 = get_xmax();
  const double ymin0 = get_ymin();
  const double ymax0 = get_ymax();
  const double xpit0 = get_x_interval();
  const double ypit0 = get_y_interval();
  fprintf(fout,"%d %.3lf %.3lf %.3lf\n",nbinx,xmin0,xmax0,xpit0);
  fprintf(fout,"%d %.3lf %.3lf %.3lf\n",nbiny,ymin0,ymax0,ypit0);
}

// Function to get std::vector< std::array<ix,iy> > from xmin xmax ymin ymax
std::vector< Ixiy > Grid2d::get_vec_ixiy(
  const double xmin_in, const double xmax_in
, const double ymin_in, const double ymax_in ) const
{
  std::vector<Ixiy> vec_ixiy;
  constexpr double small = 1000.0*std::numeric_limits<double>::epsilon();
 
  const int ix_min = get_ix( xmin_in + small );
  const int ix_max = get_ix( xmax_in - small );
  const int iy_min = get_iy( ymin_in + small );
  const int iy_max = get_iy( ymax_in - small );

  for(int iy=iy_min;iy<=iy_max;iy++){
    for(int ix=ix_min;ix<=ix_max;ix++){
      vec_ixiy.push_back( {ix,iy} );
    }
  }
  return vec_ixiy;
}

/// @brief get the index pair in the circle area
std::vector< Ixiy > Grid2d::get_vec_ixiy_in_circle(
  const double x_center, const double y_center, const double radius ) const
{
  std::vector<Grid2d::Ixiy> vec_ixiy;

  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  const double r2 = radius * radius;

  const int nxlen = static_cast<int>(2.0 * radius / get_x_interval()) + 2;
  const int nylen = static_cast<int>(2.0 * radius / get_y_interval()) + 2;

  // Pre-allocate reserve (optional)
  vec_ixiy.reserve(
    static_cast<std::size_t>(nxlen) * static_cast<std::size_t>(nylen));

  const double x_interval = get_x_interval();
  const double r_eps =  x_interval * std::numeric_limits<double>::epsilon();

  // 1) First create a bounding box containing the circle in physical coordinates
  double xmin = x_center - radius - r_eps;
  double xmax = x_center + radius + r_eps;
  double ymin = y_center - radius - r_eps;
  double ymax = y_center + radius + r_eps;

  // Clip to grid range
  xmin = std::max(xmin, get_xmin());
  xmax = std::min(xmax, get_xmax());
  ymin = std::max(ymin, get_ymin());
  ymax = std::min(ymax, get_ymax());

  // Convert y → iy, x → ix
  int ix_min = get_ix(xmin);
  int ix_max = get_ix(xmax);
  int iy_min = get_iy(ymin);
  int iy_max = get_iy(ymax);

  // ====== iy loop is outer from here ======
  for (int iy = iy_min; iy <= iy_max; ++iy) {
    const double y  = get_ycnt(iy);
    const double dy = y - y_center;

    // Scan only inner ix
    for (int ix = ix_min; ix <= ix_max; ++ix) {
      const double x  = get_xcnt(ix);
      const double dx = x - x_center;

      const double dist2 = dx * dx + dy * dy;

      if (dist2 <= r2) {  // Include boundary
        vec_ixiy.push_back(Ixiy{ix, iy});
      }
    }
  }

  return vec_ixiy;
}

/// @brief Convert a 2D Ixiy coordinate list to a 1D Ixiy coordinate list
///
/// @details
/// Flattens a 2D vector (outer index: y, inner index: x) of Ixiy coordinates
/// into a 1D vector in row-major order.
///
/// Ixiy is defined as:
/// @code
/// using Ixiy = std::array<int, 2>;
/// @endcode
///
/// @param[in] vec_vec_ixiy 2D Ixiy coordinate list
/// @return std::vector<Ixiy> 1D Ixiy coordinate list (row-major order)
///
/// @note
/// Example:
/// @code
/// std::vector<std::vector<Ixiy>> grid = {
///   { {1, 1}, {2, 1} },
///   { {1, 2}, {2, 2} }
/// };
/// auto flat = obj.get_vec_ixiy(grid);
/// // flat = { {1, 1}, {2, 1}, {1, 2}, {2, 2} }
/// @endcode
std::vector<Ixiy> Grid2d::get_vec_ixiy( 
    const std::vector<std::vector<Ixiy>> &vec_vec_ixiy ) const
{
  std::vector<Ixiy> vec_ixiy;
  for( const auto &vec_ixiy2 : vec_vec_ixiy ){
    for( const auto &ixiy : vec_ixiy2 ){
      vec_ixiy.push_back(ixiy);
    }
  }
  return vec_ixiy;
}

// set functions
// set x_axis and y_axis from std::vector\<Eigen::Vector3d\>
// assuming that data order is 
// x0 y0 z00         x0 y1 z01         x0 y2 z02
// x1 y0 z10 , ... , x1 y1 z11 , ... , x1 y2 z12
// x2 y0 z20         x2 y1 z21         x2 y2 z22
// and x0 < x1 < x2, y0 < y1 < y2
// 2023-04-05 14:48:16 update, add const bool tf_shift_x, tf_shift_y
// 2024-03-25 10:37:24 change return type void to bool
bool Grid2d::set_xy_axis_from_vec_xyz(
    const std::vector<std::array<double,3>> &vec_xyz
  , const bool tf_shift_x, const bool tf_shift_y
  , const double tolerance_ratio )
{
  constexpr double large_double = std::numeric_limits<double>::max();
  // check interval of x-grid and y-grid
  int i;
  double dx = -large_double;
  double dy = -large_double;
  double xmin_tmp =  large_double;
  double xmax_tmp = -large_double;
  double ymin_tmp =  large_double;
  double ymax_tmp = -large_double;

  // initial value is false
  // if first dx or dy (!= 0) come, they becomes true.
  bool flg_first_dx = false;
  bool flg_first_dy = false;
  double dx_first = -large_double;
  double dy_first = -large_double;

  // foreach loop for getting xmin, xmax, ymin, ymax
  for (const auto& [x, y, z] : vec_xyz) {
    if (x < xmin_tmp) xmin_tmp = x;
    if (y < ymin_tmp) ymin_tmp = y;
    if (x > xmax_tmp) xmax_tmp = x;
    if (y > ymax_tmp) ymax_tmp = y;
  }
  
  // loop for checking the interval
  for(i=0;i<vec_xyz.size()-1;i++){
    // Eigen::Vector3d &v3_0 = vec_v3.at(i+0); // ERROR
    const auto [x0,y0,z0] = vec_xyz.at(i+0);
    const auto [x1,y1,z1] = vec_xyz.at(i+1);
    dx = x1 - x0;
    dy = y1 - y0;
    // to get first dx
    if( flg_first_dx == false ){
      flg_first_dx = true;
      dx_first = fabs(dx);
      continue;
    }

    // if dx is not same, exit program
    // if( flg_first_dx == true && fabs(dx)!=dx_first && dy==0.0 ){
    if( flg_first_dx == true && fabs(fabs(dx)-fabs(dx_first)) > fabs(dx)*tolerance_ratio && dy==0.0 ){
      LOG_DEBUG("dx={:.15E}",dx);
      LOG_DEBUG("dx_first={:.15E}",dx_first);
      LOG_DEBUG("fabs(fabs(dx)-fabs(dx_first))={:.15E}",fabs(fabs(dx)-fabs(dx_first)));
      LOG_DEBUG("fabs(dx)*tolerance_ratio={:.15E}",fabs(dx)*tolerance_ratio);
      LOG_ERROR("dx != dx_first.");
      return false;
    }
    
    // to get the first dy
    if( flg_first_dy==false && dy!=0.0 ){
      flg_first_dy = true;
      dy_first = fabs(dy);
      continue;
    }

    // if dy interval is not same, exit program
    // if( flg_first_dy==true && dy!=0.0 && fabs(dy) != dy_first ){
    if( flg_first_dy==true && dy!=0.0 && fabs(fabs(dy)-fabs(dy_first)) > fabs(dy)*tolerance_ratio ){
      LOG_DEBUG("dy={:.15E}",dy);
      LOG_DEBUG("dy_first={:.15E}",dy_first);
      LOG_DEBUG("fabs(fabs(dy)-fabs(dy_first))={:.15E}",fabs(fabs(dy)-fabs(dy_first)));
      LOG_DEBUG("fabs(dy)*tolerance_ratio={:.15E}",fabs(dy)*tolerance_ratio);
      LOG_ERROR("dy != dy_first.");
      return false;
    }
  }
  
  // Round to nearest before +1: dx_first/dy_first carry tiny float excess
  // (e.g. 0.01 parsed as 0.010000000000000009), so a bare (int) cast truncates
  // an exact-multiple span one bin short. std::lround mirrors the eps guard in
  // Grid1d::get_index and keeps axis bin count consistent with the data.
  const int nbinx = (int)std::lround((xmax_tmp - xmin_tmp)/dx_first)+1;
  const int nbiny = (int)std::lround((ymax_tmp - ymin_tmp)/dy_first)+1;

  // assign of x_axis, y_axis
  double xmin_in, xmax_in, ymin_in, ymax_in;
  if( tf_shift_x==true ){
    xmin_in = xmin_tmp - 0.5*dx_first;
    xmax_in = xmax_tmp + 0.5*dx_first;
  }else{
    xmin_in = xmin_tmp;
    xmax_in = xmax_tmp + 1.0*dx_first;
  }
  if( tf_shift_y==true ){
    ymin_in = ymin_tmp - 0.5*dy_first;
    ymax_in = ymax_tmp + 0.5*dy_first;
  }else{
    ymin_in = ymin_tmp;
    ymax_in = ymax_tmp + 1.0*dy_first;
  }
  x_axis.assign("x_axis", nbinx, xmin_in, xmax_in, dx_first, tolerance_ratio);
  y_axis.assign("y_axis", nbiny, ymin_in, ymax_in, dy_first, tolerance_ratio);
  rebuild_cached_aabb2d();

  return true;
}

bool Grid2d::is_adjacent(
  const Ixiy& a, const Ixiy& b, const int distance_threshold) const
{
  const int dist = std::abs(a[0] - b[0]) + std::abs(a[1] - b[1]); 
  if( dist==0 ){
    LOG_DEBUG("dist=0");
    return false;
  }
  if( dist > distance_threshold ){
    return false;
  }
  return true;
}

void Grid2d::check_ix_inside(const int ix, const srcloc loc) const
{
  if (is_ix_inside(ix)) return;

  LOG_ERROR("ix is out of range.\n file = {}, function={}, line={}\n ix={}",
    loc.file_name(), loc.function_name(), loc.line(), ix);
  LOG_ERROR("ix = {}\n", ix);
  LOG_ERROR("nbinx = {}\n", x_axis.get_nbin());
  THROW_ERROR("Grid2d::check_ix_inside : ix is out of range {}", ix);
}

void Grid2d::check_iy_inside(const int iy, const srcloc loc) const
{
  if (is_iy_inside(iy)) return;

  LOG_ERROR("iy is out of range.\n file = {}, function={}, line={}\n iy={}",
    loc.file_name(), loc.function_name(), loc.line(), iy);
  LOG_ERROR("iy = {}\n", iy);
  LOG_ERROR("nbiny = {}\n", y_axis.get_nbin());
  THROW_ERROR("Grid2d::check_iy_inside : iy is out of range {}", iy);
}

// bool is inside Grid2d for (x,y)
bool Grid2d::is_inside( const double x, const double y ) const
{
  if( x < get_xmin() ) return false;
  if( x >= get_xmax() ) return false;
  if( y < get_ymin() ) return false;
  if( y >= get_ymax() ) return false;
  return true;
}

std::vector<std::vector<Ixiy>> 
  Grid2d::get_vec_vec_ixiy( const std::vector<Ixiy> &vec_ixiy  ) const
{
  std::vector<std::vector<Ixiy>> vec_vec_ixiy;
  const auto [ixmin, ixmax, iymin, iymax] = get_ixmin_ixmax_iymin_iymax(vec_ixiy);
  const int ixlen = ixmax - ixmin + 1;
  const int iylen = iymax - iymin + 1;

  // memory allocation
  vec_vec_ixiy.resize(iylen);
  for( int iy=0;iy<iylen;iy++ ){
    vec_vec_ixiy.at(iy).resize(ixlen);
  }

  // value assignment
  for(int iy=0;iy<iylen;iy++){
    for(int ix=0;ix<ixlen;ix++){
      vec_vec_ixiy.at(iy).at(ix) = {ix+ixmin,iy+iymin};
    }
  }

  return vec_vec_ixiy;
}

// get ixmin_ixmax_iymin_iymax from vec_ixiy
std::array<int,4> Grid2d::get_ixmin_ixmax_iymin_iymax(
  const std::vector<Ixiy> &vec_ixiy ) const
{
  constexpr int ibig = std::numeric_limits<int>::max();
  int ixmin =  ibig;
  int ixmax = -ibig; 
  int iymin =  ibig;
  int iymax = -ibig;
  for(const auto [ix,iy] : vec_ixiy){
    ixmin = std::min(ixmin,ix);
    ixmax = std::max(ixmax,ix);
    iymin = std::min(iymin,iy);
    iymax = std::max(iymax,iy);
  }
  return {ixmin,ixmax,iymin,iymax};
}

Ixiy Grid2d::get_ixiylen(
  const std::vector<Ixiy> &vec_ixiy ) const 
{
  const auto [ixmin,ixmax,iymin,iymax] = get_ixmin_ixmax_iymin_iymax(vec_ixiy);
  const int ixlen = ixmax-ixmin+1;
  const int iylen = iymax-iymin+1;
  return {ixlen,iylen};
}

Grid2d Grid2d::cut(
    const double x_lower, const double x_upper
  , const double y_lower, const double y_upper
  , const double x_eps, const double y_eps ) const
{
  if (x_eps < 0.0) THROW_ERROR("Grid2d::cut : x_eps < 0.0");
  if (y_eps < 0.0) THROW_ERROR("Grid2d::cut : y_eps < 0.0");

  Grid1d x_cut = x_axis.cut(x_lower, x_upper, x_eps);
  Grid1d y_cut = y_axis.cut(y_lower, y_upper, y_eps);

  Grid2d result(*this);

  const std::string x_name_org = x_axis.get_name();
  const std::string y_name_org = y_axis.get_name();

  result.set_x_axis(x_cut);
  result.set_y_axis(y_cut);
  result.set_x_axis_name(x_name_org);
  result.set_y_axis_name(y_name_org);

  return result;
}

// 2022-09-11 10:26:57
std::vector<Ixiy> Grid2d::get_hit_boxes_index_old(
  const double x0, const double y0, const double vx, const double vy) const {
  std::vector<Ixiy> vec_pair_index;
  vec_pair_index.clear();
  if(vx == 0.0){
    if( x0 <  get_xmin() ) return vec_pair_index;
    if( x0 >= get_xmax() ) return vec_pair_index;
    const int ix = get_ix(x0);
    for(int iy=0;iy<get_nbiny(); iy++){
      vec_pair_index.push_back( {ix, iy} );
    }
    return vec_pair_index;
  }
  else{
    // y = a*x + b
    const double a = vy/vx;
    const double b = -a*x0 + y0;
    int delta_iy=1; // if iy0 <= iy1;
    for(int ix=0;ix<get_nbinx();ix++){
      const double x0 = get_xlow(ix);
      const double x1 = get_xup(ix);
      const double y0 = a * x0 + b;
      const double y1 = a * x1 + b;
      int iy0 = get_iy(y0);
      int iy1 = get_iy(y1);

       // both are less than min
      if( iy0==-1 && iy1==-1 ) continue;

       // both are more than max
      if( is_iy_below_upper_limit(iy0) && is_iy_below_upper_limit(iy1) ) continue;

      // if not in the case above
      if( is_iy_above_lower_limit(iy0) ) iy0 = 0;
      if( is_iy_above_lower_limit(iy1) ) iy1 = 0;
      if( is_iy_below_upper_limit(iy0) ) iy0 = get_nbiny()-1;
      if( is_iy_below_upper_limit(iy1) ) iy1 = get_nbiny()-1;

      if( iy0 > iy1 ) delta_iy = -1;
      
      int iy=iy0;
      while( iy != iy1 ){
        vec_pair_index.push_back( {ix,iy} );
        iy += delta_iy;
      }
      vec_pair_index.push_back( {ix,iy} ); // case for iy=iy1
    }
    return vec_pair_index;
  }
}

/// @brief get ix, iy pair index array
/// @details with deistance in each box \n
///  https://gamedev.stackexchange.com/questions/81267/how-do-i-generalise-bresenhams-line-algorithm-to-floating-point-endpoints \n
///  @return std::vector<std::pair<ix,iy>> \n
// Build DDA traversal parameters from a valid ray-AABB intersection result.
// Extracts the common ~60 lines shared by get_hit_boxes_index() variants
// and traverse_ray_2d(), eliminating code duplication.
Grid2d::VoxelTraversalParameters2d Grid2d::build_traversal_params_2d(
    const Ray2d &ray2d, double eps, double tmin, double tmax) const
{
  VoxelTraversalParameters2d prm;

  const auto [inclination_factor_dx, inclination_factor_dy] = ray2d.get_path_factor();
  const auto [vx_sign, vy_sign] = ray2d.get_sign();
  prm.vx_sign = vx_sign;
  prm.vy_sign = vy_sign;

  prm.step_x = get_x_interval() * vx_sign;
  prm.step_y = get_y_interval() * vy_sign;

  const double t_eps = std::max(get_x_interval(), get_y_interval()) * eps;

  // factor (1.0 +- d_small) is for the case if x/y,1,2_area is on the edge of the boundary of the voxel
  const Eigen::Vector2d v2_pos1 = ray2d.pos() + (tmin + t_eps) * ray2d.dir();
  const Eigen::Vector2d v2_pos2 = ray2d.pos() + (tmax - t_eps) * ray2d.dir();

  if (tmin <= 0 && tmax >= 0) { // v2_pos0 is inside the box; x/y_start is v2_pos0
    prm.x_start = ray2d.x();
    prm.y_start = ray2d.y();
  }
  if (tmin > 0 && tmax > 0) { // v2_pos0 is outside the box; x/y_start is v2_pos1
    prm.x_start = v2_pos1.x();
    prm.y_start = v2_pos1.y();
  }

  prm.t_delta_x = prm.step_x * inclination_factor_dx;
  prm.t_delta_y = prm.step_y * inclination_factor_dy;

  const double x_end = v2_pos2.x();
  const double y_end = v2_pos2.y();

  int ix_start = get_ix(prm.x_start);
  int iy_start = get_iy(prm.y_start);
  int ix_end = get_ix(x_end);
  int iy_end = get_iy(y_end);

  // Start point should be inside the grid (AABB clipping ensures this)
  check_ix_inside(ix_start);
  check_iy_inside(iy_start);

  // End point might be on the boundary (OUT_OF_RANGE_UPPER/LOWER) due to floating-point precision
  // or ray exiting through DEM boundary. Clamp to valid range instead of throwing error.
  bool tf_hit_boundary = false;
  if (ix_end == Grid1d::OUT_OF_RANGE_UPPER) { ix_end = get_nbinx() - 1; tf_hit_boundary = true; }
  if (ix_end == Grid1d::OUT_OF_RANGE_LOWER) { ix_end = 0;               tf_hit_boundary = true; }
  if (iy_end == Grid1d::OUT_OF_RANGE_UPPER) { iy_end = get_nbiny() - 1; tf_hit_boundary = true; }
  if (iy_end == Grid1d::OUT_OF_RANGE_LOWER) { iy_end = 0;               tf_hit_boundary = true; }

  if (tf_hit_boundary) {
    LOG_WARN_ND("Ray hit DEM boundary. Consider using larger DEM. Ray2d={}", ray2d.to_string());
  }

  // x/y distance between current position to next x/y line
  double x_offset = 0.0;
  if (prm.step_x > 0) x_offset = get_xup(prm.x_start) - prm.x_start;
  if (prm.step_x < 0) x_offset = prm.x_start - get_xlow(prm.x_start);

  double y_offset = 0.0;
  if (prm.step_y > 0) y_offset = get_yup(prm.y_start) - prm.y_start;
  if (prm.step_y < 0) y_offset = prm.y_start - get_ylow(prm.y_start);

  prm.t_max_x0 = x_offset * inclination_factor_dx * vx_sign;
  prm.t_max_y0 = y_offset * inclination_factor_dy * vy_sign;
  prm.manhattan_distance = abs(ix_end - ix_start) + abs(iy_end - iy_start);

  return prm;
}

std::vector<Ixiy>
  Grid2d::get_hit_boxes_index( const Ray2d &ray2d, const double eps ) const
{
  auto [tf_hit_volume, tmin, tmax] = ray2d.is_intersect(cached_aabb2d_);
  if (!tf_hit_volume) return {};
  if (tmin < 0 && tmax < 0) return {};

  const auto prm = build_traversal_params_2d(ray2d, eps, tmin, tmax);
  return get_hit_index_by_ray_tracing_algorithm(prm);
}

/// @brief Get vector of hit box indices with beam length warning (optimization #4/#5)
std::vector<Ixiy>
  Grid2d::get_hit_boxes_index( const Ray2d &ray2d, const double BL_max, const double eps ) const
{
  // Single intersection test using cached AABB (optimizations #4 and #5)
  auto [tf_hit_volume, tmin, tmax] = ray2d.is_intersect(cached_aabb2d_);
  if (!tf_hit_volume) return {};
  if (tmin < 0 && tmax < 0) return {};

  // BL_max warning check using the same intersection result (no double intersection)
  if (BL_max > 0.0 && tmax > BL_max) {
    LOG_WARN_ND("Ray beam length ({:.1f}m) exceeds BL_max ({:.1f}m). Ray2d={}", tmax, BL_max, ray2d.to_string());
  }

  const auto prm = build_traversal_params_2d(ray2d, eps, tmin, tmax);
  return get_hit_index_by_ray_tracing_algorithm(prm);
}

/// @brief [Legacy] Trace a 2D ray — constructs a local AABB each call.
/// @deprecated Retained for regression testing; prefer get_hit_boxes_index.
std::vector<Ixiy>
  Grid2d::get_hit_boxes_index_legacy( const Ray2d &ray2d, const double eps ) const
{
  const Eigen::Vector2d v2_AABB_min(get_xmin(),get_ymin());
  const Eigen::Vector2d v2_AABB_max(get_xmax(),get_ymax());
  const AABB2d aabb2d(v2_AABB_min,v2_AABB_max);
  auto [tf_hit_volume,tmin,tmax] = ray2d.is_intersect(aabb2d);
  if( tf_hit_volume == false ) return {};
  if( tmin < 0 && tmax < 0 ) return {};

  const auto prm = build_traversal_params_2d(ray2d, eps, tmin, tmax);
  return get_hit_index_by_ray_tracing_algorithm(prm);
}

/// @brief [Legacy] Trace with BL_max — performs double ray-AABB intersection.
/// @deprecated Retained for regression testing; prefer get_hit_boxes_index.
std::vector<Ixiy>
  Grid2d::get_hit_boxes_index_legacy( const Ray2d &ray2d, const double BL_max, const double eps ) const
{
  if( BL_max <= 0.0 ){
    return get_hit_boxes_index_legacy(ray2d, eps);
  }

  const Eigen::Vector2d v2_AABB_min(get_xmin(),get_ymin());
  const Eigen::Vector2d v2_AABB_max(get_xmax(),get_ymax());
  const AABB2d aabb2d(v2_AABB_min,v2_AABB_max);
  auto [tf_hit_volume,tmin,tmax] = ray2d.is_intersect(aabb2d);

  if( tf_hit_volume && tmax > BL_max ){
    LOG_WARN_ND("Ray beam length ({:.1f}m) exceeds BL_max ({:.1f}m). Ray2d={}", tmax, BL_max, ray2d.to_string());
  }

  return get_hit_boxes_index_legacy(ray2d, eps);
}

// Get the std::vector<std::array<ix,iy>> in the original non-merged Grid2d corresponding to a cell in the merged Grid2d.
std::vector<Ixiy>
  Grid2d::get_non_merged_vec_ixiy( const Grid2d &g2d_merged, const int ix_merged, const int iy_merged ) const
{
  // Grid1d of *this
  Grid1d x_axis_org = this->get_x_axis();
  Grid1d y_axis_org = this->get_y_axis();

  // Grid1d after merged
  Grid1d x_axis_merged = g2d_merged.get_x_axis();
  Grid1d y_axis_merged = g2d_merged.get_y_axis();

  const auto [ixmin,ixmax] = x_axis_merged.get_original_index_min_max(x_axis_org, x_axis_merged, ix_merged);
  const auto [iymin,iymax] = y_axis_merged.get_original_index_min_max(y_axis_org, y_axis_merged, iy_merged);

  std::vector<Ixiy> vec_ixiy;
  for(int iy=iymin;iy<=iymax;iy++){
    for(int ix=ixmin;ix<=ixmax;ix++){
      vec_ixiy.push_back( {ix,iy} );
    }
  }
  if(vec_ixiy.size()==0) THROW_ERROR("vec_ixiy.size()==0");
  return vec_ixiy;
}

// return vector of the index of hit voxel 
std::vector<Ixiy> 
  Grid2d::get_hit_index_by_ray_tracing_algorithm(const Grid2d::VoxelTraversalParameters2d &prm) const
{
  std::vector<Ixiy> vec_int2;
  vec_int2.reserve(prm.manhattan_distance + 1);

  // get number of voxel in each direction
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();

  int ix = get_ix(prm.x_start);
  int iy = get_iy(prm.y_start);

  double t_max_x = prm.t_max_x0;
  double t_max_y = prm.t_max_y0;

  vec_int2.push_back( {ix,iy} );

  for(int t=0; t<=prm.manhattan_distance; t++){
    if( fabs(t_max_x) < fabs(t_max_y) ){
      t_max_x += prm.t_delta_x;
      // x += step_x;
      ix += (int)prm.vx_sign;
      // if out of volume, stop function
      if( ix < 0 || ix > nbinx-1 ) return vec_int2;
    } else { // |t_max_x| >= |t_max_y|
      t_max_y += prm.t_delta_y;
      // y += step_y;
      iy += (int)prm.vy_sign;
      // if out of volume, stop function
      if( iy < 0 || iy > nbiny-1 ) return vec_int2;
    }
    vec_int2.push_back( {ix,iy} );
  }
  return vec_int2;
}


// for bug check
void Grid2d::out_grid(const std::filesystem::path& pathout) const
{
  FILE *fout = myapp::get_fout(pathout);
  for(int iy=0;iy<get_nbiny();iy++){
    for(int ix=0;ix<get_nbinx();ix++){
      const double xlow = get_xlow(ix);
      const double xup  = get_xup(ix);
      const double ylow = get_ylow(iy);
      const double yup  = get_yup(iy);
      fprintf(fout,"%d %d %E %E %E %E\n"
      ,ix,iy,xlow,xup,ylow,yup);
    }
  }
  myapp::close(fout,pathout);
}

// Rebuild cached AABB2d from current axis min/max values.
// Called after any axis change so ray tracing always uses up-to-date bounds.
void Grid2d::rebuild_cached_aabb2d()
{
  cached_aabb2d_ = AABB2d(
    Eigen::Vector2d(x_axis.get_min(), y_axis.get_min()),
    Eigen::Vector2d(x_axis.get_max(), y_axis.get_max()));
}

//######################################################################
// Constructor
//######################################################################
Grid2d::Grid2d(const Grid1d &x_axis_in, const Grid1d &y_axis_in)
  : name(x_axis_in.get_name() + "_" + y_axis_in.get_name()),
    x_axis(x_axis_in),
    y_axis(y_axis_in)
{
  rebuild_cached_aabb2d();
}

//######################################################################
// Operators
//######################################################################
bool Grid2d::operator==(const Grid2d& other) const
{
  return !(*this != other);
}

//######################################################################
// Check functions
//######################################################################
bool Grid2d::is_inside(const int ix, const int iy) const
{
  return is_ix_inside(ix) && is_iy_inside(iy);
}

//######################################################################
// Getter functions
//######################################################################
double Grid2d::get_xcnt(const double value) const
{
  const int index = get_ix(value);
  check_ix_inside(index);
  return get_xcnt(index);
}

double Grid2d::get_ycnt(const double value) const
{
  const int index = get_iy(value);
  check_iy_inside(index);
  return get_ycnt(index);
}

double Grid2d::get_xlow(const double value) const
{
  const int index = get_ix(value);
  check_ix_inside(index);
  return get_xlow(index);
}

double Grid2d::get_ylow(const double value) const
{
  const int index = get_iy(value);
  check_iy_inside(index);
  return get_ylow(index);
}

double Grid2d::get_xup(const double value) const
{
  const int index = get_ix(value);
  check_ix_inside(index);
  return get_xup(index);
}

double Grid2d::get_yup(const double value) const
{
  const int index = get_iy(value);
  check_iy_inside(index);
  return get_yup(index);
}

double Grid2d::get_xlen(const int ix) const
{
  return get_xup(ix) - get_xlow(ix);
}

double Grid2d::get_ylen(const int iy) const
{
  return get_yup(iy) - get_ylow(iy);
}

std::array<double,2> Grid2d::get_xcnt_ycnt(const Ixiy ixiy)
{
  const auto [ix,iy] = ixiy;
  return {get_xcnt(ix), get_ycnt(iy)};
}

//######################################################################
// Setter functions
//######################################################################
void Grid2d::set_x_axis(const std::string x_name_in, const int nbinx_in,
  const double xmin_in, const double xmax_in, const double xpit_in)
{
  x_axis.assign(x_name_in, nbinx_in, xmin_in, xmax_in, xpit_in);
  rebuild_cached_aabb2d();
}

void Grid2d::set_y_axis(const std::string y_name_in, const int nbiny_in,
  const double ymin_in, const double ymax_in, const double ypit_in)
{
  y_axis.assign(y_name_in, nbiny_in, ymin_in, ymax_in, ypit_in);
  rebuild_cached_aabb2d();
}

//######################################################################
// Binary IO functions
//######################################################################
void Grid2d::save(std::ofstream &ofs) const
{
  io_binary::write_string(ofs, name);
  x_axis.save(ofs);
  y_axis.save(ofs);
}

void Grid2d::load(std::ifstream &ifs, double tolerance_ratio)
{
  name = io_binary::read_string(ifs);
  x_axis.load(ifs, tolerance_ratio);
  y_axis.load(ifs, tolerance_ratio);
  rebuild_cached_aabb2d();
}
