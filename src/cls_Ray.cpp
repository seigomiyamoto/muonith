/// @file cls_Ray.cpp
/// @brief Implementation of Ray2d and Ray3d classes
#include "cls_Ray.hpp"
#include "ns_mymacro.hpp"
#include <limits>
#include <cmath>

//###########################################
// class Ray2d
//###########################################

// Compute perpendicular distance from a point to this ray
double Ray2d::distanceTo(const Eigen::Vector2d& point) const
{
  Eigen::Vector2d diff = point - v2_pos;
  Eigen::Vector2d projection = diff.dot(v2_dir) * v2_dir;
  Eigen::Vector2d rejection = diff - projection;
  return rejection.norm();
}

// Project a point onto this ray
Eigen::Vector2d Ray2d::projectOnto(const Eigen::Vector2d& point) const
{
  Eigen::Vector2d diff = point - v2_pos;
  double t = diff.dot(v2_dir);
  return v2_pos + t * v2_dir;
}

// Get position along ray at parameter t
Eigen::Vector2d Ray2d::getPointAtParameter(const double t) const
{
  return v2_pos + t * v2_dir;
}

// Validate that direction vector is unit length
void Ray2d::checkDirection(const double eps) const
{
  if (std::abs(v2_dir.norm() - 1.0) > eps) {
    THROW_ERROR("Ray2d::checkDirection: Direction vector must be a unit vector. norm={}",
                v2_dir.norm());
  }
}

// Compute path length factors for each axis
std::array<double,2> Ray2d::get_path_factor() const
{
  const double dydx = this->vy()/this->vx();
  const double dxdy = 1.0/dydx;
  const double path_factor_dx = std::sqrt(1.0 + dydx*dydx);
  const double path_factor_dy = path_factor_dx * std::fabs(dxdy);
  return {path_factor_dx, path_factor_dy};
}

// Get sign of direction components
std::array<double,2> Ray2d::get_sign() const
{
  const double vx_sign = (this->vx() > 0) - (this->vx() < 0);
  const double vy_sign = (this->vy() > 0) - (this->vy() < 0);
  return {vx_sign, vy_sign};
}

// Convert position to string representation
std::string Ray2d::pos_to_string(int width, int precision, char type) const {
  std::string fmt_str;
  if (width > 0) {
    fmt_str = fmt::format("[{{:{}.{}{}}}, {{:{}.{}{}}}]",
      width, precision, type, width, precision, type);
  } else {
    fmt_str = fmt::format("[{{:.{}{}}}, {{:.{}{}}}]",
      precision, type, precision, type);
  }
  return fmt::format(fmt::runtime(fmt_str), v2_pos.x(), v2_pos.y());
}

// Convert direction to string representation
std::string Ray2d::dir_to_string(int width, int precision, char type) const {
  std::string fmt_str;
  if (width > 0) {
    fmt_str = fmt::format("[{{:{}.{}{}}}, {{:{}.{}{}}}]",
      width, precision, type, width, precision, type);
  } else {
    fmt_str = fmt::format("[{{:.{}{}}}, {{:.{}{}}}]",
      precision, type, precision, type);
  }
  return fmt::format(fmt::runtime(fmt_str), v2_dir.x(), v2_dir.y());
}

// Convert ray to string representation for debugging
std::string Ray2d::to_string(int pos_width, int pos_precision, char pos_type,
                             int dir_width, int dir_precision, char dir_type) const {
  return fmt::format("Ray2d(pos={}, dir={})",
    pos_to_string(pos_width, pos_precision, pos_type),
    dir_to_string(dir_width, dir_precision, dir_type));
}

//==========================
// Ray2d-AABB intersection
//==========================
// Slab method implementation
// Reference: https://gdbooks.gitbooks.io/3dcollisions/content/Chapter3/raycast_aabb.html
bool Ray2d::is_intersect(const AABB2d &aabb2d,
  Eigen::Vector2d *v2_pos1, Eigen::Vector2d *v2_pos2) const
{
  using namespace std;
  double t1[2], t2[2];
  const double inf = std::numeric_limits<double>::infinity();
  constexpr double eps = 1.0e-12;
  for (int i = 0; i < 2; i++) {
    const double dir = this->dir(i);
    const double pos = this->pos(i);
    const double aabb_min = aabb2d.min(i);
    const double aabb_max = aabb2d.max(i);
    if (std::abs(dir) < eps) {
      // Ray is parallel to this slab
      if (pos < aabb_min || pos > aabb_max) {
        return false;  // Outside slab, no intersection
      }
      t1[i] = -inf;
      t2[i] =  inf;
    } else {
      t1[i] = (aabb_min - pos) / dir;
      t2[i] = (aabb_max - pos) / dir;
    }
  }
  double tmin = max(min(t1[0], t2[0]), min(t1[1], t2[1]));
  double tmax = min(max(t1[0], t2[0]), max(t1[1], t2[1]));

  // If tmax < 0, AABB is behind ray origin
  if (tmax < 0) return false;

  // If tmin >= tmax, ray misses the AABB
  if (tmin >= tmax) return false;

  *v2_pos1 = this->pos() + tmin * this->dir();
  *v2_pos2 = this->pos() + tmax * this->dir();

  return true;
}

// Slab method returning t-parameters
bool Ray2d::is_intersect(const AABB2d &aabb2d,
  double *tmin_ret, double *tmax_ret) const
{
  using namespace std;
  double t1[2], t2[2];
  const double inf = std::numeric_limits<double>::infinity();
  constexpr double eps = 1.0e-12;
  for (int i = 0; i < 2; i++) {
    const double dir = this->dir(i);
    const double pos = this->pos(i);
    const double aabb_min = aabb2d.min(i);
    const double aabb_max = aabb2d.max(i);
    if (std::abs(dir) < eps) {
      // Ray is parallel to this slab
      if (pos < aabb_min || pos > aabb_max) {
        return false;
      }
      t1[i] = -inf;
      t2[i] =  inf;
    } else {
      t1[i] = (aabb_min - pos) / dir;
      t2[i] = (aabb_max - pos) / dir;
    }
  }
  double tmin = max(min(t1[0], t2[0]), min(t1[1], t2[1]));
  double tmax = min(max(t1[0], t2[0]), max(t1[1], t2[1]));

  // If tmax < 0, AABB is behind ray origin
  if (tmax < 0) return false;

  // If tmin >= tmax, ray misses the AABB
  if (tmin >= tmax) return false;

  *tmin_ret = tmin;
  *tmax_ret = tmax;
  return true;
}

// Tuple return version (no special handling for parallel rays)
std::tuple<bool,double,double> Ray2d::is_intersect(const AABB2d &aabb2d) const
{
  double t1[2], t2[2];
  for (int i = 0; i < 2; i++) {
    t1[i] = (aabb2d.min(i) - this->pos(i)) / this->dir(i);
    t2[i] = (aabb2d.max(i) - this->pos(i)) / this->dir(i);
  }
  double tmin = std::max(std::min(t1[0], t2[0]), std::min(t1[1], t2[1]));
  double tmax = std::min(std::max(t1[0], t2[0]), std::max(t1[1], t2[1]));

  // If tmax < 0, AABB is behind ray origin
  if (tmax < 0) return std::make_tuple(false, tmin, tmax);

  // If tmin >= tmax, ray misses the AABB
  if (tmin >= tmax) return std::make_tuple(false, tmin, tmax);

  return std::make_tuple(true, tmin, tmax);
}

//###########################################
// class Ray3d
//###########################################

// Compute perpendicular distance from a point to this ray
double Ray3d::distanceTo(const Eigen::Vector3d& point) const
{
  Eigen::Vector3d diff = point - v3_pos;
  Eigen::Vector3d projection = diff.dot(v3_dir) * v3_dir;
  Eigen::Vector3d rejection = diff - projection;
  return rejection.norm();
}

// Project a point onto this ray
Eigen::Vector3d Ray3d::projectOnto(const Eigen::Vector3d& point) const
{
  Eigen::Vector3d diff = point - v3_pos;
  double t = diff.dot(v3_dir);
  return v3_pos + t * v3_dir;
}

// Get position along ray at parameter t
Eigen::Vector3d Ray3d::getPointAtParameter(const double t) const
{
  return v3_pos + t * v3_dir;
}

// Normalize direction vector or throw if invalid
void Ray3d::normalize_dir_or_throw(const double eps)
{
  const double n = v3_dir.norm();
  if (!std::isfinite(n) || n <= eps) {
    THROW_ERROR("Ray3d::normalize_dir_or_throw: Direction vector is zero or non-finite. norm={}",
                n);
  }
  v3_dir /= n;
}

//==========================
// Ray3d-AABB intersection
//==========================
// Slab method implementation
// Reference: https://gdbooks.gitbooks.io/3dcollisions/content/Chapter3/raycast_aabb.html
bool Ray3d::is_intersect(const AABB3d &aabb3d,
  Eigen::Vector3d *v3_pos1, Eigen::Vector3d *v3_pos2) const
{
  using namespace std;
  double t1[3], t2[3];
  const double inf = std::numeric_limits<double>::infinity();
  constexpr double eps = 1.0e-12;
  for (int i = 0; i < 3; i++) {
    const double dir = this->dir(i);
    const double pos = this->pos(i);
    const double aabb_min = aabb3d.min(i);
    const double aabb_max = aabb3d.max(i);
    // Handle near-zero direction to avoid 0/0 -> NaN
    if (std::abs(dir) < eps) {
      // Ray is parallel to this slab
      if (pos < aabb_min || pos > aabb_max) {
        return false;  // Outside slab, no intersection
      }
      t1[i] = -inf;
      t2[i] =  inf;
    } else {
      t1[i] = (aabb_min - pos) / dir;
      t2[i] = (aabb_max - pos) / dir;
    }
  }
  const double tmin = max(max(min(t1[0], t2[0]), min(t1[1], t2[1])), min(t1[2], t2[2]));
  const double tmax = min(min(max(t1[0], t2[0]), max(t1[1], t2[1])), max(t1[2], t2[2]));

  // If tmax < 0, AABB is behind ray origin
  if (tmax < 0) return false;

  // If tmin >= tmax, ray misses the AABB
  if (tmin >= tmax) return false;

  *v3_pos1 = this->pos() + tmin * this->dir();
  *v3_pos2 = this->pos() + tmax * this->dir();

  return true;
}

// Slab method returning t-parameters
bool Ray3d::is_intersect(const AABB3d &aabb3d,
  double *tmin_ret, double *tmax_ret) const
{
  using namespace std;
  double t1[3], t2[3];
  const double inf = std::numeric_limits<double>::infinity();
  constexpr double eps = 1.0e-12;
  for (int i = 0; i < 3; i++) {
    const double dir = this->dir(i);
    const double pos = this->pos(i);
    const double aabb_min = aabb3d.min(i);
    const double aabb_max = aabb3d.max(i);
    if (std::abs(dir) < eps) {
      // Ray is parallel to this slab
      if (pos < aabb_min || pos > aabb_max) {
        return false;
      }
      t1[i] = -inf;
      t2[i] =  inf;
    } else {
      t1[i] = (aabb_min - pos) / dir;
      t2[i] = (aabb_max - pos) / dir;
    }
  }
  // The largest entry t-value
  double tmin = max(
    max(min(t1[0], t2[0]), min(t1[1], t2[1])),
    min(t1[2], t2[2])
  );

  // The smallest exit t-value
  double tmax = min(
    min(max(t1[0], t2[0]), max(t1[1], t2[1])),
    max(t1[2], t2[2])
  );

  // Always set output parameters (useful for debugging even on miss)
  *tmin_ret = tmin;
  *tmax_ret = tmax;

  // If tmax < 0, AABB is behind ray origin
  if (tmax < 0) return false;

  // If tmin >= tmax, ray misses the AABB
  if (tmin >= tmax) return false;

  return true;
}

// Debug version with printf output for development
bool Ray3d::is_intersect_debug(const AABB3d &aabb3d,
  double *tmin_ret, double *tmax_ret) const
{
  using namespace std;
  double t1[3], t2[3];
  for (int i = 0; i < 3; i++) {
    t1[i] = (aabb3d.min(i) - this->pos(i)) / this->dir(i);
    t2[i] = (aabb3d.max(i) - this->pos(i)) / this->dir(i);
    std::fprintf(stderr, "debug_mode in %s : t1[%d] = %9.6lf, t2[%d] = %9.6lf\n",
                 __FUNCTION__, i, t1[i], i, t2[i]);
  }
  // The largest entry t-value
  double tmin = max(
    max(min(t1[0], t2[0]), min(t1[1], t2[1])),
    min(t1[2], t2[2])
  );

  // The smallest exit t-value
  double tmax = min(
    min(max(t1[0], t2[0]), max(t1[1], t2[1])),
    max(t1[2], t2[2])
  );
  std::fprintf(stderr, "debug_mode in %s : tmax-tmin = %E\n", __FUNCTION__, tmax - tmin);

  // If tmax < 0, AABB is behind ray origin
  if (tmax < 0) return false;

  // If tmin >= tmax, ray misses the AABB
  if (tmin >= tmax) return false;

  std::fprintf(stderr, "### intersection: tmax-tmin = %E\n", tmax - tmin);

  return true;
}

// Tuple return version with parallel-ray handling
std::tuple<bool,double,double> Ray3d::is_intersect(const AABB3d &aabb3d) const
{
  double t1[3], t2[3];
  const double inf = std::numeric_limits<double>::infinity();
  constexpr double eps = 1.0e-12;
  for (int i = 0; i < 3; i++) {
    const double dir = this->dir(i);
    const double pos = this->pos(i);
    const double aabb_min = aabb3d.min(i);
    const double aabb_max = aabb3d.max(i);
    // Handle near-zero direction to avoid 0/0 -> NaN
    if (std::abs(dir) < eps) {
      // Ray is parallel to this slab
      if (pos < aabb_min || pos > aabb_max) {
        return std::make_tuple(false, 0.0, 0.0);  // Outside slab, no intersection
      }
      t1[i] = -inf;
      t2[i] =  inf;
    } else {
      t1[i] = (aabb_min - pos) / dir;
      t2[i] = (aabb_max - pos) / dir;
    }
  }
  // The largest entry t-value
  double tmin = std::max(
    std::max(std::min(t1[0], t2[0]), std::min(t1[1], t2[1])),
    std::min(t1[2], t2[2])
  );

  // The smallest exit t-value
  double tmax = std::min(
    std::min(std::max(t1[0], t2[0]), std::max(t1[1], t2[1])),
    std::max(t1[2], t2[2])
  );

  // If tmax < 0, AABB is behind ray origin
  if (tmax < 0) return std::make_tuple(false, tmin, tmax);

  // If tmin >= tmax, ray misses the AABB
  if (tmin >= tmax) return std::make_tuple(false, tmin, tmax);

  // Origin inside AABB if tmin < 0
  return std::make_tuple(true, tmin, tmax);
}

// Compute path length factors for each axis
std::array<double,3> Ray3d::get_path_factor() const
{
  const double dydx = this->vy() / this->vx();
  const double dxdy = 1.0 / dydx;
  const double dzdx = this->vz() / this->vx();
  const double dxdz = 1.0 / dzdx;
  const double path_factor_dx = std::sqrt(1.0 + dydx * dydx);
  const double path_factor_dy = path_factor_dx * std::fabs(dxdy);
  const double path_factor_dz = path_factor_dx * std::fabs(dxdz);
  return {path_factor_dx, path_factor_dy, path_factor_dz};
}

// Get sign of direction components
std::array<double,3> Ray3d::get_sign() const
{
  const double vx_sign = (this->vx() > 0) - (this->vx() < 0);
  const double vy_sign = (this->vy() > 0) - (this->vy() < 0);
  const double vz_sign = (this->vz() > 0) - (this->vz() < 0);
  return {vx_sign, vy_sign, vz_sign};
}

// Convert position to string representation
std::string Ray3d::pos_to_string(int width, int precision, char type) const {
  std::string fmt_str;
  if (width > 0) {
    fmt_str = fmt::format("[{{:{}.{}{}}}, {{:{}.{}{}}}, {{:{}.{}{}}}]",
      width, precision, type, width, precision, type, width, precision, type);
  } else {
    fmt_str = fmt::format("[{{:.{}{}}}, {{:.{}{}}}, {{:.{}{}}}]",
      precision, type, precision, type, precision, type);
  }
  return fmt::format(fmt::runtime(fmt_str), v3_pos.x(), v3_pos.y(), v3_pos.z());
}

// Convert direction to string representation
std::string Ray3d::dir_to_string(int width, int precision, char type) const {
  std::string fmt_str;
  if (width > 0) {
    fmt_str = fmt::format("[{{:{}.{}{}}}, {{:{}.{}{}}}, {{:{}.{}{}}}]",
      width, precision, type, width, precision, type, width, precision, type);
  } else {
    fmt_str = fmt::format("[{{:.{}{}}}, {{:.{}{}}}, {{:.{}{}}}]",
      precision, type, precision, type, precision, type);
  }
  return fmt::format(fmt::runtime(fmt_str), v3_dir.x(), v3_dir.y(), v3_dir.z());
}

// Convert ray to string representation for debugging
std::string Ray3d::to_string(int pos_width, int pos_precision, char pos_type,
                             int dir_width, int dir_precision, char dir_type) const {
  return fmt::format("Ray3d(pos={}, dir={})",
    pos_to_string(pos_width, pos_precision, pos_type),
    dir_to_string(dir_width, dir_precision, dir_type));
}

//==========================
// Binary I/O functions
//==========================

// Write Eigen::Vector3d to binary stream
void Ray3d::write_vec3d(std::ofstream &ofs, const Eigen::Vector3d &vec3d) const
{
  ofs.write(reinterpret_cast<const char*>(&vec3d(0)), sizeof(double));
  ofs.write(reinterpret_cast<const char*>(&vec3d(1)), sizeof(double));
  ofs.write(reinterpret_cast<const char*>(&vec3d(2)), sizeof(double));
}

// Read Eigen::Vector3d from binary stream
Eigen::Vector3d Ray3d::read_vec3d(std::ifstream &ifs)
{
  Eigen::Vector3d vec3d;
  ifs.read(reinterpret_cast<char*>(&vec3d(0)), sizeof(double));
  ifs.read(reinterpret_cast<char*>(&vec3d(1)), sizeof(double));
  ifs.read(reinterpret_cast<char*>(&vec3d(2)), sizeof(double));
  return vec3d;
}

// Save ray to binary stream
void Ray3d::save(std::ofstream &ofs) const
{
  write_vec3d(ofs, v3_pos);
  write_vec3d(ofs, v3_dir);
}

// Load ray from binary stream
void Ray3d::load(std::ifstream &ifs)
{
  v3_pos = read_vec3d(ifs);
  v3_dir = read_vec3d(ifs);
  normalize_dir_or_throw();
}
