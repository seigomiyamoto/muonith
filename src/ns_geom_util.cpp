// src/ns_geom_util.cpp
#include "ns_geom_util.hpp"

#include <cmath>
#include <stdexcept>
#include "ns_mylogger.hpp"
#include "ns_mymacro.hpp"
#include "ns_angle_util.hpp"
#include "cls_Angle.hpp"

/// @brief Computes direction vector v3_dir=(dx,dy,0) from 2D positions v2_pos0 and v2_pos1.
/// @param v2_pos0 Start position vector (x0,y0)
/// @param v2_pos1 End position vector (x1,y1)
/// @return v3_dir Direction vector (dx,dy,0) as unit vector
Eigen::Vector3d geom_util::calc_v3_dir(
    const Eigen::Vector2d& v2_pos0, const Eigen::Vector2d& v2_pos1)
{
  Eigen::Vector2d v2_dir = v2_pos1 - v2_pos0;
  double norm = v2_dir.norm();
  if (norm < 1e-12){
    THROW_ERROR("geom_util::calc_v3_dir: zero-length direction vector");
  }
  v2_dir /= norm;
  return Eigen::Vector3d(v2_dir.x(), v2_dir.y(), 0.0);
}

/// @brief Computes rotation matrix R that rotates v3_dir_src to v3_dir_dst (z-dominant version: R = R2 * Rz).
/// @param v3_dir_src Source direction vector (unit vector)
/// @param v3_dir_dst Destination direction vector (unit vector)
/// @return Rotation matrix (3x3)
Eigen::Matrix3d geom_util::calc_rotmat_from_v3_to_v3_z_dominant(
  const Eigen::Vector3d& v3_dir_src, const Eigen::Vector3d& v3_dir_dst)
{
  const double eps = 1e-12;

  Eigen::Vector3d a = v3_dir_src.normalized();
  Eigen::Vector3d b = v3_dir_dst.normalized();

  // Projection onto XY plane
  Eigen::Vector2d a_xy(a.x(), a.y());
  Eigen::Vector2d b_xy(b.x(), b.y());
  const double na = a_xy.norm();
  const double nb = b_xy.norm();

  // Apply yaw (Z-rotation) first
  Eigen::Matrix3d Rz = Eigen::Matrix3d::Identity();
  if (na > 1e-15 && nb > 1e-15) {
    double yaw_a = std::atan2(a_xy.y(), a_xy.x());
    double yaw_b = std::atan2(b_xy.y(), b_xy.x());
    double dpsi  = yaw_b - yaw_a;
    // Normalization (optional)
    // while (dpsi >  M_PI) dpsi -= 2*M_PI;
    // while (dpsi <= -M_PI) dpsi += 2*M_PI;

    const double c = std::cos(dpsi), s = std::sin(dpsi);
    Rz <<  c, -s, 0,
           s,  c, 0,
           0,  0, 1;
  }

  // Align residual with shortest rotation
  Eigen::Vector3d a1 = Rz * a;
  double cos_theta = a1.dot(b);
  cos_theta = std::max(-1.0, std::min(1.0, cos_theta));
  if (cos_theta >= 1.0 - eps) {
    // Nearly aligned: Z-rotation alone suffices
    return Rz;
  }

  Eigen::Vector3d axis = a1.cross(b);
  double s = axis.norm();
  if (s < 1e-15) {
    // Antiparallel: choose axis deterministically
    Eigen::Vector3d e = (std::abs(a1.x()) < std::abs(a1.y()) && std::abs(a1.x()) < std::abs(a1.z()))
                          ? Eigen::Vector3d::UnitX()
                          : (std::abs(a1.y()) < std::abs(a1.z()) ? Eigen::Vector3d::UnitY()
                                                                 : Eigen::Vector3d::UnitZ());
    axis = a1.cross(e).normalized();
    return Eigen::AngleAxisd(M_PI, axis).toRotationMatrix() * Rz;
  }

  axis /= s;
  Angle ang = Angle(std::acos(cos_theta), Angle::Unit::Radian);
  Eigen::Matrix3d R2 = Eigen::AngleAxisd(ang.rad(), axis).toRotationMatrix();

  return R2 * Rz;  // Z-rotation first, then shortest rotation for residual
}

std::array<double,2> geom_util::world_to_local_2d(
    double dx, double dy, double cos_t, double sin_t)
{
  // Rz(-θ) = [ cosθ  sinθ; -sinθ  cosθ ]
  return { cos_t * dx + sin_t * dy,
          -sin_t * dx + cos_t * dy };
}

bool geom_util::intersect(const HorizontalEllipseDisk& disk, const Ray3d& ray)
{
  // Argument check
  if (disk.a <= 0.0 || disk.b <= 0.0){
    THROW_ERROR("geom_util::intersect: ellipse radii must be positive, now ({}, {}) [m]", disk.a, disk.b);
  }

  // Intersection with plane z=z_plane
  const double vz = ray.vz();
  if (std::abs(vz) < disk.eps) return false;  // Parallel
  const double t = (disk.z_plane - ray.z()) / vz;
  if (t < 0.0) return false;                  // Behind origin

  const double dx = ray.x() + ray.vx() * t - disk.xcnt;
  const double dy = ray.y() + ray.vy() * t - disk.ycnt;

  const double theta_rad = disk.yaw.rad();
  const double cos_t = std::cos(theta_rad);
  const double sin_t = std::sin(theta_rad);
  auto [dxx, dyy] = world_to_local_2d(dx, dy, cos_t, sin_t);

  const double val = (dxx*dxx)/(disk.a*disk.a) + (dyy*dyy)/(disk.b*disk.b);
  return val <= 1.0 + disk.eps;
}

bool geom_util::hits_cylinder_side(
  double z_base, double height, const Ray3d& ray, double t, double eps)
{
  if (height <= 0.0){
    THROW_ERROR("geom_util::hits_cylinder_side: height must be positive, now {} [m]", height);
  }
  if (eps <= 0.0){
    THROW_ERROR("geom_util::hits_cylinder_side: eps must be positive, now {}", eps);
  }
  if (t < 0.0) return false; // Behind half-line

  const double z = ray.z() + ray.vz() * t;
  return (z >= z_base - eps) && (z <= z_base + height + eps);
}

bool geom_util::intersect_side(const VerticalEllipticCylinderOpen& cyl, const Ray3d& ray)
{
  if (cyl.a <= 0.0 || cyl.b <= 0.0){
    THROW_ERROR("geom_util::intersect_side: ellipse radii must be positive, now ({}, {}) [m]", cyl.a, cyl.b);
  }
  if (cyl.height <= 0.0){
    THROW_ERROR("geom_util::intersect_side: height must be positive, now {} [m]", cyl.height);
  }
  const double theta = cyl.yaw.rad();
  const double cos_t = std::cos(theta);
  const double sin_t = std::sin(theta);

  // World to local coordinates
  const double dx = ray.x() - cyl.xcnt;
  const double dy = ray.y() - cyl.ycnt;
  auto [dxx, dyy] = world_to_local_2d(dx, dy, cos_t, sin_t);
  auto [vxx, vyy] = world_to_local_2d(ray.vx(), ray.vy(), cos_t, sin_t);

  // Quadratic equation satisfying (u/a)^2+(v/b)^2=1
  const double inva2 = 1.0/(cyl.a*cyl.a);
  const double invb2 = 1.0/(cyl.b*cyl.b);
  const double A = (vxx*vxx)*inva2 + (vyy*vyy)*invb2;
  const double B = 2.0*((dxx*vxx)*inva2 + (dyy*vyy)*invb2);
  const double C = (dxx*dxx)*inva2 + (dyy*dyy)*invb2 - 1.0;

  if (std::abs(A) < cyl.eps) {
    // Linear: B t + C = 0
    if (std::abs(B) < cyl.eps) {
      // C≈0: Ray origin is on side surface (tangent, generator, etc.)
      if (std::abs(C) <= cyl.eps)
        return hits_cylinder_side(cyl.z_base, cyl.height, ray, 0.0, cyl.eps);
      return false; // Otherwise no intersection
    }
    const double t_lin = -C / B;
    return hits_cylinder_side(cyl.z_base, cyl.height, ray, t_lin, cyl.eps);
  }
  double disc = B*B - 4*A*C;
  if (disc < -cyl.eps) return false;
  if (disc < 0.0) disc = 0.0;

  const double sqrt_disc = std::sqrt(disc);
  const double t1 = (-B - sqrt_disc) / (2*A);
  const double t2 = (-B + sqrt_disc) / (2*A);

  return (hits_cylinder_side(cyl.z_base, cyl.height, ray, t1, cyl.eps) ||
          hits_cylinder_side(cyl.z_base, cyl.height, ray, t2, cyl.eps));
}

bool geom_util::intersect_cap_plane(double z_plane
  , const VerticalEllipticCylinderCapped& cyl, const Ray3d& ray)
{
  if (cyl.a <= 0.0 || cyl.b <= 0.0){
    THROW_ERROR("geom_util::intersect_cap_plane: ellipse radii must be positive, now ({}, {}) [m]", cyl.a, cyl.b);
  }

  // Nearly parallel to plane: no intersection
  if (std::abs(ray.vz()) < cyl.eps) return false;

  // Intersection parameter t with plane (half-line, so t<0 is invalid)
  const double t = (z_plane - ray.z()) / ray.vz();
  if (t < 0.0) return false;

  // Convert intersection point (x,y) to offset from center
  const double x_hit = ray.x() + ray.vx() * t;
  const double y_hit = ray.y() + ray.vy() * t;
  const double dx = x_hit - cyl.xcnt;
  const double dy = y_hit - cyl.ycnt;

  // yaw is Angle class, right-hand rule with counterclockwise positive
  const double theta = cyl.yaw.rad();
  const double cos_t = std::cos(theta);
  const double sin_t = std::sin(theta);

  // Rotate world to local (ellipse principal axes)
  auto [dxx, dyy] = world_to_local_2d(dx, dy, cos_t, sin_t);

  // Check if inside ellipse
  const double a2 = cyl.a * cyl.a;
  const double b2 = cyl.b * cyl.b;
  const double val = (dxx*dxx)/a2 + (dyy*dyy)/b2;

  return val <= 1.0 + cyl.eps;
}

bool geom_util::intersect_caps(const VerticalEllipticCylinderCapped& cyl, const Ray3d& ray)
{
  if (cyl.a <= 0.0 || cyl.b <= 0.0){
    THROW_ERROR("geom_util::intersect_caps: ellipse radii must be positive, now ({}, {}) [m]", cyl.a, cyl.b);
  }
  if (cyl.height <= 0.0){
    THROW_ERROR("geom_util::intersect_caps: height must be positive, now {} [m]", cyl.height);
  }
  return intersect_cap_plane(cyl.z_base, cyl, ray)
      || intersect_cap_plane(cyl.z_base + cyl.height, cyl, ray);
}

bool geom_util::intersect(
  const VerticalEllipticCylinderCapped& cyl, const Ray3d& ray)
{
  return intersect_side(cyl, ray) || intersect_caps(cyl, ray);
}

/// @brief Anonymous namespace: local helper functions for computing path length inside cylinder.
/// @details
/// - Auxiliary functions called only by geom_util::path_length_inside().
/// - Finds ray-cylinder intersection parameters t (P(t)=P0+tV) and
///   enumerates entry/exit candidates.
/// - Placed in anonymous namespace to keep internal to this translation unit (.cpp), not exposed to public API.
namespace {

  /// @brief Checks if point is inside vertical elliptic cylinder at local coordinates (u,v,z).
  /// @param cyl Target elliptic cylinder
  /// @param ray Ray P(t)=P0+tV
  /// @param t Ray parameter
  /// @return true if within z range and (u/a)^2+(v/b)^2<=1
  ///
  /// - world→local rotation: reflects yaw only (in horizontal plane)
  /// - Also checks if z is in range [z_base, z_base+height]
  static bool _inside_vertical_elliptic_cylinder_open(
      const VerticalEllipticCylinderOpen& cyl, const Ray3d& ray, double t)
  {
    const double x = ray.x() + ray.vx() * t;
    const double y = ray.y() + ray.vy() * t;
    const double z = ray.z() + ray.vz() * t;

    // Height range check
    if (z < cyl.z_base - cyl.eps || z > cyl.z_base + cyl.height + cyl.eps)
      return false;

    // Convert offset from ellipse center to local (u,v) axes
    const double dx = x - cyl.xcnt;
    const double dy = y - cyl.ycnt;
    const double theta = cyl.yaw.rad();
    const double c = std::cos(theta), s = std::sin(theta);
    auto xxyy = geom_util::world_to_local_2d(dx, dy, c, s);

    const double xx = xxyy[0], yy = xxyy[1];
    const double val = (xx*xx)/(cyl.a*cyl.a) + (yy*yy)/(cyl.b*cyl.b);

    return val <= 1.0 + cyl.eps;
  }

  /// @brief Enumerates intersection parameters t from side surface quadratic equation.
  /// @param cyl Elliptic cylinder
  /// @param ray Ray
  /// @param ts Output: side surface intersection t candidates (stores only t>=0)
  ///
  /// - Equation: substituting ray into ((x/a)^2+(y/b)^2=1) → A t^2 + B t + C = 0
  /// - Pushes back real solutions t>=0 that are within height range.
  /// - If discriminant<0, no intersection; A≈0 is linear degenerate case.
  static void _collect_side_ts(
    const VerticalEllipticCylinderOpen& cyl, const Ray3d& ray, std::vector<double>& ts)
  {
    const double theta = cyl.yaw.rad();
    const double c = std::cos(theta), s = std::sin(theta);

    // World to local coordinate transformation (horizontal plane)
    const double dx = ray.x() - cyl.xcnt;
    const double dy = ray.y() - cyl.ycnt;
    auto dxxyy    = geom_util::world_to_local_2d(dx,     dy,     c, s);
    auto dir_xxyy = geom_util::world_to_local_2d(ray.vx(), ray.vy(), c, s);

    const double dxx = dxxyy[0], dyy = dxxyy[1];
    const double vxx = dir_xxyy[0], vyy = dir_xxyy[1];
    const double inva2 = 1.0/(cyl.a*cyl.a);
    const double invb2 = 1.0/(cyl.b*cyl.b);

    const double A = (vxx*vxx)*inva2 + (vyy*vyy)*invb2;
    const double B = 2.0*((dxx*vxx)*inva2 + (dyy*vyy)*invb2);
    const double C = (dxx*dxx)*inva2 + (dyy*dyy)*invb2 - 1.0;

    // Degenerate case (v nearly tangential, A≈0)
    if (std::abs(A) < cyl.eps) {
      if (std::abs(B) < cyl.eps) return; // No intersection or infinite contact (treat as length 0)
      const double t = -C / B;
      if (t >= 0.0 && geom_util::hits_cylinder_side(cyl.z_base, cyl.height, ray, t, cyl.eps))
        ts.push_back(t);
      return;
    }

    // Standard quadratic solution
    double disc = B*B - 4*A*C;
    if (disc < -cyl.eps) return;
    if (disc < 0.0) disc = 0.0;
    const double sd = std::sqrt(disc);

    const double t1 = (-B - sd) / (2*A);
    const double t2 = (-B + sd) / (2*A);

    // Add only valid range
    if (t1 >= 0.0 && geom_util::hits_cylinder_side(cyl.z_base, cyl.height, ray, t1, cyl.eps))
      ts.push_back(t1);
    if (t2 >= 0.0 && geom_util::hits_cylinder_side(cyl.z_base, cyl.height, ray, t2, cyl.eps))
      ts.push_back(t2);
  }

  /// @brief Evaluates ray intersection with cap plane at z=z_plane, adds t if inside ellipse.
  /// @param z_plane Cap plane z-coordinate
  /// @details
  /// - Ignores if parallel (Vz≈0).
  /// - Used to correctly handle passage through height range boundaries even for Open type.
  static void _collect_cap_t_if_inside(
      double z_plane, const VerticalEllipticCylinderOpen& cyl
    , const Ray3d& ray, std::vector<double>& ts)
  {
    if (std::abs(ray.vz()) < cyl.eps) return; // Parallel → no intersection
    const double t = (z_plane - ray.z()) / ray.vz();
    if (t < 0.0) return; // Behind half-line
    if (_inside_vertical_elliptic_cylinder_open(cyl, ray, t))
      ts.push_back(t);
  }

  /// @brief Sort and remove duplicate t values within tolerance.
  /// @param v List of t values
  /// @param eps Tolerance (e.g., 1e-12)
  static void _unique_sort(std::vector<double>& v, double eps)
  {
    std::sort(v.begin(), v.end());
    std::vector<double> out;
    for (double x : v) {
      if (out.empty() || std::abs(x - out.back()) > eps)
        out.push_back(x);
    }
    v.swap(out);
  }
}  // namespace anonymous

//======================================================================
// Main functions
//======================================================================

double geom_util::path_length_inside(const VerticalEllipticCylinderOpen& cyl, const Ray3d& ray)
{
  // ------------------ Input validation ------------------
  if (cyl.a <= 0.0 || cyl.b <= 0.0) {
    THROW_ERROR("geom_util::path_length_inside: ellipse radii must be positive, now ({}, {}) [m]", cyl.a, cyl.b);
  }
  if (cyl.height <= 0.0) {
    THROW_ERROR("geom_util::path_length_inside: height must be positive, now {} [m]", cyl.height);
  }

  // ------------------ Collect intersection candidates ------------------
  std::vector<double> ts;
  ts.reserve(4); // Roughly 2 from side + 2 from top/bottom

  _collect_side_ts(cyl, ray, ts); // Side surface intersections
  _collect_cap_t_if_inside(cyl.z_base, cyl, ray, ts); // Bottom cap z=z_base
  _collect_cap_t_if_inside(cyl.z_base + cyl.height, cyl, ray, ts); // Top cap z=z_base+height

  // Add t=0 if starting point is inside
  const bool inside0 = _inside_vertical_elliptic_cylinder_open(cyl, ray, 0.0);
  if (inside0) ts.push_back(0.0);

  if (ts.empty()) return 0.0; // No intersections

  _unique_sort(ts, std::max(1e-12, cyl.eps)); // Remove small duplicates

  // ------------------ Determine interval ------------------
  double t_enter = -1.0, t_exit = -1.0;

  if (inside0) {
    // Inside start → inside interval is up to first boundary
    t_enter = 0.0;
    for (double t : ts) {
      if (t > 0.0) { t_exit = t; break; }
    }
  } else {
    // Outside start → interval formed by first two intersections
    if (ts.size() < 2) return 0.0; // Entry/exit pair not complete
    int count = 0;
    for (double t : ts) {
      if (t >= 0.0) {
        if (count == 0) t_enter = t;
        else { t_exit = t; break; }
        ++count;
      }
    }
  }

  if (t_enter < 0.0 || t_exit < 0.0) return 0.0;

  // ------------------ Compute length ------------------
  const double vnorm = std::sqrt(ray.vx()*ray.vx() + ray.vy()*ray.vy() + ray.vz()*ray.vz());
  if (vnorm <= 0.0) {
    THROW_ERROR("geom_util::path_length_inside: ray direction must be non-zero");
  }

  const double len = (t_exit - t_enter) * vnorm;
  return (len > cyl.eps) ? len : 0.0;
}

double geom_util::path_length_inside(const VerticalEllipticCylinderCapped& cyl, const Ray3d& ray)
{
  // No copy/slice: safely upcast to base reference and forward
  const VerticalEllipticCylinderOpen& base =
    static_cast<const VerticalEllipticCylinderOpen&>(cyl);
  return path_length_inside(base, ray);
}

double geom_util::get_cylinder_diameter(const double L, const Angle dphi, const int n)
{
  if (L <= 0.0) THROW_ERROR("geom_util::get_cylinder_diameter: L must be positive, now {} [m]", L);
  if (dphi.rad() <= 0.0) THROW_ERROR("geom_util::get_cylinder_diameter: dphi must be positive, now {} [rad]", dphi.rad());
  if (n <= 0) THROW_ERROR("geom_util::get_cylinder_diameter: n must be positive, now {}", n);
  return L * std::tan(dphi.rad()) * static_cast<double>(n);
}

std::array<Angle,2> geom_util::calc_detector_delta_theta_range(
  const Ray3d& ray_det, const Eigen::Vector2d& v2_pos_obj, const double cyl_diameter)
{
  if (cyl_diameter <= 0.0){
    THROW_ERROR("cyl_diameter must be positive [m]");
  }
  const double r  = 0.5 * cyl_diameter;

  const double dx = v2_pos_obj.x() - ray_det.x();
  const double dy = v2_pos_obj.y() - ray_det.y();
  const double L2 = dx*dx + dy*dy;
  if (L2 <= r*r){
    LOG_ERROR("ray origin is inside the cylinder (L={}, r={}) [m]", std::sqrt(L2), r);
    THROW_ERROR("ray origin is inside the cylinder");
  }

  const double T = std::sqrt(L2 - r*r);
  const double alpha = std::atan2(r, T);          // Apparent half-angle
  const double theta_c   = std::atan2(dy, dx);    // Direction to center
  const double theta_det = std::atan2(ray_det.vy(), ray_det.vx()); // Detector front direction

  auto wrap_pi = [](double a){
    while (a <= -M_PI) a += 2*M_PI;
    while (a >   M_PI) a -= 2*M_PI;
    return a;
  };

  // Relative angle (lower and upper bounds with detector front as 0 rad)
  const double delta_min = wrap_pi((theta_c - alpha) - theta_det);
  const double delta_max = wrap_pi((theta_c + alpha) - theta_det);

  if (delta_min < delta_max)
    return { Angle(delta_min, Angle::Unit::Radian), Angle(delta_max, Angle::Unit::Radian) };
  else
    return { Angle(delta_max, Angle::Unit::Radian), Angle(delta_min, Angle::Unit::Radian) };
}

double geom_util::calc_cylinder_radius_from_horizontal_angle_size(
  const Ray3d& ray_det, const Eigen::Vector2d& v2_pos_obj, const Angle ang_tangent)
{
  const double dx = v2_pos_obj.x() - ray_det.x();
  const double dy = v2_pos_obj.y() - ray_det.y();
  const double L  = std::sqrt(dx*dx + dy*dy);
  return L * std::sin(ang_tangent.rad());
}

Angle geom_util::calc_horizontal_angle_from_cylinder_radius(
  const Ray3d& ray_det, const Eigen::Vector2d& v2_pos_obj, const double radius)
{
  if (radius <= 0.0){
    THROW_ERROR("radius must be positive [m]");
  }
  const double dx = v2_pos_obj.x() - ray_det.x();
  const double dy = v2_pos_obj.y() - ray_det.y();
  const double L  = std::sqrt(dx*dx + dy*dy);
  if (L <= radius){
    LOG_ERROR("ray origin is inside the cylinder (L={}, r={}) [m]", L, radius);
    THROW_ERROR("ray origin is inside the cylinder");
  }
  const double sin_alpha = radius / L;
  const double alpha_rad = std::asin(sin_alpha);
  return Angle(alpha_rad, Angle::Unit::Radian);
}


// Checks if a point is inside or outside an ellipse in 2D space
bool geom_util::isInsideEllipse(
    const Eigen::Vector2d& center, const Eigen::Vector2d& length
  , const Angle &alpha
  , const Eigen::Vector2d& pos)
{
  // Rotate the position vector by -alpha to align with the x-axis.
  // First rotate pos by -alpha around center to compute pos_rotated in new coordinate system.
  // This rotates the ellipse axes to be parallel to the x-axis.
  Eigen::Rotation2Dd rot_mat(-alpha.rad());
  Eigen::Vector2d pos_rotated = rot_mat * (pos - center);

  // Check if the rotated position is inside the unit circle.
  // Scale pos_rotated by ellipse axis lengths
  // to fit position into unit circle in new coordinate system.
  double x = pos_rotated(0) / length(0);
  double y = pos_rotated(1) / length(1);

  // Finally, check if inside unit circle and return boolean result.
  return x * x + y * y <= 1.0;
}

// Checks if a point is inside or outside an ellipsoid in 3D space
bool geom_util::isInsideEllipsoid(
    const Eigen::Vector3d& v3_center, const Eigen::Vector3d& v3_length
  , const Angle& theta_x, const Angle& theta_y, const Angle& theta_z
  , const angle_util::Rotation3dType &rotation_type
  , const Eigen::Vector3d& v3_pos)
{
  // Create rotation matrix
  const Eigen::Matrix3d rotationMatrix
   = angle_util::make_rotation_3d_matrix_ZYX(theta_x, theta_y, theta_z, rotation_type);

  // Apply inverse rotation to coordinates so principal axes align with coordinate axes
  // rotationMatrix.transpose() = rotationMatrix.inverse()
  // but transpose is computationally cheaper.
  const Eigen::Vector3d v3_pos_rotated = rotationMatrix.transpose() * (v3_pos - v3_center);

  // Scale coordinates by each ellipsoid axis
  const double x = v3_pos_rotated(0) / v3_length(0);
  const double y = v3_pos_rotated(1) / v3_length(1);
  const double z = v3_pos_rotated(2) / v3_length(2);

  // Check if point is inside ellipsoid
  return x * x + y * y + z * z <= 1.0;
}

// Checks if a point is inside or outside a cylinder in 3D space
// v3_length represents three lengths in x,y,z.
// v3_length(2) is length along cylinder axis, with fabs(z)<=0.5*v3_length(2)
// v3_length(0), v3_length(1) are lengths of two axes perpendicular to cylinder axis.
bool geom_util::isInsideCylinder(
    const Eigen::Vector3d& v3_center, const Eigen::Vector3d& v3_length
  , const Angle& theta_x, const Angle& theta_y, const Angle& theta_z
  , const angle_util::Rotation3dType &rotation_type
  , const Eigen::Vector3d& v3_pos)
{
  // Create rotation matrix
  const Eigen::Matrix3d rotationMatrix
   = angle_util::make_rotation_3d_matrix_ZYX(theta_x, theta_y, theta_z, rotation_type);

  // Transform point to elliptic cylinder local coordinate system
  // rotationMatrix.transpose() = rotationMatrix.inverse()
  // but transpose is computationally cheaper.
  const Eigen::Vector3d local_pos = rotationMatrix.transpose() * (v3_pos - v3_center);

  // Check if inside elliptic cylinder
  const double x = local_pos(0) / v3_length(0);
  const double y = local_pos(1) / v3_length(1);
  const double z = local_pos(2);

  if ((x * x + y * y <= 1.0) && (fabs(z) <= 0.5 * v3_length(2))) {
    return true;
  }
  return false;
}

// Checks if a point is inside or outside a cuboid (rectangular box) in 3D space.
// v3_length holds the full side lengths along x,y,z; half-extent is 0.5*v3_length on each axis.
bool geom_util::isInsideCuboid(
    const Eigen::Vector3d& v3_center, const Eigen::Vector3d& v3_length
  , const Angle& theta_x, const Angle& theta_y, const Angle& theta_z
  , const angle_util::Rotation3dType &rotation_type
  , const Eigen::Vector3d& v3_pos)
{
  // Create rotation matrix
  const Eigen::Matrix3d rotationMatrix
   = angle_util::make_rotation_3d_matrix_ZYX(theta_x, theta_y, theta_z, rotation_type);

  // Transform point to cuboid local coordinate system
  // rotationMatrix.transpose() = rotationMatrix.inverse()
  // but transpose is computationally cheaper.
  const Eigen::Vector3d local_pos = rotationMatrix.transpose() * (v3_pos - v3_center);

  // Inside if within half-extent on every axis
  return (fabs(local_pos(0)) <= 0.5 * v3_length(0))
      && (fabs(local_pos(1)) <= 0.5 * v3_length(1))
      && (fabs(local_pos(2)) <= 0.5 * v3_length(2));
}


// Determines which sign (positive or negative) a given 2D coordinate v2_pos falls into on a checkerboard.
// Pattern size is specified by v2_pattern_size.
// v2_start_pos is the checkerboard start point coordinates.
double geom_util::get_sign_uniform_checkerboard_2d(
    const Eigen::Vector2d &v2_start_pos
  , const Eigen::Vector2d &v2_pattern_size
  , const Eigen::Vector2d &v2_pos)
{
  // Offset coordinates by pattern start point
  const Eigen::Vector2d v2_offset_pos = v2_pos - v2_start_pos;

  // Normalize coordinates by pattern size
  const Eigen::Vector2d v2_normalized_pos = v2_offset_pos.cwiseQuotient(v2_pattern_size);

  // Determine which checkerboard cell the coordinates fall into
  const int ix = static_cast<int>( floor(v2_normalized_pos.x()) );
  const int iy = static_cast<int>( floor(v2_normalized_pos.y()) );

  // Determine if cell is positive or negative
  double sign_value = (ix % 2 == 0) ^ (iy % 2 == 0) ? -1.0 : 1.0;

  return sign_value;
}

// Determines which sign (positive or negative) a given 3D coordinate v3_pos falls into on a checkerboard.
// Pattern size is specified by v3_pattern_size.
// v3_start_pos is the checkerboard start point coordinates.
double geom_util::get_sign_uniform_checkerboard_3d(
    const Eigen::Vector3d &v3_start_pos,
    const Eigen::Vector3d &v3_pattern_size,
    const Eigen::Vector3d &v3_pos)
{
  // Offset coordinates by pattern start point
  const Eigen::Vector3d v3_offset_pos = v3_pos - v3_start_pos;

  // Normalize coordinates by pattern size
  const Eigen::Vector3d v3_normalized_pos = v3_offset_pos.cwiseQuotient(v3_pattern_size);

  // Determine which checkerboard cell the coordinates fall into
  const int ix = static_cast<int>( floor(v3_normalized_pos.x()) );
  const int iy = static_cast<int>( floor(v3_normalized_pos.y()) );
  const int iz = static_cast<int>( floor(v3_normalized_pos.z()) );

  // Determine if cell is positive or negative
  const double sign_value = (ix % 2 == 0) ^ (iy % 2 == 0) ^ (iz % 2 == 0) ? 1.0 : -1.0;

  return sign_value;
}
