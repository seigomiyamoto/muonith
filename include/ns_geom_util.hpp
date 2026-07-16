/// @file ns_geom_util.hpp
/// @brief Geometric utility functions
/// @details Namespace providing geometric calculations for intersections, distances, and transformations.
#pragma once
#include <array>
#include <Eigen/Dense>
#include "cls_Ray.hpp"
#include "cls_Angle.hpp"
#include "st_shapes.hpp"
#include "ns_angle_util.hpp"

//#########################################################################
//#########################################################################
/// @namespace geom_util
/// @brief Geometric utility functions
//#########################################################################
//#########################################################################
namespace geom_util {

  /// @brief Computes direction vector v3_dir=(dx,dy,0) from 2D positions v2_pos0 and v2_pos1.
  /// @param v2_pos0 Start position vector (x0,y0)
  /// @param v2_pos1 End position vector (x1,y1)
  /// @return v3_dir Direction vector (dx,dy,0) as unit vector
  Eigen::Vector3d calc_v3_dir(
      const Eigen::Vector2d& v2_pos0, const Eigen::Vector2d& v2_pos1);

  /// @brief Computes rotation matrix R that rotates v3_dir_src to v3_dir_dst (z-dominant version: R = R2 * Rz).
  /// @param v3_dir_src Source direction vector (unit vector)
  /// @param v3_dir_dst Destination direction vector (unit vector)
  /// @return Rotation matrix (3x3)
  Eigen::Matrix3d calc_rotmat_from_v3_to_v3_z_dominant(
    const Eigen::Vector3d& v3_dir_src, const Eigen::Vector3d& v3_dir_dst);

  /// @brief Transforms world coordinates (dx,dy) to local coordinates (dxx,dyy) using rotation angle (cos_t, sin_t).
  /// @param dx x - xcnt
  /// @param dy y - ycnt
  /// @param cos_t cos(theta)
  /// @param sin_t sin(theta)
  /// @returns (dxx,dyy) Coordinates in local frame
  /// @details
  /// Local frame is the world frame rotated by theta around the z-axis.
  /// Rotation matrix Rz(θ) follows right-hand rule (counterclockwise positive):
  std::array<double,2> world_to_local_2d(
      double dx, double dy, double cos_t, double sin_t);

  /// @brief Checks intersection with a horizontal ellipse disk.
  /// @param disk Horizontal ellipse disk parameters
  /// @param ray Ray
  /// @return true if ray intersects inside the disk
  bool intersect(const HorizontalEllipseDisk& disk, const Ray3d& ray);

  /// @brief Checks if z(t) of ray at parameter t lies within cylinder height interval [z_base, z_base+height].
  /// @param z_base Cylinder base z-coordinate [m]
  /// @param height Cylinder height [m]
  /// @param ray Ray (origin P0 and direction V)
  /// @param t Parameter of ray equation P(t)=P0+tV (typically from side surface quadratic solution)
  /// @param eps Boundary tolerance (default 1e-12)
  /// @return true if within interval (t>=0 is also required)
  /// @details
  ///   Returns true if z(t) = ray.z() + ray.vz() * t is in [z_base, z_base+height].
  ///   Returns false if t<0 (behind the half-line origin).
  bool hits_cylinder_side(double z_base, double height
        , const Ray3d& ray, double t, double eps = 1e-12);

  /// @brief Checks if ray intersects the side surface of a vertical elliptic cylinder (open, no caps).
  /// @param cyl Vertical elliptic cylinder parameters (open)
  /// @param ray Ray (P(t)=P0+tV, t>=0)
  /// @return true if ray intersects the side surface
  /// @details
  ///   - Transforms world coordinates to cylinder local coordinates (u,v,z), solves quadratic equation,
  ///     and returns true if real solution t>=0 and z_base<=z<=z_base+height.
  bool intersect_side(const VerticalEllipticCylinderOpen& cyl, const Ray3d& ray);

  /// @brief Checks if ray intersects horizontal ellipse cap plane at z=z_plane (uses cylinder yaw).
  /// @param z_plane Plane z-coordinate [m]
  /// @param cyl Vertical elliptic cylinder parameters (capped; uses center (xc,yc), radii a,b, yaw, eps)
  /// @param ray Ray (P(t)=P0+tV, t>=0)
  /// @return true if ray intersects inside ellipse on cap plane, false otherwise
  /// @details
  /// - Ray equation: \f$ \mathbf{P}(t) = \mathbf{P}_0 + t\,\mathbf{V}, \; t \ge 0 \f$
  /// - Cap plane: \f$ z = z_{\rm plane} \f$
  /// - Intersection parameter: \f$ t = (z_{\rm plane}-P_{0z})/V_z \f$
  /// - Intersection point \f$ (dx,dy) \f$ is rotated to ellipse local coordinates \f$ (dxx,dyy) \f$ (right-hand rule, yaw>0 counterclockwise)
  ///   \f[
  ///     \begin{pmatrix}dxx\\dyy\end{pmatrix}
  ///     =
  ///     R_z(-\theta)\,
  ///     \begin{pmatrix}x-x_c\\y-y_c\end{pmatrix},\quad
  ///     R_z(-\theta)=\begin{pmatrix}\cos\theta&\sin\theta\\-\sin\theta&\cos\theta\end{pmatrix}
  ///   \f]
  /// - Inside ellipse test: \f$ (dxx/a)^2 + (dyy/b)^2 \le 1 \f$
  /// - Numerical stability: eps used for semi-open boundary tolerance
  bool intersect_cap_plane(double z_plane
  , const VerticalEllipticCylinderCapped& cyl, const Ray3d& ray);

  /// @brief Checks if ray intersects top or bottom caps (z=z_base or z=z_base+height) of vertical elliptic cylinder.
  /// @param cyl Vertical elliptic cylinder parameters (capped)
  /// @param ray Ray (P(t)=P0+tV, t>=0)
  /// @return true if ray intersects either cap plane
  /// @details
  ///   - Finds intersection with plane z=z_plane and checks if inside ellipse: (dxx/a)^2+(dyy/b)^2<=1.
  bool intersect_caps(const VerticalEllipticCylinderCapped& cyl, const Ray3d& ray);

  /// @brief Checks overall intersection with vertical elliptic cylinder (capped, with top and bottom).
  /// @param cyl Vertical elliptic cylinder parameters (capped)
  /// @param ray Ray
  /// @return true if ray intersects side surface or caps
  bool intersect(const VerticalEllipticCylinderCapped& cyl, const Ray3d& ray);

  /// @brief Computes total distance traveled by ray inside vertical elliptic cylinder.
  /// @param cyl Vertical elliptic cylinder parameters
  /// @param ray Ray (P(t)=P0+tV, t>0)
  /// @return Path length inside [m]
  ///
  /// Algorithm overview:
  ///  1. Collect all intersection candidate t values with side surface and top/bottom caps
  ///  2. Extract only t>=0 range and sort in ascending order
  ///  3. If starting point is inside, treat [0, first boundary] as inside interval;
  ///     if outside, treat [first intersection, next intersection] as inside interval
  ///  4. Return interval length (t_exit - t_enter)*|V|
  ///
  /// Notes:
  ///  - Considers case where ray starts from inside
  ///  - t and V can be non-normalized (actual length corrected by V norm)
  ///  - Tangent and backward intersections are automatically excluded
  double path_length_inside(const VerticalEllipticCylinderOpen& cyl, const Ray3d& ray);

  /// @brief (Overload) Same logic for capped cylinder.
  double path_length_inside(const VerticalEllipticCylinderCapped& cyl, const Ray3d& ray);

  /// @brief Gets cylinder diameter.
  /// @param L Distance from detector to cylinder center
  /// @param dphi Detector minimum angular element
  double get_cylinder_diameter(const double L, const Angle dphi, const int n);

  /// @brief Computes relative horizontal angle range (lower and upper bounds) from detector position/direction vector and cylinder horizontal center/diameter, \n
  ///        treating detector orientation as 0 rad.
  /// @returns std::array<Angle,2> {delta_min, delta_max}  [rad]
  /// @param ray_det Detector position and direction vector
  /// @param v2_pos_obj Cylinder horizontal center (x,y)
  /// @param cyl_diameter Cylinder diameter
  std::array<Angle,2> calc_detector_delta_theta_range(
    const Ray3d& ray_det, const Eigen::Vector2d& v2_pos_obj, const double cyl_diameter);

  /// @brief Computes circle radius such that the angle extension lines are tangent to the circle.
  /// @param ray_det Detector position and direction vector
  /// @param v2_pos_obj Circle center position (x,y)
  /// @param ang_tangent Horizontal angular size (half-angle) [rad]
  /// @return Circle radius [m]
  double calc_cylinder_radius_from_horizontal_angle_size(
    const Ray3d& ray_det, const Eigen::Vector2d& v2_pos_obj, const Angle ang_tangent);

  /// @brief Computes the angle between two lines emanating from the same point that are tangent to a circle of radius.
  /// @param ray_det Detector position and direction vector
  /// @param v2_pos_obj Circle center position (x,y)
  /// @param radius Circle radius [m]
  /// @return Angle between the two lines (half-angle)
  Angle calc_horizontal_angle_from_cylinder_radius(
    const Ray3d& ray_det, const Eigen::Vector2d& v2_pos_obj, const double radius);

  //=================================================================
  /// @name structure_myapp
  /// @details Geometric shape functions
  /// @{

  /// @brief Checks if a point is inside or outside an ellipse in 2D space.
  /// @param[in] center Ellipse center coordinates
  /// @param[in] length Lengths in x and y directions
  /// @param[in] alpha Ellipse rotation angle
  /// @param[in] pos Point coordinates to test
  /// @return true: inside ellipse, false: outside ellipse
  bool isInsideEllipse(
      const Eigen::Vector2d& center, const Eigen::Vector2d& length
    , const Angle &alpha
    , const Eigen::Vector2d& pos);

  /// @brief Checks if a point is inside or outside an ellipsoid in 3D space.
  /// @details Rotation axis order is z-axis, y-axis, x-axis
  /// @param[in] v3_center Ellipsoid center coordinates (x,y,z)
  /// @param[in] v3_length Lengths in x, y, z directions
  /// @param[in] theta_x Rotation angle around x-axis (counterclockwise)
  /// @param[in] theta_y Rotation angle around y-axis (counterclockwise)
  /// @param[in] theta_z Rotation angle around z-axis (counterclockwise)
  /// @param[in] rotation_type LOCAL or GLOBAL
  /// @param[in] v3_pos Point coordinates to test
  /// @return true: inside ellipsoid, false: outside ellipsoid
  bool isInsideEllipsoid(
          const Eigen::Vector3d& v3_center, const Eigen::Vector3d& v3_length
        , const Angle& theta_x, const Angle& theta_y, const Angle& theta_z
        , const angle_util::Rotation3dType &rotation_type
        , const Eigen::Vector3d& v3_pos);

  /// @brief Checks if a point is inside or outside a cylinder in 3D space.
  /// @details Rotation axis order is z-axis, y-axis, x-axis
  /// @param[in] v3_center Cylinder center coordinates (x,y,z)
  /// @param[in] v3_length v3_length(0), v3_length(1) are lengths of two axes perpendicular to cylinder axis, \n
  /// v3_length(2) is length along cylinder axis, with fabs(z)<=0.5*v3_length(2)
  /// @param[in] theta_x Rotation angle around x-axis (counterclockwise)
  /// @param[in] theta_y Rotation angle around y-axis (counterclockwise)
  /// @param[in] theta_z Rotation angle around z-axis (counterclockwise)
  /// @param[in] rotation_type LOCAL or GLOBAL
  /// @param[in] v3_pos Point coordinates to test
  /// @return true: inside cylinder, false: outside cylinder
  bool isInsideCylinder(
          const Eigen::Vector3d& v3_center, const Eigen::Vector3d& v3_length
        , const Angle& theta_x, const Angle& theta_y, const Angle& theta_z
        , const angle_util::Rotation3dType &rotation_type
        , const Eigen::Vector3d& v3_pos);

  /// @brief Checks if a point is inside or outside a cuboid (rectangular box) in 3D space.
  /// @details Rotation axis order is z-axis, y-axis, x-axis.
  /// @param[in] v3_center Cuboid center coordinates (x,y,z)
  /// @param[in] v3_length Full side lengths along x, y, z; half-extent is 0.5*v3_length on each axis
  /// @param[in] theta_x Rotation angle around x-axis (counterclockwise)
  /// @param[in] theta_y Rotation angle around y-axis (counterclockwise)
  /// @param[in] theta_z Rotation angle around z-axis (counterclockwise)
  /// @param[in] rotation_type LOCAL or GLOBAL
  /// @param[in] v3_pos Point coordinates to test
  /// @return true: inside cuboid, false: outside cuboid
  bool isInsideCuboid(
          const Eigen::Vector3d& v3_center, const Eigen::Vector3d& v3_length
        , const Angle& theta_x, const Angle& theta_y, const Angle& theta_z
        , const angle_util::Rotation3dType &rotation_type
        , const Eigen::Vector3d& v3_pos);

  /// @brief Determines which sign (positive or negative) a given 2D coordinate v2_pos falls into on a checkerboard.
  /// @param[in] v2_start_pos Checkerboard start point coordinates
  /// @param[in] v2_pattern_size Checkerboard pattern size
  /// @param[in] v2_pos Point coordinates to test
  /// @return 1: positive, -1: negative
  double get_sign_uniform_checkerboard_2d(
      const Eigen::Vector2d &v2_start_pos
    , const Eigen::Vector2d &v2_pattern_size
    , const Eigen::Vector2d &v2_pos );

  /// @brief Determines which sign (positive or negative) a given 3D coordinate v3_pos falls into on a checkerboard.
  /// @param[in] v3_start_pos Checkerboard start point coordinates
  /// @param[in] v3_pattern_size Checkerboard pattern size
  /// @param[in] v3_pos Point coordinates to test
  /// @return 1: positive, -1: negative
  double get_sign_uniform_checkerboard_3d(
      const Eigen::Vector3d &v3_start_pos
    , const Eigen::Vector3d &v3_pattern_size
    , const Eigen::Vector3d &v3_pos );

  ///@} ------------------------------------------------------------------

};
