/// @file st_shapes.hpp
/// @brief Geometric shape parameter structures for elliptical geometries
/// @details Defines parameter structures for horizontal ellipse disks and vertical elliptical cylinders.
///
/// This module provides simple POD-like structures representing:
/// - HorizontalEllipseDisk: A horizontal elliptical disk at a fixed z-plane
/// - VerticalEllipticCylinderOpen: A vertical elliptical cylinder without end caps
/// - VerticalEllipticCylinderCapped: A vertical elliptical cylinder with end caps
///
/// @note Coordinate system: Right-handed with z-axis vertical (z-up)
/// @note Units: All distances in meters, angles in radians, density in kg/m³
/// @note Thread-safety: All structures are trivially copyable and thread-safe for read operations
#pragma once
#include <array>
#include <Eigen/Dense>
#include "cls_Angle.hpp"

/// @brief Horizontal elliptical disk parameter structure
/// @details Represents an elliptical disk lying in a horizontal plane at z = z_plane.
///
/// The ellipse is defined in the local coordinate frame with semi-axes a and b,
/// then rotated by yaw around the z-axis, and finally translated to (xcnt, ycnt, z_plane).
///
/// @note Coordinate system: Right-handed, z-up
/// @note Units: xcnt, ycnt, z_plane, a, b in meters; yaw in radians
/// @note Thread-safety: Trivially copyable, safe for concurrent reads
struct HorizontalEllipseDisk {
  double xcnt, ycnt, z_plane;                 ///< Center position (x, y) and plane z-coordinate [meters]
  double a, b;                                 ///< Semi-axes: a (local x-axis), b (local y-axis) [meters]
  Angle  yaw = Angle(0.0, Angle::Unit::Radian);///< Rotation angle around z-axis (right-handed) [radians]
  double eps = 1e-12;                          ///< Numerical tolerance for geometric operations [dimensionless]

  //==================================================================
  /// @name Constructors
  ///@{

  /// @brief Default constructor
  /// @details Initializes a unit circle at the origin in the z=0 plane.
  ///
  /// Default values: center=(0,0,0), a=1.0m, b=1.0m, yaw=0rad, eps=1e-12
  inline HorizontalEllipseDisk() noexcept
    : xcnt(0.0), ycnt(0.0), z_plane(0.0), a(1.0), b(1.0)
    , yaw(Angle(0.0, Angle::Unit::Radian)), eps(1e-12) {}

  /// @brief Parameterized constructor
  /// @param[in] xcnt_in Center x-coordinate [meters]
  /// @param[in] ycnt_in Center y-coordinate [meters]
  /// @param[in] z_plane_in Plane z-coordinate [meters]
  /// @param[in] a_in Semi-axis along local x-axis [meters]
  /// @param[in] b_in Semi-axis along local y-axis [meters]
  /// @param[in] yaw_in Rotation angle around z-axis [radians] (default: 0)
  /// @param[in] eps_in Numerical tolerance (default: 1e-12)
  inline HorizontalEllipseDisk(
      double xcnt_in, double ycnt_in, double z_plane_in
    , double a_in, double b_in
    , Angle yaw_in = Angle(0.0, Angle::Unit::Radian)
    , double eps_in = 1e-12) noexcept
    : xcnt(xcnt_in), ycnt(ycnt_in), z_plane(z_plane_in), a(a_in), b(b_in)
      , yaw(yaw_in), eps(eps_in) {}

  /// @brief Convenience constructor with 2D center vector
  /// @param[in] v2_pos_cnt 2D center position (x, y) [meters]
  /// @param[in] z_plane_in Plane z-coordinate [meters]
  /// @param[in] a_in Semi-axis along local x-axis [meters]
  /// @param[in] b_in Semi-axis along local y-axis [meters]
  /// @param[in] yaw_in Rotation angle around z-axis [radians] (default: 0)
  /// @param[in] eps_in Numerical tolerance (default: 1e-12)
  inline HorizontalEllipseDisk(
      const Eigen::Vector2d& v2_pos_cnt
    , double z_plane_in
    , double a_in, double b_in
    , Angle yaw_in = Angle(0.0, Angle::Unit::Radian)
    , double eps_in = 1e-12) noexcept
  : xcnt(v2_pos_cnt.x()), ycnt(v2_pos_cnt.y()), z_plane(z_plane_in),
    a(a_in), b(b_in), yaw(yaw_in), eps(eps_in) {}

  /// @brief Copy constructor (defaulted)
  inline HorizontalEllipseDisk(const HorizontalEllipseDisk& rhs) noexcept = default;
  ///@}
  //==================================================================
};

/// @brief Vertical elliptical cylinder parameter structure (open, without end caps)
/// @details Represents a vertical elliptical cylinder extending from z_base to z_base+height.
///
/// The ellipse is defined in the local horizontal coordinate frame with semi-axes a and b,
/// then rotated by yaw around the z-axis, and finally translated to (xcnt, ycnt).
/// The cylinder extends vertically over z ∈ [z_base, z_base+height] but does not include
/// top or bottom caps.
///
/// @note Coordinate system: Right-handed, z-up
/// @note Units: xcnt, ycnt, z_base, height, a, b in meters; yaw in radians; density in kg/m³
/// @note Thread-safety: Trivially copyable, safe for concurrent reads
struct VerticalEllipticCylinderOpen {
  double xcnt, ycnt;            ///< Center position (x, y) [meters]
  double z_base, height;        ///< Vertical extent: z ∈ [z_base, z_base+height] [meters]
  double a, b;                  ///< Semi-axes: a (local x-axis), b (local y-axis) [meters]
  Angle  yaw = Angle(0.0, Angle::Unit::Radian); ///< Rotation angle around z-axis (right-handed) [radians]
  double density = 0.0;         ///< Material density [kg/m³]
  double eps = 1e-12;           ///< Numerical tolerance for geometric operations [dimensionless]

  /// @name Constructors
  ///@{

  /// @brief Default constructor
  /// @details Initializes a unit circular cylinder from z=0 to z=1 at the origin, with zero density.
  ///
  /// Default values: center=(0,0), z_base=0, height=1.0m, a=1.0m, b=1.0m, yaw=0rad, density=0, eps=1e-12
  inline VerticalEllipticCylinderOpen() noexcept
  : xcnt(0.0), ycnt(0.0),
    z_base(0.0), height(1.0),
    a(1.0), b(1.0),
    yaw(Angle(0.0, Angle::Unit::Radian)),
    density(0.0), eps(1e-12) {}

  /// @brief Parameterized constructor
  /// @param[in] xcnt_in Center x-coordinate [meters]
  /// @param[in] ycnt_in Center y-coordinate [meters]
  /// @param[in] z_base_in Base z-coordinate [meters]
  /// @param[in] height_in Cylinder height [meters]
  /// @param[in] a_in Semi-axis along local x-axis [meters]
  /// @param[in] b_in Semi-axis along local y-axis [meters]
  /// @param[in] yaw_in Rotation angle around z-axis [radians] (default: 0)
  /// @param[in] density_in Material density [kg/m³] (default: 0)
  /// @param[in] eps_in Numerical tolerance (default: 1e-12)
  inline VerticalEllipticCylinderOpen(
      double xcnt_in, double ycnt_in
    , double z_base_in, double height_in
    , double a_in, double b_in
    , Angle yaw_in = Angle(0.0, Angle::Unit::Radian)
    , double density_in = 0.0
    , double eps_in = 1e-12) noexcept
  : xcnt(xcnt_in), ycnt(ycnt_in),
    z_base(z_base_in), height(height_in),
    a(a_in), b(b_in), yaw(yaw_in),
    density(density_in), eps(eps_in) {}

  /// @brief Convenience constructor with 2D center vector
  /// @param[in] v2_pos_cnt 2D center position (x, y) [meters]
  /// @param[in] z_base_in Base z-coordinate [meters]
  /// @param[in] height_in Cylinder height [meters]
  /// @param[in] a_in Semi-axis along local x-axis [meters]
  /// @param[in] b_in Semi-axis along local y-axis [meters]
  /// @param[in] yaw_in Rotation angle around z-axis [radians] (default: 0)
  /// @param[in] density_in Material density [kg/m³] (default: 0)
  /// @param[in] eps_in Numerical tolerance (default: 1e-12)
  inline VerticalEllipticCylinderOpen(
      const Eigen::Vector2d& v2_pos_cnt
    , double z_base_in, double height_in
    , double a_in, double b_in
    , Angle yaw_in = Angle(0.0, Angle::Unit::Radian)
    , double density_in = 0.0
    , double eps_in = 1e-12) noexcept
  : xcnt(v2_pos_cnt.x()), ycnt(v2_pos_cnt.y()),
    z_base(z_base_in), height(height_in),
    a(a_in), b(b_in), yaw(yaw_in),
    density(density_in), eps(eps_in) {}

  /// @brief Copy constructor (defaulted)
  inline VerticalEllipticCylinderOpen(const VerticalEllipticCylinderOpen& rhs) noexcept = default;
  ///@}
};

/// @brief Vertical elliptical cylinder parameter structure (capped, with end caps)
/// @details Represents a vertical elliptical cylinder with closed top and bottom caps.
///
/// Inherits from VerticalEllipticCylinderOpen and adds the semantics that the cylinder
/// includes elliptical caps at z = z_base and z = z_base + height.
///
/// @note Coordinate system: Right-handed, z-up
/// @note Units: All dimensions in meters, angles in radians, density in kg/m³
/// @note Thread-safety: Trivially copyable, safe for concurrent reads
struct VerticalEllipticCylinderCapped : VerticalEllipticCylinderOpen {
  /// @name Constructors
  ///@{

  /// @brief Default constructor
  /// @details Inherits default values from VerticalEllipticCylinderOpen.
  inline VerticalEllipticCylinderCapped() noexcept = default;

  /// @brief Parameterized constructor
  /// @param[in] xcnt_in Center x-coordinate [meters]
  /// @param[in] ycnt_in Center y-coordinate [meters]
  /// @param[in] z_base_in Base z-coordinate [meters]
  /// @param[in] height_in Cylinder height [meters]
  /// @param[in] a_in Semi-axis along local x-axis [meters]
  /// @param[in] b_in Semi-axis along local y-axis [meters]
  /// @param[in] yaw_in Rotation angle around z-axis [radians] (default: 0)
  /// @param[in] density_in Material density [kg/m³] (default: 0)
  /// @param[in] eps_in Numerical tolerance (default: 1e-12)
  inline VerticalEllipticCylinderCapped(
      double xcnt_in, double ycnt_in
    , double z_base_in, double height_in
    , double a_in, double b_in
    , Angle yaw_in = Angle(0.0, Angle::Unit::Radian)
    , double density_in = 0.0
    , double eps_in = 1e-12) noexcept {
    this->xcnt = xcnt_in; this->ycnt = ycnt_in;
    this->z_base = z_base_in; this->height = height_in;
    this->a = a_in; this->b = b_in;
    this->yaw = yaw_in;
    this->density = density_in;
    this->eps = eps_in;
  }

  /// @brief Convenience constructor with 2D center vector
  /// @param[in] v2_pos_cnt 2D center position (x, y) [meters]
  /// @param[in] z_base_in Base z-coordinate [meters]
  /// @param[in] height_in Cylinder height [meters]
  /// @param[in] a_in Semi-axis along local x-axis [meters]
  /// @param[in] b_in Semi-axis along local y-axis [meters]
  /// @param[in] yaw_in Rotation angle around z-axis [radians] (default: 0)
  /// @param[in] density_in Material density [kg/m³] (default: 0)
  /// @param[in] eps_in Numerical tolerance (default: 1e-12)
  inline VerticalEllipticCylinderCapped(
      const Eigen::Vector2d& v2_pos_cnt
    , double z_base_in, double height_in
    , double a_in, double b_in
    , Angle yaw_in = Angle(0.0, Angle::Unit::Radian)
    , double density_in = 0.0
    , double eps_in = 1e-12) noexcept {
    this->xcnt = v2_pos_cnt.x(); this->ycnt = v2_pos_cnt.y();
    this->z_base = z_base_in; this->height = height_in;
    this->a = a_in; this->b = b_in;
    this->yaw = yaw_in;
    this->density = density_in;
    this->eps = eps_in;
  }

  /// @brief Copy constructor (defaulted)
  /// @note Also allows implicit conversion from VerticalEllipticCylinderOpen via inherited copy constructor
  inline VerticalEllipticCylinderCapped(const VerticalEllipticCylinderCapped& rhs) noexcept = default;
  ///@}
};
