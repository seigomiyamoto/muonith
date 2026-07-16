/// @file cls_Ray.hpp
/// @brief Ray representation for ray tracing
/// @details
/// Defines Ray2d and Ray3d classes representing rays with origin and direction vectors.
///
/// ## Overview
/// This module provides two ray classes for ray tracing operations:
/// - Ray2d: 2D ray with position and direction in Eigen::Vector2d
/// - Ray3d: 3D ray with position and direction in Eigen::Vector3d
///
/// ## Coordinate System
/// - Right-handed coordinate system
/// - Units: meters (assumed, application-dependent)
/// - Direction vectors are automatically normalized to unit length
///
/// ## Typical Workflow
/// 1. Create a ray with origin position and direction
/// 2. Test intersection with AABB (Axis-Aligned Bounding Box)
/// 3. Retrieve intersection parameters (t-values) or intersection points
/// 4. Use path factors for ray marching/traversal algorithms
///
/// ## Thread Safety
/// - All const methods are thread-safe for concurrent reads
/// - Non-const methods require external synchronization
///
/// @see cls_AABB.hpp for AABB2d and AABB3d classes
#pragma once
#include <Eigen/Dense>

#include <array>
#include <fstream>
#include <string>
#include <tuple>

#include "spdlog_pch.hpp"
#include "ns_mylogger.hpp"
#include "ns_mymacro.hpp"

#include "cls_AABB.hpp"

//##################################################################################
/// @class Ray2d
/// @brief Represents a 2D ray with origin position and direction vector.
/// @details
/// A Ray2d consists of:
/// - Position (origin): Eigen::Vector2d v2_pos
/// - Direction: Eigen::Vector2d v2_dir (not required to be unit length)
///
/// ## Usage Example
/// @code
/// Eigen::Vector2d origin(0.0, 0.0);
/// Eigen::Vector2d direction(1.0, 0.0);
/// Ray2d ray(origin, direction);
///
/// // Test intersection with AABB
/// AABB2d box(Eigen::Vector2d(1.0, -1.0), Eigen::Vector2d(3.0, 1.0));
/// auto [hit, tmin, tmax] = ray.is_intersect(box);
/// if (hit) {
///   Eigen::Vector2d entry = ray.getPointAtParameter(tmin);
/// }
/// @endcode
///
/// ## Thread Safety
/// - All const methods are thread-safe
/// - Non-const setters require external synchronization
///
/// @note Unlike Ray3d, direction is NOT automatically normalized
/// @ingroup geometryClasses
//##################################################################################
class Ray2d {
  private:

    /// @brief Origin position of the ray
    Eigen::Vector2d v2_pos = Eigen::Vector2d::Zero();

    /// @brief Direction vector of the ray
    Eigen::Vector2d v2_dir = Eigen::Vector2d::Zero();

  public:
  ///============================================================================
  /// @name constructor_destructor
  ///@{
  /// @brief Default constructor (origin and direction set to zero)
  Ray2d() = default;

  /// @brief Construct a ray from position and direction vectors
  /// @param[in] v2_pos_in Origin position of the ray
  /// @param[in] v2_dir_in Direction vector (not required to be normalized)
  /// @note Direction validation via checkDirection() is disabled because
  ///       Ray3d::toRay2d() may produce non-unit 2D projections
  Ray2d(const Eigen::Vector2d& v2_pos_in, const Eigen::Vector2d& v2_dir_in)
    : v2_pos(v2_pos_in), v2_dir(v2_dir_in) {
    };

  /// @brief Copy constructor
  Ray2d(const Ray2d& ray_in) = default;

  /// @brief Destructor
  ~Ray2d() = default;
  ///@} ------------------------------------------------------------------

  ///============================================================================
  /// @name operators
  ///@{

  /// @brief Assignment operator
  Ray2d& operator=(const Ray2d& other) = default;

  /// @brief Inequality operator
  /// @param[in] other Ray to compare with
  /// @return true if rays differ, false otherwise
  /// @note In debug builds (NODEBUG not defined), logs which member differs
  bool operator!=(const Ray2d& other) const {
#ifdef NODEBUG
    if (v2_pos != other.v2_pos) return true;
    if (v2_dir != other.v2_dir) return true;
#else
    if (v2_pos != other.v2_pos) { LOG_WARN("Ray2d: v2_pos differs"); return true; }
    if (v2_dir != other.v2_dir) { LOG_WARN("Ray2d: v2_dir differs"); return true; }
#endif
    return false;
  }

  /// @brief Equality operator
  /// @param[in] other Ray to compare with
  /// @return true if rays are equal, false otherwise
  bool operator==(const Ray2d& other) const {
    return !(*this != other);
  }
  ///@} ------------------------------------------------------------------

  //============================================================================
  /// @name getter_Ray2d
  /// @{
  /// @brief Get position vector
  /// @return Origin position of the ray
  Eigen::Vector2d pos() const { return v2_pos; }
  /// @brief Get direction vector
  /// @return Direction vector of the ray
  Eigen::Vector2d dir() const { return v2_dir; }
  /// @brief Get position component by index
  /// @param[in] i Component index (0=x, 1=y)
  /// @return Position component value
  double pos(const int i) const { return v2_pos(i); }
  /// @brief Get direction component by index
  /// @param[in] i Component index (0=x, 1=y)
  /// @return Direction component value
  double dir(const int i) const { return v2_dir(i); }
  /// @brief Get x component of position
  /// @return v2_pos.x()
  double  x() const { return v2_pos.x(); }
  /// @brief Get y component of position
  /// @return v2_pos.y()
  double  y() const { return v2_pos.y(); }
  /// @brief Get x component of direction
  /// @return v2_dir.x()
  double vx() const { return v2_dir.x(); }
  /// @brief Get y component of direction
  /// @return v2_dir.y()
  double vy() const { return v2_dir.y(); }

  /// @brief Convert position to string representation
  /// @param[in] width Minimum field width (0 = no width spec, default: 0)
  /// @param[in] precision Number of decimal places (default: 6)
  /// @param[in] type Format type: 'f' (fixed), 'e'/'E' (scientific), 'g'/'G' (general) (default: 'E')
  /// @return String in format "[x, y]"
  std::string pos_to_string(int width = 0, int precision = 6, char type = 'E') const;

  /// @brief Convert direction to string representation
  /// @param[in] width Minimum field width (0 = no width spec, default: 7)
  /// @param[in] precision Number of decimal places (default: 4)
  /// @param[in] type Format type: 'f', 'e', 'E', 'g', 'G' (default: 'f')
  /// @return String in format "[vx, vy]"
  std::string dir_to_string(int width = 7, int precision = 4, char type = 'f') const;

  /// @brief Convert ray to string representation for debugging
  /// @param[in] pos_width Minimum field width for position (0 = no width spec, default: 0)
  /// @param[in] pos_precision Number of decimal places for position (default: 6)
  /// @param[in] pos_type Format type for position: 'f', 'e', 'E', 'g', 'G' (default: 'E')
  /// @param[in] dir_width Minimum field width for direction (0 = no width spec, default: 7)
  /// @param[in] dir_precision Number of decimal places for direction (default: 4)
  /// @param[in] dir_type Format type for direction: 'f', 'e', 'E', 'g', 'G' (default: 'f')
  /// @return String in format "Ray2d(pos=[x, y], dir=[vx, vy])"
  std::string to_string(int pos_width = 0, int pos_precision = 6, char pos_type = 'E',
                        int dir_width = 7, int dir_precision = 4, char dir_type = 'f') const;

  ///@} ------------------------------------------------------------------

  //============================================================================
  /// @name setter_Ray2d
  /// @{
  /// @brief Set position from vector
  /// @param[in] v2_pos_in New position vector
  void set_pos(const Eigen::Vector2d& v2_pos_in) { v2_pos = v2_pos_in; }
  /// @brief Set direction from vector
  /// @param[in] v2_dir_in New direction vector
  void set_dir(const Eigen::Vector2d& v2_dir_in) { v2_dir = v2_dir_in; }
  /// @brief Set position from x,y components
  /// @param[in] x_in X component of new position
  /// @param[in] y_in Y component of new position
  void set_pos(const double& x_in, const double& y_in) {
    v2_pos << x_in, y_in;
  }
  /// @brief Set direction from vx,vy components
  /// @param[in] vx_in X component of new direction
  /// @param[in] vy_in Y component of new direction
  void set_dir(const double& vx_in, const double& vy_in) {
    v2_dir << vx_in, vy_in;
  }
  ///@} ------------------------------------------------------------------

  /// @brief Compute perpendicular distance from a point to this ray
  /// @param[in] point The point to measure distance to
  /// @return Perpendicular distance from point to ray line
  double distanceTo(const Eigen::Vector2d& point) const;

  /// @brief Project a point onto this ray
  /// @param[in] point The point to project
  /// @return Closest point on the ray to the given point
  Eigen::Vector2d projectOnto(const Eigen::Vector2d& point) const;

  /// @brief Get position along ray at parameter t
  /// @param[in] t Parameter value (position = origin + t * direction)
  /// @return Point along the ray at parameter t
  Eigen::Vector2d getPointAtParameter(const double t) const;

  //============================================================================
  /// @name path_length_Ray2d
  /// @{

  /// @brief Compute path length factors for each axis
  /// @details Computes factors to convert axis-aligned displacements to path length:
  /// - path_length = dx * factor[0]
  /// - path_length = dy * factor[1]
  /// @return Array of {path_factor_dx, path_factor_dy}
  /// @note Assumes direction is a unit vector for correct results
  /// @ingroup ray_tracing_functions
  std::array<double,2> get_path_factor() const;

  /// @brief Get sign of direction components
  /// @return Array of {sign(vx), sign(vy)} where sign is -1, 0, or +1
  /// @ingroup ray_tracing_functions
  std::array<double,2> get_sign() const;

  /// @brief Test intersection with 2D AABB and return intersection points
  /// @details Implements slab method for ray-AABB intersection test.
  /// Reference: https://gdbooks.gitbooks.io/3dcollisions/content/Chapter3/raycast_aabb.html
  /// @param[in] aabb2d The axis-aligned bounding box to test
  /// @param[out] v2_pos1 Entry intersection point (if hit)
  /// @param[out] v2_pos2 Exit intersection point (if hit)
  /// @return true if ray intersects the AABB, false otherwise
  /// @ingroup ray_tracing_functions
  bool is_intersect(
    const AABB2d &aabb2d, Eigen::Vector2d *v2_pos1, Eigen::Vector2d *v2_pos2 ) const;

  /// @brief Test intersection with 2D AABB and return t-parameters
  /// @param[in] aabb2d The axis-aligned bounding box to test
  /// @param[out] tmin_ret Entry parameter (point = pos + tmin * dir)
  /// @param[out] tmax_ret Exit parameter (point = pos + tmax * dir)
  /// @return true if ray intersects the AABB, false otherwise
  /// @ingroup ray_tracing_functions
  bool is_intersect(
    const AABB2d &aabb2d, double *tmin_ret, double *tmax_ret ) const;

  /// @brief Test intersection with 2D AABB (tuple return version)
  /// @param[in] aabb2d The axis-aligned bounding box to test
  /// @return Tuple of (hit, tmin, tmax)
  /// @ingroup ray_tracing_functions
  std::tuple<bool,double,double> is_intersect( const AABB2d &aabb2d ) const;

  ///@} ------------------------------------------------------------------

private:
  /// @brief Validate that direction vector is unit length
  /// @param[in] eps Tolerance for unit length check (default 1e-6)
  /// @throws std::invalid_argument if direction is not unit length
  void checkDirection( const double eps=1.0e-6 ) const;

};

//##################################################################################
/// @class Ray3d
/// @brief Represents a 3D ray with origin position and unit direction vector.
/// @details
/// A Ray3d consists of:
/// - Position (origin): Eigen::Vector3d v3_pos
/// - Direction: Eigen::Vector3d v3_dir (automatically normalized to unit length)
///
/// ## Usage Example
/// @code
/// Eigen::Vector3d origin(0.0, 0.0, 0.0);
/// Eigen::Vector3d direction(1.0, 1.0, 0.0);  // Will be normalized
/// Ray3d ray(origin, direction);
///
/// // Test intersection with AABB
/// AABB3d box(Eigen::Vector3d(1.0, 1.0, -1.0), Eigen::Vector3d(3.0, 3.0, 1.0));
/// auto [hit, tmin, tmax] = ray.is_intersect(box);
/// if (hit) {
///   Eigen::Vector3d entry = ray.getPointAtParameter(tmin);
///   Eigen::Vector3d exit  = ray.getPointAtParameter(tmax);
/// }
/// @endcode
///
/// ## Thread Safety
/// - All const methods are thread-safe
/// - Non-const setters require external synchronization
///
/// ## Binary I/O
/// - save()/load() methods for binary serialization
/// - Format: 6 doubles (3 for position, 3 for direction)
///
/// @note Direction is automatically normalized on construction and when set
/// @ingroup geometryClasses
//##################################################################################
class Ray3d {
  private:

    /// @brief Origin position of the ray
    Eigen::Vector3d v3_pos = Eigen::Vector3d::Zero();

    /// @brief Direction vector (always unit length)
    Eigen::Vector3d v3_dir = Eigen::Vector3d::Zero();

    /// @brief Normalize direction vector or throw if zero/non-finite
    /// @param[in] eps Minimum acceptable norm (default 1e-12)
    /// @throws std::invalid_argument if direction is zero or contains non-finite values
    void normalize_dir_or_throw(const double eps=1.0e-12);

  public:

    /// @brief Default constructor (origin and direction set to zero)
    Ray3d() = default;

    /// @brief Copy constructor
    Ray3d(const Ray3d& ray_in) = default;

    /// @brief Move constructor
    Ray3d(Ray3d&& ray_in) noexcept = default;

    /// @brief Destructor
    ~Ray3d() = default;

    /// @brief Construct a ray from position and direction vectors
    /// @param[in] v_pos_in Origin position of the ray
    /// @param[in] v3_dir_in Direction vector (will be normalized)
    /// @throws std::invalid_argument if direction is zero or non-finite
    Ray3d(const Eigen::Vector3d& v_pos_in, const Eigen::Vector3d& v3_dir_in)
      : v3_pos(v_pos_in), v3_dir(v3_dir_in) {
        normalize_dir_or_throw();
      }

    /// @brief Construct a ray from individual coordinates
    /// @param[in] x_in X coordinate of origin
    /// @param[in] y_in Y coordinate of origin
    /// @param[in] z_in Z coordinate of origin
    /// @param[in] vx_in X component of direction
    /// @param[in] vy_in Y component of direction
    /// @param[in] vz_in Z component of direction
    /// @throws std::invalid_argument if direction is zero or non-finite
    Ray3d(const double x_in, const double y_in, const double z_in,
          const double vx_in, const double vy_in, const double vz_in)
      : v3_pos(x_in, y_in, z_in), v3_dir(vx_in, vy_in, vz_in) {
        normalize_dir_or_throw();
      }

    /// @brief Copy assignment operator
    Ray3d& operator=(const Ray3d& other) {
      if (this == &other) return *this;
      v3_pos = other.v3_pos;
      v3_dir = other.v3_dir;
      return *this;
    }
  
    /// @brief Inequality operator
    /// @param[in] other Ray to compare with
    /// @return true if rays differ, false otherwise
    /// @note In debug builds (NODEBUG not defined), logs which member differs
    bool operator!=(const Ray3d& other) const {
#ifdef NODEBUG
      if (v3_pos != other.v3_pos) return true;
      if (v3_dir != other.v3_dir) return true;
#else
      if (v3_pos != other.v3_pos) { LOG_WARN("Ray3d: v3_pos differs"); return true; }
      if (v3_dir != other.v3_dir) { LOG_WARN("Ray3d: v3_dir differs"); return true; }
#endif
      return false;
    }

    /// @brief Equality operator
    /// @param[in] other Ray to compare with
    /// @return true if rays are equal, false otherwise
    bool operator==(const Ray3d& other) const {
      return !(*this != other);
    }

  public:
  //============================================================================
  /// @name getter_Ray3d
  /// @{

  /// @brief Get position vector
  /// @return Origin position of the ray
  Eigen::Vector3d pos() const { return v3_pos; }

  /// @brief Get direction vector (unit length)
  /// @return Direction vector of the ray
  Eigen::Vector3d dir() const { return v3_dir; }

  /// @brief Convert to 2D ray by ignoring z component
  /// @return Ray2d with xy components of position and direction
  /// @note The 2D direction may not be unit length after projection
  Ray2d toRay2d() const { return Ray2d(v3_pos.head(2), v3_dir.head(2)); }

  /// @brief Get 2D position (xy components)
  /// @return XY components of position
  Eigen::Vector2d pos2d() const { return v3_pos.head(2); }

  /// @brief Get 2D direction (xy components)
  /// @return XY components of direction
  Eigen::Vector2d dir2d() const { return v3_dir.head(2); }

  /// @brief Get position component by index
  /// @param[in] i Component index (0=x, 1=y, 2=z)
  /// @return Position component value
  double pos(const int i) const { return v3_pos(i); }

  /// @brief Get direction component by index
  /// @param[in] i Component index (0=x, 1=y, 2=z)
  /// @return Direction component value
  double dir(const int i) const { return v3_dir(i); }

  /// @brief Get x component of position
  /// @return v3_pos.x()
  double  x() const { return v3_pos.x(); }

  /// @brief Get y component of position
  /// @return v3_pos.y()
  double  y() const { return v3_pos.y(); }

  /// @brief Get z component of position
  /// @return v3_pos.z()
  double  z() const { return v3_pos.z(); }

  /// @brief Get x component of direction
  /// @return v3_dir.x()
  double vx() const { return v3_dir.x(); }

  /// @brief Get y component of direction
  /// @return v3_dir.y()
  double vy() const { return v3_dir.y(); }

  /// @brief Get z component of direction
  /// @return v3_dir.z()
  double vz() const { return v3_dir.z(); }

  /// @brief Convert position to string representation
  /// @param[in] width Minimum field width (0 = no width spec, default: 0)
  /// @param[in] precision Number of decimal places (default: 6)
  /// @param[in] type Format type: 'f', 'e', 'E', 'g', 'G' (default: 'E')
  /// @return String in format "[x, y, z]"
  std::string pos_to_string(int width = 0, int precision = 6, char type = 'E') const;

  /// @brief Convert direction to string representation
  /// @param[in] width Minimum field width (0 = no width spec, default: 7)
  /// @param[in] precision Number of decimal places (default: 4)
  /// @param[in] type Format type: 'f', 'e', 'E', 'g', 'G' (default: 'f')
  /// @return String in format "[vx, vy, vz]"
  std::string dir_to_string(int width = 7, int precision = 4, char type = 'f') const;

  /// @brief Convert ray to string representation for debugging
  /// @param[in] pos_width Minimum field width for position (0 = no width spec, default: 0)
  /// @param[in] pos_precision Number of decimal places for position (default: 6)
  /// @param[in] pos_type Format type for position: 'f', 'e', 'E', 'g', 'G' (default: 'E')
  /// @param[in] dir_width Minimum field width for direction (0 = no width spec, default: 7)
  /// @param[in] dir_precision Number of decimal places for direction (default: 4)
  /// @param[in] dir_type Format type for direction: 'f', 'e', 'E', 'g', 'G' (default: 'f')
  /// @return String in format "Ray3d(pos=[x, y, z], dir=[vx, vy, vz])"
  std::string to_string(int pos_width = 0, int pos_precision = 6, char pos_type = 'E',
                        int dir_width = 7, int dir_precision = 4, char dir_type = 'f') const;

  ///@} ------------------------------------------------------------------

  //============================================================================
  /// @name setter_Ray3d
  /// @{
  /// @brief Set position from vector
  /// @param[in] v3_pos_in New position vector
  void set_pos(const Eigen::Vector3d& v3_pos_in) { v3_pos = v3_pos_in; }

  /// @brief Set direction from vector (will be normalized)
  /// @param[in] v3_dir_in New direction vector
  /// @throws std::invalid_argument if direction is zero or non-finite
  void set_dir(const Eigen::Vector3d& v3_dir_in) { v3_dir = v3_dir_in; normalize_dir_or_throw(); }

  /// @brief Set position from x,y,z components
  /// @param[in] x_in X component of new position
  /// @param[in] y_in Y component of new position
  /// @param[in] z_in Z component of new position
  void set_pos(const double& x_in, const double& y_in, const double& z_in) {
    v3_pos << x_in, y_in, z_in;
  }

  /// @brief Set direction from vx,vy,vz components (will be normalized)
  /// @param[in] vx_in X component of new direction
  /// @param[in] vy_in Y component of new direction
  /// @param[in] vz_in Z component of new direction
  /// @throws std::invalid_argument if direction is zero or non-finite
  void set_dir(const double& vx_in, const double& vy_in, const double& vz_in) {
    v3_dir << vx_in, vy_in, vz_in; normalize_dir_or_throw();
  }
  ///@} ------------------------------------------------------------------

  /// @brief Compute perpendicular distance from a point to this ray
  /// @param[in] point The point to measure distance to
  /// @return Perpendicular distance from point to ray line
  double distanceTo(const Eigen::Vector3d& point) const;

  /// @brief Project a point onto this ray
  /// @param[in] point The point to project
  /// @return Closest point on the ray to the given point
  Eigen::Vector3d projectOnto(const Eigen::Vector3d& point) const;

  /// @brief Get position along ray at parameter t
  /// @param[in] t Parameter value (position = origin + t * direction)
  /// @return Point along the ray at parameter t
  Eigen::Vector3d getPointAtParameter(const double t) const;

  //============================================================================
  /// @name path_length_Ray3d
  /// @{

  /// @brief Test intersection with 3D AABB and return intersection points
  /// @details Implements slab method for ray-AABB intersection test.
  /// Reference: https://gdbooks.gitbooks.io/3dcollisions/content/Chapter3/raycast_aabb.html
  /// @param[in] aabb3d The axis-aligned bounding box to test
  /// @param[out] v3_pos1 Entry intersection point (if hit)
  /// @param[out] v3_pos2 Exit intersection point (if hit)
  /// @return true if ray intersects the AABB, false otherwise
  /// @ingroup ray_tracing_functions
  bool is_intersect(
    const AABB3d &aabb3d, Eigen::Vector3d *v3_pos1, Eigen::Vector3d *v3_pos2 ) const;

  /// @brief Test intersection with 3D AABB and return t-parameters
  /// @param[in] aabb3d The axis-aligned bounding box to test
  /// @param[out] tmin_ret Entry parameter (point = pos + tmin * dir)
  /// @param[out] tmax_ret Exit parameter (point = pos + tmax * dir)
  /// @return true if ray intersects the AABB, false otherwise
  /// @ingroup ray_tracing_functions
  bool is_intersect(
    const AABB3d &aabb3d, double *tmin_ret, double *tmax_ret ) const;

  /// @brief Debug version of intersection test with printf output
  /// @param[in] aabb3d The axis-aligned bounding box to test
  /// @param[out] tmin_ret Entry parameter
  /// @param[out] tmax_ret Exit parameter
  /// @return true if ray intersects the AABB, false otherwise
  /// @note Outputs debug information via printf (for development use only)
  /// @ingroup ray_tracing_functions
  bool is_intersect_debug( const AABB3d &aabb3d
    , double *tmin_ret, double *tmax_ret ) const;

  /// @brief Test intersection with 3D AABB (tuple return version)
  /// @param[in] aabb3d The axis-aligned bounding box to test
  /// @return Tuple of (hit, tmin, tmax)
  /// @ingroup ray_tracing_functions
  std::tuple<bool,double,double> is_intersect( const AABB3d &aabb3d ) const;

  /// @brief Compute path length factors for each axis
  /// @details Computes factors to convert axis-aligned displacements to path length:
  /// - path_length = dx * factor[0]
  /// - path_length = dy * factor[1]
  /// - path_length = dz * factor[2]
  /// @return Array of {path_factor_dx, path_factor_dy, path_factor_dz}
  /// @ingroup ray_tracing_functions
  std::array<double,3> get_path_factor() const;

  /// @brief Get sign of direction components
  /// @return Array of {sign(vx), sign(vy), sign(vz)} where sign is -1, 0, or +1
  /// @ingroup ray_tracing_functions
  std::array<double,3> get_sign() const;

  ///@} ------------------------------------------------------------------

  //============================================================================
  /// @name binary_io_functions
  /// @{

  /// @brief Write Eigen::Vector3d to binary output stream
  /// @param[in,out] ofs Output file stream (binary mode)
  /// @param[in] vec3d Vector to write
  /// @note Writes 3 doubles (24 bytes) in native byte order
  void write_vec3d( std::ofstream &ofs, const Eigen::Vector3d &vec3d ) const;

  /// @brief Read Eigen::Vector3d from binary input stream
  /// @param[in,out] ifs Input file stream (binary mode)
  /// @return Vector read from stream
  /// @note Reads 3 doubles (24 bytes) in native byte order
  Eigen::Vector3d read_vec3d( std::ifstream &ifs );

  /// @brief Save ray to binary output stream
  /// @param[in,out] ofs Output file stream (binary mode)
  /// @note Writes position then direction (48 bytes total)
  void save( std::ofstream &ofs ) const;

  /// @brief Load ray from binary input stream
  /// @param[in,out] ifs Input file stream (binary mode)
  /// @throws std::invalid_argument if loaded direction is zero or non-finite
  /// @note Reads position then direction (48 bytes total)
  void load( std::ifstream &ifs );

  ///@} ------------------------------------------------------------------
};
