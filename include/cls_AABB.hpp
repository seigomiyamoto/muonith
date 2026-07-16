/// @file cls_AABB.hpp
/// @brief Axis-Aligned Bounding Box (AABB) classes for 2D and 3D geometry
/// @details
/// This file defines AABB2d and AABB3d classes for representing axis-aligned bounding boxes.
/// These classes provide efficient spatial representations and queries for geometric calculations,
/// collision detection, and spatial partitioning.
///
/// **Key Features:**
/// - AABB2d: 2D rectangular bounding box
/// - AABB3d: 3D cuboid bounding box
/// - Point containment testing (axis-aligned and rotated)
/// - Overlap detection between AABBs
/// - Adjacency detection for neighboring AABBs
///
/// **Thread Safety:**
/// - All const methods are thread-safe
/// - Non-const methods (setters) are not thread-safe
/// - No internal synchronization is provided
///
/// **Coordinate System:**
/// - Right-handed coordinate system assumed for 3D operations
/// - No specific units enforced (user must maintain consistency)
#pragma once

#include <map>
#include <fstream>
#include <iostream>
#include <sstream> // istringstream
#include <string>
#include <cstdio>
#include <cmath>
#include <functional>  //for sorting
#include <algorithm>//for sorting
#include <vector>

#include <Eigen/Dense>
#include "cls_Angle.hpp"



//###############################################################

/// @brief definition of vertex of AABB2d
using AABB2DVertices = std::array<Eigen::Vector2d, 4>;

/// @brief definition of vertex of AABB3d
using AABB3DVertices = std::array<Eigen::Vector3d, 8>;
//###############################################################


//###############################################################
/// @class AABB2d
/// @brief Class representing a rectangle as a 2D Axis-Aligned Bounding Box.
/// @ingroup geometryClasses
///
/// @details
/// AABB2d represents a 2D rectangle aligned with the coordinate axes, defined by
/// minimum and maximum corner points. It provides efficient spatial queries including
/// point containment (for both axis-aligned and rotated rectangles), overlap detection,
/// and adjacency testing.
///
/// **Invariants:**
/// - v2_min.x() <= v2_max.x()
/// - v2_min.y() <= v2_max.y()
/// - Violations are detected in checkMinMax() and throw an error
///
/// **Primary Use Cases:**
/// - 2D collision detection
/// - Spatial partitioning and quadtree nodes
/// - Screen/viewport bounds checking
///
/// **Example Usage:**
/// @code
/// // Create AABB from corner points
/// AABB2d box(Eigen::Vector2d(0, 0), Eigen::Vector2d(10, 10));
///
/// // Check if point is inside
/// Eigen::Vector2d point(5, 5);
/// if (box.is_inside(point)) {
///   // Point is inside the box
/// }
///
/// // Check overlap with another AABB
/// AABB2d other(Eigen::Vector2d(8, 8), Eigen::Vector2d(12, 12));
/// if (box.is_overlap(other)) {
///   // Boxes overlap
/// }
/// @endcode
///
/// @note Thread-safety: const methods are thread-safe, setters are not
//###############################################################
class AABB2d {
private:
  /// @brief lower left edge of AABB2d
  Eigen::Vector2d v2_min = Eigen::Vector2d(0.0, 0.0);

  /// @brief upper right edge of AABB2d
  Eigen::Vector2d v2_max = Eigen::Vector2d(0.0, 0.0);
  
public:
  //==================================================================
  /// @name constructor_destructor
  ///@{

  /// @brief default constructor
  AABB2d() = default;

  /// @brief copy constructor
  AABB2d(const AABB2d &other) = default;

  /// @brief move constructor
  AABB2d(AABB2d &&other) noexcept = default;

  /// @brief destructor
  virtual ~AABB2d() = default;

  /// @brief constructor from two Eigen::Vector2d
  /// @param[in] v2_min_in lower edge of AABB2d
  /// @param[in] v2_max_in upper edge of AABB2d
  AABB2d(const Eigen::Vector2d& v2_min_in, const Eigen::Vector2d& v2_max_in)
    : v2_min(v2_min_in), v2_max(v2_max_in) { checkMinMax(); };

  /// @brief constructor2
  /// @param[in] x_min lower edge of x
  /// @param[in] x_max upper edge of x
  /// @param[in] y_min lower edge of y
  /// @param[in] y_max upper edge of y
  AABB2d(const double& x_min, const double& x_max, const double& y_min, const double& y_max)
    : v2_min(x_min, y_min), v2_max(x_max, y_max) { checkMinMax(); };
  
  ///@} ------------------------------------------------------------------
  //======================================================================
  /// @name operators
  ///@{

  /// @brief copy assignment operator
  AABB2d& operator=(const AABB2d& other) = default;

  /// @brief move assignment operator
  AABB2d& operator=(AABB2d&& other) noexcept = default;

  /// @brief inequality operator
  bool operator!=(const AABB2d& other) const;

  /// @brief equality operator
  bool operator==(const AABB2d& other) const {
    return !(*this != other);
  };
  ///@} ------------------------------------------------------------------
  
  //==================================================================
  /// @name getter_functions
  ///@{

  /// @brief return lower edge of AABB2d
  Eigen::Vector2d min() const { return v2_min; };

  /// @brief  return upper edge of AABB2d
  Eigen::Vector2d max() const { return v2_max; };

  /// @brief return lower edge of AABB2d
  /// @param[in] i dimension index (0=x, 1=y)
  /// @return lower edge value for the specified dimension
  /// @throws std::runtime_error if i is out of range [0, 1]
  double min( const int i ) const;

  /// @brief return upper edge of AABB2d
  /// @param[in] i dimension index (0=x, 1=y)
  /// @return upper edge value for the specified dimension
  /// @throws std::runtime_error if i is out of range [0, 1]
  double max( const int i ) const;

  /// @brief return lower edge of x
  double xmin() const { return v2_min.x(); };

  /// @brief return lower edge of y
  double ymin() const { return v2_min.y(); };

  /// @brief return upper edge of x
  double xmax() const { return v2_max.x(); };

  /// @brief return upper edge of y
  double ymax() const { return v2_max.y(); };

  /// @brief get the center of x
  double xcnt() const { return 0.5*(v2_min.x()+v2_max.x()); };

  /// @brief get the center of y
  double ycnt() const { return 0.5*(v2_min.y()+v2_max.y()); };

  /// @brief get the length of x
  double xlen() const { return v2_max.x()-v2_min.x(); };

  /// @brief get the length of y
  double ylen() const { return v2_max.y()-v2_min.y(); };

  /// @brief get the vertices of AABB2d
  /// @return std::array<Eigen::Vector2d, 4> containing 4 vertices in order: (x-,y-), (x+,y-), (x-,y+), (x+,y+)
  /// @note Vertices are ordered counter-clockwise starting from lower-left corner
  AABB2DVertices get_vertices() const;

  /// @brief get the xlow, xup, ylow, yup bounds
  /// @return std::array<double, 4> containing {xmin, xmax, ymin, ymax}
  std::array<double,4> get_xlow_xup_ylow_yup() const {
    return {v2_min.x(),v2_max.x(),v2_min.y(),v2_max.y()};
  };

  /// @brief return center coordinate of AABB2d
  Eigen::Vector2d center() const { return 0.5*(v2_min+v2_max); };

  ///@} ------------------------------------------------------------------

  //==================================================================
  /// @name setter_functions
  ///@{

  /// @brief set lower edge of x
  void set_xmin( const double x_min ) { v2_min.x() = x_min; };

  /// @brief set lower edge of y
  void set_ymin( const double y_min ) { v2_min.y() = y_min; };

  /// @brief set upper edge of x
  void set_xmax( const double x_max ) { v2_max.x() = x_max; };

  /// @brief set upper edge of y
  void set_ymax( const double y_max ) { v2_max.y() = y_max; };

  /// @brief set lower edge of AABB2d
  void set_min( const Eigen::Vector2d& v2_min_in ) { v2_min = v2_min_in; };

  /// @brief set upper edge of AABB2d
  void set_max( const Eigen::Vector2d& v2_max_in ) { v2_max = v2_max_in; };

  /// @brief set both edge of AABB2d
  void set( const AABB2d& other ) { v2_min = other.v2_min; v2_max = other.v2_max; };

  ///@} ------------------------------------------------------------------

  /// @brief Check if a point is inside the rectangle
  /// @details The rectangle is assumed to be parallel to the x-axis.
  /// @param[in] v2_pos point position to test
  /// @return true if point is inside (xmin <= x < xmax and ymin <= y < ymax), false otherwise
  /// @note Uses half-open interval: includes minimum edge, excludes maximum edge
  bool is_inside(const Eigen::Vector2d &v2_pos ) const;

  /// @brief Check if a point is inside the rectangle
  /// @param[in] x x-coordinate of the point
  /// @param[in] y y-coordinate of the point
  /// @return true if point is inside, false otherwise
  bool is_inside(const double x, const double y ) const {
    return is_inside(Eigen::Vector2d(x,y));
  };

  /// @brief Check if a point is inside a rotated rectangle
  /// @details The rectangle is assumed to be rotated counter-clockwise by theta [deg] from the x-axis.
  ///          The point is inverse-rotated and then tested against the axis-aligned bounds.
  /// @param[in] v2_pos point position to test
  /// @param[in] ang_theta rotation angle in degrees (counter-clockwise from x-axis)
  /// @return true if point is inside the rotated rectangle, false otherwise
  /// @note Rotation is performed around the rectangle's center
  bool is_inside(
    const Eigen::Vector2d &v2_pos, const Angle &ang_theta) const;

  /// @brief check if this AABB overlaps with another
  /// @param[in] other the other AABB2d to test against
  /// @return true if overlap exists, false otherwise
  /// @note Uses separating axis theorem: returns false if any axis shows separation
  bool is_overlap(const AABB2d &other) const;

  /// @brief enum class for adjacency
  enum class Adjacency {
    None,       ///< Not adjacent
    Xm,         ///< Other rectangle is adjacent on the left
    Xp,         ///< Other rectangle is adjacent on the right
    Ym,         ///< Other rectangle is adjacent below
    Yp,         ///< Other rectangle is adjacent above
    Overlap     ///< Rectangles overlap (intersect)
  };

  /// @brief check adjacency relationship with another AABB2d
  /// @param other the other AABB2d to test against
  /// @return Adjacency enum indicating the relationship (None, Xm, Xp, Ym, Yp, or Overlap)
  /// @note Adjacency requires exact edge alignment (floating-point equality)
  /// @note Returns Overlap if AABBs intersect
  Adjacency is_adjacent(const AABB2d &other) const;

private:
  /// @brief validate that min <= max for all dimensions
  /// @throws std::runtime_error if any min > max
  void checkMinMax() const;
};

//##################################################################
/// @class AABB3d
/// @brief Class representing a rectangular cuboid as a 3D Axis-Aligned Bounding Box.
/// @ingroup geometryClasses
///
/// @details
/// AABB3d represents a 3D rectangular cuboid aligned with the coordinate axes, defined by
/// minimum and maximum corner points. It provides efficient 3D spatial queries including
/// point containment, overlap detection, and adjacency testing for neighboring cuboids.
///
/// **Invariants:**
/// - v3_min.x() <= v3_max.x()
/// - v3_min.y() <= v3_max.y()
/// - v3_min.z() <= v3_max.z()
/// - Violations are detected in checkMinMax() and throw an error
///
/// **Primary Use Cases:**
/// - 3D collision detection and physics simulations
/// - Spatial partitioning and octree nodes
/// - Voxel representation and volume bounds
///
/// **Example Usage:**
/// @code
/// // Create AABB from corner points
/// AABB3d box(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(10, 10, 10));
///
/// // Check if point is inside
/// Eigen::Vector3d point(5, 5, 5);
/// if (box.is_inside(point)) {
///   // Point is inside the box
/// }
///
/// // Check overlap with another AABB
/// AABB3d other(Eigen::Vector3d(8, 8, 8), Eigen::Vector3d(12, 12, 12));
/// if (box.is_overlap(other)) {
///   // Boxes overlap
/// }
/// @endcode
///
/// @note Thread-safety: const methods are thread-safe, setters are not
/// @note Coordinate system: Right-handed coordinate system with z-up convention
//##################################################################
class AABB3d {
private:
  /// @brief lower x/y/z edge of AABB3d
  Eigen::Vector3d v3_min = Eigen::Vector3d(0.0, 0.0, 0.0); // min

  /// @brief upper x/y/z edge of AABB3d
  Eigen::Vector3d v3_max = Eigen::Vector3d(0.0, 0.0, 0.0); // max

public:
  //==================================================================
  /// @name constructor_destructor
  ///@{

  /// @brief default constructor
  AABB3d() = default;

  /// @brief constructor from two Eigen::Vector3d
  /// @param[in] v3_min_in lower edge of AABB3d
  /// @param[in] v3_max_in upper edge of AABB3d
  AABB3d(const Eigen::Vector3d& v3_min_in, const Eigen::Vector3d& v3_max_in)
    : v3_min(v3_min_in), v3_max(v3_max_in) {
      checkMinMax();
  };

  /// @brief constructor2
  /// @param[in] x_min lower edge of x
  /// @param[in] x_max upper edge of x
  /// @param[in] y_min lower edge of y
  /// @param[in] y_max upper edge of y
  /// @param[in] z_min lower edge of z
  /// @param[in] z_max upper edge of z
  AABB3d(const double& x_min, const double& x_max
       , const double& y_min, const double& y_max
       , const double& z_min, const double& z_max )
  : v3_min(x_min, y_min, z_min), v3_max(x_max, y_max, z_max) { checkMinMax(); };

  /// @brief copy constructor
  AABB3d(const AABB3d& other) = default;

  /// @brief move constructor
  AABB3d(AABB3d&& other) noexcept = default;

  /// @brief destructor
  ~AABB3d() = default;
  ///@} ------------------------------------------------------------------
  //======================================================================
  /// @name operators
  ///@{

  /// @brief copy assignment operator
  AABB3d& operator=(const AABB3d& other) = default;

  /// @brief move assignment operator
  AABB3d& operator=(AABB3d&& other) noexcept = default;

  /// @brief inequality operator
  bool operator!=(const AABB3d& other) const;

  /// @brief equality operator
  bool operator==(const AABB3d& other) const {
    return !(*this != other);
  };
  ///@} ------------------------------------------------------------------

  //==================================================================
  /// @name getter_functions
  ///@{
  
  /// @brief return lower edge of AABB3d
  Eigen::Vector3d min() const { return v3_min; };

  /// @brief return upper edge of AABB3d
  Eigen::Vector3d max() const { return v3_max; };

  /// @brief return lower edge of AABB3d
  /// @param[in] i dimension index (0=x, 1=y, 2=z)
  /// @return lower edge value for the specified dimension
  /// @throws std::runtime_error if i is out of range [0, 2]
  double min( const int i ) const;

  /// @brief return upper edge of AABB3d
  /// @param[in] i dimension index (0=x, 1=y, 2=z)
  /// @return upper edge value for the specified dimension
  /// @throws std::runtime_error if i is out of range [0, 2]
  double max( const int i ) const;

  /// @brief return lower edge of x
  double xmin() const { return v3_min.x(); };

  /// @brief return lower edge of y
  double ymin() const { return v3_min.y(); };

  /// @brief return lower edge of z
  double zmin() const { return v3_min.z(); };

  /// @brief return upper edge of x
  double xmax() const { return v3_max.x(); };

  /// @brief return upper edge of y
  double ymax() const { return v3_max.y(); };

  /// @brief return upper edge of z
  double zmax() const { return v3_max.z(); };

  /// @brief get the center of x
  double xcnt() const { return 0.5*(v3_min.x()+v3_max.x()); };

  /// @brief get the center of y
  double ycnt() const { return 0.5*(v3_min.y()+v3_max.y()); };

  /// @brief get the center of z
  double zcnt() const { return 0.5*(v3_min.z()+v3_max.z()); };

  /// @brief get the length of x
  double xlen() const { return v3_max.x()-v3_min.x(); };

  /// @brief get the length of y
  double ylen() const { return v3_max.y()-v3_min.y(); };

  /// @brief get the length of z
  double zlen() const { return v3_max.z()-v3_min.z(); };

  /// @brief return center coordinate of AABB3d
  Eigen::Vector3d center() const { return 0.5*(v3_min+v3_max); };

  /// @brief get the vertices of AABB3d as std::array of Eigen::Vector3d
  /// @return std::array containing 8 vertices in order:
  ///         (x-,y-,z-), (x+,y-,z-), (x-,y+,z-), (x+,y+,z-),
  ///         (x-,y-,z+), (x+,y-,z+), (x-,y+,z+), (x+,y+,z+)
  /// @note First 4 vertices are at z_min, last 4 at z_max
  /// @note Vertices ordered by increasing x, then y, then z
  AABB3DVertices get_vertices() const;

  ///@} ------------------------------------------------------------------

  //==================================================================
  /// @name setter_functions
  ///@{
  
  /// @brief set lower edge of x
  /// @details v3_min.x() = x_min;
  void set_xmin( const double x_min ) { v3_min.x() = x_min; };

  /// @brief set lower edge of x
  /// @details v3_max.x() = x_max;
  void set_xmax( const double x_max ) { v3_max.x() = x_max; };

  /// @brief set lower edge of y
  /// @details v3_min.y() = y_min;
  void set_ymin( const double y_min ) { v3_min.y() = y_min; };

  /// @brief set lower edge of y
  /// @details v3_max.y() = y_max;
  void set_ymax( const double y_max ) { v3_max.y() = y_max; };

  /// @brief set lower edge of z
  /// @details v3_min.z() = z_min;
  void set_zmin( const double z_min ) { v3_min.z() = z_min; };

  /// @brief set lower edge of z
  /// @details v3_max.z() = z_max;
  void set_zmax( const double z_max ) { v3_max.z() = z_max; };
  ///@} ------------------------------------------------------------------
  
  /// @brief Check if a point is inside the cuboid
  /// @details The edges of the cuboid are assumed to be parallel to the x, y, z axes.
  /// @param[in] v3_pos point position to test
  /// @return true if point is inside (xmin <= x < xmax and ymin <= y < ymax and zmin <= z < zmax), false otherwise
  /// @note Uses half-open interval: includes minimum edges, excludes maximum edges
  bool is_inside(const Eigen::Vector3d &v3_pos ) const;

  /// @brief check if this AABB overlaps with another
  /// @param[in] other the other AABB3d to test against
  /// @return true if overlap exists, false otherwise
  /// @note Uses separating axis theorem: returns false if any axis shows separation
  bool is_overlap(const AABB3d &other) const;

  /// @brief enum class for adjacency
  enum class Adjacency{
    None,      ///< Not adjacent
    Xm,        ///< Other cuboid is adjacent on the left
    Xp,        ///< Other cuboid is adjacent on the right
    Zm,        ///< Other cuboid is adjacent below
    Zp,        ///< Other cuboid is adjacent above
    Ym,        ///< Other cuboid is adjacent in front
    Yp,        ///< Other cuboid is adjacent behind
    Overlap    ///< Cuboids overlap (intersect)
  };

  /// @brief check adjacency relationship with another AABB3d
  /// @param other the other AABB3d to test against
  /// @return Adjacency enum indicating the relationship (None, Xm, Xp, Ym, Yp, Zm, Zp, or Overlap)
  /// @note Adjacency requires exact edge alignment (floating-point equality)
  /// @note Returns Overlap if AABBs intersect
  Adjacency is_adjacent(const AABB3d &other) const;

  /// @brief display AABB bounds to FILE stream
  /// @param[in] fout output FILE* stream (must be valid and open)
  /// @note Output format: "function_name | x:min-max, y:min-max, z:min-max\n"
  void disp(FILE *fout) const;

private:
  /// @brief validate that min <= max for all dimensions
  /// @throws std::runtime_error if any min > max
  void checkMinMax() const;

};
