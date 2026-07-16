/// @file cls_Grid2d.hpp
/// @brief Two-dimensional grid class for spatial discretization
/// @details
/// Provides a 2D grid structure composed of two orthogonal Grid1d axes (x_axis and y_axis).
/// This class supports spatial discretization for ray tracing, imaging, and voxel-based computations.
///
/// Key features:
/// - Grid cell indexing via (ix, iy) pairs
/// - Coordinate-to-index conversion and vice versa
/// - Ray tracing support via voxel traversal algorithms
/// - Circle and rectangular region queries
/// - Grid cutting and merging operations
/// - Binary I/O for serialization
///
/// @note Thread-safety: No. Concurrent access to mutable operations is not safe.
/// @note Units: Physical coordinates (x, y) are unitless; interpretation depends on application context.
/// @note Memory layout: Cells are accessed via [ix][iy]; iteration order is application-dependent.
#pragma once

#include <cstdio>
#include <map>
#include <fstream>
#include <iostream>
#include <sstream> // istringstream
#include <string>
#include <cmath>
#include <functional>  //for sorting
#include <algorithm>//for sorting
#include <vector>
#include <source_location>

#include <Eigen/Dense>
#include "spdlog_pch.hpp"
#include "ns_mylogger.hpp"
#include "ns_mymacro.hpp"
#include "cls_AABB.hpp"
#include "cls_Ray.hpp"
#include "cls_Grid1d.hpp"
#include "ns_io_binary.hpp"

//########################################################################
//########################################################################
/// @class Grid2d
/// @brief Two-dimensional grid class composed of Grid1d x_axis and y_axis.
///
/// @details
/// Grid2d represents a regular 2D grid with uniform or non-uniform spacing along each axis.
/// It combines two independent Grid1d instances (x_axis and y_axis) to form a Cartesian product grid.
///
/// Responsibilities:
/// - Maintain two orthogonal 1D grids (x_axis, y_axis)
/// - Convert between physical coordinates (x, y) and grid indices (ix, iy)
/// - Provide boundary checking and range queries
/// - Support ray tracing via voxel traversal algorithms
/// - Enable grid manipulation (cutting, merging)
///
/// Invariants:
/// - x_axis and y_axis must have valid configurations (nbin > 0, min < max)
/// - Index ranges: 0 <= ix < nbinx, 0 <= iy < nbiny
/// - Coordinate ranges: xmin <= x < xmax, ymin <= y < ymax
///
/// Typical usage:
/// @code
/// Grid1d x_axis("x", 10, 0.0, 10.0, 1.0);
/// Grid1d y_axis("y", 10, 0.0, 10.0, 1.0);
/// Grid2d grid(x_axis, y_axis);
/// int ix = grid.get_ix(5.5);  // Get index for x=5.5
/// int iy = grid.get_iy(7.2);  // Get index for y=7.2
/// if (grid.is_inside(ix, iy)) {
///   double xcnt = grid.get_xcnt(ix);
///   double ycnt = grid.get_ycnt(iy);
/// }
/// @endcode
///
/// @ingroup basicGridClasses
//########################################################################
//########################################################################
class Grid2d {
  private:
    /// @brief Type alias for std::source_location
    using srcloc = std::source_location;
  private:
    /// @brief the name of the Grid2d
    std::string name = "none";

    /// @brief the x-axis of Grid1d
    Grid1d x_axis;

    /// @brief the y-axis of Grid1d
    Grid1d y_axis;

    /// @brief Cached AABB2d for the grid bounding box, avoiding per-ray reconstruction.
    /// @note Updated by rebuild_cached_aabb2d(). Must be called after axis changes.
    AABB2d cached_aabb2d_;

    /// @brief Rebuild cached_aabb2d_ from current axis bounds.
    void rebuild_cached_aabb2d();

  public:
    //======================================================================
    /// @name type_definitions
    ///@{
    using Ixiy = std::array<int,2>;
    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name constructor_desctructor
    ///@{
    
    /// @brief  Default constructor
    /// @details Also rebuilds cached_aabb2d_ so that ray tracing can read
    ///          cached_aabb2d_ instead of constructing a temporary AABB per ray.
    Grid2d() : name("none"), x_axis(), y_axis() { rebuild_cached_aabb2d(); };

    /// @brief  copy constructor
    Grid2d(const Grid2d &org) = default;

    /// @brief  move constructor
    Grid2d(Grid2d &&other) noexcept = default;

    /// @brief  destructor
    virtual ~Grid2d() = default;

    /// @brief constructor using initial values
    Grid2d(const Grid1d &x_axis_in, const Grid1d &y_axis_in);

    ///@} ------------------------------------------------------------------

    //==================================================================
    /// @name operators
    ///@{

    /// @brief Assignment operator
    Grid2d& operator=(const Grid2d& other) = default;

    /// @brief Inequality operator
    bool operator!=(const Grid2d& other) const;

    /// @brief Equality operator using non-equal operator
    bool operator==(const Grid2d& other) const;
    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name check_functions
    ///@{
  
    /// @brief Check if two Ixiy coordinates are adjacent
    /// @param[in] a First grid coordinate
    /// @param[in] b Second grid coordinate
    /// @param[in] distance_threshold Manhattan distance threshold (default: 1)
    /// @return true if Manhattan distance is <= threshold and > 0, false otherwise
    /// @note Returns false if a == b (distance = 0)
    bool is_adjacent(const Ixiy& a, const Ixiy& b, const int distance_threshold = 1) const;

    /// @brief Check if ix is within valid range [0, nbinx)
    /// @param[in] ix Grid index in x direction
    /// @return true if 0 <= ix < nbinx, false otherwise
    bool is_ix_inside( const int ix ) const { return x_axis.is_index_inside(ix); };

    /// @brief Check if iy is within valid range [0, nbiny)
    /// @param[in] iy Grid index in y direction
    /// @return true if 0 <= iy < nbiny, false otherwise
    bool is_iy_inside( const int iy ) const { return y_axis.is_index_inside(iy); };

    /// @brief Check if both ix and iy are within valid ranges
    /// @param[in] ix Grid index in x direction
    /// @param[in] iy Grid index in y direction
    /// @return true if both indices are inside, false otherwise
    bool is_inside( const int ix, const int iy ) const;

    /// @brief Check if x coordinate is within valid range [xmin, xmax)
    /// @param[in] x Physical coordinate in x direction
    /// @return true if xmin <= x < xmax, false otherwise
    bool is_x_inside( const double x ) const { return x_axis.is_value_inside(x); };

    /// @brief Check if y coordinate is within valid range [ymin, ymax)
    /// @param[in] y Physical coordinate in y direction
    /// @return true if ymin <= y < ymax, false otherwise
    bool is_y_inside( const double y ) const { return y_axis.is_value_inside(y); };

    /// @brief check if ix is inside the valid range, throw error if not
    /// @param ix index to check
    /// @param loc source location (automatically captured at call site)
    /// @throw std::runtime_error if ix is out of range
    void check_ix_inside(const int ix,
      const srcloc loc = srcloc::current()) const;

    /// @brief check if iy is inside the valid range, throw error if not
    /// @param iy index to check
    /// @param loc source location (automatically captured at call site)
    /// @throw std::runtime_error if iy is out of range
    void check_iy_inside(const int iy,
      const srcloc loc = srcloc::current()) const;

    /// @brief return true if index is above lower limit. if not, return false. 
    bool is_ix_above_lower_limit( const int ix ) const { return x_axis.is_above_lower_limit(ix); };
    
    /// @brief return true if index is above lower limit. if not, return false. 
    bool is_iy_above_lower_limit( const int iy ) const { return y_axis.is_above_lower_limit(iy); };
    
    /// @brief return true if index is below upper limit. if not, return false.
    bool is_ix_below_upper_limit( const int ix ) const { return x_axis.is_below_upper_limit(ix); };
    
    /// @brief return true if index is below upper limit. if not, return false.
    bool is_iy_below_upper_limit( const int iy ) const { return y_axis.is_below_upper_limit(iy); };

    /// @brief Check if physical coordinates (x, y) are inside the grid
    /// @param[in] x Physical coordinate in x direction
    /// @param[in] y Physical coordinate in y direction
    /// @return true if xmin <= x < xmax AND ymin <= y < ymax, false otherwise
    bool is_inside( const double x, const double y ) const;

    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name getter_functions
    ///@{
    
    /// @brief get name of this instance
    std::string get_name() const { return name; };

    /// @brief get reference of x_axis
    const Grid1d& get_x_axis() const { return x_axis; };
    /// @brief get reference of y_axis
    const Grid1d& get_y_axis() const { return y_axis; };

    /// @brief get nbinx
    int get_nbinx() const { return x_axis.get_nbin(); };

    /// @brief get nbiny
    int get_nbiny() const { return y_axis.get_nbin(); };

    /// @brief get xmin
    double get_xmin() const { return x_axis.get_min(); };

    /// @brief get xmax
    double get_xmax() const { return x_axis.get_max(); };

    /// @brief get ymin
    double get_ymin() const { return y_axis.get_min(); };

    /// @brief get ymax
    double get_ymax() const { return y_axis.get_max(); };

    /// @brief get x_interval
    double get_x_interval() const { return x_axis.get_interval(); };
    
    /// @brief get y_interval
    double get_y_interval() const { return y_axis.get_interval(); };

    /// @brief  return the index of x_axis
    /// @details in the case of out of range, return -1;
    /// @param value
    /// @return x_axis.get_index(value);
    int get_ix( const double value ) const {
      return x_axis.get_index(value);
    };

    /// @brief  return the index of y_axis
    /// @details in the case of out of range, return -1;
    /// @param value
    /// @return y_axis.get_index(value);
    int get_iy( const double value ) const {
      return y_axis.get_index(value);
    };

    /// @brief  return the center value in the bin of the index on x_axis
    /// @param index
    /// @return x_axis.get_center_value(index);
    double get_xcnt( const int index ) const {
      return x_axis.get_center_value(index);
    };

    /// @brief  return the center value in the bin of the value on x_axis
    /// @param value
    /// @return index = get_ix(value); return this->get_xcnt(index);
    double get_xcnt( const double value ) const;

    /// @brief  return the center value in the bin of the index on y_axis
    /// @param index
    /// @return y_axis.get_center_value(index);
    double get_ycnt( const int index ) const {
      return y_axis.get_center_value(index);
    };

    /// @brief return the center value in the bin of the value on y_axis
    /// @param value
    /// @return index = get_iy(value); return this->get_ycnt(index);
    double get_ycnt( const double value ) const;

    /// @brief return the lower value in the bin of the index on x_axis
    /// @param index
    /// @return return x_axis.get_lower_value(index);
    double get_xlow( const int index ) const {
      return x_axis.get_lower_value(index);
    };

    /// @brief return the lower value in the bin of the index on y_axis
    /// @param index
    /// @return return y_axis.get_lower_value(index);
    double get_ylow( const int index ) const {
      return y_axis.get_lower_value(index);
    };

    /// @brief return the upper value in the bin of the index on x_axis
    /// @param index
    /// @return return x_axis.get_upper_value(index);
    double get_xup( const int index ) const {
      return x_axis.get_upper_value(index);
    };

    /// @brief return the upper value in the bin of the index on y_axis
    /// @param index
    /// @return return y_axis.get_upper_value(index);
    double get_yup( const int index ) const {
      return y_axis.get_upper_value(index);
    };

    /// @brief return the lower value in the bin of the index on x_axis
    /// @param value
    /// @return index = get_ix(value); return get_xlow(index);
    double get_xlow( const double value ) const;

    /// @brief return the lower value in the bin of the index on y_axis
    /// @param value
    /// @return index = get_iy(value); return get_ylow(index);
    double get_ylow( const double value ) const;

    /// @brief return the upper value in the bin of the index on x_axis
    /// @param value
    /// @return index = get_ix(value); return get_xup(index);
    double get_xup( const double value ) const;

    /// @brief return the upper value in the bin of the index on y_axis
    /// @param value
    /// @return index = get_ix(value); return get_xup(index);
    double get_yup( const double value ) const;

    /// @brief get get_xup() - get_xlow();
    double get_xlen( const int ix ) const;

    /// @brief get get_yup() - get_ylow();
    double get_ylen( const int iy ) const;

    /// @brief return xlower, xuppeer in index ix, and ylower, yupper in index iy
    /// @param ix 
    /// @param iy 
    /// @return {xlower, xupper, ylower, yupper}
    /// @details usage const auto bounds get_xy_lower_upper(ix,iy); xmin = bounds[0]; xmax = bounds[1]; ymin = bounds[2]; ymax = bounds[3]; \n 
    /// or const auto [xmin,xmax,ymin,ymax] = get_xy_lower_upper(ix,iy); \n
    /// @note DO NOT LIKE const auto& [xmin,xmax,ymin,ymax] = get_xy_lower_upper(ix,iy);
    std::array<double,4>
      get_xy_lower_upper( const int ix, const int iy ) const; 

    /// @brief output out_g2_header
    void out_g2_header( FILE *fout ) const;

    /// @brief Get vector of (ix, iy) pairs within specified physical coordinate range
    /// @param[in] xmin_in Minimum x coordinate
    /// @param[in] xmax_in Maximum x coordinate
    /// @param[in] ymin_in Minimum y coordinate
    /// @param[in] ymax_in Maximum y coordinate
    /// @return Vector of Ixiy coordinates within [xmin_in, xmax_in] x [ymin_in, ymax_in]
    /// @note Small epsilon is added/subtracted internally to handle boundary cases
    std::vector< Ixiy > get_vec_ixiy(
      const double xmin_in, const double xmax_in
    , const double ymin_in, const double ymax_in ) const;

    /// @brief Get vector of (ix, iy) pairs within a circular region
    /// @param[in] x_center Circle center x coordinate
    /// @param[in] y_center Circle center y coordinate
    /// @param[in] radius Circle radius
    /// @return Vector of Ixiy coordinates where cell centers satisfy (x-x_center)^2 + (y-y_center)^2 <= radius^2
    /// @note Uses cell center positions for distance calculation
    /// @note Boundary cells (dist^2 == radius^2) are included
    std::vector< Ixiy > get_vec_ixiy_in_circle(
      const double x_center, const double y_center
    , const double radius ) const;

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
    std::vector<Ixiy> get_vec_ixiy(const std::vector<std::vector<Ixiy>> &vec_vec_ixiy ) const;

    /// @brief Return the minimal bounding rectangle as a 2D grid containing all given Ixiy coordinates
    ///
    /// @details
    /// From the input std::vector<Ixiy> coordinates,
    /// extract min/max values in x and y directions, and return a 2D vector
    /// filled with all Ixiy coordinates (std::array<int, 2>) in that range.
    /// The returned vector has a nested structure with y as outer index and x as inner index.
    ///
    /// Ixiy is defined as:
    /// @code
    /// using Ixiy = std::array<int, 2>;
    /// Ixiy[0] = ix, Ixiy[1] = iy
    /// @endcode
    ///
    /// @param[in] vec_ixiy Vector of Ixiy coordinates
    /// @return std::vector<std::vector<Ixiy>>
    ///         Format: vec_vec_ixiy[iy][ix]
    ///         Range: (ixmin ≤ ix ≤ ixmax, iymin ≤ iy ≤ iymax)
    ///
    /// @note
    /// Usage example:
    /// @code
    /// std::vector<Ixiy> input = { {1, 1}, {2, 3} };
    /// auto vec_vec_ixiy = obj.get_vec_vec_ixiy(input);
    /// vec_vec_ixiy[0][0] = {1, 1},
    /// vec_vec_ixiy[1][0] = {2, 1},
    /// vec_vec_ixiy[0][1] = {1, 2},
    /// vec_vec_ixiy[1][1] = {2, 2},
    /// vec_vec_ixiy[0][2] = {1, 3},
    /// vec_vec_ixiy[1][2] = {2, 3}, returning a 2D array of Ixiy coordinates.
    /// @endcode
    std::vector<std::vector<Ixiy>> 
      get_vec_vec_ixiy( const std::vector<Ixiy> &vec_ixiy ) const;

    /// @brief get ixmin_ixmax_iymin_iymax from vec_ixiy
    std::array<int,4> get_ixmin_ixmax_iymin_iymax(
      const std::vector<Ixiy> &vec_ixiy ) const;

    /// @brief get length num of ix and iy from vec_ixiy
    Ixiy get_ixiylen(
      const std::vector<Ixiy> &vec_ixiy ) const;

    /// @brief Create a sub-grid clipped by the provided physical coordinate bounds
    /// @param[in] x_lower Lower x boundary
    /// @param[in] x_upper Upper x boundary
    /// @param[in] y_lower Lower y boundary
    /// @param[in] y_upper Upper y boundary
    /// @param[in] x_eps Epsilon for x boundary tolerance (default: 1.0e-6)
    /// @param[in] y_eps Epsilon for y boundary tolerance (default: 1.0e-6)
    /// @return New Grid2d instance containing only cells within specified bounds
    /// @throws std::runtime_error if x_eps < 0 or y_eps < 0
    Grid2d cut(const double x_lower, const double x_upper
             , const double y_lower, const double y_upper
             , const double x_eps = 1.0e-6, const double y_eps = 1.0e-6) const;

    /// @brief get tuple of xcnt,ycnt from tp_ixiy
    std::array<double,2> get_xcnt_ycnt( const Ixiy ixiy );

    /// @brief get the AABB index of the Grid2d
    /// @param x0 position x0
    /// @param y0 position y0
    /// @param vx unit direction vector
    /// @param vy unit direction vector
    /// @return vector of hit box index (ix, iy) 
    /// @ingroup ray_tracing_functions
    std::vector<Ixiy> get_hit_boxes_index_old(
      const double x0, const double y0, const double vx, const double vy) const;

    // get_hit_boxes_index declarations moved to path_length_functions section below.

    /// @brief Get the std::vector<Ixiy> in the original non-merged Grid2d corresponding to a cell in the merged Grid2d.
    std::vector<Ixiy>
      get_non_merged_vec_ixiy( const Grid2d &g2d_merged, const int ix, const int iy ) const;

    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name setter_functions
    ///@{
    
    /// @brief  set x_axis and y_axis from std::vector\<Eigen::Vector3d\>
    /** @details  assuming that data order is \n
        x0 y0 z00 \n
        x0 y1 z01 \n
        , ... ,   \n
        x0 y2 z02 \n
        x1 y0 z10 \n
        x1 y1 z11 \n
        x1 y2 z12 \n
        , ... ,   \n
        x2 y0 z20 \n
        x2 y1 z21 \n
        x2 y2 z22 \n
        , ... ,   \n
        and x0 < x1 < x2, y0 < y1 < y2 \n
    **/
    /// @param std::vector\<Eigen::Vector3d\> &vec_v3
    /// @param bintype of x_axis
    /// @param bintype of y_axis
    /// @param tolerance_ratio of same grid interval
    /// @note Unused as of 2024-12-03 12:23:28
    // void set_xy_axis_from_vec_v3( const std::vector<Eigen::Vector3d> &vec_v3
    //       , const bool tf_shift_x, const bool tf_shift_y
    //       , const double tolerance_ratio );
    
    /// @brief set x_axis and y_axis from std::vector\<std::tuple\<double,double,double\>\>
    /** @details  assuming that data order is \n
        x0 y0 z00 \n
        x0 y1 z01 \n
        , ... ,   \n
        x0 y2 z02 \n
        x1 y0 z10 \n
        x1 y1 z11 \n
        x1 y2 z12 \n
        , ... ,   \n
        x2 y0 z20 \n
        x2 y1 z21 \n
        x2 y2 z22 \n
        , ... ,   \n
        and x0 < x1 < x2, y0 < y1 < y2 \n
    **/
    /// @param std::vector<std::array<double,3>> &vec_xyz
    /// @param bintype of x_axis
    /// @param bintype of y_axis
    /// @param tolerance_ratio of same grid interval
    /// @return bool if interval of x and y is not same, return false
    bool set_xy_axis_from_vec_xyz(
        const std::vector<std::array<double,3>> &vec_xyz
      , const bool tf_shift_x, const bool tf_shift_y
      , const double tolerance_ratio );

    /// @brief x_axis.set_name(name_in);
    /// @details do x_axis.set_name(name_in);
    /// @param name_in 
    void set_x_axis_name( const std::string &name_in ){ x_axis.set_name(name_in); };

    /// @brief y_axis.set_name(name_in);
    /// @details do y_axis.set_name(name_in);
    /// @param name_in 
    void set_y_axis_name( const std::string &name_in ){ y_axis.set_name(name_in); };

    /// @brief x_axis = x_axis_in;
    /// @details Also rebuilds cached_aabb2d_. See Grid3d::set_x_axis().
    /// @param x_axis_in
    void set_x_axis( const Grid1d &x_axis_in ){ x_axis = x_axis_in; rebuild_cached_aabb2d(); };

    /// @brief y_axis = y_axis_in;
    /// @details Also rebuilds cached_aabb2d_. See Grid3d::set_y_axis().
    /// @param y_axis_in
    void set_y_axis( const Grid1d &y_axis_in ){ y_axis = y_axis_in; rebuild_cached_aabb2d(); };

    /// @brief Grid1d::x_axis.assign(x_name_in,nbinx_in,xmin_in,xmax_in,xpit_in,tf_shift_x);
    /// @param x_name_in 
    /// @param nbinx_in 
    /// @param xmin_in 
    /// @param xmax_in 
    /// @param xpit_in 
    /// @param tf_shift_x 
    void set_x_axis( const std::string x_name_in, const int nbinx_in
    , const double xmin_in, const double xmax_in, const double xpit_in);

    /// @brief Grid1d::y_axis.assign(y_name_in,nbiny_in,ymin_in,ymax_in,ypit_in,tf_shift_y);
    /// @param y_name_in
    /// @param nbiny_in
    /// @param ymin_in
    /// @param ymax_in
    /// @param ypit_in
    /// @param tf_shift_y
    void set_y_axis( const std::string y_name_in, const int nbiny_in
    , const double ymin_in, const double ymax_in, const double ypit_in );
    ///@} ------------------------------------------------------------------
    
    /// @brief Get the cached AABB2d of the grid bounding box.
    /// @return Const reference to the cached AABB2d.
    const AABB2d& get_cached_aabb2d() const { return cached_aabb2d_; }

    //======================================================================
    /// @name path_length_functions
    ///@{

    /// @brief Parameter structure for passing arguments to get_hit_boxes_index
    struct VoxelTraversalParameters2d {
      double x_start, y_start;
      int manhattan_distance;
      double step_x, step_y;
      double t_delta_x, t_delta_y;
      double t_max_x0, t_max_y0;
      double vx_sign, vy_sign;
    };

    /// @brief Build DDA traversal parameters from a valid ray-AABB intersection result.
    /// @param ray2d Ray2d with origin and direction.
    /// @param eps Small scalar to avoid missing boundary cells (relative to grid spacing).
    /// @param tmin Intersection entry parameter from ray2d.is_intersect().
    /// @param tmax Intersection exit parameter from ray2d.is_intersect().
    /// @return Filled VoxelTraversalParameters2d, ready for get_hit_index_by_ray_tracing_algorithm.
    /// @note Extracts the common ~60 lines shared by get_hit_boxes_index() variants.
    VoxelTraversalParameters2d build_traversal_params_2d(
        const Ray2d &ray2d, double eps, double tmin, double tmax) const;

    /// @brief Trace a 2D ray through the grid to find hit cell indices.
    /// @param ray2d Ray2d with origin and direction.
    /// @param eps Small scalar to avoid missing boundary cells (relative to grid spacing).
    /// @return Vector of hit box indices (ix, iy), ordered along the ray from near to far.
    /// @note Uses cached_aabb2d_ to skip per-ray AABB construction (optimization #5).
    /// @note Thread-safe: Yes.
    /// @ingroup ray_tracing_functions
    std::vector<Ixiy>
      get_hit_boxes_index( const Ray2d &ray2d, const double eps = 1.0e-6 ) const;

    /// @brief Trace a 2D ray with beam length warning.
    /// @param ray2d Ray2d with origin and direction.
    /// @param BL_max Maximum beam length threshold [meters]. If <= 0, warning is disabled.
    /// @param eps Small scalar to avoid missing boundary cells.
    /// @return Vector of hit box indices, ordered along the ray.
    /// @note Uses cached_aabb2d_ and a single intersection test (optimization #4/#5).
    /// @note Thread-safe: Yes.
    /// @ingroup ray_tracing_functions
    std::vector<Ixiy>
      get_hit_boxes_index( const Ray2d &ray2d, const double BL_max, const double eps ) const;

    /// @brief [Legacy] Trace a 2D ray — constructs a local AABB each call.
    /// @deprecated Retained for regression testing; prefer get_hit_boxes_index.
    std::vector<Ixiy>
      get_hit_boxes_index_legacy( const Ray2d &ray2d, const double eps ) const;

    /// @brief [Legacy] Trace with BL_max — performs double ray-AABB intersection.
    /// @deprecated Retained for regression testing; prefer get_hit_boxes_index.
    std::vector<Ixiy>
      get_hit_boxes_index_legacy( const Ray2d &ray2d, const double BL_max, const double eps ) const;

    /// @brief Traverse 2D ray through grid, calling visitor for each hit cell.
    /// @tparam Visitor Callable with signature void(const Ixiy&).
    /// @param ray2d Ray2d with origin and direction.
    /// @param eps Small scalar to avoid missing boundary cells.
    /// @param visitor Callback invoked for each hit cell index.
    /// @note Eliminates per-ray std::vector heap allocation (optimization #3).
    /// @note Uses cached_aabb2d_ (optimization #5).
    /// @note Thread-safe: Yes (no internal mutation).
    /// @note Implementation in header — required for template instantiation.
    template<typename Visitor>
    void traverse_ray_2d(
        const Ray2d &ray2d, double eps, Visitor&& visitor) const;

    /// @brief helper function of get_hit_boxes_index( ray2d )
    /// @return the index of hit 2d boxes from near to far (returns voxel indices in order of collision from nearest to farthest).
    /// @ingroup ray_tracing_functions
    std::vector<Ixiy>
      get_hit_index_by_ray_tracing_algorithm(
        const Grid2d::VoxelTraversalParameters2d &prm) const;
    ///@} ------------------------------------------------------------------

    /// @brief for bug check, output ix iy xlow xup ylow yup
    void out_grid(const std::filesystem::path& pathout) const;

    //======================================================================
    /// @name binary_io_functions
    ///@{
    
    /// @brief Write Grid2d to binary stream (canonical coordinates via Grid1d::save).
    void save( std::ofstream &ofs ) const;

    /// @brief Read Grid2d from binary stream (pass-through via Grid1d::load).
    /// @param[in] tolerance_ratio  Tolerance for Grid1d consistency checks
    void load( std::ifstream &ifs, double tolerance_ratio = 1.0e-4 );

    ///@} ------------------------------------------------------------------

};

//----------------------------------------------------------------------
// Template implementation: traverse_ray_2d
// Must reside in header for template instantiation.
//----------------------------------------------------------------------
template<typename Visitor>
void Grid2d::traverse_ray_2d(
    const Ray2d &ray2d, double eps, Visitor&& visitor) const
{
  auto [tf_hit_volume, tmin, tmax] = ray2d.is_intersect(cached_aabb2d_);
  if (!tf_hit_volume) return;
  if (tmin < 0 && tmax < 0) return;

  const auto prm = build_traversal_params_2d(ray2d, eps, tmin, tmax);

  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();

  int ix = get_ix(prm.x_start);
  int iy = get_iy(prm.y_start);

  double t_max_x = prm.t_max_x0;
  double t_max_y = prm.t_max_y0;

  visitor(Ixiy{ix, iy});

  for (int t = 0; t <= prm.manhattan_distance; t++) {
    if (fabs(t_max_x) < fabs(t_max_y)) {
      t_max_x += prm.t_delta_x;
      ix += (int)prm.vx_sign;
      if (ix < 0 || ix > nbinx - 1) return;
    } else {
      t_max_y += prm.t_delta_y;
      iy += (int)prm.vy_sign;
      if (iy < 0 || iy > nbiny - 1) return;
    }
    visitor(Ixiy{ix, iy});
  }
}
