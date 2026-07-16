/// @file cls_Grid3d.hpp
/// @brief 3D uniform grid defined by three orthogonal Grid1d axes.
///
/// @details This file defines Grid3d, a lightweight geometry class that stores
///          three Grid1d axes and provides index/coordinate conversions, AABB
///          queries, and ray traversal utilities.
///
///          Typical workflow:
///          1. Build Grid1d axes for x/y/z (min, max, interval).
///          2. Construct Grid3d from the axes.
///          3. Query voxel indices or bounds for spatial points.
///          4. Trace rays through voxels with get_hit_voxel_index.
///          5. Export grid ranges with out_grid_2d/out_grid_3d.
///
///          Coordinate system:
///          - Right-handed, z-axis is vertical (positive upward)
///          - Units follow the Grid1d axis values (typically meters)
///
///          Memory layout:
///          - Grid3d stores only axis metadata (no voxel payload)
///          - out_grid_2d iterates iy -> ix
///          - out_grid_3d iterates iz -> iy -> ix
///
///          Thread safety:
///          - const methods are thread-safe
///          - mutating methods are not thread-safe
///
///          I/O formats:
///          - out_grid_2d/out_grid_3d write ASCII columns to file
#pragma once

#include <array>
#include <filesystem>
#include <map>
#include <unordered_map>
#include <fstream>
#include <iostream>
#include <sstream> // istringstream
#include <string>
#include <cstdio>
#include <stdexcept>
#include <functional>  //for sorting
#include <algorithm>//for sorting
#include <vector>

#include "spdlog_pch.hpp"
#include "ns_mylogger.hpp"
#include <Eigen/Dense>

#include "cls_AABB.hpp"
#include "cls_Ray.hpp"
#include "cls_Grid1d.hpp"
#include "ns_type_definitions.hpp"
//======================================================================
/// @name check macro for Grid3d
///@{

/// @brief Check ix against x_axis and throw on out-of-range.
#define G3_CHECK_IX_INSIDE(instance,index) (instance).check_ix_inside(index, __FILE__, __PRETTY_FUNCTION__, __LINE__)
/// @brief Check iy against y_axis and throw on out-of-range.
#define G3_CHECK_IY_INSIDE(instance,index) (instance).check_iy_inside(index, __FILE__, __PRETTY_FUNCTION__, __LINE__)
/// @brief Check iz against z_axis and throw on out-of-range.
#define G3_CHECK_IZ_INSIDE(instance,index) (instance).check_iz_inside(index, __FILE__, __PRETTY_FUNCTION__, __LINE__)
///@} ------------------------------------------------------------------

//#######################################################################################
//#######################################################################################
/// @class Grid3d
/// @brief 3D uniform grid composed of three Grid1d axes.
///
/// @details Grid3d provides index/coordinate conversions and AABB queries
///          without storing voxel payloads. It is commonly used as a geometry
///          base class for voxelized data containers (e.g., Grid3dVoxel).
///
///          Primary use cases:
///          - Map continuous coordinates to discrete voxel indices
///          - Query voxel AABB bounds for spatial computations
///          - Trace rays through voxels using a 3D DDA traversal
///
///          Example usage:
///          @code
///          Grid1d gx("x", 100, 0.0, 1000.0, 10.0);
///          Grid1d gy("y", 100, 0.0, 1000.0, 10.0);
///          Grid1d gz("z", 50, 0.0, 500.0, 10.0);
///          Grid3d grid(gx, gy, gz);
///          Grid3d::Ixiyiz idx = grid.get_ixiyiz(12.3, 45.6, 78.9);
///          AABB3d cell = grid.get_AABB3d(idx[0], idx[1], idx[2]);
///          @endcode
///
/// @ingroup basicGridClasses
//#######################################################################################
//#######################################################################################
class Grid3d {
  private:
    /// @brief name of this instance
    std::string name;

    /// @brief Grid1d x_axis
    Grid1d x_axis;
    /// @brief Grid1d y_axis
    Grid1d y_axis;
    /// @brief Grid1d z_axis
    Grid1d z_axis;

    /// @brief Cached AABB3d for the entire grid.
    /// @details get_hit_voxel_index() previously constructed a temporary AABB3d
    ///          from 6 getter calls on every ray. Since the grid bounds never
    ///          change after construction, caching eliminates that per-ray cost.
    /// @note Updated by rebuild_cached_aabb3d(). Must be called after axis changes.
    AABB3d cached_aabb3d_;

    /// @brief Rebuild cached_aabb3d_ from current axis bounds.
    void rebuild_cached_aabb3d();

  public:
    //======================================================================
    /// @name type_definitions
    ///@{
    /// @brief Index triplet (ix, iy, iz).
    using Ixiyiz = std::array<int, 3>;

    /// @brief Unique voxel index type.
    using Uqiv = int;

    /// @brief Sentinel value for "not found" (ix, iy, iz) = (-1, -1, -1).
    static constexpr Ixiyiz IxiyizNotFound = { -1, -1, -1 };

    /// @brief Sentinel value for "not found" unique index = -1.
    static constexpr Uqiv UqivNotFound = -1;

    /// @brief Hash function for Ixiyiz.
    struct IxiyizHash {
      std::size_t operator()(const Ixiyiz& key) const noexcept {
        std::size_t h1 = std::hash<int>{}(key[0]);
        std::size_t h2 = std::hash<int>{}(key[1]);
        std::size_t h3 = std::hash<int>{}(key[2]);
        return h1 ^ (h2 << 1) ^ (h3 << 2);
      }
    };

    /// @brief Equality comparison for Ixiyiz.
    struct IxiyizEqual {
      bool operator()(const Ixiyiz& a, const Ixiyiz& b) const noexcept {
        return a == b;
      }
    };

    /// @brief Hash map from unique index to (ix, iy, iz).
    using UmpUqivIxiyiz = std::unordered_map<Uqiv, Ixiyiz>;

    /// @brief Hash map from (ix, iy, iz) to unique index.
    using UmpIxiyizUqiv = std::unordered_map<Ixiyiz, Uqiv, IxiyizHash, IxiyizEqual>;
    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name constructor_destructor
    ///@{

    /// @brief Default constructor (name = "none").
    /// @details Calls rebuild_cached_aabb3d() so that ray tracing can read
    ///          cached_aabb3d_ instead of constructing a temporary AABB per ray.
    /// @note Thread-safe: No (construction).
    Grid3d() :name("none"), x_axis(), y_axis(), z_axis() { rebuild_cached_aabb3d(); };

    /// @brief Copy constructor.
    /// @param org Source Grid3d.
    /// @note Thread-safe: Yes.
    Grid3d(const Grid3d &org) = default;

    /// @brief Move constructor.
    /// @param other Source Grid3d.
    /// @note Thread-safe: Yes.
    Grid3d(Grid3d &&other) noexcept = default;

    /// @brief Destructor.
    /// @note Thread-safe: Yes.
    virtual ~Grid3d() = default;

    /// @brief Construct from x/y/z Grid1d axes.
    /// @param x_axis_in X-axis Grid1d.
    /// @param y_axis_in Y-axis Grid1d.
    /// @param z_axis_in Z-axis Grid1d.
    /// @note Thread-safe: No (construction).
    Grid3d( const Grid1d &x_axis_in, const Grid1d &y_axis_in, const Grid1d &z_axis_in );

    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name operator
    ///@{
    
    /// @brief Copy assignment operator.
    /// @param other Source Grid3d.
    /// @return Reference to this instance.
    /// @note Thread-safe: No.
    Grid3d& operator=(const Grid3d& other);

    /// @brief Inequality operator.
    /// @param other Grid3d to compare.
    /// @return True if axes differ (name is ignored).
    /// @note Thread-safe: Yes.
    bool operator!=(const Grid3d& other) const;

    /// @brief Equality operator (defined via operator!=).
    /// @param other Grid3d to compare.
    /// @return True if axes match (name is ignored).
    /// @note Thread-safe: Yes.
    bool operator==(const Grid3d& other) const {
      return !(*this != other);
    };
    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name checker_functions
    /// @{

    /// @brief Check whether ix is inside the valid range.
    /// @param ix Grid index in x-direction.
    /// @return True if ix is within [0, nbinx-1].
    /// @note Thread-safe: Yes. Time complexity: O(1).
    bool is_ix_inside( const int ix ) const { return x_axis.is_index_inside(ix); };
    /// @brief Check whether iy is inside the valid range.
    /// @param iy Grid index in y-direction.
    /// @return True if iy is within [0, nbiny-1].
    /// @note Thread-safe: Yes. Time complexity: O(1).
    bool is_iy_inside( const int iy ) const { return y_axis.is_index_inside(iy); };
    /// @brief Check whether iz is inside the valid range.
    /// @param iz Grid index in z-direction.
    /// @return True if iz is within [0, nbinz-1].
    /// @note Thread-safe: Yes. Time complexity: O(1).
    bool is_iz_inside( const int iz ) const { return z_axis.is_index_inside(iz); };

    /// @brief Check whether ix is above the lower limit.
    /// @param ix Grid index in x-direction.
    /// @return True if ix is at or above the lower limit.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    bool is_ix_above_lower_limit( const int ix ) const { return x_axis.is_above_lower_limit(ix); };

    /// @brief Check whether iy is above the lower limit.
    /// @param iy Grid index in y-direction.
    /// @return True if iy is at or above the lower limit.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    bool is_iy_above_lower_limit( const int iy ) const { return y_axis.is_above_lower_limit(iy); };

    /// @brief Check whether iz is above the lower limit.
    /// @param iz Grid index in z-direction.
    /// @return True if iz is at or above the lower limit.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    bool is_iz_above_lower_limit( const int iz ) const { return z_axis.is_above_lower_limit(iz); };

    /// @brief Check whether ix is below the upper limit.
    /// @param ix Grid index in x-direction.
    /// @return True if ix is at or below the upper limit.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    bool is_ix_below_upper_limit( const int ix ) const { return x_axis.is_below_upper_limit(ix); };

    /// @brief Check whether iy is below the upper limit.
    /// @param iy Grid index in y-direction.
    /// @return True if iy is at or below the upper limit.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    bool is_iy_below_upper_limit( const int iy ) const { return y_axis.is_below_upper_limit(iy); };

    /// @brief Check whether iz is below the upper limit.
    /// @param iz Grid index in z-direction.
    /// @return True if iz is at or below the upper limit.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    bool is_iz_below_upper_limit( const int iz ) const { return z_axis.is_below_upper_limit(iz); };

    /// @brief Throw if ix is outside the valid range.
    /// @param ix Grid index in x-direction.
    /// @throws std::runtime_error if ix is out of range.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    void check_ix_inside( const int ix ) const;

    /// @brief Throw if iy is outside the valid range.
    /// @param iy Grid index in y-direction.
    /// @throws std::runtime_error if iy is out of range.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    void check_iy_inside( const int iy ) const;

    /// @brief Throw if iz is outside the valid range.
    /// @param iz Grid index in z-direction.
    /// @throws std::runtime_error if iz is out of range.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    void check_iz_inside( const int iz ) const;

    /// @brief Throw if ix is outside the valid range (macro helper).
    /// @param ix Grid index in x-direction.
    /// @param filename Source file (macro-provided).
    /// @param function_name Function name (macro-provided).
    /// @param line Source line (macro-provided).
    /// @throws std::runtime_error if ix is out of range.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    void check_ix_inside( const int ix, const std::string &filename
      , const std::string &function_name, const int line ) const;

    /// @brief Throw if iy is outside the valid range (macro helper).
    /// @param iy Grid index in y-direction.
    /// @param filename Source file (macro-provided).
    /// @param function_name Function name (macro-provided).
    /// @param line Source line (macro-provided).
    /// @throws std::runtime_error if iy is out of range.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    void check_iy_inside( const int iy, const std::string &filename
      , const std::string &function_name, const int line ) const;

    /// @brief Throw if iz is outside the valid range (macro helper).
    /// @param iz Grid index in z-direction.
    /// @param filename Source file (macro-provided).
    /// @param function_name Function name (macro-provided).
    /// @param line Source line (macro-provided).
    /// @throws std::runtime_error if iz is out of range.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    void check_iz_inside( const int iz, const std::string &filename
      , const std::string &function_name, const int line ) const;
    ///@} ------------------------------------------------------------------

    //========================================================
    /// @name getter_Grid3d
    /// @{
    
    /// @brief Get name of this instance.
    /// @return Instance name string.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    std::string get_name() const { return name; };

    /// @brief Get reference to the x-axis.
    /// @return Const reference to Grid1d x-axis.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    const Grid1d& get_x_axis() const { return x_axis; };

    /// @brief Get reference to the y-axis.
    /// @return Const reference to Grid1d y-axis.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    const Grid1d& get_y_axis() const { return y_axis; };

    /// @brief Get reference to the z-axis.
    /// @return Const reference to Grid1d z-axis.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    const Grid1d& get_z_axis() const { return z_axis; };

    /// @brief Get number of bins along x.
    /// @return Number of bins in x-axis.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    int get_nbinx() const { return x_axis.get_nbin(); }

    /// @brief Get number of bins along y.
    /// @return Number of bins in y-axis.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    int get_nbiny() const { return y_axis.get_nbin(); }

    /// @brief Get number of bins along z.
    /// @return Number of bins in z-axis.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    int get_nbinz() const { return z_axis.get_nbin(); }

    /// @brief Get bin interval along x.
    /// @return Bin interval in x-axis units.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    double get_x_interval() const { return x_axis.get_interval(); }

    /// @brief Get bin interval along y.
    /// @return Bin interval in y-axis units.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    double get_y_interval() const { return y_axis.get_interval(); }

    /// @brief Get bin interval along z.
    /// @return Bin interval in z-axis units.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    double get_z_interval() const { return z_axis.get_interval(); }

    /// @brief Get minimum x value.
    /// @return Minimum x-axis coordinate [axis units].
    /// @note Thread-safe: Yes. Time complexity: O(1).
    double get_xmin() const { return x_axis.get_min(); };

    /// @brief Get maximum x value.
    /// @return Maximum x-axis coordinate [axis units].
    /// @note Thread-safe: Yes. Time complexity: O(1).
    double get_xmax() const { return x_axis.get_max(); };

    /// @brief Get minimum y value.
    /// @return Minimum y-axis coordinate [axis units].
    /// @note Thread-safe: Yes. Time complexity: O(1).
    double get_ymin() const { return y_axis.get_min(); };

    /// @brief Get maximum y value.
    /// @return Maximum y-axis coordinate [axis units].
    /// @note Thread-safe: Yes. Time complexity: O(1).
    double get_ymax() const { return y_axis.get_max(); };

    /// @brief Get minimum z value.
    /// @return Minimum z-axis coordinate [axis units].
    /// @note Thread-safe: Yes. Time complexity: O(1).
    double get_zmin() const { return z_axis.get_min(); };

    /// @brief Get maximum z value.
    /// @return Maximum z-axis coordinate [axis units].
    /// @note Thread-safe: Yes. Time complexity: O(1).
    double get_zmax() const { return z_axis.get_max(); };

    /// @brief Convert coordinates to (ix, iy, iz).
    /// @param x_in X-coordinate [axis units].
    /// @param y_in Y-coordinate [axis units].
    /// @param z_in Z-coordinate [axis units].
    /// @return Index triplet (ix, iy, iz).
    /// @throws std::runtime_error if any coordinate is outside the grid.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    Ixiyiz get_ixiyiz(
      const double x_in, const double y_in, const double z_in ) const;

    /// @brief Get AABB minimum corner for voxel (ix, iy, iz).
    /// @param ix Grid index in x-direction.
    /// @param iy Grid index in y-direction.
    /// @param iz Grid index in z-direction.
    /// @return AABB minimum corner (xlow, ylow, zlow) [axis units].
    /// @throws std::runtime_error if any index is out of range.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    Eigen::Vector3d get_v3_AABB_min( const int ix, const int iy, const int iz ) const;

    /// @brief Get AABB maximum corner for voxel (ix, iy, iz).
    /// @param ix Grid index in x-direction.
    /// @param iy Grid index in y-direction.
    /// @param iz Grid index in z-direction.
    /// @return AABB maximum corner (xup, yup, zup) [axis units].
    /// @throws std::runtime_error if any index is out of range.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    Eigen::Vector3d get_v3_AABB_max( const int ix, const int iy, const int iz ) const;

    /// @brief Get AABB3d for voxel (ix, iy, iz).
    /// @param ix Grid index in x-direction.
    /// @param iy Grid index in y-direction.
    /// @param iz Grid index in z-direction.
    /// @return AABB3d with min/max corners [axis units].
    /// @throws std::runtime_error if any index is out of range.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    AABB3d get_AABB3d( const int ix, const int iy, const int iz ) const;

    /// @brief Get x-axis bin index for a coordinate.
    /// @param value X-coordinate [axis units].
    /// @return Bin index in x-axis, or -1 if out of range.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    int get_ix( const double value ) const { return x_axis.get_index(value); };

    /// @brief Get y-axis bin index for a coordinate.
    /// @param value Y-coordinate [axis units].
    /// @return Bin index in y-axis, or -1 if out of range.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    int get_iy( const double value ) const { return y_axis.get_index(value); };

    /// @brief Get z-axis bin index for a coordinate.
    /// @param value Z-coordinate [axis units].
    /// @return Bin index in z-axis, or -1 if out of range.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    int get_iz( const double value ) const { return z_axis.get_index(value); };

    /// @brief Get center coordinate for x-axis bin.
    /// @param index Bin index in x-axis.
    /// @return Center coordinate [axis units].
    /// @throws std::runtime_error if index is out of range.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    double get_xcnt( const int index ) const { return x_axis.get_center_value(index); };

    /// @brief Get center coordinate for y-axis bin.
    /// @param index Bin index in y-axis.
    /// @return Center coordinate [axis units].
    /// @throws std::runtime_error if index is out of range.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    double get_ycnt( const int index ) const { return y_axis.get_center_value(index); };

    /// @brief Get center coordinate for z-axis bin.
    /// @param index Bin index in z-axis.
    /// @return Center coordinate [axis units].
    /// @throws std::runtime_error if index is out of range.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    double get_zcnt( const int index ) const { return z_axis.get_center_value(index); };

    /// @brief Get center coordinate for x-axis bin containing value.
    /// @param value X-coordinate [axis units].
    /// @return Center coordinate [axis units].
    /// @throws std::runtime_error if value is out of range.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    double get_xcnt( const double value ) const;

    /// @brief Get center coordinate for y-axis bin containing value.
    /// @param value Y-coordinate [axis units].
    /// @return Center coordinate [axis units].
    /// @throws std::runtime_error if value is out of range.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    double get_ycnt( const double value ) const;

    /// @brief Get center coordinate for z-axis bin containing value.
    /// @param value Z-coordinate [axis units].
    /// @return Center coordinate [axis units].
    /// @throws std::runtime_error if value is out of range.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    double get_zcnt( const double value ) const;

    /// @brief Get upper boundary coordinate for x-axis bin.
    /// @param index Bin index in x-axis.
    /// @return Upper boundary [axis units].
    /// @throws std::runtime_error if index is out of range.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    double get_xup( const int index ) const { return x_axis.get_upper_value(index); };

    /// @brief Get upper boundary coordinate for y-axis bin.
    /// @param index Bin index in y-axis.
    /// @return Upper boundary [axis units].
    /// @throws std::runtime_error if index is out of range.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    double get_yup( const int index ) const { return y_axis.get_upper_value(index); };

    /// @brief Get upper boundary coordinate for z-axis bin.
    /// @param index Bin index in z-axis.
    /// @return Upper boundary [axis units].
    /// @throws std::runtime_error if index is out of range.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    double get_zup( const int index ) const { return z_axis.get_upper_value(index); };

    /// @brief Get lower boundary coordinate for x-axis bin.
    /// @param index Bin index in x-axis.
    /// @return Lower boundary [axis units].
    /// @throws std::runtime_error if index is out of range.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    double get_xlow( const int index ) const { return x_axis.get_lower_value(index); };

    /// @brief Get lower boundary coordinate for y-axis bin.
    /// @param index Bin index in y-axis.
    /// @return Lower boundary [axis units].
    /// @throws std::runtime_error if index is out of range.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    double get_ylow( const int index ) const { return y_axis.get_lower_value(index); };

    /// @brief Get lower boundary coordinate for z-axis bin.
    /// @param index Bin index in z-axis.
    /// @return Lower boundary [axis units].
    /// @throws std::runtime_error if index is out of range.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    double get_zlow( const int index ) const { return z_axis.get_lower_value(index); };

    /// @brief Get upper boundary coordinate for x-axis bin containing value.
    /// @param value X-coordinate [axis units].
    /// @return Upper boundary [axis units].
    /// @throws std::runtime_error if value is out of range.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    double get_xup( const double value ) const;

    /// @brief Get upper boundary coordinate for y-axis bin containing value.
    /// @param value Y-coordinate [axis units].
    /// @return Upper boundary [axis units].
    /// @throws std::runtime_error if value is out of range.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    double get_yup( const double value ) const;

    /// @brief Get upper boundary coordinate for z-axis bin containing value.
    /// @param value Z-coordinate [axis units].
    /// @return Upper boundary [axis units].
    /// @throws std::runtime_error if value is out of range.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    double get_zup( const double value ) const;

    /// @brief Get lower boundary coordinate for x-axis bin containing value.
    /// @param value X-coordinate [axis units].
    /// @return Lower boundary [axis units].
    /// @throws std::runtime_error if value is out of range.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    
    double get_xlow( const double value ) const;
    /// @brief Get lower boundary coordinate for y-axis bin containing value.
    /// @param value Y-coordinate [axis units].
    /// @return Lower boundary [axis units].
    /// @throws std::runtime_error if value is out of range.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    double get_ylow( const double value ) const;

    /// @brief Get lower boundary coordinate for z-axis bin containing value.
    /// @param value Z-coordinate [axis units].
    /// @return Lower boundary [axis units].
    /// @throws std::runtime_error if value is out of range.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    double get_zlow( const double value ) const;

    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name setter_functions
    ///@{
    
    /// @brief Set instance name.
    /// @param name_in New name string.
    /// @note Thread-safe: No.
    void set_name( const std::string &name_in ){ name = name_in; };

    /// @brief Set x-axis name.
    /// @param name_in New name string.
    /// @note Thread-safe: No.
    void set_x_axis_name( const std::string &name_in ){ x_axis.set_name(name_in); };

    /// @brief Set y-axis name.
    /// @param name_in New name string.
    /// @note Thread-safe: No.
    void set_y_axis_name( const std::string &name_in ){ y_axis.set_name(name_in); };

    /// @brief Set z-axis name.
    /// @param name_in New name string.
    /// @note Thread-safe: No.
    void set_z_axis_name( const std::string &name_in ){ z_axis.set_name(name_in); };

    /// @brief Replace x-axis definition.
    /// @param x_axis_in New x-axis Grid1d.
    /// @details Also rebuilds cached_aabb3d_ so that ray tracing
    ///          always uses up-to-date grid bounds.
    /// @note Thread-safe: No.
    void set_x_axis( const Grid1d &x_axis_in ){ x_axis = x_axis_in; rebuild_cached_aabb3d(); };

    /// @brief Replace y-axis definition.
    /// @param y_axis_in New y-axis Grid1d.
    /// @details Also rebuilds cached_aabb3d_. See set_x_axis().
    /// @note Thread-safe: No.
    void set_y_axis( const Grid1d &y_axis_in ){ y_axis = y_axis_in; rebuild_cached_aabb3d(); };

    /// @brief Replace z-axis definition.
    /// @param z_axis_in New z-axis Grid1d.
    /// @details Also rebuilds cached_aabb3d_. See set_x_axis().
    /// @note Thread-safe: No.
    void set_z_axis( const Grid1d &z_axis_in ){ z_axis = z_axis_in; rebuild_cached_aabb3d(); };

    /// @brief Get cached AABB3d for the entire grid.
    /// @return Const reference to the cached AABB3d.
    /// @note Thread-safe: Yes. Time complexity: O(1).
    const AABB3d& get_cached_aabb3d() const { return cached_aabb3d_; }

    ///@} ------------------------------------------------------------------

    //===================================================
    /// @name path_length_Grid3d
    /// @{
    
    /// @brief Parameter bundle for 3D voxel traversal.
    /// @details Used internally by the DDA ray traversal helper.
    struct VoxelTraversalParameters3d {
      double x_start, y_start, z_start;
      int manhattan_distance;
      double step_x, step_y, step_z;
      double t_delta_x, t_delta_y, t_delta_z;
      double t_max_x0, t_max_y0, t_max_z0;
      double vx_sign, vy_sign, vz_sign;
    };

    /// @brief Build DDA traversal parameters from a valid ray-AABB intersection.
    /// @param ray3d Ray with origin and direction.
    /// @param eps Boundary tolerance scalar.
    /// @param tmin Entry parameter from ray-AABB intersection.
    /// @param tmax Exit parameter from ray-AABB intersection.
    /// @return Populated VoxelTraversalParameters3d ready for DDA loop.
    /// @throws std::runtime_error if start indices are out of range.
    /// @note Caller must verify tf_hit_volume and tmin/tmax signs beforehand.
    VoxelTraversalParameters3d build_traversal_params(
        const Ray3d &ray3d, double eps, double tmin, double tmax) const;

    /// @brief Helper function for DDA voxel traversal.
    /// @param prm Precomputed traversal parameters.
    /// @return Vector of hit voxel indices (ix, iy, iz).
    /// @note Thread-safe: Yes. Time complexity: O(N) with N = path length voxels.
    /// @ingroup ray_tracing_functions
    std::vector<Ixiyiz>
      get_hit_index_by_ray_tracing_algorithm(const Grid3d::VoxelTraversalParameters3d &prm) const;

    /// @brief Trace a ray through the grid and return hit voxel indices.
    /// @param ray3d Ray3d with origin and direction (direction is normalized in Ray3d).
    /// @param eps Small scalar to avoid missing boundary voxels (relative to grid spacing).
    /// @return Vector of hit voxel indices (ix, iy, iz), ordered along the ray.
    /// @throws std::runtime_error if the entry/exit indices are out of range.
    /// @note Uses cached_aabb3d_ to skip per-ray AABB construction (optimization #4/#5).
    /// @note Thread-safe: Yes. Time complexity: O(N) with N = path length voxels.
    /// @note Coordinate system: right-handed, z-up; units follow the axis units.
    /// @ingroup ray_tracing_functions
    std::vector<Ixiyiz>
      get_hit_voxel_index( const Ray3d &ray3d, const double eps ) const;

    /// @brief Trace a ray through the grid with beam length warning.
    /// @param ray3d Ray3d with origin and direction.
    /// @param BL_max Maximum beam length threshold [meters]. If <= 0, warning is disabled.
    /// @param eps Small scalar to avoid missing boundary voxels.
    /// @return Vector of hit voxel indices, ordered along the ray.
    /// @note Uses cached_aabb3d_ and a single intersection test to avoid the
    ///       double ray-AABB intersection that the legacy version performed
    ///       (optimization #4).
    /// @note If tmax > BL_max, logs a warning via LOG_WARN_ND but continues processing.
    /// @note Thread-safe: Yes.
    /// @ingroup ray_tracing_functions
    std::vector<Ixiyiz>
      get_hit_voxel_index( const Ray3d &ray3d, const double BL_max, const double eps ) const;

    /// @brief [Legacy] Trace a ray — constructs a local AABB each call.
    /// @deprecated Retained for regression testing; prefer get_hit_voxel_index.
    std::vector<Ixiyiz>
      get_hit_voxel_index_legacy( const Ray3d &ray3d, const double eps ) const;

    /// @brief [Legacy] Trace with BL_max — performs double ray-AABB intersection.
    /// @deprecated Retained for regression testing; prefer get_hit_voxel_index.
    std::vector<Ixiyiz>
      get_hit_voxel_index_legacy( const Ray3d &ray3d, const double BL_max, const double eps ) const;

    /// @brief Traverse ray through grid, calling visitor for each hit voxel
    ///        with its path length computed inside the DDA loop.
    /// @tparam Visitor Callable with signature void(const Ixiyiz&, double delta_path).
    /// @param ray3d Ray3d with origin and direction.
    /// @param eps Small scalar to avoid missing boundary voxels.
    /// @param visitor Callback invoked for each hit voxel with its path length.
    /// @note Combines DDA traversal and path length computation in a single pass,
    ///       eliminating per-voxel ray-AABB intersection (optimization #1)
    ///       and per-ray std::vector heap allocation (optimization #3).
    /// @note Uses cached_aabb3d_ (optimization #5).
    /// @note Thread-safe: Yes (no internal mutation).
    /// @note Implementation in header — required for template instantiation.
    template<typename Visitor>
    void traverse_ray_with_pathlength(
        const Ray3d &ray3d, double eps, Visitor&& visitor) const;

    // end group path_length_Grid3d
    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name output_functions
    ///@{

    /// @brief Write a single-line header describing grid extents.
    /// @param fout Output file handle (must be non-null).
    /// @throws std::runtime_error if fout is null.
    /// @note Output format: "# Grid3d name=... nbinx=... xmin=... xmax=... ..."
    /// @note Thread-safe: Yes (no internal mutation). Time complexity: O(1).
    void out_header( FILE *fout );

    /// @brief Write 2D grid bounds (x/y only) to file.
    /// @param pathout Output path for ASCII file.
    /// @throws std::runtime_error if the output file cannot be opened.
    /// @note Output columns: ix iy xlow xup ylow yup (all in axis units).
    /// @note Thread-safe: Yes (const). Time complexity: O(nbinx*nbiny).
    void out_grid_2d(const std::filesystem::path& pathout) const;

    /// @brief Write 3D grid bounds (x/y/z) to file.
    /// @param pathout Output path for ASCII file.
    /// @throws std::runtime_error if the output file cannot be opened.
    /// @note Output columns: ix iy iz xlow xup ylow yup zlow zup (axis units).
    /// @note Thread-safe: Yes (const). Time complexity: O(nbinx*nbiny*nbinz).
    void out_grid_3d(const std::filesystem::path& pathout) const;

    /// @brief Print voxel index and coordinate bounds for a position.
    /// @param fout Output file handle (must be non-null).
    /// @param x X-coordinate [axis units].
    /// @param y Y-coordinate [axis units].
    /// @param z Z-coordinate [axis units].
    /// @throws std::runtime_error if (x,y,z) is outside the grid.
    /// @note Output format: "ix,iy,iz | x:[low,up] y:[low,up] z:[low,up]".
    /// @note Thread-safe: Yes (const). Time complexity: O(1).
    void display_voxel_index( FILE *fout,const double &x, const double &y, const double &z) const;

    /// @brief Print voxel coordinate bounds for an index triplet.
    /// @param fout Output file handle (must be non-null).
    /// @param ix Grid index in x-direction.
    /// @param iy Grid index in y-direction.
    /// @param iz Grid index in z-direction.
    /// @throws std::runtime_error if any index is out of range.
    /// @note Output format: "ix,iy,iz | x:[low,up] y:[low,up] z:[low,up]".
    /// @note Thread-safe: Yes (const). Time complexity: O(1).
    void display_voxel_coordinate( FILE *fout, const int &ix, const int &iy, const int &iz) const ;
    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name calculation_functions
    ///@{
    
    /// @brief Compute 3D Euclidean distance between two voxels.
    /// @param ixiyiz0 Index triplet of the first voxel.
    /// @param ixiyiz1 Index triplet of the second voxel.
    /// @return Distance in axis units (typically meters).
    /// @note Thread-safe: Yes. Time complexity: O(1).
    double calc_dist_dxyz(
      const Ixiyiz &ixiyiz0, const Ixiyiz &ixiyiz1 ) const;

    /// @brief Compute horizontal (x/y) distance between two voxels.
    /// @param ixiyiz0 Index triplet of the first voxel.
    /// @param ixiyiz1 Index triplet of the second voxel.
    /// @return Horizontal distance in axis units (typically meters).
    /// @note Thread-safe: Yes. Time complexity: O(1).
    double calc_dist_dxy(
      const Ixiyiz &ixiyiz0, const Ixiyiz &ixiyiz1 ) const;

    /// @brief Compute vertical (z) distance between two voxels.
    /// @param ixiyiz0 Index triplet of the first voxel.
    /// @param ixiyiz1 Index triplet of the second voxel.
    /// @return Vertical distance in axis units (typically meters).
    /// @note Thread-safe: Yes. Time complexity: O(1).
    double calc_dist_dz(
      const Ixiyiz &ixiyiz0, const Ixiyiz &ixiyiz1 ) const;

    ///@} ------------------------------------------------------------------

    //===============================================
    /// @name binary_io_functions
    /// @{

    /// @brief Save Grid3d to a binary stream.
    /// @param ofs Output stream (already opened).
    /// @note Thread-safe: Yes (const). Time complexity: O(1).
    void save( std::ofstream &ofs ) const;

    /// @brief Load Grid3d from a binary stream (pass-through via Grid1d::load).
    /// @param ifs Input stream (already opened).
    /// @param[in] tolerance_ratio  Tolerance for Grid1d consistency checks
    /// @note Thread-safe: No (mutating). Time complexity: O(1).
    void load( std::ifstream &ifs, double tolerance_ratio = 1.0e-4 );
    ///@} ------------------------------------------------------------------
};

//======================================================================
// Template implementation — must remain in header for instantiation.
//======================================================================

template<typename Visitor>
void Grid3d::traverse_ray_with_pathlength(
    const Ray3d &ray3d, double eps, Visitor&& visitor) const
{
  // --- ray-AABB intersection using cached grid bounds (opt #5) ---
  auto [tf_hit_volume, tmin, tmax] = ray3d.is_intersect(cached_aabb3d_);
  if (!tf_hit_volume) return;
  if (tmin < 0 && tmax < 0) return;

  // delegate DDA setup to shared helper — eliminates duplicated ~80-line block
  const auto prm = build_traversal_params(ray3d, eps, tmin, tmax);

  int ix = get_ix(prm.x_start);
  int iy = get_iy(prm.y_start);
  int iz = get_iz(prm.z_start);

  double t_max_x = prm.t_max_x0;
  double t_max_y = prm.t_max_y0;
  double t_max_z = prm.t_max_z0;

  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  const int nbinz = get_nbinz();

  // --- DDA traversal with inline path-length computation (opt #1) ---
  // t_prev tracks the cumulative path at the entry of the current voxel.
  // t_next = min(|t_max|) gives the boundary crossing that exits the voxel.
  // delta_path = t_next - t_prev is the path through the current voxel,
  // computed from DDA state — no per-voxel ray-AABB intersection needed.
  double t_prev = 0.0;

  for (int t = 0; t <= prm.manhattan_distance; t++) {
    // next boundary crossing distance (absolute value since t_max can be negative)
    const double t_next =
        std::min({std::fabs(t_max_x), std::fabs(t_max_y), std::fabs(t_max_z)});

    // emit current voxel with path length computed from DDA state (opt #1)
    visitor(Ixiyiz{ix, iy, iz}, t_next - t_prev);
    t_prev = t_next;

    // advance to next voxel (same branching as original DDA)
    if (std::fabs(t_max_x) < std::fabs(t_max_y)) {
      if (std::fabs(t_max_x) < std::fabs(t_max_z)) {
        t_max_x += prm.t_delta_x;
        ix += static_cast<int>(prm.vx_sign);
        if (ix < 0 || ix > nbinx - 1) return;
      } else {
        t_max_z += prm.t_delta_z;
        iz += static_cast<int>(prm.vz_sign);
        if (iz < 0 || iz > nbinz - 1) return;
      }
    } else {
      if (std::fabs(t_max_y) < std::fabs(t_max_z)) {
        t_max_y += prm.t_delta_y;
        iy += static_cast<int>(prm.vy_sign);
        if (iy < 0 || iy > nbiny - 1) return;
      } else {
        t_max_z += prm.t_delta_z;
        iz += static_cast<int>(prm.vz_sign);
        if (iz < 0 || iz > nbinz - 1) return;
      }
    }
  }

  // final voxel after all manhattan_distance steps — compute its exit path
  const double t_final =
      std::min({std::fabs(t_max_x), std::fabs(t_max_y), std::fabs(t_max_z)});
  visitor(Ixiyiz{ix, iy, iz}, t_final - t_prev);
}
