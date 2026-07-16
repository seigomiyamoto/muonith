/// @file cls_Grid2dVoxel.hpp
/// @brief 2D grid with per-cell voxel data and a fixed vertical extent.
/// @details Defines Grid2dVoxel, a Grid2d-derived container that assigns a Voxel
///          to each (ix, iy) cell with a uniform z-range [zmin, zmax].
///
///          Typical workflow:
///          1. Build Grid1d x/y axes and construct Grid2dVoxel with zmin/zmax.
///          2. Allocate voxel storage with vec_vec_memory_allocate(n_detector).
///          3. Set voxel flags/densities using callVoxel/getVoxel.
///          4. Export with out_voxel/out_voxel_all(_binary) or save/load.
///
///          Coordinate system and units:
///          - Right-handed, z-up; zmin/zmax are vertical bounds [m].
///          - x/y coordinates use the Grid1d axis units [m].
///          - Density is stored in kg/m^3 (SI units).
///
///          Memory layout:
///          - vec_vec_Voxel.at(iy).at(ix)
///          - Outer vector: y-direction (size = nbiny)
///          - Inner vector: x-direction (size = nbinx)
///
///          Thread safety:
///          - Const accessors: thread-safe (read-only)
///          - Mutating methods: not thread-safe
///          - vec_vec_memory_allocate uses OpenMP when enabled
///
///          I/O formats:
///          - out_voxel/out_voxel_all: ASCII lines
///            (tf_exist, xcnt[m], ycnt[m], zcnt[m], density[kg/m^3],
///             n_det, tf_hit..., n_hit_det)
///          - out_voxel_all_binary/save/load: native-endian binary
#pragma once

#include <map>
#include <fstream>
#include <iostream>
#include <sstream> // istringstream
#include <string>
#include <cstdio>
#include <cstdint>
#include <cmath>
#include <functional>  // for sorting
#include <algorithm> // for sorting
#include <vector>
#include <filesystem>

#include <Eigen/Dense>
#include "cls_Ray.hpp"
#include "cls_AABB.hpp"

#include "ns_tuple_int.hpp"
#include "cls_Voxel.hpp"
#include "cls_Grid2d.hpp"

namespace fs = std::filesystem;

//########################################################################
//########################################################################
/// @class Grid2dVoxel
/// @brief A 2D grid that associates a Voxel with each cell and a fixed z-range
///
/// @details This class inherits from Grid2d and stores a Voxel per (ix, iy)
///          grid cell. All voxels share the same vertical span [zmin, zmax].
///
///          Primary use cases:
///          - 2D density maps with per-cell hit flags
///          - Exporting voxelized density fields for ray/path calculations
///          - A lightweight counterpart to Grid3dVoxel in 2D setups
///
///          Example usage:
///          @code
///          Grid1d x_axis("x", 10, 0.0, 100.0, 10.0);
///          Grid1d y_axis("y", 20, 0.0, 200.0, 10.0);
///          Grid2dVoxel g2vox(x_axis, y_axis, 0.0, 50.0, 4);
///
///          g2vox.callVoxel(0, 0).set_tf_exist(true);
///          g2vox.callVoxel(0, 0).set_density(2500.0); // kg/m^3
///          @endcode
///
/// @ingroup basicGridClasses
/// @ingroup terrainClasses
//########################################################################
//########################################################################
class Grid2dVoxel : public Grid2d {
  private:
    /// @brief name of this instance
    std::string name = "g2vox";

    /// @brief vertical range [m]
    double zmin = 999.0;
    double zmax = -999.0;

    /// @brief contains Voxels as vector of vector
    /// @details access order is vec_vec_Voxel.at(iy).at(ix)
    std::vector<std::vector<Voxel>> vec_vec_Voxel;

  public:
    //======================================================================
    /// @name constructor_destructor
    ///@{

    /// @brief Default constructor
    Grid2dVoxel() = default;

    /// @brief Copy constructor
    Grid2dVoxel( const Grid2dVoxel &org ) = default;

    /// @brief Destructor
    ~Grid2dVoxel() = default;

    /// @brief Build from Grid1d x/y axes and vertical bounds
    /// @param x_axis_in x-axis definition [m]
    /// @param y_axis_in y-axis definition [m]
    /// @param zmin_in Lower z bound for all voxels [m]
    /// @param zmax_in Upper z bound for all voxels [m]
    /// @param n_detector Number of detector channels per voxel (>= 0)
    /// @throws std::runtime_error if n_detector is negative
    /// @note Thread-safe: No (allocates and mutates internal storage)
    Grid2dVoxel(
      const Grid1d &x_axis_in, const Grid1d &y_axis_in,
      const double zmin_in, const double zmax_in, const int n_detector)
    : Grid2d(x_axis_in,y_axis_in), zmin(zmin_in), zmax(zmax_in)
    { vec_vec_memory_allocate(n_detector); };

    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name build functions
    ///@{

    /// @brief Allocate voxel storage and initialize defaults
    /// @param n_det Number of detector channels per voxel (>= 0)
    /// @throws std::runtime_error if n_det is negative
    /// @note Thread-safe: No (mutates internal storage)
    /// @note Time complexity: O(nbinx × nbiny)
    /// @note Uses OpenMP parallelization when enabled
    void vec_vec_memory_allocate(const int n_det);
    
    /// @brief OpenMP-enabled allocation wrapper
    /// @details Currently delegates to vec_vec_memory_allocate
    /// @param n_det Number of detector channels per voxel (>= 0)
    /// @throws std::runtime_error if n_det is negative
    /// @note Thread-safe: No (mutates internal storage)
    void mp_vec_vec_memory_allocate(const int n_det);

    ///@} ------------------------------------------------------------------


    //======================================================================
    /// @name operator
    ///@{

    /// @brief Compare voxel storage content between two grids
    /// @param g2vox_in Other grid to compare against
    /// @return true if sizes and voxel contents match
    /// @note Thread-safe: Yes (read-only)
    bool is_vec_vec_Voxel_same( const Grid2dVoxel &g2vox_in ) const;

    /// @brief Non-equality operator
    /// @param other Other grid to compare against
    /// @return true if any voxel content or z-range differs
    /// @note Instance names are intentionally ignored
    /// @note Thread-safe: Yes (read-only)
    bool operator!=(const Grid2dVoxel& other) const;

    /// @brief Equality operator (defined via operator!=)
    /// @param other Other grid to compare against
    /// @return true if voxel content and z-range are identical
    /// @note Instance names are intentionally ignored
    /// @note Thread-safe: Yes (read-only)
    bool operator==(const Grid2dVoxel& other) const {
      return !(*this != other);
    }
    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name getter_functions
    ///@{

    /// @brief Get name
    /// @return Instance name
    /// @note Thread-safe: Yes (read-only)
    std::string get_name() const { return name; };

    /// @brief Get voxel by grid indices
    /// @param ix Grid index in x-direction
    /// @param iy Grid index in y-direction
    /// @return Const reference to voxel at (ix, iy)
    /// @throws std::runtime_error if indices are out of bounds
    /// @note Thread-safe: Yes (read-only)
    const Voxel& getVoxel( const int ix, const int iy ) const;

    /// @brief Get voxel by coordinates
    /// @param x X coordinate [m]
    /// @param y Y coordinate [m]
    /// @return Const reference to voxel at (x, y)
    /// @throws std::runtime_error if (x, y) is outside grid bounds
    /// @note Thread-safe: Yes (read-only)
    const Voxel& getVoxel( const double x, const double y ) const;

    /// @brief Return a flattened density vector for all voxels
    /// @details Flatten order is (iy, ix): index = iy * nbinx + ix
    /// @return Vector of densities [kg/m^3] (stored as float)
    /// @note Thread-safe: Yes (read-only)
    /// @note Time complexity: O(nbinx × nbiny)
    Eigen::VectorXf get_vecxf_density() const;

    /// @brief Get voxel density by grid indices
    /// @param ix Grid index in x-direction
    /// @param iy Grid index in y-direction
    /// @return Density [kg/m^3]
    /// @throws std::runtime_error if indices are out of bounds
    /// @note Thread-safe: Yes (read-only)
    double get_density( const int ix, const int iy ) const {
      return getVoxel(ix,iy).get_density();
    };
    
    /// @brief Get AABB2d for (ix, iy)
    /// @param ix Grid index in x-direction
    /// @param iy Grid index in y-direction
    /// @return 2D axis-aligned bounding box in XY [m]
    /// @throws std::runtime_error if indices are out of bounds
    /// @note Thread-safe: Yes (read-only)
    AABB2d get_AABB2d(const int ix, const int iy) const;
    
    /// @brief Get AABB3d for (ix, iy) with z in [zmin, zmax]
    /// @param ix Grid index in x-direction
    /// @param iy Grid index in y-direction
    /// @return 3D axis-aligned bounding box in XYZ [m]
    /// @throws std::runtime_error if indices are out of bounds
    /// @note Thread-safe: Yes (read-only)
    AABB3d get_AABB3d(const int ix, const int iy) const;
    
    /// @brief Get (density, AABB2d) for (ix, iy)
    /// @param ix Grid index in x-direction
    /// @param iy Grid index in y-direction
    /// @return Tuple of density [kg/m^3] and AABB2d [m]
    /// @throws std::runtime_error if indices are out of bounds
    /// @note Thread-safe: Yes (read-only)
    std::tuple< double, AABB2d >
      get_density_AABB2d(const int ix, const int iy) const;

    /// @brief Get (density, AABB3d) for (ix, iy)
    /// @param ix Grid index in x-direction
    /// @param iy Grid index in y-direction
    /// @return Tuple of density [kg/m^3] and AABB3d [m]
    /// @throws std::runtime_error if indices are out of bounds
    /// @note Thread-safe: Yes (read-only)
    std::tuple< double, AABB3d >
      get_density_AABB3d(const int ix, const int iy) const;

    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name call_functions
    /// @brief Get mutable references to voxels
    ///@{
    
    /// @brief Get mutable voxel by grid indices
    /// @param ix Grid index in x-direction
    /// @param iy Grid index in y-direction
    /// @return Mutable reference to voxel at (ix, iy)
    /// @throws std::runtime_error if indices are out of bounds
    /// @note Thread-safe: No (non-const method)
    Voxel& callVoxel( const int ix, const int iy);

    /// @brief Get mutable voxel by coordinates
    /// @param x X coordinate [m]
    /// @param y Y coordinate [m]
    /// @return Mutable reference to voxel at (x, y)
    /// @throws std::runtime_error if (x, y) is outside grid bounds
    /// @note Thread-safe: No (non-const method)
    Voxel& callVoxel( const double x, const double y);
    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name check_functions
    ///@{

    /// @brief Check whether (ix, iy) is outside grid bounds
    /// @param ix Grid index in x-direction
    /// @param iy Grid index in y-direction
    /// @return true if out of range; false otherwise
    /// @note Thread-safe: Yes (read-only)
    bool is_ixiy_out_of_range( const int ix, const int iy ) const {
      if( ix < 0 ) return true;
      if( iy < 0 ) return true;
      if( ix > get_nbinx()-1 ) return true;
      if( iy > get_nbiny()-1 ) return true;
      return false;
    };
    ///@} ------------------------------------------------------------------
    
    //======================================================================
    /// @name read_write_functions
    ///@{
      
    /// @brief Output voxel info in ASCII format
    /// @details Format:
    ///          tf_exist xcnt ycnt zcnt density n_det [tf_hit...] n_hit_det
    /// @param fout Output file pointer (must be valid)
    /// @param ix Grid index in x-direction
    /// @param iy Grid index in y-direction
    /// @throws std::runtime_error if fout is null or indices are out of bounds
    /// @note Thread-safe: Yes (read-only)
    void out_voxel( FILE *fout, const int ix, const int iy) const;
    
    /// @brief Output all voxels in ASCII format
    /// @param fout Output file pointer (must be valid)
    /// @throws std::runtime_error if fout is null
    /// @note Thread-safe: Yes (read-only)
    /// @note Time complexity: O(nbinx × nbiny)
    void out_voxel_all( FILE *fout ) const;

    /// @brief Output all voxels in binary format
    /// @details Native-endian binary layout in (iy, ix) order:
    ///          tf_exist(uint8), xcnt(double), ycnt(double), zcnt(double),
    ///          density(double), n_det(int32), tf_hit(uint8) * n_det,
    ///          n_hit_det(int32)
    /// @param fout Output file pointer (must be valid)
    /// @throws std::runtime_error if fout is null or write fails
    /// @note Thread-safe: Yes (read-only)
    /// @note Time complexity: O(nbinx × nbiny)
    void out_voxel_all_binary( FILE *fout ) const;
    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name setter_functions
    ///@{

    /// @brief Copy from another instance
    /// @param org Source instance
    /// @note Thread-safe: No (mutates internal state)
    void set( const Grid2dVoxel &org ){ *this = org; };

    /// @brief Set name of this instance
    /// @param name_in Name to assign
    /// @note Thread-safe: No (mutates internal state)
    void set_name( const std::string name_in ){ name=name_in; };

    /// @brief Set lower z bound
    /// @param zmin_in Lower z bound [m]
    /// @note Thread-safe: No (mutates internal state)
    void set_zmin( const double zmin_in ){ zmin = zmin_in; };

    /// @brief Set upper z bound
    /// @param zmax_in Upper z bound [m]
    /// @note Thread-safe: No (mutates internal state)
    void set_zmax( const double zmax_in ){ zmax = zmax_in; };

    /// @brief Replace voxel storage with input data
    /// @param vec_vec_Voxel_in Source voxel grid (size nbiny × nbinx)
    /// @throws std::runtime_error if size mismatch or grid is non-rectangular
    /// @note Thread-safe: No (mutates internal state)
    void set_vec_vec_Voxel(
      const std::vector<std::vector<Voxel>> &vec_vec_Voxel_in );

    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name binary_io_Grid2dVoxel
    ///@{

    /// @brief Save voxel storage to a binary stream
    /// @param ofs Output stream (binary)
    /// @throws std::runtime_error if write fails
    /// @note Thread-safe: Yes (read-only)
    void save_vec_vec_Voxel( std::ofstream &ofs ) const;
    
    /// @brief Load voxel storage from a binary stream
    /// @param ifs Input stream (binary)
    /// @throws std::runtime_error if read fails
    /// @note Thread-safe: No (mutates internal state)
    void load_vec_vec_Voxel(std::ifstream& ifs);

    /// @brief Save Grid2dVoxel to a binary stream
    /// @param ofs Output stream (binary)
    /// @throws std::runtime_error if write fails
    /// @note Thread-safe: Yes (read-only)
    void save( std::ofstream &ofs ) const;

    /// @brief Save Grid2dVoxel to a binary file
    /// @param pathout Output file path
    /// @throws std::runtime_error if file cannot be opened or write fails
    /// @note Thread-safe: Yes (read-only)
    void save( const fs::path& pathout) const;

    /// @brief Load Grid2dVoxel from a binary stream
    /// @param ifs Input stream (binary)
    /// @throws std::runtime_error if read fails
    /// @note Thread-safe: No (mutates internal state)
    void load( std::ifstream &ifs );

    /// @brief Load Grid2dVoxel from a binary file
    /// @param path_in Input file path
    /// @throws std::runtime_error if file cannot be opened or read fails
    /// @note Thread-safe: No (mutates internal state)
    void load( const fs::path &path_in );
    ///@} ------------------------------------------------------------------    

};
