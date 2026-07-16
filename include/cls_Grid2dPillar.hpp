/// @file cls_Grid2dPillar.hpp
/// @brief 2D grid mapped onto a 3D cuboid geometry for terrain and density modeling
///
/// @details This file defines the Grid2dPillar class, which extends Grid2d to represent
///          3D terrain or density structures by associating a vertical cuboid with each
///          2D grid cell.
///
///          Typical workflow:
///          1. Load DEM (Digital Elevation Model) data from file
///          2. Grid structure is automatically inferred from point cloud
///          3. Add density structures (cylinders, dikes, checkerboards) if needed
///          4. Use for muon radiography path calculations or cross-section analysis
///
///          Coordinate system:
///          - z-axis is vertical (positive upward)
///          - Elevations in meters above sea level (m a.s.l.)
///          - Densities in kg/m³ (SI units)
///
///          Memory layout:
///          - vec_vec_Pillar is column-major: access as vec_vec_Pillar.at(iy).at(ix)
///          - Outer vector: y-direction (size = nbiny)
///          - Inner vector: x-direction (size = nbinx)
///
///          Thread safety:
///          - Read operations (const methods): thread-safe
///          - Write operations: not thread-safe, external synchronization required
///          - Some methods use OpenMP parallelization internally (documented per-method)
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
#include "cls_AABB.hpp"
#include "cls_Angle.hpp"
#include "cls_Pillar.hpp"
#include "cls_Grid2d.hpp"
#include "ns_type_definitions.hpp"
#include "ns_iodir.hpp"

class Grid2dVoxel;
class Ray3d;

//################################################################################
//################################################################################
/// @class Grid2dPillar
/// @brief A class that places a Pillar at every grid point of a 2D grid (Grid2d)
///
/// @details This class inherits from Grid2d and extends it with 3D capabilities
///          by associating a vertical cuboid (Pillar) with each 2D grid cell.
///
///          Primary use cases:
///          - Digital Elevation Model (DEM) representation for terrain analysis
///          - Muon radiography path length calculations through geological structures
///          - Ray tracing through heterogeneous density distributions
///          - Cross-sectional density analysis at various elevations
///
///          Example usage:
///          @code
///          // Load terrain from DEM file
///          Grid2dPillar terrain("terrain.xyz", zmin=0.0, density=2500.0,
///                               tf_shift_x=false, tf_shift_y=false, tolerance=0.01);
///
///          // Add a cylindrical density anomaly
///          Grid2dPillar::VerticalCylinderParameters cyl;
///          cyl.v2_pos_cnt = Eigen::Vector2d(100.0, 200.0);
///          cyl.v2_length = Eigen::Vector2d(50.0, 50.0);  // radii
///          cyl.delta_density = 500.0;  // kg/m³
///          terrain.add_density_structure(cyl);
///
///          // Calculate muon path length through terrain
///          Ray3d muon_ray(pos, dir);
///          Ixiy hit_cell = {10, 20};
///          double path_length = terrain.get_delta_path(hit_cell, muon_ray);
///          @endcode
///
/// @ingroup basicGridClasses
/// @ingroup terrainClasses
//################################################################################
//################################################################################
class Grid2dPillar : public Grid2d {
  public:
    //======================================================================
    /// @name forward declaration of parameter class for build instance
    ///@{
    class Parameters;
    class VerticalCylinderParameters;
    class VerticalDikeParameters;
    class VerticalCheckerBoardParameters;
    ///@} ------------------------------------------------------------------

  private:
    /// @brief name of this instance
    std::string name = "non";

    /// @brief data storage of every Pillar
    /// @details be careful of index order, z.at(iy).at(ix) (column major)
    std::vector<std::vector<Pillar>> vec_vec_Pillar;

  public:
    //==================================================================
    /// @name const_values
    ///@{
    
    /// @brief if Pillar is not hit, return PATH_NO_HIT
    static constexpr double PATH_NO_HIT = 0.0;
    ///@} ------------------------------------------------------------------

    //==================================================================
    /// @name constructor_destructor
    ///@{
    
    /// @brief default constructor
    Grid2dPillar() = default;
    
    /// @brief copy constructor
    Grid2dPillar(const Grid2dPillar &org) = default;

    /// @brief move constructor
    Grid2dPillar(Grid2dPillar &&other) noexcept = default;

    /// @brief destructor    
    virtual ~Grid2dPillar() = default;
    
    /// @brief assignment operator
    Grid2dPillar& operator=(const Grid2dPillar&) = default;

    /// @brief not equal operator
    bool operator!=(const Grid2dPillar &other) const;

    /// @brief Construct Grid2dPillar from DEM file
    /// @param path_in Path to DEM file (ASCII or binary format with x,y,z data)
    /// @param zmin Base elevation for all cuboids [m a.s.l.]
    /// @param density_in Initial uniform density for all cuboids [kg/m³]
    /// @param tf_shift_x If true, shift x-axis by half bin width
    /// @param tf_shift_y If true, shift y-axis by half bin width
    /// @param tolerance_ratio Tolerance for grid regularity detection (0.0-1.0)
    /// @throws std::runtime_error if file cannot be read or grid structure cannot be inferred
    /// @note Calls read_dem internally
    Grid2dPillar(const std::filesystem::path &path_in,
                 const double zmin,
                 const double density_in,
                 const bool tf_shift_x,
                 const bool tf_shift_y,
                 const double tolerance_ratio);

    /// @brief Construct Grid2dPillar from DEM file with specified name
    /// @param name_in Name for this instance
    /// @param path_in Path to DEM file (ASCII or binary format with x,y,z data)
    /// @param zmin Base elevation for all cuboids [m a.s.l.]
    /// @param density_in Initial uniform density for all cuboids [kg/m³]
    /// @param tf_shift_x If true, shift x-axis by half bin width
    /// @param tf_shift_y If true, shift y-axis by half bin width
    /// @param tolerance_ratio Tolerance for grid regularity detection (0.0-1.0)
    /// @throws std::runtime_error if file cannot be read or grid structure cannot be inferred
    /// @note Calls read_dem and set_name internally
    Grid2dPillar(const std::string &name_in,
                 const std::filesystem::path &path_in,
                 const double zmin,
                 const double density_in,
                 const bool tf_shift_x,
                 const bool tf_shift_y,
                 const double tolerance_ratio);

    /// @brief Construct Grid2dPillar from Parameters structure
    /// @param prm Parameters structure containing all configuration
    /// @throws std::runtime_error if file cannot be read or grid structure cannot be inferred
    /// @note Calls read_dem, set_name, and add_density_structure (for cylinder, dike, checkerboard) internally
    Grid2dPillar(const Grid2dPillar::Parameters &prm);
    
    ///@} ------------------------------------------------------------------

    //==================================================================
    /// @name call functions
    /// @brief Get mutable references to Pillar objects
    /// @remark Pillar is not const, so you can modify the returned object
    ///@{

    /// @brief Get mutable reference to Pillar at specified position
    /// @param x X-coordinate [m]
    /// @param y Y-coordinate [m]
    /// @return Mutable reference to Pillar at (x, y)
    /// @throws std::runtime_error if (x, y) is outside grid bounds
    /// @note Thread-safe: No (non-const method)
    Pillar& callPillar( const double x, const double y );

    /// @brief Get mutable reference to Pillar at specified grid indices
    /// @param ix Grid index in x-direction
    /// @param iy Grid index in y-direction
    /// @return Mutable reference to Pillar at (ix, iy)
    /// @throws std::runtime_error if indices are out of bounds
    /// @note Thread-safe: No (non-const method)
    Pillar& callPillar( const int ix, const int iy );

    /// @brief Get mutable reference to Pillar at specified grid indices
    /// @param tpl_hit_cub Grid indices (ix, iy)
    /// @return Mutable reference to Pillar at (ix, iy)
    /// @throws std::runtime_error if indices are out of bounds
    /// @note Thread-safe: No (non-const method)
    Pillar& callPillar( const Ixiy &tpl_hit_cub ){
      auto [ix,iy] = tpl_hit_cub;  return callPillar(ix,iy);
    };
    ///@} ------------------------------------------------------------------

    /// @brief Set name of this instance
    /// @param name_in Name to assign
    void set_name( const std::string name_in ){ name = name_in; };

    /// @brief Set uniform density to all Pillar objects
    /// @param density_in Density value [kg/m³]
    /// @note Thread-safe: uses OpenMP parallelization internally
    /// @note Time complexity: O(nbinx × nbiny)
    void set_uniform_density( const double density_in );

    /// @brief Allocate memory for vec_vec_Pillar
    /// @note Allocates nbiny × nbinx Pillar objects
    /// @note Memory layout is column-major: vec_vec_Pillar.at(iy).at(ix)
    void vec_vec_memory_allocate();

    /// @brief Reconstruct Grid2dPillar from point cloud data
    /// @param vec_xyz Vector of (x, y, z) points representing terrain surface
    ///                Will be sorted if not already in (y, x, z) order
    /// @param zmin Base elevation for all cuboids [m a.s.l.]
    /// @param density_in Initial uniform density for all cuboids [kg/m³]
    /// @param tf_shift_x If true, shift x-axis by half bin width
    /// @param tf_shift_y If true, shift y-axis by half bin width
    /// @param tolerance_ratio Tolerance for grid regularity detection (0.0-1.0)
    /// @throws std::runtime_error if grid structure cannot be inferred even after sorting
    /// @note Thread-safe: No (modifies internal state)
    /// @note Time complexity: O(N log N) for sorting + O(N) for assignment, where N = vec_xyz.size()
    /// @note Points outside computed grid bounds are silently skipped
    void convert_from_vec_xyz(
        std::vector<std::array<double,3>> &vec_xyz
      , const double zmin, const double density_in
      , const bool tf_shift_x, const bool tf_shift_y
      , const double tolerance_ratio );

    //==================================================================
    /// @name getter_functions
    ///@{

    /// @brief Get immutable reference to Pillar at specified position
    /// @param x X-coordinate [m]
    /// @param y Y-coordinate [m]
    /// @return Const reference to Pillar at (x, y)
    /// @throws std::runtime_error if (x, y) is outside grid bounds
    /// @note Thread-safe: Yes (const method, read-only access)
    /// @see callPillar for mutable access
    const Pillar& getPillar( const double x, const double y ) const;

    /// @brief Get immutable reference to Pillar at specified grid indices
    /// @param ix Grid index in x-direction
    /// @param iy Grid index in y-direction
    /// @return Const reference to Pillar at (ix, iy)
    /// @throws std::runtime_error if indices are out of bounds
    /// @note Thread-safe: Yes (const method, read-only access)
    /// @see callPillar for mutable access
    const Pillar& getPillar( const int ix, const int iy ) const;

    /// @brief Get immutable reference to Pillar at specified grid indices
    /// @param ixiy Grid indices (ix, iy)
    /// @return Const reference to Pillar at (ix, iy)
    /// @throws std::runtime_error if indices are out of bounds
    /// @note Thread-safe: Yes (const method, read-only access)
    const Pillar& getPillar( const Ixiy &ixiy ) const{
      const auto [ix,iy] = ixiy;
      return getPillar(ix,iy);
    };

    /// @brief Get the name of this instance
    /// @return Name string
    std::string get_name() const { return name; };

    /// @brief Get voxel density at specified grid indices
    /// @param ix Grid index in x-direction
    /// @param iy Grid index in y-direction
    /// @return Density value [kg/m³]
    double get_density( const int ix, const int iy) const {
      return getPillar(ix,iy).get_density();
    };

    /// @brief Get voxel density at specified grid indices
    /// @param ixiy Grid indices (ix, iy)
    /// @return Density value [kg/m³]
    double get_density( const Ixiy &ixiy ) const {
      return getPillar(ixiy).get_density();
    };

    /// @brief Get minimum zmax of Pillar within circular region
    /// @param x_center Center x-coordinate [m]
    /// @param y_center Center y-coordinate [m]
    /// @param radius Search radius [m]
    /// @return Minimum zmax value within specified circle [m a.s.l.]
    double getMinimumZmaxCircle(
      const double x_center, const double y_center, const double radius ) const;
    
    /// @brief Get AABB3d (bounding box) for specified grid cell
    /// @param ixiy_hit_cub Grid indices (ix, iy)
    /// @return AABB3d representing the cuboid's bounding box
    AABB3d get_AABB3d(const Ixiy &ixiy_hit_cub) const;

    /// @brief Get density and AABB3d (bounding box) for specified grid cell
    /// @param ix Grid index in x-direction
    /// @param iy Grid index in y-direction
    /// @return Tuple of (density [kg/m³], AABB3d)
    std::tuple< double, AABB3d >
      get_density_AABB3d(const int ix, const int iy) const;

    /// @brief Get elevation difference (terrain_zmax - point_z) at point location
    /// @param xyz Point coordinates (x, y, z) [m]
    /// @return Elevation difference [m], or const_dz_no_hit if point is outside grid
    double get_dz( const std::array<double,3> &xyz ) const;

    /// @brief Get elevation difference (terrain_zmax - point_z) at colored point location
    /// @param tp_xyzrgb Colored point coordinates (x, y, z, r, g, b)
    /// @return Elevation difference [m], or const_dz_no_hit if point is outside grid
    double get_dz( const std::tuple<double,double,double,int,int,int> &tp_xyzrgb ) const;

    /// @brief Constant returned by get_dz when point is outside grid bounds
    static constexpr double const_dz_no_hit = -987654321.0;

    /// @brief Get vector of dz (elevation differences) from point cloud
    /// @param vec_xyz Vector of (x, y, z) points
    /// @return Vector of dz values (terrain_zmax - point_z)
    /// @note Created for comparing LiDAR point cloud data with DEM
    std::vector<double> get_vec_dz(
      const std::vector<std::array<double,3>> &vec_xyz ) const;

    /// @brief Calculate depth from zmax at specified location along detector ray's elevation angle
    /// @param v2_pos_obj 2D position of object (x, y)
    /// @param ray_det Detector ray
    /// @param elev_ang_center_rad Central elevation angle [radians]
    /// @return Depth value [m]
    double calc_depth_from_zmax(
      const Eigen::Vector2d& v2_pos_obj,
      const Ray3d& ray_det,
      double elev_ang_center_rad) const;

    /// @brief Get vector of (x, y, z, dz) from point cloud
    /// @param vec_xyz Vector of (x, y, z) points
    /// @return Vector of (x, y, z, dz) where dz = terrain_zmax - point_z
    /// @note Created for comparing LiDAR point cloud data with DEM
    std::vector<std::array<double,4>> get_vec_xyzdz(
      const std::vector<std::array<double,3>> &vec_xyz ) const;

    /// @brief Get vector of (x, y, z, dz, r, g, b) from colored point cloud
    /// @param vec_xyzrgb Vector of (x, y, z, r, g, b) colored points
    /// @return Vector of (x, y, z, dz, r, g, b) with elevation differences
    /// @note Created for comparing LiDAR point cloud data with DEM
    std::vector<std::tuple<double,double,double,double,int,int,int>> get_vec_xyzdzrgb(
      const std::vector<std::tuple<double,double,double,int,int,int>> &vec_xyzrgb ) const;

    /// @brief Get terrain elevation profile (zmax - ray_z) along Ray3d
    /// @param ray3d Ray along which to sample terrain
    /// @param max_distance Maximum distance along ray to sample [m]
    /// @return Vector of (distance, dzmax) pairs
    std::vector<std::tuple<double,double>>
      get_dzmax_profile( const Ray3d &ray3d, const double max_distance ) const;
    
    /// @brief Save terrain elevation profile along ray to file
    /// @param outputpath Output file path
    /// @param ray3d Ray along which to sample terrain
    /// @param max_distance Maximum distance along ray to sample [m]
    /// @note Output format: ASCII text with (distance, dzmax) pairs
    /// @note Calls get_dzmax_profile internally and writes results to file
    void out_dzmax_profile(
        std::filesystem::path &outputpath
      , const Ray3d &ray3d , const double max_distance ) const;
    
    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name read_write_functions
    ///@{

    /// @brief Read DEM (Digital Elevation Model) from file and construct grid
    /// @param path_in Path to input file (ASCII xyz, or .g2zbin binary grid)
    /// @param zmin Base elevation for all cuboids [m a.s.l.]
    /// @param density Initial uniform density for all cuboids [kg/m³]
    /// @param tf_shift_x If true, shift x-axis by half bin width
    /// @param tf_shift_y If true, shift y-axis by half bin width
    /// @param tolerance_ratio Tolerance for grid regularity detection (0.0-1.0)
    /// @throws std::runtime_error if path_in is empty
    /// @throws std::runtime_error if file cannot be read
    /// @throws std::runtime_error if grid structure cannot be inferred
    /// @note Thread-safe: No (modifies internal state)
    /// @note Terrain heights below zmin are clamped to zmin using OpenMP parallelization
    /// @note For .g2zbin, uses the Grid2dXYZ constructor; for ASCII xyz, calls myapp::read_vec_xyz and convert_from_vec_xyz
    void read_dem( const std::filesystem::path &path_in
      , const double zmin, const double density
      , const bool tf_shift_x, const bool tf_shift_y
      , const double tolerance_ratio );

    /// @brief Output Grid2dPillar data to file
    /// @param pathout Output file path
    /// @throws std::runtime_error if file cannot be opened for writing
    void out( const std::filesystem::path& pathout ) const;

    /// @brief Output Grid2dPillar data to file using instance name
    /// @note Filename is automatically generated as get_name() + ".tmp"
    void out() const {
      const fs::path pathout = iodir::make_pathout( get_name() + ".tmp");
      out(pathout);
    };

    /// @brief Display information about nearest grid point
    /// @param x_in X-coordinate [m]
    /// @param y_in Y-coordinate [m]
    /// @note Displays grid indices and zmin/zmax values for debugging
    void disp_nearest_point( const double x_in, const double y_in ) const;

    //----------------------------------------------------------------------
    /// @brief Parameter structure for out_cross_section_z methods
    /// @details Defines spatial range, step size, and output format for cross-section output
    struct CrossSectionZParameters {
      double xmin = 999.0; ///< Minimum x coordinate for output [m]
      double xmax =-999.0; ///< Maximum x coordinate for output [m]
      double xstep=-1.0; ///< x-direction step size [m]
      double ymin = 999.0; ///< Minimum y coordinate for output [m]
      double ymax =-999.0; ///< Maximum y coordinate for output [m]
      double ystep=-1.0; ///< y-direction step size [m]
      double zmin = 999.0; ///< Minimum z coordinate (elevation) for output [m a.s.l.]
      double zmax =-999.0; ///< Maximum z coordinate (elevation) for output [m a.s.l.]
      double zstep=-1.0; ///< z-direction step size [m]
      int n_detector = 0; ///< Number of detectors
      bool output_binary = false; ///< If true, output in binary format; otherwise ASCII
    };

    /// @brief Output header information for cross-section output
    /// @param fout File pointer for output
    /// @param prm_zcross Cross-section parameters
    /// @note Writes metadata including grid dimensions, z-range, and number of detectors
    void out_cross_section_z_header(
      FILE *fout, const CrossSectionZParameters& prm_zcross) const;

    /// @brief Output multiple cross-sections at different elevations (ASCII format)
    /// @param path_out Output file path
    /// @param prm_zcross Cross-section parameters defining z-range and step size
    /// @throws std::runtime_error if file cannot be opened for writing
    /// @note Generates cross-sections from prm_zcross.zmin to zmax with zstep increments
    /// @note Output format: ASCII text with header followed by (x, y, density) data for each z-level
    void out_cross_section_z_all(
      const std::filesystem::path &path_out,
      const CrossSectionZParameters &prm_zcross) const;

    /// @brief Output multiple cross-sections at different elevations (binary format)
    /// @param path_out Output file path
    /// @param prm_zcross Cross-section parameters defining z-range and step size
    /// @throws std::runtime_error if file cannot be opened for writing
    /// @note Generates cross-sections from prm_zcross.zmin to zmax with zstep increments
    /// @note Output format: binary with header followed by double-precision (x, y, density) data
    void out_cross_section_z_all_binary(
      const std::filesystem::path &path_out,
      const CrossSectionZParameters &prm_zcross) const;

    /// @brief Output single cross-section at specified elevation
    /// @param pathout Output file path
    /// @param z_in Elevation at which to extract cross-section [m a.s.l.]
    /// @throws std::runtime_error if file cannot be opened for writing
    /// @note Output format: ASCII text with (x, y, density) for each grid cell
    void out_cross_section_z(const std::filesystem::path& pathout, const double z_in) const;

    /// @brief Generate cross-section of Grid2dPillar as Grid2dVoxel at specified elevation
    /// @param z Elevation for cross-section [m a.s.l.]
    /// @param n_detector Number of detectors (metadata for Grid2dVoxel)
    /// @return Grid2dVoxel representing the cross-section at elevation z
    /// @note Used for visualization or further analysis of density distribution at fixed elevation
    Grid2dVoxel make_cross_section_voxel(const double z, const int n_detector) const;

    ///@} ------------------------------------------------------------------

    //==================================================================
    /// @name add_structure_functions
    ///@{

    /// @brief Add vertical cylindrical density structure to existing Grid2dPillar
    /// @param prm Cylinder parameters (position, radius, density change, vertical range)
    /// @note Thread-safe: uses OpenMP parallelization internally
    /// @note Modifies density of cuboids within elliptical cylinder region
    /// @note If cuboid center is inside cylinder, delta_density is added to existing density
    /// @see VerticalCylinderParameters for parameter details
    void add_density_structure( const VerticalCylinderParameters &prm );

    /// @brief Add vertical dike (rectangular block) density structure
    /// @param prm Dike parameters (position, dimensions, rotation angle, density change, vertical range)
    /// @note Thread-safe: uses OpenMP parallelization internally
    /// @note Modifies density of cuboids within rotated rectangular region
    /// @note If cuboid center is inside dike rectangle, delta_density is added to existing density
    /// @see VerticalDikeParameters for parameter details
    void add_density_structure( const VerticalDikeParameters &prm );

    /// @brief Add vertical checkerboard density pattern
    /// @param prm Checkerboard parameters (origin, cell size, density variations, vertical range)
    /// @note Thread-safe: uses OpenMP parallelization internally
    /// @note Applies alternating density pattern based on grid cell position
    /// @note Useful for creating test patterns or synthetic density models
    /// @see VerticalCheckerBoardParameters for parameter details
    void add_density_structure( const VerticalCheckerBoardParameters &prm );
    ///@} ------------------------------------------------------------------

    //==================================================================
    /// @name path length calc function
    ///@{

    /// @brief Compute path length through cuboid intersected by ray
    /// @param tpl_hit_cub Grid indices (ix, iy) of the cuboid to test
    /// @param ray3d Ray to intersect with cuboid
    /// @return Path length through cuboid [m], or PATH_NO_HIT (0.0) if no intersection
    /// @details Primarily called during muon radiography path calculations.
    ///          Finds intersection between Ray3d and Pillar at (ix, iy).
    ///          Returns PATH_NO_HIT if no intersection, otherwise returns delta_path.
    /// @note Thread-safe: Yes (const method, read-only access)
    double get_delta_path(
      const Ixiy &tpl_hit_cub, const Ray3d &ray3d ) const;
    ///@} ------------------------------------------------------------------

    //==================================================================
    /// @name binary_io_functions
    ///@{

    /// @brief Save Grid2dPillar to binary stream.
    /// @param[out] ofs Output stream (must be opened in binary mode).
    void save(std::ofstream& ofs) const;

    /// @brief Load Grid2dPillar from binary stream.
    /// @param[in] ifs Input stream (must be opened in binary mode).
    void load(std::ifstream& ifs);

    /// @brief Save Grid2dPillar to file.
    /// @param[in] pathout Output file path.
    void save(const std::filesystem::path& pathout) const;

    /// @brief Load Grid2dPillar from file.
    /// @param[in] path_in Input file path.
    void load(const std::filesystem::path& path_in);

    ///@} ------------------------------------------------------------------
};
