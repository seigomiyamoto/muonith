/// @file cls_Grid3dVoxelParameters.hpp
/// @brief Parameter configuration classes for Grid3dVoxel
/// @details
/// Defines parameter classes for configuring 3D voxel grids and density structures:
/// - Grid3dVoxel::Parameters: Main configuration for voxel grid dimensions and hit detection
/// - Grid3dVoxel::CheckerBoard3dParameters: Checkerboard density pattern configuration
/// - Grid3dVoxel::EllipsoidParameters: Ellipsoidal density structure configuration
/// - Grid3dVoxel::CylinderParameters: Cylindrical density structure configuration
/// - Grid3dVoxel::ReconstVoxelsParameters: Reconstruction-region voxel masking configuration
/// - Grid3dVoxel::MergeParameters: Voxel merging operation configuration
///
/// ## Typical Usage Workflow
/// 1. Load Grid3dVoxel::Parameters from JSON configuration
/// 2. Optionally load structure parameters (checkerboard, ellipsoid, cylinder)
/// 3. Configure merging parameters if voxel coarsening is needed
/// 4. Pass parameters to Grid3dVoxel construction and methods
///
/// ## Coordinate System
/// - Right-handed coordinate system with z-up convention
/// - All spatial units are in meters (m) unless otherwise specified
/// - Elevation (z-axis) typically represents altitude in meters above sea level (m a.s.l.)
///
/// ## Thread Safety
/// - Read-only access is thread-safe
/// - Modification requires external synchronization
/// - assign_parameters() methods are NOT thread-safe
///
/// ## JSON Configuration Format
/// Parameters are loaded from JSON with hierarchical sections:
/// - Main section: "GRID3D_VOXEL_PARAMETERS" (customizable)
/// - Subsections: "checkerboard_3d_params", "ellipsoid_params", "cylinder_params", "merge_params", "reconst_voxels"
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
#include <nlohmann/json.hpp>

#include "cls_Ray.hpp"
#include "cls_AABB.hpp"

#include "ns_tuple_int.hpp"
#include "cls_Voxel.hpp"
#include "cls_Grid3d.hpp"
#include "cls_Grid3dVoxel.hpp"
#include "cls_Grid2dPillar.hpp"
#include "cls_VoxelUniqueIndexMapContainer.hpp"
#include "ns_angle_util.hpp"

//#####################################################################
/// @class Grid3dVoxel::CheckerBoard3dParameters
/// @brief Parameter container for passing checkerboard structure configuration to @ref Grid3dVoxel::add_density_structure(const CheckerBoard3dParameters &)
/// @details
/// Configures a 3D checkerboard density pattern within the voxel grid. The pattern alternates
/// density values in a regular grid pattern, useful for testing and validation scenarios.
///
/// ## Usage Example
/// @code
/// nlohmann::json js = /* load from file */;
/// Grid3dVoxel::CheckerBoard3dParameters prm(js);
/// grid3d_voxel.add_density_structure(prm);
/// @endcode
///
/// ## Density Application
/// For each voxel in the checkerboard pattern:
/// - density += delta_density_offset + delta_density
///
/// ## Thread Safety
/// - assign_parameters() is NOT thread-safe
/// - Read-only access after initialization is thread-safe
///
/// @ingroup parameterClasses
//#####################################################################
class Grid3dVoxel::CheckerBoard3dParameters {
  private:
  public:
    //======================================================================
    /// @name public member variables
    ///@{
    /// @brief Name identifier for this checkerboard pattern
    std::string name = "none";

    /// @brief Enable flag: if true, execute Grid3dVoxel::add_structure
    bool tf_exec = false;

    /// @brief Density offset applied before delta_density (kg/m^3)
    double delta_density_offset = 0.0;

    /// @brief Primary density change applied to pattern voxels (kg/m^3)
    /// @details Final density = density + delta_density_offset + delta_density
    double delta_density = 0.0;

    /// @brief Axis-aligned bounding box defining the checkerboard region (m)
    AABB3d aabb3d = AABB3d();

    /// @brief Center position of the checkerboard structure (m)
    Eigen::Vector3d v3_pos_cnt = Eigen::Vector3d(0.0, 0.0, 0.0);

    /// @brief Checkerboard cell size as interval multipliers (must be positive integers)
    /// @details Actual cell dimensions are computed at runtime as:
    ///   v3_length_computed.x() = grid.get_x_interval() * v3_len_interval_mult.x()
    /// @note Units: dimensionless (grid interval count)
    Eigen::Vector3i v3_len_interval_mult = Eigen::Vector3i(4, 4, 4);

    /// @brief Computed cell dimensions in meters (populated by add_density_structure_all)
    /// @note Internal use only; computed from v3_len_interval_mult and grid intervals
    Eigen::Vector3d v3_length_computed = Eigen::Vector3d(0.0, 0.0, 0.0);

    /// @brief If true, snap aabb3d boundaries and v3_pos_cnt to grid boundaries/centers
    /// @details When enabled, coordinates are automatically aligned to grid for better pattern integrity
    bool tf_snap_to_grid = true;

    /// @brief Total number of checkerboard cells per axis (alternative to AABB)
    /// @details When specified, AABB is computed symmetrically from center:
    ///   xmin = xcnt_snapped - xlen_cells * 0.5 * cell_size
    ///   xmax = xcnt_snapped + xlen_cells * 0.5 * cell_size
    ///   where cell_size = get_x_interval() * xlen_interval_mult
    ///   Total cells per axis = xlen_cells
    /// @note If set to (0,0,0), AABB from JSON is used directly (backward compatible)
    /// @note Odd values are supported (e.g. xlen_cells=5 gives half-extent of 2.5 cells)
    Eigen::Vector3i v3_len_cells = Eigen::Vector3i(0, 0, 0);

    /// @brief Region type for checkerboard pattern application
    /// @details "aabb": axis-aligned bounding box (default),
    ///          "cylinder": vertical cylinder (XY circle x Z range)
    std::string region_type = "aabb";

    /// @brief Elliptic cylinder radii in meters
    /// @details Only used when region_type == "cylinder".
    /// @note Units: meters
    double radius_x_meters = 0.0;
    double radius_y_meters = 0.0;

    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name constructor_destructor
    ///@{
    /// @brief default constructor
    CheckerBoard3dParameters() = default;

    /// @brief copy constructor
    CheckerBoard3dParameters(const CheckerBoard3dParameters& org) = default;

    /// @brief move constructor
    CheckerBoard3dParameters(CheckerBoard3dParameters&& other) noexcept = default;

    /// @brief build from nlohmann::json
    CheckerBoard3dParameters(const nlohmann::json& js )
      : CheckerBoard3dParameters() { assign_parameters(js); };

    /// @brief destructor
    ~CheckerBoard3dParameters() = default;
    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name operators
    ///@{

    /// @brief assign operator
    CheckerBoard3dParameters& operator=(const CheckerBoard3dParameters& other) = default;

    /// @brief not equal operator
    /// @note Does not compare name field.
    bool operator!=(const CheckerBoard3dParameters& other) const;

    /// @brief equal operator
    /// @note Does not compare name field. Defined using inequality operator.
    bool operator==(const CheckerBoard3dParameters& other) const {
      return !(*this != other);
    };
    ///@} ------------------------------------------------------------------

    /// @brief Get the name identifier
    /// @return Name string of this checkerboard pattern
    std::string get_name() const { return name; };

    /// @brief Load parameters from JSON configuration
    /// @param js JSON object containing checkerboard configuration
    /// @throws std::exception if required keys are missing or values are invalid
    /// @note Thread-safe: No
    void assign_parameters(const nlohmann::json& js);
};

//#####################################################################
/// @class Grid3dVoxel::EllipsoidParameters
/// @brief Parameter container for passing ellipsoid structure configuration to @ref Grid3dVoxel::add_density_structure(const EllipsoidParameters&)
/// @details
/// Configures an ellipsoidal density anomaly within the voxel grid. The ellipsoid can be
/// rotated using Euler angles with either local (intrinsic) or global (extrinsic) rotation conventions.
///
/// ## Usage Example
/// @code
/// nlohmann::json js = /* load from file */;
/// Grid3dVoxel::EllipsoidParameters prm(js);
/// grid3d_voxel.add_density_structure(prm);
/// @endcode
///
/// ## Density Application
/// For voxels inside the ellipsoid: density += delta_density
///
/// ## Rotation Convention
/// - LOCAL: Intrinsic rotations (rotations applied in body frame)
/// - GLOBAL: Extrinsic rotations (rotations applied in fixed frame)
///
/// ## Thread Safety
/// - assign_parameters() is NOT thread-safe
/// - Read-only access after initialization is thread-safe
///
/// @ingroup parameterClasses
//#####################################################################
class Grid3dVoxel::EllipsoidParameters {
  private:
  public:
    //======================================================================
    /// @name public member variables
    ///@{
    /// @brief Enable flag: if true, execute Grid3dVoxel::add_structure
    bool tf_exec = false;

    /// @brief Name identifier for this ellipsoid structure
    std::string name = "none";

    /// @brief Density change applied to ellipsoid interior voxels (kg/m^3)
    /// @details For voxels inside: density += delta_density
    double delta_density = 0.0;

    /// @brief Center position of the ellipsoid (m)
    Eigen::Vector3d v3_pos_cnt = Eigen::Vector3d(0.0, 0.0, 0.0);

    /// @brief Semi-axis lengths in x, y, z directions before rotation (m)
    Eigen::Vector3d v3_length = Eigen::Vector3d(0.0, 0.0, 0.0);

    /// @brief Rotation angles around x, y, z axes (degrees)
    /// @details Applied according to rotation_type (LOCAL or GLOBAL)
    Angle theta_x = Angle(0.0, Angle::Unit::Degree),
          theta_y = Angle(0.0, Angle::Unit::Degree),
          theta_z = Angle(0.0, Angle::Unit::Degree);

    /// @brief Rotation convention: LOCAL (intrinsic) or GLOBAL (extrinsic)
    angle_util::Rotation3dType rotation_type = angle_util::Rotation3dType::LOCAL;

    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name constructor_destructor
    ///@{
    /// @brief default constructor
    EllipsoidParameters() = default;
  
    /// @brief copy constructor
    EllipsoidParameters(const EllipsoidParameters& org) = default;
  
    /// @brief build from nlohmann::json
    EllipsoidParameters(const nlohmann::json& js )
      : EllipsoidParameters() { assign_parameters(js); }
  
    /// @brief destructor
    ~EllipsoidParameters() = default;
    ///@} ------------------------------------------------------------------
  
    //======================================================================
    /// @name operators
    ///@{
    /// @brief assignment operator
    EllipsoidParameters& operator=(const EllipsoidParameters& other) = default;

    /// @brief not equal operator
    /// @note Does not compare name field.
    bool operator!=(const EllipsoidParameters& other) const;

    /// @brief equal operator
    /// @note Does not compare name field. Defined using inequality operator.
    bool operator==(const EllipsoidParameters& other) const {
      return !(*this != other);
    };
    ///@} ------------------------------------------------------------------
  
    /// @brief Get the name identifier
    /// @return Name string of this ellipsoid structure
    std::string get_name() const { return name; };

    /// @brief Load parameters from JSON configuration
    /// @param js JSON object containing ellipsoid configuration
    /// @throws std::exception if required keys are missing or rotation_type is invalid
    /// @note Thread-safe: No
    void assign_parameters(const nlohmann::json& js );
};

//#####################################################################
  /// @class Grid3dVoxel::CylinderParameters
  /// @brief Parameter container for passing cylinder structure configuration to @ref Grid3dVoxel::add_density_structure(const CylinderParameters &)
  /// @details
  /// Configures a cylindrical (or elliptical cylinder) density anomaly within the voxel grid.
  /// The cylinder has an elliptical cross-section defined by x_length and y_length, with height z_length.
  /// Rotation is supported using Euler angles with LOCAL or GLOBAL conventions.
  ///
  /// ## Usage Example
  /// @code
  /// nlohmann::json js = /* load from file */;
  /// Grid3dVoxel::CylinderParameters prm(js);
  /// grid3d_voxel.add_density_structure(prm);
  /// @endcode
  ///
  /// ## Geometry
  /// - Cross-section: Ellipse with semi-axes x_length and y_length in the xy-plane
  /// - Height: z_length, centered at v3_pos_cnt.z() ± z_length/2
  ///
  /// ## Density Application
  /// For voxels inside the cylinder: density += delta_density
  ///
  /// ## Thread Safety
  /// - assign_parameters() is NOT thread-safe
  /// - Read-only access after initialization is thread-safe
  ///
  /// @ingroup parameterClasses
//####################################################################
class Grid3dVoxel::CylinderParameters {
  private:
  public:
    //======================================================================
    /// @name public member variables
    ///@{
    /// @brief Enable flag: if true, execute Grid3dVoxel::add_structure
    bool tf_exec = false;

    /// @brief Name identifier for this cylinder structure
    std::string name = "none";

    /// @brief Density change applied to cylinder interior voxels (kg/m^3)
    /// @details For voxels inside: density += delta_density
    double delta_density = 0.0;

    /// @brief Center position of the cylinder (m)
    Eigen::Vector3d v3_pos_cnt = Eigen::Vector3d(0.0, 0.0, 0.0);

    /// @brief Cylinder dimensions: x_length and y_length are elliptical cross-section semi-axes, z_length is height (m)
    /// @details The cylinder extends from v3_pos_cnt.z() ± v3_length.z() * 0.5 in the z-direction.
    Eigen::Vector3d v3_length = Eigen::Vector3d(0.0, 0.0, 0.0);

    /// @brief Rotation angles around x, y, z axes (degrees)
    /// @details Applied according to rotation_type (LOCAL or GLOBAL)
    Angle theta_x = Angle(0.0, Angle::Unit::Degree),
          theta_y = Angle(0.0, Angle::Unit::Degree),
          theta_z = Angle(0.0, Angle::Unit::Degree);

    /// @brief Rotation convention: LOCAL (intrinsic) or GLOBAL (extrinsic)
    angle_util::Rotation3dType rotation_type = angle_util::Rotation3dType::LOCAL;

    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name constructor_destructor
    ///@{
    /// @brief default constructor
    CylinderParameters() = default;

    /// @brief copy constructor
    CylinderParameters(const CylinderParameters& org) = default;

    /// @brief build from nlohmann::json
    CylinderParameters(const nlohmann::json& js)
      : CylinderParameters() { assign_parameters(js); };

    /// @brief destructor
    ~CylinderParameters() = default;
    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name operators
    ///@{
    /// @brief assignment operator
    CylinderParameters& operator=(const CylinderParameters& other) = default;

    /// @brief not equal operator
    /// @note Does not compare name field.
    bool operator!=(const CylinderParameters& other) const;

    /// @brief equal operator
    /// @note Does not compare name field. Defined using inequality operator.
    bool operator==(const CylinderParameters& other) const {
      return !(*this != other);
    };
    ///@} ------------------------------------------------------------------

    /// @brief Get the name identifier
    /// @return Name string of this cylinder structure
    std::string get_name() const { return name; };

    /// @brief Load parameters from JSON configuration
    /// @param js JSON object containing cylinder configuration
    /// @throws std::exception if required keys are missing, lengths are non-positive, or rotation_type is invalid
    /// @note Thread-safe: No
    void assign_parameters(const nlohmann::json& js );
};

//#####################################################################
  /// @class Grid3dVoxel::CuboidParameters
  /// @brief Parameter container for passing cuboid (rectangular box) structure configuration to @ref Grid3dVoxel::add_density_structure(const CuboidParameters &)
  /// @details
  /// Configures a rectangular-box (cuboid) density anomaly within the voxel grid.
  /// The box has full side lengths x_length, y_length, z_length centered at v3_pos_cnt.
  /// Rotation is supported using Euler angles with LOCAL or GLOBAL conventions.
  ///
  /// ## Density Application
  /// For voxels inside the cuboid: density += delta_density
  ///
  /// ## Thread Safety
  /// - assign_parameters() is NOT thread-safe
  /// - Read-only access after initialization is thread-safe
  ///
  /// @ingroup parameterClasses
//####################################################################
class Grid3dVoxel::CuboidParameters {
  private:
  public:
    //======================================================================
    /// @name public member variables
    ///@{
    /// @brief Enable flag: if true, execute Grid3dVoxel::add_structure
    bool tf_exec = false;

    /// @brief Name identifier for this cuboid structure
    std::string name = "none";

    /// @brief Density change applied to cuboid interior voxels (kg/m^3)
    /// @details For voxels inside: density += delta_density
    double delta_density = 0.0;

    /// @brief Center position of the cuboid (m)
    Eigen::Vector3d v3_pos_cnt = Eigen::Vector3d(0.0, 0.0, 0.0);

    /// @brief Cuboid full side lengths along x, y, z (m)
    /// @details Half-extent along each axis is 0.5 * v3_length.
    Eigen::Vector3d v3_length = Eigen::Vector3d(0.0, 0.0, 0.0);

    /// @brief Rotation angles around x, y, z axes (degrees)
    /// @details Applied according to rotation_type (LOCAL or GLOBAL)
    Angle theta_x = Angle(0.0, Angle::Unit::Degree),
          theta_y = Angle(0.0, Angle::Unit::Degree),
          theta_z = Angle(0.0, Angle::Unit::Degree);

    /// @brief Rotation convention: LOCAL (intrinsic) or GLOBAL (extrinsic)
    angle_util::Rotation3dType rotation_type = angle_util::Rotation3dType::LOCAL;
    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name constructor_destructor
    ///@{
    /// @brief default constructor
    CuboidParameters() = default;

    /// @brief copy constructor
    CuboidParameters(const CuboidParameters& org) = default;

    /// @brief build from nlohmann::json
    CuboidParameters(const nlohmann::json& js)
      : CuboidParameters() { assign_parameters(js); };

    /// @brief destructor
    ~CuboidParameters() = default;
    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name operators
    ///@{
    /// @brief assignment operator
    CuboidParameters& operator=(const CuboidParameters& other) = default;

    /// @brief not equal operator
    /// @note Does not compare name field.
    bool operator!=(const CuboidParameters& other) const;

    /// @brief equal operator
    /// @note Does not compare name field. Defined using inequality operator.
    bool operator==(const CuboidParameters& other) const {
      return !(*this != other);
    };
    ///@} ------------------------------------------------------------------

    /// @brief Get the name identifier
    /// @return Name string of this cuboid structure
    std::string get_name() const { return name; };

    /// @brief Load parameters from JSON configuration
    /// @param js JSON object containing cuboid configuration
    /// @throws std::exception if required keys are missing, lengths are non-positive, or rotation_type is invalid
    /// @note Thread-safe: No
    void assign_parameters(const nlohmann::json& js );
};

//#####################################################################
  /// @class Grid3dVoxel::ReconstVoxelsParameters
  /// @brief Parameter container for reconstruction-region masking of post-merge voxels
  /// @details
  /// Defines the spatial mask applied to merged voxels to select the reconstruction region.
  /// Two mask shapes are supported: axis-aligned bounding box (AABB) and vertical elliptic
  /// cylinder. Both are mutually exclusive and specified in post-merge voxel-cell counts
  /// centered at the merge center.
  ///
  /// ## Usage Example
  /// @code
  /// nlohmann::json js = /* load from file */;
  /// Grid3dVoxel::ReconstVoxelsParameters prm(js);
  /// @endcode
  ///
  /// ## Thread Safety
  /// - assign_parameters() is NOT thread-safe
  /// - Read-only access after initialization is thread-safe
  ///
  /// @ingroup parameterClasses
//#####################################################################
class Grid3dVoxel::ReconstVoxelsParameters {
  public:
    //======================================================================
    /// @name public member variables
    ///@{
    /// @brief Enable AABB mask after merge.
    /// @details If true, voxels outside the specified AABB are set to
    ///          tf_exist=false after the merge operation. The AABB is defined
    ///          by cell counts in post-merge voxel units, centered at merge_center.
    bool tf_aabb = false;

    /// @brief X center of the AABB mask (m).
    /// @note Units: meters (same coordinate system as merge_center)
    double x_aabb_cnt = 0.0;

    /// @brief Y center of the AABB mask (m).
    /// @note Units: meters
    double y_aabb_cnt = 0.0;

    /// @brief Full X extent of the AABB mask (m).
    /// @details AABB half-extent = x_aabb_meters * 0.5. Center is x_aabb_cnt.
    /// @note Units: meters
    double x_aabb_meters = 0.0;

    /// @brief Full Y extent of the AABB mask (m).
    /// @details AABB half-extent = y_aabb_meters * 0.5. Center is y_aabb_cnt.
    /// @note Units: meters
    double y_aabb_meters = 0.0;

    /// @brief Z-min determination mode for AABB: "g3vox_zmin" or "manual".
    /// @details "g3vox_zmin": use merged grid zmin (no lower shell gap).
    ///          "manual": use aabb_zmin_value (Phase 2 shell fix required).
    std::string aabb_zmin_mode = "g3vox_zmin";

    /// @brief Explicit z-min value (m). Used only when aabb_zmin_mode == "manual".
    /// @note Units: meters
    double aabb_zmin_value = 0.0;

    /// @brief Z-max of the AABB (m). Voxels above this are masked.
    /// @details If terrain surface is lower than aabb_zmax, terrain prevails
    ///          (tf_exist from DEM is already false above terrain).
    /// @note Units: meters
    double aabb_zmax = std::numeric_limits<double>::max();

    /// @brief Enable vertical elliptic cylinder mask after merge.
    /// @details If true, voxels outside the cylinder are set to tf_exist=false.
    ///          The cylinder is centered at merge_center in the XY plane.
    bool tf_cylinder = false;

    /// @brief X center of the cylinder mask (m).
    /// @note Units: meters
    double x_cyl_cnt = 0.0;

    /// @brief Y center of the cylinder mask (m).
    /// @note Units: meters
    double y_cyl_cnt = 0.0;

    /// @brief Cylinder radius in x-direction (m).
    /// @details Defines the semi-axis of the elliptic cylinder in X.
    /// @note Units: meters. Must be > 0 when tf_cylinder is true.
    double cylinder_radius_x_meters = 0.0;

    /// @brief Cylinder radius in y-direction (m).
    /// @details Defines the semi-axis of the elliptic cylinder in Y.
    /// @note Units: meters. Must be > 0 when tf_cylinder is true.
    double cylinder_radius_y_meters = 0.0;
    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name constructor_destructor
    ///@{
    /// @brief default constructor
    ReconstVoxelsParameters() = default;

    /// @brief copy constructor
    ReconstVoxelsParameters(const ReconstVoxelsParameters& org) = default;

    /// @brief build from nlohmann::json
    ReconstVoxelsParameters(const nlohmann::json& js)
      : ReconstVoxelsParameters() { assign_parameters(js); }

    /// @brief destructor
    ~ReconstVoxelsParameters() = default;
    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name operators
    ///@{
    /// @brief assignment operator
    ReconstVoxelsParameters& operator=(const ReconstVoxelsParameters& other) = default;

    /// @brief not equal operator
    bool operator!=(const ReconstVoxelsParameters& other) const;

    /// @brief equal operator. Defined using inequality operator.
    bool operator==(const ReconstVoxelsParameters& other) const {
      return !(*this != other);
    };
    ///@} ------------------------------------------------------------------

    /// @brief Load parameters from JSON configuration
    /// @param js JSON object containing reconst_voxels configuration
    /// @throws std::runtime_error if validation fails (e.g. cylinder radii <= 0 when enabled)
    /// @note Thread-safe: No
    void assign_parameters(const nlohmann::json& js);
};

//#####################################################################
  /// @class Grid3dVoxel::MergeParameters
  /// @brief Parameter container for voxel merging operations passed to Grid3dVoxel::merge
  /// @details
  /// Configures how voxels are merged (coarsened) in the grid. Merging combines multiple
  /// adjacent voxels into a single larger voxel, reducing memory and computation costs.
  ///
  /// ## Usage Example
  /// @code
  /// nlohmann::json js = /* load from file */;
  /// Grid3dVoxel::MergeParameters prm(js);
  /// grid3d_voxel.merge(prm);
  /// @endcode
  ///
  /// ## Merge Operation
  /// - Merge factors define how many voxels to combine in each direction
  /// - Merge center defines the alignment point for the merging grid
  /// - All merge factors must be ≥ 2
  ///
  /// ## Thread Safety
  /// - assign_parameters() is NOT thread-safe
  /// - Read-only access after initialization is thread-safe
  ///
  /// @ingroup parameterClasses
//#####################################################################
class Grid3dVoxel::MergeParameters {
  public:
    //======================================================================
    /// @name public member variables
    ///@{
    /// @brief Enable flag: if true, execute Grid3dVoxel::merge
    bool tf_exec = false;

    /// @brief Name identifier for this merge configuration
    std::string name = "GRID3D_VOXEL_MERGE_PARAMETERS";

    /// @brief Merge grid alignment center coordinates (m)
    double x_merge_center = 0.0, y_merge_center = 0.0, z_merge_center = 0.0;

    /// @brief Merge factors: number of voxels to combine in each direction (must be ≥ 2)
    int x_merge_factor = 1, y_merge_factor = 1, z_merge_factor = 1;
    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name constructor_destructor
    ///@{
    /// @brief default constructor
    MergeParameters() = default;

    /// @brief copy constructor
    MergeParameters(const MergeParameters& org) = default;

    /// @brief build from nlohmann::json
    MergeParameters(const nlohmann::json& js)
      : MergeParameters() { assign_parameters(js); }

    /// @brief destructor
    ~MergeParameters() = default;
    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name operators
    ///@{
    /// @brief assignment operator
    MergeParameters& operator=(const MergeParameters& other) = default;

    /// @brief not equal operator
    /// @note Does not compare name field.
    bool operator!=(const MergeParameters& other) const;

    /// @brief equal operator
    /// @note Does not compare name field. Defined using inequality operator.
    bool operator==(const MergeParameters& other) const {
      return !(*this != other);
    };
    ///@} ------------------------------------------------------------------

    /// @brief Get the name identifier
    /// @return Name string of this merge configuration
    std::string get_name() const { return name; };

    /// @brief Load parameters from JSON configuration
    /// @param js JSON object containing merge configuration
    /// @throws std::exception if merge factors are < 2
    /// @note Thread-safe: No
    void assign_parameters(const nlohmann::json& js );
};

//#####################################################################
/// @class Grid3dVoxel::Parameters
/// @brief Main parameter configuration for Grid3dVoxel
/// @details
/// Configures the 3D voxel grid dimensions, hit detection criteria, and embedded structure parameters.
/// Supports loading from JSON with customizable section names to allow multiple configurations.
///
/// ## Usage Example
/// @code
/// nlohmann::json js = /* load from file */;
/// Grid3dVoxel::Parameters prm(js, "GRID3D_VOXEL_PARAMETERS");
/// Grid3dVoxel grid(prm);
/// @endcode
///
/// ## Grid Configuration
/// - x-axis: [xmin, xmax] with pitch x_pitch
/// - y-axis: [ymin, ymax] with pitch y_pitch
/// - z-axis: [zmin, zmax] with pitch z_pitch
/// - Can optionally inherit x/y settings from Grid2dPillar via use_grid2d_of_g2pil flag
///
/// ## Hit Detection Filtering
/// - n_hit_det_min/max: Valid range for detector panel hit counts
/// - n_hit_ele_min/max: Valid range for detector element hit counts
///
/// ## JSON Section Names
/// - Default main section: "GRID3D_VOXEL_PARAMETERS" (customizable via constructor)
/// - Subsections for structures are fixed (see section_name_* constants)
///
/// ## Thread Safety
/// - assign_parameters() is NOT thread-safe
/// - Read-only access after initialization is thread-safe
///
/// @note "GRID3D_VOXEL_PARAMETERS" is the default JSON key, but you can specify
///       custom section names in the constructor to create multiple configurations.
/// @ingroup parameterClasses
//#####################################################################
class Grid3dVoxel::Parameters {
  private:
  public:
    //======================================================================
    /// @name public member variables
    ///@{
    /// @brief Large integer constant for unbounded maximum values
    static constexpr int INT_LARGE = std::numeric_limits<int>::max();

    /// @brief Name identifier for this parameter instance
    std::string name = "Grid3dVoxel::Parameters";

    /// @brief X-axis grid boundaries and spacing (m)
    double xmin = 0, xmax = 100, x_pitch = 10;

    /// @brief Y-axis grid boundaries and spacing (m)
    double ymin = 0, ymax = 100, y_pitch = 10;

    /// @brief Z-axis grid boundaries and spacing (m a.s.l.)
    double zmin = 0, zmax = 100, z_pitch = 10;
    
    /// @brief Valid range for n_hit_det (number of DetectorPanels that hit a voxel)
    int n_hit_det_min = 0, n_hit_det_max = INT_LARGE;

    /// @brief Valid range for n_hit_ele (number of DetectorElements that hit a voxel)
    int n_hit_ele_min = 0, n_hit_ele_max = INT_LARGE;
  
    /// @brief If true, terminate program after building merged g3vox_input
    bool tf_end_after_merged = false;

    /// @brief If true, inherit x/y grid settings from Grid2dPillar instead of JSON
    bool use_grid2d_of_g2pil = false;

    /// @brief If false, skip Grid3dVoxel construction in execution flow
    bool tf_build_g3vox = true;

    //============================
    /// @brief JSON section name for CheckerBoard3dParameters
    static constexpr const char*
      section_name_checkerboard = "checkerboard_3d_params";

    /// @brief JSON section name for EllipsoidParameters
    static constexpr const char*
      section_name_ellipsoid = "ellipsoid_params";

    /// @brief JSON section name for CylinderParameters
    static constexpr const char*
      section_name_cylinder = "cylinder_params";

    /// @brief JSON section name for CuboidParameters
    static constexpr const char*
      section_name_cuboid = "cuboid_params";

    /// @brief JSON section name for MergeParameters
    static constexpr const char*
      section_name_merge = "merge_params";

    /// @brief JSON section name for ReconstVoxelsParameters
    static constexpr const char*
      section_name_reconst_voxels = "reconst_voxels";

    //============================
    /// @brief Parameter storage for checkerboard structures
    std::vector<CheckerBoard3dParameters> vec_checkerboard_prm = {};

    /// @brief Parameter storage for ellipsoid structures
    std::vector<EllipsoidParameters> vec_ellipsoid_prm = {};

    /// @brief Parameter storage for cylinder structures
    std::vector<CylinderParameters> vec_cylinder_prm = {};

    /// @brief Parameter storage for cuboid structures
    std::vector<CuboidParameters> vec_cuboid_prm = {};

    //============================
    /// @brief Parameter storage for voxel merging configuration
    Grid3dVoxel::MergeParameters prm_merge;

    /// @brief Parameter storage for reconstruction-region voxel masking
    Grid3dVoxel::ReconstVoxelsParameters prm_reconst_voxels;
    ///@} ------------------------------------------------------------------
    
    //======================================================================
    /// @name constructor_destructor
    ///@{
    /// @brief default constructor
    Parameters() = default;
    
    /// @brief copy constructor
    Parameters(const Parameters& org) = default;
    
    /// @brief constructor from nlohmann::json
    Parameters(const nlohmann::json& js, const std::string &section_name = "GRID3D_VOXEL_PARAMETERS")
      : Parameters() { assign_parameters(js, section_name); };
    
    /// @brief destructor
    ~Parameters() = default;
    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name operators
    ///@{
    /// @brief assignment operator
    Parameters& operator=(const Parameters& other) = default;

    /// @brief not equal operator
    /// @note Does not compare name field.
    bool operator!=(const Parameters& other) const;

    /// @brief equal operator
    /// @note Does not compare name field. Defined using inequality operator.
    bool operator==(const Parameters& other) const {
      return !(*this != other);
    };
    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name getter_functions
    ///@{
    /// @brief Get the name identifier
    /// @return Name string of this parameter instance
    std::string get_name() const { return name; };

    /// @brief Get number of voxels in x-direction
    /// @return Number of bins along x-axis
    int get_nbinx() const { return static_cast<int>((xmax - xmin) / x_pitch); };

    /// @brief Get number of voxels in y-direction
    /// @return Number of bins along y-axis
    int get_nbiny() const { return static_cast<int>((ymax - ymin) / y_pitch); };

    /// @brief Get number of voxels in z-direction
    /// @return Number of bins along z-axis
    int get_nbinz() const { return static_cast<int>((zmax - zmin) / z_pitch); };
    ///@} ------------------------------------------------------------------

    /// @brief Load parameters from JSON configuration
    /// @param js JSON object containing the full configuration
    /// @param section_name Section name within JSON (default: "GRID3D_VOXEL_PARAMETERS")
    /// @throws std::exception if required keys are missing or values violate constraints
    /// @note Thread-safe: No
    /// @note Validates: min < max for all axes, pitch > 0, n_hit_min < n_hit_max
    void assign_parameters(const nlohmann::json& js, const std::string &section_name);
};
