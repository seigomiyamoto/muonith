/// @file cls_Grid2dPillarParameters.hpp
/// @brief Parameter configuration for Grid2dPillar density structures
/// @details
/// Defines parameter classes for configuring Grid2dPillar objects and their associated
/// vertical density structures (cylinders, dikes, checkerboards).
///
/// ## Typical workflow
/// 1. Create a `Grid2dPillar::Parameters` instance from JSON configuration
/// 2. Set DEM file path and initial uniform density (kg/m³)
/// 3. Optionally add vertical structure parameters (cylinder, dike, checkerboard)
/// 4. Pass parameters to Grid2dPillar constructor
/// 5. Grid2dPillar applies DEM topography and adds density structures
///
/// ## Coordinate system
/// - Right-handed coordinate system with z-up (vertical axis points upward)
/// - x, y: horizontal coordinates (m)
/// - z: elevation above sea level (m a.s.l.)
///
/// ## Units
/// - Spatial dimensions: meters (m)
/// - Density: kilograms per cubic meter (kg/m³)
/// - Angles: degrees (input), converted internally to radians
///
/// ## Thread safety
/// - Parameter objects are NOT thread-safe during modification (assign_parameters)
/// - Read-only access after construction is safe for concurrent reads
/// - JSON parsing is single-threaded
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

#include "cls_AABB.hpp"
#include "cls_Angle.hpp"
#include "cls_Pillar.hpp"
#include "cls_Grid2d.hpp"
#include "cls_Grid2dPillar.hpp"

//##################################################################
/// @class Grid2dPillar::Parameters
/// @brief parameters to build Grid2dPillar
/// @ingroup parameterClasses
//##################################################################

//################################################################################
//################################################################################
/// @class Grid2dPillar::VerticalCylinderParameters
/// @brief Parameters for vertical cylinder structure used in Grid2dPillar::add_density_structure
/// @details
/// Configures a vertical cylinder (elliptical cross-section) with delta density anomaly.
/// The cylinder extends infinitely in the vertical (z) direction and is defined by:
/// - Center position (x, y)
/// - Semi-axes lengths (x_length, y_length)
/// - Rotation angle around z-axis
/// - Delta density (kg/m³) added to the background density
///
/// ### Typical use case
/// @code
/// Grid2dPillar::VerticalCylinderParameters cyl_prm;
/// cyl_prm.name = "magma_conduit";
/// cyl_prm.delta_density = -200.0;  // kg/m³ (less dense than background)
/// cyl_prm.v2_pos_cnt = Eigen::Vector2d(500.0, 600.0);  // m
/// cyl_prm.v2_length = Eigen::Vector2d(50.0, 50.0);     // m (semi-axes)
/// cyl_prm.angle_rot.setDegree(45.0);                  // rotation angle
/// cyl_prm.tf_exec = true;
/// @endcode
///
/// @note Thread-safe: No. Modification is not thread-safe.
/// @note Coordinate system: Right-handed, z-up
/// @note Units: meters (m) for lengths, kg/m³ for density
//################################################################################
//################################################################################
class Grid2dPillar::VerticalCylinderParameters {
  private:
  public:
    //==================================================================
    /// @name public member variables
    //@{
    /// @brief if true, Grid2dPillar::add_structure will apply this cylinder structure
    bool tf_exec = false;

    /// @brief name identifier for this cylinder structure
    std::string name = "none";

    /// @brief delta density anomaly added to background density (kg/m³)
    /// @note Positive values increase density, negative values decrease density
    double delta_density = 0.0;

    /// @brief center position (x, y) of the cylinder in horizontal plane (m)
    Eigen::Vector2d v2_pos_cnt = Eigen::Vector2d(0.0, 0.0);

    /// @brief semi-axes lengths (x_length, y_length) defining elliptical cross-section (m)
    /// @note For circular cylinder, set x_length == y_length
    Eigen::Vector2d v2_length = Eigen::Vector2d(0.0, 0.0);

    /// @brief rotation angle of the cylinder around z-axis (counterclockwise from x-axis)
    /// @note Input in degrees, stored internally. Positive rotation follows right-hand rule around z-up.
    Angle angle_rot = Angle(0.0, Angle::Unit::Degree);
    //@} ------------------------------------------------------------------

    //==================================================================
    /// @name constructor_destructor
    //@{
    /// @brief default constructor
    /// @note All members initialized to default values (tf_exec=false)
    VerticalCylinderParameters() = default;

    /// @brief copy constructor
    /// @param org source object to copy from
    VerticalCylinderParameters(const VerticalCylinderParameters& org) = default;

    /// @brief constructor from nlohmann::json configuration
    /// @param js JSON object containing configuration
    /// @param section_name section key within JSON to read from
    /// @note Calls assign_parameters internally
    VerticalCylinderParameters(const nlohmann::json& js, const std::string &section_name)
      : VerticalCylinderParameters() { assign_parameters(js, section_name); };

    /// @brief constructor from explicit parameters
    /// @param name_in identifier name for this cylinder
    /// @param delta_density_in density anomaly (kg/m³)
    /// @param v2_pos_cnt_in center position (x, y) in m
    /// @param v2_length_in semi-axes lengths (x_length, y_length) in m
    /// @param angle_rot_in rotation angle around z-axis
    /// @note Automatically sets tf_exec to true
    VerticalCylinderParameters(const std::string &name_in,
                               const double delta_density_in,
                               const Eigen::Vector2d &v2_pos_cnt_in,
                               const Eigen::Vector2d &v2_length_in,
                               const Angle &angle_rot_in)
      : tf_exec(true),
        name(name_in),
        delta_density(delta_density_in),
        v2_pos_cnt(v2_pos_cnt_in),
        v2_length(v2_length_in),
        angle_rot(angle_rot_in) {};

    /// @brief destructor
    ~VerticalCylinderParameters() = default;
    //@} ------------------------------------------------------------------

    //==================================================================
    /// @name operators
    //@{
    /// @brief assignment operator
    /// @param other source object to assign from
    /// @return reference to this object
    VerticalCylinderParameters& operator=(const VerticalCylinderParameters& other) = default;
    //@} ------------------------------------------------------------------

    /// @brief get name of this instance
    /// @return name identifier string
    std::string get_name() const { return name; };

    /// @brief load and assign parameters from JSON configuration
    /// @param js JSON object containing configuration
    /// @param section_name section key within JSON (used for logging only in nested objects)
    /// @note Missing keys in JSON will trigger warnings and use default values
    /// @note Thread-safe: No
    void assign_parameters(const nlohmann::json& js, const std::string &section_name);
};

//################################################################################
//################################################################################
/// @class Grid2dPillar::VerticalDikeParameters
/// @brief Parameters for vertical dike structure used in Grid2dPillar::add_density_structure
/// @details
/// Configures a vertical dike (rectangular prism extending infinitely in z-direction).
/// Defined by an axis-aligned bounding box (AABB) in the horizontal plane and rotation angle.
///
/// @note Thread-safe: No. Modification is not thread-safe.
/// @note Coordinate system: Right-handed, z-up
/// @note Units: meters (m) for lengths, kg/m³ for density
//################################################################################
//################################################################################
class Grid2dPillar::VerticalDikeParameters {
  private:
  public:
    //==================================================================
    /// @name public member variables
    //@{
    /// @brief if true, Grid2dPillar::add_structure will apply this dike structure
    bool tf_exec = false;

    /// @brief name identifier for this dike structure
    std::string name = "none";

    /// @brief delta density anomaly added to background density (kg/m³)
    double delta_density = 0.0;

    /// @brief axis-aligned bounding box defining dike extent in horizontal plane (m)
    /// @note AABB is defined before rotation is applied
    AABB2d aabb2d = AABB2d();

    /// @brief rotation angle of the dike around z-axis (counterclockwise from x-axis)
    Angle angle_rot = Angle(0.0, Angle::Unit::Degree);
    //@} ------------------------------------------------------------------

    //==================================================================
    /// @name constructor_destructor
    //@{
    /// @brief default constructor
    VerticalDikeParameters() = default;

    /// @brief copy constructor
    /// @param org source object to copy from
    VerticalDikeParameters(const VerticalDikeParameters& org) = default;

    /// @brief constructor from nlohmann::json configuration
    /// @param js JSON object containing configuration
    /// @param section_name section key within JSON
    VerticalDikeParameters(const nlohmann::json& js, const std::string &section_name)
      : VerticalDikeParameters() { assign_parameters(js, section_name); };

    /// @brief destructor
    ~VerticalDikeParameters() = default;
    //@} ------------------------------------------------------------------

    //==================================================================
    /// @name operators
    //@{
    /// @brief assignment operator
    /// @param other source object to assign from
    /// @return reference to this object
    VerticalDikeParameters& operator=(const VerticalDikeParameters& other) = default;
    //@} ------------------------------------------------------------------

    /// @brief get name of this instance
    /// @return name identifier string
    std::string get_name() const { return name; };

    /// @brief load and assign parameters from JSON configuration
    /// @param js JSON object containing configuration
    /// @param section_name section key within JSON (used for logging)
    /// @note Missing keys will trigger warnings and use defaults
    void assign_parameters(const nlohmann::json& js, const std::string &section_name);
};

//################################################################################
//################################################################################
/// @class Grid2dPillar::VerticalCheckerBoardParameters
/// @brief Parameters for vertical checkerboard structure used in Grid2dPillar::add_density_structure
/// @details
/// Configures a vertical checkerboard pattern (alternating density cells extending infinitely in z).
/// Each cell has density: `background_density + delta_density_offset ± delta_density`.
///
/// @note Thread-safe: No. Modification is not thread-safe.
/// @note Coordinate system: Right-handed, z-up
/// @note Units: meters (m) for lengths, kg/m³ for density
//################################################################################
//################################################################################
class Grid2dPillar::VerticalCheckerBoardParameters {
  private:
  public:
    //==================================================================
    /// @name public member variables
    //@{
    /// @brief if true, Grid2dPillar::add_structure will apply this checkerboard structure
    bool tf_exec = false;

    /// @brief name identifier for this checkerboard structure
    std::string name = "none";

    /// @brief constant offset added to all checkerboard cells (kg/m³)
    /// @note Final density = background + delta_density_offset ± delta_density
    double delta_density_offset = 0.0;

    /// @brief amplitude of alternating density anomaly (kg/m³)
    /// @note Cells alternate between +delta_density and -delta_density
    double delta_density = 0.0;

    /// @brief axis-aligned bounding box defining checkerboard extent in horizontal plane (m)
    AABB2d aabb2d = AABB2d();

    /// @brief center position (x, y) of the checkerboard pattern (m)
    Eigen::Vector2d v2_pos_cnt = Eigen::Vector2d(0.0, 0.0);

    /// @brief cell dimensions (x_length, y_length) for each checkerboard square (m)
    Eigen::Vector2d v2_length = Eigen::Vector2d(0.0, 0.0);
    //@} ------------------------------------------------------------------

    //==================================================================
    /// @name constructor_destructor
    //@{
    /// @brief default constructor
    VerticalCheckerBoardParameters() = default;

    /// @brief copy constructor
    /// @param org source object to copy from
    VerticalCheckerBoardParameters(const VerticalCheckerBoardParameters& org) = default;

    /// @brief constructor from nlohmann::json configuration
    /// @param js JSON object containing configuration
    /// @param section_name section key within JSON
    VerticalCheckerBoardParameters(const nlohmann::json& js, const std::string &section_name)
      : VerticalCheckerBoardParameters() { assign_parameters(js, section_name); };

    /// @brief destructor
    ~VerticalCheckerBoardParameters() = default;
    //@} ------------------------------------------------------------------

    //==================================================================
    /// @name operators
    //@{
    /// @brief assignment operator
    /// @param other source object to assign from
    /// @return reference to this object
    VerticalCheckerBoardParameters& operator=(const VerticalCheckerBoardParameters& other) = default;
    //@} ------------------------------------------------------------------

    /// @brief get name of this instance
    /// @return name identifier string
    std::string get_name() const { return name; };

    /// @brief load and assign parameters from JSON configuration
    /// @param j JSON object containing configuration
    /// @param section_name section key within JSON (used for logging)
    /// @note Missing keys will trigger warnings and use defaults
    void assign_parameters(const nlohmann::json& j, const std::string &section_name);
};

//###################################################################
/// @class Grid2dPillar::Parameters
/// @brief Main parameter configuration for Grid2dPillar construction
/// @details
/// Configures the overall Grid2dPillar object including:
/// - Initial uniform density (kg/m³)
/// - DEM (Digital Elevation Model) file path for topography
/// - Grid alignment (shift flags for x, y)
/// - Vertical structures (cylinders, dikes, checkerboards)
///
/// ### Typical usage
/// @code
/// nlohmann::json config = load_json("config.json");
/// Grid2dPillar::Parameters params(config, "GRID2D_PILLAR_PARAMETERS");
/// Grid2dPillar grid(params);
/// @endcode
///
/// ### JSON configuration
/// The default JSON section key is `GRID2D_PILLAR_PARAMETERS`, but you can
/// define multiple sections with different keys to create multiple Grid2dPillar instances.
///
/// @note Thread-safe: No. Parameter modification is not thread-safe.
/// @note Units: kg/m³ for density, m for spatial coordinates
/// @ingroup parameterClasses
//###################################################################
class Grid2dPillar::Parameters {
  private:
  public:
    //===========================================================================
    /// @name public member variables
    ///@{
    /// @brief name identifier for this Grid2dPillar parameter set
    std::string name = "Grid2dPillarParameters";

    /// @brief initial uniform density in SI units (kg/m³)
    /// @note Must be positive. Typically 1000-3000 kg/m³ for geological materials.
    double initial_uniform_density = 1000.0;

    /// @brief file path to DEM (Digital Elevation Model) data
    /// @note DEM defines the surface topography. Format must be compatible with Grid2dPillar::read_dem.
    std::filesystem::path path_dem = "none";

    /// @brief tolerance ratio for Grid2dPillar::convert_from_vec_xyz coordinate matching
    /// @note Relative tolerance for floating-point coordinate comparisons. Must be positive.
    double tolerance_ratio = 1.0e-3;

    /// @brief minimum z elevation threshold (m a.s.l.)
    /// @note Grid cells below this elevation may be excluded or handled specially
    double zmin = 0.0;

    /// @brief if true, shift x-direction bin boundaries by half-bin width
    /// @note Affects grid cell alignment in x-direction
    bool tf_shift_x = false;

    /// @brief if true, shift y-direction bin boundaries by half-bin width
    /// @note Affects grid cell alignment in y-direction
    bool tf_shift_y = false;

    //================================
    /// @brief JSON section name for vertical cylinder structure parameters array
    static constexpr const char*
      section_name_cylinder = "vertical_cylinder_params";

    /// @brief JSON section name for vertical dike structure parameters array
    static constexpr const char*
      section_name_dike = "vertical_dike_params";

    /// @brief JSON section name for vertical checkerboard structure parameters array
    static constexpr const char*
      section_name_checkerboard = "vertical_checkerboard_params";
    //================================
    /// @brief collection of vertical cylinder density structures to apply
    /// @note Each element defines one cylinder to be added to the density model
    std::vector<VerticalCylinderParameters> vec_vertical_cylinder_parameters;

    /// @brief collection of vertical dike density structures to apply
    /// @note Each element defines one dike to be added to the density model
    std::vector<VerticalDikeParameters> vec_vertical_dike_parameters;

    /// @brief collection of vertical checkerboard density structures to apply
    /// @note Each element defines one checkerboard pattern to be added to the density model
    std::vector<VerticalCheckerBoardParameters> vec_vertical_checkerboard_parameters;
    ///@} ------------------------------------------------------------------

    //==================================================================
    /// @name constructor_destructor
    ///@{

    /// @brief default constructor
    /// @note All members initialized to default values
    Parameters() = default;

    /// @brief copy constructor
    /// @param org source object to copy from
    Parameters(const Parameters& org) = default;

    /// @brief constructor from nlohmann::json configuration
    /// @param js JSON object containing configuration
    /// @param section_name JSON section key (default: "GRID2D_PILLAR_PARAMETERS")
    /// @note Calls assign_parameters internally
    Parameters(const nlohmann::json& js, const std::string &section_name)
      : Parameters() { assign_parameters(js,section_name); };

    /// @brief destructor
    ~Parameters() = default;
    ///@} ------------------------------------------------------------------

    /// @brief get name of this parameter set
    /// @return name identifier string
    std::string get_name() const { return name; };

    /// @brief load and assign all parameters from JSON configuration
    /// @param js JSON object containing configuration
    /// @param section_name JSON section key (default: "GRID2D_PILLAR_PARAMETERS")
    /// @throws std::runtime_error if validation fails (e.g., negative density, empty path_dem)
    /// @note Reads main parameters and all nested structure parameter arrays
    /// @note Missing optional keys will use default values; required keys trigger errors
    /// @note Thread-safe: No
    void assign_parameters(const nlohmann::json& js
      , const std::string &section_name = "GRID2D_PILLAR_PARAMETERS");

    /// @brief assignment operator
    /// @param other source object to assign from
    /// @return reference to this object
    Parameters& operator=(const Parameters& other) = default;
};
