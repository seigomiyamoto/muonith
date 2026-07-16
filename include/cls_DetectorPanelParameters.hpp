/// @file cls_DetectorPanelParameters.hpp
/// @brief Parameter configuration classes for detector panels
/// @details
/// This file defines two nested classes within DetectorPanel:
/// - DetectorPanel::Parameters: Configuration for a single detector panel, including geometry,
///   orientation (yaw/pitch/roll), angular binning, and efficiency table paths.
/// - DetectorPanel::ParameterLists: Configuration for multiple detector panels, typically loaded
///   from a collection of JSON parameter files.
///
/// Typical workflow:
/// 1. Load JSON configuration file(s)
/// 2. Construct Parameters or ParameterLists from JSON sections
/// 3. Use these parameters to initialize DetectorPanel objects
///
/// @note Units: angles in degrees (converted internally), positions/lengths in meters.
/// @note Coordinate system: Right-handed, z-up (yaw=0 is North, 90 is East, clockwise from above).
/// @note Thread-safety: No. Each instance should be constructed/modified by a single thread.
#pragma once
#include "cls_DetectorPanel.hpp"
#include "cls_EfficiencyModel.hpp"
#include <nlohmann/json.hpp>
#include "ns_angle_util.hpp"

//###########################################################################
/// @class DetectorPanel::Parameters
/// @brief Parameter set for configuring a single DetectorPanel instance
/// @details
/// This class holds all configuration parameters needed to construct and initialize a DetectorPanel.
/// Parameters include:
/// - Geometric properties (detector size, position, orientation via yaw/pitch/roll)
/// - Angular binning configuration (nbinx, nbiny, angle ranges, angle unit)
/// - Physical properties (exposure time, number of detector units)
/// - File paths (bin list, efficiency table)
///
/// Invariants:
/// - nbinx, nbiny > 0
/// - txmin < txmax, tymin < tymax
/// - length_hori, length_vert > 0; length_dept >= 0
/// - n_unit > 0
///
/// Typical usage:
/// @code
/// nlohmann::json js = /* load JSON */;
/// DetectorPanel::Parameters params(js, "DETECTOR_PARAMETERS");
/// DetectorPanel panel;
/// panel.setup_from_parameters(params);
/// @endcode
///
/// @note Angles are stored internally as Angle objects but read from JSON as degrees.
/// @note Rotation order: yaw -> pitch -> roll (applied according to rotation_type).
//###########################################################################
class DetectorPanel::Parameters {
  private:
  public:
    //=========================================================================
    /// @name public member variables
    ///@{
    /// @brief Name identifier for this detector panel instance
    std::string name = "none";

    /// @brief If true, read angular binning configuration from external file
    bool tf_read_bin_list = false;

    /// @brief Path to bin list file (used when tf_read_bin_list is true)
    std::filesystem::path filepath_bin_list = "none";

    /// @brief Number of angular bins in x (horizontal) and y (vertical) directions
    /// @note Must be > 0
    int nbinx = 0, nbiny = 0; 

    /// @brief Angular unit for binning (Tangent, Degree, or Radian)
    /// @note This determines interpretation of txmin/txmax/tymin/tymax ranges
    DetectorElement::AngleUnit angle_unit = DetectorElement::AngleUnit::Tangent;

    /// @brief Angular bin ranges: [txmin, txmax] x [tymin, tymax]
    /// @note Units depend on angle_unit (tangent, degrees, or radians)
    /// @note Must satisfy: txmin < txmax, tymin < tymax
    /// @note For Degree: -90 <= tymin, tymax <= 90; for Radian: -π/2 <= tymin, tymax <= π/2
    double txmin = 0.0, txmax = 0.0, tymin = 0.0, tymax = 0.0;

    /// @brief Detector physical dimensions: (horizontal, vertical, depth)
    /// @note Units: meters
    /// @note horizontal and vertical must be > 0; depth must be >= 0
    Eigen::Vector3d v3_det_length = Eigen::Vector3d(0.0, 0.0, 0.0);

    /// @brief Number of detector units (e.g., scintillator modules or pixel count)
    /// @note Must be > 0
    double n_unit = 0.0;

    /// @brief Detector orientation angles: yaw, pitch, roll
    /// @details
    /// Definition: \image html def_yaw_pitch_roll.dio.png
    /// - Yaw: rotation around vertical (z) axis
    ///   - 0° = North, 90° = East, 180° = South, -90° = West (clockwise from above)
    /// - Pitch: rotation around horizontal transverse axis
    /// - Roll: rotation around horizontal longitudinal axis
    /// @note Rotation order: yaw -> pitch -> roll (composition depends on rotation_type)
    /// @note Units: degrees (stored internally as Angle objects)
    /// @note Coordinate system: right-handed, z-up
    Angle yaw = Angle(0.0, Angle::Unit::Degree),
          pitch = Angle(0.0, Angle::Unit::Degree),
          roll = Angle(0.0, Angle::Unit::Degree);

    /// @brief Rotation composition type (LOCAL or GLOBAL frame)
    /// @details
    /// - LOCAL: Intrinsic rotations - each rotation is applied in the current (rotated) frame.
    ///   The reference axes rotate with the object.
    /// - GLOBAL: Extrinsic rotations - all rotations are applied in the fixed (original) frame.
    ///   The reference axes remain fixed.
    /// @note Default is LOCAL
    angle_util::Rotation3dType rotation_type = angle_util::Rotation3dType::LOCAL;

    /// @brief Detector center position in global coordinates (x, y, z)
    /// @note Units: meters
    /// @note Coordinate system: right-handed, z-up
    Eigen::Vector3d v3_position = Eigen::Vector3d(0.0, 0.0, 0.0);

    /// @brief Detector exposure time
    /// @note Units: days
    double days = 0.0;

    /// @brief Number of elements to pre-allocate in vec_tf_in_PL (performance hint)
    /// @note Must be >= 0. Default is 10.
    int n_reserve_vec_tf_in_PL = 10;

    /// @brief Path to detector efficiency table file
    /// @note If set to "none", no efficiency correction is applied
    /// @note Legacy path, kept for migration; path_eff_model takes precedence when set
    std::filesystem::path path_eff_table = "none";

    /// @brief Path to the efficiency-model JSON5 file (eff_table_lab/configs format)
    /// @note If set to "none" (default), the legacy path_eff_table is used instead
    std::filesystem::path path_eff_model = "none";

    /// @brief Efficiency model loaded from path_eff_model
    /// @note Empty (get_tf_loaded() == false) when path_eff_model is "none"
    EfficiencyModel eff_model = EfficiencyModel();

    ///@} ---------------------------------------------------------------------

    //=========================================================================
    /// @name constructor_destructor
    ///@{
    
    /// @brief default constructor
    Parameters() = default;
    
    /// @brief copy constructor
    Parameters(const Parameters& org) = default;

    /// @brief move constructor
    Parameters(Parameters&& other) noexcept = default;

    /// @brief destructor
    ~Parameters() = default;

    /// @brief constructor from json
    /// @param js json object
    /// @param section_name section name in json
    Parameters(const nlohmann::json& js,
      const std::string &section_name = "DETECTOR_PARAMETERS")
      : Parameters() { assign_parameters(js, section_name); }
    ///@} ---------------------------------------------------------------------

    /// @brief Get the name identifier of this parameter set
    /// @return Name string
    std::string get_name() const { return name; };

    /// @brief Load and assign parameters from JSON configuration
    /// @param[in] js JSON object containing detector configuration
    /// @param[in] section_name Section name within JSON to read from (default: "DETECTOR_PARAMETERS")
    /// @throws std::runtime_error If required parameters are missing or invalid
    /// @throws std::runtime_error If angle_unit string is not recognized
    /// @throws std::runtime_error If rotation_type string is not "LOCAL" or "GLOBAL"
    /// @note Validates all parameters (bin counts > 0, ranges, dimensions) after loading
    void assign_parameters(const nlohmann::json& js,
      const std::string& section_name = "DETECTOR_PARAMETERS" );

    /// @brief assignment operator
    Parameters& operator=(const Parameters& other) = default;
};

//##################################################################################
/// @class DetectorPanel::ParameterLists
/// @brief Container for multiple detector panel parameter file paths
/// @details
/// This class manages configuration for multiple detector panels, typically used when
/// setting up a detector array. It stores:
/// - A list of parameter file paths, each describing one detector panel
/// - A global efficiency application flag
///
/// Typical usage:
/// @code
/// nlohmann::json js = /* load JSON */;
/// DetectorPanel::ParameterLists param_lists(js, "DETECTOR_ARRAY");
/// for(const auto& path : param_lists.vec_parameter_file_path) {
///   // Load each detector parameter file and construct DetectorPanel
/// }
/// @endcode
///
/// @note Thread-safety: No. Construct/modify in a single thread.
//##################################################################################
class DetectorPanel::ParameterLists {
  public:
    /// @brief Name identifier for this parameter list
    std::string name = "DetectorPanel::ParameterLists";

    /// @brief List of parameter file paths, one per detector panel
    /// @note Each path should point to a valid JSON file containing detector parameters
    std::vector<std::filesystem::path> vec_parameter_file_path;

    /// @brief If true, apply detector efficiency corrections to signals
    bool tf_apply_eff = false;

    /// @brief If true, write the per-element txty ASCII dumps under tmp/
    /// @details Controls out_txtyPL(), out_txtyDens() and out_txtySignal(),
    ///          which write tmp/(name)_txty_(PL|dens|signal)_det(NN).tmp.
    /// @note Default true reproduces the behavior of releases before this flag existed.
    /// @note Setting it to false only stops the ASCII dumps under tmp/. With the current
    ///       swp001 settings the per-element figures are plotted from the saved binary,
    ///       so no figure is lost.
    bool tf_out_txty_ascii = true;

    /// @brief If true, write the g2bg ASCII dumps under tmp/
    /// @details Controls out_g2bg_all(), which writes tmp/g2bg_(name)_det(NN).tmp.
    /// @note Default true reproduces the behavior of releases before this flag existed.
    /// @note Setting it to false only stops the ASCII dumps under tmp/. With the current
    ///       swp001 settings the efficiency figures are plotted from the saved binary,
    ///       so no figure is lost.
    bool tf_out_g2bg_ascii = true;

    /// @brief If true, save the prior DetectorPanelArray as det/arrdet_g3vox_prior(suffix).bin
    /// @details Controls the save() call at the end of
    ///          exemdl::build_prior::rebuild_prior_from_shell_PL. Each file is about
    ///          160 MB, and tf_prior_error = true writes three of them (lower, center, upper).
    /// @note Default false. Enable it only when the prior array must be plotted from
    ///       the binary, e.g. with tf_out_txty_ascii = false.
    bool tf_save_arrdet_prior = false;

    //=========================================================================
    /// @name angle_bin_override
    /// @brief Members for overriding angular binning parameters
    /// @details When tf_override_angle_bin is true, the values below override
    ///          the angular binning settings in each detector's parameter file.
    ///@{

    /// @brief If true, override angular binning for all detectors
    bool tf_override_angle_bin = false;

    /// @brief Override values for angular binning (used when tf_override_angle_bin is true)
    int nbinx_override = 0;
    double txmin_override = 0.0;
    double txmax_override = 0.0;
    int nbiny_override = 0;
    double tymin_override = 0.0;
    double tymax_override = 0.0;

    /// @brief Set angular binning override values
    /// @param[in] nbinx Number of bins in x (horizontal) direction
    /// @param[in] txmin Minimum x angle
    /// @param[in] txmax Maximum x angle
    /// @param[in] nbiny Number of bins in y (vertical) direction
    /// @param[in] tymin Minimum y angle
    /// @param[in] tymax Maximum y angle
    /// @note This also sets tf_override_angle_bin to true
    void set_angle_bin_override(int nbinx, double txmin, double txmax,
                                int nbiny, double tymin, double tymax);
    ///@}

    /// @brief Get the number of detector panels
    /// @return Number of parameter files (i.e., number of detectors)
    int get_n_det() const { return vec_parameter_file_path.size(); }

    /// @brief default constructor
    ParameterLists() = default;
  
    /// @brief copy constructor
    ParameterLists(const ParameterLists &org) = default;
  
    /// @brief build from json
    /// @param js json object
    /// @param section_name section name in json
    ParameterLists(const nlohmann::json& js, const std::string &section_name)
      : ParameterLists() { assign_parameters(js, section_name); }
  
    /// @brief assignment operator
    ParameterLists& operator=(const ParameterLists& other) = default;
  
    /// @brief Load and assign parameter file paths from JSON configuration
    /// @param[in] js JSON object containing detector array configuration
    /// @param[in] section_name Section name within JSON to read from
    /// @throws std::runtime_error If any file in vec_parameter_file_path does not exist
    /// @note Reads "det_files" array from the specified section
    /// @note Validates file existence using myapp::filecheck()
    void assign_parameters(const nlohmann::json& js, const std::string &section_name);
};
