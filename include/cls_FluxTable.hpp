/// @file cls_FluxTable.hpp
/// @brief Flux table class for muon flux data management
/// @details
/// This file defines the FluxTable class, which manages muon flux data in tabular format
/// for ray tracing calculations. The class loads and stores two types of 2D flux tables:
/// - Integrated penetrating muon flux (log10 scale) as a function of cos(θz) and Range
/// - Differential flux dF/dR as a function of cos(θz) and Range
///
/// The tables are loaded from file paths specified in JSON configuration and stored as
/// Grid2dXYZ objects. All flux values use bilinear interpolation for arbitrary query points.
///
/// @note Units: Range in kg/m², flux in (m² sec sr)⁻¹
/// @note Thread-safety: Read-only after construction. Safe for concurrent reads.
#pragma once
#include <optional>
#include <string>
#include <filesystem>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include "ns_mymacro.hpp"
#include "ns_mylogger.hpp"
#include "cls_Grid2dXYZ.hpp"
#include "cls_Grid1dXZ.hpp"
#include "cls_DetectorElement.hpp"
#include "cls_FluxTable.hpp"

namespace fs = std::filesystem;

//##############################################
//##############################################
/// @class FluxTable
/// @brief Read-only container for muon flux tables loaded from files
/// @ingroup flux_tools
///
/// @details
/// FluxTable manages two types of muon flux data tables:
/// 1. Integrated penetrating muon flux (log10 scale): g2_log_penef_R_costhz
/// 2. Differential flux dF/dR: g2_dFdR_R_costhz
///
/// Both tables are 2D grids with:
/// - x-axis: cos(θz) where θz is the zenith angle (vertical direction)
/// - y-axis: Range R in kg/m² (path length through material)
/// - z-axis: flux value (log10 scale for peneflux, linear for dFdR)
///
/// **Workflow:**
/// 1. Construct with JSON config containing file paths and bin type flags
/// 2. Tables are automatically loaded via load_tables()
/// 3. Query flux values using get_peneflux() or calc_dFdR()
///
/// **Usage example:**
/// @code
/// json config = {{"FluxTableSection", {
///   {"pathin_log_peneflux", "path/to/peneflux.dat"},
///   {"pathin_dFdR_R_costhz", "path/to/dFdR.dat"}
/// }}};
/// FluxTable ft(config, "FluxTableSection");
/// double flux = ft.get_peneflux(0.8, 1000.0); // cos(θz)=0.8, DL=1000 kg/m²
/// @endcode
///
/// @note Units: Range in kg/m², flux in (m² sec sr)⁻¹, angles via cos(θz)
/// @note Thread-safety: Read-only after construction. Safe for concurrent reads.
/// @note The 'name' member is not compared in operator!=
//##############################################
//##############################################
class FluxTable {
private:
  /// @brief name of this instance
  std::string name;

  /// @brief Integrated penetrating flux table (log10 scale)
  /// @note x=linear costhz, y=linear Range
  std::optional<Grid2dXYZ> g2_log_penef_R_costhz;

  /// @brief dFdR table (Range × costhz)
  /// @note x=linear costhz, y=linear Range
  std::optional<Grid2dXYZ> g2_dFdR_R_costhz;

  /// @brief Path to penetrating muon flux data table
  fs::path pathin_log_peneflux{"none"};

  /// @brief Bin types for x/y axes of pene muflux table
  bool tf_xcnt_peneflux{true};
  bool tf_ycnt_peneflux{false};

  /// @brief Path to g2_dFdR_R_costhz table
  fs::path pathin_dFdR_R_costhz{"none"};

  /// @brief Bin types for x/y axes of g2_dFdR_R_costhz table
  bool tf_xcnt_dFdR{true};
  bool tf_ycnt_dFdR{false};

  /// @brief Enable the dF/dR divergence (numerical-noise) check after loading
  bool tf_check_dFdR_divergence{true};
  /// @brief Sign-flip ratio threshold per costhz slice for the divergence check
  /// @note Default 0.05; see @ref check_dFdR_divergence() for the threshold calibration.
  double dFdR_divergence_threshold{0.05};
  /// @brief If true, abort on a divergent dF/dR table; if false, only warn
  bool tf_dFdR_divergence_fatal{true};


public:
  //==================================================================
  /// @name constructor_destructor
  ///@{

  /// @brief default constructor
  FluxTable() = default;

  /// @brief Copy constructor
  FluxTable(const FluxTable &other) = default;

  /// @brief Move constructor
  FluxTable(FluxTable &&other) noexcept = default;

  /// @brief destructor
  ~FluxTable() = default;

  /// @brief Constructor that loads flux tables from files according to JSON parameters
  /// @param[in] js JSON object containing configuration
  /// @param[in] section_name Section name in JSON to read parameters from
  /// @throws std::runtime_error If section_name not found in JSON
  /// @throws std::runtime_error If file paths are invalid or files cannot be loaded
  /// @note Automatically calls assign_parameters() and load_tables()
  explicit FluxTable(const json &js, const std::string &section_name){
    assign_parameters(js, section_name);
    load_tables();
  };

  /// @brief Assign parameters from JSON configuration
  /// @param[in] js JSON object containing configuration
  /// @param[in] section_name Section name in JSON to read parameters from
  /// @throws std::runtime_error If section_name not found in JSON
  /// @details Reads the following optional keys from JSON:
  /// - pathin_log_peneflux: file path to penetrating flux table
  /// - tf_xcnt_peneflux: x-axis bin type flag for peneflux table
  /// - tf_ycnt_peneflux: y-axis bin type flag for peneflux table
  /// - pathin_dFdR_R_costhz: file path to dFdR table
  /// - tf_xcnt_dFdR: x-axis bin type flag for dFdR table
  /// - tf_ycnt_dFdR: y-axis bin type flag for dFdR table
  void assign_parameters(
    const nlohmann::json& js, const std::string& section_name);

  /// @brief Load Grid2dXYZ tables from file paths specified in assign_parameters()
  /// @throws std::runtime_error If file paths are invalid or files cannot be loaded
  /// @note Skips loading if path is set to "none" (logs warning instead)
  /// @note Tables are stored in std::optional members and only allocated if loaded
  void load_tables();

  ///@}---------------------------------------------------------------
    //==================================================================
  /// @name operators
  ///@{

  /// @brief Assignment operator
  FluxTable& operator=(const FluxTable& other) = default;

  /// @brief Inequality operator for comparing two FluxTable instances
  /// @param[in] other FluxTable instance to compare with
  /// @return true if any member (except name) differs, false otherwise
  /// @note The 'name' member is intentionally excluded from comparison
  /// @note In non-NODEBUG builds, logs a warning when differences are found
  bool operator!=(const FluxTable& other) const;

  /// @brief Equality operator defined via the inequality operator
  /// @param[in] other FluxTable instance to compare with
  /// @return true if all members (except name) are equal, false otherwise
  /// @note Implemented as negation of operator!=
  bool operator==(const FluxTable& other) const{
    return !(*this != other);
  }

  ///@}-------------------------------------------------------------------
  //==================================================================
  /// @name setter_functions
  ///@{

    /// @brief Set name of this FluxTable instance
    /// @param[in] name_in Name to assign
    void set_name( const std::string &name_in ){ name = name_in; };

    /// @brief Set integrated penetrating muon flux table
    /// @param[in] g2_log_peneflux_in Grid2dXYZ containing log10(flux) data
    /// @note x-axis: cos(θz), y-axis: Range (kg/m²), z-axis: log10(flux) in (m² sec sr)⁻¹
    void set_g2_log_penef( const Grid2dXYZ &g2_log_peneflux_in ){
      g2_log_penef_R_costhz = g2_log_peneflux_in;
    };

    /// @brief Set differential flux dF/dR table
    /// @param[in] g2_dFdR_R_costhz_in Grid2dXYZ containing dF/dR data
    /// @note x-axis: cos(θz), y-axis: Range (kg/m²), z-axis: dF/dR
    void set_g2_dFdR( const Grid2dXYZ &g2_dFdR_R_costhz_in ){
      g2_dFdR_R_costhz = g2_dFdR_R_costhz_in;
    };

  ///@}-------------------------------------------------------------------

  //==================================================================
  /// @name getter_functions
  ///@{

  /// @brief Get name of this FluxTable instance
  /// @return Name string
  std::string get_name() const { return name; };

  /// @brief Get reference to integrated penetrating flux table
  /// @return Const reference to Grid2dXYZ containing log10(flux) data
  /// @throws std::runtime_error If g2_log_penef_R_costhz is not loaded
  /// @note x-axis: cos(θz) (linear), y-axis: Range (kg/m², linear), z-axis: log10(flux)
  /// @note Call load_tables() before using this function
  const Grid2dXYZ& get_g2_log_peneflux() const;

  /// @brief Get reference to differential flux dF/dR table
  /// @return Const reference to Grid2dXYZ containing dF/dR data
  /// @throws std::runtime_error If g2_dFdR_R_costhz is not loaded
  /// @note x-axis: cos(θz) (linear), y-axis: Range (kg/m², linear), z-axis: dF/dR
  /// @note Call load_tables() before using this function
  const Grid2dXYZ& get_g2_dFdR_R_costhz() const;

  /// @brief Get penetrating muon flux from integrated table using bilinear interpolation
  /// @param[in] costhz Value of cos(θz) where θz is zenith angle
  /// @param[in] DL Path length / Range in kg/m²
  /// @return Interpolated log10(flux) value
  /// @throws std::runtime_error If g2_log_penef_R_costhz is not loaded
  /// @note Call load_tables() before using this function
  double get_peneflux(const double costhz, const double DL) const;

  /// @brief Calculate differential flux dF/dR from DetectorElement and path length
  /// @param[in] ele Detector element containing ray direction information
  /// @param[in] DL_in Path length / Range in kg/m²
  /// @return Interpolated dF/dR value using bilinear interpolation
  /// @throws std::runtime_error If g2_dFdR_R_costhz is not loaded
  /// @note Extracts cos(θz) from ele.get_ray3d().vz()
  /// @note DL_in is automatically clamped to valid Range bounds [Rmin, Rmax]
  /// @note Call load_tables() before using this function
  double calc_dFdR(const DetectorElement& ele, const double DL_in) const;
  ///@}-------------------------------------------------------------------

  //==================================================================
  /// @name checker_functions
  ///@{

  /// @brief Check whether integrated penetrating flux table is loaded
  /// @return true if g2_log_penef_R_costhz is loaded, false otherwise
  bool is_loaded_log_peneflux() const {
    return g2_log_penef_R_costhz.has_value();
  };

  /// @brief Check whether differential flux dF/dR table is loaded
  /// @return true if g2_dFdR_R_costhz is loaded, false otherwise
  bool is_loaded_dFdR() const {
    return g2_dFdR_R_costhz.has_value();
  };
  ///@}-------------------------------------------------------------------

};

//======================================================================
/// @brief Inspect a dF/dR table for divergence (numerical noise) and warn or abort.
/// @param[in] g2_dFdR_R_costhz Finished dF/dR table (x = costhz, y = range R, z = dF/dR).
/// @param[in] threshold Maximum tolerated sign-flip ratio per costhz slice (e.g. 0.05).
/// @param[in] fatal If true, throw on a divergent table; if false, only LOG_WARN and continue.
/// @param[in] context Caller label used in log and exception messages.
/// @throws std::runtime_error If the table has non-finite values, or if it looks
///         divergent and fatal == true.
/// @note A healthy dF/dR table is finite and strictly negative everywhere
///       (F(R) decreases monotonically with range R), so positive or
///       sign-flipping values are the numerical-noise signature seen in
///       issue 90gy73. Non-finite values are always fatal; the divergence
///       signal (positive-value ratio and per-costhz sign-flip ratio) aborts
///       only when fatal == true.
/// @note Threshold calibration. The sign-flip ratio is flips / (points - 1)
///       per costhz slice, measured on ~4995 range points per slice:
///         - Healthy (daemon, fixed): 0 flips   -> ratio ~ 0
///         - Healthy (Honda):         2 flips   -> ratio ~ 4e-4
///         - Divergent (Honda, misconfig): 902 flips -> ratio ~ 0.18
///         - Divergent (daemon, broken):   905 flips -> ratio ~ 0.18
///       The healthy and divergent populations are ~2 orders of magnitude
///       apart, so any value in the 0.05-0.10 band separates them with margin;
///       the default 0.05 is the conservative lower end. This default is
///       provisional: confirm zero false positives on several healthy tables
///       before treating it as final. The all-negative invariant (positive-value
///       ratio: healthy ~0 vs divergent ~0.5) is the primary separator, so this
///       ratio is a complementary check.
/// @note Complexity: O(nbinx * nbiny), single pass. Thread-safety: Yes (read-only).
void check_dFdR_divergence(
  const Grid2dXYZ& g2_dFdR_R_costhz
, const double threshold
, const bool fatal
, const std::string& context );
