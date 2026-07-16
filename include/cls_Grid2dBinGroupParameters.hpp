/// @file cls_Grid2dBinGroupParameters.hpp
/// @brief Parameter configuration for 2D bin grouping
/// @details
/// This file defines the Grid2dBinGroup::Parameters class, which configures
/// bin grouping and subdivision behavior for 2D grid-based detector analysis.
///
/// **Workflow:**
/// 1. Create Parameters instance from JSON runcard
/// 2. Configure grouping mode (initial/automatic/manual)
/// 3. Set thresholds (PL, DL, signal/noise ratio)
/// 4. Specify subdivision limits (minimum size, maximum iterations)
/// 5. Optionally load manual grouping from external files
///
/// **Units:**
/// - PathLength (PL) and DensityLength (DL): detector-specific units
/// - signal_noise_group_trig: dimensionless ratio
/// - Spatial dimensions (ixlen_min, iylen_min): bin counts
///
/// **Thread-safety:**
/// This class is not thread-safe. Parameters should be configured before
/// parallel processing begins. Read-only access during parallel sections is safe.
///
/// @ingroup parameterClasses
#pragma once

#include <map>
#include <set>
#include <fstream>
#include <iostream>
#include <sstream> // istringstream
#include <string>
#include <cstdio>
#include <cmath>
#include <functional>  //for sorting
#include <algorithm>//for sorting
#include <vector>
#include <filesystem> // for std::filesystem::path
#include <nlohmann/json.hpp> // for nlohmann::json

//###################################################################
//###################################################################
/// @class Grid2dBinGroup::Parameters
/// @brief Configuration parameters for 2D bin grouping and subdivision
///
/// **Responsibilities:**
/// - Configure initial and automatic bin grouping modes
/// - Define quality thresholds (PathLength, DensityLength, signal/noise)
/// - Specify subdivision constraints (minimum size, iteration limits)
/// - Manage manual grouping file paths for multi-detector setups
///
/// **Invariants:**
/// - PL_thres >= 0.0 and DL_thres >= 0.0
/// - nx_div_init >= 1 and ny_div_init >= 1
/// - ixlen_min >= 1 and iylen_min >= 1
/// - nloop_limit >= 1
/// - igroup_start >= 0 and n_detector_grouping_manual >= 0
///
/// **Usage example:**
/// @code
/// nlohmann::json js = /* load from file */;
/// Grid2dBinGroup::Parameters params(js, "BINGROUP_PARAMETERS");
/// // Use params to configure Grid2dBinGroup instance
/// @endcode
///
/// @note The default JSON key is "BINGROUP_PARAMETERS", but custom keys can be
/// specified to create multiple parameter configurations in the same runcard.
///
/// **Thread-safety:** Not thread-safe for modification. Safe for concurrent
/// read-only access after initialization.
///
/// @ingroup parameterClasses
//###################################################################
//###################################################################
class Grid2dBinGroup::Parameters {
  private:
  public:
    //==================================================================
    // @name public member variables
    //@{
    /// @brief name of this instance
    std::string name = "Grid2dBinGroup::Parameters";

    /// @brief initial value of signal and noise
    double signal_init = 0.0;
    double noise_init = 0.0;

    /// @brief initial value of is_avail
    bool is_avail_init = true;

    /// @brief PathLength threshold
    double PL_thres = 0.0;

    /// @brief DensityLength threshold
    double DL_thres = 0.0;

    /// @brief Signal value for bins where PL < PL_thres or DL < DL_thres
    double signal_under_thres = 0.0;

    /// @brief Noise value for bins where PL < PL_thres or DL < DL_thres
    double noise_under_thres = 0.0;

    /// @brief Availability flag for bins where PL < PL_thres or DL < DL_thres
    bool is_avail_under_thres = false;

    /// @brief Whether to run initial bin grouping
    bool tf_run_1st_grouping = false;

    /// @brief Whether to run automatic bin grouping
    bool tf_run_auto_grouping = false;

    /// @brief Starting group ID
    int igroup_start = 0;

    /// @brief Initial number of divisions in X and Y directions
    int nx_div_init = 1;
    int ny_div_init = 1;

    /// @brief Signal-to-noise threshold for bin subdivision
    /// @details Bins with signal/noise ratio >= this value will be subdivided
    double signal_noise_group_trig = 20.0;

    /// @brief Minimum bin size in X and Y directions
    /// @details Bins smaller than these dimensions will not be subdivided
    int ixlen_min = 1;
    int iylen_min = 1;

    /// @brief Preferred split direction when xlen == ylen
    /// @details If true, prefer splitting in Y direction; if false, prefer X direction
    bool tf_prefer_split_x = false;

    /// @brief Maximum number of subdivision iterations
    int nloop_limit = 10000;

    /// @brief Number of detectors for manual grouping
    int n_detector_grouping_manual = 0;

    /// @brief Flags indicating whether to read bin group list from external file for each detector
    std::vector<bool> vec_tf_read_bin_group_list;

    /// @brief File paths for bin group list files
    std::vector<std::filesystem::path> vec_file_path_bin_group_list;
    //@} ------------------------------------------------------------------

    //==================================================================
    /// @name constructor_destructor
    ///@{
    
    /// @brief default constructor
    Parameters() = default;

    /// @brief copy constructor
    Parameters(const Grid2dBinGroup::Parameters& org) = default;

    /// @brief move constructor
    Parameters(Grid2dBinGroup::Parameters&& other) = default;

    /// @brief destructor
    ~Parameters() = default;

    /// @brief Construct from JSON configuration
    /// @param[in] js JSON object containing configuration
    /// @param[in] section_name Section key in JSON (default: "BINGROUP_PARAMETERS")
    /// @throws std::runtime_error If parameter validation fails (via assert in assign_parameters)
    Parameters(const nlohmann::json& js, const std::string &section_name = "BINGROUP_PARAMETERS")
      : Parameters() { assign_parameters(js, section_name); };
    ///@} ------------------------------------------------------------------

    //==================================================================
    /// @name operators
    ///@{
    
    /// @brief Copy assignment operator
    Parameters& operator=(const Grid2dBinGroup::Parameters& other) = default;

    /// @brief Inequality operator
    /// @details Name field is not compared
    bool operator!=(const Grid2dBinGroup::Parameters& other) const;

    /// @brief Equality operator
    /// @details Name field is not compared. Defined in terms of operator!=
    bool operator==(const Grid2dBinGroup::Parameters& other) const {
      return !(*this != other);
    };
    ///@} ------------------------------------------------------------------

    /// @brief Get instance name
    /// @return Copy of the name string
    std::string get_name() const { return name; };

    /// @brief Load and validate parameters from JSON
    /// @param[in] js JSON object containing configuration
    /// @param[in] section_name Section key to read from
    /// @throws std::runtime_error If validation fails (non-negative checks, minimum bounds)
    /// @note Validates: PL_thres >= 0, DL_thres >= 0, nx_div_init >= 1, ny_div_init >= 1,
    ///       ixlen_min >= 1, iylen_min >= 1, nloop_limit >= 1, igroup_start >= 0
    void assign_parameters(const nlohmann::json& js, const std::string &section_name);

    //==================================================================
    /// @name binary_IO_functions
    ///@{

    /// @brief Serialize all member variables to binary stream
    /// @param[in,out] ofs Output file stream
    /// @throws std::runtime_error If stream operation fails
    void save( std::ofstream& ofs ) const;

    /// @brief Deserialize all member variables from binary stream
    /// @param[in,out] ifs Input file stream
    /// @throws std::runtime_error If stream operation fails
    /// @note Must match save() format exactly for proper deserialization
    void load( std::ifstream& ifs );
    ///@}
};
