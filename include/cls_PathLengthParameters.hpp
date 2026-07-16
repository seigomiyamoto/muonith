/// @file cls_PathLengthParameters.hpp
/// @brief Path length calculation parameter configuration
/// @details
/// Defines pathcalc::Parameters for controlling path length computation settings.
///
/// ## Typical workflow
/// 1. Create a JSON runcard with [PATH_LENGTH_PARAMETERS] section
/// 2. Load JSON with nlohmann::json
/// 3. Construct pathcalc::Parameters from JSON
/// 4. Access path length range (PL_min, PL_max, PL_pit) and I/O flags
///
/// ## Parameter groups
/// - **naive_calc_parameters**: PL_min, PL_max, PL_pit for path length range
/// - **g2pil parameters**: tf_load/save_arrdet_g2pil, path_arrdet_g2pil_bin
/// - **g3vox parameters**: tf_load/save_arrdet_g3vox, shell addition, sparse matrix settings
/// - **data_load_save_parameters**: Binary matrix I/O settings for dN/dD observation matrices
///
/// ## Units
/// - Path lengths (PL_min, PL_max, PL_pit): meters
/// - Sparse matrix tolerance: dimensionless ratio
///
/// ## Thread safety
/// - Thread-safe: No (not designed for concurrent modification)
/// - Instances should be constructed once and used as const references
///
/// @see pathcalc::Parameters
#pragma once

#include <string>
#include <filesystem>
#include <nlohmann/json.hpp>

/// @brief Namespace for path length calculation utilities
namespace pathcalc {
  class Parameters;
}

/// @class pathcalc::Parameters
/// @brief Configuration parameters for path length calculations
///
/// @details
/// This class holds all parameters needed for path length computation in the
/// raytracing pipeline. Parameters are typically loaded from a JSON runcard
/// under the [PATH_LENGTH_PARAMETERS] section.
///
/// ## Responsibilities
/// - Define path length range (min, max, pitch) in meters
/// - Control loading/saving of precomputed detector arrays (g2pil, g3vox)
/// - Configure sparse matrix tolerance for path length matrices
/// - Manage binary I/O paths for observation matrices
///
/// ## Usage example
/// @code
/// nlohmann::json js = nlohmann::json::parse(runcard_stream);
/// pathcalc::Parameters params(js, "PATH_LENGTH_PARAMETERS");
/// // Access parameters
/// double range = params.PL_max - params.PL_min;
/// int n_bins = static_cast<int>((params.PL_max - params.PL_min) / params.PL_pit);
/// @endcode
///
/// ## Multiple instances
/// The default JSON key is "PATH_LENGTH_PARAMETERS", but you can specify
/// different keys to create multiple parameter instances from the same runcard.
///
/// @note Name field is not compared in equality operators (for instance tracking)
/// @ingroup parameterClasses
class pathcalc::Parameters {
  public:
    /// @brief Instance name for identification/logging purposes
    std::string name = "none";

    //======================================================================
    /// @name naive_calc_parameters
    /// @brief Basic path length range parameters
    ///@{

    /// @brief Minimum path length [meters]
    double PL_min = 0.0;

    /// @brief Maximum path length [meters]
    double PL_max = 0.0;

    /// @brief Path length pitch (bin width) [meters]
    double PL_pit = 1.0;

    /// @brief Maximum beam length [meters] for ray tracing warning
    /// @details If BL_max > 0, rays whose total beam length (tmax from DEM AABB intersection)
    /// exceeds this value will log a warning (but still be processed normally).
    /// If BL_max <= 0, no beam length check is performed.
    /// This is useful for identifying low-elevation rays that extend beyond the DEM boundary.
    /// @note Default value is param_constants::BL_max_default() (set during JSON loading).
    double BL_max = 0.0;

    ///@}

    //======================================================================
    /// @name g2pil_parameters
    /// @brief Path length I/O parameters for g2pil (2D cubic grid) detectors
    ///@{

    /// @brief If true, load precomputed arrdet_g2pil from binary file
    bool tf_load_arrdet_g2pil = false;

    /// @brief If true, save computed arrdet_g2pil to binary file
    bool tf_save_arrdet_g2pil = false;

    /// @brief Binary file path for arrdet_g2pil data
    std::filesystem::path path_arrdet_g2pil_bin = "arrdet_g2pil.tmp.bin";

    ///@}
    
    //======================================================================
    /// @name g3vox_io_parameters
    /// @brief Path length I/O parameters for g3vox (3D voxel grid) and sparse matrix
    ///@{

    /// @brief If true, load precomputed arrdet_g3vox from binary file
    bool tf_load_arrdet_g3vox = false;

    /// @brief If true, save computed arrdet_g3vox to binary file
    bool tf_save_arrdet_g3vox = false;

    /// @brief If true, also compute path lengths for g2pil_shell along with g3vox
    bool tf_add_shell = true;

    /// @brief Binary file path for arrdet_g3vox data
    std::filesystem::path path_arrdet_g3vox_bin = "arrdet_g3vox.tmp.bin";

    /// @brief Binary file path for vec_spmat_PL (sparse path length matrix vector)
    std::filesystem::path path_vec_spmat_PL_bin = "vec_spmat_PL.tmp.bin";

    ///@}
    
    //======================================================================
    /// @name g3vox_calc_parameters
    /// @brief Calculation control parameters for g3vox path length computation
    ///@{

    /// @brief If true, add computed path length to detector accumulator
    bool tf_add_PLDL = true;

    /// @brief If true, increment nhit_det counter for g3vox hits
    bool tf_incr_nhit_det = true;

    /// @brief If true, increment nhit_ele counter for element hits
    bool tf_incr_nhit_ele = true;

    /// @brief Reference order of magnitude for sparse matrix values
    /// @details Used with epsilon_matPL_sparse to compute tolerance threshold
    float reference_matPL_sparse = 1.0f;

    /// @brief Relative precision for treating sparse matrix entries as zero
    /// @details Values below reference_matPL_sparse * epsilon_matPL_sparse are zeroed
    float epsilon_matPL_sparse = 1.0E-9f;

    /// @brief Compute absolute tolerance for sparse matrix zero-detection
    /// @return reference_matPL_sparse * epsilon_matPL_sparse
    float get_tolerance() const { return reference_matPL_sparse * epsilon_matPL_sparse; }

    ///@}
    
    //======================================================================
    /// @name data_load_save_parameters
    /// @brief Binary I/O parameters for dN/dD observation matrices
    ///@{

    /// @brief If true, load obs_mat_dNdD from binary instead of computing
    /// @details Controls whether calc_dNdD::make_grouped_alldet_mat_dNdD runs
    bool tf_load_bin_obs_mat_dNdD = false;

    /// @brief If true, save computed obs_mat_dNdD to binary file
    bool tf_save_bin_obs_mat_dNdD = false;

    /// @brief Binary file path for obs_mat_dNdD data
    std::filesystem::path path_bin_obs_mat_dNdD = "none";

    ///@}
    
    /// @brief Default constructor with default parameter values
    Parameters() = default;

    /// @brief Copy constructor
    Parameters(const pathcalc::Parameters& org) = default;

    /// @brief Move constructor
    Parameters(pathcalc::Parameters&& org) noexcept = default;

    /// @brief Destructor
    ~Parameters() = default;

    /// @brief Construct from JSON configuration
    /// @param[in] js JSON object containing the runcard
    /// @param[in] section_name Section key in JSON (default: "PATH_LENGTH_PARAMETERS")
    /// @throws std::runtime_error if required keys are missing from section
    Parameters(const nlohmann::json& js, const std::string &section_name = "PATH_LENGTH_PARAMETERS")
      : Parameters() { assign_parameters(js, section_name); }

    /// @brief Copy assignment operator
    Parameters& operator=(const pathcalc::Parameters& other) = default;

    /// @brief Move assignment operator
    Parameters& operator=(pathcalc::Parameters&& other) noexcept = default;

    /// @brief Set instance name for identification
    /// @param[in] name_in New name string
    void set_name(const std::string &name_in) { name = name_in; }

    /// @brief Get instance name
    /// @return Current name string
    std::string get_name() const { return name; }

    /// @brief Load parameters from JSON configuration
    /// @param[in] js JSON object containing the runcard
    /// @param[in] section_name Section key in JSON to read from
    /// @throws std::runtime_error if required keys are missing
    void assign_parameters(const nlohmann::json& js, const std::string &section_name);

    /// @brief Inequality comparison operator
    /// @param[in] other Parameters to compare against
    /// @return true if any parameter value differs (name excluded)
    /// @note Name field is intentionally excluded from comparison
    bool operator!=(const Parameters& other) const;

    /// @brief Equality comparison operator
    /// @param[in] other Parameters to compare against
    /// @return true if all parameter values match (name excluded)
    /// @note Name field is intentionally excluded from comparison
    bool operator==(const Parameters& other) const { return !(*this != other); }
};
