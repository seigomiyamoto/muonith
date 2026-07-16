/// @file cls_NagaInvManager.hpp
/// @brief Manager class for multiple NagaInv inversion parameter sets
/// @details
/// NagaInvManager provides a container for multiple NagaInvParameters configurations
/// and a factory method to create NagaInv instances from those configurations.
///
/// ## Responsibilities
/// - Load multiple inversion parameter sets from JSON configuration
/// - Store and provide indexed access to parameter sets
/// - Factory creation of NagaInv instances with specified parameter sets
///
/// ## Typical Workflow
/// 1. Create NagaInvManager from JSON configuration file
/// 2. Retrieve parameter count via vec_params.size()
/// 3. Access individual parameters via get_param(idx)
/// 4. Create NagaInv instances via create_nagainv_instance()
///
/// ## JSON Structure
/// The expected JSON structure is:
/// @code{.json}
/// {
///   "NAGAINV_PARAMETERS": [
///     { ...parameter set 0... },
///     { ...parameter set 1... }
///   ]
/// }
/// @endcode
///
/// ## Thread Safety
/// - Not thread-safe. External synchronization required for concurrent access.
///
/// @see NagaInvParameters for individual parameter set details
/// @see NagaInv for the inversion solver
#pragma once

#include <string>
#include <vector>
#include <memory>
#include <nlohmann/json.hpp>
#include "cls_NagaInvParameters.hpp"
#include "cls_NagaInv.hpp"

/// @brief Manager class for multiple NagaInvParameters configurations
/// @ingroup matrixClasses
///
/// @details
/// This class manages a collection of NagaInvParameters objects and provides
/// factory methods to create NagaInv solver instances.
///
/// ## Usage Example
/// @code
/// // Load from JSON file
/// nlohmann::json config = nlohmann::json::parse(std::ifstream("config.json"));
/// NagaInvManager manager(config);
///
/// // Create solver instance
/// auto solver = manager.create_nagainv_instance(
///     0,                    // parameter index
///     "inversion_run_1",    // name
///     prior_grid,           // Grid3dVoxel prior model
///     sensitivity_matrix,   // dN/dD matrix
///     observed_muons,       // observed muon counts
///     prior_muons           // prior muon counts
/// );
/// @endcode
class NagaInvManager {
public:
  /// @brief Default JSON section name for parameter arrays
  static constexpr const char* SECTION_NAME = "NAGAINV_PARAMETERS";

  /// @brief Collection of inversion parameter sets
  /// @note Direct access provided; use get_param() for bounds-checked access
  std::vector<NagaInvParameters> vec_params;

  /// @brief Construct from existing parameter vector
  /// @param[in] vec_params_in Source parameter vector to clone
  /// @note Creates deep copies via NagaInvParameters::clone() since copy is disabled
  NagaInvManager(const std::vector<NagaInvParameters>& vec_params_in);

  /// @brief Construct from JSON configuration
  /// @details calls \ref load_from_json(const nlohmann::json&, const std::string&)
  /// @param[in] js JSON object containing parameter array
  /// @param[in] section_name JSON key for the parameter array (default: "NAGAINV_PARAMETERS")
  /// @throws std::runtime_error if section_name key is missing or value is not an array
  NagaInvManager(const nlohmann::json& js, const std::string& section_name = SECTION_NAME);

  /// @brief Load parameters from JSON object
  /// @param[in] js JSON object containing parameter array
  /// @param[in] section_name JSON key for the parameter array (default: "NAGAINV_PARAMETERS")
  /// @throws std::runtime_error if section_name key is missing
  /// @throws std::runtime_error if JSON value is not an array
  /// @note Clears any existing parameters before loading
  void load_from_json(const nlohmann::json& js,
        const std::string& section_name = SECTION_NAME);

  /// @brief Get parameter set by index with bounds checking
  /// @param[in] idx Index of parameter set [0, vec_params.size())
  /// @return Immutable reference to NagaInvParameters
  /// @throws std::runtime_error if idx is out of range
  const NagaInvParameters& get_param(const int idx) const;

  /// @brief Factory method to create NagaInv solver instance
  /// @param[in] idx Index of parameter set to use [0, vec_params.size())
  /// @param[in] name Identifier name for the solver instance
  /// @param[in] g3vox_prior Prior 3D voxel grid model
  /// @param[in] mat_dNdD Sensitivity matrix (dN/dDensity)
  /// @param[in] vecxf_nmuon_obs Observed muon counts per path
  /// @param[in] vecxf_nmuon_prior Prior muon counts per path
  /// @return Unique pointer to newly created NagaInv instance
  /// @throws std::runtime_error if idx is out of range
  std::unique_ptr<NagaInv> create_nagainv_instance(
      const int idx
    , const std::string& name
    , const Grid3dVoxel& g3vox_prior
    , const Eigen::MatrixXf& mat_dNdD
    , const Eigen::VectorXf& vecxf_nmuon_obs
    , const Eigen::VectorXf& vecxf_nmuon_prior);
};
