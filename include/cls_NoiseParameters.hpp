/// @file cls_NoiseParameters.hpp
/// @brief Noise parameter configuration for muon detectors
/// @details Defines NoiseParameters class for configuring various noise components
///          in the muon detection system. Noise parameters are loaded from JSON
///          configuration files.
///
/// @par Typical workflow:
/// 1. Construct NoiseParameters from JSON config (section-based)
/// 2. Query execution flag via get_tf_exec()
/// 3. Access noise ratio components (flux-proportional, SOT-proportional, user-defined)
/// 4. Optionally load user-defined noise flux tables via file paths
///
/// @par Thread-safety:
/// Not thread-safe. External synchronization required for concurrent access.
#pragma once
#include <string>
#include <filesystem>
#include <vector>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include "ns_mymacro.hpp"
#include "ns_mylogger.hpp"

//##############################################
//##############################################
/// @class NoiseParameters
/// @brief Configuration parameters for noise modeling in muon detection
/// @details This class holds noise configuration parameters read from JSON files.
///          Noise is modeled by two sources, each split into a deterministic
///          floor part and a Poisson-fluctuated part:
///          - Flux-proportional noise (scales with signal flux): floor + poisson
///          - SOT-proportional noise (angle-independent background): floor + poisson
///          Total noise = (flux_floor + SOT_floor) + Poisson(flux_poisson + SOT_poisson).
///          A separate user-defined noise source (from external flux tables) is reserved.
///
/// @par Usage example:
/// @code
/// nlohmann::json config = nlohmann::json::parse(config_file);
/// NoiseParameters noise_params(config, "noise_section");
/// if (noise_params.get_tf_exec()) {
///     double flux_poisson = noise_params.get_flux_proport_ratio_poisson();
///     // Apply noise modeling...
/// }
/// @endcode
///
/// @note The four ratio values default to 0.0, indicating "no contribution".
//##############################################
//##############################################
class NoiseParameters {
private:
  /// @brief Identifier name for this parameter set
  std::string name;

  /// @brief Execution flag (true = apply noise modeling, false = skip)
  bool tf_exec{false};

  /// @brief Flux-proportional noise ratio, deterministic floor part.
  ///        Added without Poisson fluctuation. Dimensionless, 0.0 = no contribution.
  double flux_proport_ratio_floor{0.0};

  /// @brief Flux-proportional noise ratio, Poisson-fluctuated part.
  ///        Summed into the Poisson bucket. Dimensionless, 0.0 = no contribution.
  double flux_proport_ratio_poisson{0.0};

  /// @brief SOT-proportional (angle-independent) noise ratio, deterministic floor part.
  ///        Added without Poisson fluctuation. Dimensionless, 0.0 = no contribution.
  double SOT_proport_noise_ratio_floor{0.0};

  /// @brief SOT-proportional (angle-independent) noise ratio, Poisson-fluctuated part.
  ///        Summed into the Poisson bucket. Dimensionless, 0.0 = no contribution.
  double SOT_proport_noise_ratio_poisson{0.0};

  /// @brief User-defined noise flux ratio (dimensionless, 0.0 = no contribution)
  double user_defined_noise_flux_ratio{0.0};

  /// @brief File paths to user-defined noise flux tables
  std::vector<std::filesystem::path> vec_path_user_defined_noise_flux;


public:
  //==================================================================
  /// @name Constructors and Destructor
  ///@{

  /// @brief Default constructor
  NoiseParameters() = default;

  /// @brief Copy constructor
  NoiseParameters(const NoiseParameters &other) = default;

  /// @brief Move constructor
  NoiseParameters(NoiseParameters &&other) noexcept = default;

  /// @brief Destructor
  ~NoiseParameters() = default;

  /// @brief Construct from JSON configuration
  /// @param[in] js JSON object containing the configuration
  /// @param[in] section_name Name of the section to read parameters from
  /// @throws std::runtime_error If required fields have invalid types
  NoiseParameters(const nlohmann::json& js, const std::string& section_name) {
    assign_parameters(js, section_name);
  }

  /// @brief Load parameters from JSON configuration
  /// @param[in] js JSON object containing the configuration
  /// @param[in] section_name Name of the section to read parameters from
  /// @throws std::runtime_error If required fields have invalid types
  /// @note Missing optional fields are logged as warnings but do not throw.
  void assign_parameters(const nlohmann::json& js, const std::string& section_name);

  ///@}

  //==================================================================
  /// @name Operators
  ///@{

  /// @brief Copy assignment operator
  NoiseParameters& operator=(const NoiseParameters& other) = default;

  /// @brief Inequality operator (compares all members except name)
  /// @param[in] other The instance to compare against
  /// @return true if any member differs, false otherwise
  bool operator!=(const NoiseParameters& other) const;

  /// @brief Equality operator (defined in terms of operator!=)
  /// @param[in] other The instance to compare against
  /// @return true if all members are equal, false otherwise
  bool operator==(const NoiseParameters& other) const {
    return !(*this != other);
  }

  ///@}

  //==================================================================
  /// @name Setters
  ///@{

  /// @brief Set the identifier name
  /// @param[in] name_in The name to set
  void set_name(const std::string& name_in) { name = name_in; }

  /// @brief Set the flux-proportional deterministic floor ratio
  /// @param[in] flux_proport_ratio_floor_in The ratio value (dimensionless)
  void set_flux_proport_ratio_floor(const double flux_proport_ratio_floor_in) {
    flux_proport_ratio_floor = flux_proport_ratio_floor_in;
  }

  /// @brief Set the flux-proportional Poisson-fluctuated ratio
  /// @param[in] flux_proport_ratio_poisson_in The ratio value (dimensionless)
  void set_flux_proport_ratio_poisson(const double flux_proport_ratio_poisson_in) {
    flux_proport_ratio_poisson = flux_proport_ratio_poisson_in;
  }

  /// @brief Set the SOT-proportional deterministic floor ratio
  /// @param[in] SOT_proport_noise_ratio_floor_in The ratio value (dimensionless)
  void set_SOT_proport_noise_ratio_floor(const double SOT_proport_noise_ratio_floor_in) {
    SOT_proport_noise_ratio_floor = SOT_proport_noise_ratio_floor_in;
  }

  /// @brief Set the SOT-proportional Poisson-fluctuated ratio
  /// @param[in] SOT_proport_noise_ratio_poisson_in The ratio value (dimensionless)
  void set_SOT_proport_noise_ratio_poisson(const double SOT_proport_noise_ratio_poisson_in) {
    SOT_proport_noise_ratio_poisson = SOT_proport_noise_ratio_poisson_in;
  }

  ///@}

  //==================================================================
  /// @name Getters
  ///@{

  /// @brief Get the identifier name
  /// @return The name of this parameter set
  std::string get_name() const { return name; }

  /// @brief Get the execution flag
  /// @return true if noise modeling should be applied, false to skip
  bool get_tf_exec() const { return tf_exec; }

  /// @brief Get the flux-proportional deterministic floor ratio
  /// @return The ratio value (dimensionless), 0.0 if not set
  double get_flux_proport_ratio_floor() const { return flux_proport_ratio_floor; }

  /// @brief Get the flux-proportional Poisson-fluctuated ratio
  /// @return The ratio value (dimensionless), 0.0 if not set
  double get_flux_proport_ratio_poisson() const { return flux_proport_ratio_poisson; }

  /// @brief Get the SOT-proportional deterministic floor ratio
  /// @return The ratio value (dimensionless), 0.0 if not set
  double get_SOT_proport_noise_ratio_floor() const { return SOT_proport_noise_ratio_floor; }

  /// @brief Get the SOT-proportional Poisson-fluctuated ratio
  /// @return The ratio value (dimensionless), 0.0 if not set
  double get_SOT_proport_noise_ratio_poisson() const { return SOT_proport_noise_ratio_poisson; }

  /// @brief Get the user-defined noise flux ratio
  /// @return The ratio value (dimensionless), or -1.0 if not set
  double get_user_defined_noise_flux_ratio() const { return user_defined_noise_flux_ratio; }

  /// @brief Get the file paths for user-defined noise flux tables
  /// @return Const reference to the vector of file paths
  const std::vector<std::filesystem::path>& get_vec_path_user_defined_noise_flux() const {
    return vec_path_user_defined_noise_flux;
  }

  ///@}
};