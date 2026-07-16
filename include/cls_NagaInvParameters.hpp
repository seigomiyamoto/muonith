/// @file cls_NagaInvParameters.hpp
/// @brief NagaInv inversion parameter configuration
/// @details
/// Defines the NagaInvParameters class for configuring all inversion solver
/// settings and regularization parameters. This class serves as the central
/// configuration container for muon tomography inversion solvers.
///
/// ## Typical workflow
/// 1. Load JSON configuration file containing inversion parameters
/// 2. Construct NagaInvParameters from the JSON with section name
/// 3. Query individual parameters via accessor methods (get_corr_length(), etc.)
/// 4. Pass to NagaInv solver for inversion computation
///
/// ## Key parameters
/// - **corr_length**: Spatial correlation length for regularization [meters]
/// - **sigma_rho**: Standard deviation of density prior [kg/m^3]
/// - **uniform_prior_density**: Initial uniform density assumption [kg/m^3]
/// - **nmuon_thres**: Minimum muon count threshold for data acceptance
///
/// ## Units
/// - Density: kg/m^3
/// - Length: meters
/// - Muon counts: dimensionless (particle counts)
///
/// ## Thread safety
/// - Thread-safe: Read-only accessors (get_*() methods) after construction
/// - Not thread-safe: assign_parameters() must be called before concurrent access
/// - Recommendation: Complete all parameter loading in single thread before sharing
///
/// @ingroup parameterClasses
#pragma once

#include <string>
#include <cstdio>
#include <nlohmann/json.hpp>
#include "ns_myapp.hpp"

//----------------------------------------------------------------------------------------------------
/// @class NagaInvParameters
/// @brief Configuration parameters for NagaInv muon tomography inversion solver
/// @ingroup parameterClasses
///
/// @details
/// This class stores and provides access to all configuration parameters needed
/// for running the NagaInv inversion algorithm. Parameters are loaded from a
/// JSON configuration file and accessed through type-safe getter methods.
///
/// ## Responsibilities
/// - Load and validate inversion parameters from JSON configuration
/// - Provide type-safe access to parameters with default fallback values
/// - Support explicit deep copying via clone() method
/// - Maintain parameter consistency throughout inversion workflow
///
/// ## Class invariants
/// - json_params contains valid JSON after construction from JSON or after assign_parameters()
/// - Default values are returned when keys are not present in JSON
/// - All getter methods are const and do not modify internal state
/// - Type mismatches in get_param<T>() throw std::runtime_error
///
/// ## Copy semantics
/// Copy constructor and copy assignment are deleted to prevent accidental deep copies.
/// Use clone() for explicit copying when needed.
///
/// ## Usage example
/// @code
/// nlohmann::json config = nlohmann::json::parse(config_file);
/// NagaInvParameters params(config, "inversion");
/// double corr_len = params.get_corr_length();  // returns value or default
/// bool exec = params.get_tf_exec();            // type-safe boolean access
/// @endcode
///
/// @note Thread-safe for concurrent read access after parameter assignment completes
//----------------------------------------------------------------------------------------------------
class NagaInvParameters {
  public:

  /// @name Default parameter values
  /// @brief Compile-time default constants for inversion parameters
  /// @{

  /// @brief Default muon count threshold (minimum muons for valid data)
  static constexpr double nmuon_thres_init = 0.0;

  /// @brief Default under-threshold muon count value
  static constexpr double nmuon_under_thres_init = 0.0;

  /// @brief Default uniform prior density [kg/m^3]
  static constexpr double uniform_prior_density_init = 1000.0;

  /// @brief Default spatial correlation length [meters]
  static constexpr double corr_length_init = 100.0;

  /// @brief Default standard deviation of density prior [kg/m^3]
  static constexpr double sigma_rho_init = 300.0;

  /// @brief Default diagonal self-error standard deviation
  /// @note Negative value indicates auto-computation
  static constexpr double sigma_rho_diag_init = -1.0;

  /// @brief Default horizontal (xy) correlation length [meters]
  /// @note Negative value indicates fallback to corr_length (isotropic)
  static constexpr double corr_length_xy_init = -1.0;

  /// @brief Default vertical (z) correlation length [meters]
  /// @note Negative value indicates fallback to corr_length (isotropic)
  static constexpr double corr_length_z_init = -1.0;

  /// @}

  /// @brief JSON object storing all loaded parameters
  nlohmann::json json_params;

  /// @brief Default constructor
  /// @details Creates an empty parameter set; call assign_parameters() to load values
  NagaInvParameters() {};

  /// @brief Copy constructor (deleted)
  /// @details Use clone() for explicit copying
  NagaInvParameters(const NagaInvParameters& org) = delete;

  /// @brief Create a deep copy of this parameter set
  /// @return New NagaInvParameters with copied json_params
  NagaInvParameters clone() const;

  /// @brief Copy assignment operator (deleted)
  /// @details Use clone() for explicit copying
  NagaInvParameters& operator=(const NagaInvParameters& org) = delete;

  /// @brief Move constructor
  /// @param[in] other Source object to move from
  NagaInvParameters(NagaInvParameters&& other) noexcept
    : json_params(std::move(other.json_params)) {};

  /// @brief Move assignment operator
  /// @param[in] other Source object to move from
  /// @return Reference to this object
  NagaInvParameters& operator=(NagaInvParameters&& other) noexcept {
    if (this != &other) {
      json_params = std::move(other.json_params);
    }
    return *this;
  };

  /// @brief Construct from JSON configuration
  /// @param[in] js JSON object containing configuration
  /// @param[in] section_name Key to look up in js for this parameter section
  /// @throws std::runtime_error If section_name is not found in js
  NagaInvParameters(const nlohmann::json& js, const std::string &section_name)
  { assign_parameters(js, section_name); };

  /// @brief Destructor
  ~NagaInvParameters() = default;

  /// @brief Print json_params to FILE stream
  /// @param[in] fout Output file pointer; if nullptr, defaults to stdout
  /// @note Output format: "NagaInvParameters::json_params:" header followed by JSON dump
  void print_json_params(FILE* fout) const;

  /// @brief Print json_params via spdlog at specified level
  /// @param[in] log_level spdlog level for output; spdlog::level::off suppresses output
  /// @note Output format: "NagaInvParameters::json_params:" header followed by JSON dump
  void print_json_params(const spdlog::level::level_enum& log_level) const;

  /// @brief Print json_params to stderr
  /// @note Convenience wrapper for print_json_params(stderr)
  void print_json_params() const {
    print_json_params(stderr);
  };

  /// @brief Load parameters from JSON configuration
  /// @param[in] js JSON object containing configuration
  /// @param[in] section_name Key to look up in js; if empty string, uses js directly
  /// @throws std::runtime_error If section_name is not empty and not found in js
  /// @note Thread-unsafe: do not call concurrently with any other methods
  /// @note Overwrites existing json_params completely; previous values are lost
  void assign_parameters(
    const nlohmann::json& js, const std::string &section_name);

  /// @brief Get the name identifier of this parameter set
  /// @return Name string from JSON, or "NagaInvParameters" if not set
  std::string get_name() const {
    return json_params.value("name", "NagaInvParameters");
  };

  /// @brief Get execution flag
  /// @return True if this inversion should be executed
  bool get_tf_exec() const {
    return json_params.value("tf_exec", false);
  };

  /// @brief Get signal Poisson fluctuation flag
  /// @return True if signal counts should carry Poisson error
  bool get_tf_signal_poisson() const {
    return json_params.value("tf_signal_poisson", false);
  };

  /// @brief Get efficiency-uncertainty C_N diagonal flag
  /// @return True if efficiency uncertainty is added as analytic variance to the C_N diagonal
  bool get_tf_eff_cn_diag() const {
    return json_params.value("tf_eff_cn_diag", false);
  };

  /// @brief Get efficiency-uncertainty independence flag (within a merge group)
  /// @return True: per-element uncertainties are independent (sum of squares).
  ///         False (default): common within the group (square of the sum).
  bool get_tf_eff_cn_diag_independent() const {
    return json_params.value("tf_eff_cn_diag_independent", false);
  };

  /// @brief Get chi-square / p_eff self-evaluation flag
  /// @return True if chi2/ndf and the effective number of parameters p_eff = trace(R)
  ///         are computed and reported. Default false (no extra computation).
  bool get_tf_calc_chi2ndf() const {
    return json_params.value("tf_calc_chi2ndf", false);
  };

  /// @brief Get muon count threshold
  /// @return Minimum muon count for valid data acceptance
  double get_nmuon_thres() const {
    return json_params.value("nmuon_thres", nmuon_thres_init);
  };

  /// @brief Get under-threshold muon count value
  /// @return Value to use when muon count is below threshold
  double get_nmuon_under_thres() const {
    return json_params.value("nmuon_under_thres", nmuon_under_thres_init);
  };

  /// @brief Get uniform prior density
  /// @return Initial uniform density assumption [kg/m^3]
  double get_uniform_prior_density() const {
    return json_params.value("uniform_prior_density", uniform_prior_density_init);
  };

  /// @brief Get shell density for upper region
  /// @param[in] fallback Value used when key is absent in JSON [kg/m^3]
  /// @return Shell density for upper region [kg/m^3]
  /// @note Falls back to the provided fallback (typically the swept uniform_prior_density).
  double get_shell_density_upper(double fallback) const {
    return json_params.value("shell_density_upper", fallback);
  };

  /// @brief Get shell density for lower region
  /// @param[in] fallback Value used when key is absent in JSON [kg/m^3]
  /// @return Shell density for lower region [kg/m^3]
  /// @note Falls back to the provided fallback (typically the swept uniform_prior_density).
  double get_shell_density_lower(double fallback) const {
    return json_params.value("shell_density_lower", fallback);
  };

  /// @brief Get shell density for lateral region
  /// @param[in] fallback Value used when key is absent in JSON [kg/m^3]
  /// @return Shell density for lateral region [kg/m^3]
  /// @note Falls back to the provided fallback (typically the swept uniform_prior_density).
  double get_shell_density_lateral(double fallback) const {
    return json_params.value("shell_density_lateral", fallback);
  };

  /// @brief Get spatial correlation length
  /// @return Correlation length for regularization [meters]
  double get_corr_length() const {
    return json_params.value("corr_length", corr_length_init);
  };

  /// @brief Get density prior standard deviation
  /// @return Standard deviation of density prior [kg/m^3]
  double get_sigma_rho() const {
    return json_params.value("sigma_rho", sigma_rho_init);
  };

  /// @brief Get diagonal self-error standard deviation
  /// @return Diagonal error term; negative value indicates auto-computation
  double get_sigma_rho_diag() const {
    return json_params.value("sigma_rho_diag", sigma_rho_diag_init);
  };

  /// @brief Get horizontal (xy) correlation length [meters]
  /// @return corr_length_xy if set (>0), otherwise falls back to corr_length
  double get_corr_length_xy() const {
    const double v = json_params.value("corr_length_xy", corr_length_xy_init);
    return (v > 0.0) ? v : get_corr_length();
  };

  /// @brief Get vertical (z) correlation length [meters]
  /// @return corr_length_z if set (>0), otherwise falls back to corr_length
  double get_corr_length_z() const {
    const double v = json_params.value("corr_length_z", corr_length_z_init);
    return (v > 0.0) ? v : get_corr_length();
  };

  /// @brief Get anisotropic covariance formula type
  /// @return "separable" or "ellipsoidal" (default: "separable")
  std::string get_aniso_cov_type() const {
    return json_params.value("aniso_cov_type", std::string("separable"));
  };

  /// @brief Check if anisotropic covariance mode is enabled
  /// @return True if tf_aniso is explicitly set to true in JSON,
  ///         or if corr_length_xy/corr_length_z is set to a positive value
  bool is_anisotropic() const {
    if (json_params.contains("tf_aniso")) {
      return json_params.value("tf_aniso", false);
    }
    return json_params.value("corr_length_xy", corr_length_xy_init) > 0.0
        || json_params.value("corr_length_z",  corr_length_z_init)  > 0.0;
  };

  /// @brief Get arbitrary parameter by key with type checking
  /// @tparam T Expected type (double, bool, int, or std::string)
  /// @param[in] key Parameter key name
  /// @param[in] default_value Value to return if key not found
  /// @return Parameter value if found, otherwise default_value
  /// @throws std::runtime_error If key exists but type does not match T
  /// @note Supported types: double (any number), bool, int (integer numbers only), std::string
  /// @note Type checking: double accepts any JSON number, int requires integer, bool requires boolean
  /// @note Example: `double val = params.get_param<double>("my_key", 1.0);`
  template<typename T>
  T get_param(const std::string &key, const T &default_value) const {
    if (!json_params.contains(key)) {
      return default_value;
    }
    const auto& value = json_params[key];

    // Type validation for common types
    if constexpr (std::is_same_v<T, double>) {
      if (!value.is_number()) {
        THROW_ERROR("NagaInvParameters::get_param: invalid type for key '{}' (expected number)", key);
      }
    } else if constexpr (std::is_same_v<T, bool>) {
      if (!value.is_boolean()) {
        THROW_ERROR("NagaInvParameters::get_param: invalid type for key '{}' (expected boolean)", key);
      }
    } else if constexpr (std::is_same_v<T, int>) {
      if (!value.is_number_integer()) {
        THROW_ERROR("NagaInvParameters::get_param: invalid type for key '{}' (expected integer)", key);
      }
    } else if constexpr (std::is_same_v<T, std::string>) {
      if (!value.is_string()) {
        THROW_ERROR("NagaInvParameters::get_param: invalid type for key '{}' (expected string)", key);
      }
    }
    return value.get<T>();
  };

};
