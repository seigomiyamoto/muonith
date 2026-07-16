/// @file ns_param_constants.hpp
/// @brief Parameter constants shared across multiple classes.
/// @details
/// This namespace provides centralized management of parameter constants
/// such as density length (DL) and cosine of zenith angle (costhz) bounds.
/// Constants can be initialized from JSON configuration or use default values.
///
/// @note Thread-safety: Not thread-safe. Call init() once at startup
///       before any multi-threaded access to getters.
#pragma once

#include <nlohmann/json.hpp>
#include <string>

/// @namespace param_constants
/// @brief Parameter constants for muography calculations.
namespace param_constants {

/// @brief Initialize constants from JSON configuration.
/// @param[in] js JSON object containing configuration.
/// @param[in] section Section name in JSON (default: "PARAM_CONSTANTS").
/// @note If section is not found, default values are used with a warning.
/// @note Call this function once at startup, before other processing.
void init(const nlohmann::json& js, const std::string& section = "PARAM_CONSTANTS");

/// @brief Get minimum density length.
/// @return DL_min [kg/m^2]. Default: 5000.0 if init() not called.
/// @note 1 m.w.e. = 100 g/cm^2 = 1000 kg/m^2.
double DL_min();

/// @brief Get maximum density length.
/// @return DL_max [kg/m^2]. Default: 4999000.0 if init() not called.
double DL_max();

/// @brief Check whether DL_min was explicitly provided to init().
/// @return true if the JSON section contained "DL_min".
bool DL_min_is_set();

/// @brief Check whether DL_max was explicitly provided to init().
/// @return true if the JSON section contained "DL_max".
bool DL_max_is_set();

/// @brief Get minimum cosine of zenith angle.
/// @return costhz_min (dimensionless). Default: 0.0.
double costhz_min();

/// @brief Get maximum cosine of zenith angle.
/// @return costhz_max (dimensionless). Default: 1.0.
double costhz_max();

/// @brief Get maximum beam length (default).
/// @return BL_max_default [m]. Default: 5000.0 if init() not called.
double BL_max_default();

} // namespace param_constants
