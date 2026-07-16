/// @file ns_param_util.hpp
/// @brief Parameter loading utilities for JSON configuration
/// @details
/// Provides utility functions for reading parameters from JSON configuration files.
/// These functions are used by various *Parameters classes to simplify JSON parsing
/// and reduce code duplication.
///
/// ## Typical usage
/// @code
/// #include "ns_param_util.hpp"
///
/// void SomeParameters::assign_parameters(const nlohmann::json& js,
///                                        const std::string& section_name) {
///   param_util::read_json_value(js, section_name, TOSTRING(param_name), param_name);
///   param_util::read_json_path(js, section_name, TOSTRING(path_param), path_param);
/// }
/// @endcode
///
/// @namespace param_util
/// @brief Utilities for JSON parameter loading
#pragma once

#include <string>
#include <filesystem>
#include <sstream>
#include <nlohmann/json.hpp>
#include <Eigen/Dense>
#include "ns_mymacro.hpp"
#include "ns_mylogger.hpp"

namespace fs = std::filesystem;

/// @namespace param_util
/// @brief Utilities for JSON parameter loading
namespace param_util {

  /// @brief Read a value from JSON configuration
  /// @tparam T The type of the value to read (e.g., int, double, bool, std::string)
  /// @param js The JSON object containing the configuration
  /// @param section_name The section name within the JSON object
  /// @param key_str The key string for the parameter
  /// @param value Reference to the variable where the value will be stored
  /// @details
  /// If the section and key exist in the JSON, the value is read and assigned.
  /// If they don't exist, the value remains unchanged (preserving default values).
  /// @note This is a template function defined in the header for inline expansion
  template<typename T>
  void read_json_value(const nlohmann::json& js,
                       const std::string& section_name,
                       const char* key_str,
                       T& value) {
    if(js.contains(section_name) && js.at(section_name).contains(key_str)) {
      value = js.at(section_name).at(key_str).get<T>();
    }
  }

  /// @brief Read a required value from JSON configuration
  /// @tparam T The type of the value to read (e.g., int, double, bool, std::string)
  /// @param js The JSON object containing the configuration
  /// @param section_name The section name within the JSON object
  /// @param key_str The key string for the parameter
  /// @param value Reference to the variable where the value will be stored
  /// @details
  /// If the section or key is missing, throws with a formatted error message.
  template<typename T>
  void read_json_value_required(const nlohmann::json& js,
                                const std::string& section_name,
                                const char* key_str,
                                T& value) {
    if(!js.contains(section_name)) {
      THROW_ERROR("param_util::read_json_value_required: section '{}' not found", section_name);
    }
    const auto& section = js.at(section_name);
    if(!section.contains(key_str)) {
      THROW_ERROR("param_util::read_json_value_required: key '{}' missing in section '{}'",
                  key_str, section_name);
    }
    value = section.at(key_str).get<T>();
  }

  /// @brief Read a value from JSON with warning log if key is missing
  /// @tparam T The type of the value to read (e.g., int, double, bool, std::string)
  /// @param js The JSON object containing the configuration
  /// @param key_str The key string for the parameter
  /// @param value Reference to the variable where the value will be stored
  /// @param class_name Name of the class for logging (e.g., "VerticalCylinderParameters")
  /// @param entry_name Name of the entry instance for logging
  /// @details
  /// If the key exists in the JSON, the value is read and assigned.
  /// If the key doesn't exist, a warning is logged and the value remains unchanged.
  /// @note This version does NOT use section_name, only direct key lookup
  template<typename T>
  void read_json_value_warn(const nlohmann::json& js,
                            const char* key_str,
                            T& value,
                            const std::string& class_name,
                            const std::string& entry_name) {
    if(js.contains(key_str)) {
      value = js.at(key_str).get<T>();
    } else {
      std::ostringstream oss;
      oss << value;
      LOG_WARN("{}: key '{}' missing for '{}'; using default {}.",
               class_name, key_str, entry_name, oss.str());
    }
  }

  /// @brief Read a filesystem path from JSON configuration
  /// @param js The JSON object containing the configuration
  /// @param section_name The section name within the JSON object
  /// @param key_str The key string for the path parameter
  /// @param value Reference to the fs::path where the path will be stored
  /// @details
  /// Reads a string from JSON and converts it to fs::path.
  /// If the section and key don't exist, the value remains unchanged.
  inline void read_json_path(const nlohmann::json& js,
                             const std::string& section_name,
                             const char* key_str,
                             fs::path& value) {
    if(js.contains(section_name) && js.at(section_name).contains(key_str)) {
      value = fs::path(js.at(section_name).at(key_str).get<std::string>());
    }
  }

  /// @brief Read a required filesystem path from JSON configuration
  /// @param js The JSON object containing the configuration
  /// @param section_name The section name within the JSON object
  /// @param key_str The key string for the path parameter
  /// @param value Reference to the fs::path where the path will be stored
  /// @details
  /// If the section or key is missing, throws with a formatted error message.
  inline void read_json_path_required(const nlohmann::json& js,
                                      const std::string& section_name,
                                      const char* key_str,
                                      fs::path& value) {
    if(!js.contains(section_name)) {
      THROW_ERROR("param_util::read_json_path_required: section '{}' not found", section_name);
    }
    const auto& section = js.at(section_name);
    if(!section.contains(key_str)) {
      THROW_ERROR("param_util::read_json_path_required: key '{}' missing in section '{}'",
                  key_str, section_name);
    }
    value = fs::path(section.at(key_str).get<std::string>());
  }

  /// @brief Read a value from JSON without section_name (direct key lookup)
  /// @tparam T The type of the value to read (e.g., int, double, bool, std::string)
  /// @param js The JSON object containing the configuration
  /// @param key_str The key string for the parameter
  /// @param value Reference to the variable where the value will be stored
  /// @details
  /// If the key exists in the JSON, the value is read and assigned.
  /// If it doesn't exist, the value remains unchanged (preserving default values).
  /// This version does NOT use section_name, only direct key lookup.
  template<typename T>
  void read_json_direct(const nlohmann::json& js,
                        const char* key_str,
                        T& value) {
    if(js.contains(key_str)) {
      value = js.at(key_str).get<T>();
    }
  }

  /// @brief Read Eigen::Vector3d components from JSON (direct key lookup)
  /// @param js JSON object containing the configuration
  /// @param x_key Key for x component
  /// @param y_key Key for y component
  /// @param z_key Key for z component
  /// @param vec Reference to Eigen::Vector3d to store the result
  /// @details Reads three scalar values and assigns them to vec.x(), vec.y(), vec.z()
  inline void read_vector3d_direct(const nlohmann::json& js,
                                   const char* x_key,
                                   const char* y_key,
                                   const char* z_key,
                                   Eigen::Vector3d& vec) {
    double x = 0.0, y = 0.0, z = 0.0;
    read_json_direct(js, x_key, x);
    read_json_direct(js, y_key, y);
    read_json_direct(js, z_key, z);
    vec.x() = x;
    vec.y() = y;
    vec.z() = z;
  }

} // namespace param_util
