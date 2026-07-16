/// @file cls_NoiseParameters.cpp
/// @brief Implementation of NoiseParameters class
/// @details Contains parameter loading from JSON and comparison operators.

#include "cls_NoiseParameters.hpp"
#include "ns_myapp.hpp"

bool NoiseParameters::operator!=(const NoiseParameters& other) const {
  #ifdef NODEBUG
    if (tf_exec != other.tf_exec) return true;
    if (flux_proport_ratio_floor != other.flux_proport_ratio_floor) return true;
    if (flux_proport_ratio_poisson != other.flux_proport_ratio_poisson) return true;
    if (SOT_proport_noise_ratio_floor != other.SOT_proport_noise_ratio_floor) return true;
    if (SOT_proport_noise_ratio_poisson != other.SOT_proport_noise_ratio_poisson) return true;
    if (user_defined_noise_flux_ratio != other.user_defined_noise_flux_ratio) return true;
    if (vec_path_user_defined_noise_flux.size() != other.vec_path_user_defined_noise_flux.size()) return true;
    for (size_t i = 0; i < vec_path_user_defined_noise_flux.size(); ++i) {
      if (vec_path_user_defined_noise_flux[i] != other.vec_path_user_defined_noise_flux[i]) return true;
    }
  #else
    if (tf_exec != other.tf_exec) { LOG_DEBUG_TAG("NoiseParameters", "tf_exec differs: {} vs {}", tf_exec, other.tf_exec); return true; }
    if (flux_proport_ratio_floor != other.flux_proport_ratio_floor) { LOG_DEBUG_TAG("NoiseParameters", "flux_proport_ratio_floor differs: {} vs {}", flux_proport_ratio_floor, other.flux_proport_ratio_floor); return true; }
    if (flux_proport_ratio_poisson != other.flux_proport_ratio_poisson) { LOG_DEBUG_TAG("NoiseParameters", "flux_proport_ratio_poisson differs: {} vs {}", flux_proport_ratio_poisson, other.flux_proport_ratio_poisson); return true; }
    if (SOT_proport_noise_ratio_floor != other.SOT_proport_noise_ratio_floor) { LOG_DEBUG_TAG("NoiseParameters", "SOT_proport_noise_ratio_floor differs: {} vs {}", SOT_proport_noise_ratio_floor, other.SOT_proport_noise_ratio_floor); return true; }
    if (SOT_proport_noise_ratio_poisson != other.SOT_proport_noise_ratio_poisson) { LOG_DEBUG_TAG("NoiseParameters", "SOT_proport_noise_ratio_poisson differs: {} vs {}", SOT_proport_noise_ratio_poisson, other.SOT_proport_noise_ratio_poisson); return true; }
    if (user_defined_noise_flux_ratio != other.user_defined_noise_flux_ratio) { LOG_DEBUG_TAG("NoiseParameters", "user_defined_noise_flux_ratio differs: {} vs {}", user_defined_noise_flux_ratio, other.user_defined_noise_flux_ratio); return true; }
    if (vec_path_user_defined_noise_flux != other.vec_path_user_defined_noise_flux) { LOG_DEBUG_TAG("NoiseParameters", "vec_path_user_defined_noise_flux differs (size: {} vs {})", vec_path_user_defined_noise_flux.size(), other.vec_path_user_defined_noise_flux.size()); return true; }
    for (size_t i = 0; i < vec_path_user_defined_noise_flux.size(); ++i) {
      if (vec_path_user_defined_noise_flux[i] != other.vec_path_user_defined_noise_flux[i]) { LOG_DEBUG_TAG("NoiseParameters", "vec_path_user_defined_noise_flux[{}] differs: {} vs {}", i, vec_path_user_defined_noise_flux[i], other.vec_path_user_defined_noise_flux[i]); return true; }
    }
  #endif
  return false;
}



void NoiseParameters::assign_parameters(
  const nlohmann::json& js, const std::string& section_name)
{
  if (!js.contains(section_name)) {
    LOG_WARN("Section '{}' not found in JSON config.", section_name);
    return;
  }

  const auto& sec = js.at(section_name);
  std::string key;

  // name
  key = TOSTRING(name);
  if (sec.contains(key)) {
    name = sec.at(key).get<std::string>();
    LOG_INFO_TAG("NoiseParameters", "name: {}", name);
  } else {
    LOG_WARN_TAG("NoiseParameters", "name not found in section '{}'", section_name);
  }

  // tf_exec
  key = TOSTRING(tf_exec);
  if (sec.contains(key)) {
    tf_exec = sec.at(key).get<bool>();
    LOG_INFO_TAG("NoiseParameters", "tf_exec: {}", tf_exec);
  } else {
    LOG_WARN_TAG("NoiseParameters", "tf_exec not found in section '{}'", section_name);
  }

  // flux_proport_ratio_floor
  key = TOSTRING(flux_proport_ratio_floor);
  if (sec.contains(key)) {
    flux_proport_ratio_floor = sec.at(key).get<double>();
  } else {
    LOG_WARN_TAG("NoiseParameters", "flux_proport_ratio_floor not found in section '{}'", section_name);
  }

  // flux_proport_ratio_poisson
  key = TOSTRING(flux_proport_ratio_poisson);
  if (sec.contains(key)) {
    flux_proport_ratio_poisson = sec.at(key).get<double>();
  } else {
    LOG_WARN_TAG("NoiseParameters", "flux_proport_ratio_poisson not found in section '{}'", section_name);
  }

  // SOT_proport_noise_ratio_floor
  key = TOSTRING(SOT_proport_noise_ratio_floor);
  if (sec.contains(key)) {
    SOT_proport_noise_ratio_floor = sec.at(key).get<double>();
  } else {
    LOG_WARN_TAG("NoiseParameters", "SOT_proport_noise_ratio_floor not found in section '{}'", section_name);
  }

  // SOT_proport_noise_ratio_poisson
  key = TOSTRING(SOT_proport_noise_ratio_poisson);
  if (sec.contains(key)) {
    SOT_proport_noise_ratio_poisson = sec.at(key).get<double>();
  } else {
    LOG_WARN_TAG("NoiseParameters", "SOT_proport_noise_ratio_poisson not found in section '{}'", section_name);
  }

  // user_defined_noise_flux_ratio
  key = TOSTRING(user_defined_noise_flux_ratio);
  if (sec.contains(key)) {
    user_defined_noise_flux_ratio = sec.at(key).get<double>();
  } else {
    LOG_WARN_TAG("NoiseParameters", "user_defined_noise_flux_ratio not found in section '{}'", section_name);
  }

  // pathin_user_defined_noise_distribution
  key = TOSTRING(pathin_user_defined_noise_distribution);
  if (sec.contains(key)) {
    const auto& arr = sec.at(key);
    if (!arr.is_array()) {
      THROW_ERROR("NoiseParameters::assign_parameters: '{}' must be an array in section '{}'", key, section_name);
    }
    vec_path_user_defined_noise_flux.clear();
    for (const auto& item : arr) {
      if (!item.is_string()) {
        THROW_ERROR("NoiseParameters::assign_parameters: All elements of '{}' must be strings in section '{}'", key, section_name);
      }
      vec_path_user_defined_noise_flux.emplace_back(item.get<std::string>());
    }
    LOG_INFO_TAG("NoiseParameters", "Loaded {} user-defined noise file paths.", vec_path_user_defined_noise_flux.size());
  } else {
    LOG_INFO_TAG("NoiseParameters", "'{}' not found in section '{}'", key, section_name);
  }
}