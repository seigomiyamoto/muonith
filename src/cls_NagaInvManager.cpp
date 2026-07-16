/// @file cls_NagaInvManager.cpp
/// @brief Implementation of NagaInvManager class

#include "cls_NagaInvManager.hpp"
#include <stdexcept>
#include "spdlog_pch.hpp"
#include "ns_mylogger.hpp"

/// @brief Construct from existing parameter vector
/// @details Creates deep copies of all parameters via clone() method
NagaInvManager::NagaInvManager(const std::vector<NagaInvParameters>& vec_params_in) {
  vec_params.reserve(vec_params_in.size());
  for (const auto& prm : vec_params_in) {
    vec_params.emplace_back(prm.clone());  // deep copy via clone()
  }
}

NagaInvManager::NagaInvManager(const nlohmann::json& js, const std::string& section_name)
{
  load_from_json(js, section_name);
}

/// @brief Get parameter set by index with bounds checking
const NagaInvParameters& NagaInvManager::get_param(const int idx) const {
  if (idx < 0 || idx >= static_cast<int>(vec_params.size())) {
    THROW_ERROR("NagaInvManager::get_param: index {} out of range [0, {})", idx, vec_params.size());
  }
  return vec_params.at(idx);
}

/// @brief Load parameters from JSON object
void NagaInvManager::load_from_json(
  const nlohmann::json& js, const std::string& section_name)
{
  LOG_INFO("called.");

  // Verify section key exists
  if (!js.contains(section_name)) {
    THROW_ERROR("NagaInvManager::load_from_json: section key '{}' not found in JSON", section_name);
  }

  // Retrieve the JSON array
  const auto& paramsArray = js.at(section_name);

  // Verify the value is an array
  if (!paramsArray.is_array()) {
    THROW_ERROR("NagaInvManager::load_from_json: JSON value for key '{}' is not an array", section_name);
  }

  // Clear existing parameters
  vec_params.clear();

  // Load each parameter set from the array (empty section_name for individual items)
  for (const auto& js_param : paramsArray) {
    NagaInvParameters param(js_param, "");
    vec_params.emplace_back(std::move(param));
  }

  LOG_INFO("{} parameters loaded", vec_params.size());
}

/// @brief Factory method to create NagaInv solver instance
std::unique_ptr<NagaInv> NagaInvManager::create_nagainv_instance(
    const int idx
  , const std::string &name
  , const Grid3dVoxel& g3vox_prior
  , const Eigen::MatrixXf& mat_dNdD
  , const Eigen::VectorXf& vecxf_nmuon_obs
  , const Eigen::VectorXf& vecxf_nmuon_prior)
{
  LOG_INFO("Creating NagaInv instance with index: {}", idx);
  if (idx < 0 || idx >= static_cast<int>(vec_params.size())) {
    THROW_ERROR("NagaInvManager::create_nagainv_instance: index {} out of range [0, {})", idx, vec_params.size());
  }
  return std::make_unique<NagaInv>(
      name
    , g3vox_prior
    , vec_params.at(idx)
    , mat_dNdD
    , vecxf_nmuon_obs
    , vecxf_nmuon_prior
  );
}
