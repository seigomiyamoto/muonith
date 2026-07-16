/// @file cls_NagaInvParameters.cpp
/// @brief Implementation of NagaInvParameters class
/// @details Provides parameter loading, cloning, and output functionality for NagaInv
/// inversion solver configuration. All implementations preserve thread-safety guarantees
/// documented in the header.

#include "cls_NagaInvParameters.hpp"
#include "spdlog_pch.hpp"
#include "ns_mylogger.hpp"
#include "ns_mymacro.hpp"

void NagaInvParameters::assign_parameters(
  const nlohmann::json& js, const std::string &section_name)
{
  LOG_INFO("called.");
  try {
    // If section_name is provided, extract that subsection from JSON
    if (!section_name.empty()) {
      if (js.find(section_name) == js.end()) {
        THROW_ERROR("NagaInvParameters::assign_parameters: section '{}' not found in JSON",
                    section_name);
      }
      json_params = js.at(section_name);
    } else {
      json_params = js;
    }

  } catch (const std::exception &e) {
    LOG_ERROR("error: {}", e.what());
    throw;
  }
}

NagaInvParameters NagaInvParameters::clone() const {
  NagaInvParameters copy;

  // Deep copy json_params
  copy.json_params = this->json_params;

  return copy;
}

void NagaInvParameters::print_json_params(FILE* fout) const {
  if (fout == nullptr) {
    fout = stdout;
  }
  fprintf(fout, "NagaInvParameters::json_params:\n");
  myapp::dump_json(json_params, {}, "", fout);
}

void NagaInvParameters::print_json_params(const spdlog::level::level_enum& log_level) const {
  if (log_level == spdlog::level::off) {
    return;
  }
  mylogger::g_logger->log(log_level, "NagaInvParameters::json_params:");
  myapp::dump_json(json_params, {}, "", log_level);
}
