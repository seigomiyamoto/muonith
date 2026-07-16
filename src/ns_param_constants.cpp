/// @file ns_param_constants.cpp
/// @brief Implementation of parameter constants management.
#include "ns_param_constants.hpp"
#include "ns_mylogger.hpp"

namespace {
// Private static variables with default values
double s_DL_min         = 5000.0;      // [kg/m^2]
double s_DL_max         = 4999000.0;   // [kg/m^2]
bool   s_DL_min_set     = false;       // true when DL_min came from JSON
bool   s_DL_max_set     = false;       // true when DL_max came from JSON
double s_costhz_min     = 0.0;         // dimensionless
double s_costhz_max     = 1.0;         // dimensionless
double s_BL_max_default = 5000.0;      // [m] beam length
} // anonymous namespace

namespace param_constants {

void init(const nlohmann::json& js, const std::string& section) {
  if (!js.contains(section)) {
    LOG_INFO("param_constants::init: Section '{}' not found. Using defaults.", section);
    return;
  }

  const auto& sec = js.at(section);

  if (sec.contains("DL_min")) {
    s_DL_min = sec.at("DL_min").get<double>();
    s_DL_min_set = true;
    LOG_INFO("param_constants: DL_min = {} [kg/m^2]", s_DL_min);
  }

  if (sec.contains("DL_max")) {
    s_DL_max = sec.at("DL_max").get<double>();
    s_DL_max_set = true;
    LOG_INFO("param_constants: DL_max = {} [kg/m^2]", s_DL_max);
  }

  if (sec.contains("costhz_min")) {
    s_costhz_min = sec.at("costhz_min").get<double>();
    LOG_INFO("param_constants: costhz_min = {}", s_costhz_min);
  }

  if (sec.contains("costhz_max")) {
    s_costhz_max = sec.at("costhz_max").get<double>();
    LOG_INFO("param_constants: costhz_max = {}", s_costhz_max);
  }

  if (sec.contains("BL_max_default")) {
    s_BL_max_default = sec.at("BL_max_default").get<double>();
    LOG_INFO("param_constants: BL_max_default = {} [m]", s_BL_max_default);
  }
}

double DL_min()         { return s_DL_min; }
double DL_max()         { return s_DL_max; }
bool DL_min_is_set()    { return s_DL_min_set; }
bool DL_max_is_set()    { return s_DL_max_set; }
double costhz_min()     { return s_costhz_min; }
double costhz_max()     { return s_costhz_max; }
double BL_max_default() { return s_BL_max_default; }

} // namespace param_constants
