// exec/make_peneflux.cpp
#include <iostream>
#include <filesystem>
#include <vector>
#include <cmath>
#include <limits>
#include <omp.h>
#include <nlohmann/json.hpp>
using json = nlohmann::json;
#include <regex>
#include <stdexcept>

#include "cls_Grid2dXYZ.hpp"
#include "cls_Grid1dXZ.hpp"
#include "cls_Grid1d.hpp"
#include "ns_myapp.hpp"
#include "ns_mymacro.hpp"
#include "ns_mylogger.hpp"
#include "cls_FluxTableBuilder.hpp"

/// @brief Build the penetrating flux table according to the config file.
/// @param js JSON object.
/// @param section_name Section name inside the JSON.
/// @throws std::exception If an error occurs while loading or writing.
void do_build_peneflux_table(const json &js, const std::string &section_name)
{
  START_FUNC();

  // Create the builder.
  FluxTableBuilder builder(js, section_name);

  // Report the range information.
  LOG_INFO("dlogE = {}", builder.get_dlogE());
  LOG_INFO("Rmin  = {}", builder.get_Rmin());
  LOG_INFO("Rmax  = {}", builder.get_Rmax());
  LOG_INFO("dR    = {}", builder.get_dR());

  // Build the cutoff table.
  LOG_INFO(" building log10Ecut table.");
  Grid1d g1_linR = builder.make_lin_range_axis();
  builder.build_log10Ecut_table(g1_linR);

  // Build the penetrating flux.
  LOG_INFO(" building log10(peneflux) table.");
  Grid2dXYZ g2_log_peneflux = builder.build_log_peneflux();

  // Build dFdR.
  LOG_INFO(" building dF/dR table from log10(peneflux) table.");
  Grid2dXYZ g2_dFdR = builder.build_dFdR_from_peneflux(g2_log_peneflux);

  // Write the output.
  builder.out_built_tables(g2_log_peneflux, g2_dFdR);

  END_FUNC();
}

int main(int argc, char** argv)
{
  START_FUNC();
  if (argc != 3) {
    throw std::runtime_error("Usage: " + std::string(argv[0]) + " <config.json> <SECTION_NAME>");
  }
hoge(1000000);
  const std::filesystem::path path_config(argv[1]);
  const json js = myapp::load_json(path_config);
  const std::string section_name = argv[2];

hoge(1100000);
  const auto [stdout_level, stderr_level, file_level] = mylogger::get_log_levels(js, "LOG_FILE");
  std::string basename = std::string("main_") + std::filesystem::path(path_config).stem().string();

hoge(1200000);
  std::filesystem::path log_dir("logs");
  std::filesystem::create_directories(log_dir);
  auto log_filepath = mylogger::generate_unique_log_path(log_dir, basename, ".log");
  bool archive_existing = mylogger::get_log_archive_existing(js, "LOG_FILE");

hoge(1300000);
  mylogger::g_logger = mylogger::create_logger(stdout_level, stderr_level, file_level, log_filepath, archive_existing);
  mylogger::configure_rate_limit(js, "LOG_FILE");
  LOG_INFO("Logging started to {}", log_filepath.string());
  LOG_INFO("Default max count for rate-limited logs: {}", mylogger::get_default_max_count());

  try {
    // Run.
    do_build_peneflux_table(js, section_name);
    mylogger::g_logger->flush();
  } catch (const std::exception& e) {
    LOG_CRITICAL("Unhandled exception: {}", e.what());
    mylogger::g_logger->flush();
    END_FUNC();
    return 1;
  }

  END_FUNC();
  return 0;
}