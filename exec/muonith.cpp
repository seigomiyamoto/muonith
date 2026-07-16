/// @file muonith.cpp
/// @brief Parameter sweep executable for 3D density reconstruction.
/// @details
/// Uses exemdl::pipeline for modular execution and exemdl::sweep for
/// parameter sweep orchestration. Reads NAGAINV_PARAM_SWEEP section
/// from JSON5 configuration.

#include <string>
#include <filesystem>
#include <nlohmann/json.hpp>
#include <getopt.h>
#include <cstdlib>

#include "exemdl_sweep.hpp"
#include "ns_myapp.hpp"
#include "ns_mylogger.hpp"
#include "ns_seed.hpp"
#include "ns_mymacro.hpp"
#include "ns_param_constants.hpp"
#include "spdlog_pch.hpp"

namespace fs = std::filesystem;
using json = nlohmann::json;

static std::string fname_json;
static bool seed_given = false;
static unsigned seed_value = 42;
static fs::path resume_dir;
static int end_stage = -1;  // sentinel: resolve from JSON5 or default(8)

void parse_args(int argc, char** argv)
{
  const struct option longopts[] = {
      {"json",      required_argument, nullptr, 'j'}
    , {"seed",      required_argument, nullptr, 's'}
    , {"resume",    required_argument, nullptr, 'r'}
    , {"end-stage", required_argument, nullptr, 'E'}
    , {nullptr, 0, nullptr, 0}
  };

  int opt;
  while ((opt = getopt_long(argc, argv, "j:s:r:E:", longopts, nullptr)) != -1) {
    switch (opt) {
      case 'j': fname_json = optarg; break;
      case 's': seed_given = true; seed_value = std::stoul(optarg); break;
      case 'r': resume_dir = optarg; break;
      case 'E': end_stage = exemdl::sweep::resolve_end_stage(optarg); break;
      default:
        std::cerr << "Usage: " << argv[0]
                  << " --json <file> [--seed <val>] [--resume <checkpoint_dir>]"
                  << " [--end-stage <3-8 | build_geometry|trace_path_lengths|compute_prior|build_observation_matrix|invert_density|analyze_errors>]" << std::endl;
        std::exit(1);
    }
  }

  if (fname_json.empty()) {
    std::cerr << "Error: --json must be specified\n";
    std::exit(1);
  }
}

int main(int argc, char** argv)
{
  START_FUNC();

  parse_args(argc, argv);

  // Load JSON for logger setup
  const json js = myapp::load_json(fname_json);

  if (!js.contains("LOG_FILE")) {
    THROW_ERROR("JSON does not contain LOG_FILE section.");
  }

  // Logger setup
  const auto [stdout_level, stderr_level, file_level] = mylogger::get_log_levels(js, "LOG_FILE");
  std::string basename = std::string("sweep_") + fs::path(fname_json).stem().string();
  std::string log_dir_str = js.at("LOG_FILE").at("path_log_dir").get<std::string>();
  bool give_new_num = js.at("LOG_FILE").at("give_new_number").get<bool>();

  if (log_dir_str.empty()) log_dir_str = "logs";
  fs::path log_dir(log_dir_str);
  fs::create_directories(log_dir);
  auto log_filepath = mylogger::generate_unique_log_path(log_dir, basename, ".log", give_new_num);
  bool archive_existing = mylogger::get_log_archive_existing(js, "LOG_FILE");

  auto logger = mylogger::create_logger(stdout_level, stderr_level, file_level, log_filepath, archive_existing);
  spdlog::set_default_logger(logger);
  mylogger::configure_rate_limit(js, "LOG_FILE");

  LOG_INFO("Logging started to {}", log_filepath.string());

  // Seed setup: CLI > JSON > default
  if (seed_given) {
    seed::set_global_seed(seed_value);
    LOG_INFO("Seed from CLI: {}", seed_value);
  } else if (js.contains("seed")) {
    unsigned json_seed = js["seed"];
    seed::set_global_seed(json_seed);
    LOG_INFO("Seed from JSON: {}", json_seed);
  } else {
    seed::set_global_seed(42);
    LOG_WARN("No seed given. Using default seed 42.");
  }

  // Initialize parameter constants from JSON
  param_constants::init(js);

  // Execute sweep
  try {
    auto results = exemdl::sweep::run_sweep(fname_json, seed_given, seed_value, resume_dir, end_stage);
    LOG_INFO("Sweep completed: {} configurations executed", results.size());
    spdlog::default_logger()->flush();
  } catch (...) {
    spdlog::default_logger()->flush();
    LOG_CRITICAL("Unhandled exception in sweep");
    return 1;
  }

  END_FUNC();
  LOG_INFO("Seed was {}", seed::get_global_seed());
  spdlog::shutdown();
  return 0;
}
