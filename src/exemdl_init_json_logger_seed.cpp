// src/exemdl_init_json_logger_seed.cpp
#include "exemdl_init_json_logger_seed.hpp"
#include "ns_myapp.hpp"
#include "ns_seed.hpp"
#include "ns_mylogger.hpp"
#include <spdlog/spdlog.h>
#include <thread>
#include <omp.h>
#include <Eigen/Core>

nlohmann::json exemdl::init_json_logger_seed::load_and_validate_json(
  const std::filesystem::path& path_json)
{
  // Load JSON
  auto js = myapp::load_json(path_json.string());
  // Check required sections
  if (!js.contains("LOG_FILE")) {
    THROW_ERROR("load_and_validate_json: LOG_FILE section missing. path={}", path_json.string());
  }
  return js;
}

void exemdl::init_json_logger_seed::configure_logging_and_seed(
    const nlohmann::json& js, const std::filesystem::path& path_json
  , const bool seed_given , const unsigned seed_value )
{
  // Configure logging
  const auto [stdout_level, stderr_level, file_level]
   = mylogger::get_log_levels(js, "LOG_FILE");

  const std::string basename =
    "main_" + path_json.stem().string();

  std::string log_dir_str =
    js.at("LOG_FILE").at("path_log_dir").get<std::string>();

  if (log_dir_str.empty()) log_dir_str = "logs";

  // Create log directory
  std::filesystem::path log_dir(log_dir_str);
  std::filesystem::create_directories(log_dir);

  // Generate log file path
  const bool give_new_num =
    js.at("LOG_FILE").at("give_new_number").get<bool>();

  const auto log_filepath
    = mylogger::generate_unique_log_path( log_dir, basename, ".log", give_new_num );

  // Create logger
  auto logger = mylogger::create_logger(
    stdout_level, stderr_level, file_level, log_filepath );

  // Set default logger
  spdlog::set_default_logger(logger);

  LOG_INFO("Logging started to {}", log_filepath.string());

  // Configure seed
  if (seed_given) {
    seed::set_global_seed(seed_value);
    LOG_INFO("Seed from CLI: {}", seed_value);
  }
  else if (js.contains("seed")) {
    unsigned json_seed = js["seed"].get<unsigned>();
    seed::set_global_seed(json_seed);
    LOG_INFO("Seed from JSON: {}", json_seed);
  } else {
    seed::set_global_seed(42);
    LOG_WARN("No seed provided. Using default seed: 42");
  }
}

nlohmann::json exemdl::init_json_logger_seed::init_all(
  const InitArgs& args)
{
  // 1) Load JSON and validate required keys
  nlohmann::json js = load_and_validate_json(args.path_json);

  // 2) Configure logging and seed
  configure_logging_and_seed(
      js
    , args.path_json
    , args.seed_given
    , args.seed_value
  );

  // 3) Check thread count and configure Eigen
  unsigned int n_max = std::thread::hardware_concurrency();
  LOG_INFO("Maximum threads on this PC: {}", n_max);

  int n_omp = omp_get_max_threads();
  LOG_INFO("Maximum threads for this program: {}", n_omp);
  SLEEP_MSEC(500);

  int n_eigen = n_omp / 2;
  LOG_INFO("Eigen threads: {}", n_eigen);
  Eigen::setNbThreads(n_eigen);

  return js;
}

