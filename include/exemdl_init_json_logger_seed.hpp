/// @file exemdl_init_json_logger_seed.hpp
/// @brief JSON configuration, logger and random seed initialization
/// @details
/// This module provides the initialization workflow for applications that require:
/// - JSON configuration file loading and validation
/// - Logger setup with configurable output levels (stdout, stderr, file)
/// - Random seed initialization from CLI arguments or JSON configuration
/// - Thread configuration for OpenMP and Eigen
///
/// Typical usage:
/// @code
/// exemdl::init_json_logger_seed::InitArgs args{path_json, seed_given, seed_value};
/// nlohmann::json config = exemdl::init_json_logger_seed::init_all(args);
/// @endcode
///
/// @note Thread-safety: This module uses OpenMP (omp_get_max_threads) and Eigen thread
/// configuration. Call init_all() once during application startup before spawning threads.
#pragma once

#include <filesystem>
#include <nlohmann/json.hpp>

/// @namespace exemdl::init_json_logger_seed
/// @brief Namespace for initializing JSON logger and seed
/// @details
/// This namespace encapsulates the initialization sequence required at application startup:
/// JSON loading, logger configuration, seed initialization, and thread setup.
/// @ingroup ExecModule
namespace exemdl::init_json_logger_seed {

  /// @brief Arguments required for module initialization
  /// @details
  /// This structure aggregates all inputs needed for the init_all() workflow.
  /// Typically constructed from command-line arguments parsed by the main function.
  struct InitArgs {
    const std::filesystem::path& path_json; ///< Path to JSON configuration file
    const bool                   seed_given; ///< Whether seed is explicitly provided via CLI
    const unsigned               seed_value; ///< Seed value (only used if seed_given is true)
  };

  /// @brief Load JSON file and validate existence of required sections
  /// @param[in] path_json Path to JSON configuration file
  /// @return Loaded JSON object
  /// @throws std::runtime_error If JSON file cannot be loaded
  /// @throws std::runtime_error If required section "LOG_FILE" is missing
  /// @note This function verifies that the LOG_FILE section exists in the JSON configuration.
  /// Other sections may be validated by the caller as needed.
  nlohmann::json load_and_validate_json(
    const std::filesystem::path& path_json );

  /// @brief Initialize seed and logging output based on CLI and JSON configuration
  /// @param[in] js Loaded JSON object containing LOG_FILE section and optional seed field
  /// @param[in] path_json Path to JSON file (used for generating log name prefix)
  /// @param[in] seed_given Whether seed was specified via CLI
  /// @param[in] seed_value Seed value specified via CLI
  /// @details
  /// This function performs the following steps:
  /// 1. Extracts log levels (stdout, stderr, file) from js["LOG_FILE"]
  /// 2. Creates log directory from js["LOG_FILE"]["path_log_dir"] (defaults to "logs")
  /// 3. Generates unique log file path with optional numbering
  /// 4. Creates and sets the default spdlog logger
  /// 5. Initializes global random seed (priority: CLI > JSON > default=42)
  ///
  /// @note Seed priority: CLI argument takes precedence over JSON configuration.
  /// If neither is provided, uses default seed value 42.
  /// @note Thread-safety: Must be called before spawning threads. Logger creation is not thread-safe.
  void configure_logging_and_seed(
      const nlohmann::json& js
    , const std::filesystem::path& path_json
    , const bool seed_given
    , const unsigned seed_value );

  /// @brief Perform JSON loading, logging/seed configuration, and thread initialization, then return JSON
  /// @param[in] args Input arguments containing JSON path, seed configuration
  /// @return Loaded JSON object for further application use
  /// @throws std::runtime_error If JSON loading or validation fails
  /// @details
  /// This is the main entry point for application initialization. It orchestrates:
  /// 1. JSON loading and validation (via load_and_validate_json)
  /// 2. Logger and seed configuration (via configure_logging_and_seed)
  /// 3. Thread detection and Eigen thread configuration
  ///
  /// Thread configuration:
  /// - Detects hardware concurrency (std::thread::hardware_concurrency)
  /// - Queries OpenMP max threads (omp_get_max_threads)
  /// - Sets Eigen threads to half of OpenMP threads (n_omp / 2)
  ///
  /// @note Thread-safety: Must be called once at application startup before spawning threads.
  /// Uses OpenMP and Eigen thread configuration which are not thread-safe during initialization.
  nlohmann::json init_all(const InitArgs& args);

}; // namespace exemdl::init_json_logger_seed
