/// @file ns_iodir.hpp
/// @brief Directory and file path utilities
/// @details
/// This module provides utilities for output directory management and path manipulation.
/// Key features:
/// - Global default output directory management (default: ./tmp)
/// - Path combining with automatic directory creation
/// - Multiple overloads for const char*, std::string, and fs::path inputs
///
/// Workflow:
/// 1. Optionally set a default output directory using set_default_output_dir()
/// 2. Use make_pathout() to combine file paths with output directories
/// 3. If the target directory does not exist, it is created automatically
/// 4. Absolute paths are returned unchanged (bypass directory combining)
///
/// Thread-safety: No. The global default output directory is not thread-safe.
/// Concurrent calls to set_default_output_dir() and make_pathout() may race.
///
/// @note Dependencies: Requires C++17 std::filesystem
#pragma once
#include <filesystem>
namespace fs = std::filesystem;

/// @brief Namespace for output directory management
/// @details The global default output directory is defined as a static variable in an anonymous namespace.
namespace iodir {
  
  /// @brief Set the default output directory
  /// @param[in] dir New default output directory. If empty, resets to "./tmp".
  /// @note Thread-safety: No. Not safe for concurrent calls.
  void set_default_output_dir(const fs::path &dir);
  
  /// @brief Get the current default output directory
  /// @return Reference to the default output directory (default: "./tmp")
  /// @note Thread-safety: No. Concurrent modifications via set_default_output_dir() may cause race conditions.
  const fs::path &get_default_output_dir();

  /// @brief Create a path by combining file_path with the specified output directory
  /// @param[in] file_path Input file path (relative or absolute)
  /// @param[in] output_dir Output directory to combine with
  /// @return Combined path (output_dir / file_path), or file_path unchanged if absolute
  /// @throws std::runtime_error If directory creation fails
  /// @note Creates intermediate directories automatically if they do not exist.
  ///       Absolute paths bypass combination and are returned unchanged.
  fs::path make_pathout( const fs::path &file_path, const fs::path& output_dir);

  /// @brief Create a path by combining cfname with the specified output directory
  /// @param[in] cfname C-style string file name (relative or absolute)
  /// @param[in] output_dir Output directory to combine with
  /// @return Combined path (output_dir / cfname), or cfname unchanged if absolute
  /// @throws std::runtime_error If directory creation fails
  /// @note Thin wrapper over make_pathout(fs::path, fs::path).
  fs::path make_pathout( const char* cfname, const fs::path& output_dir);

  /// @brief Create a path by combining fname_str with the specified output directory
  /// @param[in] fname_str String file name (relative or absolute)
  /// @param[in] output_dir Output directory to combine with
  /// @return Combined path (output_dir / fname_str), or fname_str unchanged if absolute
  /// @throws std::runtime_error If directory creation fails
  /// @note Thin wrapper over make_pathout(fs::path, fs::path).
  fs::path make_pathout( const std::string& fname_str, const fs::path& output_dir);

  /// @brief Create a path using the globally-set default output directory (./tmp)
  /// @param[in] file_path Input file path (relative or absolute)
  /// @return Combined path using default output directory, or file_path unchanged if absolute
  /// @throws std::runtime_error If directory creation fails
  /// @note Uses the default output directory set via set_default_output_dir() (default: "./tmp").
  fs::path make_pathout( const fs::path& file_path);

  /// @brief Create a path using the globally-set default output directory (./tmp)
  /// @param[in] cfname C-style string file name (relative or absolute)
  /// @return Combined path using default output directory, or cfname unchanged if absolute
  /// @throws std::runtime_error If directory creation fails
  /// @note Thin wrapper over make_pathout(fs::path).
  fs::path make_pathout( const char* cfname );

  /// @brief Create a path using the globally-set default output directory (./tmp)
  /// @param[in] fname_str String file name (relative or absolute)
  /// @return Combined path using default output directory, or fname_str unchanged if absolute
  /// @throws std::runtime_error If directory creation fails
  /// @note Thin wrapper over make_pathout(fs::path).
  fs::path make_pathout( const std::string& fname_str );

}; // namespace iodir