/// @file ns_mymacro.hpp
/// @brief Common macro definitions and debug utilities
/// @details
/// This file provides the foundational macro and function infrastructure for error handling,
/// debugging, and diagnostic output across the muonith project. It includes:
/// - Modern error throwing functions using std::source_location and fmt::format
/// - Debug output macros (PRINT_MAT, PRINTF, PRINT_INT, etc.)
/// - Function entry/exit tracing (START_FUNC, END_FUNC)
/// - Legacy error handling functions (deprecated, kept for backward compatibility)
/// - Compile-time debug mode control (IF_DEBUG, PRAGMA)
///
/// Thread-safety: Most macros output to stderr/stdout without synchronization. Error throwing
/// functions are thread-safe for throwing exceptions, but logging output may interleave in
/// multithreaded contexts.
///
/// Dependencies: Requires spdlog, fmt, Eigen, nlohmann::json, and OpenBLAS/Accelerate.
#pragma once

#include <cstdio>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <chrono>
#include <thread>
#include <source_location>

#define _USE_MATH_DEFINES // Must be defined before including cmath to enable M_PI

// std::map<>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <cassert>
#include <vector>
#include <filesystem> // for std::filesystem::exists
#include <limits>
#include <tuple>
#include <stdexcept>  // for std::runtime_error
#include "spdlog_pch.hpp"
#include "ns_mylogger.hpp"  // spdlog logger
#include <Eigen/Dense>
#include <nlohmann/json.hpp>  // for nlohmann::json
#include <fmt/format.h>

using srcloc = std::source_location;
//==============================================================
/// @name macro_debug "the group of macros mainly for debug"
/// @{
//==============================================================

/// @brief print Eigen::Matrix, Vectors, so on
/// @details https://gist.github.com/AtsushiSakai/5175898
#define PRINT_MAT(X) std::cout << #X << ":\n" << X << "\n" << std::endl

/// @brief print float or double to stdout
#define PRINTF(X)    printf("DEBUG: file=%s, line=%d, func=%s : %s = %E\n",__FILE__,__LINE__,__PRETTY_FUNCTION__,#X,X)

/// @brief print float or double to stderr
#define PRINTF_E(X) fprintf(stderr,"DEBUG: file=%s, line=%d, func=%s : %s = %E\n",__FILE__,__LINE__,__PRETTY_FUNCTION__,#X,X)

/// @brief print int to stdout
#define PRINT_INT(X)   std::cout << "DEBUG: line=" << __LINE__ << ", func=" << __PRETTY_FUNCTION__ << " : " << #X << " = " << (X) << std::endl

/// @brief print int to stderr
#define PRINT_INT_E(X) std::cerr << "DEBUG: line=" << __LINE__ << ", func=" << __PRETTY_FUNCTION__ << " : " << #X << " = " << (X) << std::endl

/// @brief print matrix size to stderr
#define MAT_SIZE2(X) fprintf(stderr,"MAT_SIZE : %s.rows()=%ld, cols()=%ld\n",#X,X.rows(),X.cols())

/// @brief print variable name and its value to stderr
#define hogee(X)   std::cerr << "func=" << __PRETTY_FUNCTION__ << " : " << #X << " = " << X << std::endl

/// @brief print variable name and its value to std::cerr
#define hoge(X) { \
  std::ostringstream oss; \
  oss << X; \
  std::cerr << "hoge @ " << __PRETTY_FUNCTION__ << " : " << #X << " = " << oss.str() << std::endl; \
}


/// @brief print function name and start to std::cerr or mylogger
#define START_FUNC() \
  do { \
    if (mylogger::g_logger) { \
      LOG_DEBUG("{} start.", __PRETTY_FUNCTION__); \
    } else { \
      std::cerr << "[DEBUG] " << __PRETTY_FUNCTION__ << " start." << std::endl; \
    } \
  } while (0)

/// @brief print function name and end to std::cerr or mylogger
#define END_FUNC() \
  do { \
    if (mylogger::g_logger) { \
      LOG_DEBUG("{} end.", __PRETTY_FUNCTION__); \
    } else { \
      std::cerr << "[DEBUG] " << __PRETTY_FUNCTION__ << " end." << std::endl; \
    } \
  } while (0)

/// @brief sleep for X seconds
#define SLEEP_SEC(X)  fprintf(stderr,"sleep for %d seconds.\n",X); std::this_thread::sleep_for(std::chrono::seconds(X))

/// @brief sleep for X milliseconds
#define SLEEP_MSEC(X) fprintf(stderr,"sleep for %d milliseconds.\n",X); std::this_thread::sleep_for(std::chrono::milliseconds(X))

/// @brief get the value of current time
#define time_now std::chrono::system_clock::now()

// Historical THROW_ERROR implementations (kept for reference):
// Original attempt: throw std::runtime_error(__FUNCTION__ " : " X) - failed with gcc
// Alternative with sleep: caused incorrect behavior
// Simple version: throw std::runtime_error(X) - replaced by fmt-based version below

/// @brief helper function for THROW_ERROR
#define TOSTRING(x) #x
// #define CREATE_INSTANCE(T, x, ...) T x(__VA_ARGS__); x.set_name(TOSTRING(x))

/// @brief helper function for THROW_ERROR
#define SET_NAME(x) x.set_name(TOSTRING(x))

///@} ------------------------------------------------------------------

/// @brief Overload stream operator for tuple<int,int> to enable direct output
/// @details Enables direct output such as: std::cout << "it->second=" << it->second << std::endl;
inline std::ostream& operator<<(std::ostream& os, const std::tuple<int, int>& t) {
  os << "(" << std::get<0>(t) << "," << std::get<1>(t) << ")";
  return os;
}

/// @namespace mymacro
/// @brief Functions and macros primarily for debugging
namespace mymacro {

  /// @brief Configuration struct for debug macros
  struct DebugConfig {
    #ifdef NODEBUG
        static constexpr bool enabled = false;
    #else
        static constexpr bool enabled = true;
    #endif
  };  

  /// @brief Convert any streamable type to string
  /// @tparam T Type that supports operator<<(std::ostream&)
  /// @param[in] value The value to convert to string
  /// @return String representation of the value
  /// @note Used internally by deprecated THROW_ERROR2/3 macros
  template <typename T>
  std::string to_string_t(T value)
  {
      std::ostringstream oss;
      oss << value;
      return oss.str();
  };

  // ===== Modern error functions using source_location and fmt =====

  /// @brief Internal implementation: formats and throws error message
  /// @param[in] msg Error message to throw
  /// @param[in] loc Source location information (file, line, function)
  /// @throws std::runtime_error Always throws with formatted message including source location
  /// @note Logs to mylogger::g_logger at critical level if available
  [[noreturn]] inline void throw_error_impl(
    const std::string& msg,
    srcloc loc)
  {
    const std::string fullmsg = fmt::format(
      "{} at {}:{} - {}",
      loc.function_name(), loc.file_name(), loc.line(), msg);

    if (mylogger::g_logger) {
      mylogger::g_logger->critical("{}", fullmsg);
    }
    throw std::runtime_error(fullmsg);
  }

  /// @brief General error throwing function (runtime string version)
  /// @param[in] msg Error message string
  /// @param[in] loc Source location (defaults to call site)
  /// @throws std::runtime_error Always throws with formatted message
  /// @note Prefer using THROW_ERROR() macro instead of calling directly
  [[noreturn]] inline void throw_error(
    const std::string& msg,
    srcloc loc = srcloc::current())
  {
    throw_error_impl(msg, loc);
  }

  /// @brief General error throwing function (fmt::format compatible, compile-time string version)
  /// @tparam Args Variadic template parameters for format arguments
  /// @param[in] fmt_str Format string (compile-time checked)
  /// @param[in] args Arguments to format into the string
  /// @throws std::runtime_error Always throws with formatted message
  /// @note Prefer using THROW_ERROR("format {}", value) macro. Provides compile-time format string validation.
  template<typename... Args>
    requires (sizeof...(Args) > 0)
  [[noreturn]] inline void throw_error(
    fmt::format_string<Args...> fmt_str,
    Args&&... args)
  {
    const std::string msg = fmt::format(fmt_str, std::forward<Args>(args)...);
    throw_error_impl(msg, srcloc::current());
  }

  /// @brief Error throwing with instance name (runtime string version)
  /// @tparam T Type that has get_name() method returning std::string
  /// @param[in] instance Pointer to instance (must have get_name() method)
  /// @param[in] msg Error message string
  /// @param[in] loc Source location (defaults to call site)
  /// @throws std::runtime_error Always throws with formatted message including instance name
  /// @note Prefer using THROW_ERROR_NAME() macro from member functions. Requires instance->get_name() to be valid.
  template<typename T>
  [[noreturn]] inline void throw_error_with_name(
    const T* instance,
    const std::string& msg,
    srcloc loc = srcloc::current())
  {
    const std::string fullmsg = fmt::format(
      "{} at {}:{} - Instance: {} - {}",
      loc.function_name(), loc.file_name(), loc.line(),
      instance->get_name(), msg);

    if (mylogger::g_logger) {
      mylogger::g_logger->critical("{}", fullmsg);
    }
    throw std::runtime_error(fullmsg);
  }

  /// @brief Error throwing with instance name (fmt::format compatible, compile-time string version)
  /// @tparam T Type that has get_name() method returning std::string
  /// @tparam Args Variadic template parameters for format arguments
  /// @param[in] instance Pointer to instance (must have get_name() method)
  /// @param[in] fmt_str Format string (compile-time checked)
  /// @param[in] args Arguments to format into the string
  /// @throws std::runtime_error Always throws with formatted message including instance name
  /// @note Prefer using THROW_ERROR_NAME("format {}", value) macro. Provides compile-time format string validation.
  template<typename T, typename... Args>
    requires (sizeof...(Args) > 0)
  [[noreturn]] inline void throw_error_with_name(
    const T* instance,
    fmt::format_string<Args...> fmt_str,
    Args&&... args)
  {
    const std::string msg = fmt::format(fmt_str, std::forward<Args>(args)...);
    const auto loc = srcloc::current();
    const std::string fullmsg = fmt::format(
      "{} at {}:{} - Instance: {} - {}",
      loc.function_name(), loc.file_name(), loc.line(),
      instance->get_name(), msg);

    if (mylogger::g_logger) {
      mylogger::g_logger->critical("{}", fullmsg);
    }
    throw std::runtime_error(fullmsg);
  }

  // ===== Legacy functions (kept for compatibility, deprecated) =====

  /// @brief throw error (deprecated)
  void throw_error(const std::string& msg
  , const std::string& var_name, const std::string& var_value
  , const std::string& function, const std::string& file, const int line);

  /// @brief throw error (deprecated)
  void throw_error(const std::string& msg
  , const std::string& var1_name, const std::string& var1_value
  , const std::string& var2_name, const std::string& var2_value
  , const std::string& function, const std::string& file, const int line);

  /// @brief throw_error with instance name (deprecated)
  void throw_error_name(const std::string& msg
  , const std::string& var_name, const std::string& var_value
  , const std::string& function, const std::string& file
  , const int line, const std::string &instance_name);

  /// @brief throw_error with instance name (deprecated)
  void throw_error_name(const std::string& msg
  , const std::string& var1_name, const std::string& var1_value
  , const std::string& var2_name, const std::string& var2_value
  , const std::string& function, const std::string& file
  , const int line, const std::string &instance_name);

  /// @brief Display matrix dimensions via logger
  /// @param[in] mat_name Name of the matrix (for display)
  /// @param[in] mat Eigen matrix to display dimensions for
  /// @param[in] log_level spdlog log level (debug, info, warn, etc.)
  /// @note Logs to mylogger::g_logger. Outputs "mat_name=<name> rows=<R> cols=<C>"
  void disp_mat_size(
          const std::string &mat_name
        , const Eigen::MatrixXf &mat
        , const spdlog::level::level_enum log_level);
};

/// @brief helper function for THROW_ERROR
/// @ingroup macro_debug
#define STRINGIFY_HELPER(X) #X
/// @brief helper function for THROW_ERROR
/// @ingroup macro_debug
#define STRINGIFY(X) STRINGIFY_HELPER(X)

// ===== Modern macros using source_location and fmt =====

/// @brief throw error with fmt::format syntax and source_location
/// @ingroup macro_debug
/// @details Usage example: THROW_ERROR("Invalid value: {}", value);
#define THROW_ERROR(...) mymacro::throw_error(__VA_ARGS__)

/// @brief throw error with instance name, fmt::format syntax and source_location
/// @ingroup macro_debug
/// @details Usage example: THROW_ERROR_NAME("Invalid value: {}", value);
#define THROW_ERROR_NAME(...) mymacro::throw_error_with_name(this, __VA_ARGS__)

// ===== Legacy macros (kept for compatibility, deprecated) =====

/// @brief throw error (deprecated - use THROW_ERROR with fmt syntax instead)
/// @ingroup macro_debug
#define THROW_ERROR2(msg, var) mymacro::throw_error(msg, #var, mymacro::to_string_t(var), __PRETTY_FUNCTION__, __FILE__, __LINE__)

/// @brief throw error (deprecated - use THROW_ERROR with fmt syntax instead)
/// @ingroup macro_debug
#define THROW_ERROR3(msg, var1, var2) mymacro::throw_error(msg, #var1, mymacro::to_string_t(var1), #var2, mymacro::to_string_t(var2), __PRETTY_FUNCTION__, __FILE__, __LINE__)

/// @brief throw error (deprecated - use THROW_ERROR_NAME with fmt syntax instead)
/// @ingroup macro_debug
#define THROW_ERROR_NAME2(msg, var) mymacro::throw_error_name(msg, #var, mymacro::to_string_t(var), __PRETTY_FUNCTION__, __FILE__, __LINE__, this->get_name() )

/// @brief throw error (deprecated - use THROW_ERROR_NAME with fmt syntax instead)
/// @ingroup macro_debug
#define THROW_ERROR_NAME3(msg, var1, var2) mymacro::throw_error_name(msg, #var1, mymacro::to_string_t(var1), #var2, mymacro::to_string_t(var2), __PRETTY_FUNCTION__, __FILE__, __LINE__, this->get_name() )

/// @brief display matrix size
#define MAT_SIZE(X) mymacro::disp_mat_size(#X, X, spdlog::level::debug)

/// @brief Display a title message with decorative borders
/// @param[in] msg The message to display
/// @details Outputs a title message surrounded by decorative borders (######).
///          - If mylogger::g_logger is initialized: outputs via logger at INFO level
///          - If mylogger::g_logger is nullptr: falls back to std::cerr
///          This ensures the function never crashes even when the logger is uninitialized.
/// @note Thread-safety: Uses logger if available (thread-safe), otherwise stderr (not synchronized)
inline void infotitle(const std::string& msg) {
  if (mylogger::g_logger) {
    mylogger::g_logger->info("######################################################");
    mylogger::g_logger->info("{}", msg);
    mylogger::g_logger->info("######################################################");
  } else {
    std::cerr << "######################################################\n";
    std::cerr << msg << "\n";
    std::cerr << "######################################################" << std::endl;
  }
}

/// @brief Convert boolean to string literal "true" or "false"
/// @param[in] b Boolean value to convert
/// @return Pointer to static string literal "true" or "false"
/// @note Returns pointer to static storage, safe to use without copying
inline const char* tf(const bool b) {
  return b ? "true" : "false";
};

/// @brief Macro to manage DEBUG mode in a single line
#define IF_DEBUG if constexpr (mymacro::DebugConfig::enabled)

/// @brief Macro to toggle OPENMP pragmas between DEBUG and Release builds
#ifdef NODEBUG
  #define PRAGMA(options)
#else
  #define PRAGMA(options) _Pragma(#options)
#endif
