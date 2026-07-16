/// @file ns_mylogger.hpp
/// @brief spdlog-based logging system for the application
/// @details
/// This module provides a multi-sink logger based on spdlog, supporting
/// simultaneous output to stdout, stderr, and file with independent log levels.
///
/// ## Basic Usage
/// 1. Initialize the global logger (g_logger) using create_logger()
/// 2. Use logging macros (LOG_INFO(), LOG_ERROR(), etc.) to output logs
///
/// ## Behavior When Logger is Uninitialized
/// - If g_logger is nullptr, all logging functions automatically fall back to spdlog's default logger
/// - The program will NOT crash and continues to operate safely
/// - safe_log() outputs "[fallback]" prefix when falling back
///
/// ## Thread Safety
/// - Thread-safe: Yes. The logger uses multi-threaded sinks (_mt variants) with internal synchronization.
/// - All logging functions (LOG_*) are safe to call concurrently from multiple threads.
/// - Logger initialization (create_logger) should be called only once, preferably before spawning threads.
///
/// ## Example Usage
/// @code
/// // Initialize logger
/// auto logger = mylogger::create_logger(
///     spdlog::level::info,  // stdout level
///     spdlog::level::warn,  // stderr level
///     spdlog::level::trace, // file level
///     "logs/app.log");
///
/// // Log output
/// LOG_INFO("Application started");
/// LOG_ERROR("Error code: {}", error_code);
/// LOG_DEBUG_TAG("NETWORK", "Connected to {}", host);
/// @endcode
///
/// ## Available Macros
/// - LOG_TRACE(...): Trace level (with source location)
/// - LOG_INFO(...): Info level (with source location)
/// - LOG_DEBUG(...): Debug level (with source location)
/// - LOG_WARN(...): Warning level (with source location)
/// - LOG_ERROR(...): Error level (with source location)
/// - LOG_CRITICAL(...): Critical level (with source location)
/// - LOG_*_TAG(tag, ...): Tagged versions of the above macros
#pragma once

#include <atomic>
#include <cstdio>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <filesystem> // for std::filesystem::exists
#include <chrono>
#include <thread>
#include <nlohmann/json.hpp>  // for json
#include <fmt/format.h>
#include <fmt/std.h>   // Enable fmt support for std::filesystem::path
#include <fmt/ranges.h>  // For join and parse_context
#include "spdlog_pch.hpp"

/// @namespace mylogger
/// @brief Logging-related functions and utilities
namespace mylogger {

  /// @brief Global logger pointer
  /// @details This shared pointer holds the application's primary logger instance.
  ///          Initialize it using create_logger() before logging.
  ///          If nullptr, logging functions fall back to spdlog's default logger.
  extern std::shared_ptr<class spdlog::logger> g_logger;

  /// @brief Generate a unique log file path
  /// @param[in] log_dir Directory for log files
  /// @param[in] base_name Base name of the log file
  /// @param[in] extension File extension (including the dot)
  /// @param[in] give_new_number If true, append a sequential number if the file exists
  /// @return Unique log file path
  /// @details Creates the directory if it doesn't exist. If give_new_number is true
  ///          and the file already exists, appends "_0000", "_0001", etc. until a
  ///          unique filename is found.
  std::filesystem::path
    generate_unique_log_path(
        const std::filesystem::path log_dir = "logs",
        const std::string& base_name = "log",
        const std::string& extension = ".log",
        const bool give_new_number = false);

  /// @brief Create a multi-sink logger (stdout, stderr, file)
  /// @param[in] stdout_level Log level for standard output
  /// @param[in] stderr_level Log level for standard error output
  /// @param[in] file_level   Log level for file output
  /// @param[in] log_file_path Path to the log file
  /// @param[in] archive_existing If true, rename existing file to {stem}_{0001}.{ext} and create new log.
  ///            If false (default), overwrite existing log file.
  /// @return Shared pointer to the created logger
  /// @note This function sets the created logger as g_logger (globally accessible).
  ///       It also registers the logger as spdlog's default logger.
  ///       Thread-safety: Should be called only once during initialization, before multi-threaded execution.
  /// @throws std::exception If the file path is invalid or cannot be opened
  /// @details The logger uses the pattern: "[%Y-%m-%d %H:%M:%S] [%^%l%$] [%s:%# %!] %v"
  ///          where %s=file, %#=line, %!=function, %v=message.
  ///          Uses multi-threaded sinks (_mt) for thread-safe concurrent logging.
  std::shared_ptr<spdlog::logger>
    create_logger(
        spdlog::level::level_enum stdout_level,
        spdlog::level::level_enum stderr_level,
        spdlog::level::level_enum file_level,
        const std::filesystem::path& log_file_path,
        bool archive_existing = false);

  /// @brief Parse log levels from JSON configuration
  /// @param[in] config JSON configuration object
  /// @param[in] section_name Target section name in the JSON
  /// @return Tuple of (stdout_level, stderr_level, file_level)
  /// @note Expected JSON keys: "stdout_level", "stderr_level", "file_level"
  ///       Valid values: "off", "trace", "debug", "info", "warn", "warning", "err", "error", "critical"
  /// @details Default values: stdout="info", stderr="warn", file="trace"
  std::tuple<spdlog::level::level_enum, spdlog::level::level_enum, spdlog::level::level_enum>
    get_log_levels(const nlohmann::json &config, const std::string &section_name);

  /// @brief Parse archive_existing setting from JSON configuration
  /// @param[in] config JSON configuration object
  /// @param[in] section_name Target section name in the JSON
  /// @return true to archive existing log file with numbered suffix, false to overwrite
  /// @details Default: false (overwrite). JSON key: "archive_existing"
  bool get_log_archive_existing(const nlohmann::json& config, const std::string& section_name);

  /// @brief Safely output a log message (without source location)
  /// @param[in] level Log level (spdlog::level::level_enum)
  /// @param[in] msg Message string to output
  /// @details If g_logger is initialized, outputs via g_logger.
  ///          If g_logger is nullptr, falls back to spdlog::log() with "[fallback]" prefix.
  ///          This function never crashes even if the logger is uninitialized.
  void safe_log(spdlog::level::level_enum level, const std::string &msg);

  /// @brief Check if logging should proceed based on rate limit counter
  /// @param[in,out] counter Atomic counter for tracking call count (static per call site)
  /// @param[in] max_count Maximum number of times to allow logging
  /// @return true if logging should proceed, false if rate limit exceeded
  /// @details Thread-safe. Uses optimized double-check pattern:
  ///          first load() for fast path, then fetch_add() for actual counting.
  /// @note This function is called by LOG_*_N macros. Not intended for direct use.
  bool should_log_n(std::atomic<int>& counter, int max_count);

  // ========= Rate-limit default configuration =========

  /// @brief Default max count for rate-limited logging (configurable via JSON)
  /// @details Initial value is 100. Can be changed via set_default_max_count() or configure_rate_limit().
  extern int g_default_max_count;

  /// @brief Get current default max count for rate-limited logging
  /// @return Current default max count value
  [[nodiscard]] int get_default_max_count() noexcept;

  /// @brief Set default max count for rate-limited logging
  /// @param[in] value New default max count value
  /// @note Typically called during initialization
  void set_default_max_count(int value) noexcept;

  /// @brief Parse rate_limit config from JSON and set defaults
  /// @param[in] config JSON configuration object
  /// @param[in] section_name Target section name in the JSON (e.g., "LOG_FILE")
  /// @details Expected JSON structure:
  ///          { "section_name": { "rate_limit": { "default_max_count": 10 } } }
  ///          If the section or rate_limit key is missing, defaults are unchanged.
  void configure_rate_limit(const nlohmann::json& config, const std::string& section_name);

  // ========= Internal implementation functions (non-template, defined in .cpp) =========

  /// @brief TRACE level logging with source location
  /// @param[in] file Source file name (__FILE__)
  /// @param[in] line Line number (__LINE__)
  /// @param[in] func Function name (__PRETTY_FUNCTION__)
  /// @param[in] tf_show_detail If true, include source location in output
  /// @param[in] msg Message to log
  /// @details Falls back to spdlog::log() if g_logger is nullptr. Never crashes.
  void log_trace_msg(const char* file, int line, const char* func, bool tf_show_detail, const std::string& msg);

  /// @brief INFO level logging with source location
  /// @param[in] file Source file name (__FILE__)
  /// @param[in] line Line number (__LINE__)
  /// @param[in] func Function name (__PRETTY_FUNCTION__)
  /// @param[in] tf_show_detail If true, include source location in output
  /// @param[in] msg Message to log
  /// @details Falls back to spdlog::log() if g_logger is nullptr. Never crashes.
  void log_info_msg (const char* file, int line, const char* func, bool tf_show_detail, const std::string& msg);

  /// @brief INFO level logging with tag and source location
  /// @param[in] file Source file name (__FILE__)
  /// @param[in] line Line number (__LINE__)
  /// @param[in] func Function name (__PRETTY_FUNCTION__)
  /// @param[in] tf_show_detail If true, include source location in output
  /// @param[in] tag Tag prefix for the message (e.g., "NETWORK", "DATABASE")
  /// @param[in] msg Message to log
  /// @details Falls back to spdlog::log() if g_logger is nullptr. Never crashes.
  void log_info_tag_msg(const char* file, int line, const char* func, bool tf_show_detail, const std::string& tag, const std::string& msg);

  /// @brief DEBUG level logging with source location
  /// @param[in] file Source file name (__FILE__)
  /// @param[in] line Line number (__LINE__)
  /// @param[in] func Function name (__PRETTY_FUNCTION__)
  /// @param[in] tf_show_detail If true, include source location in output
  /// @param[in] msg Message to log
  /// @details Falls back to spdlog::log() if g_logger is nullptr. Never crashes.
  void log_debug_loc_msg   (const char* file, int line, const char* func, bool tf_show_detail, const std::string& msg);

  /// @brief WARN level logging with source location
  /// @param[in] file Source file name (__FILE__)
  /// @param[in] line Line number (__LINE__)
  /// @param[in] func Function name (__PRETTY_FUNCTION__)
  /// @param[in] tf_show_detail If true, include source location in output
  /// @param[in] msg Message to log
  /// @details Falls back to spdlog::log() if g_logger is nullptr. Never crashes.
  void log_warn_loc_msg    (const char* file, int line, const char* func, bool tf_show_detail, const std::string& msg);

  /// @brief ERROR level logging with source location
  /// @param[in] file Source file name (__FILE__)
  /// @param[in] line Line number (__LINE__)
  /// @param[in] func Function name (__PRETTY_FUNCTION__)
  /// @param[in] tf_show_detail If true, include source location in output
  /// @param[in] msg Message to log
  /// @details Falls back to spdlog::log() if g_logger is nullptr. Never crashes.
  void log_error_loc_msg   (const char* file, int line, const char* func, bool tf_show_detail, const std::string& msg);

  /// @brief CRITICAL level logging with source location
  /// @param[in] file Source file name (__FILE__)
  /// @param[in] line Line number (__LINE__)
  /// @param[in] func Function name (__PRETTY_FUNCTION__)
  /// @param[in] tf_show_detail If true, include source location in output
  /// @param[in] msg Message to log
  /// @details Falls back to spdlog::log() if g_logger is nullptr. Never crashes.
  void log_critical_loc_msg(const char* file, int line, const char* func, bool tf_show_detail, const std::string& msg);

  /// @brief DEBUG level logging with tag and source location
  /// @param[in] file Source file name (__FILE__)
  /// @param[in] line Line number (__LINE__)
  /// @param[in] func Function name (__PRETTY_FUNCTION__)
  /// @param[in] tf_show_detail If true, include source location in output
  /// @param[in] tag Tag prefix for the message
  /// @param[in] msg Message to log
  /// @details Falls back to spdlog::log() if g_logger is nullptr. Never crashes.
  void log_debug_tag_loc_msg   (const char* file, int line, const char* func, bool tf_show_detail, const std::string& tag, const std::string& msg);

  /// @brief WARN level logging with tag and source location
  /// @param[in] file Source file name (__FILE__)
  /// @param[in] line Line number (__LINE__)
  /// @param[in] func Function name (__PRETTY_FUNCTION__)
  /// @param[in] tf_show_detail If true, include source location in output
  /// @param[in] tag Tag prefix for the message
  /// @param[in] msg Message to log
  /// @details Falls back to spdlog::log() if g_logger is nullptr. Never crashes.
  void log_warn_tag_loc_msg    (const char* file, int line, const char* func, bool tf_show_detail, const std::string& tag, const std::string& msg);

  /// @brief ERROR level logging with tag and source location
  /// @param[in] file Source file name (__FILE__)
  /// @param[in] line Line number (__LINE__)
  /// @param[in] func Function name (__PRETTY_FUNCTION__)
  /// @param[in] tf_show_detail If true, include source location in output
  /// @param[in] tag Tag prefix for the message
  /// @param[in] msg Message to log
  /// @details Falls back to spdlog::log() if g_logger is nullptr. Never crashes.
  void log_error_tag_loc_msg   (const char* file, int line, const char* func, bool tf_show_detail, const std::string& tag, const std::string& msg);

  /// @brief CRITICAL level logging with tag and source location
  /// @param[in] file Source file name (__FILE__)
  /// @param[in] line Line number (__LINE__)
  /// @param[in] func Function name (__PRETTY_FUNCTION__)
  /// @param[in] tf_show_detail If true, include source location in output
  /// @param[in] tag Tag prefix for the message
  /// @param[in] msg Message to log
  /// @details Falls back to spdlog::log() if g_logger is nullptr. Never crashes.
  void log_critical_tag_loc_msg(const char* file, int line, const char* func, bool tf_show_detail, const std::string& tag, const std::string& msg);

  /// @brief Extract basename from file path at compile time
  /// @param[in] path Full file path (typically __FILE__)
  /// @return Pointer to the basename portion of the path
  /// @details No memory allocation. Constexpr for compile-time evaluation.
  static inline constexpr const char* basename_cstr(const char* path) noexcept {
    const char* last = path;
    for (const char* p = path; *p; ++p) {
      if (*p == '/' || *p == '\\') last = p + 1;
    }
    return last;
  }

  /// @brief Extract ClassName::method from __PRETTY_FUNCTION__ or __PRETTY_FUNCTION__
  /// @param[in] pretty_func Full function signature (e.g., "void ClassName::method(int)")
  /// @return Extracted "ClassName::method" string, or the original if extraction fails
  /// @details Handles GCC/Clang __PRETTY_FUNCTION__ format.
  ///          For free functions, returns just the function name.
  std::string extract_class_method(const char* pretty_func);

  // ========= Lightweight inline wrapper layer for header size reduction =========
  /// @namespace mylogger::detail
  /// @brief Internal implementation details. Do not call these functions directly.
  ///        Use the LOG_* macros instead.
  namespace detail {

    // --- TRACE/INFO ---
    inline void info_here(const char* f,int l,const char* fn,bool tf_detail,const std::string& msg){
      mylogger::log_info_msg(f,l,fn,tf_detail,msg);
    }
    template<class... Args>
    inline void info_here(const char* f,int l,const char* fn,bool tf_detail, fmt::format_string<Args...> s, Args&&... args){
      mylogger::log_info_msg(f,l,fn,tf_detail, fmt::format(s, std::forward<Args>(args)...));
    }

    inline void trace_here(const char* f,int l,const char* fn,bool tf_detail,const std::string& msg){
      mylogger::log_trace_msg(f,l,fn,tf_detail,msg);
    }
    template<class... Args>
    inline void trace_here(const char* f,int l,const char* fn,bool tf_detail, fmt::format_string<Args...> s, Args&&... args){
      mylogger::log_trace_msg(f,l,fn,tf_detail, fmt::format(s, std::forward<Args>(args)...));
    }

    // --- DEBUG/WARN/ERROR/CRITICAL ---
    inline void debug_here(const char* f,int l,const char* fn,bool tf_detail,const std::string& msg){
      mylogger::log_debug_loc_msg(f,l,fn,tf_detail,msg);
    }
    template<class... Args>
    inline void debug_here(const char* f,int l,const char* fn,bool tf_detail, fmt::format_string<Args...> s, Args&&... args){
      mylogger::log_debug_loc_msg(f,l,fn,tf_detail, fmt::format(s, std::forward<Args>(args)...));
    }

    inline void warn_here(const char* f,int l,const char* fn,bool tf_detail,const std::string& msg){
      mylogger::log_warn_loc_msg(f,l,fn,tf_detail,msg);
    }
    template<class... Args>
    inline void warn_here(const char* f,int l,const char* fn,bool tf_detail, fmt::format_string<Args...> s, Args&&... args){
      mylogger::log_warn_loc_msg(f,l,fn,tf_detail, fmt::format(s, std::forward<Args>(args)...));
    }

    inline void error_here(const char* f,int l,const char* fn,bool tf_detail,const std::string& msg){
      mylogger::log_error_loc_msg(f,l,fn,tf_detail,msg);
    }
    template<class... Args>
    inline void error_here(const char* f,int l,const char* fn,bool tf_detail, fmt::format_string<Args...> s, Args&&... args){
      mylogger::log_error_loc_msg(f,l,fn,tf_detail, fmt::format(s, std::forward<Args>(args)...));
    }

    inline void critical_here(const char* f,int l,const char* fn,bool tf_detail,const std::string& msg){
      mylogger::log_critical_loc_msg(f,l,fn,tf_detail,msg);
    }
    template<class... Args>
    inline void critical_here(const char* f,int l,const char* fn,bool tf_detail, fmt::format_string<Args...> s, Args&&... args){
      mylogger::log_critical_loc_msg(f,l,fn,tf_detail, fmt::format(s, std::forward<Args>(args)...));
    }

    // --- Tagged variants ---
    inline void info_tag_here(const char* f,int l,const char* fn,bool tf_detail,const std::string& tag,const std::string& msg){
      mylogger::log_info_tag_msg(f,l,fn,tf_detail,tag,msg);
    }
    template<class... Args>
    inline void info_tag_here(const char* f,int l,const char* fn,bool tf_detail,const std::string& tag, fmt::format_string<Args...> s, Args&&... args){
      mylogger::log_info_tag_msg(f,l,fn,tf_detail,tag, fmt::format(s, std::forward<Args>(args)...));
    }

    inline void debug_tag_here(const char* f,int l,const char* fn,bool tf_detail,const std::string& tag,const std::string& msg){
      mylogger::log_debug_tag_loc_msg(f,l,fn,tf_detail,tag,msg);
    }
    template<class... Args>
    inline void debug_tag_here(const char* f,int l,const char* fn,bool tf_detail,const std::string& tag, fmt::format_string<Args...> s, Args&&... args){
      mylogger::log_debug_tag_loc_msg(f,l,fn,tf_detail,tag, fmt::format(s, std::forward<Args>(args)...));
    }

    inline void warn_tag_here(const char* f,int l,const char* fn,bool tf_detail,const std::string& tag,const std::string& msg){
      mylogger::log_warn_tag_loc_msg(f,l,fn,tf_detail,tag,msg);
    }
    template<class... Args>
    inline void warn_tag_here(const char* f,int l,const char* fn,bool tf_detail,const std::string& tag, fmt::format_string<Args...> s, Args&&... args){
      mylogger::log_warn_tag_loc_msg(f,l,fn,tf_detail,tag, fmt::format(s, std::forward<Args>(args)...));
    }

    inline void error_tag_here(const char* f,int l,const char* fn,bool tf_detail,const std::string& tag,const std::string& msg){
      mylogger::log_error_tag_loc_msg(f,l,fn,tf_detail,tag,msg);
    }
    template<class... Args>
    inline void error_tag_here(const char* f,int l,const char* fn,bool tf_detail,const std::string& tag, fmt::format_string<Args...> s, Args&&... args){
      mylogger::log_error_tag_loc_msg(f,l,fn,tf_detail,tag, fmt::format(s, std::forward<Args>(args)...));
    }

    inline void critical_tag_here(const char* f,int l,const char* fn,bool tf_detail,const std::string& tag,const std::string& msg){
      mylogger::log_critical_tag_loc_msg(f,l,fn,tf_detail,tag,msg);
    }
    template<class... Args>
    inline void critical_tag_here(const char* f,int l,const char* fn,bool tf_detail,const std::string& tag, fmt::format_string<Args...> s, Args&&... args){
      mylogger::log_critical_tag_loc_msg(f,l,fn,tf_detail,tag, fmt::format(s, std::forward<Args>(args)...));
    }

  } // namespace detail
} // namespace mylogger


// ========================================================================
// Logging Macros - Main Interface
// ========================================================================

/// @defgroup logging_macros Logging Macros
/// @brief Primary interface for logging throughout the application
/// @details These macros automatically capture source location (__FILE__, __LINE__, function name)
///          and support fmt::format syntax for message formatting.
///
///          All macros safely fall back to spdlog's default logger if g_logger is nullptr.
///          The program will never crash due to uninitialized logger.
///
///          Usage examples:
///          @code
///          LOG_INFO("Application started");
///          LOG_ERROR("Failed to open file: {}", filename);
///          LOG_DEBUG_TAG("NETWORK", "Received {} bytes from {}", size, host);
///          @endcode
/// @{

/// @brief Log a TRACE level message with source location
/// @param ... Message (string or fmt::format string with arguments)
#define LOG_TRACE(...)      do { mylogger::detail::trace_here   (mylogger::basename_cstr(__FILE__), __LINE__, __PRETTY_FUNCTION__, true, __VA_ARGS__); } while(0)

/// @brief Log an INFO level message with source location
/// @param ... Message (string or fmt::format string with arguments)
#define LOG_INFO(...)       do { mylogger::detail::info_here    (mylogger::basename_cstr(__FILE__), __LINE__, __PRETTY_FUNCTION__, true, __VA_ARGS__); } while(0)

/// @brief Log a DEBUG level message with source location
/// @param ... Message (string or fmt::format string with arguments)
#define LOG_DEBUG(...)      do { mylogger::detail::debug_here   (mylogger::basename_cstr(__FILE__), __LINE__, __PRETTY_FUNCTION__, true, __VA_ARGS__); } while(0)

/// @brief Log a WARN level message with source location
/// @param ... Message (string or fmt::format string with arguments)
#define LOG_WARN(...)       do { mylogger::detail::warn_here    (mylogger::basename_cstr(__FILE__), __LINE__, __PRETTY_FUNCTION__, true, __VA_ARGS__); } while(0)

/// @brief Log an ERROR level message with source location
/// @param ... Message (string or fmt::format string with arguments)
#define LOG_ERROR(...)      do { mylogger::detail::error_here   (mylogger::basename_cstr(__FILE__), __LINE__, __PRETTY_FUNCTION__, true, __VA_ARGS__); } while(0)

/// @brief Log a CRITICAL level message with source location
/// @param ... Message (string or fmt::format string with arguments)
#define LOG_CRITICAL(...)   do { mylogger::detail::critical_here(mylogger::basename_cstr(__FILE__), __LINE__, __PRETTY_FUNCTION__, true, __VA_ARGS__); } while(0)

/// @brief Log an INFO level message with tag and source location
/// @param tag Tag prefix string (e.g., "NETWORK", "DATABASE")
/// @param ... Message (string or fmt::format string with arguments)
#define LOG_INFO_TAG(tag, ...)     do { mylogger::detail::info_tag_here    (mylogger::basename_cstr(__FILE__), __LINE__, __PRETTY_FUNCTION__, true, tag, __VA_ARGS__); } while(0)

/// @brief Log a DEBUG level message with tag and source location
/// @param tag Tag prefix string (e.g., "NETWORK", "DATABASE")
/// @param ... Message (string or fmt::format string with arguments)
#define LOG_DEBUG_TAG(tag, ...)    do { mylogger::detail::debug_tag_here   (mylogger::basename_cstr(__FILE__), __LINE__, __PRETTY_FUNCTION__, true, tag, __VA_ARGS__); } while(0)

/// @brief Log a WARN level message with tag and source location
/// @param tag Tag prefix string (e.g., "NETWORK", "DATABASE")
/// @param ... Message (string or fmt::format string with arguments)
#define LOG_WARN_TAG(tag, ...)     do { mylogger::detail::warn_tag_here    (mylogger::basename_cstr(__FILE__), __LINE__, __PRETTY_FUNCTION__, true, tag, __VA_ARGS__); } while(0)

/// @brief Log an ERROR level message with tag and source location
/// @param tag Tag prefix string (e.g., "NETWORK", "DATABASE")
/// @param ... Message (string or fmt::format string with arguments)
#define LOG_ERROR_TAG(tag, ...)    do { mylogger::detail::error_tag_here   (mylogger::basename_cstr(__FILE__), __LINE__, __PRETTY_FUNCTION__, true, tag, __VA_ARGS__); } while(0)

/// @brief Log a CRITICAL level message with tag and source location
/// @param tag Tag prefix string (e.g., "NETWORK", "DATABASE")
/// @param ... Message (string or fmt::format string with arguments)
#define LOG_CRITICAL_TAG(tag, ...) do { mylogger::detail::critical_tag_here(mylogger::basename_cstr(__FILE__), __LINE__, __PRETTY_FUNCTION__, true, tag, __VA_ARGS__); } while(0)

/// @}  // end of logging_macros group

// ========================================================================
// Rate-Limited Logging Macros (LOG_*_N)
// ========================================================================

/// @defgroup logging_macros_n Rate-Limited Logging Macros
/// @brief Logging macros with call-site based rate limiting
/// @details These macros limit the number of log outputs from each call site.
///          Uses static atomic counter per call site for thread-safe counting.
///          Once the limit is reached, further calls from that location are silently ignored.
///
///          Usage examples:
///          @code
///          for (int i = 0; i < 10000; ++i) {
///              LOG_WARN_N(5, "Loop iteration: {}", i);  // Only first 5 times
///          }
///          LOG_ERROR_N(1, "This error appears only once");  // Same as LOG_*_ONCE
///          @endcode
///
/// @note Thread-safety: Yes. Uses std::atomic for counter operations.
/// @note The counter persists for the lifetime of the program (static storage).
/// @{

/// @brief Log a TRACE level message, limited to first N calls from this location
/// @param max_count Maximum number of times to output from this call site
/// @param ... Message (string or fmt::format string with arguments)
#define LOG_TRACE_N(max_count, ...) do { \
  static std::atomic<int> _log_n_count_{0}; \
  if (mylogger::should_log_n(_log_n_count_, max_count)) \
    mylogger::detail::trace_here(mylogger::basename_cstr(__FILE__), __LINE__, __PRETTY_FUNCTION__, true, __VA_ARGS__); \
} while(0)

/// @brief Log an INFO level message, limited to first N calls from this location
/// @param max_count Maximum number of times to output from this call site
/// @param ... Message (string or fmt::format string with arguments)
#define LOG_INFO_N(max_count, ...) do { \
  static std::atomic<int> _log_n_count_{0}; \
  if (mylogger::should_log_n(_log_n_count_, max_count)) \
    mylogger::detail::info_here(mylogger::basename_cstr(__FILE__), __LINE__, __PRETTY_FUNCTION__, true, __VA_ARGS__); \
} while(0)

/// @brief Log a DEBUG level message, limited to first N calls from this location
/// @param max_count Maximum number of times to output from this call site
/// @param ... Message (string or fmt::format string with arguments)
#define LOG_DEBUG_N(max_count, ...) do { \
  static std::atomic<int> _log_n_count_{0}; \
  if (mylogger::should_log_n(_log_n_count_, max_count)) \
    mylogger::detail::debug_here(mylogger::basename_cstr(__FILE__), __LINE__, __PRETTY_FUNCTION__, true, __VA_ARGS__); \
} while(0)

/// @brief Log a WARN level message, limited to first N calls from this location
/// @param max_count Maximum number of times to output from this call site
/// @param ... Message (string or fmt::format string with arguments)
#define LOG_WARN_N(max_count, ...) do { \
  static std::atomic<int> _log_n_count_{0}; \
  if (mylogger::should_log_n(_log_n_count_, max_count)) \
    mylogger::detail::warn_here(mylogger::basename_cstr(__FILE__), __LINE__, __PRETTY_FUNCTION__, true, __VA_ARGS__); \
} while(0)

/// @brief Log an ERROR level message, limited to first N calls from this location
/// @param max_count Maximum number of times to output from this call site
/// @param ... Message (string or fmt::format string with arguments)
#define LOG_ERROR_N(max_count, ...) do { \
  static std::atomic<int> _log_n_count_{0}; \
  if (mylogger::should_log_n(_log_n_count_, max_count)) \
    mylogger::detail::error_here(mylogger::basename_cstr(__FILE__), __LINE__, __PRETTY_FUNCTION__, true, __VA_ARGS__); \
} while(0)

/// @brief Log a CRITICAL level message, limited to first N calls from this location
/// @param max_count Maximum number of times to output from this call site
/// @param ... Message (string or fmt::format string with arguments)
#define LOG_CRITICAL_N(max_count, ...) do { \
  static std::atomic<int> _log_n_count_{0}; \
  if (mylogger::should_log_n(_log_n_count_, max_count)) \
    mylogger::detail::critical_here(mylogger::basename_cstr(__FILE__), __LINE__, __PRETTY_FUNCTION__, true, __VA_ARGS__); \
} while(0)

/// @brief Log an INFO level message with tag, limited to first N calls from this location
/// @param max_count Maximum number of times to output from this call site
/// @param tag Tag prefix string (e.g., "NETWORK", "DATABASE")
/// @param ... Message (string or fmt::format string with arguments)
#define LOG_INFO_TAG_N(max_count, tag, ...) do { \
  static std::atomic<int> _log_n_count_{0}; \
  if (mylogger::should_log_n(_log_n_count_, max_count)) \
    mylogger::detail::info_tag_here(mylogger::basename_cstr(__FILE__), __LINE__, __PRETTY_FUNCTION__, true, tag, __VA_ARGS__); \
} while(0)

/// @brief Log a DEBUG level message with tag, limited to first N calls from this location
/// @param max_count Maximum number of times to output from this call site
/// @param tag Tag prefix string (e.g., "NETWORK", "DATABASE")
/// @param ... Message (string or fmt::format string with arguments)
#define LOG_DEBUG_TAG_N(max_count, tag, ...) do { \
  static std::atomic<int> _log_n_count_{0}; \
  if (mylogger::should_log_n(_log_n_count_, max_count)) \
    mylogger::detail::debug_tag_here(mylogger::basename_cstr(__FILE__), __LINE__, __PRETTY_FUNCTION__, true, tag, __VA_ARGS__); \
} while(0)

/// @brief Log a WARN level message with tag, limited to first N calls from this location
/// @param max_count Maximum number of times to output from this call site
/// @param tag Tag prefix string (e.g., "NETWORK", "DATABASE")
/// @param ... Message (string or fmt::format string with arguments)
#define LOG_WARN_TAG_N(max_count, tag, ...) do { \
  static std::atomic<int> _log_n_count_{0}; \
  if (mylogger::should_log_n(_log_n_count_, max_count)) \
    mylogger::detail::warn_tag_here(mylogger::basename_cstr(__FILE__), __LINE__, __PRETTY_FUNCTION__, true, tag, __VA_ARGS__); \
} while(0)

/// @brief Log an ERROR level message with tag, limited to first N calls from this location
/// @param max_count Maximum number of times to output from this call site
/// @param tag Tag prefix string (e.g., "NETWORK", "DATABASE")
/// @param ... Message (string or fmt::format string with arguments)
#define LOG_ERROR_TAG_N(max_count, tag, ...) do { \
  static std::atomic<int> _log_n_count_{0}; \
  if (mylogger::should_log_n(_log_n_count_, max_count)) \
    mylogger::detail::error_tag_here(mylogger::basename_cstr(__FILE__), __LINE__, __PRETTY_FUNCTION__, true, tag, __VA_ARGS__); \
} while(0)

/// @brief Log a CRITICAL level message with tag, limited to first N calls from this location
/// @param max_count Maximum number of times to output from this call site
/// @param tag Tag prefix string (e.g., "NETWORK", "DATABASE")
/// @param ... Message (string or fmt::format string with arguments)
#define LOG_CRITICAL_TAG_N(max_count, tag, ...) do { \
  static std::atomic<int> _log_n_count_{0}; \
  if (mylogger::should_log_n(_log_n_count_, max_count)) \
    mylogger::detail::critical_tag_here(mylogger::basename_cstr(__FILE__), __LINE__, __PRETTY_FUNCTION__, true, tag, __VA_ARGS__); \
} while(0)

/// @}  // end of logging_macros_n group

// ========================================================================
// Rate-Limited Logging Macros with Default Count (LOG_*_ND)
// ========================================================================

/// @defgroup logging_macros_nd Rate-Limited Logging Macros with Default Count
/// @brief Logging macros using config-based default rate limit
/// @details These macros use g_default_max_count (configurable via JSON or set_default_max_count())
///          instead of requiring an explicit count parameter.
///
///          Usage examples:
///          @code
///          // Initialize default (optional, defaults to 10)
///          mylogger::configure_rate_limit(config, "LOG_FILE");
///
///          for (int i = 0; i < 10000; ++i) {
///              LOG_WARN_ND("Loop warning: {}", i);  // Uses g_default_max_count
///          }
///          @endcode
///
/// @note Thread-safety: Yes. Uses std::atomic for counter operations.
/// @note ND = N Default (uses g_default_max_count instead of explicit N)
/// @{

/// @brief Log a TRACE level message, limited by g_default_max_count
/// @param ... Message (string or fmt::format string with arguments)
#define LOG_TRACE_ND(...) do { \
  static std::atomic<int> _log_n_count_{0}; \
  if (mylogger::should_log_n(_log_n_count_, mylogger::g_default_max_count)) \
    mylogger::detail::trace_here(mylogger::basename_cstr(__FILE__), __LINE__, __PRETTY_FUNCTION__, true, __VA_ARGS__); \
} while(0)

/// @brief Log an INFO level message, limited by g_default_max_count
/// @param ... Message (string or fmt::format string with arguments)
#define LOG_INFO_ND(...) do { \
  static std::atomic<int> _log_n_count_{0}; \
  if (mylogger::should_log_n(_log_n_count_, mylogger::g_default_max_count)) \
    mylogger::detail::info_here(mylogger::basename_cstr(__FILE__), __LINE__, __PRETTY_FUNCTION__, true, __VA_ARGS__); \
} while(0)

/// @brief Log a DEBUG level message, limited by g_default_max_count
/// @param ... Message (string or fmt::format string with arguments)
#define LOG_DEBUG_ND(...) do { \
  static std::atomic<int> _log_n_count_{0}; \
  if (mylogger::should_log_n(_log_n_count_, mylogger::g_default_max_count)) \
    mylogger::detail::debug_here(mylogger::basename_cstr(__FILE__), __LINE__, __PRETTY_FUNCTION__, true, __VA_ARGS__); \
} while(0)

/// @brief Log a WARN level message, limited by g_default_max_count
/// @param ... Message (string or fmt::format string with arguments)
#define LOG_WARN_ND(...) do { \
  static std::atomic<int> _log_n_count_{0}; \
  if (mylogger::should_log_n(_log_n_count_, mylogger::g_default_max_count)) \
    mylogger::detail::warn_here(mylogger::basename_cstr(__FILE__), __LINE__, __PRETTY_FUNCTION__, true, __VA_ARGS__); \
} while(0)

/// @brief Log an ERROR level message, limited by g_default_max_count
/// @param ... Message (string or fmt::format string with arguments)
#define LOG_ERROR_ND(...) do { \
  static std::atomic<int> _log_n_count_{0}; \
  if (mylogger::should_log_n(_log_n_count_, mylogger::g_default_max_count)) \
    mylogger::detail::error_here(mylogger::basename_cstr(__FILE__), __LINE__, __PRETTY_FUNCTION__, true, __VA_ARGS__); \
} while(0)

/// @brief Log a CRITICAL level message, limited by g_default_max_count
/// @param ... Message (string or fmt::format string with arguments)
#define LOG_CRITICAL_ND(...) do { \
  static std::atomic<int> _log_n_count_{0}; \
  if (mylogger::should_log_n(_log_n_count_, mylogger::g_default_max_count)) \
    mylogger::detail::critical_here(mylogger::basename_cstr(__FILE__), __LINE__, __PRETTY_FUNCTION__, true, __VA_ARGS__); \
} while(0)

/// @brief Log an INFO level message with tag, limited by g_default_max_count
/// @param tag Tag prefix string (e.g., "NETWORK", "DATABASE")
/// @param ... Message (string or fmt::format string with arguments)
#define LOG_INFO_TAG_ND(tag, ...) do { \
  static std::atomic<int> _log_n_count_{0}; \
  if (mylogger::should_log_n(_log_n_count_, mylogger::g_default_max_count)) \
    mylogger::detail::info_tag_here(mylogger::basename_cstr(__FILE__), __LINE__, __PRETTY_FUNCTION__, true, tag, __VA_ARGS__); \
} while(0)

/// @brief Log a DEBUG level message with tag, limited by g_default_max_count
/// @param tag Tag prefix string (e.g., "NETWORK", "DATABASE")
/// @param ... Message (string or fmt::format string with arguments)
#define LOG_DEBUG_TAG_ND(tag, ...) do { \
  static std::atomic<int> _log_n_count_{0}; \
  if (mylogger::should_log_n(_log_n_count_, mylogger::g_default_max_count)) \
    mylogger::detail::debug_tag_here(mylogger::basename_cstr(__FILE__), __LINE__, __PRETTY_FUNCTION__, true, tag, __VA_ARGS__); \
} while(0)

/// @brief Log a WARN level message with tag, limited by g_default_max_count
/// @param tag Tag prefix string (e.g., "NETWORK", "DATABASE")
/// @param ... Message (string or fmt::format string with arguments)
#define LOG_WARN_TAG_ND(tag, ...) do { \
  static std::atomic<int> _log_n_count_{0}; \
  if (mylogger::should_log_n(_log_n_count_, mylogger::g_default_max_count)) \
    mylogger::detail::warn_tag_here(mylogger::basename_cstr(__FILE__), __LINE__, __PRETTY_FUNCTION__, true, tag, __VA_ARGS__); \
} while(0)

/// @brief Log an ERROR level message with tag, limited by g_default_max_count
/// @param tag Tag prefix string (e.g., "NETWORK", "DATABASE")
/// @param ... Message (string or fmt::format string with arguments)
#define LOG_ERROR_TAG_ND(tag, ...) do { \
  static std::atomic<int> _log_n_count_{0}; \
  if (mylogger::should_log_n(_log_n_count_, mylogger::g_default_max_count)) \
    mylogger::detail::error_tag_here(mylogger::basename_cstr(__FILE__), __LINE__, __PRETTY_FUNCTION__, true, tag, __VA_ARGS__); \
} while(0)

/// @brief Log a CRITICAL level message with tag, limited by g_default_max_count
/// @param tag Tag prefix string (e.g., "NETWORK", "DATABASE")
/// @param ... Message (string or fmt::format string with arguments)
#define LOG_CRITICAL_TAG_ND(tag, ...) do { \
  static std::atomic<int> _log_n_count_{0}; \
  if (mylogger::should_log_n(_log_n_count_, mylogger::g_default_max_count)) \
    mylogger::detail::critical_tag_here(mylogger::basename_cstr(__FILE__), __LINE__, __PRETTY_FUNCTION__, true, tag, __VA_ARGS__); \
} while(0)

/// @}  // end of logging_macros_nd group
