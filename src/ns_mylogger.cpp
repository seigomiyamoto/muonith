// src/mylogger.cpp
#include "ns_mylogger.hpp"
#include <iostream>
#include <tuple>
#include <stdexcept>
#include <fmt/format.h>
#include <fmt/std.h>   // Enable fmt support for std::filesystem::path
#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>

namespace mylogger {
  std::shared_ptr<spdlog::logger> g_logger = nullptr;
}

std::string mylogger::extract_class_method(const char* pretty_func) {
  // Input examples (GCC/Clang):
  //   "void ClassName::method(int)"
  //   "static void ClassName::method()"
  //   "void ClassName::method() const"
  //   "void free_function(int)"
  //   "auto ClassName::method(int) -> ReturnType"
  //
  // Goal: Extract "ClassName::method" or "free_function"

  std::string_view sv(pretty_func);

  // Find '(' to get the end of function name
  auto paren_pos = sv.find('(');
  if (paren_pos == std::string_view::npos) {
    return std::string(pretty_func);  // Fallback
  }

  // Work backwards from '(' to find function name start
  // Skip any template parameters (< >)
  size_t end_pos = paren_pos;
  int angle_depth = 0;
  while (end_pos > 0) {
    --end_pos;
    char c = sv[end_pos];
    if (c == '>') {
      ++angle_depth;
    } else if (c == '<') {
      --angle_depth;
    } else if (angle_depth == 0 && (c == ' ' || c == '*' || c == '&')) {
      ++end_pos;  // Move back to first char of name
      break;
    }
  }
  if (end_pos == 0 && angle_depth == 0) {
    // No space found, name starts at beginning
    end_pos = 0;
  }

  // Extract from end_pos to paren_pos
  std::string_view name_part = sv.substr(end_pos, paren_pos - end_pos);

  // Remove leading spaces if any
  while (!name_part.empty() && name_part.front() == ' ') {
    name_part.remove_prefix(1);
  }

  return std::string(name_part);
}

// ===== Log message implementation (non-template, defined in .cpp) =====

namespace {
  // Log with [file:line ClassName::method] prefix (for DEBUG/WARN/ERROR/TRACE/CRITICAL)
  inline void log_with_file_line_func(
      spdlog::level::level_enum level
    , const char* file
    , int line
    , const char* func
    , bool tf_show_detail
    , const std::string& formatted_message)
  {
    std::string final_msg;
    if (tf_show_detail) {
      std::string class_method = mylogger::extract_class_method(func);
      final_msg = fmt::format("[{}:{} {}] {}", file, line, class_method, formatted_message);
    } else {
      final_msg = formatted_message;
    }

    if (mylogger::g_logger) {
      mylogger::g_logger->log(level, "{}", final_msg);
    } else {
      spdlog::log(level, "{}", final_msg);
    }
  }

  // Log with [ClassName::method] prefix only (for INFO)
  inline void log_with_func_only(
      spdlog::level::level_enum level
    , const char* /*file*/
    , int /*line*/
    , const char* func
    , bool tf_show_detail
    , const std::string& formatted_message)
  {
    std::string final_msg;
    if (tf_show_detail) {
      std::string class_method = mylogger::extract_class_method(func);
      final_msg = fmt::format("[{}] {}", class_method, formatted_message);
    } else {
      final_msg = formatted_message;
    }

    if (mylogger::g_logger) {
      mylogger::g_logger->log(level, "{}", final_msg);
    } else {
      spdlog::log(level, "{}", final_msg);
    }
  }

  // Log with [file:line ClassName::method] and [tag] prefix
  inline void log_with_file_line_func_tagged(
      spdlog::level::level_enum level
    , const char* file
    , int line
    , const char* func
    , bool tf_show_detail
    , const std::string& tag
    , const std::string& formatted_message)
  {
    std::string final_msg;
    if (tf_show_detail) {
      std::string class_method = mylogger::extract_class_method(func);
      final_msg = fmt::format("[{}:{} {}] [{}] {}", file, line, class_method, tag, formatted_message);
    } else {
      final_msg = fmt::format("[{}] {}", tag, formatted_message);
    }

    if (mylogger::g_logger) {
      mylogger::g_logger->log(level, "{}", final_msg);
    } else {
      spdlog::log(level, "{}", final_msg);
    }
  }

  // Log with [ClassName::method] and [tag] prefix (for INFO_TAG)
  inline void log_with_func_only_tagged(
      spdlog::level::level_enum level
    , const char* /*file*/
    , int /*line*/
    , const char* func
    , bool tf_show_detail
    , const std::string& tag
    , const std::string& formatted_message)
  {
    std::string final_msg;
    if (tf_show_detail) {
      std::string class_method = mylogger::extract_class_method(func);
      final_msg = fmt::format("[{}] [{}] {}", class_method, tag, formatted_message);
    } else {
      final_msg = fmt::format("[{}] {}", tag, formatted_message);
    }

    if (mylogger::g_logger) {
      mylogger::g_logger->log(level, "{}", final_msg);
    } else {
      spdlog::log(level, "{}", final_msg);
    }
  }
} // namespace

void mylogger::log_trace_msg(const char* file, int line, const char* func, bool tf_show_detail, const std::string& msg) {
  log_with_file_line_func(spdlog::level::trace, file, line, func, tf_show_detail, msg);
}

void mylogger::log_info_msg(const char* file, int line, const char* func, bool tf_show_detail, const std::string& msg) {
  log_with_func_only(spdlog::level::info, file, line, func, tf_show_detail, msg);
}

void mylogger::log_info_tag_msg(const char* file, int line, const char* func, bool tf_show_detail, const std::string& tag, const std::string& msg) {
  log_with_func_only_tagged(spdlog::level::info, file, line, func, tf_show_detail, tag, msg);
}

void mylogger::log_debug_loc_msg(const char* file, int line, const char* func, bool tf_show_detail, const std::string& msg) {
  log_with_file_line_func(spdlog::level::debug, file, line, func, tf_show_detail, msg);
}

void mylogger::log_warn_loc_msg(const char* file, int line, const char* func, bool tf_show_detail, const std::string& msg) {
  log_with_file_line_func(spdlog::level::warn, file, line, func, tf_show_detail, msg);
}

void mylogger::log_error_loc_msg(const char* file, int line, const char* func, bool tf_show_detail, const std::string& msg) {
  log_with_file_line_func(spdlog::level::err, file, line, func, tf_show_detail, msg);
}

void mylogger::log_critical_loc_msg(const char* file, int line, const char* func, bool tf_show_detail, const std::string& msg) {
  log_with_file_line_func(spdlog::level::critical, file, line, func, tf_show_detail, msg);
}

void mylogger::log_debug_tag_loc_msg(const char* file, int line, const char* func, bool tf_show_detail, const std::string& tag, const std::string& msg) {
  log_with_file_line_func_tagged(spdlog::level::debug, file, line, func, tf_show_detail, tag, msg);
}

void mylogger::log_warn_tag_loc_msg(const char* file, int line, const char* func, bool tf_show_detail, const std::string& tag, const std::string& msg) {
  log_with_file_line_func_tagged(spdlog::level::warn, file, line, func, tf_show_detail, tag, msg);
}

void mylogger::log_error_tag_loc_msg(const char* file, int line, const char* func, bool tf_show_detail, const std::string& tag, const std::string& msg) {
  log_with_file_line_func_tagged(spdlog::level::err, file, line, func, tf_show_detail, tag, msg);
}

void mylogger::log_critical_tag_loc_msg(const char* file, int line, const char* func, bool tf_show_detail, const std::string& tag, const std::string& msg) {
  log_with_file_line_func_tagged(spdlog::level::critical, file, line, func, tf_show_detail, tag, msg);
}

// ===== Public API functions =====

std::filesystem::path 
mylogger::generate_unique_log_path(
    const std::filesystem::path log_dir,
    const std::string& base_name, const std::string& extension,
    const bool give_new_number)
{
  if (!std::filesystem::exists(log_dir)) {
    std::filesystem::create_directories(log_dir);
  }
  std::filesystem::path log_path = log_dir / (base_name + extension);
  if (!give_new_number) return log_path;

  int counter = 0;
  while (std::filesystem::exists(log_path)) {
    log_path = log_dir / fmt::format("{}_{:04d}{}", base_name, counter++, extension);
  }
  return log_path;
}

std::shared_ptr<spdlog::logger>
mylogger::create_logger(
    spdlog::level::level_enum stdout_level,
    spdlog::level::level_enum stderr_level,
    spdlog::level::level_enum file_level,
    const std::filesystem::path& log_file_path,
    bool archive_existing)
{
  // If archive_existing=true and file exists, rename existing file with numbered suffix
  if (archive_existing && std::filesystem::exists(log_file_path)) {
    auto parent = log_file_path.parent_path();
    auto stem = log_file_path.stem().string();
    auto ext = log_file_path.extension().string();

    // Find next available number
    int counter = 1;
    std::filesystem::path archived_path;
    do {
      archived_path = parent / fmt::format("{}_{:04d}{}", stem, counter++, ext);
    } while (std::filesystem::exists(archived_path));

    std::filesystem::rename(log_file_path, archived_path);
  }

  auto stdout_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
  stdout_sink->set_level(stdout_level);

  // Use dedicated stderr sink to separate from stdout
  auto stderr_sink = std::make_shared<spdlog::sinks::stderr_color_sink_mt>();
  stderr_sink->set_level(stderr_level);

  // Always create new file (truncate=true) since we renamed existing file above
  auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(log_file_path.string(), /*truncate=*/true);
  file_sink->set_level(file_level);

  std::vector<spdlog::sink_ptr> sinks{stdout_sink, stderr_sink, file_sink};
  auto logger = std::make_shared<spdlog::logger>("multi_sink", sinks.begin(), sinks.end());

  const auto overall_level = std::min({stdout_level, stderr_level, file_level});
  logger->set_level(overall_level);

  // Pattern: [datetime] [level] message
  // Source location (file:line, function) is formatted by log_*_msg functions
  logger->set_pattern("[%Y-%m-%d %H:%M:%S] [%^%l%$] %v");
  logger->flush_on(spdlog::level::trace);

  mylogger::g_logger = logger;
  spdlog::set_default_logger(logger);

  LOG_INFO("Logger setup complete. default_logger={} sinks={}, stdout_level={}, stderr_level={}, file_level={}",
           static_cast<void*>(spdlog::default_logger().get()),
           sinks.size(),
           static_cast<int>(stdout_level),
           static_cast<int>(stderr_level),
           static_cast<int>(file_level));

  return logger;
}

std::tuple<spdlog::level::level_enum, spdlog::level::level_enum, spdlog::level::level_enum>
mylogger::get_log_levels(const nlohmann::json &config, const std::string &section_name)
{
  const auto& logConfig = config.at(section_name).at("log_level");
  auto parse = [](const std::string &s)->spdlog::level::level_enum{
    if (s=="trace") return spdlog::level::trace;
    if (s=="debug") return spdlog::level::debug;
    if (s=="info")  return spdlog::level::info;
    if (s=="warn"||s=="warning") return spdlog::level::warn;
    if (s=="err" ||s=="error")   return spdlog::level::err;
    if (s=="critical")           return spdlog::level::critical;
    if (s=="off") return spdlog::level::off;
    return spdlog::level::info;
  };
  return {
    parse(logConfig.value("stdout_level","info")),
    parse(logConfig.value("stderr_level","warn")),
    parse(logConfig.value("file_level",  "trace"))
  };
}

bool mylogger::get_log_archive_existing(const nlohmann::json& config, const std::string& section_name) {
  if (!config.contains(section_name)) return false;
  const auto& section = config.at(section_name);
  return section.value("archive_existing", false);
}

void mylogger::safe_log(spdlog::level::level_enum level, const std::string &msg) {
  if (g_logger) {
    g_logger->log(level, "{}", msg);
  } else {
    spdlog::log(level, "[fallback] logger is nullptr. fallback log:\n{}", msg);
  }
}

bool mylogger::should_log_n(std::atomic<int>& counter, int max_count) {
  // Fast path: already exceeded limit
  if (counter.load(std::memory_order_relaxed) >= max_count) [[likely]] {
    return false;
  }
  // Actual increment and check
  return counter.fetch_add(1, std::memory_order_relaxed) < max_count;
}

// ===== Rate-limit default configuration =====

// C++20: constinit ensures zero-initialization at compile time
constinit int mylogger::g_default_max_count = 100;

[[nodiscard]] int mylogger::get_default_max_count() noexcept {
  return g_default_max_count;
}

void mylogger::set_default_max_count(int value) noexcept {
  g_default_max_count = value;
}

void mylogger::configure_rate_limit(const nlohmann::json& config, const std::string& section_name) {
  if (!config.contains(section_name)) return;
  const auto& section = config.at(section_name);
  g_default_max_count = section.value("default_max_count", 100);
}
