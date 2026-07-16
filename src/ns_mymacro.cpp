// mymacro.cpp
#include <stdexcept>
#include <sstream>
#include <thread>
#include <tuple>
#include "ns_mymacro.hpp"
#include "ns_mylogger.hpp"
#include <spdlog/sinks/basic_file_sink.h>  // for file sink
#include <spdlog/sinks/stdout_color_sinks.h>
#include <nlohmann/json.hpp>  // for json
#include "spdlog_pch.hpp"

// throw_error
void mymacro::throw_error(const std::string& msg
, const std::string& var_name, const std::string& var_value
, const std::string& function, const std::string& file, const int line)
{
  std::ostringstream oss;
  oss << function << " at " << file << ":" << line << "\n"
      << msg << "\n"
      << "current value of " << var_name << " is " << var_value << "\n";
  if (mylogger::g_logger) {
    mylogger::g_logger->error(oss.str());
  }
  throw std::runtime_error(oss.str());
}

void mymacro::throw_error(const std::string& msg
, const std::string& var1_name, const std::string& var1_value
, const std::string& var2_name, const std::string& var2_value
, const std::string& function, const std::string& file, const int line)
{
  std::ostringstream oss;
  oss << function << " at " << file << ":" << line << "\n"
      << msg << "\n"
      << "Current value of " << var1_name << " is " << var1_value << "\n"
      << "Current value of " << var2_name << " is " << var2_value << "\n";
  if (mylogger::g_logger) {
    mylogger::g_logger->error(oss.str());
  }
  throw std::runtime_error(oss.str());
}

// throw_error with instance name
void mymacro::throw_error_name(const std::string& msg
, const std::string& var_name, const std::string& var_value
, const std::string& function, const std::string& file
, const int line, const std::string &instance_name)
{
  std::ostringstream oss;
  oss << function << " at " << file << ":" << line << "\n"
      << msg << "\n"
      << "Instance name: " << instance_name << "\n"
      << "current value of " << var_name << " is " << var_value << "\n";
  if (mylogger::g_logger) {
    mylogger::g_logger->error(oss.str());
  }
  throw std::runtime_error(oss.str());
}

// throw_error with instance name
void mymacro::throw_error_name(const std::string& msg
, const std::string& var1_name, const std::string& var1_value
, const std::string& var2_name, const std::string& var2_value
, const std::string& function, const std::string& file
, const int line, const std::string &instance_name)
{
  std::ostringstream oss;
  oss << function << " at " << file << ":" << line << "\n"
      << msg << "\n"
      << "Instance name: " << instance_name << "\n"
      << "Current value of " << var1_name << " is " << var1_value << "\n"
      << "Current value of " << var2_name << " is " << var2_value << "\n";
  if (mylogger::g_logger) {
    mylogger::g_logger->error(oss.str());
  }
  throw std::runtime_error(oss.str());
}

void mymacro::disp_mat_size(
        const std::string &mat_name
      , const Eigen::MatrixXf &mat
      , const spdlog::level::level_enum log_level)
{
  mylogger::g_logger->log(log_level, "mat_name={} rows={} cols={}",mat_name,mat.rows(),mat.cols());
}
