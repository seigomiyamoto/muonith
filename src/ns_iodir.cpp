// src/ns_iodir.cpp
#include "ns_iodir.hpp"
#include "ns_mymacro.hpp"

/// @brief Anonymous namespace for internal state
namespace {
  fs::path s_default_output_dir{"./tmp"};  // Default value
}

void iodir::set_default_output_dir(const fs::path &dir)
{
  if (dir.empty()) {
    s_default_output_dir = "./tmp";
  } else {
    s_default_output_dir = dir;
  }
}

const fs::path &iodir::get_default_output_dir() {
  return s_default_output_dir;
}

fs::path iodir::make_pathout(
  const fs::path& file_path, const fs::path& output_dir)
{
  if (file_path.is_absolute()) {
    return file_path;
  }
  const auto full_path = output_dir / file_path;
  const auto dir = full_path.parent_path();
  if (!dir.empty() && !std::filesystem::exists(dir)) {
    std::error_code ec;
    std::filesystem::create_directories(dir, ec);
    if (ec) {
      THROW_ERROR("iodir::make_pathout: Failed to create directory. path={} error={}", dir.string(), ec.message());
    }
  }
  return full_path;
}

fs::path iodir::make_pathout( const char* cfname, const fs::path& output_dir)
{
  return make_pathout( fs::path(cfname), output_dir);
}

fs::path iodir::make_pathout( const std::string& fname_str, const fs::path& output_dir)
{
  return make_pathout( fs::path(fname_str), output_dir);
}

fs::path iodir::make_pathout( const fs::path &file_path)
{
  return make_pathout(file_path, s_default_output_dir);
}

fs::path iodir::make_pathout( const char* cfname )
{
  return make_pathout( fs::path(cfname) );
}

fs::path iodir::make_pathout( const std::string& fname_str )
{
  return make_pathout( fs::path(fname_str) );
}