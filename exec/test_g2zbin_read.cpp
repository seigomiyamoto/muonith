/// @file test_g2zbin_read.cpp
/// @brief Verify that Grid2dXYZ reads g2zbin v2 files in both z precisions.
/// @details
/// Background: the GeoTIFF conversion tools (relocated to the muonith-gsi-dem
/// repository) write .g2zbin version 2 files with either float32 (precision
/// byte 4) or float64 (precision byte 8) z values.  This test guards that
/// Grid2dXYZ::load_g2zbin() accepts both precisions, preserves NaN cells,
/// and returns values matching what save_g2zbin() wrote.
/// Usage: test_g2zbin_read.exe [external.g2zbin]
///   With an argument, additionally loads the given file (e.g. one produced
///   by the Python writer scripts/g2zbin_io.py in muonith-gsi-dem) and prints
///   its grid summary for cross-tool verification.

#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include "cls_Grid1d.hpp"
#include "cls_Grid2dXYZ.hpp"
#include "ns_mylogger.hpp"
#include "spdlog_pch.hpp"

namespace fs = std::filesystem;

namespace {

// Explicit check (not assert) so the test also fails under NDEBUG builds.
// Raw std::runtime_error (not THROW_ERROR) is intentional: this is a
// standalone test harness and main() catches it to return exit code 1.
void check(bool cond, const std::string& msg)
{
  if (!cond) {
    throw std::runtime_error("check failed: " + msg);
  }
}

bool bit_equal(double a, double b)
{
  std::uint64_t ua = 0;
  std::uint64_t ub = 0;
  std::memcpy(&ua, &a, sizeof(double));
  std::memcpy(&ub, &b, sizeof(double));
  return ua == ub;
}

// Byte 10 of the fixed header is the z precision (4=float32, 8=float64).
std::uint8_t read_precision_byte(const fs::path& path)
{
  std::ifstream ifs(path, std::ios::binary);
  check(static_cast<bool>(ifs), "cannot open: " + path.string());
  ifs.seekg(10);
  char precision = 0;
  ifs.read(&precision, 1);
  check(static_cast<bool>(ifs), "cannot read precision byte: " + path.string());
  return static_cast<std::uint8_t>(precision);
}

// Sample grid: x centers 101,103,105 / y centers 197,199 (v2 outer edges),
// z = 1..6 row-major with one NaN cell at (ix=1, iy=0).
Grid2dXYZ make_sample_grid()
{
  const Grid1d x_axis("x", 3, 100.0, 106.0, 2.0);
  const Grid1d y_axis("y", 2, 196.0, 200.0, 2.0);
  Grid2dXYZ org(x_axis, y_axis);
  double z = 1.0;
  for (int iy = 0; iy < 2; ++iy) {
    for (int ix = 0; ix < 3; ++ix) {
      org.set_z(ix, iy, z);
      z += 1.0;
    }
  }
  org.set_z(1, 0, std::nan(""));
  return org;
}

void check_axes_match(const Grid2dXYZ& org, const Grid2dXYZ& in)
{
  check(in.get_nbinx() == org.get_nbinx(), "nbinx mismatch");
  check(in.get_nbiny() == org.get_nbiny(), "nbiny mismatch");
  check(bit_equal(in.get_x_axis().get_min(), org.get_x_axis().get_min()),
        "x min not bit-exact");
  check(bit_equal(in.get_x_axis().get_max(), org.get_x_axis().get_max()),
        "x max not bit-exact");
  check(bit_equal(in.get_y_axis().get_min(), org.get_y_axis().get_min()),
        "y min not bit-exact");
  check(bit_equal(in.get_y_axis().get_max(), org.get_y_axis().get_max()),
        "y max not bit-exact");
}

void test_read_float32(const fs::path& dir)
{
  std::cout << "[test_read_float32] -----------------------------------------" << std::endl;
  const Grid2dXYZ org = make_sample_grid();
  const fs::path path = dir / "sample_f32.g2zbin";
  org.save_g2zbin(path, false);
  check(read_precision_byte(path) == 4, "precision byte is not 4 (float32)");

  const Grid2dXYZ in(path, false, false);
  check_axes_match(org, in);
  for (int iy = 0; iy < org.get_nbiny(); ++iy) {
    for (int ix = 0; ix < org.get_nbinx(); ++ix) {
      const double z_org = org.get_z(ix, iy);
      const double z_in = in.get_z(ix, iy);
      if (std::isnan(z_org)) {
        check(std::isnan(z_in), "NaN cell not preserved (float32)");
      } else {
        // float32 storage: expect the value after a float round-trip.
        const double z_expected = static_cast<double>(static_cast<float>(z_org));
        check(bit_equal(z_in, z_expected), "z value mismatch (float32)");
      }
    }
  }
  std::cout << "  -> ok" << std::endl;
}

void test_read_float64(const fs::path& dir)
{
  std::cout << "[test_read_float64] -----------------------------------------" << std::endl;
  const Grid2dXYZ org = make_sample_grid();
  const fs::path path = dir / "sample_f64.g2zbin";
  org.save_g2zbin(path, true);
  check(read_precision_byte(path) == 8, "precision byte is not 8 (float64)");

  const Grid2dXYZ in(path, false, false);
  check_axes_match(org, in);
  for (int iy = 0; iy < org.get_nbiny(); ++iy) {
    for (int ix = 0; ix < org.get_nbinx(); ++ix) {
      const double z_org = org.get_z(ix, iy);
      const double z_in = in.get_z(ix, iy);
      if (std::isnan(z_org)) {
        check(std::isnan(z_in), "NaN cell not preserved (float64)");
      } else {
        check(bit_equal(z_in, z_org), "z value not bit-exact (float64)");
      }
    }
  }
  std::cout << "  -> ok" << std::endl;
}

// Load an externally produced .g2zbin (e.g. from the Python writer in
// muonith-gsi-dem) and print its summary for cross-tool verification.
void load_external(const fs::path& path)
{
  std::cout << "[load_external] " << path.string() << std::endl;
  const std::uint8_t precision = read_precision_byte(path);
  const Grid2dXYZ in(path, false, false);
  std::cout << "  precision=" << static_cast<int>(precision)
            << " nbinx=" << in.get_nbinx()
            << " nbiny=" << in.get_nbiny() << std::endl;
  std::cout << "  x: [" << in.get_x_axis().get_min() << ", "
            << in.get_x_axis().get_max() << "]"
            << " y: [" << in.get_y_axis().get_min() << ", "
            << in.get_y_axis().get_max() << "]" << std::endl;
  double z_min = std::nan("");
  double z_max = std::nan("");
  long n_nan = 0;
  for (int iy = 0; iy < in.get_nbiny(); ++iy) {
    for (int ix = 0; ix < in.get_nbinx(); ++ix) {
      const double z = in.get_z(ix, iy);
      if (std::isnan(z)) {
        ++n_nan;
      } else {
        z_min = std::isnan(z_min) ? z : std::min(z_min, z);
        z_max = std::isnan(z_max) ? z : std::max(z_max, z);
      }
    }
  }
  std::cout << "  z: [" << z_min << ", " << z_max << "]"
            << " nan_cells=" << n_nan << std::endl;
  std::cout << "  -> ok" << std::endl;
}

} // namespace

int main(int argc, char* argv[])
{
  auto logger = mylogger::create_logger(
      spdlog::level::info,
      spdlog::level::err,
      spdlog::level::info,
      "logs/test_g2zbin_read.log");
  (void)logger;

  const fs::path dir = fs::temp_directory_path() / "test_g2zbin_read";
  fs::create_directories(dir);

  try {
    test_read_float32(dir);
    test_read_float64(dir);
    if (argc > 1) {
      load_external(argv[1]);
    }
  } catch (const std::exception& ex) {
    std::cerr << "Test failed: " << ex.what() << std::endl;
    return 1;
  }

  std::cout << "All Grid2dXYZ g2zbin read tests passed." << std::endl;
  return 0;
}
