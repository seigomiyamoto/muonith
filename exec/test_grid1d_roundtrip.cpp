/// @file test_grid1d_roundtrip.cpp
/// @brief Verify that Grid1d/Grid2d save+load is bit-exact pass-through.
/// @details
/// Background: Grid1d used to hold a RAM-only tf_half_shift flag; save()
/// reverse-transformed min/max by +/-0.5*interval and load() re-applied the
/// shift from external JSON5 flags.  The two authorities could disagree
/// (checkpoint_m3 resume bug), and the transform itself could drift by
/// ~1 ULP per roundtrip.  The flag was removed: canonical coordinates are
/// now written and read as-is.  This test guards that contract.

#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include "cls_Grid1d.hpp"
#include "cls_Grid2d.hpp"
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

void check_grid1d_bit_equal(const Grid1d& org, const Grid1d& in)
{
  check(in.get_name() == org.get_name(), "name mismatch");
  check(in.get_nbin() == org.get_nbin(), "nbin mismatch");
  check(bit_equal(in.get_min(), org.get_min()), "min not bit-exact");
  check(bit_equal(in.get_max(), org.get_max()), "max not bit-exact");
  check(bit_equal(in.get_interval(), org.get_interval()), "interval not bit-exact");
}

Grid1d roundtrip(const Grid1d& org, const fs::path& path)
{
  {
    std::ofstream ofs(path, std::ios::binary);
    check(static_cast<bool>(ofs), "cannot open output: " + path.string());
    org.save(ofs);
  }
  Grid1d in;
  {
    std::ifstream ifs(path, std::ios::binary);
    check(static_cast<bool>(ifs), "cannot open input: " + path.string());
    in.load(ifs);
  }
  return in;
}

void test_grid1d_roundtrip_simple(const fs::path& dir)
{
  std::cout << "[test_grid1d_roundtrip_simple] ------------------------------" << std::endl;
  const Grid1d org("x_axis", 250, -2500.0, 2500.0, 20.0);
  const Grid1d in = roundtrip(org, dir / "g1d_simple.bin");
  check_grid1d_bit_equal(org, in);
  std::cout << "  -> ok" << std::endl;
}

void test_grid1d_roundtrip_shifted_edges(const fs::path& dir)
{
  std::cout << "[test_grid1d_roundtrip_shifted_edges] -----------------------" << std::endl;
  // Canonical (already half-interval shifted) edges as produced by the
  // bake-in construction paths: min = center_min - 0.5*interval.
  const double interval = 20.0;
  const double min = -2510.0;  // -2500.0 - 0.5*20.0
  const double max =  2490.0;  //  2500.0 - 0.5*20.0
  const Grid1d org("x_axis_shifted", 250, min, max, interval);
  const Grid1d in = roundtrip(org, dir / "g1d_shifted.bin");
  check_grid1d_bit_equal(org, in);
  std::cout << "  -> ok" << std::endl;
}

void test_grid1d_repeated_roundtrip_no_drift(const fs::path& dir)
{
  std::cout << "[test_grid1d_repeated_roundtrip_no_drift] -------------------" << std::endl;
  // Values with non-terminating binary fractions: the old +/-0.5*interval
  // save/load transform could drift ~1 ULP per cycle on such values.
  const int nbin = 137;
  const double min = -100.1;
  const double max = 1271.3;
  const Grid1d org("drift_axis", nbin, min, max, (max - min) / nbin);
  Grid1d cur = org;
  for (int i = 0; i < 10; ++i) {
    cur = roundtrip(cur, dir / "g1d_drift.bin");
  }
  check_grid1d_bit_equal(org, cur);
  std::cout << "  -> ok (10 cycles, bit-exact)" << std::endl;
}

void test_grid2d_roundtrip(const fs::path& dir)
{
  std::cout << "[test_grid2d_roundtrip] -------------------------------------" << std::endl;
  const Grid1d x_axis("x_axis", 250, -2510.0, 2490.0, 20.0);
  const Grid1d y_axis("y_axis", 300, -3010.0, 2990.0, 20.0);
  const Grid2d org(x_axis, y_axis);

  const fs::path path = dir / "g2d.bin";
  {
    std::ofstream ofs(path, std::ios::binary);
    check(static_cast<bool>(ofs), "cannot open output: " + path.string());
    org.save(ofs);
  }
  Grid2d in;
  {
    std::ifstream ifs(path, std::ios::binary);
    check(static_cast<bool>(ifs), "cannot open input: " + path.string());
    in.load(ifs);
  }
  check(in.get_name() == org.get_name(), "Grid2d name mismatch");
  check_grid1d_bit_equal(org.get_x_axis(), in.get_x_axis());
  check_grid1d_bit_equal(org.get_y_axis(), in.get_y_axis());
  std::cout << "  -> ok" << std::endl;
}

} // namespace

int main()
{
  auto logger = mylogger::create_logger(
      spdlog::level::info,
      spdlog::level::err,
      spdlog::level::info,
      "logs/test_grid1d_roundtrip.log");
  (void)logger;

  const fs::path dir = fs::temp_directory_path() / "test_grid1d_roundtrip";
  fs::create_directories(dir);

  try {
    test_grid1d_roundtrip_simple(dir);
    test_grid1d_roundtrip_shifted_edges(dir);
    test_grid1d_repeated_roundtrip_no_drift(dir);
    test_grid2d_roundtrip(dir);
  } catch (const std::exception& ex) {
    std::cerr << "Test failed: " << ex.what() << std::endl;
    return 1;
  }

  std::cout << "All Grid1d/Grid2d save-load roundtrip tests passed." << std::endl;
  return 0;
}
