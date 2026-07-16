// ns_stats_util.cpp
#include "ns_stats_util.hpp"
#include <algorithm>
#include <cmath>
#include <random>
#include <omp.h>
#include "ns_mylogger.hpp"
#include "ns_mymacro.hpp"
#include "ns_seed.hpp"

//-------------------------------------
// Internal utilities
//-------------------------------------
namespace {
  constexpr double INV_SQRT_2 = 0.7071067811865475;  // 1.0 / std::sqrt(2.0)

  inline void check_alpha(double alpha) {
    if (alpha <= 0.0 || alpha >= 1.0) {
      THROW_ERROR("check_alpha: alpha must be within 0.0 < alpha < 1.0, got alpha={}", alpha);
    }
  }
  inline void check_nonneg(double x, const char* name) {
    if (x < 0.0) {
      THROW_ERROR("check_nonneg: {} must be non-negative, got {}={}", name, name, x);
    }
  }
} // anonymous

//-------------------------------------
// Public functions
//-------------------------------------

double stats_util::normal_two_sided_p_value(double z_abs)
{
  if (std::isnan(z_abs) || z_abs < 0.0) {
    LOG_WARN("normal_two_sided_p_value: invalid z_abs={}", z_abs);
    return 1.0;
  }
  return std::erfc(z_abs * INV_SQRT_2);
}

std::string stats_util::to_string(const Distribution& d)
{
  switch (d) {
    case Distribution::Gaussian:  return "Gaussian";
    case Distribution::Poisson:   return "Poisson";
    case Distribution::ChiSquare: return "ChiSquare";
    default:                      return "Unknown";
  }
}

bool stats_util::isStatisticallySignificant(
  double value0, double value1, Distribution dist, double alpha)
{
  check_alpha(alpha);
  check_nonneg(value0, "value0");
  check_nonneg(value1, "value1");

  const double diff = std::fabs(value0 - value1);
  double z_abs = 0.0;

  switch (dist) {
    case Distribution::Gaussian: {
      // Approximate variance of difference as var0 + var1 (using observed values as variance estimates)
      const double var0 = std::max(value0, 1e-12);
      const double var1 = std::max(value1, 1e-12);
      z_abs = diff / std::sqrt(var0 + var1);
      break;
    }
    case Distribution::Poisson: {
      // Skellam/normal approximation: Var(n1 - n2) ≈ n1 + n2
      const double variance = std::max(value0 + value1, 1e-12);
      z_abs = diff / std::sqrt(variance);
      break;
    }
    case Distribution::ChiSquare: {
      // Simplified approximation: var ≈ 2k where k ≈ value0
      const double denom = std::max(value0, 1e-12);
      z_abs = diff / std::sqrt(2.0 * denom);
      break;
    }
    default: {
      THROW_ERROR("stats_util::isStatisticallySignificant: Unsupported distribution={}", static_cast<int>(dist));
    }
  }

  const double p_value = normal_two_sided_p_value(z_abs);
  return (p_value < alpha);
}

//========================
// Poisson: Skellam/normal approximation
//========================

stats_util::StatResult stats_util::poisson_count_ztest(double n1, double n2, double alpha, bool both_side)
{
  check_alpha(alpha);
  check_nonneg(n1, "n1");
  check_nonneg(n2, "n2");

  const double denom = std::sqrt(std::max(1e-12, n1 + n2));
  if (both_side) {
    // Two-sided test: z = |n2 - n1| / denom, p = erfc(z / sqrt(2))
    const double z = std::fabs(n2 - n1) / denom;
    const double p = normal_two_sided_p_value(z);
    return {z, p, alpha, p < alpha};
  } else {
    // One-sided test: sign preserved, p = erfc(|z| / sqrt(2)) / 2
    const double z = (n2 - n1) / denom;
    const double p = 0.5 * std::erfc(std::fabs(z) * INV_SQRT_2);
    return {z, p, alpha, p < alpha};
  }
}

constexpr double SQRT_2PI   = 2.5066282746310002;

double stats_util::asymmetric_gaussian_pdf(
  const double x, const double mu, const double sigma_low, const double sigma_upp)
{
  double sigma = (x < mu) ? sigma_low : sigma_upp;
  double norm = SQRT_2PI * 0.5 * (sigma_low + sigma_upp); // Normalization constant (simplified form)
  double exponent = -0.5 * (x - mu)*(x - mu) / (sigma*sigma);
  return std::exp(exponent) / norm;
}

double stats_util::sample_asymmetric_gaussian(
  const double val_low, const double val_cnt, const double val_upp)
{
  // Left and right standard deviations
  double sigma_low = val_cnt - val_low;
  double sigma_upp = val_upp - val_cnt;

  // Thread-local RNG initialization (reproducible global_seed + thread_id)
  thread_local std::mt19937 gen = []() {
    int tid = omp_get_thread_num();
    std::seed_seq seq{seed::get_global_seed(), static_cast<unsigned>(tid)};
    return std::mt19937(seq);
  }();

  // Choose left/right with uniform random number
  std::uniform_real_distribution<double> uni(0.0, 1.0);
  double u = uni(gen);

  if (u < sigma_low / (sigma_low + sigma_upp)) {
    // Left side: val_cnt - abs(N(0, sigma_low))
    std::normal_distribution<double> norm(0.0, sigma_low);
    return val_cnt - std::abs(norm(gen));
  } else {
    // Right side: val_cnt + abs(N(0, sigma_upp))
    std::normal_distribution<double> norm(0.0, sigma_upp);
    return val_cnt + std::abs(norm(gen));
  }
}
