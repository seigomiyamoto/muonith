/// @file ns_stats_util.hpp
/// @brief Statistical utility functions for hypothesis testing and distribution calculations
///
/// @details
/// This module provides statistical hypothesis testing utilities, including:
/// - Two-sided p-value calculation from z-scores
/// - Statistical significance evaluation for multiple distribution types
/// - Poisson count comparison (z-test)
///
/// Supported distributions: Gaussian (normal), Poisson, Chi-square
///
/// Thread-safety: All functions are thread-safe (pure computations, no shared state).
///
/// Complexity: All functions have O(1) time complexity (constant time operations).
///
/// Example usage:
/// @code
/// // Compare two Poisson counts
/// auto result = stats_util::poisson_count_ztest(100.0, 120.0, 0.05);
/// if (result.is_significant) {
///   std::cout << "Significant difference: p = " << result.p_value << std::endl;
/// }
///
/// // Check general statistical significance
/// bool sig = stats_util::isStatisticallySignificant(50.0, 70.0,
///                                                    stats_util::Distribution::Poisson, 0.05);
/// @endcode
#pragma once

#include <string>
#include <cmath>

//#########################################################################
//#########################################################################
/// @namespace stats_util
/// @brief Statistical utility functions for hypothesis testing and distribution calculations
//#########################################################################
//#########################################################################
namespace stats_util {

  /// @brief Enumeration of statistical distribution types
  enum class Distribution {
      Gaussian     ///< Gaussian (normal) distribution
    , Poisson      ///< Poisson distribution
    , ChiSquare    ///< Chi-square distribution
  };

  /// @brief Compute two-sided p-value from z-score for standard normal distribution.
  /// @param[in] z_abs Absolute value of the z-score (|z|)
  /// @return Two-sided p-value P(|Z| > z_abs). Returns 1.0 if z_abs is NaN or negative.
  ///
  /// @details
  /// For standard normal distribution \(Z \sim N(0,1)\), the two-sided p-value is:
  /// \f[
  ///   p = P(|Z| > z) = 2 \cdot (1 - \Phi(z)) = \mathrm{erfc}\!\left(\frac{z}{\sqrt{2}}\right)
  /// \f]
  /// where \(\Phi(z)\) is the cumulative distribution function (CDF) of the standard normal,
  /// and \(\mathrm{erfc}(x)\) is the complementary error function.
  ///
  /// The `std::erfc` implementation in C++ standard library `<cmath>` uses polynomial
  /// approximation or asymptotic expansion with O(1) complexity (constant time).
  /// A single call typically takes tens to hundreds of nanoseconds and is rarely
  /// a performance bottleneck in statistical computations.
  ///
  /// @note If z_abs is NaN or negative, logs a warning and returns 1.0 (no significance).
  /// @note Complexity: O(1)
  /// @note Thread-safety: Yes (pure function, no shared state)
  ///
  /// Example: z_abs = 1.96 yields p ≈ 0.05.
  double normal_two_sided_p_value(const double z_abs);

  /// @brief Convert Distribution enum to string representation
  /// @param[in] d Distribution enum value
  /// @return String representation ("Gaussian", "Poisson", "ChiSquare", or "Unknown")
  /// @note Complexity: O(1)
  /// @note Thread-safety: Yes
  std::string to_string(const Distribution& d);

  /// @brief Evaluate statistical significance between two datasets
  /// @param[in] value0 Count/observation data 0 (must be non-negative)
  /// @param[in] value1 Count/observation data 1 (must be non-negative)
  /// @param[in] dist Distribution type (Gaussian, Poisson, ChiSquare)
  /// @param[in] alpha Significance level (must satisfy 0.0 < alpha < 1.0, e.g., 0.05 for 95% confidence)
  /// @return true if statistically significant (p < alpha), false otherwise
  ///
  /// @throws std::runtime_error (via THROW_ERROR) if alpha not in (0,1), value0 < 0, value1 < 0, or unsupported distribution
  ///
  /// @details
  /// Approximation methods by distribution type:
  /// - Gaussian: Variance approximation var(diff) ≈ var0 + var1 where var ≈ observed value
  /// - Poisson: Skellam/normal approximation var(n1 - n2) ≈ n1 + n2
  /// - ChiSquare: Simplified approximation var ≈ 2k where k ≈ value0
  ///
  /// All methods use normal approximation and compute z-score, then two-sided p-value.
  ///
  /// @note Floor values at 1e-12 to avoid division by zero
  /// @note Complexity: O(1)
  /// @note Thread-safety: Yes
  bool isStatisticallySignificant(
    const double value0, const double value1, const Distribution dist, const double alpha);

  /// @brief Structure to hold statistical test results
  struct StatResult {
    double z_value = 0.0;   ///< z-score: test statistic in standard normal distribution (e.g., |z| > 1.96 for 95% significance)
    double p_value = 1.0;   ///< p-value: probability of observing the difference by chance
    double alpha = -1.0;    ///< Significance level (stored from input for debugging purposes)
    bool   is_significant = false; ///< true if statistically significant (p < alpha)
  };

  /// @brief Compare two Poisson counts using z-test. Null hypothesis: the two Poisson processes have equal mean rates.
  /// @param[in] n1 Count 1 (must be non-negative)
  /// @param[in] n2 Count 2 (must be non-negative)
  /// @param[in] alpha Significance level (default: 0.10 means 90% confidence level, must satisfy 0.0 < alpha < 1.0)
  /// @param[in] both_side If true, use two-sided test (p = erfc(|z|/sqrt(2))).
  ///            If false (default), use one-sided test (p = erfc(|z|/sqrt(2)) / 2).
  ///            When one-sided, z preserves the sign of (n2 - n1).
  /// @return StatResult containing z_value, p_value, alpha, and is_significant flag
  ///
  /// @throws std::runtime_error (via THROW_ERROR) if alpha not in (0,1), n1 < 0, or n2 < 0
  ///
  /// @details
  /// Uses Skellam/normal approximation.
  /// - Two-sided (both_side=true):  z = |n2 - n1| / sqrt(n1 + n2), p = erfc(z / sqrt(2))
  /// - One-sided (both_side=false): z = (n2 - n1) / sqrt(n1 + n2),  p = erfc(|z| / sqrt(2)) / 2
  /// Variance is floored at 1e-12 to avoid division by zero when both counts are zero.
  ///
  /// @note Complexity: O(1)
  /// @note Thread-safety: Yes
  StatResult poisson_count_ztest(double n1, double n2, double alpha = 0.10,
                                 bool both_side = false);

  // ======================================================================
  /// @name pdf_gaussian_functions
  /// @details Probability density related functions
  /// @{

  /// @brief PDF of asymmetric Gaussian distribution
  /// @param x Evaluation point
  /// @param mu Center value val_cnt
  /// @param sigma_low Left side σ (val_cnt - val_low)
  /// @param sigma_upp Right side σ (val_upp - val_cnt)
  /// @return Probability density
  double asymmetric_gaussian_pdf(
    const double x, const double mu, const double sigma_low, const double sigma_upp);

  /// @brief Sampling from asymmetric Gaussian distribution
  /// @param val_cnt Center value
  /// @param val_low Left side σ (val_cnt - val_low)
  /// @param val_upp Right side σ (val_upp - val_cnt)
  /// @return Sampled value
  /// @note Function to sample from asymmetric Gaussian distribution with different σ on left and right
  /// @note Uses OpenMP (omp_get_thread_num for thread-local RNG seeding).
  double sample_asymmetric_gaussian(
     const double val_low, const double val_cnt, const double val_upp);

  ///@} ------------------------------------------------------------------

};