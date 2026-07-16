/// @file cls_PrefixSum2D.hpp
/// @brief 2D prefix sum (summed-area table) for O(1) rectangular range queries
/// @details
/// This class builds a 2D cumulative sum table from a DetectorPanel's signal
/// or noise values.  Once built, any axis-aligned rectangular sub-region sum
/// can be answered in O(1) using the inclusion–exclusion identity:
///
///   sum(x1..x2, y1..y2) = P[y2][x2] - P[y1-1][x2]
///                        - P[y2][x1-1] + P[y1-1][x1-1]
///
/// Designed to replace the repeated O(M) summation inside
/// DetectorPanel::mp_get_signal_noise_sum_range_xy() when called many times
/// with different rectangular windows (e.g. the y-sweep loop).
///
/// Memory layout: prefix_[iy][ix], 1-indexed internally (0-th row/col = 0).
///
/// @note Thread-safety: Build is single-threaded; query is thread-safe (const).
/// @note Units: values are raw signal or noise counts (double).
#pragma once

#include <vector>

class DetectorPanel;   // forward declaration

/// @brief 2D prefix sum table for O(1) rectangular range sum queries
///
/// Build once per DetectorPanel state, then query arbitrary sub-rectangles.
///
/// Minimal usage example:
/// @code
///   PrefixSum2D ps;
///   ps.build_signal(panel);                // O(nbinx * nbiny)
///   double s = ps.query(ix0, ix1, iy0, iy1); // O(1)
/// @endcode
class PrefixSum2D {
public:
  PrefixSum2D() = default;

  // -----------------------------------------------------------------
  //  Build methods  (call exactly one before querying)
  // -----------------------------------------------------------------

  /// @brief Build prefix sum from DetectorPanel's **signal** values
  /// @param[in] panel  Source detector panel (must be fully populated)
  /// @note Complexity: O(nbinx * nbiny)
  void build_signal(const DetectorPanel& panel);

  /// @brief Build prefix sum from DetectorPanel's **noise** values
  /// @param[in] panel  Source detector panel (must be fully populated)
  /// @note Complexity: O(nbinx * nbiny)
  void build_noise(const DetectorPanel& panel);

  // -----------------------------------------------------------------
  //  Query
  // -----------------------------------------------------------------

  /// @brief Sum of values in the closed rectangle [ix_lo..ix_hi, iy_lo..iy_hi]
  /// @param[in] ix_lo  Left   column index (0-based, inclusive)
  /// @param[in] ix_hi  Right  column index (0-based, inclusive)
  /// @param[in] iy_lo  Bottom row    index (0-based, inclusive)
  /// @param[in] iy_hi  Top    row    index (0-based, inclusive)
  /// @return Sum of stored values within the rectangle
  /// @note Complexity: O(1)
  /// @pre  build_signal() or build_noise() must have been called
  double query(int ix_lo, int ix_hi, int iy_lo, int iy_hi) const;

  // -----------------------------------------------------------------
  //  Accessors
  // -----------------------------------------------------------------

  int nx() const { return nx_; }
  int ny() const { return ny_; }
  bool is_built() const { return !prefix_.empty(); }

private:
  /// @brief Internal: allocate and fill prefix table from raw 2D values
  /// @param[in] raw  2D array [iy][ix] of source values  (0-indexed, size ny x nx)
  void build_from_raw(const std::vector<std::vector<double>>& raw);

  int nx_ = 0;  ///< number of x bins (columns)
  int ny_ = 0;  ///< number of y bins (rows)

  /// @brief 1-indexed prefix table, size (ny_+1) x (nx_+1).
  ///        prefix_[0][*] = prefix_[*][0] = 0  (sentinel row/col).
  std::vector<std::vector<double>> prefix_;
};
