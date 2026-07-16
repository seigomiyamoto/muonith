/// @file cls_PrefixSum2D.cpp
/// @brief Implementation of PrefixSum2D — 2D summed-area table
///
/// # What this does (plain language)
///
///   Given a 2D grid of numbers ("raw"), this builds a helper table ("prefix")
///   so that we can instantly answer the question:
///
///     "What is the sum of all values in the rectangle (x0..x1, y0..y1)?"
///
///   Without this table: we must loop over every cell in the rectangle → O(W×H).
///   With    this table: 4 lookups + 3 additions                      → O(1).
///
/// # Concrete numeric example (3×3 grid)
///
///   Suppose we have a 3-column × 3-row grid of raw values:
///
///       raw[iy][ix]:       ix=0  ix=1  ix=2
///                   iy=0 [  1     2     3  ]
///                   iy=1 [  4     5     6  ]
///                   iy=2 [  7     8     9  ]
///
///   The prefix table is one row and one column larger (4×4),
///   with the 0-th row and 0-th column filled with zeros:
///
///       prefix_[iy][ix]:   ix=0  ix=1  ix=2  ix=3
///                   iy=0 [  0     0     0     0  ]   ← sentinel row (all zeros)
///                   iy=1 [  0     1     3     6  ]
///                   iy=2 [  0     5    12    21  ]
///                   iy=3 [  0    12    27    45  ]
///                          ↑
///                     sentinel col
///
///   Each cell prefix_[iy+1][ix+1] = sum of all raw values at (0..iy, 0..ix).
///
///   Example: prefix_[2][3] = 21 = sum of raw[0..1][0..2] = 1+2+3+4+5+6 = 21.
///
/// # How build works (step by step)
///
///   To fill prefix_[iy+1][ix+1], we use three already-computed neighbours:
///
///       prefix_[iy+1][ix+1] = raw[iy][ix]            ... the new cell itself
///                            + prefix_[iy  ][ix+1]    ... sum of all rows above
///                            + prefix_[iy+1][ix  ]    ... sum of all cols to the left
///                            - prefix_[iy  ][ix  ]    ... subtract double-counted corner
///
///   Visual explanation for cell (iy=1, ix=1), i.e. raw value = 5:
///
///       prefix_[2][2] = raw[1][1]          = 5       ... the cell itself
///                      + prefix_[1][2]     = 3       ... sum above  (1+2)
///                      + prefix_[2][1]     = 5       ... sum left   (1+4)
///                      - prefix_[1][1]     = 1       ... corner     (1)
///                      = 5 + 3 + 5 - 1    = 12  ✓   (= 1+2+4+5)
///
/// # How query works
///
///   To get sum of raw values in rectangle [ix_lo..ix_hi] × [iy_lo..iy_hi]:
///
///       result = prefix_[iy_hi+1][ix_hi+1]   ... whole rectangle from origin
///              - prefix_[iy_lo  ][ix_hi+1]   ... subtract rows above the window
///              - prefix_[iy_hi+1][ix_lo  ]   ... subtract cols left of the window
///              + prefix_[iy_lo  ][ix_lo  ]   ... add back doubly-subtracted corner
///
///   Example: sum of raw[1..2][1..2] (= 5+6+8+9 = 28)
///
///       = prefix_[3][3] - prefix_[1][3] - prefix_[3][1] + prefix_[1][1]
///       = 45            - 6             - 12             + 1
///       = 28 ✓

#include "cls_PrefixSum2D.hpp"
#include "cls_DetectorPanel.hpp"   // getDetectorElement, get_nbinx, get_nbiny
#include "ns_mymacro.hpp"          // THROW_ERROR

#include <cassert>

// =====================================================================
//  build_signal / build_noise  —  extract raw grid then delegate
// =====================================================================

void PrefixSum2D::build_signal(const DetectorPanel& panel)
{
  const int nx = panel.get_nbinx();
  const int ny = panel.get_nbiny();

  // --- Extract raw signal values into a plain 2D array [iy][ix] ---
  std::vector<std::vector<double>> raw(ny, std::vector<double>(nx));
  for (int iy = 0; iy < ny; ++iy)
    for (int ix = 0; ix < nx; ++ix)
      raw.at(iy).at(ix) = panel.getDetectorElement(ix, iy).get_signal();

  build_from_raw(raw);
}

void PrefixSum2D::build_noise(const DetectorPanel& panel)
{
  const int nx = panel.get_nbinx();
  const int ny = panel.get_nbiny();

  // --- Extract raw noise values into a plain 2D array [iy][ix] ---
  std::vector<std::vector<double>> raw(ny, std::vector<double>(nx));
  for (int iy = 0; iy < ny; ++iy)
    for (int ix = 0; ix < nx; ++ix)
      raw.at(iy).at(ix) = panel.getDetectorElement(ix, iy).get_noise();

  build_from_raw(raw);
}

// =====================================================================
//  build_from_raw  —  the core prefix-sum construction
// =====================================================================
//
//   Builds a "summed-area table" from the raw 2D grid.
//
//   After this function completes, every cell in the prefix table holds
//   the sum of ALL raw values in the rectangle from the top-left corner
//   (0,0) down to that cell:
//
//       prefix_[iy+1][ix+1] = Σ raw[r][c]   for r in 0..iy, c in 0..ix
//
//   The +1 offset exists because row 0 and column 0 of prefix_ are
//   filled with zeros ("sentinels"), which eliminates boundary checks.
//
//   The recurrence uses three already-filled neighbours:
//
//       ┌──────────────┬──────────┐
//       │  P[iy][ix]   │ P[iy][ix+1]    ← already computed (previous row)
//       │  (corner)    │ (above)  │
//       ├──────────────┼──────────┤
//       │ P[iy+1][ix]  │ iy+1,ix+1│    ← the cell we are computing NOW
//       │  (left)      │ = ?      │
//       └──────────────┴──────────┘
//
//       P[iy+1][ix+1] = raw[iy][ix]        ... the new cell's own value
//                      + P[iy  ][ix+1]     ... everything above (already summed)
//                      + P[iy+1][ix  ]     ... everything left  (already summed)
//                      - P[iy  ][ix  ]     ... undo double-count of upper-left
//
//   Why subtract P[iy][ix]?  Because "above" and "left" both include
//   the upper-left rectangle, so it gets counted twice → subtract once.
//
void PrefixSum2D::build_from_raw(const std::vector<std::vector<double>>& raw)
{
  ny_ = static_cast<int>(raw.size());
  if (ny_ == 0) { nx_ = 0; prefix_.clear(); return; }
  nx_ = static_cast<int>(raw.at(0).size());

  // Allocate (ny_+1) × (nx_+1), all zeros.
  // The extra row/col of zeros act as sentinels:
  //   prefix_[0][*] = 0   →  "sum of zero rows = 0"
  //   prefix_[*][0] = 0   →  "sum of zero columns = 0"
  // This means the recurrence works without any if-checks at the boundary.
  prefix_.assign(ny_ + 1, std::vector<double>(nx_ + 1, 0.0));

  // Scan raw grid from top-left to bottom-right.
  // Each cell only reads neighbours above and to the left, which are
  // already filled, so one pass in row-major order is sufficient.
  for (int iy = 0; iy < ny_; ++iy) {
    for (int ix = 0; ix < nx_; ++ix) {
      prefix_.at(iy + 1).at(ix + 1)
        = raw.at(iy).at(ix)             // this cell's own value
        + prefix_.at(iy    ).at(ix + 1) // sum of all rows above this cell
        + prefix_.at(iy + 1).at(ix    ) // sum of all cols left of this cell
        - prefix_.at(iy    ).at(ix    );// subtract upper-left (counted twice)
    }
  }
}

// =====================================================================
//  query  —  O(1) rectangular range sum
// =====================================================================
//
//   Returns the sum of raw values inside the rectangle
//   [ix_lo..ix_hi] × [iy_lo..iy_hi]  (0-based, inclusive on both ends).
//
//   The idea: prefix_[iy+1][ix+1] already holds the sum from (0,0) to (iy,ix).
//   To extract a sub-rectangle, we use inclusion–exclusion:
//
//       ┌─────────────────────────────────┐
//       │  A = P[iy_lo][ix_lo]            │ B = P[iy_lo][ix_hi+1]
//       │  (upper-left corner)            │ (upper-right corner)
//       ├────────────┬────────────────────┤
//       │            │ ← WANTED AREA →    │
//       │            │ ix_lo..ix_hi        │
//       │            │ iy_lo..iy_hi        │
//       ├────────────┴────────────────────┤
//       │  C = P[iy_hi+1][ix_lo]          │ D = P[iy_hi+1][ix_hi+1]
//       └─────────────────────────────────┘
//
//       WANTED = D - B - C + A
//
//       D = total from origin to bottom-right of the window
//       B = strip above the window  (must subtract)
//       C = strip left of the window (must subtract)
//       A = upper-left corner (subtracted twice by B and C → add back once)
//
//   Example from the file header (raw[1..2][1..2] = 5+6+8+9 = 28):
//       D=45, B=6, C=12, A=1  →  45 - 6 - 12 + 1 = 28 ✓
//
double PrefixSum2D::query(int ix_lo, int ix_hi, int iy_lo, int iy_hi) const
{
  // Clamp indices to valid range so callers don't need boundary checks.
  // If the clamped range is empty (lo > hi), return 0.
  if (ix_lo < 0)   ix_lo = 0;
  if (iy_lo < 0)   iy_lo = 0;
  if (ix_hi >= nx_) ix_hi = nx_ - 1;
  if (iy_hi >= ny_) iy_hi = ny_ - 1;
  if (ix_lo > ix_hi || iy_lo > iy_hi) return 0.0;

  // D - B - C + A
  return prefix_.at(iy_hi + 1).at(ix_hi + 1)
       - prefix_.at(iy_lo    ).at(ix_hi + 1)
       - prefix_.at(iy_hi + 1).at(ix_lo    )
       + prefix_.at(iy_lo    ).at(ix_lo    );
}
