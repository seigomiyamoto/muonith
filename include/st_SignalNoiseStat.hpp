/// @file st_SignalNoiseStat.hpp
/// @brief Signal and noise statistics structures for muon detection analysis
/// @details This module defines data structures and utility functions for analyzing
///          signal-to-noise ratios in muon detector measurements. It provides:
///          - SignalNoiseSum: Stores signal/noise sums with 2D spatial range (tx/ty)
///          - SignalNoiseDepth: Extends SignalNoiseSum with detector ID and depth range
///          - SignalNoiseStatResult: Stores statistical test results with z-values, p-values
///          - Utilities for statistical significance testing using Poisson count z-tests
///          - Output functions for ASCII and CSV formats
///          - Sorting utilities for analysis results
///
/// @note Units: depth values are in meters [m], density differences in kg/m3
/// @note Thread-safety: All output functions are not thread-safe due to shared FILE* usage
#pragma once

#include <vector>
#include <utility> // for std::pair
#include <iterator> // for std::begin(), std::end()
#include <filesystem> // for std::filesystem::path
#include <Eigen/Dense>
#include "cls_AABB.hpp"
#include "ns_type_definitions.hpp"

using namespace index_type_definitions;

/// @brief Stores signal and noise sums with 2D spatial range information
/// @details This structure holds accumulated signal and noise measurements
///          along with the spatial extent (tx/ty range) over which they were collected.
///          tx and ty represent detector-relative coordinates or angular ranges.
struct SignalNoiseSum {
  AABB2d txty_range; ///< 2D spatial range (tx/ty coordinates)
  double signal = 0.0; ///< accumulated signal sum
  double noise = 0.0; ///< accumulated noise sum

  /// @brief Get lower bound of tx range
  /// @return xmin of txty_range
  double tx_lower() const { return txty_range.xmin(); }

  /// @brief Get upper bound of tx range
  /// @return xmax of txty_range
  double tx_upper() const { return txty_range.xmax(); }

  /// @brief Get lower bound of ty range
  /// @return ymin of txty_range
  double ty_lower() const { return txty_range.ymin(); }

  /// @brief Get upper bound of ty range
  /// @return ymax of txty_range
  double ty_upper() const { return txty_range.ymax(); }
};

/// @brief Extends SignalNoiseSum with detector ID and elevation-from-detector range
/// @details This structure adds detector identification and elevation bounds to the
///          basic signal/noise measurements. Elevation from detector represents the
///          vertical distance from the detector to the measurement region in meters.
/// @note Units: elev_from_det_lower and elev_from_det_upper are in meters [m]
struct SignalNoiseDepth {
  Detid detid = DetidNotAssigned; ///< detector identifier (DetidNotAssigned = unassigned)
  SignalNoiseSum sn;  ///< signal/noise sum with tx/ty spatial range
  double elev_from_det_lower = 0.0; ///< lower elevation bound from detector [m]
  double elev_from_det_upper = 0.0; ///< upper elevation bound from detector [m]

  /// @brief Get lower bound of tx range
  /// @return xmin of sn.txty_range
  double tx_lower() const { return sn.tx_lower(); }

  /// @brief Get upper bound of tx range
  /// @return xmax of sn.txty_range
  double tx_upper() const { return sn.tx_upper(); }

  /// @brief Get lower bound of ty range
  /// @return ymin of sn.txty_range
  double ty_lower() const { return sn.ty_lower(); }

  /// @brief Get upper bound of ty range
  /// @return ymax of sn.txty_range
  double ty_upper() const { return sn.ty_upper(); }

  /// @brief Get accumulated signal sum
  /// @return signal value
  double signal() const { return sn.signal; };

  /// @brief Get accumulated noise sum
  /// @return noise value
  double noise() const { return sn.noise; };
};

/// @brief Stores complete statistical test results for signal/noise comparison
/// @details This structure holds measurements for both base and modified cases
///          along with statistical test results (z-test, p-value). Used for
///          analyzing the significance of signal differences when a structure is
///          introduced in the measurement region.
/// @note Units: obj_size and z_value in meters [m], delta_dens in kg/m3
struct SignalNoiseStatResult{
  SignalNoiseDepth snd_base; ///< base case measurements (without structure)
  SignalNoiseDepth snd_modi; ///< modified case measurements (with structure)
  double obj_size = 0.0; ///< size of the introduced object [m]
  double delta_dens = 0.0; ///< density difference of the object [kg/m3]
  double z_value = 0.0; ///< z-statistic from hypothesis test
  double p_value = 0.0; ///< p-value from hypothesis test
  double stat_alpha = 0.0; ///< significance level threshold (e.g., 0.05)
  double signal_amp = 1.0; ///< signal amplification factor applied
  bool is_significant = false; ///< true if p_value < stat_alpha

  /// @brief Get detector ID from base case
  /// @return detector identifier
  Detid detid() const { return snd_base.detid; };

  /// @brief Get lower bound of tx range (base case)
  /// @return xmin of base txty_range
  double tx_lower() const { return snd_base.tx_lower(); }

  /// @brief Get upper bound of tx range (base case)
  /// @return xmax of base txty_range
  double tx_upper() const { return snd_base.tx_upper(); }

  /// @brief Get lower bound of ty range (base case)
  /// @return ymin of base txty_range
  double ty_lower() const { return snd_base.ty_lower(); }

  /// @brief Get upper bound of ty range (base case)
  /// @return ymax of base txty_range
  double ty_upper() const { return snd_base.ty_upper(); }

  /// @brief Get lower elevation-from-detector bound (base case)
  /// @return elev_from_det_lower in meters [m]
  double elev_from_det_lower() const { return snd_base.elev_from_det_lower; };

  /// @brief Get upper elevation-from-detector bound (base case)
  /// @return elev_from_det_upper in meters [m]
  double elev_from_det_upper() const { return snd_base.elev_from_det_upper; };

  /// @brief Get signal sum for base case
  /// @return base signal value
  double signal_base() const { return snd_base.signal(); };

  /// @brief Get signal sum for modified case
  /// @return modified signal value
  double signal_modi() const { return snd_modi.signal(); };

  /// @brief Get noise sum for base case
  /// @return base noise value
  double noise_base() const { return snd_base.noise(); };

  /// @brief Get noise sum for modified case
  /// @return modified noise value
  double noise_modi() const { return snd_modi.noise(); };
};

/// @brief Sort keys for SignalNoiseStatResult
/// @details Enumeration of all sortable fields in SignalNoiseStatResult.
///          Used with SignalNoiseStatSortCondition to specify multi-level sorting.
enum class SignalNoiseStatSortKey {
  Detid          ///< Sort by detector ID
, Tx_lower       ///< Sort by lower tx bound
, Tx_upper       ///< Sort by upper tx bound
, Ty_lower       ///< Sort by lower ty bound
, Ty_upper       ///< Sort by upper ty bound
, Depth_lower    ///< Sort by lower depth bound
, Depth_upper    ///< Sort by upper depth bound
, SignalBase     ///< Sort by base case signal
, SignalModi     ///< Sort by modified case signal
, ObjSize        ///< Sort by object size
, DeltaDens      ///< Sort by density difference
, ZValue         ///< Sort by z-statistic
, PValue         ///< Sort by p-value
, StatAlpha      ///< Sort by significance level
, SignalAmp      ///< Sort by signal amplification factor
, IsSignificant  ///< Sort by significance flag
};

/// @brief Sorting condition specifying key and order
/// @details Defines a single sorting criterion with a field key and direction.
///          Multiple conditions can be combined for multi-level sorting.
struct SignalNoiseStatSortCondition {
  SignalNoiseStatSortKey key; ///< field to sort by
  bool ascending;             ///< true for ascending, false for descending
};

/// @brief Namespace for signal/noise statistical analysis utilities
/// @details Provides functions for:
///          - Statistical significance testing (Poisson z-test)
///          - Output formatting (ASCII, CSV)
///          - Sorting and filtering results
namespace signal_noise_stat_util {

  /// @brief Output std::vector<SignalNoiseDepth> as ASCII format
  /// @param vec_snd vector of signal/noise/depth data to output
  /// @param fout output file pointer (default: stdout)
  /// @note Not thread-safe due to shared FILE* usage
  void out_signal_noise_depth_info_vector(
    const std::vector<SignalNoiseDepth> &vec_snd, FILE *fout = stdout );

  /// @brief Output std::vector<SignalNoiseDepth> as ASCII format to specified path
  /// @param vec_snd vector of signal/noise/depth data to output
  /// @param path_out output file path
  /// @throws std::runtime_error if file cannot be opened
  void out_signal_noise_depth_info_vector(
      const std::vector<SignalNoiseDepth> &vec_snd
    , const std::filesystem::path &path_out );

  /// @brief Evaluate statistical significance of signal differences and store results
  /// @details Performs Poisson count z-tests comparing base and modified cases
  ///          for each detector position. Appends results to vec_stat_result.
  /// @param vec_snd_base signal/noise sum info vector for base case (without structure)
  /// @param vec_snd_modi signal/noise sum info vector for modified case (with structure)
  /// @param obj_size size of introduced object [m]
  /// @param delta_dens density difference of object [kg/m3]
  /// @param stat_alpha significance level threshold (e.g., 0.05)
  /// @param signal_amp signal amplification factor applied to counts
  /// @param vec_stat_result output vector where results are appended
  /// @throws std::runtime_error if vec_snd_base.size() != vec_snd_modi.size()
  /// @note Uses stats_util::poisson_count_ztest internally
  void eval_signal_significance(
      const std::vector<SignalNoiseDepth>& vec_snd_base
    , const std::vector<SignalNoiseDepth>& vec_snd_modi
    , const double obj_size
    , const double delta_dens
    , const double stat_alpha
    , const double signal_amp
    , std::vector<SignalNoiseStatResult>& vec_stat_result
    , const bool both_side = false );

  /// @brief Output std::vector<SignalNoiseStatResult> as ASCII format with header
  /// @param vec_stat_result vector to output
  /// @param fout output file pointer (default: stdout)
  /// @note Skips entries with negative elev_from_det_lower
  /// @note Not thread-safe due to shared FILE* usage
  void out_signal_stat_result_vector(
      const std::vector<SignalNoiseStatResult>& vec_stat_result, FILE *fout = stdout );

  /// @brief Output std::vector<SignalNoiseStatResult> as ASCII format to specified path
  /// @param vec_stat_result vector to output
  /// @param outpath output file path
  /// @throws std::runtime_error if file cannot be opened
  void out_signal_stat_result_vector(
      const std::vector<SignalNoiseStatResult>& vec_stat_result
    , const std::filesystem::path& outpath );

  /// @brief Output std::vector<SignalNoiseStatResult> as CSV format with header
  /// @param vec_stat_result vector to output
  /// @param fout output file pointer (if nullptr, defaults to stdout)
  /// @note Skips entries with negative elev_from_det_lower
  /// @note Not thread-safe due to shared FILE* usage
  void out_signal_stat_result_vector_csv(
    const std::vector<SignalNoiseStatResult>& vec_stat_result, FILE *fout);

  /// @brief Output std::vector<SignalNoiseStatResult> as CSV format to specified path
  /// @param vec_stat_result vector to output
  /// @param outpath output file path
  /// @throws std::runtime_error if file cannot be opened
  void out_signal_stat_result_vector_csv(
      const std::vector<SignalNoiseStatResult>& vec_stat_result
    , const std::filesystem::path& outpath );

  /// @brief Evaluate statistical significance and output results to file
  /// @param vec_snd_base base case signal/noise/depth data
  /// @param vec_snd_modi modified case signal/noise/depth data
  /// @param stat_alpha significance level threshold
  /// @param fout output file pointer
  /// @throws std::runtime_error if vec_snd_base.size() != vec_snd_modi.size()
  /// @note Not thread-safe due to shared FILE* usage
  void out_signal_sum_significance_vector(
      const std::vector<SignalNoiseDepth> &vec_snd_base
    , const std::vector<SignalNoiseDepth> &vec_snd_modi
    , const double stat_alpha, FILE *fout );

  /// @brief Evaluate statistical significance and output results to specified path
  /// @param vec_snd_base base case signal/noise/depth data
  /// @param vec_snd_modi modified case signal/noise/depth data
  /// @param stat_alpha significance level threshold
  /// @param path_out output file path
  /// @throws std::runtime_error if file cannot be opened or vector size mismatch
  void out_signal_sum_significance_vector(
      const std::vector<SignalNoiseDepth> &vec_snd_base
    , const std::vector<SignalNoiseDepth> &vec_snd_modi
    , const double stat_alpha
    , const std::filesystem::path &path_out );

  /// @brief Sort SignalNoiseStatResult vector based on multiple conditions
  /// @details Performs stable multi-level sorting using std::sort with a lambda comparator.
  ///          Conditions are applied in order: first condition has highest priority.
  ///          For boolean fields (IsSignificant), false < true numerically.
  /// @param vec_result results to be sorted (modified in-place)
  /// @param vec_condition sorting conditions in priority order
  /// @note Complexity: O(N log N) comparisons
  /// @details Usage example:
  /// @code{.cpp}
  ///  std::vector<SignalNoiseStatResult> results;
  ///  std::vector<SignalNoiseStatSortCondition> conds = {
  ///    {SignalNoiseStatSortKey::Detid, false},  // 1st priority: significant (true) first
  ///    {SignalNoiseStatSortKey::PValue, true},  // 2nd priority: ascending p-value
  ///    {SignalNoiseStatSortKey::ObjSize, true}  // 3rd priority: ascending object size
  ///  };
  ///  signal_noise_stat_util::sort_signal_stat_results(results, conds);
  /// @endcode
  void sort_signal_stat_results(
      std::vector<SignalNoiseStatResult>& vec_result
    , const std::vector<SignalNoiseStatSortCondition>& vec_condition);

} // namespace signal_noise_stat_util