/// @file ns_type_definitions.hpp
/// @brief Global type definitions and constants for detector indexing
/// @details
/// This header defines the fundamental index types, IDs, and sentinel values (NotFound, NotAssigned)
/// used throughout the detector geometry and tracking codebase.
///
/// Key types defined:
/// - Detid: Detector ID
/// - Uqid: Unique index for DetectorElement
/// - Inthis: ID within a detector
/// - Igroup: Group index
/// - UqigAvail: Unique index for available groups
/// - Composite types: DetIxiy (det_id, ix, iy), DetInthis, DetIgroup
///
/// Thread-safety: All types are POD (Plain Old Data) types and constexpr values; safe for concurrent read access.
/// No mutable global state. Hash functors and comparators are stateless and thread-safe.
///
/// Usage pattern: Include this header wherever detector indexing types are needed.
#pragma once

#include <string>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <set>

#include <chrono>
#include <thread>

#define _USE_MATH_DEFINES // Must be defined before including <cmath> to access M_PI

// std::map<>
#include <map>
#include <unordered_map>
#include <algorithm>

#include <cassert>
#include <vector>
#include <filesystem> // for std::filesystem::path

//#########################################################################
//#########################################################################
/// @namespace index_type_definitions
/// @brief Type definitions for indices (IDs) related to Detector Elements
///
/// @details
/// This namespace provides a comprehensive set of type aliases and sentinel values
/// for detector element indexing. All types use `int` as the underlying type.
///
/// Sentinel values:
/// - `IntNotFound = -1`: Indicates an index was not found during lookup
/// - `IntNotAssigned = -999`: Indicates an index has not yet been assigned
///
/// Each type (Detid, Uqid, Inthis, etc.) has corresponding NotFound and NotAssigned constants.
///
/// Hash and comparison functors are provided for composite types to enable their use
/// in `std::unordered_map` and `std::set` containers.
//#########################################################################
//#########################################################################
namespace index_type_definitions {

  /// @brief int value for NotFound
  static constexpr int IntNotFound = -1;

  /// @brief int value for NotAssigned
  static constexpr int IntNotAssigned = -999;

  /// @brief pair of (ix, iy) of DetectorElement
  using Ixiy = std::array<int,2>;

  /// @brief Ixiy value when not found
  static constexpr Ixiy IxiyNotFound = { IntNotFound, IntNotFound };

  /// @brief Ixiy value when not assigned
  static constexpr Ixiy IxiyNotAssigned = { IntNotAssigned, IntNotAssigned };

  /// @brief Comparison function with iy-major ordering
  /// @details
  /// Compares std::array<int,2> (typically Ixiy) with iy (index [1]) as primary key,
  /// and ix (index [0]) as secondary key.
  /// Use case: Sorting detector elements by row-first (iy-major) order.
  /// Thread-safety: Yes (stateless, const operator).
  struct IyMajorCompare {
    bool operator()(const std::array<int,2>& a, const std::array<int,2>& b) const {
      if (a[1] != b[1]) return a[1] < b[1];  // iy priority
      return a[0] < b[0];                   // if iy equal, compare by ix
    }
  };

  /// @brief Detector ID
  using Detid = int;

  /// @brief Detid value when not found
  static constexpr Detid DetidNotFound = IntNotFound;

  /// @brief Detid value when not assigned
  static constexpr Detid DetidNotAssigned = IntNotAssigned;

  /// @brief ID in a detector
  using Inthis = int;

  /// @brief Inthis value when not found
  static constexpr Inthis InthisNotFound = IntNotFound;

  /// @brief Inthis value when not assigned
  static constexpr Inthis InthisNotAssigned = IntNotAssigned;

  /// @brief unique_index for DetectorElement
  using Uqid = int;

  /// @brief Uqid value when not found
  static constexpr Uqid UqidNotFound = IntNotFound;

  /// @brief Uqid value when not assigned
  static constexpr Uqid UqidNotAssigned = IntNotAssigned;

  /// @brief igroup
  using Igroup = int;

  /// @brief Igroup value when not found
  static constexpr Igroup IgroupNotFound = IntNotFound;

  /// @brief Igroup value when not assigned
  static constexpr Igroup IgroupNotAssigned = IntNotAssigned;

  /// @brief unique_index for Group ID
  using Uqig = int;

  /// @brief Uqig value when not found
  static constexpr Uqig UqigNotFound = IntNotFound;

  /// @brief Uqig value when not assigned
  static constexpr Uqig UqigNotAssigned = IntNotAssigned;

  /// @brief unique_index for available Group ID
  /// @note Uqig and UqigAvail are not necessarily equal
  using UqigAvail = int;

  /// @brief UqigAvail value when not found
  static constexpr UqigAvail UqigAvailNotFound = IntNotFound;

  /// @brief UqigAvail value when not assigned
  static constexpr UqigAvail UqigAvailNotAssigned = IntNotAssigned;

  /// @brief new unique_index for available Group ID when we remove UqigAvail by disabled detid
  using UqigAvailDisabled = int;

  /// @brief UqigAvailDisabled value when not found
  static constexpr UqigAvailDisabled UqigAvailDisabledNotFound = IntNotFound;

  /// @brief UqigAvailDisabled value when not assigned
  static constexpr UqigAvailDisabled UqigAvailDisabledNotAssigned = IntNotAssigned;

  /// @brief Array combining (det_id, ix, iy)
  using DetIxiy = std::array<int, 3>;

  /// @brief DetIxiy value when not found
  static constexpr DetIxiy DetIxiyNotFound = { IntNotFound, IntNotFound, IntNotFound };

  /// @brief DetIxiy value when not assigned
  static constexpr DetIxiy DetIxiyNotAssigned = { IntNotAssigned, IntNotAssigned, IntNotAssigned };

  /// @brief Hash function for DetIxiy with Det → Iy → Ix priority
  /// @details
  /// Combines hashes of detector ID, iy, and ix using XOR and bit-shifting.
  /// Priority order: Det (lowest bits) → Iy (middle bits) → Ix (high bits).
  /// Enables use of DetIxiy as key in std::unordered_map.
  /// Thread-safety: Yes (stateless, const operator).
  struct DetIxiyHash {
    size_t operator()(const DetIxiy& x) const {
      size_t h_det = std::hash<int>()(x[0]);  // Det
      size_t h_iy  = std::hash<int>()(x[2]);  // Iy
      size_t h_ix  = std::hash<int>()(x[1]);  // Ix

      return h_det
          ^ (h_iy << 16)  // Iy to middle bits
          ^ (h_ix << 32); // Ix to high bits
    }
  };

  /// @brief Equality function for DetIxiy
  /// @details
  /// Compares all three components (det_id at [0], ix at [1], iy at [2]).
  /// Required for std::unordered_map with DetIxiy keys.
  /// Thread-safety: Yes (stateless, const operator).
  struct DetIxiyEq {
    bool operator()(const DetIxiy& a, const DetIxiy& b) const {
      return (a[0] == b[0] &&
              a[1] == b[1] &&
              a[2] == b[2]);
    }
  };

  /// @brief Array combining (det_id, id_in_this_detector)
  using DetInthis = std::array<int, 2>;

  /// @brief DetInthis value when not found
  static constexpr DetInthis DetInthisNotFound = { IntNotFound, IntNotFound };

  /// @brief DetInthis value when not assigned
  static constexpr DetInthis DetInthisNotAssigned = { IntNotAssigned, IntNotAssigned };

  /// @brief Array combining (det_id, igroup)
  using DetIgroup = std::array<int, 2>;

  /// @brief DetIgroup value when not found
  static constexpr DetIgroup DetIgroupNotFound = { IntNotFound, IntNotFound };

  /// @brief DetIgroup value when not assigned
  static constexpr DetIgroup DetIgroupNotAssigned = { IntNotAssigned, IntNotAssigned };

  /// @brief Hash function of DetIgroup
  /// @details
  /// Combines hashes of det_id (key[0]) and igroup (key[1]) using XOR and bit-shifting.
  /// Enables use of DetIgroup as key in std::unordered_map.
  /// Thread-safety: Yes (stateless, const operator).
  struct DetIgroupHash {
    std::size_t operator()(const std::array<int, 2>& key) const noexcept {
      std::size_t h1 = std::hash<int>{}(key[0]);
      std::size_t h2 = std::hash<int>{}(key[1]);
      return h1 ^ (h2 << 1);  // Common XOR + shift combination
    }
  };

  /// @brief Equality function of DetIgroup
  /// @details
  /// Compares both components (det_id at [0], igroup at [1]).
  /// Required for std::unordered_map with DetIgroup keys.
  /// Thread-safety: Yes (stateless, const operator).
  struct DetIgroupEq {
    bool operator()(const std::array<int, 2>& lhs, const std::array<int, 2>& rhs) const noexcept {
      return lhs[0] == rhs[0] && lhs[1] == rhs[1];
    }
  };

  /// @brief Hash function for DetInthis
  /// @details
  /// Combines hashes of det_id (x[0]) and id_in_this_detector (x[1]) using XOR and bit-shifting.
  /// Enables use of DetInthis as key in std::unordered_map.
  /// Thread-safety: Yes (stateless, const operator).
  struct DetInthisHash {
    size_t operator()(const DetInthis& x) const {
      size_t h1 = std::hash<int>()(x[0]);
      size_t h2 = std::hash<int>()(x[1]);
      return h1 ^ (h2 << 1);
    }
  };

  /// @brief Equality function for DetInthis
  /// @details
  /// Compares both components (det_id at [0], id_in_this_detector at [1]).
  /// Required for std::unordered_map with DetInthis keys.
  /// Thread-safety: Yes (stateless, const operator).
  struct DetInthisEq {
    bool operator()(const DetInthis& a, const DetInthis& b) const {
      return (a[0] == b[0] &&
              a[1] == b[1]);
    }
  };

  /// @brief Comparison function for SortedDetidUqigSet, sorting Detid and UqigAvail in ascending order
  /// @details
  /// Primary key: Detid (ascending), secondary key: UqigAvail (ascending).
  /// Enables use in std::set for sorted storage of (Detid, UqigAvail) pairs.
  /// Thread-safety: Yes (stateless, const operator).
  struct DetidUqigPairCompare {
    bool operator()(const std::pair<Detid, UqigAvail>& a,
                    const std::pair<Detid, UqigAvail>& b) const {
      if (a.first != b.first)
        return a.first < b.first; // Detid ascending
      return a.second < b.second; // UqigAvail ascending
    }
  };

  /// @brief Index definition
  /// @note Used within UqidManager
  using Index = int;

  /// @brief Index value when not found
  static constexpr Index IndexNotFound = IntNotFound;

  /// @brief Index value when not assigned
  static constexpr Index IndexNotAssigned = IntNotAssigned;

  /// @brief Predicate functor returning true if Index is NotFound
  /// @note May not have been necessary.
  struct AllowIndexNotFound {
    bool operator()(const Index& idx) const noexcept {
      return idx == IndexNotFound;
    }
  };

  /// @brief Set holding pairs of Detid and UqigAvail, sorted in ascending order
  using SortedDetidUqigSet = std::set<std::pair<Detid, UqigAvail>, DetidUqigPairCompare>;

  /// @brief Mapping table for UqigAvailDisabled when any detid is disabled
  using MapUqigAvailDisabled = std::map<UqigAvailDisabled, UqigAvail>;

}; // end of namespace index_type_definitions

/// @brief Stream output operator for std::array<int, 2>
/// @param[in,out] os Output stream
/// @param[in] arr Array to output (e.g., Ixiy, DetInthis, DetIgroup)
/// @return Reference to output stream for chaining
/// @details Formats as "[elem0, elem1]". Useful for logging and debugging.
inline std::ostream& operator<<(std::ostream& os, const std::array<int, 2>& arr) {
  return os << "[" << arr[0] << ", " << arr[1] << "]";
}

/// @brief Stream output operator for std::array<int, 3>
/// @param[in,out] os Output stream
/// @param[in] arr Array to output (e.g., DetIxiy)
/// @return Reference to output stream for chaining
/// @details Formats as "[elem0, elem1, elem2]". Useful for logging and debugging.
inline std::ostream& operator<<(std::ostream& os, const std::array<int, 3>& arr) {
  return os << "[" << arr[0] << ", " << arr[1] << ", " << arr[2] << "]";
}

