/// @file st_GroupInfo.hpp
/// @brief Grouped angular element information structure
/// @details
/// Defines GroupInfo and related utilities for managing metadata about grouped detector angular elements.
/// This header provides a lightweight data structure for representing detector groups along with
/// hash/equality functors for use in standard containers.
///
/// Thread-safety: All types are trivially copyable or stateless functors (thread-safe).
#pragma once
#include <array>
#include <functional>
#include <cstddef>  // for size_t
#include <stdexcept>
#include <fmt/format.h>
#include <fmt/ostream.h>  // for fmt::ostream
#include "ns_type_definitions.hpp"  // for Uqid, Detid, Ixiy,

using namespace index_type_definitions;

/// @brief Metadata for a grouped set of detector angular elements
///
/// @details
/// GroupInfo aggregates identification and availability metadata for a detector group.
/// Each group represents a logical collection of angular detector elements that share
/// common processing or analysis characteristics.
///
/// Invariants:
/// - If is_avail == true, then uqig_avail must hold a valid assigned ID (not NotAssigned).
/// - If is_avail == false, uqig_avail should be UqigAvailNotAssigned.
///
/// Thread-safety: Trivially copyable, safe to copy across threads.
struct GroupInfo {
  Uqig uqig = UqigNotAssigned; ///< Unique identifier for this group (global scope)
  UqigAvail uqig_avail = UqigAvailNotAssigned; ///< Unique ID for analysis-available groups (subset of all groups)
  Detid detid = DetidNotAssigned; ///< Detector ID to which this group belongs
  Igroup igroup = IgroupNotAssigned; ///< Group index within the detector (local scope)
  bool is_avail = false; ///< Whether this group is available for analysis

  /// @brief Equality comparison operator
  friend bool operator==(const GroupInfo& lhs, const GroupInfo& rhs) {
    return lhs.uqig == rhs.uqig &&
           lhs.uqig_avail == rhs.uqig_avail &&
           lhs.detid == rhs.detid &&
           lhs.igroup == rhs.igroup &&
           lhs.is_avail == rhs.is_avail;
  };

  /// @brief Inequality comparison operator
  friend bool operator!=(const GroupInfo& lhs, const GroupInfo& rhs) {
    return !(lhs == rhs);
  };
};

/// @brief Default "not found" sentinel value for GroupInfo
///
/// @details
/// This constant represents a GroupInfo instance that indicates "no group found" or "invalid group".
/// All ID fields are set to their respective NotFound sentinel values, and is_avail is false.
inline static constexpr GroupInfo GroupInfoNotFound{
  UqigNotFound         // uqig
, UqigAvailNotFound    // uqig_avail
, DetidNotFound        // detid
, IgroupNotFound       // igroup
, false                // is_avail
};

/// @brief Hash functor for GroupInfo (for use in std::unordered_map, std::unordered_set)
///
/// @details
/// Hashes GroupInfo based solely on the uqig field, which is the global unique identifier.
/// This assumes uqig uniquely identifies a group; collisions in other fields are considered
/// equivalent if uqig matches.
///
/// @note Thread-safety: Stateless functor, safe to use concurrently.
struct GroupInfoHash {
  size_t operator()(const GroupInfo& g) const noexcept {
    return std::hash<Uqig>()(g.uqig);
  }
};

/// @brief Equality functor for GroupInfo (for use in std::unordered_map, std::unordered_set)
///
/// @details
/// Delegates to GroupInfo::operator==, which compares all fields.
///
/// @note Thread-safety: Stateless functor, safe to use concurrently.
struct GroupInfoEq {
  bool operator()(const GroupInfo& a, const GroupInfo& b) const noexcept {
    return a == b;
  }
};

/// @brief Hash functor for (Detid, Igroup) pairs (for use in std::unordered_map)
///
/// @details
/// Combines hashes of Detid and Igroup using XOR and bit-shifting.
/// This is a simple mixing strategy suitable for typical use cases.
///
/// @note Thread-safety: Stateless functor, safe to use concurrently.
struct DetidIgroupPairHash {
  size_t operator()(const std::pair<Detid, Igroup>& p) const noexcept {
    // Mix Detid and Igroup hashes via XOR and shift (common simple hash combination).
    return std::hash<Detid>()(p.first) ^ (std::hash<Igroup>()(p.second) << 1);
  }
};

/// @brief Equality functor for (Detid, Igroup) pairs (for use in std::unordered_map)
///
/// @details
/// Compares both components of the pair for equality.
///
/// @note Thread-safety: Stateless functor, safe to use concurrently.
struct DetidIgroupPairEqual {
  bool operator()(const std::pair<Detid, Igroup>& a, const std::pair<Detid, Igroup>& b) const noexcept {
    return (a.first == b.first) && (a.second == b.second);
  }
};

/// @brief Stream insertion operator for GroupInfo
///
/// @param os Output stream
/// @param g GroupInfo instance to output
/// @return Reference to the output stream (for chaining)
///
/// @details
/// Outputs all fields in a structured format: { uqig=..., uqig_avail=..., detid=..., igroup=..., is_avail=... }
inline std::ostream& operator<<(std::ostream& os, const GroupInfo& g) {
  return os << "{ uqig=" << g.uqig
            << ", uqig_avail=" << g.uqig_avail
            << ", detid=" << g.detid
            << ", igroup=" << g.igroup
            << ", is_avail=" << std::boolalpha << g.is_avail << " }";
};

/// @brief fmt formatter specialization for GroupInfo
///
/// @details
/// Allows GroupInfo to be formatted with fmt::format() and fmt::print().
/// Formatting is delegated to the existing operator<< by inheriting
/// from fmt::ostream_formatter.
///
/// @note
/// fmt does not automatically use operator<< for user-defined types;
/// this specialization is required to avoid compile-time errors.
template <>
struct fmt::formatter<GroupInfo> : fmt::ostream_formatter {};
