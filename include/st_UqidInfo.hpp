/// @file st_UqidInfo.hpp
/// @brief Unique detector element ID information structure and related utilities
///
/// @details
/// This header defines UqidInfo, a Plain Old Data (POD) structure that maps globally
/// unique detector element IDs (Uqid) to detector-specific identifiers and group
/// memberships. It also provides hash and equality functors for use with STL containers,
/// as well as stream output and fmt formatting support.
///
/// Key components:
/// - UqidInfo: Core structure holding ID mappings
/// - UqidInfoNotFound: Sentinel constant for failed lookups
/// - AllowUqidInfoPartialDefault: Predicate for partial initialization validation
/// - UqidInfoHash/UqidInfoEq: Functors for hash-based containers
/// - Stream and fmt formatting support
///
/// @note Thread-safety: UqidInfo is POD and thread-safe for concurrent reads.
///       Concurrent writes require external synchronization.
/// @note Coordinate system: Not applicable (ID mapping structure)
/// @note Units: Dimensionless integer indices
#pragma once
#include <array>
#include <functional>
#include <cstddef>  // for size_t
#include <stdexcept>
#include <ostream>
#include <fmt/format.h>
#include <fmt/std.h>   // Enable fmt support for std::filesystem::path
#include <fmt/ostream.h>  // for fmt::ostream
#include "ns_type_definitions.hpp"  // for Uqid, Detid, Ixiy,

using namespace index_type_definitions;

/// @brief Unique ID information for all detector elements
///
/// @details
/// UqidInfo maps unique detector element IDs (Uqid) to detector-specific identifiers
/// and group memberships. Each detector element has a globally unique ID (uqid) that
/// maps to its detector ID (detid), spatial index (ixiy), and various group assignments.
///
/// @note Thread-safety: This structure is Plain Old Data (POD) and inherently thread-safe
///       for read operations. Concurrent writes require external synchronization.
struct UqidInfo {
  Uqid uqid = UqidNotAssigned;  ///< Globally unique detector element ID
  Detid detid = DetidNotAssigned;  ///< Detector identifier
  Ixiy ixiy = IxiyNotAssigned;  ///< Spatial index (ix, iy) within detector
  Inthis inthis = InthisNotAssigned;  ///< Local element ID within this detector
  Igroup igroup = IgroupNotAssigned;  ///< Group identifier
  Uqig uqig = UqigNotAssigned;  ///< Unique ID of the group this element belongs to
  UqigAvail uqig_avail = UqigAvailNotAssigned; ///< Unique ID of the analysis-available group this element belongs to
  bool is_avail = false; ///< Whether this element is available for analysis

  /// @brief Equality operator
  friend bool operator==(const UqidInfo& lhs, const UqidInfo& rhs) {
    return  lhs.uqid == rhs.uqid
         && lhs.detid == rhs.detid
         && lhs.ixiy == rhs.ixiy
         && lhs.inthis == rhs.inthis
         && lhs.igroup == rhs.igroup
         && lhs.uqig == rhs.uqig
         && lhs.uqig_avail == rhs.uqig_avail
         && lhs.is_avail == rhs.is_avail;
  }

  /// @brief Inequality operator
  friend bool operator!=(const UqidInfo& lhs, const UqidInfo& rhs) {
    return !(lhs == rhs);
  }
};

/// @brief Default "not found" sentinel value for UqidInfo
///
/// @details
/// This constant represents a detector element that was not found in a lookup operation.
/// All ID fields are set to their respective "NotFound" sentinel values, and is_avail
/// is false.
///
/// @note This is a compile-time constant that can be used for comparison or as a
///       return value when an element lookup fails.
inline static constexpr UqidInfo UqidInfoNotFound{
    UqidNotFound          // uqid
  , DetidNotFound         // detid
  , IxiyNotFound          // ixiy
  , InthisNotFound        // inthis
  , IgroupNotFound        // igroup
  , UqigNotFound          // uqig
  , UqigAvailNotFound     // uqig_avail
  , false                 // is_avail
};


/// @brief Predicate functor to check if specific UqidInfo fields are at default values
///
/// @details
/// Returns true if the group-related fields (igroup, uqig_avail, is_avail) are
/// unassigned. This allows partial initialization where core ID fields (uqid, detid,
/// ixiy, inthis) are set but group memberships are deferred.
///
/// @note Originally intended for UqidManager::insertInit to validate partially initialized
///       entries. May be redundant depending on actual usage patterns.
struct AllowUqidInfoPartialDefault {
  bool operator()(const UqidInfo& info) const noexcept {
    return info.igroup     == IgroupNotAssigned
        && info.uqig_avail == UqigAvailNotAssigned
        && info.is_avail   == false;
  }
};

/// @brief Hash functor for UqidInfo
///
/// @details
/// Computes hash value using only the uqid field as the primary key. This is
/// appropriate for hash-based containers (std::unordered_map, std::unordered_set)
/// where uqid uniquely identifies each detector element.
///
/// @note Only hashes the uqid field; other fields are ignored for hashing purposes.
struct UqidInfoHash {
  size_t operator()(const UqidInfo& info) const noexcept {
    // Hash only uqid as the primary key
    return std::hash<Uqid>()(info.uqid);
  }
};

/// @brief Equality functor for UqidInfo
///
/// @details
/// Delegates to UqidInfo's operator== for field-by-field comparison. Used in
/// hash-based containers to resolve hash collisions.
///
/// @note Compares all fields of UqidInfo, not just uqid.
struct UqidInfoEq {
  bool operator()(const UqidInfo& lhs, const UqidInfo& rhs) const noexcept {
    // Delegate to UqidInfo's operator==
    return lhs == rhs;
  }
};

/// @brief Stream insertion operator for UqidInfo
///
/// @param[in,out] os Output stream
/// @param[in] info UqidInfo instance to output
/// @return Reference to the output stream (for chaining)
///
/// @details
/// Formats UqidInfo as a brace-enclosed, comma-separated list of field=value pairs.
/// The is_avail field is printed as "true" or "false" (not 0/1).
///
/// @note Output format: "{ uqid=X, detid=Y, ixiy=Z, inthis=W, igroup=A, uqig=B, uqig_avail=C, is_avail=D }"
inline std::ostream& operator<<(std::ostream& os, const UqidInfo& info) {
  return os << "{ uqid=" << info.uqid
            << ", detid=" << info.detid
            << ", ixiy=" << info.ixiy
            << ", inthis=" << info.inthis
            << ", igroup=" << info.igroup
            << ", uqig=" << info.uqig
            << ", uqig_avail=" << info.uqig_avail
            << ", is_avail=" << std::boolalpha << info.is_avail
            << " }";
};

/// @brief fmt formatter specialization for UqidInfo
///
/// @details
/// Allows UqidInfo to be formatted with fmt::format() and fmt::print().
/// Formatting is delegated to the existing operator<< by inheriting
/// from fmt::ostream_formatter.
///
/// @note
/// fmt does not automatically use operator<< for user-defined types;
/// this specialization is required to avoid compile-time errors.
template <>
struct fmt::formatter<UqidInfo> : fmt::ostream_formatter {};
