/// @file cls_OneToManyUOBimap.hpp
/// @brief Unordered bidirectional one-to-many map implementation
/// @details Header-only template class providing hash-based bidirectional mapping
///          between unique keys (One) and multiple values (Many).
///
/// ## Overview
/// OneToManyUOBimap maintains a bidirectional relationship where:
/// - One key can map to multiple Many values (one-to-many relationship)
/// - Each Many value maps to exactly one One key (many-to-one uniqueness)
///
/// ## Typical Workflow
/// 1. Create an instance: `OneToManyUOBimap<int, std::string> bimap;`
/// 2. Insert mappings: `bimap.insert(1, "a"); bimap.insert(1, "b");`
/// 3. Query by One: `auto vec = bimap.get_vecMany(1);` returns {"a", "b"}
/// 4. Query by Many: `auto one = bimap.getOne("a");` returns 1
/// 5. Erase mappings: `bimap.eraseMany("a");` or `bimap.eraseOne(1);`
///
/// ## Thread Safety
/// This class is NOT thread-safe. External synchronization is required for
/// concurrent access from multiple threads.
///
/// ## Memory Layout
/// Internally uses two hash maps:
/// - `map_One_Many_`: unordered_multimap<One, Many> for one-to-many lookups
/// - `map_Many_One_`: unordered_map<Many, One> for many-to-one lookups
///
/// ## Complexity
/// - insert(): O(1) average
/// - getOne(): O(1) average
/// - get_vecMany(): O(k) where k is the number of Many values for the given One
/// - eraseOne(): O(k) where k is the number of Many values for the given One
/// - eraseMany(): O(n) worst case due to std::erase_if on multimap
#pragma once

#include <unordered_map>
#include <vector>
#include <stdexcept>
#include <set>
#include <fmt/ostream.h>
#include "ns_mylogger.hpp"
#include "ns_mymacro.hpp"

/// @brief Bidirectional one-to-many map using unordered hash containers
/// @tparam One      Key type on the "One" side (unique per mapping set)
/// @tparam Many     Key type on the "Many" side (can have multiple per One)
/// @tparam HashOne  Hash function for One type (default: std::hash<One>)
/// @tparam EqOne    Equality comparator for One type (default: std::equal_to<One>)
/// @tparam HashMany Hash function for Many type (default: std::hash<Many>)
/// @tparam EqMany   Equality comparator for Many type (default: std::equal_to<Many>)
///
/// @note The One type must be an integral type for OneNotAssigned sentinel to work correctly.
///
/// @code
/// // Example usage:
/// OneToManyUOBimap<int, std::string> bimap;
/// bimap.insert(1, "apple");
/// bimap.insert(1, "banana");
/// bimap.insert(2, "cherry");
///
/// auto fruits = bimap.get_vecMany(1);  // {"apple", "banana"}
/// int id = bimap.getOne("cherry");     // 2
/// @endcode
template<
    typename One, typename Many
  , typename HashOne  = std::hash<One>
  , typename EqOne    = std::equal_to<One>
  , typename HashMany = std::hash<Many>
  , typename EqMany   = std::equal_to<Many>
>
class OneToManyUOBimap {
public:
  /// @brief Sentinel value representing an unassigned One (-1)
  /// @note One type must be an integral type for this to work correctly
  static constexpr One OneNotAssigned = static_cast<One>(-1);

private:
  /// @brief Type alias for One->Many multimap (one key can have multiple values)
  using MapOneMany = std::unordered_multimap<One, Many, HashOne, EqOne>;
  /// @brief Type alias for Many->One map (each Many maps to exactly one One)
  using MapManyOne = std::unordered_map<Many, One, HashMany, EqMany>;

  /// @brief Internal storage for One->Many mappings
  MapOneMany map_One_Many_;

  /// @brief Internal storage for Many->One mappings
  MapManyOne map_Many_One_;

public:
  /// @brief Default constructor
  OneToManyUOBimap() = default;

  /// @brief Default destructor
  ~OneToManyUOBimap() = default;

  /// @brief Copy constructor
  OneToManyUOBimap(const OneToManyUOBimap&) = default;

  /// @brief Copy assignment operator
  OneToManyUOBimap& operator=(const OneToManyUOBimap&) = default;

  /// @brief Move constructor
  OneToManyUOBimap(OneToManyUOBimap&&) noexcept = default;

  /// @brief Move assignment operator
  OneToManyUOBimap& operator=(OneToManyUOBimap&&) noexcept = default;

  /// @brief Get const reference to the entire One->Many map
  /// @return Const reference to internal unordered_multimap
  const MapOneMany& getMapOneMany() const { return map_One_Many_; }

  /// @brief Get mutable reference to the entire One->Many map
  /// @return Mutable reference to internal unordered_multimap
  /// @warning Direct modification may break bidirectional consistency
  MapOneMany& callMapOneMany() { return map_One_Many_; }

  /// @brief Get const reference to the entire Many->One map
  /// @return Const reference to internal unordered_map
  const MapManyOne& getMapManyOne() const { return map_Many_One_; }

  /// @brief Get mutable reference to the entire Many->One map
  /// @return Mutable reference to internal unordered_map
  /// @warning Direct modification may break bidirectional consistency
  MapManyOne& callMapManyOne() { return map_Many_One_; }

  /// @brief Insert a One->Many mapping
  /// @param[in] one  The key on the One side (can have multiple Many values)
  /// @param[in] many The key on the Many side (must be unique across all One keys)
  /// @throws std::runtime_error If many is already registered to another One
  /// @note Many->One mapping must be unique; throws if duplicate many is found
  void insert(const One& one, const Many& many) {
    // If many is already registered, the uniqueness constraint of
    // Many->One mapping is violated, so throw an error.
    if (map_Many_One_.contains(many)) {
      THROW_ERROR("OneToManyUOBimap::insert: many={} already registered to another one", many);
    }

    // Add new one->many pair to multimap (One side can hold multiple Many values)
    map_One_Many_.emplace(one, many);

    // Register the unique many->one mapping
    map_Many_One_.emplace(many, one);
  }

  /// @brief Insert a One->Many mapping, overwriting if many already exists
  /// @param[in] one  The key on the One side
  /// @param[in] many The key on the Many side (will be reassigned if exists)
  /// @note If many is already mapped to another One, the old mapping is removed first
  void insertOrOverwrite(const One& one, const Many& many) {
    if (hasMany(many)) eraseMany(many);  // Many must be unique
    insert(one, many);                    // One side can have multiple Many
  }

  /// @brief Get all Many values associated with a given One key
  /// @param[in] one The One key to query
  /// @return Vector of Many values (empty if one is not registered)
  /// @note Logs a warning if one has no associated Many values
  std::vector<Many> get_vecMany(const One& one) const {
    std::vector<Many> result;
    auto range = map_One_Many_.equal_range(one);
    for (auto it = range.first; it != range.second; ++it)
      result.push_back(it->second);

    if (result.empty()) {
      LOG_WARN("OneToManyUOBimap::get_vecMany: one={} not assigned to any many", one);
    }

    return result;
  }

  /// @brief Get all Many values associated with a given One key as a sorted set
  /// @param[in] one The One key to query
  /// @return Set of Many values (empty if one is not registered)
  /// @note Uses get_vecMany() internally
  std::set<Many> get_setMany(const One& one) const {
    std::vector<Many> vecMany = get_vecMany(one);
    return std::set<Many>(vecMany.begin(), vecMany.end());
  }

  /// @brief Get all registered One keys (with duplicates for each Many)
  /// @return Vector of One keys (unordered, may contain duplicates)
  /// @note For unique One keys, use get_setOne() instead
  std::vector<One> get_vecOne() const {
    std::vector<One> result;
    result.reserve(map_One_Many_.size());
    for (const auto& pair : map_One_Many_)
      result.push_back(pair.first);
    return result;
  }

  /// @brief Get all unique registered One keys as a sorted set
  /// @return Set of unique One keys (sorted)
  /// @note Uses get_vecOne() internally and converts to set
  std::set<One> get_setOne() const {
    std::vector<One> vec = get_vecOne();
    return std::set<One>(vec.begin(), vec.end());
  }

  /// @brief Get the One key associated with a given Many value
  /// @param[in] many The Many value to query
  /// @return The associated One key, or OneNotAssigned if not found
  /// @note Logs a warning if many is not registered
  One getOne(const Many& many) const {
    // Search in the many->one map
    auto it = map_Many_One_.find(many);

    // Return sentinel if not found
    if (it == map_Many_One_.end()) {
      LOG_WARN("OneToManyUOBimap::getOne: many={} not assigned to any one", many);
      return OneNotAssigned;
    }

    return it->second;
  }

  /// @brief Alias for get_vecMany()
  /// @param[in] one The One key to query
  /// @return Vector of Many values associated with the One key
  std::vector<Many> getMany(const One& one) const { return get_vecMany(one); }

  /// @brief Erase all mappings for a given One key
  /// @param[in] one The One key to erase
  /// @details Removes all Many values associated with the One key and maintains
  ///          bidirectional consistency by also removing corresponding Many->One entries.
  /// @note Logs a warning if one is not registered
  void eraseOne(const One& one) {
    // Get all Many values associated with this One
    auto range = map_One_Many_.equal_range(one);

    // If no mappings exist, log warning and return
    if (range.first == range.second) {
      LOG_WARN("OneToManyUOBimap::eraseOne: one={} is not registered", one);
      return;
    }

    // Remove corresponding entries from Many->One map
    for (auto it = range.first; it != range.second; ++it)
      map_Many_One_.erase(it->second);

    // Remove all One->Many entries for this One
    map_One_Many_.erase(range.first, range.second);
  }

  /// @brief Erase the mapping for a given Many value
  /// @param[in] many The Many value to erase
  /// @details Removes the Many->One mapping and the corresponding entry from
  ///          One->Many multimap to maintain bidirectional consistency.
  /// @note Logs a warning if many is not registered
  void eraseMany(const Many& many) {
    // Find the Many->One entry
    auto it = map_Many_One_.find(many);
    if (it == map_Many_One_.end()) {
      LOG_WARN("OneToManyUOBimap::eraseMany: many={} is not registered", many);
      return;
    }

    const One& one = it->second;

    // Remove from One->Many multimap (only the specific many value)
    std::erase_if(map_One_Many_, [&](const auto& pair) {
      return pair.first == one && EqMany{}(pair.second, many);
    });

    // Remove from Many->One map
    map_Many_One_.erase(it);
  }

  /// @brief Clear all mappings from both maps
  void clear() {
    map_One_Many_.clear();
    map_Many_One_.clear();
  }

  /// @brief Get the number of Many->One mappings (total pair count)
  /// @return Number of unique Many values registered
  size_t countMany() const { return map_Many_One_.size(); }

  /// @brief Get the number of unique One keys
  /// @return Number of unique One keys registered
  /// @note Complexity: O(n) due to set construction
  size_t countOne() const { return get_setOne().size(); }

  /// @brief Check if a One key exists in the map
  /// @param[in] one The One key to check
  /// @return true if one is registered, false otherwise
  bool hasOne(const One& one) const { return map_One_Many_.contains(one); }

  /// @brief Check if a Many value exists in the map
  /// @param[in] many The Many value to check
  /// @return true if many is registered, false otherwise
  bool hasMany(const Many& many) const { return map_Many_One_.contains(many); }

  /// @brief Reserve additional capacity for One->Many map
  /// @param[in] n_plus Additional capacity to reserve
  void rsvplusOneMany(const size_t n_plus) {
    map_One_Many_.reserve(map_One_Many_.size() + n_plus);
  }

  /// @brief Reserve additional capacity for Many->One map
  /// @param[in] n_plus Additional capacity to reserve
  void rsvplusManyOne(const size_t n_plus) {
    map_Many_One_.reserve(map_Many_One_.size() + n_plus);
  }

  /// @brief Equality comparison operator
  /// @param[in] o The other bimap to compare
  /// @return true if both maps are equal, false otherwise
  bool operator==(const OneToManyUOBimap& o) const {
    return map_One_Many_ == o.map_One_Many_
        && map_Many_One_ == o.map_Many_One_;
  }

  /// @brief Inequality comparison operator
  /// @param[in] o The other bimap to compare
  /// @return true if not equal, false otherwise
  bool operator!=(const OneToManyUOBimap& o) const {
#ifdef NODEBUG
    if (map_One_Many_ != o.map_One_Many_) return true;
    if (map_Many_One_ != o.map_Many_One_) return true;
#else
    if (map_One_Many_ != o.map_One_Many_) { LOG_WARN("OneToManyUOBimap: map_One_Many_ differs"); return true; }
    if (map_Many_One_ != o.map_Many_One_) { LOG_WARN("OneToManyUOBimap: map_Many_One_ differs"); return true; }
#endif
    return false;
  }

  /// @brief Check bidirectional consistency of the map
  /// @return true if consistent, false if internal state is corrupted
  /// @details Verifies that:
  ///          1. Every One->Many entry has a corresponding Many->One entry
  ///          2. Every Many->One entry has a corresponding One->Many entry
  /// @note Complexity: O(n) where n is the total number of mappings
  bool isConsistent() const {
    // Step 1: Verify each One->Many pair has correct Many->One reverse mapping
    for (const auto& [one, many] : map_One_Many_) {
      auto it = map_Many_One_.find(many);
      if (it == map_Many_One_.end()) return false;
      if (it->second != one) return false;
    }

    // Step 2: Verify each Many->One pair exists in One->Many multimap
    for (const auto& [many, one] : map_Many_One_) {
      auto range = map_One_Many_.equal_range(one);
      bool found = false;
      for (auto it = range.first; it != range.second; ++it) {
        if (EqMany{}(it->second, many)) {
          found = true;
          break;
        }
      }
      if (!found) return false;
    }

    return true;
  }
};
