/// @file cls_OneToManyBimap.hpp
/// @brief Bidirectional one-to-many map implementation.
///
/// @details
/// This header-only template provides an ordered bidirectional mapping between
/// unique "One" keys and multiple "Many" values.
///
/// **Invariants:**
/// - Each "Many" value maps to exactly one "One" key.
/// - Each "One" key may map to zero or more "Many" values.
/// - Internal consistency: map_One_Many_ and map_Many_One_ are always synchronized.
///
/// **Typical workflow:**
/// 1. Create an empty OneToManyBimap instance.
/// 2. Use insert() to register One->Many mappings.
/// 3. Query with getOne(), get_vecMany(), or get_setMany().
/// 4. Remove mappings with eraseOne() or eraseMany().
/// 5. Optionally verify consistency with isConsistent().
///
/// **Performance characteristics:**
/// - Most operations are O(log n) due to underlying std::map/std::multimap.
/// - Memory overhead: Two maps are maintained for bidirectional access.
/// - Iteration order: Sorted by comparison functors (CompareOne/CompareMany).
///
/// **Thread-safety:** Not thread-safe. External synchronization is required
/// for concurrent access.
///
/// @note The "One" type must be an integral or integral-like type that can be
/// cast from -1 for the sentinel value OneNotAssigned.
#pragma once

#include <map>
#include <vector>
#include <set>
#include <stdexcept>
#include <algorithm>
#include <fmt/ostream.h>
#include "ns_mylogger.hpp"

/// @brief Bidirectional one-to-many ordered map (One <-> Many).
///
/// @details
/// A template class that maintains a bidirectional mapping where:
/// - One "One" key can map to multiple "Many" values (one-to-many).
/// - Each "Many" value maps back to exactly one "One" key.
///
/// Both internal maps use std::map/std::multimap for ordered traversal.
///
/// **Usage example:**
/// @code
/// OneToManyBimap<int, std::string> bimap;
/// bimap.insert(1, "alpha");
/// bimap.insert(1, "beta");
/// bimap.insert(2, "gamma");
///
/// // Get all Many values for One=1 -> {"alpha", "beta"}
/// auto values = bimap.get_vecMany(1);
///
/// // Get the One key for "gamma" -> 2
/// int key = bimap.getOne("gamma");
/// @endcode
///
/// @tparam One         Key type for the "one" side (must be integral-like for sentinel).
/// @tparam Many        Key type for the "many" side.
/// @tparam CompareOne  Comparison functor for One (default: std::less<One>).
/// @tparam CompareMany Comparison functor for Many (default: std::less<Many>).
template<
  typename One, typename Many
, typename CompareOne  = std::less<One>
, typename CompareMany = std::less<Many>
>
class OneToManyBimap {
public:
  /// @brief Sentinel value indicating an unassigned One (-1 cast to One type).
  static constexpr One OneNotAssigned = static_cast<One>(-1);

private:
  using MapOneMany = std::multimap<One, Many, CompareOne>;
  using MapManyOne = std::map<Many, One, CompareMany>;

  MapOneMany map_One_Many_;
  MapManyOne map_Many_One_;

public:
  /// @brief Default constructor. Creates an empty bimap.
  OneToManyBimap() = default;

  /// @brief Default destructor.
  ~OneToManyBimap() = default;

  /// @brief Copy constructor (deep copy).
  OneToManyBimap(const OneToManyBimap&) = default;

  /// @brief Copy assignment operator (deep copy).
  OneToManyBimap& operator=(const OneToManyBimap&) = default;

  /// @brief Move constructor.
  OneToManyBimap(OneToManyBimap&&) noexcept = default;

  /// @brief Move assignment operator.
  OneToManyBimap& operator=(OneToManyBimap&&) noexcept = default;

  /// @brief Get the entire One->Many multimap (const reference).
  /// @return Const reference to the internal One->Many multimap.
  const MapOneMany& getMapOneMany() const { return map_One_Many_; }

  /// @brief Get the entire One->Many multimap (mutable reference).
  /// @return Mutable reference to the internal One->Many multimap.
  /// @warning Direct modification may break internal consistency.
  MapOneMany& callMapOneMany() { return map_One_Many_; }

  /// @brief Get the entire Many->One map (const reference).
  /// @return Const reference to the internal Many->One map.
  const MapManyOne& getMapManyOne() const { return map_Many_One_; }

  /// @brief Get the entire Many->One map (mutable reference).
  /// @return Mutable reference to the internal Many->One map.
  /// @warning Direct modification may break internal consistency.
  MapManyOne& callMapManyOne() { return map_Many_One_; }

  /// @brief Insert a One->Many mapping.
  /// @param[in] one  The "one" key.
  /// @param[in] many The "many" value to associate.
  /// @throws std::runtime_error if @p many is already registered to another One.
  /// @note Complexity: O(log n) for both internal maps.
  void insert(const One& one, const Many& many) {
    if (map_Many_One_.contains(many)) {
      THROW_ERROR("OneToManyBimap::insert: many={} already registered to another one", many);
    }
    map_One_Many_.emplace(one, many);
    map_Many_One_.insert(std::make_pair(many, one));
  }

  /// @brief Insert a One->Many mapping, overwriting if Many already exists.
  /// @param[in] one  The "one" key.
  /// @param[in] many The "many" value to associate.
  /// @note If @p many is already registered, its previous mapping is removed first.
  /// @note Complexity: O(log n) if many is new, O(2 log n) if many exists and is reassigned.
  void insertOrOverwrite(const One& one, const Many& many) {
    if (hasMany(many)) eraseMany(many);
    insert(one, many);
  }

  /// @brief Get all Many values associated with a given One as a vector.
  /// @param[in] one The "one" key to query.
  /// @return Vector of Many values (empty if One is not registered).
  /// @note Logs a warning if @p one has no associated Many values.
  /// @note Complexity: O(k + log n) where k is the count of Many values for this One.
  std::vector<Many> get_vecMany(const One& one) const {
    std::vector<Many> result;
    auto range = map_One_Many_.equal_range(one);
    for (auto it = range.first; it != range.second; ++it)
      result.push_back(it->second);

    if (result.empty()) {
      LOG_WARN("OneToManyBimap::get_vecMany: one={} not assigned to any many", one);
    }
    return result;
  }

  /// @brief Get all Many values associated with a given One as a set.
  /// @param[in] one The "one" key to query.
  /// @return Set of Many values (empty if One is not registered).
  /// @note Complexity: O(k log k) where k is the count of Many values for this One.
  std::set<Many> get_setMany(const One& one) const {
    auto vec = get_vecMany(one);
    return std::set<Many>(vec.begin(), vec.end());
  }

  /// @brief Get all registered One keys (with duplicates for each Many).
  /// @return Vector of One keys in map order (may contain duplicates).
  /// @note Complexity: O(n) where n is the total number of One->Many pairs.
  std::vector<One> get_vecOne() const {
    std::vector<One> result;
    result.reserve(map_One_Many_.size());
    for (const auto& p : map_One_Many_)
      result.push_back(p.first);
    return result;
  }

  /// @brief Get all unique registered One keys as a set.
  /// @return Set of unique One keys.
  /// @note Complexity: O(n log n) where n is the total number of One->Many pairs.
  std::set<One> get_setOne() const {
    auto vec = get_vecOne();
    return std::set<One>(vec.begin(), vec.end());
  }

  /// @brief Get the One key associated with a given Many.
  /// @param[in] many The "many" value to query.
  /// @return The associated One key, or OneNotAssigned if not found.
  /// @note Logs a warning if @p many is not registered.
  /// @note Complexity: O(log n).
  One getOne(const Many& many) const {
    auto it = map_Many_One_.find(many);
    if (it == map_Many_One_.end()) {
      LOG_WARN("OneToManyBimap::getOne: many={} not assigned to any one", many);
      return OneNotAssigned;
    }
    return it->second;
  }

  /// @brief Erase all mappings for a given One key.
  /// @param[in] one The "one" key whose mappings should be removed.
  /// @note Logs a warning if @p one is not registered.
  /// @note Complexity: O(k log n) where k is the count of Many values for this One.
  void eraseOne(const One& one) {
    auto range = map_One_Many_.equal_range(one);

    if (range.first == range.second) {
      LOG_WARN("OneToManyBimap::eraseOne: one={} is not registered", one);
      return;
    }

    for (auto it = range.first; it != range.second; ++it)
      map_Many_One_.erase(it->second);

    map_One_Many_.erase(range.first, range.second);
  }

  /// @brief Erase the mapping for a given Many value.
  /// @param[in] many The "many" value whose mapping should be removed.
  /// @note Logs a warning if @p many is not registered.
  /// @note Complexity: O(k + log n) where k is the count of Many values for the associated One.
  void eraseMany(const Many& many) {
    auto it = map_Many_One_.find(many);
    if (it == map_Many_One_.end()) {
      LOG_WARN("OneToManyBimap::eraseMany: many={} is not registered", many);
      return;
    }

    const One& one = it->second;

    // Remove from One->Many multimap: find the specific (one, many) pair
    auto range = map_One_Many_.equal_range(one);
    for (auto itM = range.first; itM != range.second; ++itM) {
      if (itM->second == many) {
        map_One_Many_.erase(itM);
        break;  // Only one such entry exists
      }
    }

    map_Many_One_.erase(it);
  }

  /// @brief Clear all mappings.
  void clear() {
    map_One_Many_.clear();
    map_Many_One_.clear();
  }

  /// @brief Get the number of Many->One mappings.
  /// @return Total count of Many entries.
  size_t countMany() const { return map_Many_One_.size(); }

  /// @brief Get the number of unique One keys.
  /// @return Count of unique One entries.
  /// @note Complexity: O(u log n) where u is the count of unique One keys.
  size_t countOne() const {
    size_t count = 0;
    for (auto it = map_One_Many_.begin(); it != map_One_Many_.end();
         it = map_One_Many_.upper_bound(it->first))
      ++count;
    return count;
  }

  /// @brief Check if a One key exists in the map.
  /// @param[in] one The "one" key to check.
  /// @return True if @p one has at least one associated Many.
  /// @note Complexity: O(log n).
  bool hasOne(const One& one) const { return map_One_Many_.contains(one); }

  /// @brief Check if a Many value exists in the map.
  /// @param[in] many The "many" value to check.
  /// @return True if @p many is registered.
  /// @note Complexity: O(log n).
  bool hasMany(const Many& many) const { return map_Many_One_.contains(many); }

  /// @brief Equality comparison.
  /// @param[in] o The other bimap to compare.
  /// @return True if both internal maps are equal.
  bool operator==(const OneToManyBimap& o) const {
    return map_One_Many_ == o.map_One_Many_ &&
           map_Many_One_ == o.map_Many_One_;
  }

  /// @brief Inequality comparison.
  /// @param[in] o The other bimap to compare.
  /// @return True if any internal map differs.
  bool operator!=(const OneToManyBimap& o) const {
#ifdef NODEBUG
    if (map_One_Many_ != o.map_One_Many_) return true;
    if (map_Many_One_ != o.map_Many_One_) return true;
#else
    if (map_One_Many_ != o.map_One_Many_) { LOG_WARN("OneToManyBimap: map_One_Many_ differs"); return true; }
    if (map_Many_One_ != o.map_Many_One_) { LOG_WARN("OneToManyBimap: map_Many_One_ differs"); return true; }
#endif
    return false;
  }

  /// @brief Copy all mappings from another bimap (replaces current content).
  /// @param[in] other The source bimap to copy from.
  void copy_from(const OneToManyBimap& other) {
    map_One_Many_ = other.map_One_Many_;
    map_Many_One_ = other.map_Many_One_;
  }

  /// @brief Check internal consistency of the bidirectional mapping.
  /// @return True if both maps are mutually consistent, false otherwise.
  /// @note Use this for debugging or after direct map modification.
  /// @note Complexity: O(n * log n) where n is the total number of entries.
  bool isConsistent() const {
    // Verify every One->Many entry has a matching Many->One entry
    for (const auto& [one, many] : map_One_Many_) {
      auto it = map_Many_One_.find(many);
      if (it == map_Many_One_.end() || it->second != one)
        return false;
    }
    // Verify every Many->One entry has a matching One->Many entry
    for (const auto& [many, one] : map_Many_One_) {
      auto range = map_One_Many_.equal_range(one);
      bool found = false;
      for (auto it = range.first; it != range.second; ++it) {
        if (it->second == many) {
          found = true;
          break;
        }
      }
      if (!found)
        return false;
    }
    return true;
  }
};
