/// @file cls_Bimap.hpp
/// @brief Bidirectional map (Bimap) template class
/// @details
/// Header-only implementation of a bidirectional map for efficient two-way lookups.
/// Maintains two std::map instances internally to provide O(log n) lookup in both directions.
///
/// **Main workflow:**
/// 1. Construct Bimap with optional forbidden values (defaults for not-found returns)
/// 2. Insert key-value pairs using insert() or insertOrOverwrite()
/// 3. Perform bidirectional lookups using getAB()/getBA() or their variants
/// 4. Erase pairs using eraseA()/eraseB() or clear all with clear()
///
/// **Thread-safety:** No. External synchronization required for concurrent access.
///
/// **Invariants:** Both internal maps must remain synchronized. Use isConsistent() to verify.
#pragma once

#include <map>
#include <string>
#include <stdexcept>
#include <utility>
#include <fmt/ostream.h>
#include "ns_mylogger.hpp"

/// @class Bimap
/// @brief Header-only bidirectional map (Bimap) template implementation
/// @details
/// Provides efficient two-way mapping between keys of type A and B.
/// Both A and B must support equality comparison (operator==) for forbidden value checking.
/// Both A and B must be copyable and default-constructible.
///
/// **Responsibilities:**
/// - Maintain bidirectional mapping between A and B keys
/// - Ensure internal consistency between forward (A→B) and reverse (B→A) maps
/// - Prevent insertion of forbidden values (notFoundBA_, notFoundAB_)
/// - Provide exception-safe operations for insertion and deletion
///
/// **Invariants:**
/// - For every entry (a, b) in map_AB_, there exists entry (b, a) in map_BA_
/// - No key can equal its corresponding forbidden value
/// - map_AB_.size() == map_BA_.size() always
///
/// **Usage example:**
/// @code
/// Bimap<int, std::string> bimap(-1, "");
/// bimap.insert(1, "one");
/// bimap.insert(2, "two");
/// std::string s = bimap.getAB(1);        // returns "one"
/// int i = bimap.getBA("two");            // returns 2
/// int missing = bimap.getABorDefault(99); // returns -1 (notFoundBA_)
/// @endcode
///
/// @tparam A First key type (must support operator== and be copyable)
/// @tparam B Second key type (must support operator== and be copyable)
/// @tparam CompareA Comparison function for type A (default: std::less<A>)
/// @tparam CompareB Comparison function for type B (default: std::less<B>)
///
/// @note Thread-safety: No. This class is not thread-safe. External synchronization required.
/// @note Complexity: All lookup, insertion, and deletion operations are O(log n).
template<
  typename A, typename B
  , typename CompareA = std::less<A>
  , typename CompareB = std::less<B>
>
class Bimap {
  /// @brief Type definition for A->B direction map
  using MapAB = std::map<A, B, CompareA>;
  /// @brief Type definition for B->A direction map
  using MapBA = std::map<B, A, CompareB>;

  /// @brief Map for A side (A → B)
  MapAB map_AB_;
  /// @brief Map for B side (B → A)
  MapBA map_BA_;

  /// @brief Default value returned when B key is not found (used in getBAorDefault)
  A notFoundBA_;

  /// @brief Default value returned when A key is not found (used in getABorDefault)
  B notFoundAB_;

public:
  /// @brief Constructor with forbidden values
  /// @param[in] notFoundBA Default value returned by getBAorDefault when B key not found
  /// @param[in] notFoundAB Default value returned by getABorDefault when A key not found
  /// @note These forbidden values cannot be inserted as actual keys via insert()
  Bimap( const A& notFoundBA = A{}, const B& notFoundAB = B{} )
    : notFoundBA_(notFoundBA), notFoundAB_(notFoundAB) {};

  /// @brief Default destructor
  ~Bimap() = default;

  /// @brief Copy constructor
  Bimap(const Bimap&) = default;

  /// @brief Copy assignment operator
  Bimap& operator=(const Bimap&) = default;

  /// @brief Move constructor
  Bimap(Bimap&&) = default;

  /// @brief Move assignment operator
  Bimap& operator=(Bimap&&) = default;

  /// @brief Check if values match forbidden values
  /// @param[in] a Input value of type A to check
  /// @param[in] b Input value of type B to check
  /// @throws std::invalid_argument If a equals notFoundBA_ or b equals notFoundAB_
  void checkForbiddenValue( const A& a, const B& b ) const {
    if (a == notFoundBA_) {
      THROW_ERROR("Bimap::checkForbiddenValue: A value '{}' matches forbidden value notFoundBA_", a);
    }
    if (b == notFoundAB_) {
      THROW_ERROR("Bimap::checkForbiddenValue: B value '{}' matches forbidden value notFoundAB_", b);
    }
  };

  /// @brief Get A->B map (read-only)
  /// @return Const reference to internal A->B map
  /// @note Complexity: O(1)
  const MapAB& getMapAB() const { return map_AB_; };

  /// @brief Get B->A map (read-only)
  /// @return Const reference to internal B->A map
  /// @note Complexity: O(1)
  const MapBA& getMapBA() const { return map_BA_; };

  /// @brief Get A->B map (writable)
  /// @return Mutable reference to internal A->B map
  /// @warning Direct modification may break bidirectional consistency. Use with caution.
  /// @note Complexity: O(1)
  MapAB& callMapAB() { return map_AB_; };

  /// @brief Get B->A map (writable)
  /// @return Mutable reference to internal B->A map
  /// @warning Direct modification may break bidirectional consistency. Use with caution.
  /// @note Complexity: O(1)
  MapBA& callMapBA() { return map_BA_; };

  /// @brief Register pair of A and B
  /// @param[in] a Key of type A
  /// @param[in] b Key of type B
  /// @throws std::invalid_argument If a or b equals forbidden value
  /// @throws std::runtime_error If a or b is already registered
  /// @note Complexity: O(log n)
  void insert(const A& a, const B& b) {
    // check forbidden value
    checkForbiddenValue(a, b);

    // check if a or b is already registered
    if (map_AB_.contains(a)) {
      THROW_ERROR("Bimap::insert: A key a={} already exists", a);
    }
    if (map_BA_.contains(b)) {
      THROW_ERROR("Bimap::insert: B key b={} already exists", b);
    }

    // prepare two inserts transactionally, then commit them together
    auto tmp_AB = std::make_pair(a, b);
    auto tmp_BA = std::make_pair(b, a);

    // commit to the main maps if both inserts succeed
    map_AB_.insert(std::move(tmp_AB));
    map_BA_.insert(std::move(tmp_BA));
  };

  /// @brief Register or overwrite pair of A and B
  /// @param[in] a Key of type A
  /// @param[in] b Key of type B
  /// @throws std::invalid_argument If a or b equals forbidden value
  /// @note If either key exists, the old pair is erased before inserting new pair
  /// @note Complexity: O(log n)
  void insertOrOverwrite(const A& a, const B& b) {
    // check forbidden value
    checkForbiddenValue(a, b);

    // if A key already exists, erase it
    if (hasA(a)) eraseA(a);

    // if B key already exists, erase it
    if (hasB(b)) eraseB(b);

    // insert new pair
    insert(a, b);
  };

  /// @brief Get A->B value (read-only)
  /// @param[in] a Key of type A
  /// @return Const reference to corresponding B value
  /// @throws std::out_of_range If key not found
  /// @note Complexity: O(log n)
  const B& getAB(const A& a) const {
    auto it = map_AB_.find(a);
    if (it == map_AB_.end()) {
      THROW_ERROR("Bimap::getAB: key not found (a={})", a);
    }
    return it->second;
  };

  /// @brief Get A->B value (return notFoundAB_ if key not found)
  /// @param[in] a Key of type A
  /// @return Corresponding B value, or notFoundAB_ if not found
  /// @note Complexity: O(log n)
  B getABorDefault(const A& a) const noexcept {
    auto it = map_AB_.find(a);
    return it != map_AB_.end() ? it->second : notFoundAB_;
  };

  /// @brief Get A->B value (writable)
  /// @param[in] a Key of type A
  /// @return Mutable reference to corresponding B value
  /// @throws std::out_of_range If key not found
  /// @warning Modifying the returned reference does not update the reverse map
  /// @note Complexity: O(log n)
  B& callAB(const A& a) {
    auto it = map_AB_.find(a);
    if (it == map_AB_.end()) {
      THROW_ERROR("Bimap::callAB: key not found (a={})", a);
    }
    return it->second;
  };

  /// @brief Get B->A value (read-only)
  /// @param[in] b Key of type B
  /// @return Const reference to corresponding A value
  /// @throws std::out_of_range If key not found
  /// @note Complexity: O(log n)
  const A& getBA(const B& b) const {
    auto it = map_BA_.find(b);
    if (it == map_BA_.end()) {
      THROW_ERROR("Bimap::getBA: key not found (b={})", b);
    }
    return it->second;
  };

  /// @brief Get B->A value (return notFoundBA_ if key not found)
  /// @param[in] b Key of type B
  /// @return Corresponding A value, or notFoundBA_ if not found
  /// @note Complexity: O(log n)
  A getBAorDefault(const B& b) const noexcept {
    auto it = map_BA_.find(b);
    return it != map_BA_.end() ? it->second : notFoundBA_;
  };

  /// @brief Get B->A value (writable)
  /// @param[in] b Key of type B
  /// @return Mutable reference to corresponding A value
  /// @throws std::out_of_range If key not found
  /// @warning Modifying the returned reference does not update the forward map
  /// @note Complexity: O(log n)
  A& callBA(const B& b) {
    auto it = map_BA_.find(b);
    if (it == map_BA_.end()) {
      THROW_ERROR("Bimap::callBA: key not found (b={})", b);
    }
    return it->second;
  };

  /// @brief Erase A->B pair by A key
  /// @param[in] a Key of type A to erase
  /// @note If key does not exist, logs warning and returns without error
  /// @note Complexity: O(log n)
  void eraseA(const A& a) {
    auto it = map_AB_.find(a);
    if (it == map_AB_.end()) {
      LOG_WARN("Bimap::eraseA: A key a={} does not exist", a);
      return;
    }
    map_BA_.erase(it->second);
    map_AB_.erase(it);
  };

  /// @brief Erase B->A pair by B key
  /// @param[in] b Key of type B to erase
  /// @note If key does not exist, logs warning and returns without error
  /// @note Complexity: O(log n)
  void eraseB(const B& b) {
    auto it = map_BA_.find(b);
    if (it == map_BA_.end()) {
      LOG_WARN("Bimap::eraseB: B key b={} does not exist", b);
      return;
    }
    map_AB_.erase(it->second);
    map_BA_.erase(it);
  };

  /// @brief Erase all mappings
  /// @note Complexity: O(n)
  void clear() {
    map_AB_.clear();
    map_BA_.clear();
  };

  /// @brief Get number of registered pairs
  /// @return Number of key pairs in the bimap
  /// @note Complexity: O(1)
  size_t size() const { return map_AB_.size(); };

  /// @brief Check existence of A key
  /// @param[in] a Input value of type A
  /// @return true if A key exists, false otherwise
  /// @note Complexity: O(log n)
  bool hasA(const A& a) const { return map_AB_.contains(a); };

  /// @brief Check existence of B key
  /// @param[in] b Input value of type B
  /// @return true if B key exists, false otherwise
  /// @note Complexity: O(log n)
  bool hasB(const B& b) const { return map_BA_.contains(b); };

  /// @brief Equality comparison operator
  /// @param[in] other Another Bimap to compare with
  /// @return true if both bimaps contain identical mappings, false otherwise
  /// @note Complexity: O(n)
  bool operator==(const Bimap& other) const {
    return map_AB_ == other.map_AB_ && map_BA_ == other.map_BA_;
  };

  /// @brief Inequality comparison operator
  /// @param[in] other Another Bimap to compare with
  /// @return true if bimaps differ, false otherwise
  /// @note Complexity: O(n)
  bool operator!=(const Bimap& other) const {
    return !(*this == other);
  };

  /// @brief Check consistency of the bidirectional map
  /// @return true if consistent, false if internal maps are desynchronized
  /// @note Complexity: O(n)
  /// @note Useful for debugging after direct map modifications via callMapAB()/callMapBA()
  bool isConsistent() const {
    for (const auto& [a, b] : map_AB_) {
      auto it = map_BA_.find(b);
      if (it == map_BA_.end() || it->second != a) return false;
    }
    for (const auto& [b, a] : map_BA_) {
      auto it = map_AB_.find(a);
      if (it == map_AB_.end() || it->second != b) return false;
    }
    return true;
  };
};
