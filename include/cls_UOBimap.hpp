/// @file cls_UOBimap.hpp
/// @brief Unordered bidirectional map implementation
/// @details
/// This header-only template provides hash-based bidirectional mapping between two key types.
///
/// **Main Features:**
/// - O(1) average-case lookup in both directions (A→B and B→A)
/// - Automatic consistency maintenance between forward and reverse mappings
/// - Forbidden value support for "not found" sentinel values
/// - Conditional insertion with predicate-based exceptions
/// - Thread-safety: No. User must provide external synchronization for concurrent access.
/// - Memory layout: Two independent unordered_map instances maintained in sync
///
/// **Typical Workflow:**
/// 1. Construct UOBimap with optional forbidden values and expected size
/// 2. Insert bidirectional pairs using insert() or specialized variants
/// 3. Query in either direction using getAB()/getBA() or safe variants
/// 4. Maintain consistency using isConsistent() in debug builds
///
/// **Performance Characteristics:**
/// - Insertion: O(1) average, O(n) worst case
/// - Lookup: O(1) average, O(n) worst case
/// - Deletion: O(1) average, O(n) worst case
/// - Memory overhead: 2× standard unordered_map (one for each direction)
///
/// @note This class is header-only. No separate compilation required.
#pragma once

#include <vector>
#include <set>
#include <unordered_map>
#include <stdexcept>
#include <utility>
#include <string>
#include <fmt/ostream.h>
#include <concepts> // c++20, for std::equality_comparable
#include "ns_mylogger.hpp"

/// @class UOBimap
/// @brief Header-only template implementation of unordered bidirectional map
///
/// @details
/// UOBimap (Unordered Bidirectional Map) maintains two synchronized hash maps
/// for constant-time lookups in both A→B and B→A directions.
///
/// **Invariants:**
/// - map_AB_.size() == map_BA_.size() at all times
/// - For every (a, b) in map_AB_, (b, a) exists in map_BA_
/// - No key equals its corresponding forbidden value (notFoundBA_ or notFoundAB_)
/// - All operations maintain bidirectional consistency
///
/// **Primary Use Cases:**
/// - ID↔Name mappings requiring fast reverse lookup
/// - Entity↔Handle bidirectional associations
/// - Symbol table with bidirectional resolution
///
/// **Example Usage:**
/// @code
/// UOBimap<int, std::string> idNameMap(-1, "");
/// idNameMap.insert(100, "Alice");
/// idNameMap.insert(200, "Bob");
///
/// std::string name = idNameMap.getAB(100);        // "Alice"
/// int id = idNameMap.getBA("Bob");                // 200
/// std::string safe = idNameMap.getABorDefault(999); // "" (not found)
/// @endcode
///
/// @tparam A First key type (must satisfy std::equality_comparable)
/// @tparam B Second key type (must satisfy std::equality_comparable)
/// @tparam HashA Hash function for type A (default: std::hash<A>)
/// @tparam EqA Equality comparison function for type A (default: std::equal_to<A>)
/// @tparam HashB Hash function for type B (default: std::hash<B>)
/// @tparam EqB Equality comparison function for type B (default: std::equal_to<B>)
///
/// @note Thread-safety: No. This class is not thread-safe. External synchronization required for concurrent access.
/// @note Complexity: All basic operations (insert/erase/lookup) are O(1) average case, O(n) worst case.
template<
    std::equality_comparable A
  , std::equality_comparable B
  , typename HashA = std::hash<A>
  , typename EqA   = std::equal_to<A>
  , typename HashB = std::hash<B>
  , typename EqB   = std::equal_to<B>
>
class UOBimap {
private:
  /// @brief Type alias for A->B direction map
  using MapAB = std::unordered_map<A, B, HashA, EqA>;

  /// @brief Type alias for B->A direction map
  using MapBA = std::unordered_map<B, A, HashB, EqB>;

  /// @brief Default B value returned when A key is not found
  B notFoundAB_;

  /// @brief Default A value returned when B key is not found
  A notFoundBA_;

  /// @brief A-side map (A → B)
  MapAB map_AB_;

  /// @brief B-side map (B → A)
  MapBA map_BA_;

  /// @brief Reserve capacity for A-side map
  /// @param n Expected number of elements (bucket count)
  void reserveA(const size_t n) { map_AB_.reserve(n); };

  /// @brief Reserve capacity for B-side map
  /// @param n Expected number of elements (bucket count)
  void reserveB(const size_t n) { map_BA_.reserve(n); };

  /// @brief Check for forbidden values (must be called before insertion)
  void checkForbiddenValue(const A& a, const B& b) const {
    if (a == notFoundBA_) {
      std::string msg = fmt::format("UOBimap::checkForbiddenValue: A-side value '{}' is forbidden (notFoundBA_)", a);
      LOG_ERROR(msg);
      throw std::invalid_argument("UOBimap::checkForbiddenValue: Cannot register forbidden value on A-side");
    }
    if (b == notFoundAB_) {
      std::string msg = fmt::format("UOBimap::checkForbiddenValue: B-side value '{}' is forbidden (notFoundAB_)", b);
      LOG_ERROR(msg);
      throw std::invalid_argument("UOBimap::checkForbiddenValue: Cannot register forbidden value on B-side");
    }
  };

public:
  /// @brief Constructor
  /// @param notFoundBA Forbidden value for A-side (default is A's default value)
  /// @param notFoundAB Forbidden value for B-side (default is B's default value)
  /// @param expectedSize Expected number of elements to register (0 means no reservation)
  explicit UOBimap(
    const A notFoundBA = A{}, const B notFoundAB = B{}, const size_t expectedSize = 0)
      : notFoundAB_(notFoundAB), notFoundBA_(notFoundBA) {
        if (expectedSize > 0) reserve(expectedSize);
      };

  /// @brief Default destructor
  ~UOBimap() = default;

  /// @brief Copy constructor
  UOBimap(const UOBimap&) = default;

  /// @brief Copy assignment operator
  UOBimap& operator=(const UOBimap&) = default;

  /// @brief Move constructor
  UOBimap(UOBimap&&) = default;

  /// @brief Move assignment operator
  UOBimap& operator=(UOBimap&&) = default;

  /// @brief Get const reference to A-side map
  /// @return Const reference to A-side map
  const MapAB& getMapAB() const { return map_AB_; };

  /// @brief Get const reference to B-side map
  /// @return Const reference to B-side map
  const MapBA& getMapBA() const { return map_BA_; };

  /// @brief Get mutable reference to A-side map
  /// @return Mutable reference to A-side map
  MapAB& callMapAB() { return map_AB_; };

  /// @brief Get mutable reference to B-side map
  /// @return Mutable reference to B-side map
  MapBA& callMapBA() { return map_BA_; };

  /// @brief Register a bidirectional pair
  /// @param a A-side key
  /// @param b B-side key
  /// @throws std::invalid_argument If A-side key or B-side key already exists
  /// @note Both A-side and B-side keys must be unique.
  void insert(const A& a, const B& b) {
    checkForbiddenValue(a, b);

    if (map_AB_.contains(a)) {
      std::string msg = fmt::format("UOBimap::insert: A-side key a = {} already exists", a);
      LOG_ERROR(msg);
      throw std::invalid_argument("UOBimap::insert: A-side key already exists");
    }
    if (map_BA_.contains(b)) {
      std::string msg = fmt::format("UOBimap::insert: B-side key b = {} already exists", b);
      LOG_ERROR(msg);
      throw std::invalid_argument("UOBimap::insert: B-side key already exists");
    }

    auto tmp_AB = std::make_pair(a, b);
    auto tmp_BA = std::make_pair(b, a);

    map_AB_.insert(std::move(tmp_AB));
    map_BA_.insert(std::move(tmp_BA));
  };

  /// @brief Insert allowing duplicate if A-side key matches exception value
  /// @param a A-side key
  /// @param b B-side key
  /// @param a_allow_value Exception value for A-side key duplication
  void insertWithAException(const A& a, const B& b, const A& a_allow_value) {
    checkForbiddenValue(a, b);

    const bool is_exception = (a == a_allow_value);

    if (map_AB_.contains(a) && !is_exception) {
      LOG_ERROR("UOBimap::insertWithAException: A-side key a = {} already exists", a);
      throw std::invalid_argument("UOBimap::insertWithAException: A-side key already exists");
    }
    if (map_BA_.contains(b)) {
      LOG_ERROR("UOBimap::insertWithAException: B-side key b = {} already exists", b);
      throw std::invalid_argument("UOBimap::insertWithAException: B-side key already exists");
    }

    auto tmp_AB = std::make_pair(a, b);
    auto tmp_BA = std::make_pair(b, a);

    map_AB_.insert(std::move(tmp_AB));
    map_BA_.insert(std::move(tmp_BA));
  };

  /// @brief Insert allowing duplicate if B-side key matches exception value
  /// @param a A-side key
  /// @param b B-side key
  /// @param b_allow_value Exception value for B-side key duplication
  void insertWithBException(const A& a, const B& b, const B& b_allow_value) {
    checkForbiddenValue(a, b);

    const bool is_exception = (b == b_allow_value);

    if (map_AB_.contains(a)) {
      LOG_ERROR("UOBimap::insertWithBException: A-side key a = {} already exists", a);
      throw std::invalid_argument("UOBimap::insertWithBException: A-side key already exists");
    }
    if (map_BA_.contains(b) && !is_exception) {
      LOG_ERROR("UOBimap::insertWithBException: B-side key b = {} already exists", b);
      throw std::invalid_argument("UOBimap::insertWithBException: B-side key already exists");
    }

    auto tmp_AB = std::make_pair(a, b);
    auto tmp_BA = std::make_pair(b, a);

    map_AB_.insert(std::move(tmp_AB));
    map_BA_.insert(std::move(tmp_BA));
  };


  /// @brief Insert allowing duplicate if A-side key satisfies predicate condition
  /// @tparam Pred Predicate function type for A-side key condition
  /// @param a A-side key
  /// @param b B-side key
  /// @param allowA Predicate lambda to judge if A-side key duplication is allowed
  /// @throws std::invalid_argument If uniqueness constraint is violated
  /// @note This method may not be necessary. Consider removal if unused in codebase.
  template<typename Pred>
  void insertAllowingACondition(const A& a, const B& b, Pred allowA) {
    checkForbiddenValue(a, b);

    const bool is_exception = allowA(a);

    if (map_AB_.contains(a) && !is_exception) {
      LOG_ERROR("UOBimap::insertAllowingACondition: A-side key a = {} already exists", a);
      throw std::invalid_argument("UOBimap::insertAllowingACondition: A-side key already exists");
    }
    if (map_BA_.contains(b)) {
      LOG_ERROR("UOBimap::insertAllowingACondition: B-side key b = {} already exists", b);
      throw std::invalid_argument("UOBimap::insertAllowingACondition: B-side key already exists");
    }

    auto tmp_AB = std::make_pair(a, b);
    auto tmp_BA = std::make_pair(b, a);

    map_AB_.insert(std::move(tmp_AB));
    map_BA_.insert(std::move(tmp_BA));
  };

  /// @brief Insert allowing duplicate if B-side key satisfies predicate condition
  /// @tparam Pred Predicate function type for B-side key condition
  /// @param a A-side key
  /// @param b B-side key
  /// @param allowB Predicate lambda to judge if B-side key duplication is allowed
  /// @throws std::invalid_argument If uniqueness constraint is violated
  /// @note This method may not be necessary. Consider removal if unused in codebase.
  template<typename Pred>
  void insertAllowingBCondition(const A& a, const B& b, Pred allowB) {
    checkForbiddenValue(a, b);

    const bool is_exception = allowB(b);

    if (map_AB_.contains(a)) {
      LOG_ERROR("UOBimap::insertAllowingBCondition: A-side key a = {} already exists", a);
      throw std::invalid_argument("UOBimap::insertAllowingBCondition: A-side key already exists");
    }
    if (map_BA_.contains(b) && !is_exception) {
      LOG_ERROR("UOBimap::insertAllowingBCondition: B-side key b = {} already exists", b);
      throw std::invalid_argument("UOBimap::insertAllowingBCondition: B-side key already exists");
    }

    auto tmp_AB = std::make_pair(a, b);
    auto tmp_BA = std::make_pair(b, a);

    map_AB_.insert(std::move(tmp_AB));
    map_BA_.insert(std::move(tmp_BA));
  };

  /// @brief Register A→B pair (overwrite if existing pair exists)
  /// @param a A-side key
  /// @param b B-side key
  void insertOrOverwrite(const A& a, const B& b) {
    checkForbiddenValue(a, b);

    if (hasA(a)) eraseA(a);
    if (hasB(b)) eraseB(b);

    insert(a, b);
  };

  /// @brief Get A→B value (const reference)
  /// @param a A-side key
  /// @return Const reference to corresponding B value
  /// @throws std::out_of_range If key does not exist
  const B& getAB(const A& a) const {
    auto it = map_AB_.find(a);
    if (it == map_AB_.end()) {
      throw std::out_of_range("UOBimap::getAB: key not found");
    }
    return it->second;
  }

  /// @brief Get A→B value (returns notFoundAB_ if not registered)
  B getABorDefault(const A& a) const noexcept {
    auto it = map_AB_.find(a);
    return it != map_AB_.end() ? it->second : notFoundAB_;
  };

  /// @brief Get A→B value (mutable reference)
  /// @param a A-side key
  /// @return Mutable reference to corresponding B value
  /// @throws std::out_of_range If key does not exist
  B& callAB(const A& a) {
    auto it = map_AB_.find(a);
    if (it == map_AB_.end()) {
      throw std::out_of_range("UOBimap::callAB: key not found");
    }
    return it->second;
  }

  /// @brief Get B→A value (const reference)
  /// @param b B-side key
  /// @return Const reference to corresponding A value
  /// @throws std::out_of_range If key does not exist
  const A& getBA(const B& b) const {
    auto it = map_BA_.find(b);
    if (it == map_BA_.end()) {
      throw std::out_of_range("UOBimap::getBA: key not found");
    }
    return it->second;
  }

  /// @brief Get B→A value (returns notFoundBA_ if not registered)
  A getBAorDefault(const B& b) const noexcept {
    auto it = map_BA_.find(b);
    return it != map_BA_.end() ? it->second : notFoundBA_;
  };

  /// @brief Get B→A value (mutable reference)
  /// @param b B-side key
  /// @return Mutable reference to corresponding A value
  /// @throws std::out_of_range If key does not exist
  A& callBA(const B& b) {
    auto it = map_BA_.find(b);
    if (it == map_BA_.end()) {
      throw std::out_of_range("UOBimap::callBA: key not found");
    }
    return it->second;
  }
  /// @brief Erase A→B pair (does nothing if not exists)
  /// @param a A-side key
  void eraseA(const A& a) {
      auto it = map_AB_.find(a);
      if (it == map_AB_.end()){
        LOG_WARN("UOBimap::eraseA: A-side key a = {} does not exist", a);
        return;
      }
      map_BA_.erase(it->second);
      map_AB_.erase(it);
  };

  /// @brief Erase B→A pair (does nothing if not exists)
  /// @param b B-side key
  void eraseB(const B& b) {
    auto it = map_BA_.find(b);
    if (it == map_BA_.end()){
      LOG_WARN("UOBimap::eraseB: B-side key b = {} does not exist", b);
      return;
    }
    map_AB_.erase(it->second);
    map_BA_.erase(it);
  };

  /// @brief Clear all mappings
  void clear() { map_AB_.clear(); map_BA_.clear(); };

  /// @brief Get number of registered elements (number of pairs)
  /// @return Number of registered pairs
  size_t size() const { return map_AB_.size(); };

  /// @brief Check if A-side key exists
  /// @param a A-side key
  /// @return true if registered, false if not registered
  bool hasA(const A& a) const { return map_AB_.contains(a); };

  /// @brief Check if B-side key exists
  /// @param b B-side key
  /// @return true if registered, false if not registered
  bool hasB(const B& b) const { return map_BA_.contains(b); };

  /// @brief Reserve capacity for both A/B sides simultaneously
  void reserve(const size_t n){ reserveA(n); reserveB(n); };

  /// @brief Reserve additional capacity on top of existing map size
  /// @param n_plus Number of additional elements to reserve
  void reserveAdditional(const size_t n_plus){
    reserveA(map_AB_.size() + n_plus);
    reserveB(map_BA_.size() + n_plus);
  };

  /// @brief Equality comparison between bidirectional maps
  /// @param other UOBimap to compare against
  /// @return true if both directions of maps are identical
  bool operator==(const UOBimap& other) const {
    return map_AB_ == other.map_AB_
        && map_BA_ == other.map_BA_;
  };

  /// @brief Inequality comparison between bidirectional maps
  /// @param other UOBimap to compare against
  /// @return true if any mismatch exists
  bool operator!=(const UOBimap& other) const {
    return !(*this == other);
  };

  /// @brief Check internal consistency of bidirectional map
  /// @return true if consistency is maintained, false if corrupted
  bool isConsistent() const {
    if (map_AB_.size() != map_BA_.size()) return false;
    for (const auto& [a, b] : map_AB_) {
      auto it = map_BA_.find(b);
      if (it == map_BA_.end() || it->second != a) return false;
    }
    return true;
  };

  /// @brief Get vector of A-side keys
  /// @return Vector of A-side keys
  std::vector<A> get_vec_A() const {
    std::vector<A> vec_a;
    vec_a.reserve(map_AB_.size());
    for (const auto& kv : map_AB_) {
      vec_a.push_back(kv.first);
    }
    return vec_a;
  };

  /// @brief Get set of A-side keys
  /// @return Set of A-side keys
  std::set<A> get_set_A() const {
    std::vector<A> vec_a = get_vec_A();
    return std::set<A>(vec_a.begin(), vec_a.end());
  };

  /// @brief Get vector of B-side keys
  /// @return Vector of B-side keys
  std::vector<B> get_vec_B() const {
    std::vector<B> vec_b;
    vec_b.reserve(map_BA_.size());
    for (const auto& kv : map_BA_) {
      vec_b.push_back(kv.first);
    }
    return vec_b;
  };

  /// @brief Get set of B-side keys
  /// @return Set of B-side keys
  std::set<B> get_set_B() const {
    std::vector<B> vec_b = get_vec_B();
    return std::set<B>(vec_b.begin(), vec_b.end());
  };

};
