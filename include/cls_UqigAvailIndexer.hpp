/// @file cls_UqigAvailIndexer.hpp
/// @brief Bidirectional mapping between UqigAvail identifiers and dense matrix indices.
/// @details Provides efficient bidirectional conversion between sparse UqigAvail values
///          and 0-based consecutive indices suitable for Eigen matrix/vector operations.
///          This indexer is primarily used to map sparse identifiers to dense matrix row/column
///          indices, enabling efficient matrix-based computations while maintaining meaningful
///          identifier-based access patterns.
///
/// @par Typical Workflow:
/// 1. Construct an empty UqigAvailIndexer object.
/// 2. Call build() with a vector of UqigAvail values to establish the mapping.
/// 3. Use getIndex() to convert UqigAvail to dense index for matrix operations.
/// 4. Use getUqigAvail() to convert dense index back to UqigAvail.
/// 5. Use save()/load() for binary serialization.
///
/// @par Thread Safety:
/// - Not thread-safe. External synchronization required for concurrent access.
/// - Concurrent reads after build() are safe only if no modifications occur.
///
/// @par Complexity:
/// - build(): O(n log n) where n = number of unique UqigAvail values.
/// - getIndex(), getUqigAvail(): O(1) average (hash table lookup).
///
/// @par Memory Layout:
/// - Internal storage uses UOBimap with hash-based bidirectional lookup.
/// - UqigAvail values are stored in sorted ascending order internally.
///
/// @par Sentinel Values:
/// - UqigAvailNotFound: Returned when a dense index has no corresponding UqigAvail.
/// - IndexNotFound: Returned when a UqigAvail has no corresponding dense index.
#pragma once
#include <vector>
#include <unordered_map>
#include "cls_UOBimap.hpp"
#include "ns_type_definitions.hpp"

using namespace index_type_definitions;

/// @class UqigAvailIndexer
/// @brief Maps sparse UqigAvail identifiers to 0-based consecutive indices.
/// @details Enables Eigen Vector/Matrix rows to be indexed by UqigAvail values.
///          Even when UqigAvail values are sparse integers, this class provides
///          a compact 0-based index mapping suitable for dense matrix storage.
///
/// @par Usage Example:
/// @code
/// UqigAvailIndexer indexer;
/// std::vector<UqigAvail> avails = {100, 200, 150};
/// indexer.build(avails);
/// // After build: 100->0, 150->1, 200->2 (sorted order)
/// Index idx = indexer.getIndex(150);  // returns 1
/// UqigAvail ua = indexer.getUqigAvail(2);  // returns 200
/// @endcode
///
/// @invariant After build(), the mapping is immutable until the next build() call.
/// @invariant Dense indices are always consecutive starting from 0.
class UqigAvailIndexer {
private:
  /// @brief Bidirectional map: dense index <-> original UqigAvail.
  /// @note UqigAvailNotFound and IndexNotFound are sentinel/forbidden values.
  ///       - If key=UqigAvail is not found, returns IndexNotFound.
  ///       - If key=Index is not found, returns UqigAvailNotFound.
  UOBimap<UqigAvail, Index> bimap_UqigAvail_Index{ UqigAvailNotFound, IndexNotFound };

public:
  /// @brief Default constructor. Creates an empty indexer.
  UqigAvailIndexer() = default;

  /// @brief Copy constructor.
  UqigAvailIndexer(const UqigAvailIndexer& other) = default;

  /// @brief Move constructor.
  UqigAvailIndexer(UqigAvailIndexer&& other) noexcept = default;

  /// @brief Destructor.
  ~UqigAvailIndexer() = default;

  /// @brief Copy assignment operator.
  UqigAvailIndexer& operator=(const UqigAvailIndexer& other) = default;

  /// @brief Move assignment operator.
  UqigAvailIndexer& operator=(UqigAvailIndexer&& other) noexcept = default;

  /// @brief Inequality operator.
  /// @param[in] other The indexer to compare against.
  /// @return true if the mappings differ, false otherwise.
  bool operator!=(const UqigAvailIndexer& other) const;

  /// @brief Equality operator (negation of inequality).
  /// @param[in] other The indexer to compare against.
  /// @return true if the mappings are identical, false otherwise.
  inline bool operator==(const UqigAvailIndexer& other) const {
    return !(*this != other);
  }


  /// @brief Build the indexer from a list of UqigAvail values.
  /// @param[in] vec_avail Input list of UqigAvail values (order and duplicates allowed).
  /// @details Internally sorts in ascending order and assigns consecutive 0-based indices.
  ///          Duplicate values in the input will result in only one index mapping.
  /// @note After build(), bidirectional conversion via getIndex() and getUqigAvail() is available.
  /// @note Complexity: O(n log n) where n = vec_avail.size().
  void build(const std::vector<UqigAvail>& vec_avail);

  /// @brief Convert a UqigAvail to its corresponding dense index.
  /// @param[in] avail The UqigAvail value to look up.
  /// @return The 0-based dense index, or IndexNotFound if the value is not in the mapping.
  /// @note Complexity: O(1) average (hash lookup).
  Index getIndex(const UqigAvail avail) const noexcept {
    return bimap_UqigAvail_Index.getABorDefault(avail);
  }

  /// @brief Get the total count of unique UqigAvail entries.
  /// @return The number of entries (also the upper bound of dense indices).
  size_t size() const noexcept { return bimap_UqigAvail_Index.size(); }

  /// @brief Get all UqigAvail values as a vector.
  /// @return Vector of UqigAvail values in ascending order.
  std::vector<UqigAvail> get_vec_UqigAvail() const noexcept {
    return bimap_UqigAvail_Index.get_vec_A();
  }

  /// @brief Get all UqigAvail values as a set.
  /// @return Set of UqigAvail values.
  std::set<UqigAvail> get_set_UqigAvail() const noexcept {
    return bimap_UqigAvail_Index.get_set_A();
  }

  /// @brief Get all dense indices as a vector.
  /// @return Vector of dense indices (0, 1, 2, ..., size()-1).
  std::vector<Index> get_vecIndex() const noexcept {
    return bimap_UqigAvail_Index.get_vec_B();
  }

  /// @brief Get all dense indices as a set.
  /// @return Set of dense indices.
  std::set<Index> get_setIndex() const noexcept {
    return bimap_UqigAvail_Index.get_set_B();
  }

  /// @brief Convert a dense index back to its corresponding UqigAvail.
  /// @param[in] index The 0-based dense index.
  /// @return The corresponding UqigAvail, or UqigAvailNotFound if index is out of range.
  /// @note Complexity: O(1) average (hash lookup).
  UqigAvail getUqigAvail(const Index index) const noexcept {
    return bimap_UqigAvail_Index.getBAorDefault(index);
  }

  /// @brief Serialize the indexer to a binary output stream.
  /// @param[in,out] ofs Output file stream (must be opened in binary mode).
  /// @throws std::runtime_error If write operation fails.
  /// @note Format: [int32:count][UqigAvail[0]]...[UqigAvail[count-1]] in index order.
  void save(std::ofstream& ofs) const;

  /// @brief Deserialize the indexer from a binary input stream.
  /// @param[in,out] ifs Input file stream (must be opened in binary mode).
  /// @throws std::runtime_error If read operation fails or element count is negative.
  /// @note Replaces the current mapping with the loaded data.
  /// @note Format: [int32:count][UqigAvail[0]]...[UqigAvail[count-1]] in index order.
  void load(std::ifstream& ifs);
};
