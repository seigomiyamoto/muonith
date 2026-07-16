/// @file cls_UqidManager.hpp
/// @brief Unique detector element ID manager
/// @details Manages UqidInfo mappings and provides fast lookup by various keys (uqid, detid+ixiy, etc.).
///
/// @par Workflow:
/// 1. Construct UqidManager instance
/// 2. Call initialize() with starting uqid value
/// 3. Register detector elements via insert() or insertInit()
/// 4. Query by uqid, detid/ix/iy, or group identifiers
/// 5. Optional: serialize with save()/load() for persistence
///
/// @note Thread-Safety: No. External synchronization required for concurrent access.
/// @note Memory Layout: Uses three bidirectional maps (UOBimap) internally. No specific access pattern optimization required.
#pragma once
#include <vector>
#include <unordered_map>
#include <set>
#include <functional>      // for std::hash, std::equal_to
#include "cls_UOBimap.hpp"
#include "st_UqidInfo.hpp"   // for Uqid, UqidInfo, DetIxiy, DetIxiyHash, DetIxiyEq
#include "ns_type_definitions.hpp"  // for UqidNotFound, DetidNotFound, etc.
#include "ns_mylogger.hpp"  // for LOG_ERROR, LOG_WARN, LOG_DEBUG, LOG_INFO
#include "ns_mymacro.hpp"  // for THROW_ERROR
#include <filesystem>  // for std::filesystem::path
namespace fs = std::filesystem;
using namespace index_type_definitions;


/// @class UqidManager
/// @brief Manages unique detector element IDs (Uqid) and their associated information.
/// @details Provides bidirectional mappings between Uqid, Index, DetIxiy, and UqidInfo.
///          This class is the central registry for all detector element identification
///          in the raytracing engine. Maintains three bidirectional maps for fast lookups:
///          (1) Uqid <-> Index, (2) DetIxiy <-> Uqid, (3) Uqid <-> UqidInfo.
///
/// @par Typical Usage:
/// @code
/// UqidManager mgr;
/// mgr.initialize(0);  // Start uqid from 0
/// mgr.insert(uqid_info);  // Register detector element
/// auto info = mgr.getInfo(uqid);  // Lookup by uqid (O(1) average)
/// @endcode
///
/// @par Invariants:
/// - uqid_min_ <= uqid_max_ for any registered data
/// - All three bimaps remain synchronized (same uqid set)
/// - DetIxiy uniqueness: each (detid, ix, iy) maps to exactly one uqid
///
/// @par Thread-Safety: No
/// @note Not thread-safe. External synchronization required for concurrent access.
/// @note Complexity: Most lookups are O(1) average case (hash-based). Iteration-based queries (e.g., get_vecUqid_by_detid) are O(n).
class UqidManager {
private:
  //============================================================================
  /// @name type_definitions
  /// @brief Internal type definitions for UqidManager
  ///@{

  /// @brief Bidirectional map between Uqid and Index (int supports std::hash)
  using MapUqidIndex   = UOBimap<Uqid, Index>;

  /// @brief Bidirectional map between DetIxiy and Uqid
  /// @note Requires custom hash/equality functors for DetIxiy (array<int,3>)
  using MapDetIxiyUqid = UOBimap< DetIxiy, Uqid,DetIxiyHash, DetIxiyEq >;

  /// @brief Bidirectional map between Uqid and UqidInfo
  using MapUqidInfo = UOBimap<
      Uqid                  ///< Key A: unique ID
    , UqidInfo              ///< Key B: info structure
    , std::hash<Uqid>       ///< Hash for A
    , std::equal_to<Uqid>   ///< Equality for A
    , UqidInfoHash          ///< Hash for B
    , UqidInfoEq            ///< Equality for B
  >;

  ///@}
  //----------------------------------------------------------------------------

private:

  //============================================================================
  /// @name private_members
  /// @brief Internal member variables for UqidManager
  ///@{

  /// @brief Minimum registered uqid
  Uqid uqid_min_ = UqidNotAssigned;

  /// @brief Maximum registered uqid
  Uqid uqid_max_ = UqidNotAssigned;

  /// @brief Bidirectional map between Uqid and Index
  MapUqidIndex   bimap_Uqid_Index{ UqidNotFound, IndexNotFound };

  /// @brief Bidirectional map between DetIxiy and Uqid
  MapDetIxiyUqid bimap_DetIxiy_Uqid{ DetIxiyNotFound, UqidNotFound };

  /// @brief Bidirectional map between Uqid and UqidInfo
  MapUqidInfo    bimap_Uqid_Info{UqidNotFound, UqidInfoNotFound};
  ///@}--------------------------------------------------------------------------

public:
  //============================================================================
  /// @name constructor_and_operators
  /// @brief Constructors and operators for UqidManager
  ///@{

  /// @brief Default constructor
  UqidManager() = default;

  /// @brief Copy constructor
  UqidManager(const UqidManager&) = default;

  /// @brief Move constructor
  UqidManager(UqidManager&&) noexcept = default;

  /// @brief Destructor
  ~UqidManager() = default;

  /// @brief Copy assignment operator
  UqidManager& operator=(const UqidManager&) = default;

  /// @brief Move assignment operator
  UqidManager& operator=(UqidManager&&) noexcept = default;

  /// @brief Inequality operator
  /// @param[in] other Another UqidManager to compare
  /// @return true if not equal, false otherwise
  bool operator!=(const UqidManager& other) const;

  /// @brief Equality operator (negation of inequality)
  /// @param[in] other Another UqidManager to compare
  /// @return true if equal, false otherwise
  inline bool operator==(const UqidManager& other) const {
    return !(*this != other); }

  /// @brief Initialize all members
  /// @param[in] uqid_min Initial minimum uqid value
  void initialize(const Uqid uqid_min );

  /// @brief Reserve bimap capacity for current size + n_plus elements
  /// @param[in] n_plus Additional capacity to reserve
  void reserveAdditional(const size_t n_plus);
  ///@}
  //----------------------------------------------------------------------------

  //============================================================================
  /// @name getter_functions
  /// @brief Getter functions for UqidManager
  ///@{

  /// @brief Get the minimum and maximum registered uqid
  /// @return Minimum uqid value
  Uqid get_uqid_min() const { return uqid_min_; }
  /// @brief Get the maximum registered uqid
  /// @return Maximum uqid value
  Uqid get_uqid_max() const { return uqid_max_; }

  /// @brief Get the range of registered uqig_avail values
  /// @return std::array<UqigAvail,2> {min, max}
  /// @throws std::runtime_error If no data is registered or no valid uqig_avail exists
  std::array<UqigAvail,2> getUqigAvailRange() const;

  /// @brief Get UqidInfo for a given uqid
  /// @param[in] uqid Target uqid
  /// @return Corresponding UqidInfo (by value)
  /// @throws std::runtime_error If uqid is not found
  UqidInfo getInfo(const Uqid uqid) const;

  /// @brief Get UqidInfo for a given detid, ix, iy
  /// @param[in] detid Detector ID
  /// @param[in] ix X index
  /// @param[in] iy Y index
  /// @return Corresponding UqidInfo (by value)
  /// @throws std::runtime_error If DetIxiy is not found
  UqidInfo getInfo(const Detid detid, const int ix, const int iy) const;

  /// @brief Get DetIxiy (detid, ix, iy) from uqid
  /// @param[in] uqid Target uqid
  /// @return Corresponding DetIxiy array
  /// @throws std::runtime_error If uqid is not found
  DetIxiy getDetIxiy(const Uqid uqid) const;

  /// @brief Get mutable reference to UqidInfo for a given uqid
  /// @param[in] uqid Target uqid
  /// @return Mutable reference to corresponding UqidInfo
  /// @throws std::runtime_error If uqid is not found
  UqidInfo& callInfo(const Uqid uqid);

  /// @brief Get mutable reference to UqidInfo for a given detid, ix, iy
  /// @param[in] detid Detector ID
  /// @param[in] ix X index
  /// @param[in] iy Y index
  /// @return Mutable reference to corresponding UqidInfo
  /// @throws std::runtime_error If DetIxiy is not found
  UqidInfo& callInfo(const Detid detid, const int ix, const int iy);

  /// @brief Get all registered uqids as a vector
  /// @return Vector of all uqids
  std::vector<Uqid> get_vecUqid_all() const;

  /// @brief Get all registered uqids as a sorted set (unique, ascending)
  /// @return Set of all uqids
  std::set<Uqid> get_setUqid_all() const;

  /// @brief Get all uqids belonging to a specific detid
  /// @param[in] detid Target detector ID
  /// @return Vector of uqids belonging to detid
  std::vector<Uqid> get_vecUqid_by_detid(const Detid detid) const;

  /// @brief Get all uqids belonging to a specific detid as a sorted set
  /// @param[in] detid Target detector ID
  /// @return Set of uqids (unique, ascending)
  std::set<Uqid> get_setUqid_by_detid(const Detid detid) const;

  /// @brief Get all uqids belonging to a specific uqig
  /// @param[in] uqig Target unique group ID
  /// @return Vector of uqids belonging to uqig
  std::vector<Uqid> get_vecUqid_by_uqig(const Uqig uqig) const;

  /// @brief Get all uqids belonging to a specific uqig as a sorted set
  /// @param[in] uqig Target unique group ID
  /// @return Set of uqids (unique, ascending)
  std::set<Uqid> get_setUqid_by_uqig(const Uqig uqig) const;

  /// @brief Get all uqids belonging to a specific uqig_avail
  /// @param[in] uqig_avail Target available unique group ID
  /// @return Vector of uqids belonging to uqig_avail
  std::vector<Uqid> get_vecUqid_by_uqigAvail(const UqigAvail uqig_avail) const;

  /// @brief Get all uqids belonging to a specific uqig_avail as a sorted set
  /// @param[in] uqig_avail Target available unique group ID
  /// @return Set of uqids (unique, ascending)
  std::set<Uqid> get_setUqid_by_uqigAvail(const UqigAvail uqig_avail) const;

  /// @brief Get available uqig_avail values for a specific detid
  /// @param[in] detid Target detector ID
  /// @return Vector of available uqig_avail values
  std::vector<UqigAvail> get_vecAvail_by_detid(Detid detid) const;

  /// @brief Get available uqig_avail values for a specific detid as a sorted set
  /// @param[in] detid Target detector ID
  /// @return Set of available uqig_avail values
  std::set<UqigAvail> get_setAvails_by_detid(Detid detid) const;

  /// @brief Get Index by Uqid (const reference)
  /// @param[in] uqid Target uqid
  /// @return Const reference to corresponding Index
  /// @throws std::runtime_error If uqid is not found
  const Index& getIndexByUqid(const Uqid& uqid) const;

  /// @brief Get Uqid by Index (const reference)
  /// @param[in] index Target index
  /// @return Const reference to corresponding Uqid
  /// @throws std::runtime_error If index is not found
  const Uqid& getUqidByIndex(const Index& index) const;

  /// @brief Get Index by Uqid (mutable reference)
  /// @param[in] uqid Target uqid
  /// @return Mutable reference to corresponding Index
  /// @throws std::runtime_error If uqid is not found
  Index& callIndexByUqid(const Uqid& uqid);

  /// @brief Get Uqid by Index (mutable reference)
  /// @param[in] index Target index
  /// @return Mutable reference to corresponding Uqid
  /// @throws std::runtime_error If index is not found
  Uqid& callUqidByIndex(const Index& index);

  /// @brief Get Uqid by DetIxiy (const reference)
  /// @param[in] detixiy Target DetIxiy (detid, ix, iy)
  /// @return Const reference to corresponding Uqid
  /// @throws std::runtime_error If DetIxiy is not found
  const Uqid& getUqidByDetIxiy(const DetIxiy& detixiy) const;

  /// @brief Get Uqid by DetIxiy (mutable reference)
  /// @param[in] detixiy Target DetIxiy (detid, ix, iy)
  /// @return Mutable reference to corresponding Uqid
  /// @throws std::runtime_error If DetIxiy is not found
  Uqid& callUqidByDetIxiy(const DetIxiy& detixiy);

  /// @brief Count uqids belonging to a specific detid and igroup
  /// @param[in] detid_in Target detector ID
  /// @param[in] igroup Target group index
  /// @return Number of uqids matching the criteria
  int get_n_uqid_by_detid_igroup(const Detid detid_in, const Igroup igroup) const;

  /// @brief Check if a uqid is available (enabled)
  /// @param[in] uqid Target uqid
  /// @return true if available, false otherwise
  bool is_avail(const Uqid uqid) const { return getInfo(uqid).is_avail; }

  ///@} ---------------------------------------------------------------------------

  //============================================================================
  /// @name setter_functions
  /// @brief Setter and modification functions for UqidManager
  ///@{

  /// @brief Register UqidInfo to all bimaps (Uqid_Index, DetIxiy_Uqid, Uqid_Info)
  /// @param[in] uqid_info UqidInfo to register
  void insert( const UqidInfo & uqid_info );

  /// @brief Register UqidInfo allowing duplicate initial keys
  /// @param[in] uqid_info UqidInfo to register (allows IndexNotAssigned duplicates)
  void insertInit( const UqidInfo & uqid_info );

  /// @brief Set the minimum uqid value
  /// @param[in] min Minimum uqid value to set
  void set_uqid_min(Uqid min) { uqid_min_ = min; }
  /// @brief Set the maximum uqid value
  /// @param[in] max Maximum uqid value to set
  void set_uqid_max(Uqid max) { uqid_max_ = max; }

  /// @brief Set availability flag for a uqid
  /// @param[in] uqid Target uqid
  /// @param[in] tf true to enable, false to disable
  void set_is_avail(const Uqid uqid, const bool tf) {
    callInfo(uqid).is_avail = tf;
  }

  /// @brief Register a Uqid-Index mapping
  /// @param[in] uqid Uqid to register
  /// @param[in] index Index to associate
  void regUqidIndex(const Uqid& uqid, const Index& index){
    bimap_Uqid_Index.insert(uqid, index);
  }

  /// @brief Remove a uqid and all associated information
  /// @param[in] uqid Target uqid to remove
  /// @throws std::runtime_error If uqid is out of range or not found
  void remove(const Uqid uqid);

  /// @brief Remove Uqid->Index mapping by uqid
  /// @param[in] uqid Uqid to remove from mapping
  void rmUqidIndex(const Uqid& uqid){ bimap_Uqid_Index.eraseA(uqid); }

  /// @brief Remove Index->Uqid mapping by index
  /// @param[in] index Index to remove from mapping
  void rmIndexUqid(const Index& index){ bimap_Uqid_Index.eraseB(index); }

  /// @brief Register a DetIxiy->Uqid mapping
  /// @param[in] detixiy DetIxiy key (detid, ix, iy)
  /// @param[in] uqid Uqid value to associate
  void regDetIxiyUqid(const DetIxiy& detixiy, const Uqid& uqid){
    bimap_DetIxiy_Uqid.insert(detixiy, uqid);
  }

  /// @brief Remove DetIxiy->Uqid mapping by DetIxiy
  /// @param[in] detixiy DetIxiy to remove from mapping
  void rmDetIxiyUqid(const DetIxiy& detixiy){ bimap_DetIxiy_Uqid.eraseA(detixiy); }

  /// @brief Remove Uqid->DetIxiy mapping by uqid
  /// @param[in] uqid Uqid to remove from mapping
  void rmUqidDetIxiy(const Uqid& uqid){ bimap_DetIxiy_Uqid.eraseB(uqid); }

  ///@} ----------------------------------------------------------------------------

  //============================================================================
  /// @name io_functions
  /// @brief I/O functions for serialization and output
  ///@{

  /// @brief Output all UqidInfo to an ASCII file
  /// @param[in] filepath Output file path
  /// @details Each line contains: uqid, detid, ix, iy, inthis, igroup, uqig,
  ///          uqig_avail, is_avail (tab-separated, sorted by uqid ascending).
  /// @throws std::runtime_error If file cannot be opened
  void out_UqidInfo_all(const fs::path& filepath) const;

  /// @brief Save UqidManager state to binary stream
  /// @param[out] ofs Output file stream (must be opened in binary mode)
  /// @throws std::runtime_error If write fails
  /// @note File format: uqid_min, uqid_max, count, then for each entry: uqid, detid, ix, iy, inthis, igroup, uqig, uqig_avail, is_avail
  /// @note Binary compatibility: Field order must match between save() and load()
  void save(std::ofstream& ofs) const;

  /// @brief Load UqidManager state from binary stream
  /// @param[in] ifs Input file stream (must be opened in binary mode)
  /// @throws std::runtime_error If read fails
  /// @note Clears existing data before loading. Reconstructs all three bimaps from serialized UqidInfo entries.
  /// @note Index reconstruction: index = uqid - uqid_min
  void load(std::ifstream& ifs);

  ///@} ----------------------------------------------------------------------------
};
