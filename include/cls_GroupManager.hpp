/// @file cls_GroupManager.hpp
/// @brief Group management with fast lookup capabilities
/// @details
/// Manages group (uqig) information with bidirectional maps for efficient
/// detid and uqig_avail based searches.
///
/// ## Typical Workflow
/// 1. Create a GroupManager instance
/// 2. Call `initialize()` to set minimum uqig and uqig_avail values
/// 3. Call `reserveBimapPlus()` to pre-allocate capacity (optional, for performance)
/// 4. Insert GroupInfo entries using `insert_to_Bimaps()`
/// 5. Query groups by uqig, uqig_avail, or (detid, igroup)
/// 6. Modify group availability using `set_is_avail_uqig()` or `set_is_avail_uqig_avail()`
///
/// ## Thread Safety
/// - Read operations: Thread-safe when no writes occur concurrently
/// - Write operations: NOT thread-safe, external synchronization required
/// - OpenMP: Not used in this class
///
/// ## Memory Layout
/// - Uses three bidirectional hash maps (UOBimap) for O(1) lookup:
///   - Uqig ↔ GroupInfo
///   - UqigAvail ↔ Uqig
///   - (Detid, Igroup) ↔ Uqig
/// - Capacity can be pre-allocated to avoid rehashing during bulk inserts
///
/// ## Input/Output
/// - Binary save/load: Use `save()` and `load()` with std::ofstream/std::ifstream
/// - Format: Custom binary format (not portable across architectures)
#pragma once
#include <set>
#include <unordered_map>
#include <vector>
#include <optional>
#include <stdexcept>
#include "st_GroupInfo.hpp"
#include "cls_UOBimap.hpp"
#include "cls_OneToManyUOBimap.hpp"
#include "ns_type_definitions.hpp"

using namespace index_type_definitions;

/// @brief Class to manage group (uqig) information and provide fast lookup based on detid and uqig_avail
///
/// @details
/// GroupManager maintains bidirectional mappings between different group identifiers
/// for efficient lookups in multiple directions. It tracks group metadata (GroupInfo)
/// and provides methods to query, update, and manage group availability status.
///
/// ## Key Responsibilities
/// - Store and retrieve GroupInfo by uqig, uqig_avail, or (detid, igroup)
/// - Track min/max values for uqig and uqig_avail ranges
/// - Manage group availability flags (is_avail)
/// - Provide set/vector accessors for bulk queries
/// - Support binary serialization for persistence
///
/// ## Primary Use Cases
/// - Detector grouping in muon tomography applications
/// - Fast lookup of group metadata by various identifiers
/// - Managing available vs. unavailable groups for analysis
///
/// ## Example Usage
/// @code
/// GroupManager mgr;
/// mgr.initialize(0, 0);  // Start from uqig=0, uqig_avail=0
/// mgr.reserveBimapPlus(100);  // Pre-allocate for 100 groups
///
/// // Insert a group
/// GroupInfo info{.uqig=0, .uqig_avail=0, .detid=1, .igroup=5, .is_avail=true};
/// mgr.insert_to_Bimaps(info);
///
/// // Lookup by uqig
/// const GroupInfo& retrieved = mgr.getInfo(0);
///
/// // Check availability
/// bool avail = mgr.is_avail(0);
///
/// // Get all groups for a detector
/// auto groups = mgr.get_vecGroupInfo_by_detid(1);
/// @endcode
///
/// @note All lookup operations have O(1) average time complexity
/// @note Thread-safe for concurrent reads, but NOT for concurrent writes
class GroupManager {
private:
  //============================================================================
  /// @name type_definitions
  /// @brief Type definitions used internally by GroupManager
  ///@{

  /// @brief Bidirectional map: uqig ↔ GroupInfo
  using MapUqigGroupInfo = UOBimap<
      Uqig, GroupInfo
    , std::hash<Uqig>, std::equal_to<Uqig>
    , GroupInfoHash,    GroupInfoEq>;

  /// @brief Bidirectional map: uqig_avail ↔ uqig
  using MapUqigAvailUqig = UOBimap<UqigAvail, Uqig>;

  /// @brief Bidirectional map: (detid, igroup) ↔ uqig
  using MapDetIgroupUqig = UOBimap<DetIgroup, Uqig, DetIgroupHash, DetIgroupEq>;

  /// @brief One-to-many bidirectional map: Uqig (One) ↔ Uqid (Many)
  using OtmUqigUqid = OneToManyUOBimap<Uqig, Uqid>;
  
  ///@}
  //----------------------------------------------------------------------------

private:
  //============================================================================
  /// @name private_variables
  /// @brief Member variables used internally by GroupManager
  ///@{

  /// @brief Minimum registered uqig value
  Uqig uqig_min_ = UqigNotAssigned;

  /// @brief Maximum registered uqig value
  Uqig uqig_max_ = UqigNotAssigned;

  /// @brief Minimum registered uqig_avail value
  UqigAvail uqig_avail_min_ = UqigAvailNotAssigned;

  /// @brief Maximum registered uqig_avail value
  UqigAvail uqig_avail_max_ = UqigAvailNotAssigned;

  /// @brief Main bidirectional map: Uqig ↔ GroupInfo
  /// @param notFoundAB Value returned when key not found in A→B direction
  /// @param notFoundBA Value returned when key not found in B→A direction
  MapUqigGroupInfo bimap_Uqig_GroupInfo{ UqigNotFound, GroupInfoNotFound };

  /// @brief Main bidirectional map: uqig_avail ↔ uqig
  MapUqigAvailUqig  bimap_UqigAvail_Uqig{ UqigAvailNotFound, UqigNotFound};

  /// @brief Main bidirectional map: (detid, igroup) ↔ uqig
  MapDetIgroupUqig  bimap_DetIgroup_Uqig{ DetIgroupNotFound, UqigNotFound };

  ///@} ----------------------------------------------------------------------------

public:
  //============================================================================
  /// @name constructor_and_operators
  /// @brief Constructors and operators for GroupManager
  ///@{

  /// @brief default constructor
  GroupManager() = default;

  /// @brief copy constructor
  GroupManager(const GroupManager &other) = default;

  /// @brief move constructor
  GroupManager(GroupManager &&other) noexcept = default;

  /// @brief Inequality operator
  bool operator!=(const GroupManager& other) const;

  /// @brief Equality operator (negation of inequality)
  bool operator==(const GroupManager& other) const { return !(*this != other); };

  /// @brief copy assignment operator
  GroupManager& operator=(const GroupManager &other) = default;

  /// @brief move assignment operator
  GroupManager& operator=(GroupManager&& other) noexcept = default;

  ///@} ----------------------------------------------------------------------------

  //============================================================================
  /// @name setter_functions
  /// @brief Methods to set or modify GroupManager contents
  ///@{

  /// @brief Clear all data and set uqig_min_ and uqig_avail_min_ to initial values
  void initialize(const Uqig uqig_min=0, const UqigAvail uqig_avail_min=0);

  /// @brief Insert Uqig and GroupInfo to bimap_Uqig_GroupInfo
  void insert_to_Bimaps( const GroupInfo& info );

  /// @brief Remove Uqig and GroupInfo from bimap_Uqig_GroupInfo
  void erase_from_Bimaps(const Uqig uqig);

  /// @brief Set the value of uqig_min_
  void set_uqig_min(const Uqig uqig_min) { uqig_min_ = uqig_min; };

  /// @brief Set the value of uqig_max_
  void set_uqig_max(const Uqig uqig_max) { uqig_max_ = uqig_max; };

  /// @brief Set the value of uqig_avail_min_
  void set_uqig_avail_min(const UqigAvail uqig_avail_min) { uqig_avail_min_ = uqig_avail_min; };

  /// @brief Set the value of uqig_avail_max_
  void set_uqig_avail_max(const UqigAvail uqig_avail_max) { uqig_avail_max_ = uqig_avail_max; };

  /// @brief Set the is_avail flag of the specified uqig
  void set_is_avail_uqig(const Uqig uqig, const bool tf);

  /// @brief Set the is_avail flag of the specified uqig_avail
  void set_is_avail_uqig_avail(const UqigAvail uqig_avail, const bool tf);

  ///@} ----------------------------------------------------------------------------


  //============================================================================
  /// @name reserve_methods
  /// @brief Internal structure initialization and capacity reservation
  ///@{

  /// @brief Reserve capacity for bimap_Uqig_GroupInfo, bimap_UqigAvail_Uqig, and bimap_DetIgroup_Uqig
  /// @param n_plus Number of additional elements to reserve
  void reserveBimapPlus(const size_t n_plus) {
    rsvplusBimapUqigGroupInfo(n_plus);
    rsvplusBimapUqigAvailUqig(n_plus);
    rsvplusBimapDetIgroupUqig(n_plus);
  };

  /// @brief Reserve capacity for bimap_Uqig_GroupInfo (current size + n_plus)
  void rsvplusBimapUqigGroupInfo(const size_t n_plus) {
    bimap_Uqig_GroupInfo.reserveAdditional(n_plus);
  };

  /// @brief Reserve capacity for bimap_UqigAvail_Uqig (current size + n_plus)
  void rsvplusBimapUqigAvailUqig(const size_t n_plus) {
    bimap_UqigAvail_Uqig.reserveAdditional(n_plus);
  };

  /// @brief Reserve capacity for bimap_DetIgroup_Uqig (current size + n_plus)
  void rsvplusBimapDetIgroupUqig(const size_t n_plus) {
    bimap_DetIgroup_Uqig.reserveAdditional(n_plus);
  };

  ///@}
  //----------------------------------------------------------------------------

  //============================================================================
  /// @name getter_functions
  /// @brief Getter methods for GroupManager
  ///@{

  /// @brief Get the minimum registered uqig value
  Uqig get_uqig_min() const { return uqig_min_; }

  /// @brief Get the maximum registered uqig value
  Uqig get_uqig_max() const { return uqig_max_; }

  /// @brief Get the minimum registered uqig_avail value
  UqigAvail get_uqig_avail_min() const { return uqig_avail_min_; }

  /// @brief Get the maximum registered uqig_avail value
  UqigAvail get_uqig_avail_max() const { return uqig_avail_max_; }

  /// @brief Get the range of registered uqig values (min and max)
  /// @return std::array<Uqig,2> {min, max}
  std::array<Uqig,2> getUqigRange() const;

  /// @brief Get the range of registered uqig_avail values (min and max)
  /// @return std::array<UqigAvail,2> {min, max}
  std::array<UqigAvail,2> getUqigAvailRange() const;

  /// @brief Get the GroupInfo corresponding to the specified Uqig
  const GroupInfo& getGroupInfo(const Uqig uqig) const;

  /// @brief Get a mutable reference to GroupInfo corresponding to the specified Uqig
  /// @throws std::runtime_error If not found
  GroupInfo& callInfo_by_uqig(const Uqig uqig);

  /// @brief Get a mutable reference to GroupInfo corresponding to the specified UqigAvail
  /// @throws std::runtime_error If not found
  GroupInfo& callInfo_by_uqigAvail(const UqigAvail uqigAvail);

  /// @brief Get all registered uqig values as a vector
  std::vector<Uqig> get_vecUqig_all() const;

  /// @brief Get all registered uqig values as a set
  std::set<Uqig>    get_setUqig_all() const;

  /// @brief Get all registered uqigAvail values as a vector
  std::vector<UqigAvail> get_vecUqigAvail_all() const;

  /// @brief Get all registered uqigAvail values as a set
  std::set<UqigAvail>    get_setUqigAvail_all() const;

  /// @brief Get all uqig values belonging to the specified detid as a vector
  std::vector<Uqig> get_vecUqig_by_detid(const Detid detid) const;

  /// @brief Get all uqig values belonging to the specified detid as a set
  std::set<Uqig>    get_setUqig_by_detid(const Detid detid) const;

  /// @brief Get all uqigAvail values belonging to the specified detid as a vector
  std::vector<UqigAvail> get_vecUqigAvail_by_detid(const Detid detid) const;

  /// @brief Get all uqigAvail values belonging to the specified detid as a set
  std::set<UqigAvail>    get_setUqigAvail_by_detid(const Detid detid) const;

  /// @brief Get the GroupInfo corresponding to the specified uqig
  /// @param uqig Target group index
  /// @return Reference to GroupInfo
  /// @note This is an alias for getInfo_by_uqig()
  const GroupInfo& getInfo(const Uqig uqig) const;

  /// @brief Get the GroupInfo corresponding to the specified Uqig
  /// @throws std::runtime_error If not found
  const GroupInfo& getInfo_by_uqig(const Uqig uqig) const;

  /// @brief Get the GroupInfo corresponding to the specified UqigAvail
  /// @throws std::runtime_error If not found
  const GroupInfo& getInfo_by_uqigAvail(const UqigAvail uqigAvail) const;

  /// @brief Search for uqig from uqig_avail
  /// @param uqig_avail Analysis group ID to search
  /// @return uqig if found, nullopt if not found
  std::optional<Uqig> find_uqig_by_avail(const UqigAvail uqig_avail) const;

  /// @brief Get all registered GroupInfo
  /// @return Vector of GroupInfo
  std::vector<GroupInfo> get_vecGroupInfo_all() const;

  /// @brief Get list of groups corresponding to the specified detid
  /// @param detid Target detector ID
  /// @return Vector of GroupInfo
  std::vector<GroupInfo> get_vecGroupInfo_by_detid(const Detid detid) const;

  /// @brief Get list of groups corresponding to the specified uqig_avail
  /// @param uqig_avail Target uqig_avail to search
  std::vector<GroupInfo> get_vecGroupInfo_by_avail(const UqigAvail uqig_avail) const;

  /// @brief Get the number of unique igroups registered
  int get_n_group( const Detid detid) const;

  /// @brief Get the total number of registered UqigAvail
  int get_n_uqig_avail() const;

  /// @brief Get a list of all unique igroups for the specified detid
  std::vector<Igroup> get_vec_igroup( const Detid detid ) const;

  /// @brief Get a list of all unique igroups (across all detids)
  std::vector<Igroup> get_vec_igroup_all() const;

  /// @brief Get a set of unique igroups for the specified detid
  std::set<Igroup> get_set_igroup(const Detid detid ) const;

  /// @brief Get a set of all unique igroups (across all detids)
  std::set<Igroup> get_set_igroup_all() const;

  /// @brief Get the vector of UqigAvail registered within the specified detid
  std::vector<UqigAvail> get_vec_uqigAvail(const Detid detid) const;

  /// @brief Get the set of UqigAvail registered within the specified detid
  std::set<UqigAvail> get_set_uqigAvail(const Detid detid) const;

  /// @brief Get the number of registered UqigAvail within the specified detid
  int get_n_uqigAvail(const Detid detid) const;

  /// @brief Get the maximum igroup within the specified detid
  Igroup get_igroup_max(const Detid detid) const;

  /// @brief Get uqig from (detid, igroup) pair
  Uqig get_uqig_by_detid_igroup( const Detid detid, const Igroup igroup) const;

  /// @brief Check if this uqig_avail group is enabled (on/off)
  bool is_avail( const Uqig uqig) const;

  ///@} ----------------------------------------------------------------------------

  /// @brief Remove a group by (detid, igroup) pair
  /// @param detid Detector ID
  /// @param igroup Group ID
  /// @throws std::runtime_error If the target does not exist
  void remove_by_detid_igroup(const Detid detid, const Igroup igroup);

  //============================================================================
  /// @name checker methods
  ///@{

  /// @brief Check the consistency of all maps
  bool isConsistent() const {
    return bimap_Uqig_GroupInfo.isConsistent()
        && bimap_UqigAvail_Uqig.isConsistent()
        && bimap_DetIgroup_Uqig.isConsistent();
  };

  ///@} ------------------------------------------------------------------------


  /// @brief for binary save
  void save(std::ofstream& ofs) const;

  /// @brief for binary load
  void load(std::ifstream& ifs);  
};
