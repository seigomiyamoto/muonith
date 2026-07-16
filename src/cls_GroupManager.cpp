// cls_GroupManager.cpp
#include "cls_GroupManager.hpp"
#include "ns_io_binary.hpp"
#include "ns_mymacro.hpp"

bool GroupManager::operator!=(const GroupManager& other) const
{
#ifdef NODEBUG
  if (uqig_min_ != other.uqig_min_) return true;
  if (uqig_max_ != other.uqig_max_) return true;
  if (uqig_avail_min_ != other.uqig_avail_min_) return true;
  if (uqig_avail_max_ != other.uqig_avail_max_) return true;
  if (bimap_Uqig_GroupInfo != other.bimap_Uqig_GroupInfo) return true;
  if (bimap_UqigAvail_Uqig != other.bimap_UqigAvail_Uqig) return true;
  if (bimap_DetIgroup_Uqig != other.bimap_DetIgroup_Uqig) return true;
#else
  if (uqig_min_ != other.uqig_min_) { LOG_WARN("GroupManager: uqig_min_ differs"); return true; }
  if (uqig_max_ != other.uqig_max_) { LOG_WARN("GroupManager: uqig_max_ differs"); return true; }
  if (uqig_avail_min_ != other.uqig_avail_min_) { LOG_WARN("GroupManager: uqig_avail_min_ differs"); return true; }
  if (uqig_avail_max_ != other.uqig_avail_max_) { LOG_WARN("GroupManager: uqig_avail_max_ differs"); return true; }
  if (bimap_Uqig_GroupInfo != other.bimap_Uqig_GroupInfo) { LOG_WARN("GroupManager: bimap_Uqig_GroupInfo differs"); return true; }
  if (bimap_UqigAvail_Uqig != other.bimap_UqigAvail_Uqig) { LOG_WARN("GroupManager: bimap_UqigAvail_Uqig differs"); return true; }
  if (bimap_DetIgroup_Uqig != other.bimap_DetIgroup_Uqig) { LOG_WARN("GroupManager: bimap_DetIgroup_Uqig differs"); return true; }
#endif
  return false;
}

/// @brief Clear all data and set uqig_min_ and uqig_avail_min_ to initial values
void GroupManager::initialize(const Uqig uqig_min, const UqigAvail uqig_avail_min)
{
  // Set the minimum values
  uqig_min_ = uqig_min;
  uqig_max_ = uqig_min; // Initial value is the same
  uqig_avail_min_ = uqig_avail_min;
  uqig_avail_max_ = uqig_avail_min; // Initial value is the same

  // Clear all maps
  bimap_Uqig_GroupInfo.clear();
  bimap_UqigAvail_Uqig.clear();
  bimap_DetIgroup_Uqig.clear();
}

/// @brief insert Uqig and GroupInfo to bimap_Uqig_GroupInfo
void GroupManager::insert_to_Bimaps( const GroupInfo& info )
{
  const Uqig uqig = info.uqig;
  const UqigAvail uqig_avail = info.uqig_avail;
  const Detid detid = info.detid;
  const Igroup igroup = info.igroup;

  bimap_Uqig_GroupInfo.insert(uqig, info);
  if( uqig_avail != UqigAvailNotAssigned ) 
    bimap_UqigAvail_Uqig.insert(uqig_avail, uqig);
  bimap_DetIgroup_Uqig.insert({detid, igroup}, uqig);
}

/// @brief remove Uqig from three bimaps
void GroupManager::erase_from_Bimaps(const Uqig uqig)
{
  bimap_Uqig_GroupInfo.eraseA(uqig);
  bimap_UqigAvail_Uqig.eraseB(uqig);
  bimap_DetIgroup_Uqig.eraseB(uqig);
}

/// @brief Set the is_avail flag of the specified uqig
void GroupManager::set_is_avail_uqig(const Uqig uqig, const bool tf)
{
  bimap_Uqig_GroupInfo.callAB(uqig).is_avail = tf;
}

/// @brief Set the is_avail flag of the specified uqig_avail
void GroupManager::set_is_avail_uqig_avail(const UqigAvail uqig_avail, const bool tf)
{
  const Uqig uqig = bimap_UqigAvail_Uqig.getAB(uqig_avail);
  bimap_Uqig_GroupInfo.callAB(uqig).is_avail = tf;
}

/// @brief Check if this uqig group is enabled (on/off)
bool GroupManager::is_avail(const Uqig uqig) const
{
  return bimap_Uqig_GroupInfo.getAB(uqig).is_avail;
}

std::array<Uqig,2> GroupManager::getUqigRange() const
{
  return {uqig_min_, uqig_max_};
}

std::array<UqigAvail,2> GroupManager::getUqigAvailRange() const
{
  return {uqig_avail_min_, uqig_avail_max_};
}

const GroupInfo& GroupManager::getInfo(const Uqig uqig) const
{
  return bimap_Uqig_GroupInfo.getAB(uqig);
}

/// @brief Search for uqig from uqig_avail
/// @param uqig_avail Analysis group ID to search
/// @return uqig if found, nullopt if not found
std::optional<Uqig> GroupManager::find_uqig_by_avail(const UqigAvail uqig_avail) const
{
  if (bimap_UqigAvail_Uqig.hasA(uqig_avail)) {
    return bimap_UqigAvail_Uqig.getAB(uqig_avail);
  }
  return std::nullopt;
}

std::vector<GroupInfo> GroupManager::get_vecGroupInfo_all() const
{
  std::vector<GroupInfo> result;
  result.reserve(bimap_Uqig_GroupInfo.size());
  for (const auto& [uqig, ginfo] : bimap_Uqig_GroupInfo.getMapAB()) {
    result.push_back(ginfo);
  }
  return result;
}

std::vector<GroupInfo> GroupManager::get_vecGroupInfo_by_detid(const Detid detid) const
{
  std::vector<GroupInfo> result;
  for (const auto& [uqig, group] : bimap_Uqig_GroupInfo.getMapAB()) {
    if (group.detid == detid) {
      result.push_back(group);
    }
  }
  return result;
}

std::vector<GroupInfo> GroupManager::get_vecGroupInfo_by_avail(const UqigAvail uqig_avail) const
{
  std::vector<GroupInfo> result;
  for (const auto& [uqig, group] : bimap_Uqig_GroupInfo.getMapAB()) {
    if (group.uqig_avail == uqig_avail) {
      result.push_back(group);
    }
  }
  return result;
}
/// @brief Get a list of unique igroups filtered by detid
std::vector<Igroup> GroupManager::get_vec_igroup(const Detid detid) const
{
  std::vector<Igroup> result;
  for (const auto& [uqig, ginfo] : bimap_Uqig_GroupInfo.getMapAB()) {
    if (ginfo.detid == detid) {
      result.push_back(ginfo.igroup);
    }
  }
  if (result.empty()) {
    LOG_WARN("No igroup found for detid {}", detid);
    return {}; // Return empty vector
  }
  return result;
}

/// @brief Get a set of unique igroups filtered by detid
std::set<Igroup> GroupManager::get_set_igroup(const Detid detid) const
{
  const auto vec_igroup = get_vec_igroup(detid);
  if (vec_igroup.empty()) {
    LOG_WARN("No igroup found for detid {}", detid);
    return {}; // Return empty set
  }
  std::set<Igroup> set_igroup(vec_igroup.begin(), vec_igroup.end());
  return set_igroup;
}

/// @brief Get a list of all unique igroups (across all detids)
std::vector<Igroup> GroupManager::get_vec_igroup_all() const
{
  std::vector<Igroup> result;
  for (const auto& [uqig, ginfo] : bimap_Uqig_GroupInfo.getMapAB()) {
    result.push_back(ginfo.igroup);
  }
  return result;
}

/// @brief Get a set of all unique igroups (across all detids)
std::set<Igroup> GroupManager::get_set_igroup_all() const
{
  const auto vec_igroup = get_vec_igroup_all();
  std::set<Igroup> set_igroup(vec_igroup.begin(), vec_igroup.end());
  return set_igroup;
}

/// @brief Get the number of unique igroups (without duplicates)
int GroupManager::get_n_group(const Detid detid) const
{
  std::set<Igroup> set_igroup = get_set_igroup(detid);
  return static_cast<int>(set_igroup.size());
}

/// @brief Get the vector of UqigAvail registered within the specified detid
std::vector<UqigAvail> GroupManager::get_vec_uqigAvail(const Detid detid) const
{
  std::vector<UqigAvail> vec_uqigAvail;
  for (const auto& [uqig, ginfo] : bimap_Uqig_GroupInfo.getMapAB()) {
    if (ginfo.detid != detid) continue;
    if (ginfo.uqig_avail == UqigAvailNotAssigned) {
      LOG_WARN("uqig_avail not assigned for uqig {}", uqig);
      continue;
    }
    vec_uqigAvail.push_back(ginfo.uqig_avail);
  }
  return vec_uqigAvail;
}

/// @brief Get the set of UqigAvail registered within the specified detid
std::set<UqigAvail> GroupManager::get_set_uqigAvail(const Detid detid) const
{
  std::vector<UqigAvail> vec_uqigAvail = get_vec_uqigAvail(detid);
  if (vec_uqigAvail.empty()) {
    LOG_ERROR("No uqig_avail found for detid={}", detid);
    THROW_ERROR("GroupManager::get_set_uqigAvail: No uqig_avail found for detid={}", detid);
  }
  std::set<UqigAvail> set_uqigAvail(vec_uqigAvail.begin(), vec_uqigAvail.end());
  return set_uqigAvail;
}

/// @brief Get the number of registered UqigAvail within the specified detid
int GroupManager::get_n_uqigAvail(const Detid detid) const
{
  std::set<UqigAvail> set_uqigAvail = get_set_uqigAvail(detid);
  return static_cast<int>(set_uqigAvail.size());
}

/// @brief Get the total number of registered UqigAvail
int GroupManager::get_n_uqig_avail() const
{
  return static_cast<int>(bimap_UqigAvail_Uqig.size());
}

Igroup GroupManager::get_igroup_max(const Detid detid) const
{
  Igroup igroup_map = std::numeric_limits<Igroup>::min();

  const auto& vec_group = get_vecGroupInfo_by_detid(detid);
  for(const auto& gi : vec_group){
    if(gi.igroup > igroup_map){
      igroup_map = gi.igroup;
    }
  }
  return igroup_map;
}

/// @brief Get uqig from (detid, igroup) pair
/// @return uqig if found, UqigNotFound if not found
Uqig GroupManager::get_uqig_by_detid_igroup(
  const Detid detid, const Igroup igroup) const
{
  const DetIgroup key{detid, igroup};
  if (!bimap_DetIgroup_Uqig.hasA(key)) {
    return UqigNotFound;
  }
  return bimap_DetIgroup_Uqig.getAB(key);
}


/// @brief Remove a group by (detid, igroup) pair
/// @param detid Detector ID
/// @param igroup Group ID
/// @throws std::runtime_error If the target does not exist
void GroupManager::remove_by_detid_igroup(const Detid detid, const Igroup igroup)
{
  const DetIgroup key{detid, igroup};
  if (!bimap_DetIgroup_Uqig.hasA(key)) {
    THROW_ERROR("GroupManager::remove_by_detid_igroup: Target group not found. detid={}, igroup={}", detid, igroup);
  }
  const Uqig target_uqig = bimap_DetIgroup_Uqig.getAB(key);

  bimap_Uqig_GroupInfo.eraseA(target_uqig);
  bimap_UqigAvail_Uqig.eraseB(target_uqig);
  bimap_DetIgroup_Uqig.eraseB(target_uqig);
}

/// @brief Get the GroupInfo corresponding to the specified Uqig
const GroupInfo& GroupManager::getInfo_by_uqig(const Uqig uqig) const
{
  if (!bimap_Uqig_GroupInfo.hasA(uqig)) {
    LOG_ERROR("uqig={} not found", uqig);
    THROW_ERROR("GroupManager::getInfo_by_uqig: uqig={} not found", uqig);
  }
  return bimap_Uqig_GroupInfo.getAB(uqig);
}

/// @brief Get the GroupInfo corresponding to the specified UqigAvail
const GroupInfo& GroupManager::getInfo_by_uqigAvail(const UqigAvail uqigAvail) const
{
  if (!bimap_UqigAvail_Uqig.hasA(uqigAvail)) {
    LOG_ERROR("uqigAvail={} not found", uqigAvail);
    THROW_ERROR("GroupManager::getInfo_by_uqigAvail: uqigAvail={} not found", uqigAvail);
  }
  const Uqig uqig = bimap_UqigAvail_Uqig.getAB(uqigAvail);

  if (!bimap_Uqig_GroupInfo.hasA(uqig)) {
    LOG_ERROR("uqig={} (from uqigAvail={}) not found", uqig, uqigAvail);
    THROW_ERROR("GroupManager::getInfo_by_uqigAvail: uqig={} (from uqigAvail={}) not found", uqig, uqigAvail);
  }
  return bimap_Uqig_GroupInfo.getAB(uqig);
}

/// @brief Get a mutable reference to GroupInfo corresponding to the specified Uqig
/// @throws std::runtime_error If not found
GroupInfo& GroupManager::callInfo_by_uqig(const Uqig uqig)
{
  if (!bimap_Uqig_GroupInfo.hasA(uqig)) {
    LOG_ERROR("uqig={} not found", uqig);
    THROW_ERROR("GroupManager::callInfo_by_uqig: uqig={} not found", uqig);
  }
  return bimap_Uqig_GroupInfo.callAB(uqig);
}

/// @brief Get a mutable reference to GroupInfo corresponding to the specified UqigAvail
/// @throws std::runtime_error If not found
GroupInfo& GroupManager::callInfo_by_uqigAvail(const UqigAvail uqigAvail)
{
  if (!bimap_UqigAvail_Uqig.hasA(uqigAvail)) {
    LOG_ERROR("uqigAvail={} not found", uqigAvail);
    THROW_ERROR("GroupManager::callInfo_by_uqigAvail: uqigAvail={} not found", uqigAvail);
  }
  const Uqig uqig = bimap_UqigAvail_Uqig.getAB(uqigAvail);

  if (!bimap_Uqig_GroupInfo.hasA(uqig)) {
    LOG_ERROR("uqig={} (from uqigAvail={}) not found", uqig, uqigAvail);
    THROW_ERROR("GroupManager::callInfo_by_uqigAvail: uqig={} (from uqigAvail={}) not found", uqig, uqigAvail);
  }
  return bimap_Uqig_GroupInfo.callAB(uqig);
}

/// @brief Get all registered uqig values as a vector
std::vector<Uqig> GroupManager::get_vecUqig_all() const
{
  std::vector<Uqig> result;
  result.reserve(bimap_Uqig_GroupInfo.size());
  for (const auto& [uqig, _] : bimap_Uqig_GroupInfo.getMapAB()) {
    result.push_back(uqig);
  }
  return result;
}

/// @brief Get all registered uqig values as a set
std::set<Uqig> GroupManager::get_setUqig_all() const
{
  std::vector<Uqig> vec_uqig = get_vecUqig_all();
  std::set<Uqig> set_uqig(vec_uqig.begin(), vec_uqig.end());
  return set_uqig;
}

/// @brief Get all registered uqigAvail values as a vector
std::vector<UqigAvail> GroupManager::get_vecUqigAvail_all() const
{
  std::vector<UqigAvail> result;
  result.reserve(bimap_UqigAvail_Uqig.size());
  for (const auto& [uqig_avail, uqig] : bimap_UqigAvail_Uqig.getMapAB()) {
    result.push_back(uqig_avail);
  }
  return result;
}

/// @brief Get all registered uqigAvail values as a set
std::set<UqigAvail> GroupManager::get_setUqigAvail_all() const
{
  std::vector<UqigAvail> vec_uqig_avail = get_vecUqigAvail_all();
  std::set<UqigAvail> set_uqig_avail(vec_uqig_avail.begin(), vec_uqig_avail.end());
  return set_uqig_avail;
}

/// @brief Get all uqig values belonging to the specified detid as a vector
std::vector<Uqig> GroupManager::get_vecUqig_by_detid(const Detid detid) const
{
  std::vector<Uqig> result;
  for (const auto& [uqig, group] : bimap_Uqig_GroupInfo.getMapAB()) {
    if (group.detid == detid) {
      result.push_back(uqig);
    }
  }
  return result;
}

/// @brief Get all uqig values belonging to the specified detid as a set
std::set<Uqig> GroupManager::get_setUqig_by_detid(const Detid detid) const
{
  std::vector<Uqig> vec_uqig = get_vecUqig_by_detid(detid);
  std::set<Uqig> set_uqig(vec_uqig.begin(), vec_uqig.end());
  return set_uqig;
}

/// @brief Get all uqigAvail values belonging to the specified detid as a vector
std::vector<UqigAvail> GroupManager::get_vecUqigAvail_by_detid(const Detid detid) const
{
  std::vector<UqigAvail> result;
  for (const auto& [uqig, group] : bimap_Uqig_GroupInfo.getMapAB()) {
    if (group.detid == detid) {
      result.push_back(group.uqig_avail);
    }
  }
  return result;
}

/// @brief Get all uqigAvail values belonging to the specified detid as a set
std::set<UqigAvail> GroupManager::get_setUqigAvail_by_detid(const Detid detid) const
{
  std::vector<UqigAvail> vec_avail = get_vecUqigAvail_by_detid(detid);
  std::set<UqigAvail> set_avail(vec_avail.begin(), vec_avail.end());
  return set_avail;
}


void GroupManager::save(std::ofstream& ofs) const
{
  // Save minimum and maximum values
  io_binary::write_binary(ofs, uqig_min_);
  io_binary::write_binary(ofs, uqig_max_);
  io_binary::write_binary(ofs, uqig_avail_min_);
  io_binary::write_binary(ofs, uqig_avail_max_);

  // Save bimap_Uqig_GroupInfo
  const auto& mapInfo = bimap_Uqig_GroupInfo.getMapAB();
  io_binary::write_binary(ofs, mapInfo.size());
  for (const auto& [uqig, info] : mapInfo) {
    io_binary::write_binary(ofs, uqig);
    io_binary::write_binary(ofs, info.uqig_avail);
    io_binary::write_binary(ofs, info.detid);
    io_binary::write_binary(ofs, info.igroup);
    io_binary::write_bool(ofs, info.is_avail);
  }

  // Save bimap_UqigAvail_Uqig
  const auto& mapAvail = bimap_UqigAvail_Uqig.getMapAB();
  io_binary::write_binary(ofs, mapAvail.size());
  for (const auto& [uqigAvail, uqig] : mapAvail) {
    io_binary::write_binary(ofs, uqigAvail);
    io_binary::write_binary(ofs, uqig);
  }

  if (ofs.fail()) THROW_ERROR("GroupManager::save failed.");
}

void GroupManager::load(std::ifstream& ifs)
{
  // Load minimum and maximum values
  uqig_min_       = io_binary::read_binary<Uqig>(ifs);
  uqig_max_       = io_binary::read_binary<Uqig>(ifs);
  uqig_avail_min_ = io_binary::read_binary<UqigAvail>(ifs);
  uqig_avail_max_ = io_binary::read_binary<UqigAvail>(ifs);

  // Reconstruct bimap_Uqig_GroupInfo and bimap_DetIgroup_Uqig
  const size_t nInfo = io_binary::read_binary<size_t>(ifs);
  bimap_Uqig_GroupInfo.clear();
  bimap_DetIgroup_Uqig.clear();
  for (size_t i = 0; i < nInfo; ++i) {
    Uqig      uqig       = io_binary::read_binary<Uqig>(ifs);
    GroupInfo info;
    info.uqig       = uqig;
    info.uqig_avail = io_binary::read_binary<UqigAvail>(ifs);
    info.detid      = io_binary::read_binary<Detid>(ifs);
    info.igroup     = io_binary::read_binary<Igroup>(ifs);
    info.is_avail   = io_binary::read_bool(ifs);
    bimap_Uqig_GroupInfo.insert(uqig, info);
    bimap_DetIgroup_Uqig.insert({info.detid, info.igroup}, uqig);
  }

  // Reconstruct bimap_UqigAvail_Uqig
  const size_t nAvail = io_binary::read_binary<size_t>(ifs);
  bimap_UqigAvail_Uqig.clear();
  for (size_t i = 0; i < nAvail; ++i) {
    UqigAvail uqigAvail = io_binary::read_binary<UqigAvail>(ifs);
    Uqig      uqig       = io_binary::read_binary<Uqig>(ifs);
    bimap_UqigAvail_Uqig.insert(uqigAvail, uqig);
  }

  if (ifs.fail()) THROW_ERROR("GroupManager::load failed.");
}

