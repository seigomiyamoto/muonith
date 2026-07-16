// cls_UqidManager.cpp
#include "cls_UqidManager.hpp"
#include <algorithm>
#include <limits>
#include <stdexcept>
#include "ns_io_binary.hpp"
#include "ns_iodir.hpp"
#include "ns_mylogger.hpp"
#include "ns_mymacro.hpp"

bool UqidManager::operator!=(const UqidManager& other) const
{
#ifdef NODEBUG
  if (uqid_min_ != other.uqid_min_) return true;
  if (uqid_max_ != other.uqid_max_) return true;
  if (bimap_Uqid_Index   != other.bimap_Uqid_Index) return true;
  if (bimap_DetIxiy_Uqid != other.bimap_DetIxiy_Uqid) return true;
  if (bimap_Uqid_Info    != other.bimap_Uqid_Info) return true;
#else
  if (uqid_min_ != other.uqid_min_) { LOG_WARN("UqidManager: uqid_min_ differs"); return true; }
  if (uqid_max_ != other.uqid_max_) { LOG_WARN("UqidManager: uqid_max_ differs"); return true; }
  if (bimap_Uqid_Index   != other.bimap_Uqid_Index) { LOG_WARN("UqidManager: bimap_Uqid_Index differs"); return true; }
  if (bimap_DetIxiy_Uqid != other.bimap_DetIxiy_Uqid) { LOG_WARN("UqidManager: bimap_DetIxiy_Uqid differs"); return true; }
  if (bimap_Uqid_Info    != other.bimap_Uqid_Info) { LOG_WARN("UqidManager: bimap_Uqid_Info differs"); return true; }
#endif
  return false;
}

void UqidManager::initialize(const Uqid uqid_min )
{
  uqid_min_ = uqid_min;
  uqid_max_ = uqid_min;

  bimap_Uqid_Index.clear();
  bimap_DetIxiy_Uqid.clear();
  bimap_Uqid_Info.clear();
}

void UqidManager::reserveAdditional(const size_t n_plus)
{
  bimap_Uqid_Index.reserve(bimap_Uqid_Index.size() + n_plus);
  bimap_DetIxiy_Uqid.reserve(bimap_DetIxiy_Uqid.size() + n_plus);
  bimap_Uqid_Info.reserve(bimap_Uqid_Info.size() + n_plus);
}

std::array<UqigAvail,2> UqidManager::getUqigAvailRange() const
{
  // Check if data is registered
  if (bimap_Uqid_Info.size() == 0) {
    THROW_ERROR("UqidManager::getUqigAvailRange: no data registered");
  }

  UqigAvail min_val = std::numeric_limits<UqigAvail>::max();
  UqigAvail max_val = std::numeric_limits<UqigAvail>::min();

  // Scan the A->B map of bimap_Uqid_Info
  for (const auto& kv : bimap_Uqid_Info.getMapAB()) {
    const auto& info = kv.second;
    min_val = std::min(min_val, info.uqig_avail);
    max_val = std::max(max_val, info.uqig_avail);
  }

  // Validate the result
  if (min_val == UqigAvailNotAssigned || max_val == UqigAvailNotAssigned ||
      min_val == UqigAvailNotFound || max_val == UqigAvailNotFound) {
    THROW_ERROR("UqidManager::getUqigAvailRange: no valid uqig_avail registered");
  }

  return { min_val, max_val };
}

UqidInfo UqidManager::getInfo(const Uqid uqid) const
{
  if (!bimap_Uqid_Info.hasA(uqid)) {
    THROW_ERROR2("UqidManager::getInfo: uqid not found", uqid);
  }
  return bimap_Uqid_Info.getABorDefault(uqid);
}

UqidInfo& UqidManager::callInfo(const Uqid uqid)
{
  if (!bimap_Uqid_Info.hasA(uqid)) {
    THROW_ERROR2("UqidManager::callInfo: uqid not found", uqid);
  }
  return bimap_Uqid_Info.callAB(uqid);
}

UqidInfo UqidManager::getInfo(const Detid detid, const int ix, const int iy) const
{
  DetIxiy detixiy{ detid, ix, iy };
  if (!bimap_DetIxiy_Uqid.hasA(detixiy)) {
    THROW_ERROR("UqidManager::getInfo: DetIxiy not found. detid={}, ix={}, iy={}", detid, ix, iy);
  }
  const auto& uqid = bimap_DetIxiy_Uqid.getAB(detixiy);
  return getInfo(uqid);
}

UqidInfo& UqidManager::callInfo(const Detid detid, const int ix, const int iy)
{
  DetIxiy detixiy{ detid, ix, iy };
  if (!bimap_DetIxiy_Uqid.hasA(detixiy)) {
    THROW_ERROR("UqidManager::callInfo: DetIxiy not found. detid={}, ix={}, iy={}", detid, ix, iy);
  }
  auto& uqid = bimap_DetIxiy_Uqid.callAB(detixiy);
  return callInfo(uqid);
}

DetIxiy UqidManager::getDetIxiy(const Uqid uqid) const
{
  if (!bimap_DetIxiy_Uqid.hasB(uqid)) {
    THROW_ERROR("UqidManager::getDetIxiy: Uqid not found. uqid={}", uqid);
  }
  return bimap_DetIxiy_Uqid.getBA(uqid);
}

std::vector<Uqid> UqidManager::get_vecUqid_all() const
{
  std::vector<Uqid> result;
  result.reserve(bimap_Uqid_Index.size());

  for (const auto& [uqid, _] : bimap_Uqid_Index.getMapAB()) {
    result.push_back(uqid);
  }

  return result;
}

std::set<Uqid> UqidManager::get_setUqid_all() const
{
  const auto vec = get_vecUqid_all();
  return { vec.begin(), vec.end() };
}

std::vector<Uqid> UqidManager::get_vecUqid_by_detid(const Detid detid) const
{
  std::vector<Uqid> result;

  for (const auto& [uqid, info] : bimap_Uqid_Info.getMapAB()) {
    if (info.detid == detid) {
      result.push_back(uqid);
    }
  }

  return result;
}

std::set<Uqid> UqidManager::get_setUqid_by_detid(const Detid detid) const
{
  const auto vec = get_vecUqid_by_detid(detid);
  return { vec.begin(), vec.end() };
}

std::vector<Uqid> UqidManager::get_vecUqid_by_uqig(const Uqig uqig) const
{
  std::vector<Uqid> result;
  for (const auto& [uqid, info] : bimap_Uqid_Info.getMapAB()) {
    if (info.uqig == uqig) {
      result.push_back(uqid);
    }
  }
  return result;
}

std::set<Uqid> UqidManager::get_setUqid_by_uqig(const Uqig uqig) const
{
  const auto vec = get_vecUqid_by_uqig(uqig);
  return { vec.begin(), vec.end() };
}

std::vector<Uqid> UqidManager::get_vecUqid_by_uqigAvail(const UqigAvail uqig_avail) const
{
  std::vector<Uqid> result;
  for (const auto& [uqid, info] : bimap_Uqid_Info.getMapAB()) {
    if (info.uqig_avail == uqig_avail) {
      result.push_back(uqid);
    }
  }
  return result;
}

std::set<Uqid> UqidManager::get_setUqid_by_uqigAvail(const UqigAvail uqig_avail) const
{
  const auto vec = get_vecUqid_by_uqigAvail(uqig_avail);
  return { vec.begin(), vec.end() };
}

std::vector<UqigAvail> UqidManager::get_vecAvail_by_detid(Detid detid) const
{
  std::vector<UqigAvail> vec;
  for (const auto& [uqid, info] : bimap_Uqid_Info.getMapAB()) {
    if (info.detid == detid && info.is_avail) {
      vec.push_back(info.uqig_avail);
    }
  }
  return vec;
}

std::set<UqigAvail> UqidManager::get_setAvails_by_detid(Detid detid) const
{
  const auto vec = get_vecAvail_by_detid(detid);
  return std::set<UqigAvail>(vec.begin(), vec.end());
}

const Index& UqidManager::getIndexByUqid(const Uqid& uqid) const
{
  if (!bimap_Uqid_Index.hasA(uqid)) {
    THROW_ERROR("UqidManager::getIndexByUqid: uqid not found. uqid={}", uqid);
  }
  return bimap_Uqid_Index.getAB(uqid);
}

const Uqid& UqidManager::getUqidByIndex(const Index& index) const
{
  if (!bimap_Uqid_Index.hasB(index)) {
    THROW_ERROR("UqidManager::getUqidByIndex: index not found. index={}", index);
  }
  return bimap_Uqid_Index.getBA(index);
}

Index& UqidManager::callIndexByUqid(const Uqid& uqid)
{
  if (!bimap_Uqid_Index.hasA(uqid)) {
    THROW_ERROR("UqidManager::callIndexByUqid: uqid not found. uqid={}", uqid);
  }
  return bimap_Uqid_Index.callAB(uqid);
}

Uqid& UqidManager::callUqidByIndex(const Index& index)
{
  if (!bimap_Uqid_Index.hasB(index)) {
    THROW_ERROR("UqidManager::callUqidByIndex: index not found. index={}", index);
  }
  return bimap_Uqid_Index.callBA(index);
}

const Uqid& UqidManager::getUqidByDetIxiy(const DetIxiy& detixiy) const
{
  if (!bimap_DetIxiy_Uqid.hasA(detixiy)) {
    THROW_ERROR("UqidManager::getUqidByDetIxiy: DetIxiy not found. detid={}, ix={}, iy={}",
      detixiy[0], detixiy[1], detixiy[2]);
  }
  return bimap_DetIxiy_Uqid.getAB(detixiy);
}

Uqid& UqidManager::callUqidByDetIxiy(const DetIxiy& detixiy)
{
  if (!bimap_DetIxiy_Uqid.hasA(detixiy)) {
    THROW_ERROR("UqidManager::callUqidByDetIxiy: DetIxiy not found. detid={}, ix={}, iy={}",
      detixiy[0], detixiy[1], detixiy[2]);
  }
  return bimap_DetIxiy_Uqid.callAB(detixiy);
}

int UqidManager::get_n_uqid_by_detid_igroup(const Detid detid_in, const Igroup igroup) const
{
  int count = 0;
  for (const auto& [uqid, info] : bimap_Uqid_Info.getMapAB()) {
    if (info.detid == detid_in && info.igroup == igroup) {
      count++;
    }
  }
  return count;
}

void UqidManager::insert( const UqidInfo & uqid_info )
{
  const Uqid uqid = uqid_info.uqid;
  const Detid detid = uqid_info.detid;
  const int ix = uqid_info.ixiy[0];
  const int iy = uqid_info.ixiy[1];
  const Index index = IndexNotAssigned;
  
  // insert to bimap_Uqid_Index
  bimap_Uqid_Index.insert(uqid, index);

  // insert to bimap_DetIxiy_Uqid
  DetIxiy detixiy = {detid, ix, iy};
  bimap_DetIxiy_Uqid.insert(detixiy, uqid);

  // insert to bimap_Uqid_Info
  bimap_Uqid_Info.insert(uqid, uqid_info);
}

void UqidManager::insertInit( const UqidInfo & uqid_info )
{
  constexpr Index allow_index = IndexNotAssigned;

  const Uqid  uqid = uqid_info.uqid;
  const Detid detid = uqid_info.detid;
  const int   ix = uqid_info.ixiy[0];
  const int   iy = uqid_info.ixiy[1];
  const Index index = IndexNotAssigned;  // Initial value (same as allow_index)

  DetIxiy detixiy = { detid, ix, iy };

  bimap_Uqid_Index.insertWithBException(uqid, index, allow_index);
  bimap_DetIxiy_Uqid.insert(detixiy, uqid);
  bimap_Uqid_Info.insert(uqid, uqid_info);
}

void UqidManager::remove(const Uqid uqid)
{
  // 1) Range check
  if (uqid < uqid_min_ || uqid > uqid_max_) {
    THROW_ERROR("UqidManager::remove: uqid out of range. uqid={}, range=[{}, {}]",
                uqid, uqid_min_, uqid_max_);
  }
  // 2) Remove UqidInfo
  if (!bimap_Uqid_Info.hasA(uqid)) {
    THROW_ERROR("UqidManager::remove: uqid not found. uqid={}", uqid);
  }
  bimap_Uqid_Info.eraseA(uqid);

  // 3) Remove Uqid -> Index mapping
  bimap_Uqid_Index.eraseA(uqid);

  // 4) Remove DetIxiy -> Uqid mapping
  bimap_DetIxiy_Uqid.eraseB(uqid);
}


void UqidManager::out_UqidInfo_all(const fs::path& filepath) const
{
  std::ofstream ofs(filepath);
  if (!ofs) {
    THROW_ERROR2("UqidManager::out_ascii_vec_uqid_info: failed to open file", filepath.string());
  }
  LOG_INFO("writing to {}", filepath.string());

  constexpr int width = 9;

  // Header row
  ofs << "# "
      << std::setw(width-2) << std::right << "uqid"
      << std::setw(width) << std::right << "detid"
      << std::setw(width) << std::right << "ix"
      << std::setw(width) << std::right << "iy"
      << std::setw(width) << std::right << "inthis"
      << std::setw(width) << std::right << "igroup"
      << std::setw(width) << std::right << "uqig"
      << std::setw(width) << std::right << "uqigAvl"
      << std::setw(width) << std::right << "avail"
      << "\n";

  // Copy to vector for sorting
  std::vector<std::pair<Uqid, UqidInfo>> vec_uqid_info(
    bimap_Uqid_Info.getMapAB().begin(),
    bimap_Uqid_Info.getMapAB().end()
  );

  // Sort by uqid (ascending)
  std::sort(vec_uqid_info.begin(), vec_uqid_info.end(),
    [](const auto& lhs, const auto& rhs) {
      return lhs.first < rhs.first;
    }
  );

  for (const auto& [uqid, info] : vec_uqid_info) {
    ofs << std::setw(width) << std::right << info.uqid
        << std::setw(width) << std::right << info.detid
        << std::setw(width) << std::right << info.ixiy[0]
        << std::setw(width) << std::right << info.ixiy[1]
        << std::setw(width) << std::right << info.inthis
        << std::setw(width) << std::right << info.igroup
        << std::setw(width) << std::right << info.uqig
        << std::setw(width) << std::right << info.uqig_avail
        << std::setw(width) << std::right << info.is_avail
        << "\n";
  }
}

void UqidManager::save(std::ofstream& ofs) const
{
  // 1) Range metadata
  io_binary::write_binary(ofs, uqid_min_);
  io_binary::write_binary(ofs, uqid_max_);

  // 2) Number of UqidInfo entries
  const size_t n = bimap_Uqid_Info.size();
  io_binary::write_binary(ofs, n);

  // 3) Write each UqidInfo field
  for (const auto& [uqid, info] : bimap_Uqid_Info.getMapAB()) {
    io_binary::write_binary(ofs, info.uqid);
    io_binary::write_binary(ofs, info.detid);
    io_binary::write_binary(ofs, info.ixiy[0]);
    io_binary::write_binary(ofs, info.ixiy[1]);
    io_binary::write_binary(ofs, info.inthis);
    io_binary::write_binary(ofs, info.igroup);
    io_binary::write_binary(ofs, info.uqig);
    io_binary::write_binary(ofs, info.uqig_avail);
    io_binary::write_bool  (ofs, info.is_avail);
  }

  if (ofs.fail()) {
    THROW_ERROR("UqidManager::save failed.");
  }
}

void UqidManager::load(std::ifstream& ifs)
{
  // 1) Restore range metadata
  uqid_min_ = io_binary::read_binary<Uqid>(ifs);
  uqid_max_ = io_binary::read_binary<Uqid>(ifs);

  // 2) Clear existing data
  bimap_Uqid_Info.clear();
  bimap_Uqid_Index.clear();
  bimap_DetIxiy_Uqid.clear();

  // 3) Read number of UqidInfo entries
  const size_t n = io_binary::read_binary<size_t>(ifs);

  // 4) Reconstruct the three UOBimaps from loaded entries
  for (size_t i = 0; i < n; ++i) {
    UqidInfo info;
    info.uqid       = io_binary::read_binary<Uqid>(ifs);
    info.detid      = io_binary::read_binary<Detid>(ifs);
    {
      int ix = io_binary::read_binary<int>(ifs);
      int iy = io_binary::read_binary<int>(ifs);
      info.ixiy = { ix, iy };
    }
    info.inthis = io_binary::read_binary<Inthis>(ifs);
    info.igroup = io_binary::read_binary<Igroup>(ifs);
    info.uqig       = io_binary::read_binary<Uqig>(ifs);
    info.uqig_avail = io_binary::read_binary<UqigAvail>(ifs);
    info.is_avail   = io_binary::read_bool(ifs);

    // (a) Register UqidInfo
    bimap_Uqid_Info.insert(info.uqid, info);

    // (b) Reconstruct Uqid->Index map (index = uqid - uqid_min_)
    const Index index = static_cast<Index>(info.uqid - uqid_min_);
    bimap_Uqid_Index.insert(info.uqid, index);

    // (c) Reconstruct DetIxiy->Uqid map
    DetIxiy d = { info.detid, info.ixiy[0], info.ixiy[1] };
    bimap_DetIxiy_Uqid.insert(d, info.uqid);
  }

  if (ifs.fail()) {
    THROW_ERROR("UqidManager::load failed.");
  }
}

