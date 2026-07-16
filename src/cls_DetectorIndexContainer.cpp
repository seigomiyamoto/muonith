// cls_DetectorIndexContainer.cpp
#include <algorithm>  // for std::move
#include "cls_DetectorIndexContainer.hpp"
#include "ns_myapp.hpp"
#include "ns_mylogger.hpp"
#include "ns_mymacro.hpp"
#include "ns_io_binary.hpp"

bool DetectorIndexContainer::operator!=(const DetectorIndexContainer& other) const
{
#ifdef NODEBUG
  if (tf_built_uqid != other.tf_built_uqid) return true;
  if (tf_built_uqig != other.tf_built_uqig) return true;
  if (tf_built_uqigAvail != other.tf_built_uqigAvail) return true;
  if (tf_built_indexer != other.tf_built_indexer) return true;
  if (detid_max_ != other.detid_max_) return true;
  if (group_mgr_ != other.group_mgr_) return true;
  if (uqid_mgr_ != other.uqid_mgr_) return true;
  if (avail_indexer_ != other.avail_indexer_) return true;
#else
  if (tf_built_uqid != other.tf_built_uqid) { LOG_WARN("DetectorIndexContainer: tf_built_uqid differs"); return true; }
  if (tf_built_uqig != other.tf_built_uqig) { LOG_WARN("DetectorIndexContainer: tf_built_uqig differs"); return true; }
  if (tf_built_uqigAvail != other.tf_built_uqigAvail) { LOG_WARN("DetectorIndexContainer: tf_built_uqigAvail differs"); return true; }
  if (tf_built_indexer != other.tf_built_indexer) { LOG_WARN("DetectorIndexContainer: tf_built_indexer differs"); return true; }
  if (detid_max_ != other.detid_max_) { LOG_WARN("DetectorIndexContainer: detid_max_ differs"); return true; }
  if (group_mgr_ != other.group_mgr_) { LOG_WARN("DetectorIndexContainer: group_mgr_ differs"); return true; }
  if (uqid_mgr_ != other.uqid_mgr_) { LOG_WARN("DetectorIndexContainer: uqid_mgr_ differs"); return true; }
  if (avail_indexer_ != other.avail_indexer_) { LOG_WARN("DetectorIndexContainer: avail_indexer_ differs"); return true; }
#endif
  return false;
}

void DetectorIndexContainer::initialize()
{
  // Clear all members
  tf_built_uqid = false;
  tf_built_uqig = false;
  tf_built_uqigAvail = false;
  tf_built_indexer = false;
  detid_max_ = DetidNotAssigned;
  group_mgr_ = GroupManager();
  uqid_mgr_ = UqidManager();
  avail_indexer_ = UqigAvailIndexer();
}

std::set<Ixiy> DetectorIndexContainer::get_set_ixiy_by_detid_igroup(
  const Detid detid, const Igroup igroup) const
{
  const auto vec_ixiy = get_vec_ixiy_by_detid_igroup(detid, igroup);
  return std::set<Ixiy>(vec_ixiy.begin(), vec_ixiy.end());
}

DetIxiy DetectorIndexContainer::get_detidixiy( const Uqid uqid_in ) const
{
  // Range check
  if( uqid_in < uqid_mgr_.get_uqid_min() ){
    THROW_ERROR3("uqid_in < uqid_min", uqid_in, uqid_mgr_.get_uqid_min());
  }
  if( uqid_in > uqid_mgr_.get_uqid_max() ){
    THROW_ERROR3("uqid_in > uqid_max", uqid_in, uqid_mgr_.get_uqid_max());
  }

  // Get UqidInfo
  const auto& uqid_info = uqid_mgr_.getInfo(uqid_in);

  // Assemble DetIxiy
  DetIxiy detixiy = {
      uqid_info.detid     // Detid
    , uqid_info.ixiy[0]   // ix
    , uqid_info.ixiy[1]   // iy
  };

  return detixiy;
}

std::vector<Ixiy> DetectorIndexContainer::get_vec_ixiy_by_detid_igroup(
    const Detid detid, const Igroup igroup) const
{
  std::vector<Ixiy> res;
  const std::set<Uqid> set_uqid = get_set_uqid_by_detid_igroup(detid, igroup);
  for (auto u : set_uqid) {
    res.push_back(uqid_mgr_.getInfo(u).ixiy);
  }
  return res;
}

std::set<Uqid> DetectorIndexContainer::get_set_uqid_by_detid_igroup(
  const Detid detid, const Igroup igroup) const
{
  // 1) Get uqig from (detid,igroup)
  const Uqig uqig = group_mgr_.get_uqig_by_detid_igroup(detid, igroup);
  if (uqig==UqigNotFound) return {};

  // 2) Return list of uqid belonging to that uqig
  return uqid_mgr_.get_setUqid_by_uqig(uqig);
}


Igroup DetectorIndexContainer::get_igroup(
  const Detid detid, const int ix, const int iy) const
{
  const UqidInfo& uinfo = uqid_mgr_.getInfo(detid, ix, iy);
  return uinfo.igroup;
}

/// @brief get the set of detid
std::set<Detid> DetectorIndexContainer::get_set_detid() const
{
  std::set<Detid> detid_set;
  const std::vector<GroupInfo> vec_GroupInfo_all = group_mgr_.get_vecGroupInfo_all();
  for (const auto& group : vec_GroupInfo_all) {
    detid_set.insert(group.detid);
  }
  if( detid_set.size() == 0 ){
    LOG_WARN("No detid found in GroupManager");
  }
  return detid_set;
}

std::vector<Index> DetectorIndexContainer::get_enabled_indices_excluding_detid(
  const Detid detid_disabled) const
{
  // 1. Get uqig_avail group belonging to detid_disabled and convert to row indices
  const std::vector<GroupInfo> vec_ginfo = group_mgr_.get_vecGroupInfo_by_detid(detid_disabled);
  std::set<Index> set_disabled_Index;

  for (const auto& group : vec_ginfo) {
    try {
      Index idx = avail_indexer_.getIndex(group.uqig_avail);
      set_disabled_Index.insert(idx);
    } catch (const std::out_of_range& e) {
      LOG_WARN("uqig_avail {} from detid {} not found in indexer: {}",
               group.uqig_avail, group.detid, e.what());
    }
  }

  // 2. Extract indices that are not disabled
  std::vector<Index> vec_keep_Index;
  const int total_size = avail_indexer_.size();
  for (int i = 0; i < total_size; ++i) {
    if (!set_disabled_Index.count(i)) {
      vec_keep_Index.push_back(i);
    }
  }

  return vec_keep_Index;
}

const GroupInfo& DetectorIndexContainer::getGroupInfo(
  const Detid detid, const int ix, const int iy) const
{
  const UqidInfo& uinfo = uqid_mgr_.getInfo(detid, ix, iy);
  const Uqig uqig = uinfo.uqig;
  return getGrpMgr().getInfo_by_uqig(uqig);
}

const GroupInfo& DetectorIndexContainer::getGroupInfo_by_Index(const Index index) const
{
  // Throw exception if Index is out of range
  if (index < 0 || index >= avail_indexer_.size()) {
    LOG_ERROR("Index out of range: {} (max: {})", index, avail_indexer_.size());
    THROW_ERROR("Index out of range");
  }

  // Get corresponding uqig_avail
  const UqigAvail uqig_avail = avail_indexer_.getUqigAvail(index);

  // Get GroupInfo corresponding to uqig_avail
  return group_mgr_.getInfo_by_uqigAvail(uqig_avail);
}


GroupInfo& DetectorIndexContainer::callGroupInfo(const Detid detid, const int ix, const int iy)
{
  const UqidInfo& uinfo = getUqidMgr().getInfo(detid, ix, iy);
  Uqig uqig = uinfo.uqig;
  return callGrpMgr().callInfo_by_uqig(uqig);
}

void DetectorIndexContainer::remove_group(
  const Detid detid, Igroup igroup)
{
  // Delegate group removal to GroupManager
  group_mgr_.remove_by_detid_igroup(detid, igroup);
}

Eigen::MatrixXf DetectorIndexContainer::get_disabled_mat_dNdD(
  const Detid detid, const Eigen::MatrixXf& mat_dNdD) const
{
  const auto keep_rows = get_enabled_indices_excluding_detid(detid);

  Eigen::MatrixXf result(keep_rows.size(), mat_dNdD.cols());
  for (int i = 0; i < result.rows(); ++i) {
    result.row(i) = mat_dNdD.row(keep_rows[i]);
  }
  // Return new matrix excluding rows corresponding to detid
  return result;
}

NmuonVectors DetectorIndexContainer::get_disabled_vec_nmuons(
  const Detid detid, const Eigen::VectorXf& vecxf_nmuon_obs
, const Eigen::VectorXf& vecxf_nmuon_prior) const
{
  const auto keep_indices = get_enabled_indices_excluding_detid(detid);

  Eigen::VectorXf vec_obs(keep_indices.size()), vec_prior(keep_indices.size());

  for (int i = 0; i < keep_indices.size(); ++i) {
    vec_obs(i) = vecxf_nmuon_obs(keep_indices[i]);
    vec_prior(i) = vecxf_nmuon_prior(keep_indices[i]);
  }

  return {
    .vec_obs = vec_obs,
    .vec_prior = vec_prior
  };
}

Eigen::MatrixXf DetectorIndexContainer::get_disabled_mat_cov_muon(
  Detid detid_disabled, const Eigen::MatrixXf& mat_cov_muon) const
{
  // 1. Get indices to keep (this function is already shared)
  const std::vector<Index> vec_keep_Index = get_enabled_indices_excluding_detid(detid_disabled);
  const int n = static_cast<int>(vec_keep_Index.size());

  // 2. Create square matrix for result
  Eigen::MatrixXf result(n, n);

  // 3. Extract blocks corresponding to selected indices (parallelized)
  #pragma omp parallel for collapse(2)
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      result(i, j) = mat_cov_muon(vec_keep_Index[i], vec_keep_Index[j]);
    }
  }

  return result;
}

void DetectorIndexContainer::build_uqigAvailIndexer()
{
  // Warn if already built
  if(is_built_indexer()) LOG_WARN("UqigAvailIndexer is already built.");

  // Check state of uqid_mgr_ and group_mgr_
  check_built_uqid();
  check_built_uqig();
  check_built_uqigAvail();

  std::vector<UqigAvail> vec_uqigAvail = group_mgr_.get_vecUqigAvail_all();
  avail_indexer_.build(vec_uqigAvail);

  // Set flag indicating indexer construction is complete
  set_built_indexer(true);
}

// Build UqigAvailIndexer with specified Detid disabled
// detid_disabled is the detector ID to disable
void DetectorIndexContainer::build_uqigAvailIndexer_disabled(const Detid detid_disabled)
{
  // Warn if already built
  if(is_built_indexer()) LOG_WARN("UqigAvailIndexer is already built.");

  // Check state of uqid_mgr_ and group_mgr_
  check_built_uqid();
  check_built_uqig();
  check_built_uqigAvail();

  // Range check
  if (detid_disabled < 0 || detid_disabled >= get_detid_max()) {
    LOG_ERROR("detid_disabled out of range: {} (max: {})", detid_disabled, get_detid_max());
    THROW_ERROR("detid_disabled out of range");
  }

  // Prepare vector for availindexer
  std::vector<UqigAvail> vec_uqigAvail_filtered;

  // Reserve size for vec_uqigAvail_filtered
  const Detid maxDetid = get_detid_max();
  std::vector<size_t> vec_size(maxDetid, 0);
  std::vector<size_t> vec_offset(maxDetid, 0);
  size_t totalSize = 0;

  for (Detid detid = 0; detid < maxDetid; ++detid) {
    if (detid == detid_disabled) continue;
    vec_size.at(detid) = group_mgr_.get_vecUqigAvail_by_detid(detid).size();
  }

  // Prefix sum
  for (Detid detid = 1; detid < maxDetid; ++detid) {
    vec_offset.at(detid) = vec_offset.at(detid - 1) + vec_size.at(detid - 1);
  }
  totalSize = vec_offset.at(maxDetid - 1) + vec_size.at(maxDetid - 1);

  vec_uqigAvail_filtered.resize(totalSize);

  #pragma omp parallel for
  for (Detid detid = 0; detid < maxDetid; ++detid) {
    if (detid == detid_disabled || vec_size.at(detid) == 0) continue;

    // This is thread safe
    std::vector<UqigAvail> vec = group_mgr_.get_vecUqigAvail_by_detid(detid);

    std::move(
        vec.begin()
      , vec.end()
      , vec_uqigAvail_filtered.begin() + vec_offset.at(detid)
    );
  }

  // Build UqigAvailIndexer
  avail_indexer_.build(vec_uqigAvail_filtered);

  // Set flag indicating indexer construction is complete
  set_built_indexer(true);
}

// Build UqigAvailIndexer with list of disabled Detid (vector version)
void DetectorIndexContainer::build_uqigAvailIndexer_disabled(
  const std::vector<Detid>& vec_detid_disabled)
{
  // Warn if already built
  if(is_built_indexer()) LOG_WARN("UqigAvailIndexer is already built.");

  // Check state of uqid_mgr_ and group_mgr_
  check_built_uqid();
  check_built_uqig();
  check_built_uqigAvail();

  // Get maximum detid
  const Detid maxDetid = get_detid_max();

  // 0) Validate and eliminate duplicates
  std::vector<char> vec_is_disabled(maxDetid, 0);
  for (const Detid detid : vec_detid_disabled) {
    if (detid < 0 || detid >= maxDetid)
      THROW_ERROR("detid_disabled out of range");
    vec_is_disabled.at(detid) = 1;  // Duplicates are automatically eliminated here
  }

  // 1) Collect element count for each detid (even serial O(N) << copy cost)
  std::vector<size_t> vec_size(maxDetid, 0);
  size_t totalSize = 0;
  for (Detid detid = 0; detid < maxDetid; ++detid) {
    if (vec_is_disabled.at(detid)) continue;
    vec_size.at(detid) = group_mgr_.get_vecUqigAvail_by_detid(detid).size();
    totalSize += vec_size.at(detid);
  }

  // 2) Compute offset with exclusive scan (C++20)
  std::vector<size_t> vec_offset(maxDetid, 0);
  std::exclusive_scan(vec_size.begin(), vec_size.end(),
                      vec_offset.begin(), 0);

  // 3) Reserve output buffer
  std::vector<UqigAvail> vec_uqigAvail_filtered(totalSize);

  // 4) Parallel copy: each thread modifies only its own region, so safe
  #pragma omp parallel for
  for (Detid detid = 0; detid < maxDetid; ++detid) {
    if (vec_is_disabled.at(detid) || vec_size.at(detid) == 0) continue;
    auto vec = group_mgr_.get_vecUqigAvail_by_detid(detid);
    std::move(vec.begin(), vec.end(),
              vec_uqigAvail_filtered.begin() + vec_offset.at(detid));
  }

  // 5) Empty check
  if (vec_uqigAvail_filtered.empty()) {
    LOG_ERROR("No uqig_avail found after filtering disabled detids ({})",
      fmt::join(vec_detid_disabled, ", "));
    THROW_ERROR("No valid uqig_avail found");
  }

  // 6) Build indexer
  avail_indexer_.build(vec_uqigAvail_filtered);

  // Set flag indicating indexer construction is complete
  set_built_indexer(true);
}

void DetectorIndexContainer::build_uqigAvailIndexer_enabled(
  const std::vector<Detid>& vec_detid_enabled)
{
  // Warn if already built
  if(is_built_indexer()) LOG_WARN("UqigAvailIndexer is already built.");

  // Check state of uqid_mgr_ and group_mgr_
  check_built_uqid();
  check_built_uqig();
  check_built_uqigAvail();

  // Get maximum detid
  const Detid maxDetid = get_detid_max();

  // 0) Validate and eliminate duplicates
  std::vector<char> vec_is_enabled(maxDetid, 0);
  for (const Detid detid : vec_detid_enabled) {
    if (detid < 0 || detid >= maxDetid)
      THROW_ERROR("detid_enabled out of range");
    vec_is_enabled.at(detid) = 1;  // Duplicates are eliminated here
  }

  // 1) Collect sizes
  std::vector<size_t> vec_size(maxDetid, 0);
  size_t totalSize = 0;
  for (Detid detid = 0; detid < maxDetid; ++detid) {
    if (!vec_is_enabled.at(detid)) continue;
    vec_size.at(detid) = group_mgr_.get_vecUqigAvail_by_detid(detid).size();
    totalSize += vec_size.at(detid);
  }

  // 2) Compute offset with exclusive scan
  std::vector<size_t> vec_offset(maxDetid, 0);
  std::exclusive_scan(vec_size.begin(), vec_size.end(),
                      vec_offset.begin(), 0);

  // 3) Reserve output buffer
  std::vector<UqigAvail> vec_uqigAvail_filtered(totalSize);

  // 4) Parallel copy
  #pragma omp parallel for
  for (Detid detid = 0; detid < maxDetid; ++detid) {
    if (!vec_is_enabled.at(detid) || vec_size.at(detid) == 0) continue;
    auto vec = group_mgr_.get_vecUqigAvail_by_detid(detid);
    std::move(vec.begin(), vec.end(),
              vec_uqigAvail_filtered.begin() + vec_offset.at(detid));
  }

  // 5) Empty check
  if (vec_uqigAvail_filtered.empty()) {
    LOG_ERROR("No uqig_avail found for enabled detids ({})",
      fmt::join(vec_detid_enabled, ", "));
    THROW_ERROR("No valid uqig_avail found");
  }

  // 6) Build indexer
  avail_indexer_.build(vec_uqigAvail_filtered);

  // Set flag indicating indexer construction is complete
  set_built_indexer(true);
}


void DetectorIndexContainer::set(const DetectorIndexContainer &dic_in)
{ 
  group_mgr_ = dic_in.group_mgr_;
  uqid_mgr_ = dic_in.uqid_mgr_;
  avail_indexer_ = dic_in.avail_indexer_;
  tf_built_uqid = dic_in.tf_built_uqid;
  tf_built_uqig = dic_in.tf_built_uqig;
  tf_built_uqigAvail = dic_in.tf_built_uqigAvail;
  tf_built_indexer = dic_in.tf_built_indexer;
}

void DetectorIndexContainer::check_built_uqid() const
{
  if (!tf_built_uqid) {
    THROW_ERROR("DetectorIndexContainer maps about uqid is not build.");
  }
}

void DetectorIndexContainer::check_built_uqig() const
{
  if (!tf_built_uqig) {
    THROW_ERROR("DetectorIndexContainer maps about uqig is not build.");
  }
}

void DetectorIndexContainer::check_built_uqigAvail() const
{
  if (!tf_built_uqigAvail) {
    THROW_ERROR("DetectorIndexContainer maps about uqigAvail is not build.");
  }
}

void DetectorIndexContainer::check_built_indexer() const
{
  if (!tf_built_indexer) {
    THROW_ERROR("DetectorIndexContainer maps about indexer is not built.");
  }
}

void DetectorIndexContainer::out_grouping_info(FILE* fout) const {
  const auto& uqid_mgr = getUqidMgr();
  const auto& group_mgr = getGrpMgr();
  
  for (const Uqid uqid : uqid_mgr.get_vecUqid_all()) {
    const auto& info = uqid_mgr.getInfo(uqid);
    const Igroup igroup = group_mgr.getInfo(info.uqig).igroup;

    fprintf(fout, "%3d %3d %3d %6d %5d %5d %5d\n",
        info.detid, info.ixiy[0], info.ixiy[1]
      , info.uqid, info.uqig, info.uqig_avail, igroup);
  }
}

void DetectorIndexContainer::out_grouping_info(const fs::path& pathout) const
{
  if (pathout.empty()) {
    THROW_ERROR("DetectorIndexContainer::out_grouping_info: path is empty");
  }

  FILE* fout = myapp::get_fout(pathout);
  if (fout == nullptr) {
    THROW_ERROR("DetectorIndexContainer::out_grouping_info: failed to open file: {}", pathout.string());
  }

  out_grouping_info(fout);

  myapp::close(fout, pathout);
}

void DetectorIndexContainer::out_all_index_info(FILE* fout) const
{
  const auto& uqid_mgr = getUqidMgr();
  const auto& group_mgr = getGrpMgr();
  const auto& avail_indexer = getAvailIndexer();

  fprintf(fout,
    "%6s %4s %4s %7s %6s %6s %6s %8s %11s %7s\n",
    "Detid", "Ix", "Iy", "Inthis", "Igroup",
    "Uqid", "Uqig", "is_avail", "UqigAvail", "Index"
  );

  for (const Uqid uqid : uqid_mgr.get_vecUqid_all()) {
    const UqidInfo& info = uqid_mgr.getInfo(uqid);
    const GroupInfo& group_info = group_mgr.getInfo(info.uqig);
    const int index = avail_indexer.getIndex(info.uqig_avail);

    fprintf(fout,
      "%6d %4d %4d %7d %6d %6d %6d %8d %11d %7d\n",
      info.detid,
      info.ixiy[0],
      info.ixiy[1],
      info.inthis,
      group_info.igroup,
      info.uqid,
      info.uqig,
      info.is_avail,
      info.uqig_avail,
      index
    );
  }
}

void DetectorIndexContainer::out_all_index_info(
  const fs::path& pathout) const
{
  if (pathout.empty()) {
    THROW_ERROR("DetectorIndexContainer::out_all_index_info: output file path is empty");
  }

  FILE* fout = std::fopen(pathout.string().c_str(), "wt");
  if (!fout) {
    THROW_ERROR("DetectorIndexContainer::out_all_index_info: failed to open file: {}", pathout.string());
  }

  out_all_index_info(fout);
  std::fclose(fout);
}

void DetectorIndexContainer::out_all_index_info2(
  const fs::path& pathout, const int width ) const
{
  if (pathout.empty()) {
    THROW_ERROR("DetectorIndexContainer::out_all_index_info2: output file path is empty");
  }

  std::ofstream ofs(pathout);
  if (!ofs) {
    THROW_ERROR("DetectorIndexContainer::out_all_index_info2: failed to open file: {}", pathout.string());
  }

  // Header information
  ofs << "detid_max: "       << detid_max_                            << '\n';
  ofs << "built_uqid: "      << (tf_built_uqid    ? "true" : "false") << '\n';
  ofs << "built_uqig: "      << (tf_built_uqig    ? "true" : "false") << '\n';
  ofs << "built_uqigAvail: " << (tf_built_uqigAvail ? "true" : "false") << '\n';
  ofs << "built_indexer: "   << (tf_built_indexer ? "true" : "false") << '\n';

  // Comment line: legend
  ofs << "# "
      << std::setw(width-2) << std::right << "Uqid"       << ' '
      << std::setw(width)   << std::right << "Detid"      << ' '
      << std::setw(width)   << std::right << "Ix"         << ' '
      << std::setw(width)   << std::right << "Iy"         << ' '
      << std::setw(width)   << std::right << "Inthis"     << ' '
      << std::setw(width)   << std::right << "Igroup"     << ' '
      << std::setw(width)   << std::right << "Uqig"       << ' '
      << std::setw(width)   << std::right << "UqigAvail"  << ' '
      << std::setw(width)   << std::right << "is_avail"   << ' '
      << std::setw(width)   << std::right << "Index"      << '\n';

  // Output information for each Uqid
  ofs << std::fixed;
  for (auto uqid : uqid_mgr_.get_setUqid_all()) {
    const auto info = uqid_mgr_.getInfo(uqid);
    auto index = avail_indexer_.getIndex(info.uqig_avail);
    ofs << std::setw(width) << std::right << uqid             << ' '
        << std::setw(width) << std::right << info.detid       << ' '
        << std::setw(width) << std::right << info.ixiy[0]     << ' '
        << std::setw(width) << std::right << info.ixiy[1]     << ' '
        << std::setw(width) << std::right << info.inthis  << ' '
        << std::setw(width) << std::right << info.igroup      << ' '
        << std::setw(width) << std::right << info.uqig        << ' '
        << std::setw(width) << std::right << info.uqig_avail  << ' '
        << std::setw(width) << std::right << info.is_avail    << ' '
        << std::setw(width) << std::right << index            << '\n';
  }
}

void DetectorIndexContainer::out_all_index_info_csv(
  const fs::path& pathout ) const
{
  std::ofstream ofs(pathout);
  if (!ofs) {
    THROW_ERROR("DetectorIndexContainer::out_all_index_info_csv: failed to open file: {}", pathout.string());
  }

  // Header row (comma-separated)
  ofs << "Uqid,Detid,Ix,Iy,Inthis,Igroup,Uqig,UqigAvail,is_avail,Index\n";

  // Output information for each Uqid in CSV format
  for (auto uqid : uqid_mgr_.get_setUqid_all()) {
    const auto info = uqid_mgr_.getInfo(uqid);
    auto index = avail_indexer_.getIndex(info.uqig_avail);
    ofs << uqid             << ','
        << info.detid       << ','
        << info.ixiy[0]     << ','
        << info.ixiy[1]     << ','
        << info.inthis      << ','
        << info.igroup      << ','
        << info.uqig        << ','
        << info.uqig_avail  << ','
        << info.is_avail    << ','
        << index            << '\n';
  }
}


void DetectorIndexContainer::save(std::ofstream& ofs) const
{
  // Build-state flags
  io_binary::write_binary(ofs, tf_built_uqid);
  io_binary::write_binary(ofs, tf_built_uqig);
  io_binary::write_binary(ofs, tf_built_uqigAvail);
  io_binary::write_binary(ofs, tf_built_indexer);
  io_binary::write_binary(ofs, detid_max_);
  // GroupManager
  group_mgr_.save(ofs);
  // UqidManager
  uqid_mgr_.save(ofs);
  // UqigAvailIndexer
  avail_indexer_.save(ofs);
  if (ofs.fail()) THROW_ERROR("DetectorIndexContainer::save failed.");
}

void DetectorIndexContainer::load(std::ifstream& ifs)
{
  // Build-state flags
  tf_built_uqid        = io_binary::read_binary<bool>(ifs);
  tf_built_uqig        = io_binary::read_binary<bool>(ifs);
  tf_built_uqigAvail   = io_binary::read_binary<bool>(ifs);
  tf_built_indexer      = io_binary::read_binary<bool>(ifs);
  detid_max_            = io_binary::read_binary<decltype(detid_max_)>(ifs);
  // GroupManager
  group_mgr_.load(ifs);
  // UqidManager
  uqid_mgr_.load(ifs);
  // UqigAvailIndexer
  avail_indexer_.load(ifs);
  if (ifs.fail()) THROW_ERROR("DetectorIndexContainer::load failed.");
}
