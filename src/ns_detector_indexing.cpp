// src/detector_indexing.cpp
#include "ns_detector_indexing.hpp"
#include "ns_myapp.hpp"

// @brief Extract columns specified by map_uqiv_old_new_avail values and create a new matrix
Eigen::MatrixXf detector_indexing::get_unavaible_cols_removed_matrix(
  const Eigen::MatrixXf &mat_org, const std::map<Grid3d::Uqiv,Grid3d::Uqiv> &map_uqiv_old_new_avail) 
{
  const int rows = mat_org.rows();
  const int nmap = map_uqiv_old_new_avail.size();
  
  // Create new matrix
  Eigen::MatrixXf mat_new(rows, nmap);

  // Copy only available columns to new matrix
  for (const auto &[uqiv_old, uqid_new] : map_uqiv_old_new_avail) {
    mat_new.col(uqid_new) = mat_org.col(uqiv_old);
  }

  return mat_new;
}

// @brief Extract columns specified by map_uqiv_old_new_avail values and create a new matrix
// @details Sparse matrix version
SpMatf detector_indexing::get_unavaible_cols_removed_matrix(
  const SpMatf &mat_org, const std::map<Grid3d::Uqiv,Grid3d::Uqiv> &map_uqiv_old_new_avail)
{
  const int rows = mat_org.rows();
  const int nmap = map_uqiv_old_new_avail.size();
  
  // Prepare Triplet list for creating new sparse matrix
  std::vector<Eigen::Triplet<float>> triplets;

  // For each element in map, copy non-zero elements from old column (uqiv_old) to new column (uqid_new)
  for (const auto &[uqiv_old, uqid_new] : map_uqiv_old_new_avail) {
    // Iterate over non-zero elements in column uqiv_old of the original sparse matrix
    for (SpMatf::InnerIterator it(mat_org, uqiv_old); it; ++it) {
      triplets.emplace_back(it.row(), uqid_new, it.value());
    }
  }
  
  // Build new sparse matrix from Triplet list
  SpMatf mat_new(rows, nmap);
  mat_new.setFromTriplets(triplets.begin(), triplets.end());
  mat_new.makeCompressed();
  
  return mat_new;
}

// Vector version of get_unavaible_cols_removed_matrix
std::vector<Eigen::MatrixXf> detector_indexing::get_unavaible_cols_removed_matrix(
  const std::vector<Eigen::MatrixXf> &vec_mat_old, const std::map<Grid3d::Uqiv,Grid3d::Uqiv> &map_uqiv_old_new_avail)
{
  // Get vector size
  const int n_vec_mat = vec_mat_old.size();

  LOG_DEBUG(
    ", before vec_mat_old.size() = {}"
    , n_vec_mat);

  // Create new vector
  std::vector<Eigen::MatrixXf> vec_mat_new;
  vec_mat_new.resize(n_vec_mat);
  LOG_INFO(", vec_mat_new.size() = {}", vec_mat_new.size());

  // Execute get_unavaible_cols_removed_matrix for all elements in vec_mat
  for (int i = 0; i < n_vec_mat; ++i) {
    vec_mat_new.at(i)
      = get_unavaible_cols_removed_matrix(vec_mat_old.at(i), map_uqiv_old_new_avail);
  }
  LOG_DEBUG(
    ", after vec_mat_new.size() = {}"
    , vec_mat_new.size());

  return vec_mat_new;
}


// Vector<SpMatf> version of get_unavaible_cols_removed_matrix
std::vector<SpMatf> detector_indexing::get_unavaible_cols_removed_matrix(
  const std::vector<SpMatf> &vec_mat_old, const std::map<Grid3d::Uqiv,Grid3d::Uqiv> &map_uqiv_old_new_avail)
{
  // Get vector size
  const int n_vec_mat = vec_mat_old.size();

  LOG_DEBUG(
    ", before vec_mat_old.size() = {}"
    , n_vec_mat);

  // Create new vector
  std::vector<SpMatf> vec_mat_new;
  vec_mat_new.resize(n_vec_mat);
  LOG_INFO(", vec_mat_new.size() = {}", vec_mat_new.size());

  // Execute get_unavaible_cols_removed_matrix for all elements in vec_mat
  for (int i = 0; i < n_vec_mat; ++i) {
    vec_mat_new.at(i)
      = get_unavaible_cols_removed_matrix(vec_mat_old.at(i), map_uqiv_old_new_avail);
  }
  LOG_DEBUG(
    ", after vec_mat_new.size() = {}"
    , vec_mat_new.size());

  return vec_mat_new;
}

/// @brief select rows from mat based on indices
/// @param[in] mat Source matrix
/// @param[in] indices Row indices to select
/// @return New matrix containing selected rows
Eigen::MatrixXf detector_indexing::selectRows(
  const Eigen::MatrixXf& mat, const std::vector<int>& keep_indices)
{
  const int rows = mat.rows();
  const int cols = mat.cols();
  LOG_DEBUG("mat.rows()={}, mat.cols()={}", rows, cols);

  Eigen::MatrixXf result(keep_indices.size(), cols);

  #pragma omp parallel for
  for (int i = 0; i < static_cast<int>(keep_indices.size()); ++i) {
    int idx = keep_indices.at(i);

    if (idx < 0 || idx >= rows) {
      // Within OpenMP parallel region, use critical section for safe exception throwing
      #pragma omp critical
      {
        LOG_ERROR("invalid row index={} (mat.rows()={})", idx, rows);
        THROW_ERROR("detector_indexing::selectRows: invalid row index");
      }
    }

    result.row(i) = mat.row(idx);
  }

  return result;
}

/// @brief select square block from mat based on indices
/// @param mat Source matrix
/// @param keep_indices Row and column indices to select
/// @return New matrix containing selected rows and columns 
Eigen::MatrixXf detector_indexing::selectSquareBlock(
  const Eigen::MatrixXf& mat, const std::vector<int>& keep_indices )
{
  const int n = keep_indices.size();
  Eigen::MatrixXf result(n, n);

  #pragma omp parallel for collapse(2)
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      result(i, j) = mat(keep_indices[i], keep_indices[j]);
    }
  }

  return result;
}

/// @brief get the unique, sorted Detid set from sorted_Detid_UqigAvail
std::set<Detid> detector_indexing::get_set_Detid(
  const SortedDetidUqigSet& sorted_Detid_UqigAvail)
{
  std::set<Detid> set_Detid;
  for (const auto& pair : sorted_Detid_UqigAvail) {
    set_Detid.insert(pair.first);
  }
  return set_Detid;
}

// @brief get the get_vec_uqig_avail
// @param[in] detid_in Target detector ID
// @return std::vector<UqigAvail> Vector of UqigAvail corresponding to the target detector ID
std::vector<UqigAvail> detector_indexing::get_vec_uqig_avail(
  const Detid detid_in,
  const SortedDetidUqigSet& sorted_Detid_UqigAvail)
{
  std::vector<UqigAvail> vec_uqig_avail;

  for (const auto& [detid, uqig_avail] : sorted_Detid_UqigAvail) {
    if (detid == detid_in) {
      vec_uqig_avail.push_back(uqig_avail);
    }
    else if (detid > detid_in) {
      // Detid is in ascending order, so skip remaining elements
      break;
    }
  }

  return vec_uqig_avail;
}

/// @brief get the vector of newly assigned UqigAvail_new by disabling by detid_in
/// @param[in] detid_in Target detector ID
/// @return MapUqigAvailDisabled Map of newly assigned indices after disabling specified detector ID
MapUqigAvailDisabled
  detector_indexing::get_map_uqig_disabled_avail(
    const SortedDetidUqigSet& sorted_Detid_UqigAvail
  , const Detid detid_in
  , const UqigAvailDisabled uqig_disabled_avail_start )
{
  MapUqigAvailDisabled map_uqig_disabled_avail;

  UqigAvailDisabled new_uqig_disabled_avail = uqig_disabled_avail_start;
  for (const auto& [detid, uqigAvail] : sorted_Detid_UqigAvail) {
    if (detid == detid_in) continue;
    map_uqig_disabled_avail[new_uqig_disabled_avail] = uqigAvail;
    new_uqig_disabled_avail++;
  }

  return map_uqig_disabled_avail;
}

/// @brief make and get the disabled mat_dNdD by detid
/// @param[in] detid_in Target detector ID
/// @param[in] uqig_disabled_avail_start Start index for assignment
/// @return Matrix obtained by extracting rows from mat_dNdD based on the uqig of specified detector ID
/// @note Must return a copy, not a reference
Eigen::MatrixXf detector_indexing::get_disable_mat_dNdD(
    const SortedDetidUqigSet& sorted_Detid_UqigAvail
  , const Eigen::MatrixXf& mat_dNdD_org
  , const Detid detid_in
  , const UqigAvailDisabled uqig_disabled_avail_start)
{
  // ! debug
  myapp::write_sorted_detid_uqigavail_to_file(sorted_Detid_UqigAvail, "sorted_Detid_UqigAvail.tmp");

  LOG_TRACE_ND("before detector_indexing::get_disable_mat_dNdD");
  const auto map_disabled_avail
     = detector_indexing::get_map_uqig_disabled_avail(
        sorted_Detid_UqigAvail,detid_in, uqig_disabled_avail_start);

  // ! debug
  char buf[256];
  snprintf(buf, sizeof(buf), "map_disabled_avail_detid%02d.tmp", detid_in);
  myapp::out_map(map_disabled_avail, buf);
  LOG_TRACE_ND("done");

  std::vector<UqigAvailDisabled> vec_UqigAvailDisabled;
  vec_UqigAvailDisabled.reserve(map_disabled_avail.size());
  for (const auto& [disabled, _] : map_disabled_avail) {
    vec_UqigAvailDisabled.push_back(disabled);
  }

  return detector_indexing::selectRows(mat_dNdD_org, vec_UqigAvailDisabled);
}

/// @brief make and get the disable_detid_vec_nmuon
/// @param[in] detid_in Target detector ID
/// @param[in] uqig_disabled_avail_start Start value for disabled uqig_avail
/// @return std::array<Eigen::VectorXf,2> 0:vecxf_nmuon_obs, 1:vecxf_nmuon_prior
NmuonVectors detector_indexing::get_disable_vec_nmuons(
    const SortedDetidUqigSet& sorted_Detid_UqigAvail
  , const Eigen::VectorXf &vecxf_nmuon_obs
  , const Eigen::VectorXf &vecxf_nmuon_prior
  , const Detid detid_in
  , const UqigAvailDisabled uqig_disabled_avail_start )
{
  const auto map_disabled_avail
     = detector_indexing::get_map_uqig_disabled_avail(
      sorted_Detid_UqigAvail, detid_in, uqig_disabled_avail_start );

  const int n = map_disabled_avail.size();
  Eigen::VectorXf vec_obs(n), vec_prior(n);

  // Assumes disabled indices are consecutive integers
  for (const auto& [disabled, avail] : map_disabled_avail) {
    vec_obs(static_cast<int>(disabled)) = vecxf_nmuon_obs(avail);
    vec_prior(static_cast<int>(disabled)) = vecxf_nmuon_prior(avail);
  }

  return {
    .vec_obs = vec_obs
  , .vec_prior = vec_prior
  };
}

/// @brief make and get the disabled mat_dNdD by detid
/// @param[in] detid_in Target detector ID
/// @param[in] uqig_disabled_avail_start Start index for assignment
/// @return Eigen::MatrixXf mat_dNdD excluding UqigAvail corresponding to the target detector ID
/// @note Must absolutely return a copy, not a reference
Eigen::MatrixXf detector_indexing::get_disable_mat_cov_muon(
    const SortedDetidUqigSet& sorted_Detid_UqigAvail
  , const Eigen::MatrixXf& mat_cov_muon
  , const Detid detid_in
  , const UqigAvailDisabled uqig_disabled_avail_start )
{
  const auto map_disabled_avail
   = detector_indexing::get_map_uqig_disabled_avail(
    sorted_Detid_UqigAvail, detid_in, uqig_disabled_avail_start );

  std::vector<UqigAvailDisabled> vec_UqigAvailDisabled;
  vec_UqigAvailDisabled.reserve(map_disabled_avail.size());
  for (const auto& [disabled, _] : map_disabled_avail) {
    vec_UqigAvailDisabled.push_back(disabled);
  }

  return detector_indexing::selectSquareBlock(mat_cov_muon, vec_UqigAvailDisabled);
}
