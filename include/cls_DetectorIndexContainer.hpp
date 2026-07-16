/// @file cls_DetectorIndexContainer.hpp
/// @brief Container for managing detector index information
/// @details
/// This file defines the DetectorIndexContainer class, which integrates three manager classes
/// (GroupManager, UqidManager, and UqigAvailIndexer) to provide unified access to detector indexing.
///
/// **Workflow:**
/// 1. GroupManager is built by DetectorPanelArray::build_index_container
/// 2. UqidManager is built by DetectorPanelArray::build_vec_panel
/// 3. UqigAvailIndexer is built by DetectorPanelArray::build_index_container
///
/// **Thread-safety:**
/// - Read operations are thread-safe after construction
/// - Modification operations (addUqidInfo, addGroupInfo, etc.) are NOT thread-safe
/// - build_uqigAvailIndexer_disabled/enabled use OpenMP internally for parallel construction
///
/// **Key terminology:**
/// - Detid: Detector ID
/// - Uqid: Unique ID for each (detector, ix, iy) combination
/// - Uqig: Unique ID for each (detector, igroup) combination
/// - UqigAvail: Uqig with availability flag
/// - Igroup: Group ID within a detector
#pragma once
#include "cls_UqigAvailIndexer.hpp"
#include "cls_UqidManager.hpp"
#include "cls_GroupManager.hpp"
#include <Eigen/Dense>
#include <set>
#include "st_NmuonVectors.hpp"  // for NmuonVectors definition

#include <filesystem>
namespace fs = std::filesystem;

using namespace index_type_definitions;

/// @class DetectorIndexContainer
/// @brief Unified container for managing detector index information
/// @ingroup detectorClasses
///
/// @details
/// This class integrates three manager classes to provide a unified interface for
/// detector indexing and grouping operations:
/// - GroupManager: Manages detector groups (uqig) and provides fast lookup
/// - UqidManager: Manages unique IDs (uqid) for detector pixels
/// - UqigAvailIndexer: Manages dense indices for available uqig_avail
///
/// **Responsibilities:**
/// - Provide unified access to detector index mappings
/// - Support disabling/enabling specific detectors for analysis
/// - Generate submatrices and subvectors excluding disabled detectors
/// - Serialize/deserialize index information for persistence
///
/// **Invariants:**
/// - Build flags (tf_built_*) accurately reflect construction status
/// - detid_max_ >= maximum detid in the system
/// - UqigAvailIndexer provides dense index [0, N) for N available groups
///
/// **Usage example:**
/// @code
/// DetectorIndexContainer dic;
/// // ... (populated by DetectorPanelArray)
/// dic.build_uqigAvailIndexer();  // Build dense indexer
/// auto disabled_matrix = dic.get_disabled_mat_cov_muon(detid_to_disable, mat_cov);
/// @endcode
///
/// @note This class does NOT own the detector geometry; it only manages indices.
class DetectorIndexContainer {
private:

  /// @brief Maximum Detid
  Detid detid_max_ = DetidNotAssigned;

  /// @brief Whether uqid construction is complete
  bool tf_built_uqid = false;

  /// @brief Whether uqig construction is complete
  bool tf_built_uqig = false;

  /// @brief Whether map_uqig_avail_uqig build is complete
  bool tf_built_uqigAvail = false;

  /// @brief Whether the indexer for UqAvail is built
  bool tf_built_indexer = false;

  /// @brief Instance of GroupManager
  /// @details Manages each group (uqig) and provides fast lookup based on detid and uqig_avail.
  GroupManager group_mgr_;

  /// @brief Instance of UqidManager
  /// @details Manages unique IDs (uqid) and their data, and handles registration, lookup, and updates.
  UqidManager uqid_mgr_;

  /// @brief Instance of UqigAvailIndexer
  /// @details Manages indices for available uqig_avail and provides uqig_avail conversions.
  UqigAvailIndexer avail_indexer_;
public:

  //==================================================================
  /// @name constructor_and_operators
  ///@{

  /// @brief Default constructor
  DetectorIndexContainer() = default;

  /// @brief copy constructor
  DetectorIndexContainer(const DetectorIndexContainer &other) = default;

  /// @brief move constructor
  DetectorIndexContainer(DetectorIndexContainer &&other) noexcept = default;

  /// @brief destructor
  ~DetectorIndexContainer() = default;

  /// @brief Initialize all members
  void initialize();

  /// @brief Inequality operator
  bool operator!=(const DetectorIndexContainer& other) const;

  /// @brief Equality operator
  bool operator==(const DetectorIndexContainer& other) const { return !(*this != other); };

  /// @brief copy assignment operator
  DetectorIndexContainer& operator=(const DetectorIndexContainer &other) = delete;

  /// @brief move assignment operator
  DetectorIndexContainer& operator=(DetectorIndexContainer&& other) noexcept = default;

  ///@} ------------------------------------------------------------------

  //==================================================================
  /// @name getter_functions
  ///@{

  /// @brief Get the maximum Detid
  Detid get_detid_max() const { return detid_max_; }

  /// @brief Whether uqid construction is complete
  bool is_built_uqid() const { return tf_built_uqid; }

  /// @brief Whether uqig construction is complete
  bool is_built_uqig() const { return tf_built_uqig; }

  /// @brief Whether uqigAvail construction is complete
  bool is_built_uqigAvail() const { return tf_built_uqigAvail; }

  /// @brief Whether indexer construction is complete
  bool is_built_indexer() const { return tf_built_indexer; }

  /// @brief get the minimum uqid
  Uqid get_uqid_min() const { return uqid_mgr_.get_uqid_min(); };

  /// @brief get the maximum uqid
  Uqid get_uqid_max() const { return uqid_mgr_.get_uqid_max(); };

  /// @brief get the minimum uqig
  Uqig get_uqig_min() const { return group_mgr_.get_uqig_min(); };

  /// @brief get the maximum uqig
  Uqig get_uqig_max() const { return group_mgr_.get_uqig_max(); };

  /// @brief get the minimum uqig_avail
  UqigAvail get_uqig_avail_min() const { return group_mgr_.get_uqig_avail_min(); };

  /// @brief get the maximum uqig_avail
  UqigAvail get_uqig_avail_max() const { return group_mgr_.get_uqig_avail_max(); };

  /// @brief get the mutable reference to GroupManager
  GroupManager& callGrpMgr() { return group_mgr_; }

  /// @brief get the mutable reference to UqidManager
  UqidManager& callUqidMgr() { return uqid_mgr_; }

  /// @brief get the mutable reference to UqigAvailIndexer
  UqigAvailIndexer& callAvailIndexer() { return avail_indexer_; }

  /// @brief get the const reference to GroupManager
  const GroupManager& getGrpMgr() const { return group_mgr_; }

  /// @brief get the const reference to UqidManager
  const UqidManager& getUqidMgr() const { return uqid_mgr_; }

  /// @brief get the const reference to UqigAvailIndexer
  const UqigAvailIndexer& getAvailIndexer() const { return avail_indexer_; }

  /// @brief get a copy of GroupManager
  GroupManager getGrpMgr_copy() const { return group_mgr_; }
  
  /// @brief get a copy of UqidManager
  UqidManager getUqidMgr_copy() const { return uqid_mgr_; }

  /// @brief get a copy of UqigAvailIndexer
  UqigAvailIndexer getAvailIndexer_copy() const { return avail_indexer_; }

  /// @brief Get UqidInfo for the minimum uqid
  /// @return UqidInfo corresponding to the minimum uqid
  UqidInfo get_uqid_info_min() const {
    return uqid_mgr_.getInfo(uqid_mgr_.get_uqid_min());
  };

  /// @brief Check whether the uqid range is properly set
  /// @return true if the range is set, false otherwise
  bool is_uqid_range_valid() const {
    return ( !uqid_mgr_.get_vecUqid_all().empty() ) 
        && ( uqid_mgr_.get_uqid_max() >= uqid_mgr_.get_uqid_min() );
  };
  
  /// @brief get the vector of GroupInfo from detid
  std::vector<GroupInfo> get_vecGroupInfo_by_detid(const Detid detid) const {
    return group_mgr_.get_vecGroupInfo_by_detid(detid);
  };

  /// @brief Get the set of Ixiy from detid and igroup
  /// @param detid 
  /// @param igroup 
  /// @return set of Ixiy
  std::set<Ixiy> get_set_ixiy_by_detid_igroup( const Detid detid, const Igroup igroup) const;

  /// @brief Get DetIxiy from uqid
  /// @param[in] uqid_in uqid to search for
  /// @return DetIxiy structure containing (detid, ix, iy)
  /// @throws std::runtime_error If uqid_in is out of range (via THROW_ERROR3)
  /// @note Performs range check against uqid_min and uqid_max before lookup
  DetIxiy get_detidixiy( const Uqid uqid_in ) const;

  /// @brief Get the list of (ix, iy) from (detid, igroup)
  std::vector<Ixiy> get_vec_ixiy_by_detid_igroup(
    const Detid detid, const Igroup igroup) const;

  /// @brief Get the list of uqid from (detid, igroup)
  /// @param detid detector ID
  /// @param igroup group ID
  /// @return set of uqid
  std::set<Uqid> get_set_uqid_by_detid_igroup(
    const Detid detid, const Igroup igroup) const;

  /// @brief Get igroup from (detid, ix, iy)
  /// @param detid detector ID
  /// @param ix x index
  /// @param iy y index
  /// @return corresponding igroup
  /// @throws std::out_of_range if not registered
  Igroup get_igroup(const Detid detid, const int ix, const int iy) const;

  /// @brief Get igroup from (detid, ixiy)
  Igroup get_igroup(const Detid detid, const Ixiy& ixiy) const {
    return get_igroup(detid, ixiy[0], ixiy[1]);
  };

  /// @brief Get igroup from (detid, ix, iy)
  Igroup get_igroup(const DetIxiy& detixiy) const {
    return get_igroup(detixiy[0], detixiy[1], detixiy[2]);
  };

  /// @brief Return the number of (ix, iy) that belong to the specified igroup
  int get_n_uqid(const Detid detid_in, const Igroup igroup) const{
    return uqid_mgr_.get_n_uqid_by_detid_igroup(detid_in,igroup);
  };

  /// @brief Return the number of registered igroup
  int get_n_group(const Detid detid) const { return group_mgr_.get_n_group(detid); };

  /// @brief Return the number of registered uqig_avail
  int get_n_uqig_avail() const { return group_mgr_.get_n_uqig_avail(); };

  /// @brief Return the list of registered unique igroup
  std::vector<Igroup> get_vec_igroup(const Detid detid) const {
    return group_mgr_.get_vec_igroup(detid); 
  };

  /// @brief Return the set of registered unique igroup
  std::set<Igroup> get_set_igroup(const Detid detid) const {
    return group_mgr_.get_set_igroup(detid); };

  /// @brief Wrapper that returns the list of registered unique igroup
  Igroup get_igroup_max(const Detid detid) const{
    return group_mgr_.get_igroup_max(detid);
  };

  /// @brief get the set of detid
  std::set<Detid> get_set_detid() const;

  /// @brief Get the list of enabled row indices excluding uqig_avail that belong to the specified detid
  /// @param[in] detid_disabled Detector ID to exclude
  /// @return Vector of dense indices (0-based) for all enabled groups excluding those from detid_disabled
  /// @note Warns if uqig_avail from detid_disabled is not found in indexer
  /// @note Used internally by get_disabled_mat_dNdD, get_disabled_vec_nmuons, get_disabled_mat_cov_muon
  std::vector<Index> get_enabled_indices_excluding_detid(const Detid detid_disabled) const;

  /// @brief Add UqidInfo
  void addUqidInfo( const UqidInfo& uqid_info_in){ uqid_mgr_.insert(uqid_info_in); };

  /// @brief Get const ref of UqidInfo for the specified detid, ix, iy
  const GroupInfo& getGroupInfo(const Detid detid, const int ix, const int iy) const;

  /// @brief Get const ref of GroupInfo corresponding to Index
  const GroupInfo& getGroupInfo_by_Index(const Index index) const;

  /// @brief Get mutable ref of GroupInfo for the specified detid, ix, iy
  GroupInfo& callGroupInfo(const Detid detid, const int ix, const int iy);

  /// @brief get the size of Indexer
  /// @return the size of UqigAvailIndexer
  size_t get_size_indexer() const noexcept { return avail_indexer_.size(); }

  /// @brief Return a matrix with rows (uqig_avail) belonging to the specified detid disabled
  /// @param detid detector ID to disable
  /// @param mat_dNdD original matrix (rows = number of uqig_avail)
  /// @return new matrix with rows belonging to detid removed
  Eigen::MatrixXf
    get_disabled_mat_dNdD(const Detid detid, const Eigen::MatrixXf& mat_dNdD) const;

  /// @brief Return nmuon vectors with uqig_avail belonging to the specified detid disabled
  /// @param detid detector ID to disable
  /// @param vecxf_nmuon_obs original observed vector (size = avail_indexer.size())
  /// @param vecxf_nmuon_prior original prior vector
  /// @return NmuonVectors (struct containing vec_obs and vec_prior)
  NmuonVectors
    get_disabled_vec_nmuons(const Detid detid
    , const Eigen::VectorXf& vecxf_nmuon_obs
    , const Eigen::VectorXf& vecxf_nmuon_prior) const;


  /// @brief Return a covariance matrix (copy) with UqigAvail corresponding to the specified detid excluded
  ///
  /// This function constructs a square matrix (sub-block of the covariance matrix)
  /// by removing rows and columns that correspond to all UqigAvail belonging to
  /// the specified detector ID (detid).
  ///
  /// - The target covariance matrix `mat_cov_muon` must be based on the dense index of UqigAvailIndexer.
  /// - The uqig_avail to remove are automatically obtained from GroupManager and removed based on the indexer.
  /// - The return value is a **newly created copy matrix**, and the original matrix is not modified.
  ///
  /// @param detid detector ID to disable (UqigAvail belonging to this ID are removed from rows/columns)
  /// @param mat_cov_muon original covariance matrix (size must be `avail_indexer_.size()` x `avail_indexer_.size()`)
  /// @return Eigen::MatrixXf square submatrix (copy) excluding UqigAvail belonging to the target detid
  /// @note Uses OpenMP
  Eigen::MatrixXf get_disabled_mat_cov_muon(
      const Detid detid, const Eigen::MatrixXf& mat_cov_muon) const;

  ///@} ------------------------------------------------------------------


  //==================================================================
  /// @name builder_functions
  ///@{

  /// @brief Build UqigAvailIndexer from all available uqig_avail
  /// @throws std::runtime_error If uqid, uqig, or uqigAvail is not built (via check_built_*)
  /// @note Warns if indexer is already built but rebuilds anyway
  /// @note Sets tf_built_indexer to true upon completion
  void build_uqigAvailIndexer();

  /// @brief Build UqigAvailIndexer with specified Detid disabled
  /// @param[in] detid_disabled Detector ID to disable (exclude from indexer)
  /// @throws std::runtime_error If detid_disabled is out of range or prerequisites not built
  /// @note Uses OpenMP for parallel construction
  /// @note Thread-safety: This method is NOT thread-safe during execution
  void build_uqigAvailIndexer_disabled(const Detid detid_disabled);

  /// @brief Build UqigAvailIndexer with a list of Detid disabled
  /// @param[in] vec_detid_disabled List of detector IDs to disable (exclude from indexer)
  /// @throws std::runtime_error If any detid in the list is out of range or no valid uqig_avail remains
  /// @note Duplicates in vec_detid_disabled are automatically eliminated
  /// @note Uses OpenMP for parallel construction
  /// @note Uses C++20 std::exclusive_scan for offset computation
  void build_uqigAvailIndexer_disabled(const std::vector<Detid>& vec_detid_disabled);

  /// @brief Build UqigAvailIndexer using only enabled detid
  /// @param[in] vec_detid_enabled List of enabled detector IDs (only these are included in indexer)
  /// @throws std::runtime_error If any detid in the list is out of range or no valid uqig_avail found
  /// @note Duplicates in vec_detid_enabled are automatically eliminated
  /// @note Uses OpenMP for parallel construction
  /// @note Uses C++20 std::exclusive_scan for offset computation
  void build_uqigAvailIndexer_enabled(const std::vector<Detid>& vec_detid_enabled);

  ///@} ------------------------------------------------------------------


  //==================================================================
  /// @name setter_functions
  ///@{

  /// @brief Set dic_ from another DetectorIndexContainer
  void set(const DetectorIndexContainer &dic_in);

  /// @brief Set detid_max_ of dic_
  void set_detid_max(const Detid detid_max) { detid_max_ = detid_max; }

  /// @brief Set the build status of uqid
  void set_built_uqid(const bool tf) { tf_built_uqid = tf; }

  /// @brief Set the build status of uqig
  void set_built_uqig(const bool tf) { tf_built_uqig = tf; }

  /// @brief Set the build status of uqigAvail
  void set_built_uqigAvail(const bool tf) { tf_built_uqigAvail = tf; }

  /// @brief Set the build status of indexer
  void set_built_indexer(const bool tf) { tf_built_indexer = tf; }

  /// @brief Remove the specified (detector, group) entirely
  void remove_group(const Detid detid, Igroup igroup);

  /// @brief Remove UqidInfo corresponding to the specified uqid
  /// @param uqid_in uqid to remove
  /// @throws std::out_of_range if uqid is out of range or not registered
  void remove_uqid(const Uqid uqid_in) { uqid_mgr_.remove(uqid_in);};

  /// @brief Add GroupInfo
  void addGroupInfo(const GroupInfo& group_info_in) { group_mgr_.insert_to_Bimaps(group_info_in); };

  ///@} ------------------------------------------------------------------

  //==================================================================
  /// @name checker_functions
  ///@{

  /// @brief throw error if uqid is not built
  void check_built_uqid() const;

  /// @brief throw error if uqig is not built
  void check_built_uqig() const;

  /// @brief throw error if uqigAvail is not built
  void check_built_uqigAvail() const;

  /// @brief throw error if indexer is not built
  void check_built_indexer() const;

  ///@} ------------------------------------------------------------------



  /// @brief Output grouping information to a FILE pointer
  /// @param fout output file pointer
  void out_grouping_info(FILE* fout) const;

  /// @brief Output grouping information to a specified file
  /// @param path output file path
  void out_grouping_info(const std::filesystem::path& path) const;

  /// @brief ASCII output of Detid, Ix, Iy, Inthis, Igroup, Uqid, Uqig, is_avail, UqigAvail for each Uqid
  void out_all_index_info(FILE* fout) const;

  /// @brief File-path version
  void out_all_index_info(const std::filesystem::path& path) const;

  /// @brief Output all index information to the specified file
  /// @param file_path 
  /// @param width 
  void out_all_index_info2(const std::filesystem::path& file_path, const int width=9) const;

  /// @brief File-path version
  void out_all_index_info_csv(const std::filesystem::path& file_path) const;

  /// @brief Save binary representation to output stream
  /// @param[in,out] ofs Output file stream
  /// @throws std::runtime_error If save operation fails (via THROW_ERROR)
  /// @note Saves GroupManager, UqidManager, and UqigAvailIndexer in sequence
  void save(std::ofstream& ofs) const;

  /// @brief Load binary representation from input stream
  /// @param[in,out] ifs Input file stream
  /// @throws std::runtime_error If load operation fails (via THROW_ERROR)
  /// @note Loads GroupManager, UqidManager, and UqigAvailIndexer in sequence
  void load(std::ifstream& ifs);  
};
