/// @file cls_DetectorPanelArray.hpp
/// @brief A container for managing multiple \ref DetectorPanel instances
/// @details
/// This file defines DetectorPanelArray, which manages collections of \ref DetectorPanel
/// instances for multi-panel detector configurations. It provides:
/// - Unified indexing across all detector panels via DetectorIndexContainer
/// - Factory methods for cuboid and voxel-based detector geometries
/// - Signal, noise, and path length calculations
/// - Automatic angular binning and grouping
/// - Binary I/O for detector configurations
///
/// ## Workflow
/// 1. Build panels via build_vec_panel() or factory methods (create/create_sprs)
/// 2. Build unified index container via build_index_container()
/// 3. Calculate path lengths (PL) and density lengths (DL)
/// 4. Set signal and noise via mp_calc_set_peneflux_signal_from_DL() and mp_set_noise_all()
/// 5. Perform grouping: mp_assign_1st_igroup_all() then mp_auto_grouping_by_signal_noise_group_alldet()
/// 6. Calculate group statistics via mp_calc_vec_signal_noise_group_all()
///
/// ## Coordinate System
/// Inherits coordinate system from individual DetectorPanel instances
///
/// ## Thread Safety
/// Many methods use OpenMP parallelization (prefixed with mp_). Not thread-safe for concurrent modification.
///
/// @note Units: meters for positions, radians for angles (inherited from DetectorPanel)
#pragma once

#include <map>
#include <fstream>
#include <iostream>
#include <sstream> // istringstream
#include <string>
#include <cstdio>
#include <cmath>
#include <functional>  //for sorting
#include <algorithm>//for sorting
#include <vector>
#include <tuple>
#include <filesystem> // for std::filesystem::path

// for std::map<>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <bitset>

#include <Eigen/Dense>
#include "cls_DetectorElement.hpp"
#include "cls_DetectorPanel.hpp"
#include "cls_DetectorPanelParameters.hpp"
#include "cls_FluxTable.hpp"
#include "cls_PathLengthParameters.hpp" // for PathLengthParameters
#include "cls_Grid3dVoxelParameters.hpp"
#include "ns_io_binary.hpp"
#include "cls_NoiseParameters.hpp"
#include "cls_DetectorIndexContainer.hpp"
#include "st_shapes.hpp" // for VerticalEllipticCylinderOpen
#include "ns_iodir.hpp"

namespace fs = std::filesystem;

// forward declaration
class DetectorPanelArray;

// 1) Define flag list
/// @brief Macro that defines the progress flag list
#define FLG_LIST(FLG_X)                                                  \
  FLG_X(tf_built_panels)    /**< Whether detector panel construction is complete */       \
  FLG_X(tf_calc_PL_DL)      /**< Whether PL/DL calculation is complete */                 \
  FLG_X(tf_calc_signal)     /**< Whether signal calculation is complete */                \
  FLG_X(tf_calc_noise)      /**< Whether noise calculation is complete */                 \
  FLG_X(tf_auto_grouping)   /**< Whether auto grouping is complete */                     \
  FLG_X(tf_calc_signal_group)  /**< Whether vec_signal_group calculation is complete */  \
  FLG_X(tf_calc_noise_group)  /**< Whether vec_noise_poi_group calculation is complete */    \
  FLG_X(tf_calc_proj_dens_group)  /**< Whether vec_proj_dens_group calculation is complete */  \
  FLG_X(tf_calc_volume_group)  /**< Whether vec_volume_group calculation is complete */  \
  FLG_X(tf_set_eff_group)   /**< Whether vec_eff_group setup is complete */

// 2) Enum from the above
/// @brief Enum class for progress flags
enum class FlgProg {
#define DECL_ENUM(name) name,
  FLG_LIST(DECL_ENUM)
#undef DECL_ENUM
  COUNT
};

// 3) Name array (for display)
/// @brief Array of flag names corresponding to FlgProg enum
static constexpr std::array<const char*,
  static_cast<size_t>(FlgProg::COUNT)> FlgProgNames = {
#define DECL_NAME(name) #name,
  FLG_LIST(DECL_NAME)
#undef DECL_NAME
};

//##################################################################################
//##################################################################################
/// @class DetectorPanelArray
/// @brief A collection of multiple detector panels with unified indexing across all angular elements
/// @details
/// DetectorPanelArray is a final (non-inheritable) container class that manages:
/// - A vector of \ref DetectorPanel instances (vec_panel)
/// - Unified indexing via DetectorIndexContainer (dic_) for all DetectorElement instances
/// - Progress tracking flags for workflow stages
/// - Signal/noise/path-length data aggregation
///
/// ## Responsibilities
/// - Construct and manage multiple detector panels
/// - Provide unified element access: by (detid, ix, iy) or by global unique index (uqid)
/// - Calculate path lengths (PL) and density lengths (DL) for cuboid or voxel geometries
/// - Compute signal and noise for each element
/// - Perform automatic angular binning and grouping
/// - Support volume-weighted density calculations
///
/// ## Invariants
/// - n_all_element equals the sum of all panel element counts
/// - dic_ is built before group-level operations (enforced by check_built_uqid())
/// - Grouping must complete before efficiency calculations (enforced by check_auto_grouping())
///
/// ## Typical Usage
/// @code
/// // 1. Build from parameter list
/// DetectorPanelArray arrdet(prm_list);
/// arrdet.build_vec_panel();
/// arrdet.build_index_container();
///
/// // 2. Calculate geometry and signals
/// pathcalc::g2pil::mp_add_PLDL(arrdet, g2pil);
/// arrdet.mp_calc_set_peneflux_signal_from_DL(ft, true);
/// arrdet.mp_set_noise_all(noise_prm, DL_thres);
///
/// // 3. Perform grouping
/// arrdet.mp_copy_signal_noise_to_g2bg_all(prm_bingroup);
/// arrdet.mp_assign_1st_igroup_all(prm_bingroup);
/// arrdet.mp_auto_grouping_by_signal_noise_group_alldet(prm_bingroup);
/// @endcode
///
/// @ingroup detectorClasses
//##################################################################################
//##################################################################################
class DetectorPanelArray final {
  public:
    /// @brief forward declaration of PillarBuildParams
    struct PillarBuildParams;

    /// @brief forward declaration of VoxelBuildParams
    struct VoxelBuildParams;

  private:
    /// @brief name of this instance
    std::string name = "DetectorPanelArray";

    /// @brief num of all DetectorElement
    int n_all_element = 0;

    /// @brief vector array of parameter file path \n
    /// Stores parameter file paths for each detector
    std::vector<fs::path> vec_parameter_file_path = {};

    /// @brief If true, override angular binning for all detectors
    bool tf_override_angle_bin = false;

    /// @brief Override values for angular binning (used when tf_override_angle_bin is true)
    int nbinx_override = 0;
    double txmin_override = 0.0;
    double txmax_override = 0.0;
    int nbiny_override = 0;
    double tymin_override = 0.0;
    double tymax_override = 0.0;

    /// @brief vector array of DetectorPanel
    std::vector<DetectorPanel> vec_panel = {};

    /// @brief unique index map container for all DetectorElement(det_id,ix,iy)) and group index
    // id_container::DEID uqid_maps = id_container::DEID();

    /// @brief DetectorIndexContainer
    DetectorIndexContainer dic_ = DetectorIndexContainer();

    /// @brief Bitset that manages flags defined in FLAG_LIST
    std::bitset<static_cast<size_t>(FlgProg::COUNT)> flags_;

    /// @brief set of unique voxel indices
    std::set<Grid3d::Uqiv> set_uqiv = {};

  public:
    //============================================================================
    /// @name constructor & destructor
    //============================================================================
    ///@{

    /// @brief default constructor
    DetectorPanelArray() = default;

    /// @brief copy constructor with initialization list
    DetectorPanelArray(const DetectorPanelArray &org);

    /// @brief move constructor, some members are not moved.
    DetectorPanelArray(DetectorPanelArray&& other) noexcept;
    
    /// @brief destructor
    ~DetectorPanelArray() = default;

    /// @brief Read DetectorPanel::ParameterLists and build vec_panel based on it.
    DetectorPanelArray(const DetectorPanel::ParameterLists &prm_g2det_list)
      : DetectorPanelArray() {
      vec_parameter_file_path = prm_g2det_list.vec_parameter_file_path;
      tf_override_angle_bin = prm_g2det_list.tf_override_angle_bin;
      nbinx_override = prm_g2det_list.nbinx_override;
      txmin_override = prm_g2det_list.txmin_override;
      txmax_override = prm_g2det_list.txmax_override;
      nbiny_override = prm_g2det_list.nbiny_override;
      tymin_override = prm_g2det_list.tymin_override;
      tymax_override = prm_g2det_list.tymax_override;
      build_vec_panel();
    };

    /// @brief Constructor that loads a binary file saved by save
    /// @param[in] path_in Path to the binary file
    /// @note private member function
    DetectorPanelArray(const fs::path &path_in)
      : DetectorPanelArray() {
      spdlog::info(" load DetectorPanelArray object from {}", path_in.string());
      std::this_thread::sleep_for(std::chrono::seconds(5));
      load(path_in);
    };

    ///@} ---------------------------------------------------------------------

    //============================================================================
    /// @name operators
    ///@{

    /// @brief Assignment operator
    DetectorPanelArray& operator=(const DetectorPanelArray &other) = delete;

    /// @brief Inequality operator
    /// @note Name is not compared.
    bool operator!=(const DetectorPanelArray &other) const;

    /// @brief Equality operator
    /// @details Name is not compared. Defined via the inequality operator.
    bool operator==(const DetectorPanelArray &other) const {
      return !(*this != other);
    };

    /// @brief move assignment operator
    DetectorPanelArray& operator=(DetectorPanelArray&& other) noexcept;

    ///@} ---------------------------------------------------------------------

    //============================================================================
    /// @name factory_functions
    //============================================================================
    ///@{
    
    /// @brief Factory function that selects and calls DetectorPanelArray constructors
    /// \image html images/chart_arrdet_g2pil_create.dio.png "flow chart"
    /// @param prm_arrdet_cub Struct holding parameters and data instance references
    /// @return DetectorPanelArray
    static DetectorPanelArray create( 
        const DetectorPanelArray::PillarBuildParams &prm_arrdet_cub );

    /// @brief Factory function that selects and calls DetectorPanelArray constructors
    /// \image html images/chart_arrdet_g3vox_input_vec_spmat_PL_create.dio.png "flow chart"
    /// @param prm_arrdet_vox Struct holding parameters and data instance references
    /// @return tuple of DetectorPanelArray & std::vector<SpMatf>
    static std::tuple< DetectorPanelArray, std::vector<SpMatf> >
      create( const DetectorPanelArray::VoxelBuildParams &prm_arrdet_vox );

    /// @brief Factory function that selects and calls DetectorPanelArray constructors, fully sparse version
    /// \image html images/chart_arrdet_g3vox_input_vec_spmat_PL_create.dio.png "flow chart"
    /// @param prm_arrdet_vox Struct holding parameters and data instance references
    /// @param opt_dic optional DetectorIndexContainer; can assign if already available.
    /// @return tuple of (DetectorPanelArray, vec_spmat_PL, vecxf_non_rec_vox_PL).
    ///   - vec_spmat_PL: per-detector sparse PL matrix, columns filtered by n_hit_det.
    ///   - vecxf_non_rec_vox_PL: flat per-element PL contribution from voxels that
    ///     were dropped by the n_hit_det filter (i.e. PL columns NOT in
    ///     vec_spmat_PL). Required for prior reconstruction in the compute_prior stage.
    ///     Empty vector when arrdet_g3vox was loaded from cache (in that case
    ///     the value must come from ShellPL::load instead).
    static std::tuple< DetectorPanelArray, std::vector<SpMatf>, Eigen::VectorXf >
      create_sprs( const DetectorPanelArray::VoxelBuildParams &prm_arrdet_vox
        , const std::optional<DetectorPanelArray>& opt_arrdet);

    ///@} ------------------------------------------------------------------------

    //============================================================================
    /// @name flag_functions
    /// @brief get the status flags
    ///@{

    /// @brief get the status of a flag
    void set(FlgProg f, const bool v = true) { flags_.set(static_cast<size_t>(f), v); };

    /// @brief get the status of a flag
    bool get(FlgProg f) const { return flags_.test(static_cast<size_t>(f)); };

    /// @brief Generate per-flag accessors: is_<name>() / set_<name>(bool).
    /// @details Expands FLG_LIST to avoid repetitive boilerplate when adding flags.
    #define DECL_ACCESSOR(name)                                          \
      bool is_##name() const {                                          \
        return flags_.test(static_cast<size_t>(FlgProg::name));         \
      }                                                                 \
      void set_##name(bool v) {                                         \
        flags_.set(static_cast<size_t>(FlgProg::name), v);              \
      }
      FLG_LIST(DECL_ACCESSOR)
    #undef DECL_ACCESSOR


    /// @brief display the status of this instance
    void display_status(
      const spdlog::level::level_enum spdlog_level=spdlog::level::debug) const;

    /// @brief display the status of this instance
    void display_status(FILE *fout) const;

    /// @brief Check if all detector positions are within the grid range
    /// @param xmin Minimum x coordinate of the grid [m]
    /// @param xmax Maximum x coordinate of the grid [m]
    /// @param ymin Minimum y coordinate of the grid [m]
    /// @param ymax Maximum y coordinate of the grid [m]
    /// @param margin_factor Factor to expand the grid range for checking (default: 2.0)
    /// @throws std::runtime_error If any detector is outside the expanded grid range
    /// @note This should be called before path length calculations to detect
    ///       coordinate system mismatches early
    void check_detector_positions_in_grid_xy(
      double xmin, double xmax, double ymin, double ymax,
      double margin_factor = 2.0) const;

    ///@} ------------------------------------------------------------------------


  public:
    //============================================================================
    /// @name build_functions
    //============================================================================
    ///@{

    /// @brief build vec_panel
    /// @param[in] uqid_start starting Uqid for the first DetectorElement in the first DetectorPanel
    /// @param[in] detid_start starting Detid for the first DetectorPanel
    /// @returns n_all_element_sum : total number of DetectorElement in all DetectorPanel
    /// @note use openmp
    int build_vec_panel( const Uqid uqid_start=0, const Detid detid_start=0 );

    /// @brief build DetectorIndexContainer this->dic_ from vec_panels
    std::tuple<Uqig,UqigAvail> build_index_container(
      const Uqig uqig_min=0, const UqigAvail uqig_avail_min=0);

    ///@} ------------------------------------------------------------------

    //============================================================================
    /// @name getter_functions
    /// @brief you can't change the member variables
    ///@{

    /// @brief get the name of this instance
    std::string get_name() const { return name; };

    // getter
    /// @brief get the number of DetectorPanel
    int get_n_det() const { return vec_panel.size(); };

    /// @brief get the exposure time (sec) of Detid=detid_in
    double get_exposure_time_sec_panel( const Detid detid_in ) const;

    /// @brief get the number of all DetectorElement
    int get_n_all_element() const { return n_all_element; };

    /// @brief get the minimum of unique index
    int get_uqid_min() const { return dic_.get_uqid_min(); };

    /// @brief get the maximum of unique index
    int get_uqid_max() const { return dic_.get_uqid_max(); };

    /// @brief get the immutable reference of DetectorIndexContainer
    const DetectorIndexContainer& get_dic() const { return dic_; }

    /// @brief get the mutable reference of DetectorIndexContainer
    DetectorIndexContainer& call_dic() { return dic_; }

    /// @brief get the reference of immutable DetectorElement from detid, ix, iy
    /// @details the called instance cannot be changed
    /// @param[in] detid_in Detid of the DetectorPanel
    /// @param[in] ix_in x index of the DetectorElement in the DetectorPanel
    /// @param[in] iy_in y index of the DetectorElement in the DetectorPanel
    const DetectorElement& getDetectorElement(
      const int detid_in, const int ix_in, const int iy_in) const {
      return getDetectorPanel(detid_in).getDetectorElement(ix_in,iy_in);
    };

    /// @brief get the reference of immutable DetectorElement from detid, ix, iy
    /// @details the called instance cannot be changed
    /// @param[in] detixiy_in std::array<int,3> of (detid, ix, iy)
    const DetectorElement& getDetectorElement( const DetIxiy & detixiy_in ) const {
      const auto& [detid,ix,iy] = detixiy_in;
      return getDetectorElement(detid,ix,iy);
    };

    /// @brief get the reference of immutable DetectorElement from uqid
    /// @details the called instance cannot be changed
    /// @param[in] uqid_in unique index of the DetectorElement
    const DetectorElement& getDetectorElement( const Uqid uqid_in ) const {
      return getDetectorElement(dic_.get_detidixiy(uqid_in));
    };

    /// @brief call the immutable reference of DetectorPanel from detid
    /// @details the called instance cannot be changed
    /// @param[in] detid_in Detid of the DetectorPanel
    const DetectorPanel& getDetectorPanel( const Detid detid_in ) const {
      return vec_panel.at(detid_in);
    };

    /// @brief return deep copy of DetectorPanel from detid using default copy constructor
    /// @param[in] detid_in Detid of the DetectorPanel
    DetectorPanel getDetectorPanelCopy( const Detid detid_in ) const {
      DetectorPanel panel( vec_panel.at(detid_in) );
      return panel;
    };

    /// @brief Return VectorXf of \ref Grid2dBinGroup::get_vecxf_nmuon for all UqigAvail for all detectors.
    /// @note Uses OpenMP
    Eigen::VectorXf get_vecxf_nmuon_all() const;

    /// @brief Return a vector of log10(signal_group + noise_poi_group) ordered by uqig_avail
    /// @param nmuon_group_thres  threshold below which log_value_under_thres is used
    /// @param log_value_under_thres  log value to fill when below threshold
    /// @throws std::out_of_range if uqid or group is not built
    Eigen::VectorXf get_vecxf_log_nmuon_all(
      const double nmuon_group_thres, const double log_value_under_thres ) const;

    /// @brief Return VectorXf of per-row (signal,noise) counts.
    /// @param tf_signal_poisson whether to apply Poisson error to signal counts
    /// @note noise term is always det floor + Poisson(poi bucket); no noise switch.
    /// @note Mirrors get_vecxf_nmuon_all()'s row-indexer + OpenMP structure.
    /// @note Uses OpenMP
    Eigen::VectorXf get_vecxf_nmuon_poisson_all(
        const bool tf_signal_poisson) const;

    /// @brief Return per-row efficiency-variance vector for the C_N diagonal.
    /// @param tf_independent If true, per-element efficiency uncertainties are summed
    ///        in quadrature within each merge group (independent). If false (default),
    ///        they are treated as common within the group (square of the sum).
    /// @return Vector aligned row-for-row with get_vecxf_nmuon_poisson_all().
    /// @note Per element: sigma_eff_i = (eff_upp_i - eff_low_i)/2, contribution = sigma_eff_i * b_i,
    ///       where b_i is the efficiency-free base count (DetectorElement::calc_signal()).
    /// @note Common (default): var_eff_g = (sum_i sigma_eff_i * b_i)^2.
    ///       Independent: var_eff_g = sum_i (sigma_eff_i * b_i)^2.
    /// @note Mirrors get_vecxf_nmuon_poisson_all()'s row indexer so the order matches C_N.
    /// @note Uses OpenMP.
    Eigen::VectorXf get_vecxf_var_eff_all(
        const bool tf_independent = false) const;

    /// @brief Return VectorXf of signal assigned in uqid order for all detectors
    Eigen::VectorXf get_vecxf_signal_all() const;

    /// @brief Return VectorXf by concatenating panel.get_vecxf_nmuoned() for all detectors
    Eigen::VectorXf get_vecxf_nmuoned_all() const;

    /// @brief Return an array of 3 elements (eff_low/cnt/upp_group) ordered by uqig_avail for all detectors
    /// @note Returns Eigen::VectorXf of vec_eff_low_group in uqig_avail order
    /// @throws std::out_of_range If dic_ is not constructed or grouping is incomplete
    std::array<Eigen::VectorXf,3>
      get_vecxf_eff_low_cnt_upp_group_all() const;
    
    /// @brief Return the vector of vecxf_DL for each detector.
    std::vector<Eigen::VectorXf> get_vec_vecxf_DL() const;

    /// @brief Return std::vector of vecxf_DL_prior for all detectors assuming uniform prior density dens_prior.
    std::vector<Eigen::VectorXf> get_vec_vecxf_PL_times_dens(const double dens) const;

    /// @brief get the reference of dic_
    const DetectorIndexContainer& get_dic_ref() const { return dic_; }

    /// @brief Return unique_index (uqid) from det_id, ix, iy
    /// @param det_id Detid of the DetectorPanel
    /// @param ix x index of the DetectorElement in the DetectorPanel
    /// @param iy y index of the DetectorElement in the DetectorPanel
    Uqid get_uqid( const Detid det_id, const int ix, const int iy ) const;

    /// @brief Get det_id from uqid
    /// @param uqid_in unique index of the DetectorElement
    Detid get_detid(const Uqid uqid_in) const;

    /// @brief Convert mat(irow=uqig) --> mat(irow=uqig_avail).
    /// @param mat mat(irow=uqig)
    /// @return mat(irow=uqig_avail)
    /// @note Uses OpenMP
    Eigen::MatrixXf get_matrix_uqig_avail_from_matrix_uqig(
      const Eigen::MatrixXf& mat) const;

    /// @brief get the immutable reference of the set of unique index of voxels
    const std::set<Grid3d::Uqiv>& get_set_uqiv() const { return set_uqiv; };

    /// @brief get the number of unique index of voxels
    int get_n_uqiv() const { return set_uqiv.size(); };

    /// @brief Return a new array with element-wise signal and noise subtracted from another array
    /// @details Calls \ref DetectorPanel::get_subtract() for each DetectorPanel and aggregates results into a new instance
    /// @param other The DetectorPanelArray to subtract from
    /// @return A new DetectorPanelArray containing the element-wise differences
    /// @throws std::out_of_range If dic_ is not constructed
    DetectorPanelArray get_subtract( const DetectorPanelArray &other) const;

    ///@} ------------------------------------------------------------------------

    //============================================================================
    /// @name call_functions
    /// @brief you CAN change the member variables
    ///@{

    /// @brief call the mutable reference of DetectorPanel from detid
    /// @param[in] detid_in Detid of the DetectorPanel
    DetectorPanel& callDetectorPanel( const int detid_in ){
      return vec_panel.at(detid_in);
    };

    /// @brief call the mutable reference of DetectorElement from detid, ix, iy
    /// @param[in] detid_in Detid of the DetectorPanel
    /// @param[in] ix_in x index of the DetectorElement in the DetectorPanel
    /// @param[in] iy_in y index of the DetectorElement in the DetectorPanel
    DetectorElement& callDetectorElement(
      const Detid detid_in, const int ix_in, const int iy_in){
      return callDetectorPanel(detid_in).callDetectorElement(ix_in,iy_in);
    };

    /// @brief call the mutable reference of DetectorElement from unique_index (uqid)
    /// @details Obtain detid, ix, iy via dic_ and return the mutable reference of the corresponding DetectorElement.
    /// @param uqid_in (unique ID)
    /// @return mutable reference of the corresponding DetectorElement
    DetectorElement& callDetectorElement(const Uqid uqid_in){
      const auto& uqid_info = dic_.getUqidMgr().getInfo(uqid_in);
      return callDetectorElement(uqid_info.detid, uqid_info.ixiy[0], uqid_info.ixiy[1]);
    };
    ///@} ------------------------------------------------------------------------

    //============================================================================
    /// @name setter_functions
    ///@{

    /// @brief set the name of this instance
    /// @param[in] name_in Name to set
    void set_name( const std::string name_in ){ name = name_in; };

    /// @brief set the number of all DetectorElement
    /// @param[in] n_all_element_in Number of all DetectorElement to set
    void set_n_all_element( const int n_all_element_in ){
      n_all_element = n_all_element_in;
    };

    /// @brief set the exposure time for all DetectorElement
    /// @param[in] time_sec_in Exposure time in seconds to set
    void set_exposure_time_sec( const double time_sec_in );

    /// @brief Copy panel_in Grid2dBinGroup g2bg_ into vec_panel[det_id]
    /// @param[in] det_id DetectorPanel ID
    /// @param[in] panel_in DetectorPanel to copy from
    void copy_g2bg( const Detid det_id, const DetectorPanel& panel_in );

    /// @brief Copy map_ixiy_igroup and map_igroup_ixiy for all DetectorPanel from another DetectorPanelArray
    /// @param arrdet_src Source DetectorPanelArray
    /// @note Uses OpenMP
    void mp_copy_g2bg_all( const DetectorPanelArray &arrdet_src );

    /// @brief copy bimap_ for all DetectorPanel from other DetectorPanelArray
    /// @param arrdet_src Source DetectorPanelArray
    /// @note Uses OpenMP
    void mp_copy_bimap_all( const DetectorPanelArray &arrdet_src );

    /// @brief Assign noise \n
    /// to each DetectorElement in vec_vec_DetectorElement
    /// @details openmp version
    void mp_set_noise_all( const NoiseParameters &ndist, const double DL_thres);
    
    /// @brief Assign penetrating_muonflux and signal muon counts \n
    /// to each DetectorElement in vec_vec_DetectorElement
    /// @details openmp version
    /// @param[in] tf_apply_eff_cnt If true, central count uses eff_cnt (deterministic, not dice).
    ///            Default false keeps behavior unchanged.
    void mp_calc_set_peneflux_signal_from_DL(
      const FluxTable &ft, const bool tf_apply_eff
    , const bool tf_apply_eff_cnt = false );

    /// @brief Run DetectorPanel::copy_signal_noise_to_g2bg for all DetectorPanel
    /// @details Copy DetectorElement signal/noise to Grid2dBinGroup vec_vec_signal/noise.
    /// @param signal_init initial value of signal
    /// @param noise_init initial value of noise
    /// @param is_avail_init initial value of is_avail
    /// @param PL_thres rocklength threshold
    /// @param DL_thres density length threshold
    /// @param signal_under_thres signal under threshold
    /// @param noise_under_thres noise under threshold
    /// @param is_avail_under_thres is_avail under threshold
    void mp_copy_signal_noise_to_g2bg_all(
        const double signal_init = 0.0 , const double noise_init = 0.0
      , const bool is_avail_init = true
      , const double PL_thres = 0.0
      , const double DL_thres = 0.0
      , const double signal_under_thres = 0.0
      , const double noise_under_thres = 0.0
      , const bool is_avail_under_thres = false );
    
    /// @brief Wrapper for the function above \ref mp_copy_signal_noise_to_g2bg_all(double,double,bool,double,double,double,double,bool).
    /// @param prm_bingrp 
    void mp_copy_signal_noise_to_g2bg_all(const Grid2dBinGroup::Parameters &prm_bingrp ){
        mp_copy_signal_noise_to_g2bg_all(
            prm_bingrp.signal_init, prm_bingrp.noise_init, prm_bingrp.is_avail_init
          , prm_bingrp.PL_thres, prm_bingrp.DL_thres
          , prm_bingrp.signal_under_thres, prm_bingrp.noise_under_thres
          , prm_bingrp.is_avail_under_thres );
    };

    /// @brief Copy DetectorPanel::vec_tf_in_PL for all DetectorPanel
    void mp_copy_vec_tf_in_PL_all( const DetectorPanelArray &arrdet_in );
    
    /// @brief Run DetectorPanel::grouping_by_bin_list for all DetectorPanel
    /// @param prm_bingrp
    /// @note Uses OpenMP
    void mp_grouping_by_bin_list_all(
      const Grid2dBinGroup::Parameters &prm_bingrp );

    /// @brief If tf_run_1st_grouping is true for all Grid2Detector, run assign_1st_igroup; otherwise run assign_naive_igroup
    /// @param tf_run_1st_grouping whether to run initial grouping
    /// @param PL_thres rock length threshold when binning
    /// @param igroup_start initial igroup value
    /// @param nx_div_init number of x-direction bins when binning
    /// @param ny_div_init number of y-direction bins when binning
    /// @param min_signal_noise_sum_thres minimum signal+noise count considered in initial grouping
    /// @note Uses OpenMP
    void mp_assign_1st_igroup_all(
        const bool tf_run_1st_grouping, const int igroup_start
      , const int nx_div_init, const int ny_div_init
      , const double min_signal_noise_sum_thres );
    
    /// @brief Wrapper for the function above \ref mp_assign_1st_igroup_all(bool,int,int,int,double).
    /// @param prm_bingrp Parameters for initial grouping.
    void mp_assign_1st_igroup_all( const Grid2dBinGroup::Parameters &prm_bingrp ){
        mp_assign_1st_igroup_all(
            prm_bingrp.tf_run_1st_grouping
          , prm_bingrp.igroup_start
          , prm_bingrp.nx_div_init
          , prm_bingrp.ny_div_init
          , prm_bingrp.signal_noise_group_trig );
    };

    /// @brief Run DetectorPanel::auto_divide_by_signal_noise_group_all for all DetectorPanel.
    /// @param tf_prefer_split_x : if xlen==ylen, true splits in x direction.
    /// @note Uses OpenMP
    void mp_auto_grouping_by_signal_noise_group_alldet(
        const double signal_noise_group_trig, const int nloop_limit
      , const int ixlen_min, const int iylen_min, const bool tf_prefer_split_x);
    
    /// @brief Wrapper for the function above \ref mp_auto_grouping_by_signal_noise_group_alldet(double,int,int,int,bool).
    void mp_auto_grouping_by_signal_noise_group_alldet(
      const Grid2dBinGroup::Parameters &prm_bingrp){
        mp_auto_grouping_by_signal_noise_group_alldet(
          prm_bingrp.signal_noise_group_trig, prm_bingrp.nloop_limit
        , prm_bingrp.ixlen_min, prm_bingrp.iylen_min
        , prm_bingrp.tf_prefer_split_x);
    };

    /// @brief Check that every group of every DetectorPanel satisfies
    ///        ixlen >= ixlen_min and iylen >= iylen_min.
    /// @throws std::runtime_error If any detector has a group smaller than the minimum.
    /// @note Uses OpenMP. Violations are only recorded inside the parallel loop;
    ///       the exception is thrown after the loop.
    void mp_check_group_ixiylen_min_all(
      const int ixlen_min, const int iylen_min ) const;

    /// @brief Copy signal_group to vec_signal_group for all DetectorPanel
    void mp_calc_vec_signal_noise_group_all( );

    /// @brief Assign different Poisson errors to signal/noise_poi_group for all DetectorPanel
    /// @note Uses OpenMP
    void mp_set_signal_noise_group_poisson_all();

    /// @brief Copy only bimap_ (igroup <-> ixiy map) and allocate group buffers
    /// @details
    ///   After set_dic() copies dic_, this copies only each DetectorPanel's bimap_
    ///   via mp_copy_bimap_all() (Grid2dBinGroup::copy_bimap_from).
    ///   Other g2bg_ data (e.g., signal/noise vectors) are not copied.
    ///   Instead, allocate_vec_value_group_all() reserves the group-related buffers.
    ///   Use this when you only need the mapping and will rebuild group data.
    /// @param src source DetectorPanelArray
    void set_dic_bimap_all_from( const DetectorPanelArray& src );

    /// @brief Run Grid2dBinGroup::allocate_vec_value_group for all detectors.
    /// @note Uses OpenMP
    void allocate_vec_value_group_all();

    /// @brief Set Grid2dBinGroup::Parameters for all DetectorPanel.
    /// @param prm_bingrp Parameters to set
    void set_parameters_all( const Grid2dBinGroup::Parameters &prm_bingrp ){
      for(auto& det : vec_panel) det.set_parameters(prm_bingrp);
    };

    /// @brief Run DetectorPanel::init_len_DL_pene_sig_noi for all DetectorPanel.
    /// @note Uses OpenMP
    void mp_init_PL_DL_pene_sig_noi_all( );

    /// @brief Set set_uqiv
    void set_set_uqiv( const std::set<int> &set_uqiv_in ){
      set_uqiv = set_uqiv_in;
    };

    /// @brief Set DetectorIndexContainer of DetectorPanelArray and set dic_ptr for all DetectorPanels
    void set_dic( const DetectorIndexContainer &dic_src ){
      dic_.set(dic_src);
      set_dic_ptr_all();
    };

    /// @brief Set dic_ptr for all DetectorPanels
    void set_dic_ptr_all( ){ for(auto& det : vec_panel) det.set_dic_ptr(&dic_); };

    /// @brief Run mp_calc_set_proj_dens_grouped for all DetectorPanel
    /// @note Uses OpenMP
    void mp_calc_set_proj_dens_grouped_all(
        const bool tf_signal_poisson
      , const double dens_min, const double dens_max
      , const double dens_step, const double sigma
      , const FluxTable &ft
      , const bool tf_eff = false, const bool tf_eff_independent = false );

    /// @brief Run mp_calc_set_proj_dens_grouped for all DetectorPanel
    /// @param dens_min      minimum of search range
    /// @param dens_max      maximum of search range
    /// @param vec_dens_step vector of step sizes used in each iteration (sorted in descending order)
    /// @param range_factor  scaling factor when updating the search range (e.g., expansion rate relative to initial step size)
    /// @param ft          FluxTable
    /// @note Uses OpenMP
    void mp_calc_set_proj_dens_grouped_all(
      const bool tf_signal_poisson
    , const double dens_min, const double dens_max
    , const std::vector<double> &vec_dens_step, const double range_factor
    , const double sigma, const FluxTable &ft
    , const bool tf_eff = false, const bool tf_eff_independent = false );

    /// @brief For each detector and each grouped angular bin, calculate projected mean density and store data in Grid2dBinGroup.
    /// @note If vec_dens_step.size()==1, run mp_calc_set_proj_dens_grouped_all(dens_min, dens_max, vec_dens_step.at(0), sigma, ft) \n
    /// If vec_dens_step.size()>2, run mp_calc_set_proj_dens_grouped_all(dens_min, dens_max, vec_dens_step, range_factor, sigma, ft)
    /// @param js JSON object in the section about projected density calculation PROJ_DENS_EVAL_GROUPED
    /// @param ft FluxTable
    /// @note The signal Poisson flag is read from js ("tf_signal_poisson"); the noise term
    ///       is always det floor + Poisson(poi bucket), driven by NOISE_PARAMETERS ratios.
    void mp_calc_set_proj_dens_grouped_all(
      const nlohmann::json &js, const FluxTable &ft
    , const bool tf_eff = false, const bool tf_eff_independent = false );
    
    /// @brief Run calc_set_proj_density for all minimum DetectorElement
    /// @param PL_thres rock length threshold
    /// @param DL_thres density length threshold
    /// @note Uses OpenMP
    void mp_calc_set_proj_density_all( const double PL_thres, const double DL_thres );

    /// @brief get the volume weighted average density from vec_dens_group_center
    /// @param prefix_out Prefix for output files become "prefix_out" + ".tmp"
    /// @note and set volume of each group to DetectorPanel, so it cannot be const
    double calc_volume_weighted_average_density(const fs::path& prefix_out);

    /// @brief Calculate volume-weighted average density.
    /// @details For each group, obtain volume and density, and return the lower/center/upper averages weighted by volume.
    /// @note Not const because it sets volume internally.
    /// @param prefix_out Prefix for output files become "prefix_out" + "_low/cnt/upp.tmp"
    /// @return array of 3 elements (lower, center, upper density)
    std::array<double, 3>
      calc_volume_weighted_average_density_lower_center_upper(
        const fs::path& prefix_out= "vol_wei_avr_dens_detail" );

    /// @brief get the volume sum of cone like shape for det_id
    /// @details used in \ref calc_volume_weighted_average_density
    /// @param det_id Detid of the DetectorPanel
    /// @return volume sum of cone like shape
    double get_volume_sum(const Detid det_id) const;

    /// @brief Set the group efficiency triplets of every group in every panel.
    /// @details For each group, eff_g^X = (1/n_g) * sum_{i in g} eff_i^X for X = low, cnt, upp
    ///          (simple unweighted mean of the element efficiencies, see
    ///          \ref DetectorPanel::calc_eff_mean_group). Deterministic and independent of
    ///          any density model. The former stochastic signal-ratio formula with external
    ///          error constants was replaced (it had zero callers, so the signature change
    ///          has no compatibility impact).
    /// @note should be done after grouping and DetectorPanelArray::mp_allocate_vec_value_group_all
    /// @throws std::runtime_error If auto grouping has not completed.
    void calc_set_eff_group();

    /// @brief set PL and DL for all DetectorElement in all DetectorPanel
    /// @param PL_init Path Length initial value
    /// @param DL_init Density Length initial value
    void init_PLDL_all( const double PL_init=0.0, const double DL_init=0.0);

    ///@} ------------------------------------------------------------------------

    //============================================================================
    /// @name read_write_functions
    ///@{

    /// @brief get map of Detid and opened FILE*
    /// @param prefix file name prefix
    /// @param suffix file name suffix
    /// @return map of det_id and opened FILE*
    /// @throws std::runtime_error if file cannot be opened
    /// @throws std::runtime_error if n_det <= 0
    std::map<Detid, FILE*> 
      get_outfile_map_and_out_header(
        const std::string& prefix, const std::string& suffix="") const;

    /// @brief close file map opened in get_outfile_map_and_out_header
    /// @param file_map map of det_id and opened FILE*
    void close_g2bg_outfile_map(std::map<Detid, FILE*>& file_map) const;

    /// @brief Output signal_group to file for each det_id across all uqig_avail
    void out_g2bg_all() const;

    /// @brief Output vecxf to file for each det_id across all uqig_avail
    void out_g2bg_all( const Eigen::VectorXf &vecxf
                      , const std::string &vecxf_name ) const;

    /// @brief Output (tx, ty, arbitrary data) for all DetectorElement \n
    /// basic function of \ref out_txtyDL, \ref out_txtyPL, \ref out_txtySignal, \ref out_txtyNoise, \ref out_txtySignalPlusNoise, \ref out_txtyDens
    void out_txtyData( const std::string& suffix
      , std::function<double(const DetectorElement&)> data_extractor ) const;
    
    /// @brief Output (tx, ty, Density Length) for all DetectorElement
    void out_txtyDL() const {
      out_txtyData("_DL", [](const DetectorElement& ele) {
        return ele.get_DL();
      });
    };

    /// @brief Output (tx, ty, Path Length) for all DetectorElement
    void out_txtyPL() const {
      out_txtyData("_PL", [](const DetectorElement& ele) {
        return ele.get_PL();
      });
    };

    /// @brief Output (tx, ty, Signal) for all DetectorElement
    void out_txtySignal() const {
      out_txtyData("_signal", [](const DetectorElement& ele) {
        return ele.get_signal();
      });
    };

    /// @brief Output (tx, ty, Noise) for all DetectorElement
    void out_txtyNoise() const {
      out_txtyData("_noise", [](const DetectorElement& ele) {
        return ele.get_noise();
      });
    };

    /// @brief Output (tx, ty, Signal+Noise) for all DetectorElement
    void out_txtySignalPlusNoise() const {
      out_txtyData("_signal_plus_noise", [](const DetectorElement& ele) {
        return ele.get_signal_plus_noise();
      });
    };

    /// @brief Output (tx, ty, Dens) for all DetectorElement
    void out_txtyDens() const {
      out_txtyData("_dens", [](const DetectorElement& ele) {
        return ele.get_proj_density();
      });
    };
    
    /// @brief output file names are name + "_" + panel.get_name() + ".tmp", calls \ref DetectorPanel::out_all
    void out_vec_panel( ) const;

    /// @brief Allocate vec_value_group for all DetectorPanel, calls \ref Grid2dBinGroup::allocate_vec_value_group()
    /// @note Uses OpenMP
    void mp_allocate_vec_value_group_all();

    /// @brief Output (uqig_avail, uqig, det_id, igroup, is_avail) for 
    /// 
    void out_UqidInfo_all( const fs::path& pathout ) const{
      const fs::path global_dir_path = iodir::make_pathout(pathout);
      dic_.getUqidMgr().out_UqidInfo_all(global_dir_path);
    }

    ///@} ------------------------------------------------------------------------

    //============================================================================
    /// @name other_functions
    ///@{
    
    /// @brief For all DetectorElement in *this and arrdet_in, \n
    /// add (or subtract) rock_length and DL.
    /// @remark for debug, not used now?
    /// @return arrdetay = this - arrdet_sub
    void add_sub_path_and_DL(
      const DetectorPanelArray &arrdet_in, const double sign );
    
    /// @brief
    void check_built_uqid() const {
      if( !dic_.is_built_uqid() ){
        LOG_ERROR("DetectorPanelArray::dic_ is not built");
        LOG_ERROR("Please check the dic_ is built, before using this method.");
        THROW_ERROR("DetectorPanelArray::dic_ is not built");
      }
    };

    /// @brief Throw an error if tf_auto_grouping_done is false
    void check_auto_grouping() const{
      if( !get(FlgProg::tf_auto_grouping) ){
        LOG_ERROR("DetectorPanelArray::tf_auto_grouping is false");
        LOG_ERROR("Please check the grouping process, before using this method.");
        THROW_ERROR("DetectorPanelArray::tf_auto_grouping is false");
      }
    };

    /// @brief Throw an error if tf_done_eff_group is false
    void check_set_eff_group() const{
      if( !get(FlgProg::tf_set_eff_group) ){
        LOG_ERROR("DetectorPanelArray::tf_set_eff_group is false");
        LOG_ERROR("Please check the caliculation of eff_group process is done or not, before using this method.");
        THROW_ERROR("DetectorPanelArray::tf_set_eff_group is false");
      }
    };

    /// @brief disp volume sum of all DetectorPanel
    /// @return volume sum of all DetectorPanel
    double disp_volume_sum_all(
      spdlog::level::level_enum spdlog_level = spdlog::level::info ) const;

    ///@} ------------------------------------------------------------------------

    //============================================================================
    /// @name binary_io_functions
    ///@{
    
    // write vec_parameter_file_path to std::ofstream &ofs with binary mode
    void write_vec_parameter_file_path( std::ofstream &ofs ) const;

    // read vec_parameter_file_path from std::ifstream &ifs with binary mode
    void read_vec_parameter_file_path( std::ifstream &ifs );

    // write vec_panel to std::ofstream &ofs with binary mode
    void write_vec_panel( std::ofstream &ofs ) const;

    // read vec_panel from std::ifstream &ifs with binary mode
    void read_vec_panel( std::ifstream &ifs );

    // write all private members to std::ofstream &ofs with binary mode
    void save( std::ofstream &ofs ) const;

    // read all private members from std::ifstream &ifs with binary mode
    void load( std::ifstream &ifs );

    // write all private members to fs::path &pathout with binary mode
    void save( const fs::path& pathout ) const {
      std::ofstream ofs = io_binary::open_ofstream(pathout);
      // First, write architecture information
      io_binary::ArchitectureInfo currentInfo = io_binary::get_current_architecture_info();
      io_binary::write_architecture_info(ofs, currentInfo);
      // Then, write data
      save(ofs); ofs.close();
    };

    // read all private members from fs::path &path_in with binary mode
    void load( const fs::path &path_in ){
      std::ifstream ifs = io_binary::open_ifstream(path_in);
      // 1) Read architecture information stored in the file
      io_binary::ArchitectureInfo fileInfo = io_binary::read_architecture_info(ifs);
      // 2) Check compatibility with current environment (throw exception if mismatch)
      io_binary::check_architecture_compatibility_or_throw(fileInfo);
      // 3) Load data
      load(ifs); ifs.close();
    };

    ///@} ------------------------------------------------------------------------

};

//##################################################################################
/// @struct DetectorPanelArray::PillarBuildParams
/// @brief Struct that aggregates parameters needed to compute path length and signals for Grid2dPillar in DetectorPanelArray
//##################################################################################
struct DetectorPanelArray::PillarBuildParams {
  const bool tf_load_arrdet_g2pil; ///< Whether to load DetectorPanelArray binary file
  const bool tf_save_arrdet_g2pil; ///< Whether to save DetectorPanelArray binary file
  const bool tf_run_1st_grouping; ///< Whether to run initial grouping
  const bool tf_run_auto_grouping; ///< Whether to run auto grouping
  const fs::path& path_arrdet_bin; ///< Path to DetectorPanelArray binary file
  const DetectorPanelArray& arrdet_built; ///< Built arrdet_built template
  const Grid2dBinGroup::Parameters& prm_bingroup; ///< Grid2dBinGroup::Parameters
  const Grid2dPillar& g2pil; ///< g2pil with data loaded
  const FluxTable& ft; ///< FluxRangeDataTable with data loaded
  const bool tf_apply_eff; ///< Whether to apply efficiency to signal
  const NoiseParameters& prm_noise; ///< NoiseParameters parameters

  /// @brief constructor with input values
  /// @note Cannot be defined with C++11 default because of const members.
  PillarBuildParams(
    const bool tf_load_arrdet_g2pil, const bool tf_save_arrdet_g2pil
  , const bool tf_run_1st_grouping, const bool tf_run_auto_grouping
  , const fs::path& path_arrdet_bin_in
  , const DetectorPanelArray& arrdet
  , const Grid2dBinGroup::Parameters& bgp
  , const Grid2dPillar& gc
  , const FluxTable& frdt_in
  , const bool tf_apply_eff_in
  , const NoiseParameters& prm_noise_in
  )
    : tf_load_arrdet_g2pil(tf_load_arrdet_g2pil)
    , tf_save_arrdet_g2pil(tf_save_arrdet_g2pil)
    , tf_run_1st_grouping(tf_run_1st_grouping)
    , tf_run_auto_grouping(tf_run_auto_grouping)
    , path_arrdet_bin(path_arrdet_bin_in)
    , arrdet_built(arrdet)
    , prm_bingroup(bgp)
    , g2pil(gc)
    , ft(frdt_in)
    , tf_apply_eff(tf_apply_eff_in)
    , prm_noise(prm_noise_in)
    {};

  /// @brief Inequality operator
  bool operator!=(const PillarBuildParams& other) const;

  /// @brief Equality operator is defined as the negation of the inequality operator
  bool operator==(const PillarBuildParams& other) const { return !(*this != other); };
};

//##################################################################################
/// @struct DetectorPanelArray::VoxelBuildParams
/// @brief Struct that aggregates parameters needed to compute path length and signals for Grid3dVoxel in DetectorPanelArray
//##################################################################################
struct DetectorPanelArray::VoxelBuildParams {
  const bool tf_load_arrdet_g3vox; ///< Whether to load DetectorPanelArray binary file
  const bool tf_save_arrdet_g3vox; ///< Whether to save DetectorPanelArray binary file
  const bool tf_run_1st_grouping; ///< Whether to run initial grouping
  const bool tf_run_auto_grouping; ///< Whether to run auto grouping
  const fs::path& path_arrdet_bin; ///< Path to DetectorPanelArray binary file&
  const fs::path& path_vec_spmat_PL_bin; ///< Path to std::vector<Eigen::MatrixXf> binary file&
  const DetectorPanelArray& arrdet_built; ///< Built arrdet_built template&
  const Grid2dBinGroup::Parameters& prm_bingroup; ///< Grid2dBinGroup::Parameters&
  Grid3dVoxel& g3vox; ///< g3vox with data loaded&
  const FluxTable& ft; ///< FluxRangeDataTable with data loaded&
  const pathcalc::Parameters& prm_pathcalc; ///< pathcalc::Parameters&
  const Grid3dVoxel::Parameters& prm_g3vox; ///< Grid3dVoxelParameters&
  const Grid2dPillar &g2pil_shell_upper; ///< Reference to upper g2pil_shell (terrain above g3vox)
  const Grid2dPillar &g2pil_shell_lower; ///< Reference to lower g2pil_shell (terrain below g3vox)
  const Grid2dPillar &g2pil_shell_lateral; ///< Reference to lateral g2pil_shell
  const bool has_shell_upper; ///< Whether upper shell has valid cuboids
  const bool has_shell_lower; ///< Whether lower shell has valid cuboids
  const bool has_shell_lateral; ///< Whether lateral shell has valid cuboids
  const int uqiv_start; ///< Start value of unique_index
  const bool tf_apply_eff; ///< Whether to apply efficiency to signal
  const NoiseParameters& prm_noise; ///< NoiseParameters parameters
  const bool tf_apply_eff_cnt; ///< If true, central count uses eff_cnt (deterministic). Gated by tf_eff_cn_diag.

  /// @brief constructor with input values
  /// @note Cannot be defined with C++11 default because of const members.
  VoxelBuildParams(
    const bool tf_load_arrdet_g3vox_in, const bool tf_save_arrdet_g3vox_in
  , const bool tf_run_1st_grouping_in, const bool tf_run_auto_grouping_in
  , const fs::path& path_arrdet_bin_in
  , const fs::path& path_vec_spmat_PL_bin_in
  , const DetectorPanelArray& arrdet
  , const Grid2dBinGroup::Parameters& bgp
  , Grid3dVoxel& g3vox_in
  , const FluxTable& frdt_in
  , const pathcalc::Parameters& prm_pathcalc_in
  , const Grid3dVoxel::Parameters& prm_g3vox_in
  , const Grid2dPillar &g2pil_shell_upper_in
  , const Grid2dPillar &g2pil_shell_lower_in
  , const Grid2dPillar &g2pil_shell_lateral_in
  , const bool has_shell_upper_in
  , const bool has_shell_lower_in
  , const bool has_shell_lateral_in
  , const int uqiv_start_in
  , const bool tf_apply_eff_in
  , const NoiseParameters& prm_noise_in
  , const bool tf_apply_eff_cnt_in = false)
    : tf_load_arrdet_g3vox(tf_load_arrdet_g3vox_in)
    , tf_save_arrdet_g3vox(tf_save_arrdet_g3vox_in)
    , tf_run_1st_grouping(tf_run_1st_grouping_in)
    , tf_run_auto_grouping(tf_run_auto_grouping_in)
    , path_arrdet_bin(path_arrdet_bin_in)
    , path_vec_spmat_PL_bin(path_vec_spmat_PL_bin_in)
    , arrdet_built(arrdet)
    , prm_bingroup(bgp)
    , g3vox(g3vox_in)
    , ft(frdt_in)
    , prm_pathcalc(prm_pathcalc_in)
    , prm_g3vox(prm_g3vox_in)
    , g2pil_shell_upper(g2pil_shell_upper_in)
    , g2pil_shell_lower(g2pil_shell_lower_in)
    , g2pil_shell_lateral(g2pil_shell_lateral_in)
    , has_shell_upper(has_shell_upper_in)
    , has_shell_lower(has_shell_lower_in)
    , has_shell_lateral(has_shell_lateral_in)
    , uqiv_start(uqiv_start_in)
    , tf_apply_eff(tf_apply_eff_in)
    , prm_noise(prm_noise_in)
    , tf_apply_eff_cnt(tf_apply_eff_cnt_in)
    {};

  /// @brief Inequality operator
  bool operator!=(const VoxelBuildParams& other) const;

  /// @brief Equality operator is defined as the negation of the inequality operator
  bool operator==(const VoxelBuildParams& other) const { return !(*this != other); };

};
