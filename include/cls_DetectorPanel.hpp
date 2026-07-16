/// @file cls_DetectorPanel.hpp
/// @brief Represents a muon detector panel
/// @details DetectorPanel has many \ref DetectorElement as its members, and also has \ref Grid2dBinGroup as a member.
#pragma once

#include <map>
#include <fstream>
#include <iostream>
#include <sstream> // istringstream
#include <memory>
#include <string>
#include <cstdio>
#include <cmath>
#include <functional>  //for sorting
#include <algorithm>//for sorting
#include <vector>

#include <utility> // for std::pair
#include <iterator> // for std::begin(), std::end()
#include <filesystem> // for fs::path
namespace fs = std::filesystem;

#include "cls_Angle.hpp"
#include "ns_tuple_int.hpp"

#include <Eigen/Dense>
#include "ns_myapp.hpp"
#include "cls_Ray.hpp"
#include "cls_AABB.hpp"
#include "cls_DetectorElement.hpp"
#include "cls_Grid2dBinGroup.hpp"
#include "cls_Grid2dBinGroupParameters.hpp"
#include "cls_Grid2dXYZ.hpp"
#include "cls_DetectorIndexContainer.hpp"
#include "cls_FluxTable.hpp"
#include "ns_angle_util.hpp"
#include "st_SignalNoiseStat.hpp"

// forward declaration (full definition in cls_EfficiencyModel.hpp)
class EfficiencyModel;

//###########################################################################
//###########################################################################
/// @class DetectorPanel
/// @brief Class that represents one detector with many angular elements (class DetectorElement). \n
/// (Has Grid2dBinGroup as a member) \n
/// @ingroup detectorClasses
//###########################################################################
//###########################################################################
class DetectorPanel final {
  public:
    /// @brief rotation matrix type
    /// @note LOCAL : After rotation, the reference axes rotate with the object. \n
    /// GLOBAL : After rotation, the reference axes remain fixed.
    // forward declaration
    class Parameters;
    class ParameterLists;

  private:
    /// @brief name of this instance
    std::string name = "none";

    /// @brief Grid2dBinGroup
    Grid2dBinGroup g2bg_ = Grid2dBinGroup();

    /// @brief Pointer of DetectorIndexContainer from DetectorPanelArray::dic_
    DetectorIndexContainer* pdic_ = nullptr;

    /// @brief detid is unique index of many detectors
    Detid detid_ = DetidNotAssigned;

    /// @brief detector center position and center direction
    Ray3d ray3d = Ray3d();

    /// @brief meter. [0]=horizontal, [1]=vertical, 2=[depth]
    Eigen::Vector3d v3_det_length = Eigen::Vector3d(0, 0, 0);

    /// @brief number of detector unit which has the size v3_det_length
    double n_unit = 0;

    /// @brief exposure time (days);
    double days = 0;

    /// @brief after constructor, this value will be assigned
    int n_element = 0;

    /// @brief vector of vector of DetectorElement
    /// @details be careful of index order, z.at(iy).at(ix), col-major
    std::vector<std::vector<DetectorElement>> vec_vec_DetectorElement = {};

    /// @brief Rotation matrix type
    angle_util::Rotation3dType rotation_type = angle_util::Rotation3dType::LOCAL;

    /// @brief Grid2dBinGroup::Parameters
    Grid2dBinGroup::Parameters prm_bingrp = Grid2dBinGroup::Parameters();

    /// @brief angle unit shared by all DetectorElements in this panel
    DetectorElement::AngleUnit angle_unit = DetectorElement::AngleUnit::Tangent;

  public:
    //=========================================================================
    /// @name constructor_destructor
    ///@{

    /// @brief default constructor
    DetectorPanel() = default;

    /// @brief copy constructor
    DetectorPanel(const DetectorPanel &org) = default;

    /// @brief move constructor, some members are not moved.
    DetectorPanel(DetectorPanel&& other) noexcept;

    /// @brief destructor
    ~DetectorPanel() = default;

    /// @brief Detector panel constructor
    /// @param detprm DetectorPanel::Parameters
    /// @param detid Detid of this instance
    /// @param uqid_min_tmp Uqid of the first DetectorElement in this DetectorPanel
    /// @param uqid_max_tmp Uqid of the last DetectorElement in this DetectorPanel
    /// @param vec_UqidInfo vector of UqidInfo, which is used to upper layer to build
    /// @param dic DetectorIndexContainer, just to set the pointer
    /// @note Uses OpenMP
    DetectorPanel( const DetectorPanel::Parameters &detprm
                  , const Detid detid
                  , const Uqid uqid_min_tmp, const Uqid uqid_max_tmp
                  , std::vector<UqidInfo>& vec_UqidInfo
                  , DetectorIndexContainer &dic);

    ///@} ---------------------------------------------------------------------

    /// @brief allocate memory for std::vector<std::vector<DetectorElement>> and vec_vec_XXX of g2bg_
    void vec_vec_memory_allocate();

    /// @brief reserve memory for DetectorElement::vec_tf_in_PL
    /// @note Uses OpenMP
    void reserve_vec_tf_in_PL(const int n_reserve_vec_tf_in_length);

    //==================================================================
    /// @name operators
    ///@{

    /// @brief Assignment operator operator=, copy assignment operator using copy & swap idiom
    DetectorPanel& operator=(const DetectorPanel& other) = default;

    /// @brief Inequality operator
    /// @note Name is not compared.
    bool operator!=(const DetectorPanel& other) const;

    /// @brief Equality operator
    /// @details Name is not compared. Defined via the inequality operator.
    bool operator==(const DetectorPanel& other) const {
      return !(*this != other);
    };

    /// @brief Move assignment
    DetectorPanel& operator=(DetectorPanel&& other) noexcept;

    ///@} ------------------------------------------------------------------

    //=========================================================================
    /// @name wrapper_getter_functions
    ///@{
    
    /// @brief get x axis, wrapper for Grid2d::get_x_axis()
    const Grid1d& get_x_axis() const { return g2bg_.get_x_axis(); };

    /// @brief get y axis, wrapper for Grid2d::get_y_axis()
    const Grid1d& get_y_axis() const { return g2bg_.get_y_axis(); };

    /// @brief get nbinx, wrapper for Grid2d::get_nbinx()
    int get_nbinx() const { return g2bg_.get_nbinx(); };

    /// @brief get nbiny, wrapper for Grid2d::get_nbiny()
    int get_nbiny() const { return g2bg_.get_nbiny(); };

    /// @brief get x_interval, wrapper for Grid2d::get_x_interval()
    double get_x_interval() const { return g2bg_.get_x_interval(); };

    /// @brief get y_interval, wrapper for Grid2d::get_y_interval()
    double get_y_interval() const { return g2bg_.get_y_interval(); };

    /// @brief get_ix, wrapper for Grid2d::get_ix()
    int get_ix( const double value ) const { return g2bg_.get_ix(value); };

    /// @brief get_iy, wrapper for Grid2d::get_iy()
    int get_iy( const double value ) const { return g2bg_.get_iy(value); };

    ///@} ------------------------------------------------------------------


    //=========================================================================
    /// @name getter_functions
    ///@{
    
    /// @brief return detid
    Detid get_detid() const { return detid_; };

    /// @brief get mutable reference of g2bg_
    Grid2dBinGroup& call_g2bg() { return g2bg_;};

    /// @brief get const reference of g2bg_
    const Grid2dBinGroup& get_g2bg() const { return g2bg_;};

    /// @brief get angle unit defined for this detector panel
    DetectorElement::AngleUnit get_angle_unit() const { return angle_unit; };

    /// @brief return uqid in this detector
    Inthis get_id_in_this_detector(const int ix, const int iy) const{
      return getDetectorElement(ix,iy).get_id_in_this_detector();
    };

    const DetectorIndexContainer& get_dic() const {
      if (!pdic_){
        LOG_ERROR("DetectorPanel::dic_ is not set.");
        throw std::runtime_error("DetectorPanel::dic_ is not set.");
      }
      return *pdic_;
    };

    DetectorIndexContainer& call_dic() {
      if (!pdic_){
        LOG_ERROR("DetectorPanel::dic_ is not set.");
        throw std::runtime_error("DetectorPanel::dic_ is not set.");
      }
      return *pdic_;
    };

    /// @brief get the signal sum of igroup, wrapper for Grid2dBinGroup::get_signal_group()
    double get_signal_group(const Igroup igroup) const{
      return g2bg_.get_signal_group(igroup);
    };

    /// @brief get the noise sum of igroup, wrapper for Grid2dBinGroup::get_noise_poi_group()
    double get_noise_poi_group(const Igroup igroup) const{
      return g2bg_.get_noise_poi_group(igroup);
    };

    /// @brief get the deterministic floor noise sum (det) of igroup, wrapper for Grid2dBinGroup::get_noise_det_group()
    double get_noise_det_group(const Igroup igroup) const{
      return g2bg_.get_noise_det_group(igroup);
    };

    /// @brief get the signal sum of igroup with poisson error, wrapper for Grid2dBinGroup::get_signal_group_poisson()
    double get_signal_group_poisson(const Igroup igroup) const{
      return g2bg_.get_signal_group_poisson(igroup);
    };

    /// @brief get the noise sum of igroup with poisson error, wrapper for Grid2dBinGroup::get_noise_poi_group_poisson()
    double get_noise_poi_group_poisson(const Igroup igroup) const{
      return g2bg_.get_noise_poi_group_poisson(igroup);
    };

    /// @brief get the name of this instance
    std::string get_name() const { return name; };

    /// @brief get position and direction of detector as Ray3d
    Ray3d get_ray3d() const { return ray3d; };

    /// @brief get direction of detector of x,y as Ray2d
    Ray2d get_ray2d() const { return ray3d.toRay2d(); }; 

    /// @brief get direction of detector as Eigen::Vector3d
    Eigen::Vector3d get_v3_dir() const { return ray3d.dir(); };

    /// @brief get position of detector as Eigen::Vector3d
    Eigen::Vector3d get_v3_pos() const { return ray3d.pos(); };

    /// @brief get size of detector as Eigen::Vector3d
    /// @details [0]=horizontal, [1]=vertical, [2]=depth
    Eigen::Vector3d get_v3_det_length() const { return v3_det_length; };

    /// @brief get number of detector unit
    double get_n_unit() const { return n_unit; };

    /// @brief get exposure time (days) 
    double get_days() const { return days; };

    /// @brief get exposure time (sec)
    double get_exposure_time_sec() const { return days * 86400.0; };

    /// @brief get number of detector element
    int get_n_element() const { return n_element; };

    /// @brief get reference of Grid2dBinGroup::Parameters prm_bingrp
    const Grid2dBinGroup::Parameters& get_prm_bingrp() const { return prm_bingrp; };

    /// @brief get the immutable reference of DetectorElement of ix, iy
    /// @param[in] ix index along x-axis
    /// @param[in] iy index along y-axis
    /// @details the member of DetectorElement cannot be changed
    /// @remark It is an error without const at the front.
    const DetectorElement& getDetectorElement( const int ix, const int iy ) const;

    /// @brief get the immutable reference of DetectorElement of tx, ty
    /// @param[in] tx position along x-axis
    /// @param[in] ty position along y-axis
    /// @details the member of DetectorElement cannot be changed
    /// @remark It is an error without const at the front.
    const DetectorElement& getDetectorElement( const double tx, const double ty ) const;

    /// @brief Get the summed signal over angle-bin regions defined by BinList
    Eigen::VectorXf get_vecxf_nmuoned() const;

    /// @brief get set of uqid from igroup
    std::set<Uqid> get_set_uqid_by_igroup( const Igroup igroup ) const;

    /// @brief get vector of uqid from igroup
    /// @note Uses OpenMP
    std::vector<Uqid> get_vec_uqid_by_igroup( const Igroup igroup ) const;

    /// @brief return std::vector<std::set<unique_index>>, the index of vector is igroup \n 
    /// with transforming ixiy into uqid
    std::vector<std::set<int>> get_vec_set_uqid() const;

    /// @brief get Eigen::VectorXf of PL, index is id_in_this_detector
    /// @note Uses OpenMP
    Eigen::VectorXf mp_get_vecxf_PL() const;

    /// @brief get Eigen::VectorXf of DL, index is id_in_this_detector
    /// @note Uses OpenMP
    Eigen::VectorXf mp_get_vecxf_DL() const;

    /// @brief get the signal values of all DetectorElement as Eigen::VectorXf
    Eigen::VectorXf get_vecxf_signal() const;

    /// @brief get the noise values of all DetectorElement as Eigen::VectorXf
    Eigen::VectorXf get_vecxf_noise() const;

    /// @brief get the diagonal SpMatf of eff_cnt with index of (id_inthis,id_inthis)
    SpMatf get_spmat_eff_cnt() const;

    /// @brief get the diagonal SpMatf of eff(random sampling value) with index of (id_inthis,id_inthis)
    SpMatf get_spmat_eff_sample() const;

    /// @brief get the total signal of all DetectorElement DL >= DL_thres
    /// @note Uses OpenMP
    double mp_calc_total_signal_DLthres(const double DL_thres) const;

    /// @brief get the total SOT of all DetectorElement as double
    /// @note Uses OpenMP
    double mp_calc_total_SOT_DLthres(const double DL_thres) const;

    /// @brief get_subtract of *this and other DetectorPanel
    /// @note Uses OpenMP
    DetectorPanel get_subtract(const DetectorPanel &other) const;

    ///@} ---------------------------------------------------------------------

    //=========================================================================
    // setter_functions
    ///@{

    /// @brief set x_axis, wrapper for Grid2d::set_x_axis()
    /// @param x_axis_in
    void set_x_axis( const Grid1d &x_axis_in ){ g2bg_.set_x_axis(x_axis_in); };

    /// @brief set y_axis, wrapper for Grid2d::set_y_axis()
    /// @param y_axis_in
    void set_y_axis( const Grid1d &y_axis_in ){ g2bg_.set_y_axis(y_axis_in); };

    /// @brief set DetectorIndexContainer pointer
    /// @param dic_ptr 
    void set_dic_ptr(DetectorIndexContainer* dic_ptr) { pdic_ = dic_ptr; };

    /// @brief set detid
    void set_detid( const Detid detid_in ){ detid_ = detid_in; };
    
    /// @brief set n_element
    void set_n_element( const int n_element_in ){ n_element = n_element_in; };

    /// @brief set name   
    void set_name( const std::string name_in ){ name = name_in; };

    /// @brief set direction
    void set_direction( const Eigen::Vector3d v3_direction_in ){
      ray3d.set_dir( v3_direction_in );
    };

    /// @brief set position
    void set_position( const Eigen::Vector3d v3_position_in ){
      ray3d.set_pos( v3_position_in );
    };

    /// @brief set detector size
    void set_v3_det_length( const Eigen::Vector3d v3_det_length_in ){
      v3_det_length = v3_det_length_in;
    };

    /// @brief set the exposure time (sec) for all DetectorElement
    /// @note Uses OpenMP
    void set_exposure_time_sec( const double time_sec_in );

    /// @brief copy only DetectortElement::vec_tf_in_PL
    /// @note Uses OpenMP
    void mp_set_vec_tf_in_PL( const DetectorPanel& panel );

    /// @brief Assign penetrating_muonflux and signal muon counts to each DetectorElement in vec_vec_DetectorElement
    /// @details openmp version
    /// @param[in] ft FluxTable that has penetrating muon flux table
    /// @param[in] tf_apply_eff_cnt If true, central count uses eff_cnt: signal = calc_signal()*eff_cnt
    ///            (deterministic, NOT the dice path). Default false keeps behavior unchanged.
    void mp_calc_set_peneflux_signal_from_DL(
      const FluxTable &ft, const bool tf_apply_eff
    , const bool tf_apply_eff_cnt = false );


    /// @brief calc and set the noise of each DetectorElement, split into a
    ///        deterministic floor part (noise_det) and a Poisson-bucket part (noise_poi).
    /// @details Each of the two noise sources (flux-proportional, SOT-proportional)
    ///          contributes to both the floor and the Poisson bucket:
    ///          noise_det = signal*flux_proport_ratio_floor
    ///                    + SOT_term*SOT_proport_noise_ratio_floor;
    ///          noise_poi = signal*flux_proport_ratio_poisson
    ///                    + SOT_term*SOT_proport_noise_ratio_poisson;
    ///          where SOT_term = SOT/total_SOT_DL_thres * total_signal_DL_thres.
    /// @param[in] flux_proport_ratio_floor flux-proportional ratio, deterministic floor part
    /// @param[in] flux_proport_ratio_poisson flux-proportional ratio, Poisson-bucket part
    /// @param[in] SOT_proport_noise_ratio_floor SOT-proportional (angle-independent) ratio, floor part
    /// @param[in] SOT_proport_noise_ratio_poisson SOT-proportional (angle-independent) ratio, Poisson-bucket part
    /// @param[in] user_defined_noise_flux_ratio User-defined noise flux ratio (currently not implemented; reserved for future extension to support custom noise models from external data)
    /// @param[in] vec_path_user_defined_noise_flux vector of path to user defined noise flux files
    /// @param[in] DL_thres density length threshold for SOT calculation
    /// @note Uses OpenMP
    void mp_set_noise(
      const double flux_proport_ratio_floor
    , const double flux_proport_ratio_poisson
    , const double SOT_proport_noise_ratio_floor
    , const double SOT_proport_noise_ratio_poisson
    , const double user_defined_noise_flux_ratio
    , const std::vector<fs::path> vec_path_user_defined_noise_flux
    , const double DL_thres);

    /// @brief initilize PL for all DetectorElement
    /// @note Uses OpenMP
    void mp_initPL( const double value_in=0.0 );
    
    /// @brief initilize DL for all DetectorElement
    /// @note Uses OpenMP
    void mp_initDL( const double value_in=0.0 );

    /// @brief initilize PL and PL for all DetectorElement
    /// @param PL_in rock length
    /// @param DL_in density length
    /// @note Uses OpenMP
    void mp_init_PLDL( const double PL_in=0.0, const double DL_in=0.0 );

    /// @brief initilize penetrating muon flux for all DetectorElement
    /// @note Uses OpenMP
    void mp_init_peneflux( const double value_in=0.0 );

    /// @brief initilize signal for all DetectorElement
    /// @note Uses OpenMP
    void mp_init_signal( const double value_in=0.0 );

    /// @brief initilize noise for all DetectorElement
    void mp_init_noise( const double value_in=0.0 );

    /// @brief initilize len, DL, pene, sig, noi for all DetectorElement
    /// @note Uses OpenMP
    void mp_init_PL_DL_pene_sig_noi( const double PL_in=0.0
      , const double DL_in=0.0, const double peneflux_in=0.0
      , const double sig_in=0.0, const double noi_in=0.0 );

    /// @brief set Grid2dBinGroup::Parameters
    void set_parameters( const Grid2dBinGroup::Parameters &prm_bingrp_in ){
      prm_bingrp = prm_bingrp_in;
    };

    /// @brief set signal/noise of each DetectorElement to Grid2dBinGroup vec_vec_signal/noise;
    /// @param signal_init initial value of signal
    /// @param noise_init initial value of noise
    /// @param is_avail_init initial value of is_avail
    /// @param PL_thres rocklength threshold
    /// @param signal_under_thres signal under threshold
    /// @param noise_under_thres noise under threshold
    /// @param is_avail_under_thres is_avail under threshold
    /// @return void
    /// @note Uses OpenMP
    void copy_signal_noise_to_g2bg(
        const double signal_init = 0.0, const double noise_init = 0.0
      , const bool is_avail_init = true
      , const double PL_thres = 0.0
      , const double DL_thres = 0.0
      , const double signal_under_thres = 0.0
      , const double noise_under_thres = 0.0
      , const bool is_avail_under_thres = false );

    /// @brief Determine a unique igroup directly from ix, iy
    /// @param igroup_start first group id
    /// @return int last igroup
    int assign_naive_igroup( const int igroup_start );

    /// @brief assign_1st_igroup, specify the first igroup.
    /// @param igroup_start first group id
    /// @param nx_div_init number of bin-group divisions in x direction
    /// @param ny_div_init number of bin-group divisions in y direction
    /// @return int last igroup
    int assign_1st_igroup(
        const int igroup_start
      , const int nx_div_init, const int ny_div_init );

    /// @brief Perform bin grouping as defined in an ASCII file.
    /// @note Uses OpenMP
    void grouping_by_bin_list( const fs::path &path_in );

    /// @brief copy Grid2dBinGroup, wrapper for \ref Grid2dBinGroup::copy_g2bg()
    /// @param panel_in : other DetectorPanel
    void copy_g2bg( const DetectorPanel& panel_in ){
      this->g2bg_.copy_g2bg( panel_in.g2bg_.get_g2bg() );
    };

    /// @brief copy bimap_ from other DetectorPanel
    /// @param other : other DetectorPanel
    void copy_bimap_from( const DetectorPanel &other ){
      g2bg_.copy_bimap_from( other.g2bg_ );
      g2bg_.set_done_grouping( other.g2bg_.is_done_grouping() );
    };

    ///@} ---------------------------------------------------------------------

    //=========================================================================
    /// @name call_functions
    /// @brief call the pointer of DetectorElement
    /// @remark the member of DetectorElement can be changed
    ///@{

    /// @brief  get the mutable reference of DetectorElement of ix, iy
    /// @details the member of DetectorElement can be changed
    DetectorElement& callDetectorElement(const int ix, const int iy);

    /// @brief   get the mutable reference of DetectorElement of tx, ty
    /// @details the member of DetectorElement can be changed
    DetectorElement& callDetectorElement(const double tx, const double ty);
    ///@} ---------------------------------------------------------------------

    //======================================================================
    /// @name calc_density_functions
    ///@{

    /// @brief apply calc_approx_volue to DetectorElement(ix,iy)
    double calc_approx_volume( const int ix, const int iy ) const;

    /// @brief apply calc_approx_volue to DetectorElement(Ixiy &ixiy)
    double calc_approx_volume( const Ixiy &ixiy ) const{
      return calc_approx_volume( ixiy[0], ixiy[1] );
    };

    /// @brief Calculate calc_approx_volume for {(ix,iy)} belonging to igroup and return the sum
    double calc_approx_volume_sum_grouped( const Igroup igroup ) const;

    /// @brief Calculate nmuon_grouped for a given density dens
    /// @param dens : density
    /// @param ft : FluxTable
    /// @return nmuon_grouped
    double calc_nmuon_grouped(const Igroup igroup
      , const double dens, const FluxTable &ft ) const;

    /// @brief Calculate nmuon_grouped for a given density dens
    /// @note openmp version
    /// @param dens : density
    /// @param ft : FluxTable
    /// @return nmuon_grouped
    double mp_calc_nmuon_grouped(const Igroup igroup
      , const double dens, const FluxTable &ft ) const;

    /// @brief invalid density value
    static constexpr double not_found_dens = -1.0;

    /// @brief invalid delta_nmuon value
    static constexpr double not_found_delta_nmuon = -1.0;

    /// @brief If no bin exists for the corresponding angle, return not_found_delta_nmuon. lower_delta_nmuon, center_delta_nmuon, upper_delta_nmuon
    static constexpr std::array<double,3>
      not_found_delta_nmuon_arr = {not_found_delta_nmuon, not_found_delta_nmuon, not_found_delta_nmuon};
    
    /// @brief Compute density grid-search results for nmuon_grouped
    /// @param igroup : group id
    /// @param dens_min : minimum of search range
    /// @param dens_max : maximum of search range
    /// @param dens_step : step of search range
    /// @param ft : FluxTable
    /// @return std::vector<std::array<double,2>> : {dens,nmuon_grouped}
    /// @note Uses OpenMP
    std::vector<std::array<double,2>>
      calc_vec_dens_nmuon_grouped(const Igroup igroup, const double dens_min, const double dens_max
        , const double dens_step, const FluxTable &ft ) const;

    /// @brief Perform adaptive grid search to find dens that minimizes the difference between calc_nmuon_grouped and target_value.
    /// @param igroup        group ID
    /// @param dens_min      minimum of search range
    /// @param dens_max      maximum of search range
    /// @param vec_dens_step vector of step sizes used in each iteration (sorted in descending order)
    /// @param range_factor  scaling factor when updating the search range (e.g., expansion rate relative to initial step size)
    /// @param ft          FluxTable
    /// @param target_value  target value (e.g., signal_group_center or upper/lower values)
    /// @return {best_dens, best_error}: dens that gives minimum error and that error
    std::array<double, 2> 
      calc_dens_nmuon_grouped(
        const Igroup igroup, 
        const double dens_min, const double dens_max, 
        const std::vector<double> &vec_dens_step, const double range_factor,
        const FluxTable &ft, 
        const double target_value);
    
    /// @brief Perform adaptive grid search to find dens that minimizes the difference between calc_nmuon_grouped and target_value.
    /// @param igroup        group ID
    /// @param dens_min      minimum of search range
    /// @param dens_max      maximum of search range
    /// @param initial_step  initial grid interval (coarse step)
    /// @param min_step      minimum step size that will not be refined further
    /// @param refinement_factor  step reduction factor per iteration (e.g., 2 means half)
    /// @param num_refinements    number of refinements
    /// @param ft          FluxTable
    /// @param target_value  target value (e.g., signal_group_center or upper/lower values)
    /// @return {best_dens, best_error}: dens that gives minimum error and that error
    std::array<double, 2> calc_dens_nmuon_grouped_adaptive_search(
      const Igroup igroup, 
      const double dens_min, const double dens_max, 
      const double initial_step, const double min_step,
      const double refinement_factor, 
      const int num_refinements,
      const FluxTable &ft, 
      const double target_value);

    /// @brief Calculate projected mean density for grouped bins and set vec_dens_group_lower, vec_dens_group_center, vec_dens_group_upper.
    /// @note Density is calculated from dens_min to dens_max in steps of dens_step.
    /// @param tf_signal_poisson whether to apply Poisson error to signal counts
    /// @note noise term is always det floor + Poisson(poi bucket); no noise switch.
    /// @param igroup : group id
    /// @param dens_min : minimum density
    /// @param dens_max : maximum density
    /// @param dens_step : step of density
    /// @param ft : FluxTable
    void calc_set_proj_dens_grouped(
        const bool tf_signal_poisson
      , const Igroup igroup, const double dens_min, const double dens_max
      , const double dens_step, const double sigma
      , const FluxTable &ft
      , const bool tf_eff = false, const bool tf_eff_independent = false );

    /// @brief Obtain projected density values for a group via adaptive grid search and set them to member variables.
    /// @param tf_signal_poisson whether to apply Poisson error to signal counts
    /// @note noise term is always det floor + Poisson(poi bucket); no noise switch.
    /// @param igroup        group ID
    /// @param dens_min      minimum of search range
    /// @param dens_max      maximum of search range
    /// @param vec_dens_step vector of step sizes used in each iteration (sorted in descending order)
    /// @param range_factor  scaling factor when updating the search range (e.g., expansion rate relative to initial step size)
    /// @param ft          FluxTable
    void calc_set_proj_dens_grouped(
      const bool tf_signal_poisson
    , const Igroup igroup, const double sigma
    , const double dens_min, const double dens_max
    , const std::vector<double> &vec_dens_step, const double range_factor
    , const FluxTable &ft
    , const bool tf_eff = false, const bool tf_eff_independent = false );


          /// @brief Perform adaptive grid search to find dens that minimizes the difference between calc_nmuon_grouped and target_value.
    /// @param igroup        group ID
    /// @param dens_min      minimum of search range
    /// @param dens_max      maximum of search range
    /// @param initial_step  initial grid interval (coarse step)
    /// @param min_step      minimum step size that will not be refined further
    /// @param refinement_factor  step reduction factor per iteration (e.g., 2 means half)
    /// @param num_refinements    number of refinements
    /// @param ft          FluxTable
    void calc_set_proj_dens_grouped_adaptive_search(
      const Igroup igroup
    , const double sigma , const double dens_min, const double dens_max
    , const double dens_step_init, const double dens_step_min
    , const double dens_step_refinement_factor, const int dens_step_num_refinements
    , const FluxTable &ft
    , const bool tf_eff = false, const bool tf_eff_independent = false );

    /// @brief Calculate projected mean density for all grouped bins and set vec_dens_group_lower, vec_dens_group_center, vec_dens_group_upper.
    /// @param tf_signal_poisson whether to apply Poisson error to signal counts
    /// @note noise term is always det floor + Poisson(poi bucket); no noise switch.
    /// @param dens_min : minimum density
    /// @param dens_max : maximum density
    /// @param dens_step : step of density
    /// @param ft : FluxTable
    /// @note Uses OpenMP
    void mp_calc_set_proj_dens_grouped(
        const bool tf_signal_poisson
      , const double dens_min, const double dens_max
      , const double dens_step, const double sigma
      , const FluxTable &ft
      , const bool tf_eff = false, const bool tf_eff_independent = false );

    /// @brief Calculate density for grouped bins and set vec_dens_group_lower, vec_dens_group_center, vec_dens_group_upper.
    /// @param dens_min      minimum of search range
    /// @param dens_max      maximum of search range
    /// @param vec_dens_step vector of step sizes used in each iteration (sorted in descending order)
    /// @param range_factor  scaling factor when updating the search range (e.g., expansion rate relative to initial step size)
    /// @param ft          FluxTable
    /// @note Uses OpenMP
    void mp_calc_set_proj_dens_grouped(
        const bool tf_signal_poisson
      , const double dens_min, const double dens_max
      , const std::vector<double> &vec_dens_step, const double range_factor
      , const double sigma, const FluxTable &ft
      , const bool tf_eff = false, const bool tf_eff_independent = false );

    /// @brief Calculate density for grouped bins and set vec_dens_group_lower, vec_dens_group_center, vec_dens_group_upper.
    /// @note adative search version
    /// @note Uses OpenMP
    void mp_calc_set_proj_dens_grouped_adaptive_search(
      const double dens_min, const double dens_max
    , const double dens_step_init, const double dens_step_min
    , const double dens_step_refinement_factor, const int dens_step_num_refinements
    , const double sigma, const FluxTable &ft
    , const bool tf_eff = false, const bool tf_eff_independent = false );

    /// @brief calc_signal_group without efficiency
    /// @note Uses OpenMP
    double mp_calc_signal_noeff_group(const Igroup igroup) const;

    /// @brief calc_signal_group with efficiency
    /// @note Uses OpenMP
    double mp_calc_signal_eff_group(const Igroup igroup) const;

    /// @brief Simple unweighted mean of the element efficiency triplets over one group.
    /// @details eff_g^X = (1/n_g) * sum_{i in g} eff_i^X for X = low, cnt, upp.
    ///          Deterministic (no random sampling) and independent of any density model;
    ///          it depends only on the assigned element efficiencies and the grouping.
    /// @param igroup group ID
    /// @return array of 3 elements (eff_low, eff_cnt, eff_upp)
    /// @throws std::runtime_error If igroup is out of range or the group is empty.
    std::array<double,3> calc_eff_mean_group(const Igroup igroup) const;

    /// @brief Efficiency-uncertainty variance of the grouped count, in count^2 units.
    /// @details Mirrors DetectorPanelArray::get_vecxf_var_eff_all but panel-local for one group.
    ///          Accumulates sigma_eff_i * b_i over the merge-group elements, where
    ///          sigma_eff_i = 0.5*(eff_upp - eff_low) and b_i = calc_signal() (efficiency-free
    ///          base count). Returns (sum)^2 for the common mode, or sum of squares for the
    ///          independent mode. Same quantity added to the C_N diagonal (id-6bgzem); here it is
    ///          added in quadrature to the Poisson term of the projected-density error band.
    /// @param igroup        group ID
    /// @param tf_independent true: sum of squares (independent); false: square of sum (common)
    /// @return efficiency variance var_eff_g (count^2)
    double calc_var_eff_group(const Igroup igroup, const bool tf_independent) const;

    ///@} ------------------------------------------------------------------

    //=========================================================================
    /// @name read_write_functions
    ///@{

    /// @brief out_header
    void out_header( FILE *fout = stderr ) const ;

    /// @brief out_all information to FILE *fout
    void out_all( FILE *fout = stderr ) const ;

    /// @brief output all information to file
    void out_all( const fs::path& pathout ) const ;

    /// @brief output txcnt, tycnt, signal to FILE *fout
    void out_axay_signal( FILE *fout ) const;

    /// @brief output txcnt, tycnt, signal to fs::path
    void out_axay_signal( const fs::path& pathout ) const;

    /// @brief Output (tx, ty, arbitrary data) for all DetectorElement
    void out_txtyData(
      const std::string& suffix,
      std::function<double(const DetectorElement&)> data_extractor,
      const std::string& prefix = "" ) const;

    /// @brief Output (tx, ty, Signal) for all DetectorElement
    void out_txtySignal() const;

    /// @brief Output (tx, ty, Signal) for all DetectorElement
    /// @param prefix : prefix of output file
    void out_txtySignal( const std::string& prefix ) const;

    /// @brief Output (tx, ty, PL) for all DetectorElement
    void out_txtyPL() const;

    /// @brief Output (tx, ty, PL) for all DetectorElement
    /// @param prefix : prefix of output file
    void out_txtyPL( const std::string& prefix ) const;

    /// @brief Output (tx, ty, DL) for all DetectorElement
    void out_txtyDL() const;

    /// @brief Output (tx, ty, DL) for all DetectorElement
    /// @param prefix : prefix of output file
    void out_txtyDL( const std::string& prefix ) const;

    /// @brief interface structure for efficiency table to read
    struct EffBin {
      double xlow, xup;
      double ylow, yup;
      double eff_low, eff_cnt, eff_upp;
    };

    /// @brief read efficinecy table from file
    /// @param path_in : file path
    /// @note Legacy input path; mp_assign_efficiency_model replaces it when a
    ///       model JSON5 (path_eff_model) is given.
    void read_efficiency_table( const fs::path &path_in );

    /// @brief Evaluate an efficiency model at every angular bin center and
    ///        assign eff_low / eff_cnt / eff_upp to each DetectorElement
    /// @param[in] model Loaded efficiency model (EfficiencyModel::get_tf_loaded() must be true)
    /// @throws std::runtime_error If the model is not loaded
    /// @note Uses OpenMP. Each (ix, iy) iteration writes only to its own
    ///       DetectorElement, so there is no data race.
    /// @note The angular grid comes from this panel's own Grid2dBinGroup
    ///       (g2bg_), so no grid mismatch with an external table can occur.
    void mp_assign_efficiency_model( const EfficiencyModel &model );

    ///@} ---------------------------------------------------------------------
    
    //=========================================================================
    /// @name checker_functions
    ///@{

    /// @brief wrapper for Grid2dBinGroup::check_ix_inside
    void check_ix_inside( const int ix ) const { g2bg_.check_ix_inside(ix); };

    /// @brief wrapper for Grid2dBinGroup::check_iy_inside
    void check_iy_inside( const int iy ) const { g2bg_.check_iy_inside(iy); };

    ///@} ---------------------------------------------------------------------

    //=========================================================================
    /// @name depth reso sweeper functions
    ///@{

    /// @brief get the sum of signal in the specified range. for depth vs reso sweep
    /// @param[in] tx_lower lower bound of tx
    /// @param[in] tx_upper upper bound of tx
    /// @param[in] ty_lower lower bound of ty
    /// @param[in] ty_upper upper bound of ty
    /// @return sum of signal, noise
    /// @note if tx/ty/min/max is near the edge, include the edge bins
    /// @note Uses OpenMP
    SignalNoiseSum mp_get_signal_noise_sum_range_xy(
        const double tx_lower,const double tx_upper
      , const double ty_lower,const double ty_upper ) const;

    /// @brief get the sum of signal in the specified range (y-sweep)
    /// @param[in] tx_lower lower bound of tx
    /// @param[in] tx_upper upper bound of tx
    /// @param[in] ty_lower lower bound of ty
    /// @param[in] ty_upper upper bound of ty
    /// @param[in] ty_min minimum of ty for sweeping
    /// @param[in] ty_max maximum of ty for sweeping
    /// @param[in] ty_step step of ty for sweeping
    /// @param[in] v3_pos_obj_top position of the top of the object
    /// @return vector of SignalNoiseDepth for each ty step
    std::vector<SignalNoiseDepth> mp_get_signal_noise_sum_y_sweep(
        const double tx_lower, const double tx_upper
      , const double ty_lower, const double ty_upper
      , const double ty_min, const double ty_max, const double ty_step
      , const Eigen::Vector3d& v3_pos_obj_top
      , const double diff_elev ) const;

    /// @brief Cut out a DetectorPanel within the specified range
    /// @details Cut g2bg_ with the same arguments as Grid2dBinGroup::cut,
    ///          and return a new instance with contained DetectorElement re-arranged and re-indexed.
    DetectorPanel cut(const double x_lower, const double x_upper
                     , const double y_lower, const double y_upper
                     , const double x_eps = 1.0e-6
                     , const double y_eps = 1.0e-6) const;

    /// @brief set direction from Eigen::Vector2d v2_pos
    /// @note Uses OpenMP
    void set_direction_to_v2_pos( const Eigen::Vector2d v2_pos );

    ///@} ---------------------------------------------------------------------

    //=========================================================================
    /// @name binary_io_functions
    ///@{

    /// @brief save vec_vec_DetectorElement to std::ofstream &ofs
    void save_vec_vec_DetectorElement( std::ofstream &ofs ) const;

    /// @brief load vec_vec_DetectorElement from std::ifstream &ifs
    void load_vec_vec_DetectorElement( std::ifstream &ifs );

    /// @brief save all variables of DetectorPanel to std::ofstream &ofs
    void save( std::ofstream &ofs ) const;

    /// @brief load all variables of DetectorPanel from std::ifstream &ifs
    void load( std::ifstream &ifs );

    /// @brief save all variables of DetectorPanel to fs::path
    void save( const fs::path& pathout ) const;

    /// @brief load all variables of DetectorPanel from fs::path
    void load( const fs::path &path_in );
    ///@} ---------------------------------------------------------------------
};
