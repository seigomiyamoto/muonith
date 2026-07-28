/// @file cls_Grid2dBinGroup.hpp
/// @brief Mainly used for grouping DetectorElement bins in Grid2d to increase muon signal statistics.
/// @details Extends Grid2d to support grouped bins for coarser resolution analysis.
#pragma once

#include <map>
#include <set>
#include <fstream>
#include <iostream>
#include <sstream> // istringstream
#include <string>
#include <cstdio>
#include <cmath>
#include <functional>  //for sorting
#include <algorithm>//for sorting
#include <vector>
#include <filesystem> // for std::filesystem::path
#include <iomanip>
#include <stdexcept>

#include <Eigen/Dense>
#include "cls_Grid2d.hpp"
#include "cls_Grid2dXYZ.hpp"
#include "cls_DetectorIndexContainer.hpp"
#include "ns_type_definitions.hpp"
using namespace index_type_definitions;

#include "cls_OneToManyBimap.hpp"

//#########################################################################
//#########################################################################
/// @class Grid2dBinGroup
/// @brief Class that records, in a map, which igroup each (ix, iy) belongs to, \n
/// used when grouping finely divided angular bins
/// @ingroup basicGridClasses
//#########################################################################
//#########################################################################
class Grid2dBinGroup : public Grid2d {
  public:
    class Parameters; // forward declaration

    /// @brief static constexpr members for the values below
    static constexpr bool is_avail_vecvec_init = false;
    static constexpr double signal_init_vecvec = 0.0;
    static constexpr double noise_init_vecvec = 0.0;
    static constexpr bool is_avail_group_init = false;
    static constexpr double signal_group_init = 0.0;
    static constexpr double noise_group_init = 0.0;
    static constexpr double dens_group_lower_init = 0.0;
    static constexpr double dens_group_center_init = 0.0;
    static constexpr double dens_group_upper_init = 0.0;
    static constexpr std::array<double,3> arr_dens_group_init = {
      dens_group_lower_init,
      dens_group_center_init,
      dens_group_upper_init
    };
    static constexpr double delta_nmuon_group_lower_init = 0.0;
    static constexpr double delta_nmuon_group_center_init = 0.0;
    static constexpr double delta_nmuon_group_upper_init = 0.0;
    static constexpr double volume_group_init = 0.0;
    static constexpr double eff_low_group_init = -1.0;
    static constexpr double eff_cnt_group_init = -1.0;
    static constexpr double eff_upp_group_init = -1.0;

  private:
    /// @brief name of this instance
    std::string name = "Grid2dBinGroup";

    /// @brief Map recording which group each (ix, iy) belongs to
    /// @details Container holding, for each bin of Grid2d, the group it is bundled into \n
    /// group0 : (0,0) , (0,1) , (1,0), (1,1), .... \n
    /// group1 : (2,0) , (2,1) , (3,0), (3,1), .... \n
    /// ....

    using BimapIgroupIxiy_IyMajor = OneToManyBimap<
      Igroup, Ixiy, std::less<Igroup>, IyMajorCompare
    >;

    BimapIgroupIxiy_IyMajor bimap_;

    /// @brief Detid of this instance
    Detid detid_ = DetidNotAssigned;
    
    /// @brief DetectorIndexContainer pointer from DetectorPanelArray::dic_
    DetectorIndexContainer* pdic_ = nullptr;

    /// @brief Availability flag of each bin group.
    /// @details If at least one bin in the group is false, the group is treated as not available.
    std::vector<bool> vec_is_avail_group;
    
    /// @brief Number of signal muons contained in each bin group
    /// @details Stores observed data and its pseudo data.
    std::vector<double> vec_signal_group;

    /// @brief Number of noise events in each bin group (Poisson-fluctuated side, poi)
    /// @details Stores observed data and its pseudo data.
    std::vector<double> vec_noise_poi_group;

    /// @brief Number of deterministic floor noise events in each bin group (non-fluctuated side, det)
    /// @details Consumer side adds this directly: noise = noise_det_group + Poisson(noise_poi_group).
    ///          For all configs with floor ratios = 0 this stays 0 and behavior is unchanged.
    std::vector<double> vec_noise_det_group;

    /// @brief Number of signal muons in vec_signal_group with a Poisson error applied
    std::vector<double> vec_signal_group_poisson;

    /// @brief Number of noise muons in vec_noise_poi_group with a Poisson error applied
    std::vector<double> vec_noise_poi_group_poisson;
 
    /// @brief Availability flag of each bin.
    /// @details Stored in ColMajor order
    std::vector<std::vector<bool>> vec_vec_is_avail;

    /// @brief Storage of the signal held by each finest bin.
    /// @details Stored in ColMajor order
    std::vector<std::vector<double>> vec_vec_signal;

    /// @brief Storage of the noise held by each finest bin (Poisson-fluctuated side, poi).
    /// @details Stored in ColMajor order
    std::vector<std::vector<double>> vec_vec_noise_poi;

    /// @brief Storage of the deterministic floor noise held by each finest bin (non-fluctuated side, det).
    /// @details Stored in ColMajor order. Initialized to 0; stays 0 whenever the floor ratio is 0.
    std::vector<std::vector<double>> vec_vec_noise_det;

    /// @brief is done grouping by any function
    bool done_grouping = false;
    
    /// @brief lower of average projection desntiy of each group
    std::vector<double> vec_dens_group_lower;

    /// @brief average projection desntiy of each group
    std::vector<double> vec_dens_group_center;

    /// @brief upper of average projection desntiy of each group
    std::vector<double> vec_dens_group_upper;

    /// @brief delta_nmuon when evaluate lower average projection desntiy of each group
    std::vector<double> vec_delta_nmuon_group_lower;

    /// @brief delta_nmuon when evaluate center average projection desntiy of each group
    std::vector<double> vec_delta_nmuon_group_center;

    /// @brief delta_nmuon when evaluate upper average projection desntiy of each group
    std::vector<double> vec_delta_nmuon_group_upper;

    /// @brief sum of volume of each group
    std::vector<double> vec_volume_group;

    /// @brief vec_eff_low_group
    std::vector<double> vec_eff_low_group;

    /// @brief vec_eff_cnt_group
    std::vector<double> vec_eff_cnt_group;

    /// @brief vec_eff_upp_group
    std::vector<double> vec_eff_upp_group;

    /// @brief Template helper function for get_xxx_group methods
    /// @details Reduces code duplication for getter functions with error checking
    /// @param vec Vector from which to get the value
    /// @param igroup Group index
    /// @param init_val Initial value to return if vec is empty
    /// @param name Name of the vector (for error messages)
    template<typename T>
    T get_group_value(const std::vector<T>& vec, const Igroup igroup,
                      const T init_val, const char* name) const {
      if (vec.empty()) {
        LOG_WARN("{} is empty, cannot get value for igroup {}", name, igroup);
        return init_val;
      }
      if (igroup < 0 || igroup >= static_cast<int>(vec.size())) {
        LOG_ERROR("igroup {} is out of range for {}", igroup, name);
        THROW_ERROR(std::string("igroup out of range for ") + name);
      }
      return vec.at(igroup);
    };

    /// @brief Template helper function for resize_vec_xxx_group methods
    /// @details Reduces code duplication for resize functions
    /// @param vec Vector to resize
    /// @param n_group Number of groups to resize to
    /// @param init_val Initial value to fill the vector with
    template<typename T>
    void resize_group_vector(std::vector<T>& vec, const int n_group, const T init_val) {
      vec.clear();
      vec.resize(n_group, init_val);
    };

  public:
    //======================================================================
    /// @name constructor_destructor
    ///@{
    
    /// @brief default constructors
    Grid2dBinGroup() = default;
    
    /// @brief copy constructor
    Grid2dBinGroup(const Grid2dBinGroup &org) = default;

    /// @brief move constructor
    Grid2dBinGroup(Grid2dBinGroup &&other) noexcept = default;
    
    /// @brief destructor
    ~Grid2dBinGroup() = default;

    /// @brief constructor from Grid1d x_axis & y_axis
    Grid2dBinGroup(const Grid1d &x_axis_in, const Grid1d &y_axis_in);
    
        /// @brief   constructor from Grid2dXYZ
    /// @details Assigns one single igroup to every bin
    Grid2dBinGroup(const Grid2dXYZ &g2xyz_in, const int igroup = 0);
    
    /// @brief   constructor from Grid2dXYZ & g2xyz_signal_in
    /// @details Divides the bins by nx_div, ny_div and assigns nx_div*ny_div igroups
    Grid2dBinGroup(const Grid2dXYZ &g2xyz_signal_in,
                   const int nx_div, const int ny_div);

    /// @brief it calls read_bin_list2
    Grid2dBinGroup(const std::filesystem::path &path_in);
    
    ///@} ------------------------------------------------------------------

    /// @brief memory allocation of vec_signal_group, vec_noise_poi_group,
    ///        vec_is_avail_group, vec_vec_signal, vec_vec_noise_poi,
    ///        vec_vec_is_avail, vec_dens_group_center, vec_dens_group_upper,
    ///        vec_dens_group_lower
    /// @param signal_init Initial value for signal vectors
    /// @param noise_init Initial value for noise vectors
    /// @param is_avail_init Initial value for availability flags
    /// @details Allocates and initializes all group and bin-level vectors
    ///          required for signal/noise tracking and density calculations.
    void init_vec_vec( const double signal_init = 0.0
      , const double noise_init = 0.0, const bool is_avail_init = false );

    /// @brief allocate memory for vec_value_group
    void allocate_vec_value_group();

    //======================================================================
    /// @name operators
    ///@{
    
    /// @brief assignment operator
    Grid2dBinGroup& operator=(const Grid2dBinGroup& other) = default;

    /// @brief inequality operator
    /// @details The name is not compared.
    bool operator!=(const Grid2dBinGroup& other) const;

    /// @brief equality operator
    /// @details The name is not compared.
    bool operator==(const Grid2dBinGroup& other) const {
      return !(*this != other);
    };
    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name getter_functions
    ///@{
    
    /// @brief  get name of this instance
    std::string get_name() const { return name; };

    /// @brief get the immutable ref of this class
    const Grid2dBinGroup& get_g2bg_ref() const { return *this; };
    
    /// @brief get the number of igroup in map, returns the number of unique igroup
    int get_n_group() const { return static_cast<int>(bimap_.countOne()); };

    /// @brief get the vector of igroup in this instance
    std::vector<Igroup> get_vec_igroup() const { return bimap_.get_vecOne(); };

    /// @brief get igroup from ixiy
    /// @param ix X-axis bin index
    /// @param iy Y-axis bin index
    /// @return Group index associated with the bin (ix, iy)
    Igroup get_igroup(const int ix, const int iy) const;

    /// @brief get the number of (ix,iy) in the bin group_index=igroup from map
    int get_n_uqid(const Igroup igroup) const {
      return static_cast<int>(bimap_.get_vecMany(igroup).size());
    };

    /// @brief  get the name of this instance
    Detid get_detid() const { return detid_; }

    /// @brief get the const reference of the DetectorIndexContainer
    const DetectorIndexContainer& get_dic() const {
      if (!pdic_) THROW_ERROR("Grid2dBinGroup::pdic_ is not set.");
      return *pdic_;
    };

    /// @brief get the mutable reference of the DetectorIndexContainer
    DetectorIndexContainer& call_dic() {
      if (!pdic_) THROW_ERROR("Grid2dBinGroup::pdic_ is not set.");
      return *pdic_;
    };

    /// @brief get the immutable ref of this class
    const Grid2dBinGroup& get_g2bg() const { return *this; };

    /// @brief get signal of vec_vec_signal.at(iy).at(ix)
    /// @param ix X-axis bin index
    /// @param iy Y-axis bin index
    /// @return Signal value at the specified bin
    double get_signal( const int ix, const int iy ) const;

    /// @brief get signal of vec_vec_signal.at(iy).at(ix)
    /// @param ixiy Combined bin index (ix, iy)
    /// @return Signal value at the specified bin
    double get_signal( const Ixiy& ixiy ) const;

    /// @brief get signal of vec_vec_signal.at(iy).at(ix) from x,y
    /// @param x X-coordinate in physical space
    /// @param y Y-coordinate in physical space
    /// @return Signal value at the bin containing (x, y)
    double get_signal( const double x, const double y ) const;

    /// @brief get signal vec_vec_noise_poi.at(iy).at(ix)
    /// @param ix X-axis bin index
    /// @param iy Y-axis bin index
    /// @return Noise value at the specified bin
    double get_noise( const int ix, const int iy ) const;

    /// @brief get signal vec_vec_noise_poi.at(iy).at(ix)
    /// @param ixiy Combined bin index (ix, iy)
    /// @return Noise value at the specified bin
    double get_noise( const Ixiy& ixiy ) const;

    /// @brief get signal vec_vec_noise_poi.at(iy).at(ix) from x,y
    /// @param x X-coordinate in physical space
    /// @param y Y-coordinate in physical space
    /// @return Noise value at the bin containing (x, y)
    double get_noise( const double x, const double y ) const;

    /// @brief get deterministic floor noise vec_vec_noise_det.at(iy).at(ix)
    /// @param ix X-axis bin index
    /// @param iy Y-axis bin index
    /// @return Floor (det) noise value at the specified bin
    double get_noise_det( const int ix, const int iy ) const;

    /// @brief get the signal sum in the bin group_index=id_group
    double get_signal_group( const int id_group ) const {
       return vec_signal_group.at(id_group);
    };

    /// @brief get the noise sum in the bin group_index=id_group (Poisson-bucket side, poi)
    double get_noise_poi_group( const int id_group ) const {
       return vec_noise_poi_group.at(id_group);
    };

    /// @brief get the deterministic floor noise sum in the bin group_index=id_group (det)
    double get_noise_det_group( const int id_group ) const {
       return vec_noise_det_group.at(id_group);
    };

    /// @brief get signal_group_poisson in the bin group_index=id_group
    double get_signal_group_poisson( const int id_group ) const{
      return vec_signal_group_poisson.at(id_group);
    };

    /// @brief get noise_poi_group_poisson in the bin group_index=id_group
    double get_noise_poi_group_poisson( const int id_group ) const{
      return vec_noise_poi_group_poisson.at(id_group);
    };

    /// @brief get the reference of vec_signal_group
    const std::vector<double>& get_vec_signal_group() const { return vec_signal_group; };

    /// @brief get the reference of vec_noise_poi_group
    const std::vector<double>& get_vec_noise_poi_group() const { return vec_noise_poi_group; };

    /// @brief get the vector of the number of muons in all the bin group
    Eigen::VectorXf get_vecxf_nmuon( ) const;    
    
    /// @brief push_back nmuon to vec_signal_group
    void push_back_nmuon( const double nmuon_in ){
      vec_signal_group.push_back( nmuon_in );
    };

    /// @brief clear all the signal_group
    void clear_vec_signal_group(){ vec_signal_group.clear(); };

    /// @brief clear all the noise_poi_group
    void clear_vec_noise_poi_group(){ vec_noise_poi_group.clear(); };

    /// @brief clear all the volume_group
    void clear_vec_volume_group(){ vec_volume_group.clear(); };

    /// @brief get_xmin_xmax_ymin_ymax in the group_index=igroup
    /// @param igroup Group index
    /// @return Array [xmin, xmax, ymin, ymax] defining the bounding box
    ///         of all bins in the group
    std::array<double,4> get_xmin_xmax_ymin_ymax(const int igroup) const;

    /// @brief return vec_vec_ixiy as std::vector<std::vector<Ixiy>> from igroup
    std::vector<std::vector<Ixiy>> 
      get_vec_vec_ixiy( const Igroup igroup ) const;

    /// @brief convert std::vector<std::vector<Ixiy>> vec_vec_ixiy into std::vector<Ixiy> vec_ixiy
    std::vector<Ixiy> get_vec_ixiy( 
      const std::vector<std::vector<Ixiy>> &vec_vec_ixiy ) const;

    /// @brief return the igroup with the smallest signal_group + noise_poi_group. Throws if vec_igroup is empty
    Igroup get_igroup_signal_noise_group_min( const std::vector<int> vec_igroup ) const;

    /// @brief return the igroup with the largest signal_group + noise_poi_group. Throws if vec_igroup is empty
    Igroup get_igroup_signal_noise_group_max( const std::vector<int> vec_igroup ) const;

    /// @brief get ixmin, ixmax, iymin, iymax from igroup
    /// @param igroup Group index
    /// @return Array [ixmin, ixmax, iymin, iymax] defining the bounding box
    ///         in bin index space
    std::array<int,4> get_ixmin_ixmax_iymin_iymax( const Igroup igroup ) const;

    /// @brief return the number of ix and iy in igroup_in
    Ixiy get_ixiylen( const int igroup_in ) const;

    /// @brief return the number of ix in igroup_in
    int get_ixlen( const int igroup_in ) const;

    /// @brief return the number of iy in igroup_in
    int get_iylen( const int igroup_in ) const;

    /// @brief return the vector of (ix, iy) adjacent to the right side of the vec_ixiy forming igroup
    std::vector<Ixiy> 
      get_adjacent_plus_x_vec_ixiy(
        const int igroup_in, const int ix_delta, const int iy_delta ) const;

    /// @brief return the igroup of the bin with the smallest signal_group
    Igroup get_igroup_signal_group_min() const;

    /// @brief get igroup from ixiy
    Igroup get_igroup( const Ixiy &ixiy ) const { return get_igroup(ixiy[0],ixiy[1]); };

    /// @brief get igroup set from bimap_
    std::set<Igroup> get_set_igroup() const;

    /// @brief get the vector of Ixiy in igroup
    /// @param igroup Group index
    /// @return Vector of all (ix, iy) bin indices belonging to the group
    std::vector<Ixiy> get_vec_ixiy(const Igroup igroup) const;

    /// @brief get value of vec_vec_is_avail for ix,iy
    bool get_is_avail( const int ix, const int iy ) const{
      return vec_vec_is_avail.at(iy).at(ix);
    };

    /// @brief get vector of ixiy for all igroup
    /// @details The outer vector is ordered by igroup, and the inner vector
    /// is the array of Ixiy belonging to that igroup.
    std::vector< std::vector< Ixiy > > get_vec_igroup_vec_ixiy() const;

    /// @brief get the const reference of vec_vec_signal
    const std::vector<std::vector<double>>& 
      get_vec_vec_signal() const { return vec_vec_signal; };

    /// @brief get the const reference of vec_vec_noise_poi
    const std::vector<std::vector<double>>& 
      get_vec_vec_noise_poi() const { return vec_vec_noise_poi; };

    /// @brief get the set of ix,iy
    /// @param igroup Group index
    /// @return Set of all (ix, iy) bin indices belonging to the group
    /// @details Unlike get_vec_ixiy, returns a set for fast membership testing
    std::set<Ixiy> get_set_ixiy( const Igroup igroup ) const;

    /// @brief return done_grouping
    bool is_done_grouping() const { return done_grouping; };

    /// @brief get the size of vec_signal_group
    int get_size_vec_signal_group() const {
      return vec_signal_group.size(); };

    /// @brief get the size of vec_noise_poi_group
    int get_size_vec_noise_poi_group() const {
      return vec_noise_poi_group.size(); };
    
    /// @brief get the size of vec_signal_group_poisson
    int get_size_vec_signal_group_poisson() const {
      return vec_signal_group_poisson.size(); };
    
    /// @brief get the size of vec_noise_poi_group_poisson
    int get_size_vec_noise_poi_group_poisson() const {
      return vec_noise_poi_group_poisson.size(); };

    /// @brief get the size of vec_dens_group_lower
    int get_size_vec_dens_group_lower() const {
      return vec_dens_group_lower.size(); };

    /// @brief get the size of vec_dens_group_center
    int get_size_vec_dens_group_center() const {
      return vec_dens_group_center.size(); };
    
    /// @brief get the size of vec_dens_group_upper
    int get_size_vec_dens_group_upper() const {
      return vec_dens_group_upper.size(); };
    
    /// @brief get the size of vec_delta_nmuon_group_lower
    int get_size_vec_delta_nmuon_group_lower() const {
      return vec_delta_nmuon_group_lower.size(); };
    
    /// @brief get the size of vec_delta_nmuon_group_center
    int get_size_vec_delta_nmuon_group_center() const {
      return vec_delta_nmuon_group_center.size(); };
    
    /// @brief get the size of vec_delta_nmuon_group_upper
    int get_size_vec_delta_nmuon_group_upper() const {
      return vec_delta_nmuon_group_upper.size(); };
    
    /// @brief get the size of vec_volume_group
    int get_size_vec_volume_group() const {
      return vec_volume_group.size(); };

    /// @brief get the vec_dens_group_lower of igroup
    double get_dens_group_lower( const Igroup igroup ) const;

    /// @brief get the vec_dens_group_center of igroup
    double get_dens_group_center( const Igroup igroup ) const;

    /// @brief get the vec_dens_group_upper of igroup
    double get_dens_group_upper( const Igroup igroup ) const;

    /// @brief get the vec_delta_nmuon_group_lower of igroup
    double get_delta_nmuon_group_lower( const Igroup igroup ) const;

    /// @brief get the vec_delta_nmuon_group_center of igroup
    double get_delta_nmuon_group_center( const Igroup igroup ) const;

    /// @brief get the vec_delta_nmuon_group_upper of igroup
    double get_delta_nmuon_group_upper( const Igroup igroup ) const;

    /// @brief get the vec_volume_group of igroup
    double get_volume_group( const Igroup igroup ) const;

    /// @brief get the vec_eff_low_group of igroup
    double get_eff_low_group( const Igroup igroup ) const;

    /// @brief get the vec_eff_cnt_group of igroup
    double get_eff_cnt_group( const Igroup igroup ) const;

    /// @brief get the vec_eff_upp_group of igroup
    double get_eff_upp_group( const Igroup igroup ) const;

    /// @brief get the subtracted vec_vec_signal , this - other
    /// @note Uses OpenMP
    std::vector<std::vector<double>>
      get_vec_vec_signal_subtracted( const Grid2dBinGroup &other ) const;

    /// @brief get the subtracted vec_vec_noise_poi , this - other
    /// @note Uses OpenMP
    std::vector<std::vector<double>>
      get_vec_vec_noise_poi_subtracted( const Grid2dBinGroup &other ) const;
    
    /// @brief get the tuple of subtracted vec_vec_signal, vec_vec_noise_poi
    /// @note Uses OpenMP
    std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>>
      get_vec_vec_signal_noise_subtracted( const Grid2dBinGroup &other ) const;

    /// @brief get the subtracted vec_signal_group , this - other
    /// @note Uses OpenMP
    std::vector<double>
      get_vec_signal_group_subtracted( const Grid2dBinGroup &other ) const;
    
    /// @brief get the subtracted vec_noise_poi_group , this - other
    /// @note Uses OpenMP
    std::vector<double>
      get_vec_noise_poi_group_subtracted( const Grid2dBinGroup &other ) const;

    /// @brief get the tuple of subtracted vec_signal_group, vec_noise_poi_group
    /// @note Uses OpenMP
    std::tuple<std::vector<double>, std::vector<double>>
      get_vec_signal_noise_group_subtracted( const Grid2dBinGroup &other ) const;

    /// @brief create a new Grid2dBinGroup restricted to the specified range
    /// @param x_lower Lower bound of x-range to extract
    /// @param x_upper Upper bound of x-range to extract
    /// @param y_lower Lower bound of y-range to extract
    /// @param y_upper Upper bound of y-range to extract
    /// @param x_eps Tolerance for x-boundary comparison
    /// @param y_eps Tolerance for y-boundary comparison
    /// @return New Grid2dBinGroup instance containing only bins within
    ///         specified range
    /// @details Creates a new instance by extracting bins that fall within
    ///          [x_lower, x_upper] × [y_lower, y_upper]. Epsilon parameters
    ///          handle floating-point comparison tolerance.
    Grid2dBinGroup cut(const double x_lower, const double x_upper,
                       const double y_lower, const double y_upper,
                       const double x_eps = 1.0e-6,
                       const double y_eps = 1.0e-6) const;


    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name read_write_functions
    ///@{

    /// @brief read angular bin group list from ascii file of style 2
    /// @param path_in Path to input file containing bin group definitions
    /// @details Reads a text file where each line defines a bin group with
    ///          format: igroup ix iy. Populates the bimap_ with group
    ///          assignments.
    void read_bin_list2( const std::filesystem::path &path_in );

    /// @brief disp all group range
    void disp_group_range_all() const;

    /// @brief output all igroup, ix, iy
    void out_igroup_ixiy_all(FILE* fout) const;

    /// @brief output all igroup, ix, iy to a text file
    void out_igroup_ixiy_all(const std::filesystem::path& pathout) const;

    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name setter_functions
    ///@{
    
    /// @brief set name of this instance
    void set_name( const std::string &name_in ){ name = name_in; };

    /// @brief set detid of this instance
    void set_detid(const Detid d){ detid_ = d; }

    /// @brief set the pointer to the DetectorIndexContainer
    void set_dic_ptr(DetectorIndexContainer* dic_ptr) { pdic_ = dic_ptr; };

    /// @brief set ixiy to group of *p_dic
    void set_ixiy_to_group(const Ixiy& ixiy, const Igroup igroup){
      bimap_.insert(igroup, ixiy); };

    /// @brief copy_Grid2dBinGroup
    void copy_g2bg( const Grid2dBinGroup &g2bg_in ){ *this = g2bg_in; };

    /// @brief copy_bimap_from
    /// @details copy bimap_ from another Grid2dBinGroup
    void copy_bimap_from( const Grid2dBinGroup &other ){
      bimap_.copy_from(other.bimap_);
    };

    /// @brief set vec_vec_signal
    void set_vec_vec_signal( const std::vector<std::vector<double>> &vec_vec_signal_in ){
      vec_vec_signal = vec_vec_signal_in; };

    /// @brief set vec_vec_noise_poi
    void set_vec_vec_noise_poi( const std::vector<std::vector<double>> &vec_vec_noise_poi_in ){
      vec_vec_noise_poi = vec_vec_noise_poi_in; };

    /// @brief set vec_vec_is_avail
    void set_vec_vec_is_avail( const std::vector<std::vector<bool>> &vec_vec_is_avail_in ){
      vec_vec_is_avail = vec_vec_is_avail_in; };

    /// @brief set vec_signal_group
    void set_vec_signal_group( const std::vector<double> &vec_signal_group_in ){
      vec_signal_group = vec_signal_group_in; };
    
    /// @brief set vec_noise_poi_group
    void set_vec_noise_poi_group( const std::vector<double> &vec_noise_poi_group_in ){
      vec_noise_poi_group = vec_noise_poi_group_in; };

    /// @brief set vec_noise_det_group (deterministic floor noise per group, det)
    void set_vec_noise_det_group( const std::vector<double> &vec_noise_det_group_in ){
      vec_noise_det_group = vec_noise_det_group_in; };

    /// @brief set signal=signal_in value of bin ix,iy
    void set_signal( const int ix, const int iy, const double signal_in );

    /// @brief set signal=signal_in value of bin (ix,iy)
    void set_signal( const Ixiy &tp_ixiy, const double signal_in );

    /// @brief set signal=signal_in value of bin (x,y)
    void set_signal( const double x, const double y, const double signal_in );

    /// @brief set noise=noise_in value of bin ix,iy
    void set_noise( const int ix, const int iy, const double noise_in );

    /// @brief set noise=noise_in value of bin (ix,iy)
    void set_noise( const Ixiy &tp_ixiy, const double noise_in );

    /// @brief set noise=noise_in value of bin (x,y)
    void set_noise( const double x, const double y, const double noise_in );

    /// @brief set deterministic floor noise (det) = noise_det_in value of bin ix,iy
    void set_noise_det( const int ix, const int iy, const double noise_det_in );

    /// @brief allocate memory for vec_is_avail_group
    void resize_vec_is_avail_group(const int n_group_in, const bool is_avail_init);

    /// @brief allocate memory for vec_signal_group
    void resize_vec_signal_group(const int n_group_in, const double signal_init);

    /// @brief allocate memory for vec_noise_poi_group
    void resize_vec_noise_poi_group(const int n_group_in, const double noise_init);

    /// @brief allocate memory for vec_noise_det_group (deterministic floor noise per group, det)
    void resize_vec_noise_det_group(const int n_group_in, const double noise_init);

    /// @brief allocate memory for vec_signal_group_poisson
    void resize_vec_signal_group_poisson(const int n_group_in, const double signal_init);

    /// @brief allocate memory for vec_noise_poi_group_poisson
    void resize_vec_noise_poi_group_poisson(const int n_group_in, const double noise_init);

    /// @brief allocate memory for vec_dens_group_lower
    void resize_vec_dens_group_lower(const int n_group_in, const double dens_init);

    /// @brief allocate memory for vec_dens_group_center
    void resize_vec_dens_group_center(const int n_group_in, const double dens_init);

    /// @brief allocate memory for vec_dens_group_upper
    void resize_vec_dens_group_upper(const int n_group_in, const double dens_init);

    /// @brief allocate memory for vec_delta_nmuon_group_lower
    void resize_vec_delta_nmuon_group_lower(const int n_group_in, const double delta_nmuon_init);

    /// @brief allocate memory for vec_delta_nmuon_group_center
    void resize_vec_delta_nmuon_group_center(const int n_group_in, const double delta_nmuon_init);

    /// @brief allocate memory for vec_delta_nmuon_group_upper
    void resize_vec_delta_nmuon_group_upper(const int n_group_in, const double delta_nmuon_init);

    /// @brief allocate memory for vec_volume_group
    void resize_vec_volume_group(const int n_group_in, const double volume_group_init);

    /// @brief allocate memory for vec_eff_low_group
    void resize_vec_eff_low_group(const int n_group_in, const double eff_low_init);

    /// @brief allocate memory for vec_eff_cnt_group
    void resize_vec_eff_cnt_group(const int n_group_in, const double eff_cnt_init);

    /// @brief allocate memory for vec_eff_upp_group
    void resize_vec_eff_upp_group(const int n_group_in, const double eff_upp_init);
        
    /// @brief set is_avail of bin ix,iy
    void set_is_avail( const int ix, const int iy, const bool is_avail_in );

    /// @brief set is_avail of bin (ix,iy)
    void set_is_avail( const Ixiy &tp_ixiy, const bool is_avail_in );

    /// @brief set is_avail of bin (x,y)
    void set_is_avail(const double x, const double y, const bool is_avail_in );

    /// @brief set is_avail for all bins in igroup
    void set_is_avail( const Igroup igroup, const bool is_avail_in );

    /// @brief set the is_avail_group in the bin group_index=id_group
    void set_is_avail_group( const int id_group, const bool is_avail_group_in ){
      vec_is_avail_group.at(id_group) = is_avail_group_in;
    };

    /// @brief set the signal_group in the bin group_index=id_group
    void set_signal_group( const int id_group, const double signal_group_in ){
      vec_signal_group.at(id_group) = signal_group_in;
    };

    /// @brief set the noise_poi_group in the bin group_index=id_group (Poisson-bucket side, poi)
    void set_noise_poi_group( const int id_group, const double noise_group_in ){
      vec_noise_poi_group.at(id_group) = noise_group_in;
    };

    /// @brief set the deterministic floor noise_poi_group (det) in the bin group_index=id_group
    void set_noise_det_group( const int id_group, const double noise_det_group_in ){
      vec_noise_det_group.at(id_group) = noise_det_group_in;
    };

    /// @brief set the signal_group_poisson in the bin group_index=id_group
    void set_signal_group_poisson( const int id_group, const double signal_group_poisson_in ){
      vec_signal_group_poisson.at(id_group) = signal_group_poisson_in;
    };

    /// @brief set the noise_poi_group_poisson in the bin group_index=id_group
    void set_noise_poi_group_poisson( const int id_group, const double noise_poi_group_poisson_in ){
      vec_noise_poi_group_poisson.at(id_group) = noise_poi_group_poisson_in;
    };

    /// @brief apply_calc_signal_group_to_vec_signal
    void apply_calc_signal_group_to_vec_signal();

    /// @brief set done_grouping
    void set_done_grouping( const bool done_grouping_in ){ done_grouping = done_grouping_in; };

    /// @brief set the vec_dens_group_lower of igroup
    void set_dens_group_lower( const Igroup igroup, const double dens_group_lower_in )
      { vec_dens_group_lower.at(igroup) = dens_group_lower_in; };

    /// @brief set the vec_dens_group_center of igroup
    void set_dens_group_center( const Igroup igroup, const double dens_group_center_in )
      { vec_dens_group_center.at(igroup) = dens_group_center_in; };
    
    /// @brief set the vec_dens_group_upper of igroup
    void set_dens_group_upper( const Igroup igroup, const double dens_group_upper_in )
      { vec_dens_group_upper.at(igroup) = dens_group_upper_in; };

    /// @brief set the vec_delta_nmuon_group_lower of igroup
    void set_delta_nmuon_group_lower( const Igroup igroup, const double delta_nmuon_group_lower_in )
      { vec_delta_nmuon_group_lower.at(igroup) = delta_nmuon_group_lower_in; };

    /// @brief set the vec_delta_nmuon_group_center of igroup
    void set_delta_nmuon_group_center( const Igroup igroup, const double delta_nmuon_group_center_in )
      { vec_delta_nmuon_group_center.at(igroup) = delta_nmuon_group_center_in; };
    
    /// @brief set the vec_delta_nmuon_group_upper of igroup
    void set_delta_nmuon_group_upper( const Igroup igroup, const double delta_nmuon_group_upper_in )
      { vec_delta_nmuon_group_upper.at(igroup) = delta_nmuon_group_upper_in; };

    /// @brief set the vec_volume_group of igroup
    void set_volume_group( const Igroup igroup, const double volume_group_in )
    { vec_volume_group.at(igroup) = volume_group_in; };

    /// @brief set the vec_eff_low_group of igroup
    void set_eff_low_group( const Igroup igroup, const double eff_low_group_in )
    { vec_eff_low_group.at(igroup) = eff_low_group_in; };

    /// @brief set the vec_eff_cnt_group of igroup
    void set_eff_cnt_group( const Igroup igroup, const double eff_cnt_group_in )
    { vec_eff_cnt_group.at(igroup) = eff_cnt_group_in; };

    /// @brief set the vec_eff_upp_group of igroup
    void set_eff_upp_group( const Igroup igroup, const double eff_upp_group_in )
    { vec_eff_upp_group.at(igroup) = eff_upp_group_in; };

    /// @brief renumber the registered igroup in ascending order as 0, 1, 2, ...
    void reassign_igroup_indices( const Igroup igroup_start = 0 );

    ///@} ------------------------------------------------------------------
    

    //======================================================================
    /// @name checker_functions
    ///@{

    /// @brief get is_avail of bin ix,iy
    bool is_avail( const int ix, const int iy ) const;

    /// @brief get is_avail of bin (ix,iy)
    bool is_avail( const Ixiy &tp_ixiy ) const;

    /// @brief get is_avail of bin (x,y)
    bool is_avail( const double x, const double y ) const;

    /// @brief output the ix, iy whose is_avail is false to the logger
    void log_non_avail_vec_ixiy( const spdlog::level::level_enum& level, const Igroup igroup ) const;

    /// @brief return false if is_avail is false for at least one finest bin, true if all are true
    bool is_avail_group( const Igroup igroup ) const;

    /// check that pdic_ is not nullptr
    void _ensure_dic(const DetectorIndexContainer* p_dic) const {
      if(!p_dic) THROW_ERROR("Grid2dBinGroup::pdic_ is nullptr");
    };

    /// @brief check tf_done_grouping
    void check_done_grouping() const {
      if (!done_grouping) THROW_ERROR("Grid2dBinGroup::done_grouping is false");
    };

    ///@} ------------------------------------------------------------------

    /// @brief set all is_avail to is_avail_in
    void clear_is_avail(const bool is_avail_in);

    //======================================================================
    /// @name signal_noise_group_functions
    ///@{

    /// @brief calculate the signal_group of igroup.
    /// @param igroup Group index
    /// @return Sum of signal values over all bins in the group
    /// @details Sums vec_vec_signal values for all (ix,iy) belonging to igroup
    double calc_signal_group( const Igroup igroup ) const;

    /// @brief calculate the signal_group of the group to which the bin (ix,iy) belongs.
    double calc_signal_group( const Ixiy &tp_ixiy ) const {
      int igroup = get_igroup(tp_ixiy);
      return calc_signal_group(igroup);
    };

    /// @brief calculate the signal_group of the group to which the bin ix,iy belongs.
    double calc_signal_group( const int ix, const int iy ) const {
      return calc_signal_group({ix,iy});
    };

    /// @brief Calculate sum of signal+noise over a rectangular bin region.
    /// @param ixmin Minimum x-index (inclusive)
    /// @param ixmax Maximum x-index (inclusive)
    /// @param iymin Minimum y-index (inclusive)
    /// @param iymax Maximum y-index (inclusive)
    /// @return Sum of (signal + noise) over the rectangular region
    /// @note Direct access to vec_vec_signal/vec_vec_noise_poi without bimap_ lookup.
    double calc_signal_noise_rect(int ixmin, int ixmax, int iymin, int iymax) const;

    /// @brief Result of calc_trial_split_max_half.
    struct TrialSplitResult {
      int max_half;   ///< 0-based division index with larger signal+noise (0 or 1).
      int n_halves;   ///< Number of valid halves (1 if too small to split, 2 if split succeeded).
    };

    /// @brief Determine which half has larger signal+noise after a trial 2-split.
    /// @param ixmin Minimum x-index of the group (inclusive)
    /// @param ixmax Maximum x-index of the group (inclusive)
    /// @param iymin Minimum y-index of the group (inclusive)
    /// @param iymax Maximum y-index of the group (inclusive)
    /// @param split_x If true, split along X axis; otherwise split along Y axis.
    /// @return TrialSplitResult with max_half and n_halves.
    /// @note Reproduces split2DArrayIntoGrid logic without copying Grid2dBinGroup.
    TrialSplitResult calc_trial_split_max_half(
      int ixmin, int ixmax, int iymin, int iymax, bool split_x) const;

    /// @brief calculate the noise_poi_group of igroup.
    /// @param igroup Group index
    /// @return Sum of noise values over all bins in the group
    /// @details Sums vec_vec_noise_poi values for all (ix,iy) belonging to igroup
    double calc_noise_poi_group( const Igroup igroup ) const;

    /// @brief calculate the deterministic floor noise_det_group (det) of igroup.
    /// @param igroup Group index
    /// @return Sum of floor (det) noise values over all bins in the group
    /// @details Sums vec_vec_noise_det values for all (ix,iy) belonging to igroup
    double calc_noise_det_group( const Igroup igroup ) const;

    /// @brief calculate the noise_poi_group of the group to which the bin (ix,iy) belongs.
    double calc_noise_poi_group( const Ixiy &tp_ixiy ) const {
      const Igroup igroup = get_igroup(tp_ixiy);
      return calc_noise_poi_group(igroup);
    };

    /// @brief calculate the noise_poi_group of the group to which the bin ix,iy belongs.
    double calc_noise_poi_group( const int ix, const int iy ) const {
      return calc_noise_poi_group({ix,iy});
    };

    ///@} ------------------------------------------------------------------



    /// @brief disp ixmin_ixmax_iymin_iymax
    void disp_ixmin_ixmax_iymin_iymax( const Igroup igroup ) const;

    /// @brief disp ixmin_ixmax_iymin_iymax
    void disp_ixmin_ixmax_iymin_iymax( FILE* fout, const Igroup igroup ) const;
    
    //======================================================================
    /// @name map_operation_functions
    ///@{

    //========================================================
    // map related
    //========================================================

    /// @brief clear bimap_
    void clear_map(){ bimap_.clear(); };

    /// @brief insert ixiy to map
    void insert_map( const int igroup, const Ixiy &tp_ixiy ){
      bimap_.insert(igroup, tp_ixiy);
    };
    
    /// @brief insert ixiy to map
    void insert_map( const int igroup, const int ix, const int iy ){
      insert_map(igroup,{ix,iy});
    };    
    
    /// @brief insert vector of ixiy to map
    void insert_map( const int igroup, const std::vector< Ixiy > &vec_ixiy ){
      for(const auto& ixiy : vec_ixiy) insert_map(igroup, ixiy);
    };

    //==================================================================
    // merge bin related
    //==================================================================

    /// @brief merge the single bin at the +ix_delta, +iy_delta neighbor
    /// @brief Merge bin at (ix+ix_delta, iy+iy_delta) into the group
    ///        containing tp_ixiy
    /// @param tp_ixiy Reference bin position (ix, iy)
    /// @param ix_delta X-offset to the bin being merged
    /// @param iy_delta Y-offset to the bin being merged
    /// @return true if merge succeeded, false otherwise
    /// @details Updates bimap_ after merge. Current position is tp_ixiy.
    bool merge_bin(
      const Ixiy &tp_ixiy, const int ix_delta, const int iy_delta );

    /// @brief merge the bins of ig_src into ig_dst.
    /// @param ig_src Source group index (will be merged into ig_dst)
    /// @param ig_dst Destination group index
    /// @return true if merge succeeded, false otherwise
    bool merge_group( const int ig_src, const int ig_dst );

    /// @brief search around the seed group and merge bins until z_sum reaches nmuon_thres or above
    /// @details The merged bins always form a rectangle.
    bool merge_rect_bingroup(
      const int igroup_seed, const double nmuon_thres, const int nloop_max=200 );

    /// @brief erase the elements with igroup = igroup_del from bimap_.
    void remove_group(const Igroup igroup_del){ bimap_.eraseOne(igroup_del); };

    /// @brief erase the element with ix,iy = tp_ixiy from bimap_.
    void remove_ixiy( const Ixiy &tp_ixiy ){ bimap_.eraseMany(tp_ixiy); };

    /// @brief   divide igroup_in by nx_div, ny_div.
    /// @details The igroup of the newly created bins starts at the maximum igroup + 1. \n
    /// igroup_in is not erased at the end. \n
    /// Returns the vector of the newly created igroup, ix_divide, iy_divide.
    std::vector<std::array<int,3>> 
      divide(const int igroup_in, const int nx_div, const int ny_div);

    /// @brief   new version of divide
    /// @details Divides igroup_in by nx_div, ny_div. \n
    /// The igroup of the newly created bins starts at the maximum igroup + 1. \n
    /// The remaining bins are gathered into the bin specified by ix_extra, iy_extra. \n
    /// Returns the vector of the newly created entries.
    std::vector<std::array<int,3>> 
      divide_extra( const Igroup igroup_in
        , const int nx_div, const int ny_div
        , const int ix_extra, const int iy_extra
        , const bool exec_move_ixiy);


    /// @brief   divide igroup.
    /// @details Repeats divide until signal_group becomes smaller than signal_noise_group_trigger. \n
    /// divide is applied in the order (nx_div=2, ny_div=1), then (nx_div=1, ny_div=2). \n
    /// The remainder is given to the side with the smaller signal_group. \n
    /// nloop_max is the maximum number of divide loops. \n
    /// See the image below for the procedure \n
    /// \image html images/auto_divide_by_signal_noise_group.draw.png
    /// @param tf_prefer_split_x : when xlen==ylen, split along the x direction if true.
    /// @returns tuple_bool_int | bool: true if divided, false if not divided, int: number of loops performed
    std::tuple<bool,int> auto_divide_by_signal_noise_group(
      const Igroup igroup, const double signal_noise_group_trigger, const int nloop_max,
      const int ixlen_min, const int iylen_min, const bool tf_prefer_split_x );
    
    /// @brief   divide igroup.
    /// @details Repeats divide until signal_group becomes smaller than signal_noise_group_trigger. \n
    /// divide is applied in the order (nx_div=2, ny_div=1), then (nx_div=1, ny_div=2). \n
    /// The remainder is given to the side with the smaller signal_group. \n
    /// nloop_limit is the maximum number of divide loops.
    /// @param prm_bingrp : Grid2dBinGroup::Parameters
    /// @returns tuple_bool_int | bool: true if divided, false if not divided, int: number of loops performed
    std::tuple<bool,int> auto_divide_by_signal_noise_group(
      const Igroup igroup, const Grid2dBinGroup::Parameters &prm_bingrp );

    /// @brief Legacy version of auto_divide_by_signal_noise_group (for regression testing).
    /// @details Uses full Grid2dBinGroup copy for trial splits. Retained for result comparison.
    std::tuple<bool,int> auto_divide_by_signal_noise_group_legacy(
      const Igroup igroup, const double signal_noise_group_trigger, const int nloop_max,
      const int ixlen_min, const int iylen_min, const bool tf_prefer_split_x );

    /// @brief   optimize the size of all bins with signal_group_trg, nloop_limit, ixlen_min, iylen_min.
    /// @details See the image below for an example of the final result \n
    /// \image html images/example_auto_grouping_by_signal_noise_group_all.gif "nx_div_init=4, ny_div_init=2, signal_noise_group_trig=500, ix/iylen_min=4, nloop_limit=10000" width=1000px
    /// @return the final number of groups.
    /// @param tf_prefer_split_x : when xlen==ylen, split along the x direction if true.
    int auto_divide_by_signal_noise_group_all(
      const double signal_noise_group_trigger, const int nloop_limit,
      const int ixlen_min, const int iylen_min, const bool tf_prefer_split_x );
    
    /// @brief   optimize the size of all bins with signal_group_trg, nloop_limit, ixlen_min, iylen_min.
    /// @details Internally calls @ref auto_divide_by_signal_noise_group_all(double,int,int,int).
    /// @return the final number of groups.
    int auto_divide_by_signal_noise_group_all(
      const Grid2dBinGroup::Parameters &prm_bingrp);
    
    // Optimize the size of all bins with signal_group_trg, nloop_limit, ixlen_min, iylen_min.
    // Returns the final number of groups.
    // ! ERROR
    // int mp_auto_grouping_by_signal_noise_group_all(
    //   const double signal_noise_group_trigger, const int nloop_limit,
    //   const int ixlen_min, const int iylen_min );
    
    // Optimize the size of all bins with signal_group_trg, nloop_limit, ixlen_min, iylen_min.
    // Returns the final number of groups.
    // ! ERROR
    // int mp_auto_grouping_by_signal_noise_group_all(
    //   const Grid2dBinGroup::Parameters &prm_bingrp ){
    //     return mp_auto_grouping_by_signal_noise_group_all(
    //       prm_bingrp.signal_noise_group_trig, prm_bingrp.nloop_limit,
    //       prm_bingrp.ixlen_min, prm_bingrp.iylen_min );
    // };

    /// @brief calculate signal_group / noise_poi_group and set them into vec_signal_group / vec_noise_poi_group.
    /// @note Uses OpenMP
    void calc_set_vec_signal_noise_group();

    /// @brief set vec_signal_group_poisson and vec_noise_poi_group_poisson from vec_signal_group and vec_noise_poi_group with a newly drawn Poisson error.
    /// @note Uses OpenMP
    void set_diff_poisson_again();

    /// @brief renumber the igroup of bimap_ starting from igroup_start.
    void renumbering_groups(const int igroup_start=0);
    ///@} ------------------------------------------------------------------

    //==================================================================
    /// @name check_functions
    ///@{
    
    /// @brief check with the area formula whether the bins of igroup form a rectangle. Return false if not a rectangle.
    bool is_rect( const Igroup igroup ) const;

    /// @brief check whether the signal_group of every igroup is greater than or equal to signal_thres.
    bool is_all_signal_group_larger_than_thres( const double signal_thres ) const;

    /// @brief check whether every group satisfies signal_noise_sum >= signal_noise_thres
    /// @throws std::runtime_error If a group does not satisfy the condition
    bool is_all_signal_noise_group_larger_than_thres(
      const double signal_noise_thres ) const;

    /// @brief check whether every group satisfies ixlen >= ixlen_min and iylen >= iylen_min.
    /// @details Violating groups are reported with LOG_ERROR and false is returned.
    ///          Does not throw, so it is safe to call inside an OpenMP parallel loop.
    bool is_all_group_ixiylen_larger_than_min(
      const int ixlen_min, const int iylen_min ) const;

    ///@} ------------------------------------------------------------------

    //==================================================================
    /// @name IO_functions
    ///@{
    
    /// @brief write the contents of bimap_ in ASCII to the file given by the fs::path argument. \n
    /// Rows are sorted in the order Detid, Igroup, Iy, Ix
    /// @param width Fixed width of each field (default 9).
    void write_bimap_all_to_ascii(
      const std::filesystem::path& file_path, const int width = 9) const;

    /// @brief  write the header information of bimap_ in ASCII format.
    void write_bimap_header_to_ascii(
      const std::filesystem::path& file_path, const int width=9) const;

    /// @brief save bimap_ to std::ofstream
    void save_bimap(std::ofstream& ofs) const;

    /// @brief load bimap_ from std::ifstream
    void load_bimap(std::ifstream& ifs);

    ///@brief save all member variables to std::ofstream
    void save( std::ofstream& ofs ) const;

    ///@brief load all member variables from std::ifstream
    void load( std::ifstream& ifs );

    ///@brief save all member variables to std::filesystem::path
    void save( const std::filesystem::path& pathout ) const;

    ///@brief load all member variables from std::filesystem::path
    void load( const std::filesystem::path &path_in );

    ///@} ------------------------------------------------------------------

};
