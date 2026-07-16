/// @file cls_Grid2dBinGroup.cpp
/// @brief Implementation of Grid2dBinGroup class
/// @details Extends Grid2d to support grouped bins for coarser resolution
///          analysis, mainly used for grouping DetectorElement bins to
///          increase muon signal statistics.
#include "cls_Grid2dBinGroup.hpp"
#include "cls_Grid2dBinGroupParameters.hpp"
#include "ns_myapp.hpp"
#include "ns_mymacro.hpp"
#include "spdlog_pch.hpp"
#include "ns_mylogger.hpp"

//###############################################
//###############################################
// class Grid2dBinGroup
// as derived class of the base class Grid2d
//###############################################
//###############################################

// inequality operator
bool Grid2dBinGroup::operator!=(const Grid2dBinGroup& other) const
{
  #ifdef NODEBUG
    if (Grid2d::operator!=(other)) return true;
    if (bimap_ != other.bimap_) return true;
    if (detid_ != other.detid_) return true;
    if (pdic_ != other.pdic_) return true;
    if (vec_signal_group != other.vec_signal_group) return true;
    if (vec_noise_poi_group != other.vec_noise_poi_group) return true;
    if (vec_noise_det_group != other.vec_noise_det_group) return true;
    if (vec_signal_group_poisson != other.vec_signal_group_poisson) return true;
    if (vec_noise_poi_group_poisson != other.vec_noise_poi_group_poisson) return true;
    if (vec_vec_is_avail != other.vec_vec_is_avail) return true;
    if (vec_vec_signal != other.vec_vec_signal) return true;
    if (vec_vec_noise_poi != other.vec_vec_noise_poi) return true;
    if (vec_vec_noise_det != other.vec_vec_noise_det) return true;
    if (done_grouping != other.done_grouping) return true;
    if (vec_dens_group_lower != other.vec_dens_group_lower) return true;
    if (vec_dens_group_center != other.vec_dens_group_center) return true;
    if (vec_dens_group_upper != other.vec_dens_group_upper) return true;
    if (vec_delta_nmuon_group_lower != other.vec_delta_nmuon_group_lower) return true;
    if (vec_delta_nmuon_group_center != other.vec_delta_nmuon_group_center) return true;
    if (vec_delta_nmuon_group_upper != other.vec_delta_nmuon_group_upper) return true;
    if (vec_volume_group != other.vec_volume_group) return true;
    if (vec_eff_low_group != other.vec_eff_low_group) return true;
    if (vec_eff_cnt_group != other.vec_eff_cnt_group) return true;
    if (vec_eff_upp_group != other.vec_eff_upp_group) return true;
  #else
    if (Grid2d::operator!=(other)) { LOG_WARN("Grid2dBinGroup: Base class Grid2d differs"); return true; }
    if (bimap_ != other.bimap_) { LOG_WARN("Grid2dBinGroup: bimap_ differs"); return true; }
    if (detid_ != other.detid_) { LOG_WARN("Grid2dBinGroup: detid_ differs"); return true; }
    if (pdic_ != other.pdic_) { LOG_WARN("Grid2dBinGroup: pdic_ differs"); return true; }
    if (vec_signal_group != other.vec_signal_group) { LOG_WARN("Grid2dBinGroup: vec_signal_group differs"); return true; }
    if (vec_noise_poi_group != other.vec_noise_poi_group) { LOG_WARN("Grid2dBinGroup: vec_noise_poi_group differs"); return true; }
    if (vec_noise_det_group != other.vec_noise_det_group) { LOG_WARN("Grid2dBinGroup: vec_noise_det_group differs"); return true; }
    if (vec_signal_group_poisson != other.vec_signal_group_poisson) { LOG_WARN("Grid2dBinGroup: vec_signal_group_poisson differs"); return true; }
    if (vec_noise_poi_group_poisson != other.vec_noise_poi_group_poisson) { LOG_WARN("Grid2dBinGroup: vec_noise_poi_group_poisson differs"); return true; }
    if (vec_vec_is_avail != other.vec_vec_is_avail) { LOG_WARN("Grid2dBinGroup: vec_vec_is_avail differs"); return true; }
    if (vec_vec_signal != other.vec_vec_signal) { LOG_WARN("Grid2dBinGroup: vec_vec_signal differs"); return true; }
    if (vec_vec_noise_poi != other.vec_vec_noise_poi) { LOG_WARN("Grid2dBinGroup: vec_vec_noise_poi differs"); return true; }
    if (vec_vec_noise_det != other.vec_vec_noise_det) { LOG_WARN("Grid2dBinGroup: vec_vec_noise_det differs"); return true; }
    if (done_grouping != other.done_grouping) { LOG_WARN("Grid2dBinGroup: done_grouping differs"); return true; }
    if (vec_dens_group_lower != other.vec_dens_group_lower) { LOG_WARN("Grid2dBinGroup: vec_dens_group_lower differs"); return true; }
    if (vec_dens_group_center != other.vec_dens_group_center) { LOG_WARN("Grid2dBinGroup: vec_dens_group_center differs"); return true; }
    if (vec_dens_group_upper != other.vec_dens_group_upper) { LOG_WARN("Grid2dBinGroup: vec_dens_group_upper differs"); return true; }
    if (vec_delta_nmuon_group_lower != other.vec_delta_nmuon_group_lower) { LOG_WARN("Grid2dBinGroup: vec_delta_nmuon_group_lower differs"); return true; }
    if (vec_delta_nmuon_group_center != other.vec_delta_nmuon_group_center) { LOG_WARN("Grid2dBinGroup: vec_delta_nmuon_group_center differs"); return true; }
    if (vec_delta_nmuon_group_upper != other.vec_delta_nmuon_group_upper) { LOG_WARN("Grid2dBinGroup: vec_delta_nmuon_group_upper differs"); return true; }
    if (vec_volume_group != other.vec_volume_group) { LOG_WARN("Grid2dBinGroup: vec_volume_group differs"); return true; }
    if (vec_eff_low_group != other.vec_eff_low_group) { LOG_WARN("Grid2dBinGroup: vec_eff_low_group differs"); return true; }
    if (vec_eff_cnt_group != other.vec_eff_cnt_group) { LOG_WARN("Grid2dBinGroup: vec_eff_cnt_group differs"); return true; }
    if (vec_eff_upp_group != other.vec_eff_upp_group) { LOG_WARN("Grid2dBinGroup: vec_eff_upp_group differs"); return true; }
  #endif
  return false;
}


// constructor from Grid1d x_axis & y_axis
Grid2dBinGroup::Grid2dBinGroup( const Grid1d &x_axis_in, const Grid1d &y_axis_in )
: Grid2dBinGroup()
{
  set_name("g2bg_" + x_axis_in.get_name() + "_" + y_axis_in.get_name());
  set_x_axis(x_axis_in);
  set_y_axis(y_axis_in);
  // memory allocation of vec_signal_group, vec_vec_is_avail, vec_vec_signal
  init_vec_vec();
}

// constructor from Grid2dXYZ
// assign a unique igroup to every bin
Grid2dBinGroup::Grid2dBinGroup(
  const Grid2dXYZ &g2xyz_in, const int igroup )
: Grid2dBinGroup()
{
  LOG_INFO(
    "from Grid2dXYZ {} ...\n"
    , g2xyz_in.get_name());

  // set name
  set_name( "g2bg_" + g2xyz_in.get_name() );
  
  // set x,y axis
  set_x_axis(g2xyz_in.get_x_axis());
  set_y_axis(g2xyz_in.get_y_axis());

  const auto log_level = spdlog::level::debug;

  get_x_axis().out_info(log_level);
  get_y_axis().out_info(log_level);
  
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  
  // memory allocation of vec_signal_group, vec_vec_is_avail, vec_vec_signal
  init_vec_vec();

  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      // set z for each bin
      set_signal(ix,iy,g2xyz_in.get_z(ix,iy));
      
      // set is_avail for each bin
      set_is_avail(ix,iy,false);

      // insert map
      insert_map(igroup,{ix,iy});
    }
  }
}

// constructor from Grid2dXYZ& g2xyz_signal_in
// divide the bins by nx_div, ny_div and assign nx_div*ny_div igroups
Grid2dBinGroup::Grid2dBinGroup(
  const Grid2dXYZ &g2xyz_signal_in, const int nx_div, const int ny_div )
: Grid2dBinGroup()
{
  LOG_INFO(
    "from Grid2dXYZ {} ...\n"
    , g2xyz_signal_in.get_name());

  // set name
  set_name( "g2bg_" + g2xyz_signal_in.get_name() );
  
  // set x,y axis
  set_x_axis(g2xyz_signal_in.get_x_axis());
  set_y_axis(g2xyz_signal_in.get_y_axis());
  get_x_axis().out_info(spdlog::level::debug);
  get_y_axis().out_info(spdlog::level::debug);
  
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  
  // memory allocation of vec_signal_group, vec_noise_poi_group, vec_vec_signal, , vec_vec_noise_poi, vec_vec_is_avail
  init_vec_vec();

  const int nx_per_group = std::ceil(static_cast<double>(nbinx) / nx_div);
  const int ny_per_group = std::ceil(static_cast<double>(nbiny) / ny_div);
  
  for (int iy = 0; iy < nbiny; ++iy) {
    const int igroup_y = iy / ny_per_group;
    for (int ix = 0; ix < nbinx; ++ix) {
      const int igroup_x = ix / nx_per_group;
      const int igroup = nx_div * igroup_y + igroup_x; // uniquely determine igroup here

      // set z for each bin
      set_signal(ix,iy,g2xyz_signal_in.get_z(ix,iy));

      // insert map
      insert_map(igroup,{ix,iy});
    }
  }
}

// memory allocation of vec_vec_signal, vec_vec_noise_poi, vec_vec_is_avail
void Grid2dBinGroup::init_vec_vec(
  const double signal_init, const double noise_init, const bool is_avail_init )
{
  // memory allocate of vec_vec_signal
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  vec_vec_signal.clear();
  vec_vec_signal.resize(nbiny);
  for( int iy = 0; iy < nbiny; iy++ ){
    vec_vec_signal.at(iy).resize(nbinx,signal_init);
  }

  // memory allocate of vec_vec_noise_poi (Poisson-bucket side, poi)
  vec_vec_noise_poi.clear();
  vec_vec_noise_poi.resize(nbiny);
  for( int iy = 0; iy < nbiny; iy++ ){
    vec_vec_noise_poi.at(iy).resize(nbinx,noise_init);
  }

  // memory allocate of vec_vec_noise_det (deterministic floor side, det)
  vec_vec_noise_det.clear();
  vec_vec_noise_det.resize(nbiny);
  for( int iy = 0; iy < nbiny; iy++ ){
    vec_vec_noise_det.at(iy).resize(nbinx,noise_init);
  }

  // memory allocate of vec_vec_is_avail
  vec_vec_is_avail.clear();
  vec_vec_is_avail.resize(nbiny);
  for( int iy = 0; iy < nbiny; iy++ ){
    vec_vec_is_avail.at(iy).resize(nbinx,is_avail_init);
  }
}

Grid2dBinGroup Grid2dBinGroup::cut(
  const double x_lower, const double x_upper
, const double y_lower, const double y_upper
, const double x_eps, const double y_eps) const
{
  // Get the new x / y axes with the base-class cut.
  // Its return value carries axis information only, so the BinGroup-specific data is copied manually below.
  const Grid2d base_cut = Grid2d::cut(x_lower, x_upper, y_lower, y_upper, x_eps, y_eps);

  // Get the extraction range for x and for y.
  const Grid1d::RangeIndices range_x = get_x_axis().calc_range_indices(x_lower, x_upper, x_eps);
  const Grid1d::RangeIndices range_y = get_y_axis().calc_range_indices(y_lower, y_upper, y_eps);

  // Create the result instance. Axis information, name and Detid are inherited from the source instance.
  Grid2dBinGroup result;
  result.set_name(get_name() + "_cut");
  result.set_x_axis(base_cut.get_x_axis());
  result.set_y_axis(base_cut.get_y_axis());
  result.set_x_axis_name(get_x_axis().get_name());
  result.set_y_axis_name(get_y_axis().get_name());
  result.set_detid(detid_);
  
  // Allocate the BinGroup-specific 2D buffers at the new size (all cleared to their initial values).
  result.init_vec_vec(signal_init_vecvec, noise_init_vecvec, is_avail_vecvec_init);

  const int nbinx_org = get_nbinx();
  const int nbiny_org = get_nbiny();
  const int nbinx_new = result.get_nbinx();
  const int nbiny_new = result.get_nbiny();

  // Check that the expected new size (the range length) matches the array size actually allocated.
  const int expected_nbinx = (range_x.end >= range_x.start) ? (range_x.end - range_x.start + 1) : 0;
  const int expected_nbiny = (range_y.end >= range_y.start) ? (range_y.end - range_y.start + 1) : 0;
  if (expected_nbinx != nbinx_new || expected_nbiny != nbiny_new) {
    THROW_ERROR("Grid2dBinGroup::cut : inconsistent bin counts between source and result axes");
  }

  // Using the index correspondence obtained from range_x / range_y,
  // copy the signal / noise / is_avail data into the new grid.
  // The grouping information such as bimap_ is not rebuilt here; only the physical values are transferred.
  for (int iy_new = 0; iy_new < nbiny_new; ++iy_new) {
    const int iy_org = range_y.start + iy_new;
    for (int ix_new = 0; ix_new < nbinx_new; ++ix_new) {
      const int ix_org = range_x.start + ix_new;
      result.set_signal(ix_new, iy_new, get_signal(ix_org, iy_org));
      result.set_noise(ix_new, iy_new, get_noise(ix_org, iy_org));
      result.set_noise_det(ix_new, iy_new, get_noise_det(ix_org, iy_org));
      result.set_is_avail(ix_new, iy_new, is_avail(ix_org, iy_org));
    }
  }

  // Initialize the per-group aggregation buffers as empty, matching the new structure.
  // Recalculation is required after the cut, so all of them are refilled with zero length and initial values.
  result.resize_vec_is_avail_group(0, is_avail_group_init);
  result.resize_vec_signal_group(0, signal_group_init);
  result.resize_vec_noise_poi_group(0, noise_group_init);
  result.resize_vec_noise_det_group(0, noise_group_init);
  result.resize_vec_signal_group_poisson(0, signal_group_init);
  result.resize_vec_noise_poi_group_poisson(0, noise_group_init);
  result.resize_vec_dens_group_lower(0, dens_group_lower_init);
  result.resize_vec_dens_group_center(0, dens_group_center_init);
  result.resize_vec_dens_group_upper(0, dens_group_upper_init);
  result.resize_vec_delta_nmuon_group_lower(0, delta_nmuon_group_lower_init);
  result.resize_vec_delta_nmuon_group_center(0, delta_nmuon_group_center_init);
  result.resize_vec_delta_nmuon_group_upper(0, delta_nmuon_group_upper_init);
  result.resize_vec_volume_group(0, volume_group_init);
  result.resize_vec_eff_low_group(0, eff_low_group_init);
  result.resize_vec_eff_cnt_group(0, eff_cnt_group_init);
  result.resize_vec_eff_upp_group(0, eff_upp_group_init);
  // bimap_ is left discarded, so done_grouping is reset to the not-yet-done flag.
  result.set_done_grouping(false);

  return result;
}

// allocate memory for Grid2dBinGroup::vec_value_group
void Grid2dBinGroup::allocate_vec_value_group()
{
  // check_1st_grouping
  check_done_grouping();

  // initialize vec_xxxxxxx_group
  const int ngroup = get_n_group();
  resize_vec_is_avail_group(ngroup, is_avail_group_init);
  resize_vec_signal_group(ngroup, signal_group_init);
  resize_vec_noise_poi_group(ngroup, noise_group_init);
  resize_vec_noise_det_group(ngroup, noise_group_init);
  resize_vec_signal_group_poisson(ngroup, signal_group_init);
  resize_vec_noise_poi_group_poisson(ngroup, noise_group_init);
  resize_vec_dens_group_lower(ngroup, dens_group_lower_init);
  resize_vec_dens_group_center(ngroup, dens_group_center_init);
  resize_vec_dens_group_upper(ngroup, dens_group_upper_init);
  resize_vec_delta_nmuon_group_lower(ngroup, delta_nmuon_group_lower_init);
  resize_vec_delta_nmuon_group_center(ngroup, delta_nmuon_group_center_init);
  resize_vec_delta_nmuon_group_upper(ngroup, delta_nmuon_group_upper_init);
  resize_vec_volume_group(ngroup, volume_group_init);
  resize_vec_eff_low_group(ngroup, eff_low_group_init);
  resize_vec_eff_cnt_group(ngroup, eff_cnt_group_init);
  resize_vec_eff_upp_group(ngroup, eff_upp_group_init);
}

//it calls read_bin_list2
Grid2dBinGroup::Grid2dBinGroup( const fs::path &path_in )
: Grid2dBinGroup()
{
  read_bin_list2(path_in);
}

// get the vec_dens_group_lower of igroup
double Grid2dBinGroup::get_dens_group_lower(
  const Igroup igroup ) const
{
  return get_group_value(
    vec_dens_group_lower, igroup, dens_group_lower_init, "vec_dens_group_lower");
}

// get the vec_dens_group_center of igroup
double Grid2dBinGroup::get_dens_group_center( const Igroup igroup ) const
{
  return get_group_value(
    vec_dens_group_center, igroup, dens_group_center_init, "vec_dens_group_center");
}

// get the vec_dens_group_upper of igroup
double Grid2dBinGroup::get_dens_group_upper( const Igroup igroup ) const
{
  return get_group_value(
    vec_dens_group_upper, igroup, dens_group_upper_init, "vec_dens_group_upper");
}

// get the vec_delta_nmuon_group_lower of igroup
double Grid2dBinGroup::get_delta_nmuon_group_lower( const Igroup igroup ) const
{
  return get_group_value(
    vec_delta_nmuon_group_lower, igroup, delta_nmuon_group_lower_init
    , "vec_delta_nmuon_group_lower");
}

// get the vec_delta_nmuon_group_center of igroup
double Grid2dBinGroup::get_delta_nmuon_group_center( const Igroup igroup ) const
{
  return get_group_value(
    vec_delta_nmuon_group_center, igroup, delta_nmuon_group_center_init
    , "vec_delta_nmuon_group_center");
}

// get the vec_delta_nmuon_group_upper of igroup
double Grid2dBinGroup::get_delta_nmuon_group_upper( const Igroup igroup ) const
{
  return get_group_value(
    vec_delta_nmuon_group_upper, igroup, delta_nmuon_group_upper_init
    , "vec_delta_nmuon_group_upper");
}

// get the vec_volume_group of igroup
double Grid2dBinGroup::get_volume_group( const Igroup igroup ) const
{
  return get_group_value(
    vec_volume_group, igroup, volume_group_init, "vec_volume_group");
}

// get the vec_eff_low_group of igroup
double Grid2dBinGroup::get_eff_low_group( const Igroup igroup ) const
{
  return get_group_value(
    vec_eff_low_group, igroup, eff_low_group_init, "vec_eff_low_group");
}

// get the vec_eff_cnt_group of igroup
double Grid2dBinGroup::get_eff_cnt_group( const Igroup igroup ) const
{
  return get_group_value(
    vec_eff_cnt_group, igroup, eff_cnt_group_init, "vec_eff_cnt_group");
}

// get the vec_eff_upp_group of igroup
double Grid2dBinGroup::get_eff_upp_group( const Igroup igroup ) const
{
  return get_group_value(
    vec_eff_upp_group, igroup, eff_upp_group_init, "vec_eff_upp_group");
}

// get the subtracted vec_vec_signal , this - other
std::vector<std::vector<double>>
Grid2dBinGroup::get_vec_vec_signal_subtracted( const Grid2dBinGroup &other ) const
{
  LOG_INFO("...\n");

  // check size
  if( get_nbinx() != other.get_nbinx() || get_nbiny() != other.get_nbiny() ){
    THROW_ERROR_NAME("get_vec_vec_signal_subtracted: size mismatch");
  }

  // allocate result vector
  std::vector<std::vector<double>> 
    result(get_nbiny(), std::vector<double>(get_nbinx(), 0.0));

  // subtract
  #pragma omp parallel for collapse(2)
  for(int iy = 0; iy < get_nbiny(); iy++){
    for(int ix = 0; ix < get_nbinx(); ix++){
      result.at(iy).at(ix) = this->get_signal(ix, iy) - other.get_signal(ix, iy);
    }
  }

  return result;
}

// get the subtracted vec_vec_noise_poi , this - other
std::vector<std::vector<double>>
Grid2dBinGroup::get_vec_vec_noise_poi_subtracted( const Grid2dBinGroup &other ) const
{
  LOG_INFO("...\n");

  // check size
  if( get_nbinx() != other.get_nbinx() || get_nbiny() != other.get_nbiny() ){
    THROW_ERROR_NAME("get_vec_vec_noise_poi_subtracted: size mismatch");
  }

  // allocate result vector
  std::vector<std::vector<double>> result(get_nbiny(), std::vector<double>(get_nbinx(), 0.0));

  // subtract
  #pragma omp parallel for collapse(2)
  for(int iy = 0; iy < get_nbiny(); iy++){
    for(int ix = 0; ix < get_nbinx(); ix++){
      result.at(iy).at(ix) = this->get_noise(ix, iy) - other.get_noise(ix, iy);
    }
  }

  return result;
}


// get the tuple of subtracted vec_vec_signal, vec_vec_noise_poi
std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>>
Grid2dBinGroup::get_vec_vec_signal_noise_subtracted( const Grid2dBinGroup &other ) const
{
  LOG_INFO("...\n");

  // check size
  if( get_nbinx() != other.get_nbinx() || get_nbiny() != other.get_nbiny() ){
    THROW_ERROR_NAME("get_vec_vec_signal_noise_subtracted: size mismatch");
  }

  // allocate result vectors
  std::vector<std::vector<double>>
    result_signal(get_nbiny(), std::vector<double>(get_nbinx(), 0.0));
  
  std::vector<std::vector<double>>
    result_noise(get_nbiny(), std::vector<double>(get_nbinx(), 0.0));

  // subtract
  #pragma omp parallel for collapse(2)
  for(int iy = 0; iy < get_nbiny(); iy++){
    for(int ix = 0; ix < get_nbinx(); ix++){
      result_signal.at(iy).at(ix) = get_signal(ix, iy) - other.get_signal(ix, iy);
      result_noise.at(iy).at(ix) = get_noise(ix, iy) - other.get_noise(ix, iy);
    }
  }

  return {result_signal, result_noise};
}


// get the subtracted vec_signal_group , this - other
std::vector<double>
Grid2dBinGroup::get_vec_signal_group_subtracted( const Grid2dBinGroup &other ) const
{
  LOG_INFO("...\n");

  // check size
  if( get_nbinx() != other.get_nbinx() || get_nbiny() != other.get_nbiny() ){
    THROW_ERROR_NAME("get_vec_signal_group_subtracted: size mismatch");
  }

  const std::vector<Igroup> vec_igroup = get_vec_igroup();
  const std::vector<Igroup> other_vec_igroup = other.get_vec_igroup();
  if( vec_igroup.size() != other_vec_igroup.size() ){
    THROW_ERROR_NAME("get_vec_signal_group_subtracted: vec_igroup size mismatch");
  }
  // allocate result vector
  std::vector<double> result(vec_igroup.size(), 0.0);

  // subtract
  #pragma omp parallel for
  for(size_t i = 0; i < vec_igroup.size(); i++) {
    if (vec_igroup.at(i) != other_vec_igroup.at(i)) {
      THROW_ERROR_NAME("get_vec_signal_group_subtracted: vec_igroup mismatch at index " + std::to_string(i));
    }
    result.at(i) = get_signal_group(vec_igroup.at(i)) - other.get_signal_group(other_vec_igroup.at(i));
  }
  return result;
}

// get the subtracted vec_noise_poi_group , this - other
std::vector<double>
Grid2dBinGroup::get_vec_noise_poi_group_subtracted( const Grid2dBinGroup &other ) const
{
  LOG_INFO("...\n");

  // check size
  if( get_nbinx() != other.get_nbinx() || get_nbiny() != other.get_nbiny() ){
    THROW_ERROR_NAME("get_vec_noise_poi_group_subtracted: size mismatch");
  }

  const std::vector<Igroup> vec_igroup = get_vec_igroup();
  const std::vector<Igroup> other_vec_igroup = other.get_vec_igroup();
  if( vec_igroup.size() != other_vec_igroup.size() ){
    THROW_ERROR_NAME("get_vec_noise_poi_group_subtracted: vec_igroup size mismatch");
  }
  
  // allocate result vector
  std::vector<double> result(vec_igroup.size(), 0.0);

  // subtract
  #pragma omp parallel for
  for(size_t i = 0; i < vec_igroup.size(); i++) {
    if (vec_igroup.at(i) != other_vec_igroup.at(i)) {
      THROW_ERROR_NAME("get_vec_noise_poi_group_subtracted: vec_igroup mismatch at index " + std::to_string(i));
    }
    result.at(i) = get_noise_poi_group(vec_igroup.at(i)) - other.get_noise_poi_group(other_vec_igroup.at(i));
  }

  return result;
}

// get the tuple of subtracted vec_signal_group, vec_noise_poi_group
std::tuple<std::vector<double>, std::vector<double>>
Grid2dBinGroup::get_vec_signal_noise_group_subtracted( const Grid2dBinGroup &other ) const
{
  LOG_INFO("...\n");

  // check size
  if( get_nbinx() != other.get_nbinx() || get_nbiny() != other.get_nbiny() ){
    THROW_ERROR_NAME("get_vec_signal_noise_group_subtracted: size mismatch");
  }

  const std::vector<Igroup> vec_igroup = get_vec_igroup();
  const std::vector<Igroup> other_vec_igroup = other.get_vec_igroup();
  if( vec_igroup.size() != other_vec_igroup.size() ){
    LOG_ERROR("vec_igroup.size() = {}, other_vec_igroup.size() = {}",
              vec_igroup.size(), other_vec_igroup.size());
    LOG_ERROR("get_vec_signal_noise_group_subtracted: vec_igroup size mismatch");
    THROW_ERROR_NAME("get_vec_signal_noise_group_subtracted: vec_igroup size mismatch");
  }
  
  // allocate result vectors
  std::vector<double> result_signal(vec_igroup.size(), 0.0);
  std::vector<double> result_noise(vec_igroup.size(), 0.0);

  // subtract
  #pragma omp parallel for
  for(size_t i = 0; i < vec_igroup.size(); i++) {
    if (vec_igroup.at(i) != other_vec_igroup.at(i)) {
      THROW_ERROR_NAME("get_vec_signal_noise_group_subtracted: vec_igroup mismatch at index " + std::to_string(i));
    }
    result_signal.at(i) = get_signal_group(vec_igroup.at(i)) - other.get_signal_group(other_vec_igroup.at(i));
    result_noise.at(i) = get_noise_poi_group(vec_igroup.at(i)) - other.get_noise_poi_group(other_vec_igroup.at(i));
  }

  return {result_signal, result_noise};
}

// read_bin_list2
void Grid2dBinGroup::read_bin_list2( const fs::path &path_in )
{
  LOG_INFO("reading file {} ...\n",path_in.string());

  // set name
  set_name( path_in.filename().string() );

  // check file existence
  if(!std::filesystem::exists(path_in)) THROW_ERROR_NAME( path_in.string() + " does not exist");

  // file stream open
  std::ifstream reading_file;
  reading_file.open(path_in,std::ios::in);

  // define variables for reading
  std::string reading_line_str;
  std::vector<std::string> vec_str;

  // line 1: nbinx0 xmin0 xmax0 xpit0
  std::getline(reading_file, reading_line_str);
  vec_str = myapp::split(reading_line_str,myapp::char_delim_default);
  if( vec_str.size()!=4 ) THROW_ERROR_NAME2("vec_str.size()!=4",vec_str.size());
  const int   nbinx0 = std::stoi( vec_str.at(0) );
  const double xmin0 = std::stod( vec_str.at(1) );
  const double xmax0 = std::stod( vec_str.at(2) );
  const double xpit0 = std::stod( vec_str.at(3) );
  set_x_axis("xaxis_"+get_name(),nbinx0,xmin0,xmax0,xpit0);

  // line 2: nbiny0 ymin0 ymax0 ypit0
  std::getline(reading_file, reading_line_str);
  vec_str = myapp::split(reading_line_str,myapp::char_delim_default);
  if( vec_str.size()!=4 ) THROW_ERROR_NAME2("vec_str.size()!=4",vec_str.size());
  const int   nbiny0 = std::stoi( vec_str.at(0) );
  const double ymin0 = std::stod( vec_str.at(1) );
  const double ymax0 = std::stod( vec_str.at(2) );
  const double ypit0 = std::stod( vec_str.at(3) );
  set_y_axis("yaxis_"+get_name(),nbiny0,ymin0,ymax0,ypit0);

  // line 3 onward: xmin, xmax, ymin, ymax, nmuon
  int nrow=0;
  double xmin, xmax, ymin, ymax;
  double nmuon;
  vec_signal_group.clear();
  Igroup igroup=0;
  while ( std::getline(reading_file, reading_line_str) ){
    nrow++;
    std::vector<std::string> vec_str = myapp::split(reading_line_str,myapp::char_delim_default);
    if( vec_str.size()!=5 ) THROW_ERROR_NAME2("vec_str.size()!=5",vec_str.size());
    xmin = std::stod( vec_str.at(0) );
    xmax = std::stod( vec_str.at(1) );
    ymin = std::stod( vec_str.at(2) );
    ymax = std::stod( vec_str.at(3) );
    nmuon = std::stof( vec_str.at(4) );
    std::vector<Ixiy> vec_ixiy = Grid2d::get_vec_ixiy(xmin,xmax,ymin,ymax);
    insert_map(igroup,vec_ixiy);
    push_back_nmuon(nmuon);
    igroup++;
  }
  LOG_INFO("read {} rows\n", nrow);
}

// get_xmin_xmax_ymin_ymax in the group_index=igroup_in
std::array<double,4> Grid2dBinGroup::get_xmin_xmax_ymin_ymax(const int igroup_in) const
{
  // to be returned
  constexpr double big = std::numeric_limits<double>::max();
  double xmin_result =  big;
  double xmax_result = -big;
  double ymin_result =  big;
  double ymax_result = -big;

  // get the number of (ix,iy) in the bin group_index=group_id
  const size_t n_pair = get_n_uqid(igroup_in);

  // loop for the index of pairs in the bin group_index=group_id
  const std::vector<Ixiy> vec_ixiy = get_vec_ixiy(igroup_in);
  if( n_pair==0 ) return {0,0,0,0};

  // n_pair > 0
  for(const auto [ix,iy] : vec_ixiy ){
    const auto [xmin,xmax,ymin,ymax] = Grid2d::get_xy_lower_upper(ix,iy);
    xmin_result = std::min(xmin_result,xmin);
    xmax_result = std::max(xmax_result,xmax);
    ymin_result = std::min(ymin_result,ymin);
    ymax_result = std::max(ymax_result,ymax);
  }
  return {xmin_result,xmax_result,ymin_result,ymax_result};
}

// output all igroup, ix, iy
void Grid2dBinGroup::out_igroup_ixiy_all(FILE* fout) const
{
  for (const auto& igroup : bimap_.get_setOne()) {
    for (const auto& ixiy : bimap_.get_vecMany(igroup)) {
      fprintf(fout, "%d %d %d\n", igroup, ixiy[0], ixiy[1]);
    }
  }
}

// output all igroup, ix, iy to a text file
void Grid2dBinGroup::out_igroup_ixiy_all(const fs::path& pathout) const
{
  if (std::filesystem::exists(pathout)) {
    LOG_WARN("out_igroup_ixiy_all: file already exists: {}", pathout.string());
  }
  FILE* fout = std::fopen(pathout.string().c_str(), "w");
  if (!fout) {
    THROW_ERROR("Grid2dBinGroup::out_igroup_ixiy_all: Failed to open file. path={}", pathout.string());
  }

  out_igroup_ixiy_all(fout);
  std::fclose(fout);
}



// get the vector of the number of muons in all the bin group
Eigen::VectorXf Grid2dBinGroup::get_vecxf_nmuon( ) const
{
  const int num_groups = get_n_group();
  Eigen::VectorXf vecxf_nmuon(num_groups);
  for(int i=0;i<num_groups;i++) vecxf_nmuon(i) = get_signal_group(i);
  if( vecxf_nmuon.allFinite()==false ){
    THROW_ERROR_NAME("vecxf_nmuon.allFinite()==false");
  }
  return vecxf_nmuon;
}

std::set<Igroup> Grid2dBinGroup::get_set_igroup() const
{
  std::set<Igroup> set_igroup;
  // get the set of igroup from bimap_
  const auto& map = bimap_.getMapOneMany(); // multimap<One,Many>
  bool first = true;
  Igroup last = IgroupNotAssigned; // default initial value of the type
  for (const auto& [igroup, _] : map) {
    if (first || igroup != last) {
      set_igroup.insert(igroup);
      last = igroup;
      first = false;
    }
  }
  return set_igroup;
}

std::vector<Ixiy> Grid2dBinGroup::get_vec_ixiy(const Igroup igroup) const
{
  std::vector<Ixiy> vec_ixiy = bimap_.get_vecMany(igroup);
  if (vec_ixiy.empty())
    LOG_ERROR("igroup={} has no ixiy, in detid={}"
    , igroup, get_detid());
  return vec_ixiy;
}


// get vector of ixiy for all igroup
std::vector< std::vector< Ixiy > >  
  Grid2dBinGroup::get_vec_igroup_vec_ixiy() const
{
  std::vector< std::vector< Ixiy > > vec_igroup_vec_ixiy;
  const int num_groups = get_n_group();
  for(int i_group=0;i_group<num_groups;i_group++){
    vec_igroup_vec_ixiy.push_back( get_vec_ixiy(i_group) );
  }
  return vec_igroup_vec_ixiy;
}

// return vec_vec_ixiy as std::vector<std::vector<Ixiy>>
std::vector<std::vector<Ixiy>> 
  Grid2dBinGroup::get_vec_vec_ixiy( const Igroup igroup ) const
{
  std::vector<std::vector<Ixiy>> vec_vec_ixiy;
  const auto vec_ixiy = get_vec_ixiy(igroup);
  for( const auto ixiy : vec_ixiy ){
    const auto [ix,iy] = ixiy;
    vec_vec_ixiy.at(iy).at(ix) = ixiy;
  }
  return vec_vec_ixiy;
}


// get igroup from ixiy
Igroup Grid2dBinGroup::get_igroup(const int ix, const int iy) const {
  Igroup ig = bimap_.getOne({ix, iy});
  if (ig == OneToManyBimap<Igroup, Ixiy>::OneNotAssigned){
    THROW_ERROR("Grid2dBinGroup::get_igroup: Ixiy not registered. ix={}, iy={}", ix, iy);
  }
  return ig;
}


// @brief get the set of ix,iy from igroup
// @param igroup : group index
std::set<Ixiy> Grid2dBinGroup::get_set_ixiy(const Igroup igroup) const
{
  std::vector<Ixiy> vec_ixiy = get_vec_ixiy(igroup);
  return std::set<Ixiy>(vec_ixiy.begin(), vec_ixiy.end());
}

//==============================================================
// 2023-09-21 15:28:06
//==============================================================
void Grid2dBinGroup::set_signal( const int ix, const int iy, const double signal_in )
{
  // if out of range. throw error
  check_ix_inside(ix);
  check_iy_inside(iy);
  vec_vec_signal.at(iy).at(ix) = signal_in;
}

void Grid2dBinGroup::set_signal( const Ixiy &ixiy, const double signal_in )
{
  const auto [ix,iy] = ixiy;
  set_signal(ix,iy,signal_in);
}

void Grid2dBinGroup::set_signal( const double x, const double y, const double signal_in )
{
  // if out of range. throw error
  const int ix = get_ix(x);
  const int iy = get_iy(y);
  set_signal(ix,iy,signal_in);
}

void Grid2dBinGroup::set_noise( const int ix, const int iy, const double noise_in )
{
  // if out of range. throw error
  check_ix_inside(ix);
  check_iy_inside(iy);
  vec_vec_noise_poi.at(iy).at(ix) = noise_in;
}

void Grid2dBinGroup::set_noise( const Ixiy &ixiy, const double noise_in )
{
  const auto [ix,iy] = ixiy;
  set_noise(ix,iy,noise_in);
}

void Grid2dBinGroup::set_noise( const double x, const double y, const double noise_in )
{
  // if out of range. throw error
  const int ix = get_ix(x);
  const int iy = get_iy(y);
  set_noise(ix,iy,noise_in);
}

void Grid2dBinGroup::set_noise_det( const int ix, const int iy, const double noise_det_in )
{
  // if out of range. throw error
  check_ix_inside(ix);
  check_iy_inside(iy);
  vec_vec_noise_det.at(iy).at(ix) = noise_det_in;
}

double Grid2dBinGroup::get_signal( const int ix, const int iy ) const
{
  // if out of range. throw error
  const int nsize = vec_vec_signal.size();
  if( nsize <= 0 ) THROW_ERROR_NAME("vec_vec_signal.size() <= 0");
  check_ix_inside(ix);
  check_iy_inside(iy);
  return vec_vec_signal.at(iy).at(ix);
}

double Grid2dBinGroup::get_signal( const Ixiy& ixiy ) const
{
  const auto [ix,iy] = ixiy;
  return get_signal(ix,iy);
}

double Grid2dBinGroup::get_signal( const double x, const double y ) const
{
  // if out of range. throw error
  const int ix = get_ix(x);
  const int iy = get_iy(y);
  return vec_vec_signal.at(iy).at(ix);
}

double Grid2dBinGroup::get_noise( const int ix, const int iy ) const
{
  // if out of range. throw error
  const int nsize = vec_vec_noise_poi.size();
  if( nsize <= 0 ) THROW_ERROR_NAME("vec_vec_noise_poi.size() <= 0");
  check_ix_inside(ix);
  check_iy_inside(iy);
  return vec_vec_noise_poi.at(iy).at(ix);
}

double Grid2dBinGroup::get_noise( const Ixiy& ixiy ) const
{
  const auto [ix,iy] = ixiy;
  return get_noise(ix,iy);
}

double Grid2dBinGroup::get_noise( const double x, const double y ) const
{
  // if out of range. throw error
  const int ix = get_ix(x);
  const int iy = get_iy(y);
  return vec_vec_noise_poi.at(iy).at(ix);
}

double Grid2dBinGroup::get_noise_det( const int ix, const int iy ) const
{
  // if out of range. throw error
  const int nsize = vec_vec_noise_det.size();
  if( nsize <= 0 ) THROW_ERROR_NAME("vec_vec_noise_det.size() <= 0");
  check_ix_inside(ix);
  check_iy_inside(iy);
  return vec_vec_noise_det.at(iy).at(ix);
}

void Grid2dBinGroup::set_is_avail( const int ix, const int iy, const bool is_avail_in )
{
  // if out of range. throw error
  check_ix_inside(ix);
  check_iy_inside(iy);
  vec_vec_is_avail.at(iy).at(ix) = is_avail_in;
}

void Grid2dBinGroup::set_is_avail( const Ixiy &ixiy, const bool is_avail_in )
{
  const auto [ix,iy] = ixiy;
  set_is_avail(ix,iy,is_avail_in);
}

void Grid2dBinGroup::set_is_avail(const double x, const double y, const bool is_avail_in )
{
  // if out of range. throw error
  const int ix = get_ix(x);
  const int iy = get_iy(y);
  set_is_avail(ix,iy,is_avail_in);
}

void Grid2dBinGroup::set_is_avail( const Igroup igroup, const bool is_avail_in )
{
  const auto vec_ixiy = get_vec_ixiy(igroup);
  for( const auto &ixiy : vec_ixiy ){
    set_is_avail(ixiy,is_avail_in);
  }
}

// reassign the indices of the existing igroup
void Grid2dBinGroup::reassign_igroup_indices( const Igroup igroup_start )
{
  // get the existing igroup and sort them
  auto set_igroups = get_set_igroup();
  std::vector<Igroup> sorted_igroups(set_igroups.begin(), set_igroups.end());

  // build the old-to-new mapping
  std::map<Igroup, Igroup> map_oldToNew;
  for (size_t i = 0; i < sorted_igroups.size(); ++i) {
    map_oldToNew[sorted_igroups[i]] = igroup_start + static_cast<Igroup>(i);
  }

  // copy all existing pairs
  auto all_pairs = bimap_.getMapOneMany();

  // clear and register again
  bimap_.clear();
  for (const auto& kv : all_pairs) {
    try {
      bimap_.insert(map_oldToNew[kv.first], kv.second);
    }
    catch (...) {
      THROW_ERROR("reassign_igroup_indices: insertion failed");
    }
  }

  // consistency check
  if (!bimap_.isConsistent()) {
    THROW_ERROR("reassign_igroup_indices: consistency check failed");
  }
}


bool Grid2dBinGroup::is_avail( const int ix, const int iy ) const
{
  // if out of range. throw error
  check_ix_inside(ix);
  check_iy_inside(iy);
  return vec_vec_is_avail.at(iy).at(ix);
}

bool Grid2dBinGroup::is_avail( const Ixiy &ixiy ) const
{
  const auto [ix,iy] = ixiy;
  return is_avail(ix,iy);
}

bool Grid2dBinGroup::is_avail( const double x, const double y ) const
{
  // if out of range. throw error
  const int ix = get_ix(x);
  const int iy = get_iy(y);
  return vec_vec_is_avail.at(iy).at(ix);
}

// output the ix, iy whose is_avail is false to the logger
void Grid2dBinGroup::log_non_avail_vec_ixiy( const spdlog::level::level_enum &level, const Igroup igroup ) const
{
  mylogger::g_logger->log(level,"Grid2dBinGroup::log_non_avail_vec_ixiy igroup={} is not available",igroup);
const auto [xmin,xmax,ymin,ymax] = get_xmin_xmax_ymin_ymax(igroup);
  mylogger::g_logger->log(level,"xmin,xmax,ymin,ymax={:7.4f},{:7.4f},{:7.4f},{:7.4f}",xmin,xmax,ymin,ymax);
  const auto vec_ixiy = get_vec_ixiy(igroup);
  for( const auto [ix,iy] : vec_ixiy ){
    mylogger::g_logger->log(level,"{} {}, signal(ix,iy)={:E}, is_avail(ix,iy)={}"
    , ix, iy, get_signal(ix,iy), is_avail(ix,iy));
  }
}

// return false if at least one bin is not available
bool Grid2dBinGroup::is_avail_group( const Igroup igroup ) const
{
  const auto vec_ixiy = get_vec_ixiy(igroup);
  for( const auto &ixiy : vec_ixiy ){
    if( is_avail(ixiy) == false ){
      return false;
    }
  }
  return true;
}

// set all is_avail to is_avail_in
void Grid2dBinGroup::clear_is_avail(const bool is_avail_in)
{
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  for( int iy = 0; iy < nbiny; iy++ ){
    for( int ix = 0; ix < nbinx; ix++ ){
      set_is_avail(ix,iy,is_avail_in);
    }
  }
}

// calculate the signal_group of igroup.
double Grid2dBinGroup::calc_signal_group( const Igroup igroup ) const
{
  std::vector< Ixiy > vec_ixiy = get_vec_ixiy(igroup);
  double signal_group = 0;
  for(const auto &ixiy : vec_ixiy){
    signal_group += get_signal(ixiy);
  }
  return signal_group;
}

// calculate the noise_poi_group of igroup (Poisson-bucket side, poi).
double Grid2dBinGroup::calc_noise_poi_group( const Igroup igroup ) const
{
  std::vector< Ixiy > vec_ixiy = get_vec_ixiy(igroup);
  double noise_poi_group = 0;
  for(const auto &ixiy : vec_ixiy){
    noise_poi_group += get_noise(ixiy);
  }
  return noise_poi_group;
}

// calculate the deterministic floor noise_det_group (det) of igroup.
double Grid2dBinGroup::calc_noise_det_group( const Igroup igroup ) const
{
  std::vector< Ixiy > vec_ixiy = get_vec_ixiy(igroup);
  double noise_det_group = 0;
  for(const auto &ixiy : vec_ixiy){
    const auto [ix,iy] = ixiy;
    noise_det_group += get_noise_det(ix,iy);
  }
  return noise_det_group;
}


// Sum of (signal + noise) over a rectangular bin region [ixmin..ixmax] x [iymin..iymax].
// Direct access to vec_vec_signal/vec_vec_noise_poi; no bimap_ lookup.
double Grid2dBinGroup::calc_signal_noise_rect(
  int ixmin, int ixmax, int iymin, int iymax) const
{
  double sum = 0.0;
  for (int iy = iymin; iy <= iymax; ++iy)
    for (int ix = ixmin; ix <= ixmax; ++ix)
      sum += vec_vec_signal[iy][ix] + vec_vec_noise_poi[iy][ix];
  return sum;
}

// Determine which half has larger signal+noise after a trial 2-split.
// Reproduces split2DArrayIntoGrid(arr, 2, 1, 1, 0) for split_x=true,
//                                 (arr, 1, 2, 0, 1) for split_x=false.
Grid2dBinGroup::TrialSplitResult Grid2dBinGroup::calc_trial_split_max_half(
  int ixmin, int ixmax, int iymin, int iymax, bool split_x) const
{
  if (split_x) {
    // X-split: nx_div=2, ny_div=1, ix_extra=1, iy_extra=0
    // targetGroup_x = ix_extra = 1 (0-based: group 0 gets the extra columns)
    const int ixlen = ixmax - ixmin + 1;
    const int baseCol = ixlen / 2;
    const int extraCol = ixlen % 2;
    // half 0: [ixmin, ixmin + baseCol + extraCol - 1]  (gets remainder)
    // half 1: [ixmin + baseCol + extraCol, ixmax]
    const int ix_boundary = ixmin + baseCol + extraCol;
    if (ix_boundary > ixmax) {
      // only 1 valid half (region too small to split)
      return {0, 1};
    }
    const double sn0 = calc_signal_noise_rect(ixmin, ix_boundary - 1, iymin, iymax);
    const double sn1 = calc_signal_noise_rect(ix_boundary, ixmax, iymin, iymax);
    return (sn0 >= sn1) ? TrialSplitResult{0, 2} : TrialSplitResult{1, 2};
  } else {
    // Y-split: nx_div=1, ny_div=2, ix_extra=0, iy_extra=1
    // targetGroup_y = iy_extra = 1 (0-based: group 0 gets the extra rows)
    const int iylen = iymax - iymin + 1;
    const int baseRow = iylen / 2;
    const int extraRow = iylen % 2;
    // half 0: [iymin, iymin + baseRow + extraRow - 1]  (gets remainder)
    // half 1: [iymin + baseRow + extraRow, iymax]
    const int iy_boundary = iymin + baseRow + extraRow;
    if (iy_boundary > iymax) {
      return {0, 1};
    }
    const double sn0 = calc_signal_noise_rect(ixmin, ixmax, iymin, iy_boundary - 1);
    const double sn1 = calc_signal_noise_rect(ixmin, ixmax, iy_boundary, iymax);
    return (sn0 >= sn1) ? TrialSplitResult{0, 2} : TrialSplitResult{1, 2};
  }
}

// return the igroup with the smallest signal_group + noise_poi_group
// throws if vec_igroup is empty
Igroup Grid2dBinGroup::get_igroup_signal_noise_group_min(
  const std::vector<int> vec_igroup ) const
{
  if( vec_igroup.empty() ) THROW_ERROR_NAME("vec_igroup is empty");
  constexpr double big = std::numeric_limits<double>::max();
  double signal_noise_group_min = big;
  Igroup igroup_min = IgroupNotAssigned;
  for( const auto &igroup : vec_igroup ){
    const double signal_noise_group = calc_signal_group(igroup) + calc_noise_poi_group(igroup);
    if( signal_noise_group < signal_noise_group_min ){
      signal_noise_group_min = signal_noise_group;
      igroup_min = igroup;
    }
  }
  return igroup_min;
}

// return the igroup with the largest signal_group + noise_poi_group
// throws if vec_igroup is empty
Igroup Grid2dBinGroup::get_igroup_signal_noise_group_max(
  const std::vector<int> vec_igroup ) const
{
  if( vec_igroup.empty() ) THROW_ERROR_NAME("vec_igroup is empty");
  constexpr double big = std::numeric_limits<double>::max();
  double signal_noise_group_max = -big;
  Igroup igroup_max = IgroupNotAssigned;
  for( const auto &igroup : vec_igroup ){
    const double signal_noise_group = calc_signal_group(igroup) + calc_noise_poi_group(igroup);
    if( signal_noise_group > signal_noise_group_max ){
      signal_noise_group_max = signal_noise_group;
      igroup_max = igroup;
    }
  }
  return igroup_max;
}


// apply_calc_signal_group_to_vec_signal
void Grid2dBinGroup::apply_calc_signal_group_to_vec_signal( )
{
  vec_signal_group.clear();
  const int n_group = get_n_group();
  vec_signal_group.resize(n_group);
  std::set<Igroup> set_igroup = get_set_igroup();
  for( auto igroup : set_igroup ){
    vec_signal_group.at(igroup) = calc_signal_group(igroup);
  }
}

// get ixmin, ixmax, iymin, iymax from igroup
std::array<int,4> Grid2dBinGroup::get_ixmin_ixmax_iymin_iymax( const Igroup igroup ) const
{
  const std::vector<Ixiy> vec_ixiy = get_vec_ixiy(igroup);
  return this->Grid2d::get_ixmin_ixmax_iymin_iymax(vec_ixiy);
};

// return the number of ix and iy in igroup_in
Ixiy Grid2dBinGroup::get_ixiylen( const int igroup_in ) const
{
  const auto [ixmin,ixmax,iymin,iymax] = get_ixmin_ixmax_iymin_iymax(igroup_in);
  return {ixmax-ixmin+1,iymax-iymin+1};
}

// return the number of ix in igroup_in
int Grid2dBinGroup::get_ixlen( const int igroup_in ) const
{
  return get_ixiylen(igroup_in)[0];
}

// return the number of iy in igroup_in
int Grid2dBinGroup::get_iylen( const int igroup_in ) const
{
  return get_ixiylen(igroup_in)[1];
}

// return the vector of (ix, iy) adjacent to the right side of the vec_ixiy forming igroup
std::vector<Ixiy> 
  Grid2dBinGroup::get_adjacent_plus_x_vec_ixiy(
    const int igroup_in, const int ix_delta, const int iy_delta ) const {
  if( std::abs(ix_delta)+std::abs(iy_delta)!=1 ) {
    LOG_ERROR("ix_delta+iy_delta must be 1 or -1: ix_delta={}, iy_delta={}", ix_delta, iy_delta);
    THROW_ERROR("Error! ix_delta+iy_delta must be 1 or -1");
  }
  
  // for output
  std::vector<Ixiy> vec_ixiy_out;

  // get vec_ixiy
  const std::vector<Ixiy> vec_ixiy = get_vec_ixiy(igroup_in);

  // get ixmin, ixmax, iymin iymax of igroup_in
  const auto [ixmin,ixmax,iymin,iymax] = get_ixmin_ixmax_iymin_iymax(igroup_in);

  if( ix_delta== 1 ){
    const int ixx = ixmax+1;
    if( !is_ix_inside(ixx) ) THROW_ERROR2("Error! ixx is out of range",ixx);
    for(int iy=iymin; iy<=iymax; iy++){
      const Ixiy ixiy = {ixx,iy};
      vec_ixiy_out.push_back(ixiy);
    }
  }
  if( ix_delta==-1 ){
    const int ixx = ixmin-1;
    for(int iy=iymin; iy<=iymax; iy++){
    if( !is_ix_inside(ixx) ) THROW_ERROR2("Error! ixx is out of range",ixx);
      const Ixiy ixiy = {ixx,iy};
      vec_ixiy_out.push_back(ixiy);
    }
  }
  
  return vec_ixiy_out;
}


//==================================================================
// merge bin related
//==================================================================

// among the bins that are not available, return the igroup of the bin with the smallest signal_group
int Grid2dBinGroup::get_igroup_signal_group_min() const
{
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  double z_min = 1.0e+10;
  int igroup_min = -1;
  // #pragma omp parallel for private(ix_min,iy_min,z_min)
  for( int iy = 0; iy < nbiny; iy++ ){
    for( int ix = 0; ix < nbinx; ix++ ){
      const Ixiy ixiy = {ix,iy};
      const Igroup igroup = get_igroup(ixiy);
      const double z_sum = calc_signal_group(igroup);
      if( is_avail(ixiy)==true ) continue;
      if( z_sum < z_min ){
        z_min = z_sum;
        igroup_min = igroup;
      }
    }
  }
  return igroup_min;
}

// the current position is Ixiy &ixiy
// merge the single bin at the ix_delta, iy_delta neighbor
// update the map information after the merge
bool Grid2dBinGroup::merge_bin(const Ixiy& ixiy, int ix_delta, int iy_delta)
{
  const auto [ix0, iy0] = ixiy;
  const int ix1 = ix0 + ix_delta;
  const int iy1 = iy0 + iy_delta;

  // return false if ix1, iy1 is out of range
  if(ix1 < 0 || ix1 >= get_nbinx()) return false;
  if(iy1 < 0 || iy1 >= get_nbiny()) return false;

  // build the array of ix1, iy1
  const Ixiy ixiy_new{ix1,iy1};
  // if(is_avail(ixiy_new)) return false;               // do nothing if already used

  const Igroup ig_src = get_igroup(ixiy);      // source group
  const Igroup ig_dst = get_igroup(ixiy_new);  // absorbed group

  if(ig_src == ig_dst) return false;                       // nothing to do if the same group

  // get every ixiy belonging to ig_dst and re-register it into ig_src
  const auto vec_ixiy_dst = bimap_.get_vecMany(ig_dst);
  if( vec_ixiy_dst.empty() ){
    // ! debug
    LOG_ERROR("ig_dst={} has no ixiy", ig_dst);
    return false; // nothing to do if the absorbed group is empty
  }

  for (const auto& ixiy_move : vec_ixiy_dst) {
    bimap_.insert(ig_src, ixiy_move);     // register anew
    bimap_.eraseMany(ixiy_move);          // erase the old registration (reattached to another One)
  }

  // set is_avail of ixiy_new to true
  // ! necessity is questionable
  // set_is_avail(ixiy_new,true);

  return true;
}

// merge the bins of ig_src into ig_dst.
bool Grid2dBinGroup::merge_group(const Igroup ig_src, const Igroup ig_dst)
{
  if( ig_src==ig_dst ) return false;

  // count the ixiy belonging to ig_src
  // return false if the count is 0
  // get the list of Ixiy of ig_src
  const auto vec_ixiy_src = bimap_.get_vecMany(ig_src);
  if (vec_ixiy_src.empty()){
    // ! debug
    LOG_ERROR("ig_src={} has no ixiy", ig_src);
    return false;
  }

  // update the map information
  // copy every ixiy owned by ig_src into ig_dst
  // register all Ixiy of ig_src into ig_dst (bimap_ keeps the consistency automatically)
  for (const auto& ixiy_src : vec_ixiy_src) {
    bimap_.insertOrOverwrite(ig_dst, ixiy_src); // overwrite if already registered
  }

  // erase the registration of ig_src (the Ixiy no longer needed are erased together)
  bimap_.eraseOne(ig_src);
  
#ifdef NODEBUG
#else
  assert(bimap_.isConsistent());
#endif

  return true;
}

// check whether the signal_group of every igroup is greater than or equal to signal_thres.
bool Grid2dBinGroup::is_all_signal_group_larger_than_thres(
  const double signal_thres) const
{
  for (const auto& igroup : bimap_.get_setOne()) {
    const double signal_group = calc_signal_group(igroup);
    if (signal_group < signal_thres) return false;
  }
  return true;
}

// check whether every group satisfies signal_noise_sum >= signal_noise_thres
/// @throws std::runtime_error If a group does not satisfy the condition
bool Grid2dBinGroup::is_all_signal_noise_group_larger_than_thres(
  const double signal_noise_thres ) const
{
  const Detid detid = get_detid();
  const auto set_igroup = get_set_igroup();
  for( const auto igroup : set_igroup ){
    const double signal_group = calc_signal_group(igroup);
    const double noise_poi_group = calc_noise_poi_group(igroup);
    const double signal_noise_group = signal_group + noise_poi_group;
    if( signal_group + noise_poi_group < signal_noise_thres ){
      LOG_WARN("");
      LOG_WARN(" detid={}, igroup={} has signal_group + noise_poi_group = {:E} < thres={:E}"
        , detid, igroup, signal_noise_group, signal_noise_thres);
      return false;
    }
  }
  return true;
}


// merge algorithm
// search around the seed group and merge bins until z_sum reaches thres or above
// the merged bins always form a rectangle.
bool Grid2dBinGroup::merge_rect_bingroup(
  const int igroup_seed, const double thres, const int nloop_max )
{
  // mark the input bins as merged
  set_is_avail(igroup_seed,true);
  
  // set the initial value of z_sum
  double z_sum = calc_signal_group(igroup_seed);

  // return true if z_sum is already thres or above
  if( z_sum >= thres ) return true;

  // search around the bin with the smallest z
  int n_loop=0;
  int ix_delta,iy_delta; // parameters selecting which bin to merge
  while( z_sum < thres ){
    n_loop++;
    if( n_loop%4==0 ){ ix_delta=-1; iy_delta= 0; } // search left
    if( n_loop%4==1 ){ ix_delta= 1; iy_delta= 0; } // search right
    if( n_loop%4==2 ){ ix_delta= 0; iy_delta= 1; } // search up
    if( n_loop%4==3 ){ ix_delta= 0; iy_delta=-1; } // search down

    // get the ixiy of the bins belonging to igroup_seed
    std::vector<Ixiy> vec_ixiy = get_vec_ixiy(igroup_seed);
    for( const auto& ixiy: vec_ixiy ){
      // run merge_bin
      bool is_avail = merge_bin(ixiy,ix_delta,iy_delta);

      // continue if no merge happened.
      if( is_avail == false ) continue;

      // if merged, update z_sum and break
      z_sum = calc_signal_group(igroup_seed);
      if( z_sum >= thres ) break;
    } // vec_ixiy loop end
    if( n_loop > nloop_max ) break;
  } // while loop end
  
  // return true if z_sum is thres or above
  z_sum = calc_signal_group(igroup_seed);
  if( z_sum < thres ) return false;
  return true;
}

// check whether the bins of igroup form a rectangle.
// return false if they do not form a rectangle.
bool Grid2dBinGroup::is_rect( const Igroup igroup ) const
{
  const auto vec_ixiy = get_vec_ixiy(igroup);
  const auto [ixmin,ixmax,iymin,iymax] = get_ixmin_ixmax_iymin_iymax(igroup);
  const int ixlen = ixmax - ixmin + 1;
  const int iylen = iymax - iymin + 1;
  if( vec_ixiy.size() != ixlen*iylen ) return false;
  return true;
}

// new version of divide
// divide igroup_seed by nx_div, ny_div.
// the igroup of the newly created bins starts at the maximum igroup + 1.
// the remaining bins are gathered into the bin specified by ix_extra, iy_extra.
// returns the vector of the newly created (igroup, ix_div, iy_div).
std::vector<std::array<int,3>> Grid2dBinGroup::divide_extra(
  const Igroup igroup_seed
, const int nx_div, const int ny_div
, const int ix_extra, const int iy_extra
, const bool exec_move_ixiy )
{
  // vector for the output
  std::vector<std::array<int,3>> vec_igixiy_div;

  // check nx_div, ny_div
  if( nx_div < 1 ) THROW_ERROR("nx_div must be >= 1");
  if( ny_div < 1 ) THROW_ERROR("ny_div must be >= 1");
  if( nx_div==1 && ny_div==1 ){
    LOG_WARN("nx_div=1 and ny_div=1. Nothing to do.");
    return vec_igixiy_div;
  }

  // get the ix, iy belonging to igroup_seed
  const auto vec_ixiy = get_vec_ixiy(igroup_seed);
  if (vec_ixiy.empty()) return vec_igixiy_div; // fail when the group has no bin

  if(bimap_.countMany() == 0) {
    // handling when the map is empty (log an error message)
    LOG_ERROR("bimap_ is empty. return empty result.");
    return vec_igixiy_div;
  }

  // start index of the new igroup = maximum existing igroup + 1
  const auto set_old = bimap_.get_setOne();
  const int igroup_max = *set_old.rbegin(); // get the maximum igroup
  int new_igroup_index = igroup_max + 1;
  // hoge(new_igroup_index);

  // get_vec_vec_ixiy returns, as a 2D grid, the smallest rectangle containing the Ixiy coordinates of vec_ixiy
  // that is, for std::vector<Ixiy> input = { {1, 1}, {2, 3} }; the output of get_vec_vec_ixiy is
  // vec_vec_ixiy[0][0] = {1, 1}, 
  // vec_vec_ixiy[1][0] = {2, 1},
  // vec_vec_ixiy[0][1] = {1, 2},
  // vec_vec_ixiy[1][1] = {2, 2},
  // vec_vec_ixiy[0][2] = {1, 3},
  // vec_vec_ixiy[1][2] = {2, 3}: a 2D array of Ixiy coordinates.
  const auto vec_vec_ixiy = this->Grid2d::get_vec_vec_ixiy(vec_ixiy);
  
  // vec_x4_ixiy holds the vec_vec_ixiy above inside an ny_div x nx_div 2D array.
  auto vec_x4_ixiy = myapp::split2DArrayIntoGrid( vec_vec_ixiy, nx_div, ny_div, ix_extra, iy_extra );
  
  for(int iy_div=0;iy_div<ny_div;iy_div++){
    for(int ix_div=0;ix_div<nx_div;ix_div++){
      const auto& vec_vec_ixiy_tmp = vec_x4_ixiy.at(iy_div).at(ix_div);
      const auto& vec_ixiy_tmp = this->Grid2d::get_vec_ixiy(vec_vec_ixiy_tmp);
      if( vec_ixiy_tmp.empty() ) continue;
      if( exec_move_ixiy ){
        for (const Ixiy& ixiy : vec_ixiy_tmp)
          bimap_.insertOrOverwrite(new_igroup_index, ixiy); // register into the new igroup
      }

      // store new_igroup_index, ix_div, iy_div
      if( new_igroup_index%100==0 ) fprintf(stderr,"igroup_seed=%d : ",igroup_seed);
      if( exec_move_ixiy && new_igroup_index%100==0 ) disp_ixmin_ixmax_iymin_iymax(stderr,new_igroup_index);
      vec_igixiy_div.push_back( {new_igroup_index, ix_div, iy_div} );
      ++new_igroup_index;
    }
  }
  // erase the old igroup after everything is done
  // bimap_.eraseOne(igroup_seed);
  
  return vec_igixiy_div;
}

// divide igroup.
// repeat divide until signal_group becomes smaller than signal_noise_group_trigger.
// divide is applied in the order (nx_div=2, ny_div=1), then (nx_div=1, ny_div=2).
// the remainder is given to the side with the smaller signal_group.
// nloop_limit is the maximum number of divide loops.
// see images/auto_divide_by_signal_noise_group.draw.png
// returns bool,int | bool: true if divided, false if not divided, int: number of loops performed
//
// Optimized version: uses calc_trial_split_max_half() to determine which half
// has larger signal+noise, avoiding a full Grid2dBinGroup deep copy per iteration.
std::tuple<bool,int> Grid2dBinGroup::auto_divide_by_signal_noise_group(
  const int igroup_start, const double signal_noise_group_trigger, const int nloop_limit,
  const int ixlen_min, const int iylen_min, const bool tf_prefer_split_x )
{
  const int n_group = get_n_group();
  if( n_group%100==0 ) fprintf(stderr,"n_group=%d, ",n_group);
  // igroup storage
  std::vector<Igroup> vec_igroup;
  vec_igroup.push_back(igroup_start);

  // igroup_new storage
  std::vector<Igroup> vec_igroup_new;
  std::vector<std::array<int,3>> vec_igixiy_div;

  int nloop = -1;

  // note that vec_igroup grows inside the loop, so its size changes dynamically
  for(int i = 0; i < (int)vec_igroup.size(); i++ ){
    const Igroup igroup = vec_igroup.at(i);
    nloop++;
    if( nloop > nloop_limit ) break;

    // continue if the igroup does not exist
    if(!bimap_.hasOne(igroup)) continue;

    // calc ixlen, iylen and bounding box (used for both threshold check and split)
    const auto [ixmin, ixmax, iymin, iymax] = get_ixmin_ixmax_iymin_iymax(igroup);
    const int ixlen = ixmax - ixmin + 1;
    const int iylen = iymax - iymin + 1;

    // calc signal+noise directly from rectangular region (no bimap_ lookup)
    const double signal_noise_group = calc_signal_noise_rect(ixmin, ixmax, iymin, iymax);

    // do nothing if signal_noise_group is below the threshold
    if( signal_noise_group < signal_noise_group_trigger ) continue;

    //------------------------------------------
    // decide the split direction
    //------------------------------------------
    // X-direction split: the condition is "(ixlen > iylen) or (ixlen == iylen and tf_prefer_split_x==true)" and ixlen > ixlen_min
    if( (ixlen > iylen || (ixlen == iylen && tf_prefer_split_x)) && ixlen > ixlen_min ){

      // Determine which half has larger signal+noise without copying Grid2dBinGroup.
      // Reproduces: tmp_g2bg.divide_extra(igroup, 2, 1, 1, 0, true)
      //             tmp_g2bg.get_igroup_signal_noise_group_max(vec_igroup_new)
      const auto [max_half, n_halves] = calc_trial_split_max_half(
        ixmin, ixmax, iymin, iymax, /*split_x=*/true);

      if (n_halves >= 2) {
        // Both halves are valid. max_half is the ix_div index with larger signal+noise.
        // The original code collected all (ix_div, iy_div) matching igroup_max.
        // For X-split (nx_div=2, ny_div=1), each half has exactly one (ix_div, iy_div):
        //   half 0 -> (0, 0), half 1 -> (1, 0).
        // When vec_ixiy_div.size() == 1, divide_extra is called with ix_extra = ix_div+1.
        // When vec_ixiy_div.size() >= 2, this can't happen for nx_div=2 (at most 1 match).
        // So the path is always size()==1, ix_extra = max_half + 1.
        vec_igixiy_div = divide_extra(igroup, 2, 1, max_half + 1, 0, true);
        vec_igroup_new = myapp::copyArrayElementToVector(vec_igixiy_div, 0);
        vec_igroup.insert(vec_igroup.end(), vec_igroup_new.begin(), vec_igroup_new.end());
      }
      // n_halves == 1: region too small to split into 2 halves, skip

      vec_igixiy_div.clear();
      vec_igroup_new.clear();
      continue;
    } // X-direction split end

    // Y-direction split: the condition is "(ixlen < iylen) or (ixlen == iylen and tf_prefer_split_x==false)" and iylen > iylen_min
    if( (ixlen < iylen || (ixlen == iylen && !tf_prefer_split_x)) && iylen > iylen_min ){

      // Determine which half has larger signal+noise without copying Grid2dBinGroup.
      // Reproduces: tmp_g2bg.divide_extra(igroup, 1, 2, 0, 1, true)
      //             tmp_g2bg.get_igroup_signal_noise_group_max(vec_igroup_new)
      const auto [max_half, n_halves] = calc_trial_split_max_half(
        ixmin, ixmax, iymin, iymax, /*split_x=*/false);

      if (n_halves >= 2) {
        // For Y-split (nx_div=1, ny_div=2), each half has (0, iy_div).
        // size()==1 path: iy_extra = max_half + 1.
        vec_igixiy_div = divide_extra(igroup, 1, 2, 0, max_half + 1, true);
        vec_igroup_new = myapp::copyArrayElementToVector(vec_igixiy_div, 0);
        vec_igroup.insert(vec_igroup.end(), vec_igroup_new.begin(), vec_igroup_new.end());
      }

      vec_igixiy_div.clear();
      vec_igroup_new.clear();
      continue;
    } // Y-direction split end

  } // for loop end

  if (nloop <= 0) return std::make_tuple(false, nloop);

  return std::make_tuple(true, nloop);
}

//######################################################################
// Legacy version of auto_divide_by_signal_noise_group (for regression testing).
// Uses full Grid2dBinGroup copy for trial splits.
//######################################################################
std::tuple<bool,int> Grid2dBinGroup::auto_divide_by_signal_noise_group_legacy(
  const int igroup_start, const double signal_noise_group_trigger, const int nloop_limit,
  const int ixlen_min, const int iylen_min, const bool tf_prefer_split_x )
{
  const int n_group = get_n_group();
  if( n_group%100==0 ) fprintf(stderr,"n_group=%d, ",n_group);
  std::vector<Igroup> vec_igroup;
  vec_igroup.push_back(igroup_start);

  std::vector<Igroup> vec_igroup_new;
  std::vector<std::array<int,3>> vec_igixiy_div;

  int nloop = -1;

  for(int i = 0; i < (int)vec_igroup.size(); i++ ){
    const Igroup igroup = vec_igroup.at(i);
    nloop++;
    if( nloop > nloop_limit ) break;

    if(!bimap_.hasOne(igroup)) continue;

    const auto vec_ixiy = get_vec_ixiy(igroup);
    if(vec_ixiy.empty()) continue;

    const auto [ixlen, iylen] = get_ixiylen(igroup);
    double signal_noise_group = calc_signal_group(igroup) + calc_noise_poi_group(igroup);
    if( signal_noise_group < signal_noise_group_trigger ) continue;

    // X-direction split
    if( (ixlen > iylen || (ixlen == iylen && tf_prefer_split_x)) && ixlen > ixlen_min ){
      Grid2dBinGroup tmp_g2bg = *this;
      vec_igixiy_div = tmp_g2bg.divide_extra(igroup, 2, 1, 1, 0, true);
      vec_igroup_new = myapp::copyArrayElementToVector(vec_igixiy_div, 0);
      Igroup igroup_signal_noise_group_max
        = tmp_g2bg.get_igroup_signal_noise_group_max(vec_igroup_new);

      std::vector<Ixiy> vec_ixiy_div;
      for( const auto& [ig, ix_div, iy_div] : vec_igixiy_div ){
        if( ig == igroup_signal_noise_group_max )
          vec_ixiy_div.push_back({ix_div, iy_div});
      }

      if( vec_ixiy_div.empty() ) continue;
      if( vec_ixiy_div.size() >= 2 ){
        vec_igixiy_div = divide_extra(igroup, 2, 1, 1, 0, true);
        vec_igroup_new = myapp::copyArrayElementToVector(vec_igixiy_div, 0);
        vec_igroup.insert(vec_igroup.end(), vec_igroup_new.begin(), vec_igroup_new.end());
      }
      if( vec_ixiy_div.size() == 1 ){
        const auto [ix_div, iy_div] = vec_ixiy_div.at(0);
        vec_igixiy_div = divide_extra(igroup, 2, 1, ix_div + 1, 0, true);
        vec_igroup_new = myapp::copyArrayElementToVector(vec_igixiy_div, 0);
        vec_igroup.insert(vec_igroup.end(), vec_igroup_new.begin(), vec_igroup_new.end());
      }
      vec_igixiy_div.clear();
      vec_igroup_new.clear();
      continue;
    }

    // Y-direction split
    if( (ixlen < iylen || (ixlen == iylen && !tf_prefer_split_x)) && iylen > iylen_min ){
      Grid2dBinGroup tmp_g2bg = *this;
      vec_igixiy_div = tmp_g2bg.divide_extra(igroup, 1, 2, 0, 1, true);
      vec_igroup_new = myapp::copyArrayElementToVector(vec_igixiy_div, 0);
      Igroup igroup_signal_noise_group_max
        = tmp_g2bg.get_igroup_signal_noise_group_max(vec_igroup_new);

      std::vector<Ixiy> vec_ixiy_div;
      for( const auto& [ig, ix_div, iy_div] : vec_igixiy_div ){
        if( ig == igroup_signal_noise_group_max )
          vec_ixiy_div.push_back({ix_div, iy_div});
      }
      vec_igixiy_div.clear();

      if( vec_ixiy_div.empty() ) continue;
      if( vec_ixiy_div.size() >= 2 ){
        vec_igixiy_div = divide_extra(igroup, 1, 2, 0, 1, true);
        vec_igroup_new = myapp::copyArrayElementToVector(vec_igixiy_div, 0);
        vec_igroup.insert(vec_igroup.end(), vec_igroup_new.begin(), vec_igroup_new.end());
      }
      if( vec_ixiy_div.size() == 1 ){
        const auto [ix_div, iy_div] = vec_ixiy_div.at(0);
        vec_igixiy_div = divide_extra(igroup, 1, 2, 0, iy_div + 1, true);
        vec_igroup_new = myapp::copyArrayElementToVector(vec_igixiy_div, 0);
        vec_igroup.insert(vec_igroup.end(), vec_igroup_new.begin(), vec_igroup_new.end());
      }
      vec_igixiy_div.clear();
      vec_igroup_new.clear();
      continue;
    }

  } // for loop end

  if (nloop <= 0) return std::make_tuple(false, nloop);
  return std::make_tuple(true, nloop);
}

// @brief   divide igroup.
// @details Repeats divide until signal_group becomes smaller than signal_noise_group_trigger. \n
// divide is applied in the order (nx_div=2, ny_div=1), then (nx_div=1, ny_div=2). \n
// The remainder is given to the side with the smaller signal_group. \n
// nloop_limit is the maximum number of divide loops.
std::tuple<bool,int> Grid2dBinGroup::auto_divide_by_signal_noise_group(
  const Igroup igroup, const Grid2dBinGroup::Parameters &prm_bingrp )
{
  return Grid2dBinGroup::auto_divide_by_signal_noise_group(
    igroup, prm_bingrp.signal_noise_group_trig, prm_bingrp.nloop_limit,
    prm_bingrp.ixlen_min, prm_bingrp.iylen_min, prm_bingrp.tf_prefer_split_x );
}

// optimize the size of all bins with signal_group_trg, nloop_limit, ixlen_min, iylen_min.
// returns the final number of groups.
// @ref auto_divide_by_signal_noise_group(int,double,int,int,int)
int Grid2dBinGroup::auto_divide_by_signal_noise_group_all(
  const double signal_noise_group_trigger, const int nloop_limit,
  const int ixlen_min, const int iylen_min, const bool tf_prefer_split_x )
{
  const std::set<Igroup> set_igroup_before = get_set_igroup();

  // run auto divide for each group
  for (const auto& igroup : set_igroup_before) {
    const auto [is_done, nloop] = auto_divide_by_signal_noise_group(
      igroup, signal_noise_group_trigger, nloop_limit, ixlen_min, iylen_min, tf_prefer_split_x);
    fprintf(stderr, "igroup = %d, is_done = %d\r", igroup, is_done);
  }
  std::set<Igroup> set_igroup_after = get_set_igroup();

  set_done_grouping(true);

  // renumbering igroup
  renumbering_groups();
  
  set_igroup_after = get_set_igroup();

  return static_cast<int>(set_igroup_after.size());
}
// @brief   optimize the size of all bins with signal_group_trg, nloop_limit, ixlen_min, iylen_min.
// @details Returns the final number of groups.
int Grid2dBinGroup::auto_divide_by_signal_noise_group_all(const Grid2dBinGroup::Parameters &prm_bingrp)
{
  return Grid2dBinGroup::auto_divide_by_signal_noise_group_all(
    prm_bingrp.signal_noise_group_trig, prm_bingrp.nloop_limit,
    prm_bingrp.ixlen_min, prm_bingrp.iylen_min, prm_bingrp.tf_prefer_split_x);
}

// calculate signal_group / noise_poi_group for each group and set them into vec_signal_group / vec_noise_poi_group.
void Grid2dBinGroup::calc_set_vec_signal_noise_group()
{
  const volatile int n_group = get_n_group();
  LOG_DEBUG("calc_set_vec_signal_noise_group : detid={}, n_group={}",detid_, n_group);

  vec_signal_group.clear();
  vec_signal_group.resize( n_group, 0.0 );
  vec_noise_poi_group.clear();
  vec_noise_poi_group.resize( n_group, 0.0 );
  vec_noise_det_group.clear();
  vec_noise_det_group.resize( n_group, 0.0 );

  vec_signal_group_poisson.clear();
  vec_signal_group_poisson.resize( n_group, 0.0 );
  vec_noise_poi_group_poisson.clear();
  vec_noise_poi_group_poisson.resize( n_group, 0.0 );

  // calculate noise_group and set it into vec_noise_poi_group
  // det (floor) gets no fluctuation; only poi (Poisson-bucket) is passed through poisson().
  // const std::vector<Igroup> vec_igroup = get_vec_igroup(); // ! this is wrong. it errors.
  const std::set<Igroup> set_igroup = get_set_igroup();
  #pragma omp parallel for
  for( int i =0; i<n_group; i++ ){
    const Igroup igroup = *std::next(set_igroup.begin(), i);
    const double signal_group = calc_signal_group(igroup);
    const double noise_poi_group = calc_noise_poi_group(igroup);
    const double noise_det_group = calc_noise_det_group(igroup);
    const double signal_group_poisson = static_cast<double>(myapp::poisson(signal_group));
    const double noise_poi_group_poisson = static_cast<double>(myapp::poisson(noise_poi_group));
    set_signal_group(igroup, signal_group);
    set_noise_poi_group(igroup, noise_poi_group);
    set_noise_det_group(igroup, noise_det_group);
    set_signal_group_poisson(igroup, signal_group_poisson);
    set_noise_poi_group_poisson(igroup, noise_poi_group_poisson);
  }
}

// set signal_group_poisson and noise_poi_group_poisson from signal_group and noise_poi_group with a newly drawn Poisson error.
void Grid2dBinGroup::set_diff_poisson_again()
{
  const int n_group = get_n_group();
  vec_signal_group_poisson.clear();
  vec_signal_group_poisson.resize( n_group, 0.0 );
  vec_noise_poi_group_poisson.clear();
  vec_noise_poi_group_poisson.resize( n_group, 0.0 );

  // get the signal_group and noise_poi_group values and reset vec_signal_group_poisson / vec_noise_poi_group_poisson
  const std::vector<Igroup> vec_igroup = get_vec_igroup();
  #pragma omp parallel for
  for( int i =0; i<n_group; i++ ){
    const Igroup igroup = vec_igroup.at(i);
    const double signal_group = get_signal_group(igroup);
    const double noise_poi_group = get_noise_poi_group(igroup);
    const double signal_group_poisson = static_cast<double>(myapp::poisson(signal_group));
    const double noise_poi_group_poisson = static_cast<double>(myapp::poisson(noise_poi_group));
    set_signal_group_poisson(igroup, signal_group_poisson);
    set_noise_poi_group_poisson(igroup, noise_poi_group_poisson);
  }
}

// renumber the igroup of bimap_ starting from igroup_start.
void Grid2dBinGroup::renumbering_groups(const Igroup igroup_start)
{
  BimapIgroupIxiy_IyMajor bimap_new;

  Igroup igroup_new = igroup_start;
  const auto set_old_igroups = bimap_.get_setOne();
  const int n_group = static_cast<int>(set_old_igroups.size());

  for (const auto& ig_old : set_old_igroups) {
    for (const auto& ixiy : bimap_.get_vecMany(ig_old)) {
      bimap_new.insert(igroup_new, ixiy);
    }
    ++igroup_new;
  }

  if (igroup_new - igroup_start != n_group)
    THROW_ERROR("Grid2dBinGroup::renumbering_groups: Group count mismatch. expected={}, actual={}", n_group, igroup_new - igroup_start);

  bimap_ = std::move(bimap_new);
}

//===================================================
// file I/O
//===================================================

// write the contents of bimap_ in ASCII to the file given by the fs::path argument.
// width is the fixed width of each field (default 9).
// rows are sorted in the order Detid, Igroup, Iy, Ix
void Grid2dBinGroup::write_bimap_all_to_ascii(
  const fs::path& pathout, const int width) const
{
  std::ofstream ofs(pathout);
  if (!ofs) {
    throw std::runtime_error("Failed to open file: " + pathout.string());
  }

  struct Info {
    std::array<int, 4> idx;  // {Detid, Igroup, Iy, Ix}
  };

  const auto map_tmp = bimap_.getMapOneMany();
  const size_t n = map_tmp.size();
  std::vector<Info> vec(n);

  // cache the boundary values of each Igroup
  std::unordered_map<Igroup, std::array<double, 4>> bounds_map;
  for (const Igroup ig : get_set_igroup()) {
    bounds_map[ig] = this->get_xmin_xmax_ymin_ymax(ig);
  }

  size_t index = 0;
  for (const auto& [igroup, ixiy] : map_tmp) {
    vec.at(index).idx = {detid_, igroup, ixiy[1], ixiy[0]};
    ++index;
  }

  // sort
  std::sort(vec.begin(), vec.end(),
    [](const Info& a, const Info& b) {
      return a.idx < b.idx;
    });

  // output
  ofs << std::fixed << std::setprecision(4);
  Igroup igroup_prev = -std::numeric_limits<Igroup>::max();
  for (const auto& entry : vec) {
    const Igroup igroup = entry.idx[1];
    if (igroup != igroup_prev) {

      // comment row: legend
      ofs << "# "
          << std::setw(width-2) << std::right << "Detid"  << ' '
          << std::setw(width) << std::right << "Igroup" << ' '
          << std::setw(width) << std::right << "xmin"   << ' '
          << std::setw(width) << std::right << "xmax"   << ' '
          << std::setw(width) << std::right << "ymin"   << ' '
          << std::setw(width) << std::right << "ymax"   << '\n';

      const auto& b = bounds_map.at(igroup);
      // output row: Detid, Igroup, Ix, Iy
      ofs << std::setw(width) << std::right << entry.idx[0] << ' '
          << std::setw(width) << std::right << igroup       << ' '
          << std::setw(width) << b[0] << ' '
          << std::setw(width) << b[1] << ' '
          << std::setw(width) << b[2] << ' '
          << std::setw(width) << b[3] << '\n'
          << "# "
          << std::setw(width-2) << std::right << "Ix"     << ' '
          << std::setw(width) << std::right << "Iy"     << '\n';
      igroup_prev = igroup;
    }
    // output row: Ix, Iy
    ofs << std::setw(width) << std::right << entry.idx[3] << ' '
        << std::setw(width) << std::right << entry.idx[2] << '\n';
  }
}


//  write the header information of bimap_ in ASCII format.
void Grid2dBinGroup::write_bimap_header_to_ascii(
  const fs::path& pathout, const int width) const
{
  std::ofstream ofs(pathout);
  if (!ofs) {
    throw std::runtime_error("Failed to open file: " + pathout.string());
  }

  // get set of Igroup
  const auto set_igroup = bimap_.get_setOne();
  if (set_igroup.empty()){
    LOG_WARN("write_bimap_header_to_ascii: set_igroup is empty. detid={}", detid_);
    SLEEP_MSEC(500);
    return;
  }
  
  // store to vector of <detid, igroup
  // , ixmin, ixmax, iymin, iymax, ixlen, iylen
  // , xmin, xmax, ymin, ymax, xlen, ylen>
  const int n_group = static_cast<int>(set_igroup.size());

  struct HeaderInfo {
    int detid;
    Igroup igroup;
    int ixmin, ixmax, iymin, iymax;
    int ixlen, iylen;
    double xmin, xmax, ymin, ymax;
    double xlen, ylen;
  };


  std::vector<HeaderInfo> vec_info(n_group);
  size_t index = 0;
  for (const Igroup igroup : set_igroup) {
    const auto [ixmin, ixmax, iymin, iymax] = get_ixmin_ixmax_iymin_iymax(igroup);
    const auto [xmin, xmax, ymin, ymax] = get_xmin_xmax_ymin_ymax(igroup);
    const int ixlen = ixmax - ixmin + 1;
    const int iylen = iymax - iymin + 1;
    const double xlen = xmax - xmin;
    const double ylen = ymax - ymin;
    HeaderInfo info = {
        detid_, igroup
      , ixmin, ixmax, iymin, iymax, ixlen, iylen
      , xmin, xmax, ymin, ymax, xlen, ylen
    };
    vec_info.at(index) = info;
    ++index;
  }

  // sort in ascending order of detid_, igroup
  std::sort(vec_info.begin(), vec_info.end(),
    [](const HeaderInfo& a, const HeaderInfo& b) {
      if (a.detid != b.detid) return a.detid < b.detid;
      if (a.igroup != b.igroup) return a.igroup < b.igroup;
      return a.ixmin < b.ixmin; // tie-break by ixmin
    });

  // output the header row
  ofs << "# "
      << std::setw(width-2) << std::right << "Detid"  << ' '
      << std::setw(width) << std::right << "Igroup" << ' '
      << std::setw(width) << std::right << "Ixmin"  << ' '
      << std::setw(width) << std::right << "Ixmax"  << ' '
      << std::setw(width) << std::right << "Iymin"  << ' '
      << std::setw(width) << std::right << "Iymax"  << ' '
      << std::setw(width) << std::right << "Ixlen"  << ' '
      << std::setw(width) << std::right << "Iylen"  << ' '
      << std::setw(width) << std::right << "Xmin"   << ' '
      << std::setw(width) << std::right << "Xmax"   << ' '
      << std::setw(width) << std::right << "Ymin"   << ' '
      << std::setw(width) << std::right << "Ymax"   << ' '
      << std::setw(width) << std::right << "Xlen"   << ' '
      << std::setw(width) << std::right << "Ylen"   << '\n';
  
  // output the information of each group
  ofs << std::fixed << std::setprecision(4);
  for (const auto& info : vec_info) {
    ofs << std::setw(width) << std::right << info.detid << ' '
        << std::setw(width) << std::right << info.igroup << ' '
        << std::setw(width) << std::right << info.ixmin << ' '
        << std::setw(width) << std::right << info.ixmax << ' '
        << std::setw(width) << std::right << info.iymin << ' '
        << std::setw(width) << std::right << info.iymax << ' '
        << std::setw(width) << std::right << info.ixlen << ' '
        << std::setw(width) << std::right << info.iylen << ' '
        << std::setw(width) << std::right << info.xmin << ' '
        << std::setw(width) << std::right << info.xmax << ' '
        << std::setw(width) << std::right << info.ymin << ' '
        << std::setw(width) << std::right << info.ymax << ' '
        << std::setw(width) << std::right << info.xlen << ' '
        << std::setw(width) << std::right << info.ylen << '\n';
  }
}



// save all member variables to std::ofstream
void Grid2dBinGroup::save(std::ofstream& ofs) const
{
  Grid2d::save(ofs);
  io_binary::write_string(ofs, name);
  save_bimap(ofs);
  io_binary::write_binary<int>(ofs, detid_);  
  io_binary::write_vec_bool(ofs, vec_is_avail_group );
  io_binary::write_vec<double>(ofs, vec_signal_group );
  io_binary::write_vec<double>(ofs, vec_noise_poi_group );
  io_binary::write_vec<double>(ofs, vec_noise_det_group );
  io_binary::write_vec<double>(ofs, vec_signal_group_poisson );
  io_binary::write_vec<double>(ofs, vec_noise_poi_group_poisson );
  io_binary::write_vec_vec_bool(ofs, vec_vec_is_avail );
  io_binary::write_vec_vec<double>(ofs, vec_vec_signal );
  io_binary::write_vec_vec<double>(ofs, vec_vec_noise_poi );
  io_binary::write_vec_vec<double>(ofs, vec_vec_noise_det );
  io_binary::write_bool(ofs, done_grouping );
  io_binary::write_vec<double>(ofs, vec_dens_group_lower );
  io_binary::write_vec<double>(ofs, vec_dens_group_center );
  io_binary::write_vec<double>(ofs, vec_dens_group_upper );
  io_binary::write_vec<double>(ofs, vec_delta_nmuon_group_lower );
  io_binary::write_vec<double>(ofs, vec_delta_nmuon_group_center );
  io_binary::write_vec<double>(ofs, vec_delta_nmuon_group_upper );
  io_binary::write_vec<double>(ofs, vec_volume_group );
  // group efficiencies (added in PIPELINE_VERSION 7)
  io_binary::write_vec<double>(ofs, vec_eff_low_group );
  io_binary::write_vec<double>(ofs, vec_eff_cnt_group );
  io_binary::write_vec<double>(ofs, vec_eff_upp_group );

  if (ofs.fail()) THROW_ERROR("Grid2dBinGroup::save failed");
}

// load all member variables from std::ifstream
void Grid2dBinGroup::load(std::ifstream& ifs)
{
  Grid2d::load(ifs);
  name    = io_binary::read_string(ifs);
  load_bimap(ifs);
  detid_  = io_binary::read_binary<int>(ifs);
  vec_is_avail_group       = io_binary::read_vec_bool(ifs);
  vec_signal_group         = io_binary::read_vec<double>(ifs);
  vec_noise_poi_group          = io_binary::read_vec<double>(ifs);
  vec_noise_det_group      = io_binary::read_vec<double>(ifs);
  vec_signal_group_poisson = io_binary::read_vec<double>(ifs);
  vec_noise_poi_group_poisson = io_binary::read_vec<double>(ifs);
  vec_vec_is_avail         = io_binary::read_vec_vec_bool(ifs);
  vec_vec_signal           = io_binary::read_vec_vec<double>(ifs);
  vec_vec_noise_poi            = io_binary::read_vec_vec<double>(ifs);
  vec_vec_noise_det        = io_binary::read_vec_vec<double>(ifs);
  done_grouping            = io_binary::read_bool(ifs);
  vec_dens_group_lower     = io_binary::read_vec<double>(ifs);
  vec_dens_group_center    = io_binary::read_vec<double>(ifs);
  vec_dens_group_upper     = io_binary::read_vec<double>(ifs);
  vec_delta_nmuon_group_lower  = io_binary::read_vec<double>(ifs);
  vec_delta_nmuon_group_center = io_binary::read_vec<double>(ifs);
  vec_delta_nmuon_group_upper  = io_binary::read_vec<double>(ifs);
  vec_volume_group         = io_binary::read_vec<double>(ifs);
  // group efficiencies (added in PIPELINE_VERSION 7)
  vec_eff_low_group        = io_binary::read_vec<double>(ifs);
  vec_eff_cnt_group        = io_binary::read_vec<double>(ifs);
  vec_eff_upp_group        = io_binary::read_vec<double>(ifs);

  if (ifs.fail()) THROW_ERROR("Grid2dBinGroup::load failed");
}

//==============================================================================
// Functions moved from hpp to cpp (Phase 3)
//==============================================================================

// disp all group range
void Grid2dBinGroup::disp_group_range_all() const
{
  std::set<Igroup> set_igroup = pdic_->get_set_igroup(detid_);
  for( const auto &igroup : set_igroup ){
    disp_ixmin_ixmax_iymin_iymax(stdout,igroup);
  }
}

// disp ixmin_ixmax_iymin_iymax
void Grid2dBinGroup::disp_ixmin_ixmax_iymin_iymax( const Igroup igroup ) const
{
  const auto [ixmin,ixmax,iymin,iymax] = get_ixmin_ixmax_iymin_iymax(igroup);
  std::cout << "igroup = " << igroup << ", range:"
   << ", ix :  " << ixmin << " - " << ixmax
   << ", iy :  " << iymin << " - " << iymax << std::endl;
}

// disp ixmin_ixmax_iymin_iymax
void Grid2dBinGroup::disp_ixmin_ixmax_iymin_iymax( FILE* fout, const Igroup igroup ) const
{
  const auto [ixmin,ixmax,iymin,iymax] = get_ixmin_ixmax_iymin_iymax(igroup);
  fprintf(fout, "igroup = %5d, range:"
    ", ix :  %5d - %5d"
    ", iy :  %5d - %5d\r", igroup, ixmin, ixmax, iymin, iymax );
}

// save bimap_ to std::ofstream
void Grid2dBinGroup::save_bimap(std::ofstream& ofs) const
{
  io_binary::write_multimap_int_int2(ofs, bimap_.getMapOneMany());
}

// load bimap_ from std::ifstream
void Grid2dBinGroup::load_bimap(std::ifstream& ifs)
{
  bimap_.clear();
  const auto mmap = io_binary::read_multimap_int_int2(ifs);  // multimap<int, array<int,2>>
  for (const auto& [igroup, ixiy] : mmap) bimap_.insert(igroup, ixiy); // insert from mmap<int, array<int,2>> into bimap_
}

// save all member variables to std::filesystem::path
void Grid2dBinGroup::save( const std::filesystem::path& pathout ) const
{
  std::ofstream ofs = io_binary::open_ofstream(pathout);
  // first write the architecture information
  io_binary::ArchitectureInfo currentInfo = io_binary::get_current_architecture_info();
  io_binary::write_architecture_info(ofs, currentInfo);
  // then write the data
  save(ofs); ofs.close();
}

// load all member variables from std::filesystem::path
void Grid2dBinGroup::load( const std::filesystem::path &path_in )
{
  std::ifstream ifs = io_binary::open_ifstream(path_in);
  // 1) read the architecture information stored in the file
  io_binary::ArchitectureInfo fileInfo = io_binary::read_architecture_info(ifs);
  // 2) check compatibility with the current environment (throws on mismatch)
  io_binary::check_architecture_compatibility_or_throw(fileInfo);
  // 3) read the data
  load(ifs); ifs.close();
}

// allocate memory for vec_is_avail_group
void Grid2dBinGroup::resize_vec_is_avail_group(
  const int n_group_in, const bool is_avail_init)
{
  resize_group_vector(
    vec_is_avail_group, n_group_in, is_avail_init);
}

// allocate memory for vec_signal_group
void Grid2dBinGroup::resize_vec_signal_group(
  const int n_group_in, const double signal_init)
{
  resize_group_vector(
    vec_signal_group, n_group_in, signal_init);
}

// allocate memory for vec_noise_poi_group
void Grid2dBinGroup::resize_vec_noise_poi_group(
  const int n_group_in, const double noise_init)
{
  resize_group_vector(
    vec_noise_poi_group, n_group_in, noise_init);
}

// allocate memory for vec_noise_det_group (deterministic floor noise per group, det)
void Grid2dBinGroup::resize_vec_noise_det_group(
  const int n_group_in, const double noise_init)
{
  resize_group_vector(
    vec_noise_det_group, n_group_in, noise_init);
}

// allocate memory for vec_signal_group_poisson
void Grid2dBinGroup::resize_vec_signal_group_poisson(
  const int n_group_in, const double signal_init)
{
  resize_group_vector(
    vec_signal_group_poisson, n_group_in, signal_init);
}

// allocate memory for vec_noise_poi_group_poisson
void Grid2dBinGroup::resize_vec_noise_poi_group_poisson(
  const int n_group_in, const double noise_init)
{
  resize_group_vector(
    vec_noise_poi_group_poisson, n_group_in, noise_init);
}

// allocate memory for vec_dens_group_lower
void Grid2dBinGroup::resize_vec_dens_group_lower(
  const int n_group_in, const double dens_init)
{
  resize_group_vector(
    vec_dens_group_lower, n_group_in, dens_init);
}

// allocate memory for vec_dens_group_center
void Grid2dBinGroup::resize_vec_dens_group_center(
  const int n_group_in, const double dens_init)
{
  resize_group_vector(
    vec_dens_group_center, n_group_in, dens_init);
}

// allocate memory for vec_dens_group_upper
void Grid2dBinGroup::resize_vec_dens_group_upper(
  const int n_group_in, const double dens_init)
{
  resize_group_vector(
    vec_dens_group_upper, n_group_in, dens_init);
}

// allocate memory for vec_delta_nmuon_group_lower
void Grid2dBinGroup::resize_vec_delta_nmuon_group_lower(
  const int n_group_in, const double delta_nmuon_init)
{
  resize_group_vector(
    vec_delta_nmuon_group_lower, n_group_in, delta_nmuon_init);
}

// allocate memory for vec_delta_nmuon_group_center
void Grid2dBinGroup::resize_vec_delta_nmuon_group_center(
  const int n_group_in, const double delta_nmuon_init)
{
  resize_group_vector(
    vec_delta_nmuon_group_center, n_group_in, delta_nmuon_init);
}

// allocate memory for vec_delta_nmuon_group_upper
void Grid2dBinGroup::resize_vec_delta_nmuon_group_upper(
  const int n_group_in, const double delta_nmuon_init)
{
  resize_group_vector(
    vec_delta_nmuon_group_upper, n_group_in, delta_nmuon_init);
}

// allocate memory for vec_volume_group
void Grid2dBinGroup::resize_vec_volume_group(
  const int n_group_in, const double volume_group_init)
{
  resize_group_vector(
    vec_volume_group, n_group_in, volume_group_init);
}

// allocate memory for vec_eff_low_group
void Grid2dBinGroup::resize_vec_eff_low_group(
  const int n_group_in, const double eff_low_init)
{
  resize_group_vector(
    vec_eff_low_group, n_group_in, eff_low_init);
}

// allocate memory for vec_eff_cnt_group
void Grid2dBinGroup::resize_vec_eff_cnt_group(
  const int n_group_in, const double eff_cnt_init)
{
  resize_group_vector(
    vec_eff_cnt_group, n_group_in, eff_cnt_init);
}

// allocate memory for vec_eff_upp_group
void Grid2dBinGroup::resize_vec_eff_upp_group(
  const int n_group_in, const double eff_upp_init)
{
  resize_group_vector(
    vec_eff_upp_group, n_group_in, eff_upp_init);
}
