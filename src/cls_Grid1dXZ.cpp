// cls_Grid1dXZ.cpp
#include "cls_Grid1dXZ.hpp"
#include "ns_myapp.hpp"
#include "ns_mymacro.hpp"
#include "ns_mylogger.hpp"
#include "spdlog_pch.hpp"

//########################################################
//########################################################
// class Grid1dXZ
// as derived class of the base class Grid2d
//########################################################
//########################################################

Grid1dXZ::Grid1dXZ( const std::filesystem::path &path_in, const bool tf_shift_in )
{
  build_from_ascii_xy(path_in,tf_shift_in);
}

void Grid1dXZ::build_from_ascii_xy(
    const std::filesystem::path &path_in, const bool tf_shift_in
  , const double interval_precision )
{
  myapp::filecheck(path_in);
  fprintf(stderr,"@ %s : fname_fullpath = %s\n",__PRETTY_FUNCTION__,path_in.string().c_str());
  LOG_DEBUG("@ {} : fname_fullpath = {}",__PRETTY_FUNCTION__,path_in.string());

  // read ascii file (x0,y0, x1,y1, x2,y2, ...)
  std::vector<std::pair<double,double>> vec_dpair = myapp::read_double2(path_in);
  mysort::sort_vec_double2_x_increasing(vec_dpair);
  bool is_same_interval_x = myapp::is_same_interval_x(vec_dpair,interval_precision);
  if(is_same_interval_x==false){
    THROW_ERROR("Grid1dXZ::build_from_ascii_xy: X values do not have uniform intervals. interval_precision={}", interval_precision);
  }
  
  // determine min max in Grid1d and zmin zmax in vec_z
  constexpr double big = std::numeric_limits<double>::max();
  double xmin_tmp =  big;
  double xmax_tmp = -big;
  double zmin_tmp =  big;
  double zmax_tmp = -big;
  for(const auto& [x,z] : vec_dpair){
    xmin_tmp = std::min(xmin_tmp,x);
    xmax_tmp = std::max(xmax_tmp,x);
    zmin_tmp = std::min(zmin_tmp,z);
    zmax_tmp = std::max(zmax_tmp,z);
  }


  set_zmin(zmin_tmp);
  set_zmax(zmax_tmp);

  LOG_DEBUG("xmin_tmp = {:.15E}",xmin_tmp);
  LOG_DEBUG("xmax_tmp = {:.15E}",xmax_tmp);

  // assign of Grid1d
  // stem() returns filename without extension
  const std::string name_in("xaxis_" + path_in.stem().string());
  const int nbin_in = vec_dpair.size();
  double min_in, max_in;
  const double interval_in = (xmax_tmp - xmin_tmp) / (double)(nbin_in-1);
  LOG_DEBUG("interval_in = {:.15E}",interval_in);
  if( tf_shift_in==true ){
    min_in = xmin_tmp - 0.5*interval_in;
    max_in = xmax_tmp + 0.5*interval_in;
  }
  else{
    min_in = xmin_tmp;
    max_in = xmax_tmp + 1.0*interval_in;
  }

  LOG_DEBUG("min_in = {:.15E}",min_in);
  LOG_DEBUG("max_in = {:.15E}",max_in);

  assign(name_in,nbin_in,min_in,max_in,interval_in, interval_precision);
  out_info(spdlog::level::debug);

  // memory allocate
  vec_memory_allocate();

  // assign of vec_z
  fprintf(stderr,"@%s , assigning z values ...\n",__FUNCTION__);
  LOG_DEBUG("@{} , assigning z values ...",__FUNCTION__);
  for(const auto& [x,z] : vec_dpair) set_z(x,z);
  
}

double Grid1dXZ::get_z( const int index ) const
{
  check_inside(index); // if out of range, throw error
  const double z = vec_z.at(index);
  return z;
}

/// @brief get z value from x value
double Grid1dXZ::get_z( const double value ) const
{
  const int index = get_index(value);
  return get_z(index);
}

void Grid1dXZ::set_z( const int index, const double z_in )
{
  check_inside(index);
  vec_z.at(index) = z_in;
}

void Grid1dXZ::set_z( const double value, const double z_in )
{
  const int index = get_index(value);
  check_inside(index);
  vec_z.at(index) = z_in;
}

void Grid1dXZ::multiply_z( const double factor )
{
  for( auto &z : vec_z ) z *= factor;
}

void Grid1dXZ::add_z( const double z_in )
{
  for( auto &z : vec_z ) z += z_in;
}

// allocate memory for std::vector<double> vec_z;
void Grid1dXZ::vec_memory_allocate(){
  vec_z.resize( Grid1d::get_nbin() );
}

// output x z array
void Grid1dXZ::out(const std::filesystem::path& pathout) const
{
  if(pathout.empty()) THROW_ERROR("Grid1dXZ::out: Output file path is empty");
  FILE *fout = fopen(pathout.c_str(),"wt");
  if(fout == NULL) THROW_ERROR("Grid1dXZ::out: Cannot open file for writing. path={}", pathout.string());
  double x,z;
  for(int index=0;index<get_nbin();index++){
    x = get_center_value(index);
    z = get_z(index);
    fprintf(fout,"%E %E\n",x,z);
  }
  if(fclose(fout) == EOF) THROW_ERROR("Grid1dXZ::out: Failed to close output file. path={}", pathout.string());
}

// get z_value using bilinear interpolation
double Grid1dXZ::get_linear_interpolated_z( const double value ) const
{
  if( value <  get_min() ) THROW_ERROR("Grid1dXZ::get_linear_interpolated_z: Input value is below grid minimum. value={}, min={}", value, get_min());
  if( value >= get_max() ) THROW_ERROR("Grid1dXZ::get_linear_interpolated_z: Input value is at or above grid maximum. value={}, max={}", value, get_max());
  const double value0 = get_lower_value(value);
  const double value1 = get_upper_value(value);
  const double z0 = get_z(value0);
  const double z1 = get_z(value1);
  return myapp::linear_interpolation(value0,z0,value1,z1,value);
}

// get x_value from z_value
// Note: Multiple solutions are possible depending on the function shape
std::vector<double> Grid1dXZ::get_linear_interpolated_x( const double z_value ) const
{
  std::vector<double> vec_x_value; // result container
  if( z_value < zmin ) return vec_x_value;
  if( z_value > zmax ) return vec_x_value;
  double x,x0,x1,z0,z1;
  for(int i=0;i<get_nbin()-1;i++){ // loop over x-axis points
    x0 = get_lower_value(i+0);
    x1 = get_lower_value(i+1);
    z0 = get_z(i+0);
    z1 = get_z(i+1);
    if( z0 < z1 ){ // if graph is increasing
      if( z_value < z0 ) continue; // skip if z_value is outside range
      if( z_value > z1 ) continue; // skip if z_value is outside range
      // z0 <= z_value <= z1
      x = myapp::linear_interpolation(z0,x0,z1,x1,z_value);
      vec_x_value.push_back(x);
    }else{ // z0 >= z1, if graph is decreasing
      if( z_value > z0 ) continue; // skip if z_value is outside range
      if( z_value < z1 ) continue; // skip if z_value is outside range
      // z1 <= z_value <= z0
      x = myapp::linear_interpolation(z0,x0,z1,x1,z_value);
      vec_x_value.push_back(x);
    }
  }
  return vec_x_value;
}

// make another Grid1dXZ object which is differential of this object
Grid1dXZ Grid1dXZ::make_differential() const
{
  Grid1dXZ g1d_diff(*this);
  double x0,x1,y0,y1,dx,dy,dydx;
  const int nbinx = get_nbin();
  for(int i=0;i<nbinx-1;i++){
    x0 = get_lower_value(i+0);
    x1 = get_lower_value(i+1);
    dx = x1 - x0;
    y0 = get_z(i+0);
    y1 = get_z(i+1);
    dy = y1 - y0;
    dydx = dy/dx;
    g1d_diff.set_z(i,dydx);
  }
  // Set the last point to the same value as the previous point
  g1d_diff.set_z(x1,dydx);
  g1d_diff.out("g1d_diff.tmp");
  return g1d_diff;
}

// Create vector of (x, z) pairs from Grid1dXZ
std::vector<std::pair<double,double>> Grid1dXZ::make_vec_double2_xz() const
{
  std::vector<std::pair<double,double>> vec_double2_xz;
  const int nbinx = get_nbin();
  for(int i=0;i<nbinx;i++){
    const double x = get_lower_value(i);
    const double z = get_z(i);
    vec_double2_xz.push_back(std::make_pair(x,z));
  }
  return vec_double2_xz;
}

// Create vector of (10^x, 10^z) pairs from Grid1dXZ
std::vector<std::pair<double,double>> Grid1dXZ::make_vec_double2_pow10_xz() const
{
  std::vector<std::pair<double,double>> vec_double2_pow10_xz;
  const int nbinx = get_nbin();
  for(int i=0;i<nbinx;i++){
    const double x = get_lower_value(i);
    const double z = get_z(i);
    const double pow10_x = pow(10,x);
    const double pow10_z = pow(10,z);
    vec_double2_pow10_xz.push_back(std::make_pair(pow10_x,pow10_z));
  }
  return vec_double2_pow10_xz;
}

// get the tuple of zmin and zmax
std::tuple<double,double> Grid1dXZ::get_zmin_zmax() const
{
  double zmin =  std::numeric_limits<double>::max();
  double zmax = -std::numeric_limits<double>::max();
  for(int i=0;i<get_nbin();i++){
    double z = get_z(i);
    if( z < zmin ) zmin = z;
    if( z > zmax ) zmax = z;
  }
  return std::make_tuple(zmin,zmax);
}

// get the tuple of zmin, x_zmin, zmax, x_zmax.
std::array<double,4> Grid1dXZ::get_zmin_x_zmin_zmax_x_zmax() const
{
  constexpr double big = std::numeric_limits<double>::max();
  double zmin =  big;
  double zmax = -big;
  double x_zmin=-big;
  double x_zmax= big;
  for(int i=0;i<get_nbin();i++){
    double z = get_z(i);
    if( z < zmin ){
      zmin = z;
      x_zmin = get_lower_value(i);
    }
    if( z > zmax ){
      zmax = z;
      x_zmax = get_lower_value(i);
    }
  }
  if( zmin== big ) THROW_ERROR("Grid1dXZ::get_zmin_x_zmin_zmax_x_zmax: Failed to find valid zmin (no bins processed)");
  if( zmax==-big ) THROW_ERROR("Grid1dXZ::get_zmin_x_zmin_zmax_x_zmax: Failed to find valid zmax (no bins processed)");
  return {zmin,x_zmin,zmax,x_zmax};
}
