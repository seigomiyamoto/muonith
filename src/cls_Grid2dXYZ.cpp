/// @file cls_Grid2dXYZ.cpp
/// @brief Implementation of Grid2dXYZ class
///
/// @details This file implements the Grid2dXYZ class methods including:
///          - ASCII/binary I/O
///          - Bilinear interpolation
///          - Grid transformations (log/pow10, derivatives)
///          - dF/dR table generation for muon flux calculations
///
///
/// @note Thread-safety: Read-only methods are thread-safe after construction.
///       Mutating methods require external synchronization.
#include "ns_myapp.hpp"
#include "cls_Grid2dXYZ.hpp"
#include "cls_XYArray.hpp"
#include "ns_mymacro.hpp"
#include "ns_mylogger.hpp"
#include "ns_param_constants.hpp"
#include "spdlog_pch.hpp"
#include <cstring>
#include <filesystem>

bool Grid2dXYZ::operator!=(const Grid2dXYZ& other) const
{
  #ifdef NODEBUG
    if (Grid2d::operator!=(other)) return true;
    if (vec_vec_z != other.vec_vec_z) return true;
  #else
    if (Grid2d::operator!=(other)){ LOG_WARN("Grid2dXYZ: Base Grid2d part differs"); return true; }
    if (vec_vec_z != other.vec_vec_z) {LOG_WARN("Grid2dXYZ: vec_vec_z differs"); return true; }
  #endif

  return false;
}

Grid2dXYZ::Grid2dXYZ( const fs::path &path_in
  , const bool tf_shift_x, const bool tf_shift_y
  , const double tolerance_ratio )
: Grid2dXYZ()
{
  if (path_in.extension() == ".g2zbin") {
    // g2zbin v2 stores canonical (already shifted) axes; tf_shift flags are
    // only meaningful for the ASCII path where raw points need interpretation.
    load_g2zbin(path_in, tolerance_ratio);
  } else {
    build_from_ascii_xyz(path_in, tf_shift_x, tf_shift_y, tolerance_ratio);
  }
}

void Grid2dXYZ::vec_vec_memory_allocate(const double z_ini)
{
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  vec_vec_z.resize(nbiny); // allocate y direction
  for(int iy=0;iy<nbiny;iy++){
    vec_vec_z.at(iy).resize(nbinx); // allocate x direction
    for(int ix=0;ix<nbinx;ix++){
      set_z(ix,iy,z_ini);
    }
  }
}

void Grid2dXYZ::build_from_ascii_xyz(
    const fs::path &path_in
  , const bool tf_shift_x, const bool tf_shift_y
  , const double tolerance_ratio )
{
  // Step 1: Validate input file path.
  myapp::filecheck(path_in);

  // Step 2: Set instance name from the file name.
  LOG_INFO("setting name of Grid2dXYZ from path_in: {}", path_in.string());
  set_name(path_in.filename().string());

  // Step 3: Read XYZ tuples.
  std::vector< std::array<double,3> > vec_xyz = myapp::read_vec_xyz(path_in);

  // Step 4: Sort by (y, x, z) order to stabilize axis detection.
  std::chrono::system_clock::time_point start, end;
  start = time_now;
  std::sort(vec_xyz.begin(),vec_xyz.end(),myapp::compare_array_yxz);
  end = time_now;
  myapp::cast_time_msec(stdout,"std::sort",start,end);

  // Step 5: Build x/y axes from sorted XYZ.
  LOG_INFO("...");
  set_xy_axis_from_vec_xyz(
    vec_xyz,tf_shift_x,tf_shift_y,tolerance_ratio);
  LOG_INFO("set_xy_axis_from_vec_xyz ... done.");

  // axis name assign
  const std::string name_x = get_x_axis().get_name() + "_" + path_in.stem().string();
  const std::string name_y = get_y_axis().get_name() + "_" + path_in.stem().string();
  set_x_axis_name(name_x);
  set_y_axis_name(name_y);

  // Step 6: Allocate storage and assign z values.
  vec_vec_memory_allocate();
  LOG_INFO("vec_vec_memory_allocate ... done.");

  // assignment of z values
  for(const auto& [x,y,z] : vec_xyz ){
    const int ix = get_ix(x);
    check_ix_inside(ix);
    const int iy = get_iy(y);
    check_iy_inside(iy);
    vec_vec_z.at(iy).at(ix) = z;
  }
  LOG_INFO("reading {} done.", path_in.string());
  get_x_axis().out_info(spdlog::level::debug);
  get_y_axis().out_info(spdlog::level::debug);
  LOG_INFO(", z_value assignment done.");
}

void Grid2dXYZ::out(const fs::path& pathout) const
{
  FILE *fout = myapp::get_fout(pathout);
  double x,y,z;
  for(int iy=0;iy<get_nbiny();iy++){
    y = get_ycnt(iy);
    for(int ix=0;ix<get_nbinx();ix++){
      x = get_xcnt(ix);
      z = get_z(ix,iy);
      fprintf(fout,"%E %E %E\n",x,y,z);
    }
  }
  myapp::close(fout,pathout);
}

void Grid2dXYZ::out_yx(const fs::path& pathout
  , const bool tf_xcnt, const bool tf_ycnt) const
{
  FILE *fout = myapp::get_fout(pathout);
  if (!fout) {
    LOG_ERROR("Failed to open file for writing: {}", pathout.string());
    throw std::runtime_error("Failed to open file: " + pathout.string());
  }

  for (int iy = 0; iy < get_nbiny(); ++iy) {
    double y = tf_ycnt ? get_ycnt(iy) : get_ylow(iy);
    for (int ix = 0; ix < get_nbinx(); ++ix) {
      double x = tf_xcnt ? get_xcnt(ix) : get_xlow(ix);
      double z = get_z(ix, iy);
      if (std::fprintf(fout, "%E %E %E\n", x, y, z) < 0) {
        LOG_ERROR("Failed to write to text file: {}", pathout.string());
        myapp::close(fout, pathout);
        throw std::runtime_error("Write error on: " + pathout.string());
      }
    }
  }
  myapp::close(fout, pathout);
}

void Grid2dXYZ::out_xy(const fs::path& pathout
  , const bool tf_xcnt, const bool tf_ycnt) const
{
  FILE *fout = myapp::get_fout(pathout);
  if (!fout) {
    LOG_ERROR("Failed to open file for writing: {}", pathout.string());
    throw std::runtime_error("Failed to open file: " + pathout.string());
  }

  for (int ix = 0; ix < get_nbinx(); ++ix) {
    double x = tf_xcnt ? get_xcnt(ix) : get_xlow(ix);
    for (int iy = 0; iy < get_nbiny(); ++iy) {
      double y = tf_ycnt ? get_ycnt(iy) : get_ylow(iy);
      double z = get_z(ix, iy);
      if (std::fprintf(fout, "%E %E %E\n", x, y, z) < 0) {
        LOG_ERROR("Failed to write to text file: {}", pathout.string());
        myapp::close(fout, pathout);
        throw std::runtime_error("Write error on: " + pathout.string());
      }
    }
  }
  myapp::close(fout, pathout);
}

void Grid2dXYZ::out_all_xy_order_ylower(
  const fs::path& pathout) const
{
  FILE *fout = myapp::get_fout(pathout);
  double x,y,z;
  for(int ix=0;ix<get_nbinx();ix++){
    x = get_xcnt(ix);
    for(int iy=0;iy<get_nbiny();iy++){
      y = get_ylow(iy);
      z = get_z(ix,iy);
      fprintf(fout,"%E %E %E\n",x,y,z);
    }
  }
  myapp::close(fout,pathout);
}


double Grid2dXYZ::get_z( const int ix, const int iy ) const {
  check_ix_inside(ix);
  check_iy_inside(iy);
  const double z_value = vec_vec_z.at(iy).at(ix);
  return z_value;
}

double Grid2dXYZ::get_z( const double x, const double y ) const {
  const int ix = get_ix(x);
  const int iy = get_iy(y);
  check_ix_inside(ix);
  check_iy_inside(iy);
  return get_z(ix,iy);
}

double Grid2dXYZ::get_bilinear_interpolated_z_value(const double x, const double y ) const 
{
  // Step 1: Validate input bounds.
  if( Grid2d::is_inside(x,y)==false ){
    LOG_ERROR("x,y is out of range.");
    LOG_ERROR("x = {}, y = {}",x,y);
    LOG_ERROR("xmin = {}, xmax = {}",get_xmin(),get_xmax());
    LOG_ERROR("ymin = {}, ymax = {}",get_ymin(),get_ymax());
    throw std::runtime_error("Grid2dXYZ::get_bilinear_interpolated_z_value : x,y is out of range.");
  }

  // Step 2: Clamp to avoid edge overflow on upper bounds.
  const double eps_xy = 1.0e-9;
  const double x_clamped = std::clamp(x, get_xmin(), get_xmax() - get_x_interval()*eps_xy);
  const double y_clamped = std::clamp(y, get_ymin(), get_ymax() - get_y_interval()*eps_xy);

  // Step 3: Locate the surrounding bin indices.
  const int nbinx = get_nbinx();
  const int ix1 = get_ix(x_clamped);
  check_ix_inside(ix1);
  const int ix2 = ix1 + 1;
  
  const int nbiny = get_nbiny();
  const int iy1 = get_iy(y_clamped);
  check_iy_inside(iy1);
  const int iy2 = iy1 + 1;

  double x1,x2,y1,y2;
  double z11,z12,z21,z22;

  // Step 4: Handle boundary cases with linear interpolation.
  // Handle ix == nbinx - 1 boundary.
  if( ix1==nbinx-1 ){
    if( iy1==nbiny-1 ){
      return get_z(ix1,iy1);
    }else{
      y1 = get_ylow(iy1);
      y2 = get_yup(iy1);
      z11 = get_z(ix1,iy1);
      z12 = get_z(ix1,iy2);
      return myapp::linear_interpolation(y1,z11,y2,z12,y_clamped);
    }
  }
  // Handle iy == nbiny - 1 boundary.
  if( iy1==nbiny-1 ){
    if( ix1==nbinx-1 ){
      return get_z(ix1,iy1);
    }else{
      x1 = get_xlow(ix1);
      x2 = get_xup(ix1);
      z11 = get_z(ix1,iy1);
      z21 = get_z(ix2,iy1);
      return myapp::linear_interpolation(x1,z11,x2,z21,x_clamped);
    }
  }
  // Step 5: Perform bilinear interpolation inside the grid cell.
  x1 = get_xlow(ix1);
  x2 = get_xup(ix1);
  y1 = get_ylow(iy1);
  y2 = get_yup(iy1);
  z11 = get_z(ix1,iy1);
  z21 = get_z(ix2,iy1);
  z12 = get_z(ix1,iy2);
  z22 = get_z(ix2,iy2);
  return myapp::bilinear_interpolation(x1,x2,y1,y2,z11,z21,z12,z22,x_clamped,y_clamped);
}

Grid1dXZ Grid2dXYZ::get_Grid1dXZ_when_x_fixed( const double x_in, const bool tf_xcnt ) const
{
  // Find the nearest bin center (or lower) to x_in.
  const Grid1d x_axis = get_x_axis();
  const double x_fixed = x_axis.get_nearest_value(x_in, tf_xcnt);
  const int ix_fix = x_axis.get_index(x_fixed);
  LOG_DEBUG("x_in={:.6E}, x_fixed={:.6E}, ix_fix={}", x_in, x_fixed, ix_fix);

  // Set y axis for Grid1dXZ from this grid.
  const Grid1d y_axis = get_y_axis();
  Grid1dXZ g1z_yz;
  g1z_yz.set_value(y_axis);
  const std::string name_g1z_yz = name + "_x=" + std::to_string(x_in);
  g1z_yz.set_name(name_g1z_yz);
  g1z_yz.vec_memory_allocate();

  // Assign values.
  const int nbiny = get_nbiny();
  for(int iy=0;iy<nbiny;iy++){
    double z_value = get_z(ix_fix,iy);
    g1z_yz.set_z(iy,z_value);
  }
  return g1z_yz;
}

Grid1dXZ Grid2dXYZ::get_Grid1dXZ_interp(const double x_in) const
{
  const Grid1d x_axis = get_x_axis();
  const Grid1d y_axis = get_y_axis();

  const double xmin = x_axis.get_min();
  const double xmax = x_axis.get_max();
  if (x_in < xmin || x_in >= xmax) {
    LOG_ERROR("x_in = {} is out of range. xmin = {}, xmax = {}", x_in, xmin, xmax);
    throw std::runtime_error("x_in is out of range.");
  }

  Grid1dXZ g1z_yz;
  g1z_yz.set_value(y_axis);
  g1z_yz.set_name(name + "_x=" + std::to_string(x_in));
  g1z_yz.vec_memory_allocate();

  const int nbiny = get_nbiny();
  for (int iy = 0; iy < nbiny; ++iy) {
    const double y = y_axis.get_lower_value(iy);
    const double z = get_bilinear_interpolated_z_value(x_in, y);
    g1z_yz.set_z(iy, z);
  }

  return g1z_yz;
}

Grid2dXYZ Grid2dXYZ::make_dFdR_table( const Grid1dXZ  &g1_logR_logE
                    , const Grid2dXYZ &g2_logdFdE_logE )
{
  // Step 1: Convert range table to linear scale and compute dR/dE.
  // Step 2: Determine output axis ranges (logE, costhz).
  // Step 3: Allocate output grid.
  // Step 4: For each (logE, costhz), interpolate dF/dE and compute dF/dR.
  LOG_INFO("Creating vector of (x, z) from Grid1dXZ.");
  std::vector<std::pair<double,double>> vec_logR_logE
    = g1_logR_logE.make_vec_double2_xz();
  
  // XYArray holds a set of (x, y) pairs.
  XYArray xyarr_logR_logE(vec_logR_logE);
  xyarr_logR_logE.set_name("xyarr_logR_logE");
  xyarr_logR_logE.out();
  XYArray xyarr_R_E = xyarr_logR_logE.get_converted_xy_to_pow10();
  xyarr_R_E.set_name("xyarr_R_E");
  xyarr_R_E.out();
  LOG_INFO("make_differential returns a vector of (x, dy/dx).");
  XYArray xyarr_dRdE_E = xyarr_R_E.make_differential();
  xyarr_dRdE_E.set_name("xyarr_dRdE_E");
  xyarr_dRdE_E.out();

  Grid2dXYZ g2_dFdR;

  // Find min/max logE from the differential flux table.
  const double logEmin_logdFdE = g2_logdFdE_logE.get_ymin();
  const double logEmax_logdFdE = g2_logdFdE_logE.get_ymax();
  LOG_INFO("logEmin_logdFdE={:.4E} (logGeV)",logEmin_logdFdE);
  LOG_INFO("logEmax_logdFdE={:.4E} (logGeV)",logEmax_logdFdE);
  
  // Find min/max logE from the range table.
  const double logEmin_logR = xyarr_logR_logE.get_xmin(); // GeV
  const double logEmax_logR = xyarr_logR_logE.get_xmax(); // GeV
  LOG_INFO("logEmin_logR={:.4E} (logGeV)",logEmin_logR);
  LOG_INFO("logEmax_logR={:.4E} (logGeV)",logEmax_logR);

  // Determine y-min for the output table as the larger of the two minima.
  const double logEmin_dFdR
     = (std::max)(logEmin_logdFdE,logEmin_logR);
  LOG_INFO("logEmin_dFdR={:.4E} (logGeV)",logEmin_dFdR);

  // Determine a provisional y-max for the output table (smaller of maxima).
  const double logEmax_tmp
     = (std::min)(logEmax_logdFdE,logEmax_logR);
  LOG_INFO("logEmax_tmp={:.4E} (logGeV)",logEmax_tmp);

  // y axis corresponds to logE, so use y-interval from g2_logdFdE_logE.
  const double logE_interval_dFdR = g2_logdFdE_logE.get_y_interval();
  LOG_INFO("logE_interval_dFdR={:.4E} (logGeV)",logE_interval_dFdR);

  // Compute remainder of (logEmax - logEmin) divided by interval.
  const double logE_remainder = fmod(logEmax_tmp-logEmin_dFdR,logE_interval_dFdR);
  LOG_INFO("logE_remainder=fmod(logEmax_tmp-logEmin_dFdR,logE_interval_dFdR)={:.4E} (logGeV)",logE_remainder);

  // Compute adjusted y-max aligned to the interval grid.
  const double logEmax_dFdR = logEmax_tmp - logE_remainder;
  LOG_INFO("logEmax_dFdR = logEmax_tmp - logE_remainder={:.4E} (logGeV)",logEmax_dFdR);
  
  // Compute nbiny for dFdR.
  const int nbiny_dFdR = (int)((logEmax_tmp-logEmin_dFdR)/logE_interval_dFdR);
  LOG_INFO("nbiny_dFdR=(int)((logEmax_tmp-logEmin_dFdR)/logE_interval_dFdR)={}",nbiny_dFdR);

  // Set x axis to costhz.
  g2_dFdR.set_x_axis(g2_logdFdE_logE.get_x_axis());

  // Set y axis to logE.
  const Grid1d y_axis("logE (y) axis of dFdR_table"
    ,nbiny_dFdR,logEmin_dFdR,logEmax_dFdR,logE_interval_dFdR);

  // Assign axes.
  g2_dFdR.set_y_axis(y_axis);
  g2_dFdR.set_name("g2_dFdR_costhz_logE");

  LOG_INFO("x/y_axis of {}",g2_dFdR.get_name());
  g2_dFdR.get_x_axis().out_info(spdlog::level::info);
  g2_dFdR.get_y_axis().out_info(spdlog::level::info);
  
  const int nbinx_dFdR = g2_dFdR.get_nbinx();

  // Allocate z grid.
  g2_dFdR.vec_vec_memory_allocate();

  for(int iy=0;iy<nbiny_dFdR;iy++){ // loop of logE
    // Evaluate logE at lower edge (x is costhz).
    const double logE = g2_dFdR.get_ylow(iy);
    const double E = pow(10,logE); //GeV

    const std::vector<double> vec_dRdE = xyarr_dRdE_E.get_interp_y(E);
    if( vec_dRdE.size()!=1 ) 
      THROW_ERROR("make_dFdR_table vec_dRdE.size()!=1");
    const double dRdE = vec_dRdE.at(0);

    for(int ix=0;ix<nbinx_dFdR;ix++){ // loop of costhz
      // Compute dF/dE at (costhz, logE).
      const double costhz = g2_dFdR.get_xcnt(ix);
      const double logdFdE = g2_logdFdE_logE.get_bilinear_interpolated_z_value(costhz,logE);
      const double dFdE = pow(10,logdFdE);
      mylogger::g_logger->trace("dFdE={:.4E} (m-2 sr-1 sec-1 GeV-1) at E={:.4E}(GeV), costhz={:.3f}", dFdE, E,costhz);
      mylogger::g_logger->trace("dRdE={:.4E} (kg m-2 GeV-1) at E={:.4E}(GeV), 1/dRdE={:.4E}",dRdE,E,1.0/dRdE);

      // Compute dF/dR.
      const double dFdR =  dFdE / dRdE;
      mylogger::g_logger->trace("dFdR={:.4E} (kg-1 sr-1 sec-1) at E={:.4E}(GeV), costhz={:.3f}", dFdR, E, costhz);
      
      // Assign value.
      g2_dFdR.set_z(ix,iy,dFdR);
    }
  }
  LOG_INFO("building log_dF/dR table done.");
  return g2_dFdR;
}

Grid2dXYZ Grid2dXYZ::make_dFdR_table_from_peneflux(
  const Grid2dXYZ &g2_log_peneflux_R_costhz)
{
  LOG_INFO("Create Grid2dXYZ dF/dR table from penetrating muon flux table.");
  
  // Create F, R, costhz table (linear scale).
  Grid2dXYZ g2_F_R_costhz = g2_log_peneflux_R_costhz.make_pow10_z();
  g2_F_R_costhz.set_name("g2_F_R_costhz");
  g2_F_R_costhz.get_x_axis().out_info(spdlog::level::info);
  g2_F_R_costhz.get_y_axis().out_info(spdlog::level::info);

  // Create dF/dR table.
  Grid2dXYZ g2_dFdR_costhz = g2_F_R_costhz.make_z_dzdy();

  return g2_dFdR_costhz;
}

Grid2dXYZ Grid2dXYZ::make_log_z() const
{
  // Copy Grid2dXYZ.
  Grid2dXYZ g2_logz(*this);
  g2_logz.initialize_z();
  const int nbinx = g2_logz.get_nbinx();
  const int nbiny = g2_logz.get_nbiny();
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      const double z = this->get_z(ix,iy);
      const double log_z = myapp::log10(z);
      g2_logz.set_z(ix,iy,log_z);
    }
  }
  return g2_logz;
}

Grid2dXYZ Grid2dXYZ::make_pow10_z() const
{
  // Copy Grid2dXYZ.
  Grid2dXYZ g2_pow10z(*this);
  g2_pow10z.initialize_z();
  const int nbinx = g2_pow10z.get_nbinx();
  const int nbiny = g2_pow10z.get_nbiny();
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      const double log_z = this->get_z(ix,iy);
      const double z = pow(10,log_z);
      g2_pow10z.set_z(ix,iy,z);
    }
  }
  return g2_pow10z;
}

Grid2dXYZ Grid2dXYZ::make_log_abs_z() const
{
  // Copy Grid2dXYZ.
  Grid2dXYZ g2_logz(*this);
  g2_logz.initialize_z();
  const int nbinx = g2_logz.get_nbinx();
  const int nbiny = g2_logz.get_nbiny();
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      const double z = fabs( this->get_z(ix,iy) );
      const double log_z = myapp::log10(z);
      g2_logz.set_z(ix,iy,log_z);
    }
  }
  return g2_logz;
}

void Grid2dXYZ::initialize_z( const double z_ini)
{
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      set_z(ix,iy,z_ini);
    }
  }
}

Grid2dXYZ Grid2dXYZ::make_z_dzdy() const
{
  // Copy Grid2dXYZ.
  Grid2dXYZ g2_dzdy(*this);
  g2_dzdy.initialize_z();
  const int nbinx = g2_dzdy.get_nbinx();
  const int nbiny = g2_dzdy.get_nbiny();
  const double dy = this->get_y_interval();
  for(int iy=0;iy<nbiny-1;iy++){
    for(int ix=0;ix<nbinx;ix++){
      const double z0 = this->get_z(ix,iy+0);
      const double z1 = this->get_z(ix,iy+1);
      const double dzdy = (z1-z0)/dy;
      g2_dzdy.set_z(ix,iy,dzdy);
    }
  }
  // The last row (iy=nbiny-1) reuses the previous derivative.
  const int iy = nbiny-1;
  for(int ix=0;ix<nbinx;ix++){
    const double z0 = this->get_z(ix,iy-1);
    const double z1 = this->get_z(ix,iy);
    const double dzdy = (z1-z0)/dy;
    g2_dzdy.set_z(ix,iy,dzdy);
  }
  return g2_dzdy;
}

double Grid2dXYZ::get_bilinear_interpolated_peneflux(
    const double costhz, const double DL_in ) const
{
  if (!std::isfinite(costhz) || !std::isfinite(DL_in)) {
    LOG_ERROR("non-finite input costhz={:.4f}, DL={:.1f}", costhz, DL_in);
  }
  // Effective DL clamp bounds: explicit set_DL_clamp_bounds() values win;
  // otherwise fall back to this grid's own y-axis (table-derived defaults).
  // Both are kept inside the readable axis range (y-axis max is exclusive).
  const double DL_floor = std::max(
    std::isnan(d_DL_clamp_min) ? get_ymin() : d_DL_clamp_min, get_ymin());
  const double DL_ceil = std::min(
    std::isnan(d_DL_clamp_max) ? get_ymax() : d_DL_clamp_max,
    get_ymax() * (1.0 - 1.0e-8));
  if (DL_in < 0.0) {
    LOG_WARN("DL_in is negative ({}), use DL_min={:.3E}", DL_in, DL_floor);
  }
  if( costhz < param_constants::costhz_min() ) THROW_ERROR_NAME3("costhz < costhz_min",costhz, param_constants::costhz_min());
  if( costhz > param_constants::costhz_max() ) THROW_ERROR_NAME3("costhz > costhz_max",costhz, param_constants::costhz_max());
  double DL = DL_in;
  if( DL_in > DL_ceil ){
    LOG_TRACE_ND("get_bilinear_interpolated_peneflux: DL_in({:.3E}) > DL_max({:.3E}), use DL_max", DL_in, DL_ceil);
    DL = DL_ceil;
  }
  if( DL_in <= 0.0 ){
    DL = DL_floor;
  }else if( DL_in < DL_floor ){
    LOG_TRACE_ND("get_bilinear_interpolated_peneflux: DL_in({:.3E}) < DL_min({:.3E}), use DL_min", DL_in, DL_floor);
    DL = DL_floor;
  }
  const double log10_peneflux = Grid2dXYZ::get_bilinear_interpolated_z_value(costhz,DL);
  const double peneflux = pow(10,log10_peneflux);
  return peneflux;
}

std::vector<std::array<double, 3>>
  Grid2dXYZ::make_vec_double3_xyz() const
{
  std::vector<std::array<double, 3>> vec;
  vec.reserve(get_nbinx() * get_nbiny());

  for (int ix = 0; ix < get_nbinx(); ++ix) {
    double x = get_x_axis().get_lower_value(ix);
    for (int iy = 0; iy < get_nbiny(); ++iy) {
      double y = get_y_axis().get_lower_value(iy);
      double z = get_z(ix, iy);
      vec.push_back({x, y, z});
    }
  }
  return vec;
}

//######################################################################
// g2zbin I/O
//######################################################################

void Grid2dXYZ::save_g2zbin(std::ofstream& ofs, bool use_float64) const
{
  // -- Fixed header (16 bytes) --
  constexpr char magic[8] = {'G','2','Z','B','I','N','\0','\0'};
  ofs.write(magic, 8);

  // Version 2: axes hold canonical (already shifted) min/max, written
  // pass-through by Grid1d::save().  Version 1 stored raw (pre-shift) values.
  const uint16_t version = 2;
  io_binary::write_binary(ofs, version);

  const uint8_t precision = use_float64 ? 8 : 4;
  io_binary::write_binary(ofs, precision);

  const uint8_t order = 0; // row-major (iy outer, ix inner)
  io_binary::write_binary(ofs, order);

  const uint32_t reserved = 0;
  io_binary::write_binary(ofs, reserved);

  // -- Grid2d header (name + x_axis + y_axis in raw coordinates) --
  io_binary::write_string(ofs, get_name());
  get_x_axis().save(ofs);
  get_y_axis().save(ofs);

  // -- z raster (row-major: iy outer, ix inner) --
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();

  if (use_float64) {
    for (int iy = 0; iy < nbiny; ++iy) {
      ofs.write(reinterpret_cast<const char*>(vec_vec_z.at(iy).data()),
                nbinx * static_cast<std::streamsize>(sizeof(double)));
    }
  } else {
    std::vector<float> buf(nbinx);
    for (int iy = 0; iy < nbiny; ++iy) {
      for (int ix = 0; ix < nbinx; ++ix) {
        buf[ix] = static_cast<float>(vec_vec_z.at(iy).at(ix));
      }
      ofs.write(reinterpret_cast<const char*>(buf.data()),
                nbinx * static_cast<std::streamsize>(sizeof(float)));
    }
  }

  if (ofs.fail()) {
    THROW_ERROR("Grid2dXYZ::save_g2zbin: stream write failed.");
  }
}

void Grid2dXYZ::save_g2zbin(const fs::path& pathout, bool use_float64) const
{
  std::ofstream ofs = io_binary::open_ofstream(pathout);
  save_g2zbin(ofs, use_float64);
  ofs.close();
  if (ofs.fail()) {
    THROW_ERROR("Grid2dXYZ::save_g2zbin: failed to close file. path={}",
                pathout.string());
  }
}

void Grid2dXYZ::load_g2zbin(const fs::path& path_in,
                              double tolerance_ratio)
{
  std::ifstream ifs = io_binary::open_ifstream(path_in);

  // -- Fixed header (16 bytes) --
  char magic[8];
  ifs.read(magic, 8);
  if (ifs.fail()) {
    THROW_ERROR("Grid2dXYZ::load_g2zbin: failed to read magic. path={}",
                path_in.string());
  }
  constexpr char expected[8] = {'G','2','Z','B','I','N','\0','\0'};
  if (std::memcmp(magic, expected, 8) != 0) {
    THROW_ERROR("Grid2dXYZ::load_g2zbin: invalid magic bytes. path={}",
                path_in.string());
  }

  const uint16_t version = io_binary::read_binary<uint16_t>(ifs);
  if (version != 2) {
    THROW_ERROR("Grid2dXYZ::load_g2zbin: unsupported version={} (expected 2). "
                "Version 1 files store pre-shift coordinates and must be "
                "regenerated from their source data. path={}",
                version, path_in.string());
  }

  const uint8_t precision = io_binary::read_binary<uint8_t>(ifs);
  if (precision != 8 && precision != 4) {
    THROW_ERROR("Grid2dXYZ::load_g2zbin: invalid precision={}. path={}",
                precision, path_in.string());
  }

  const uint8_t order = io_binary::read_binary<uint8_t>(ifs);
  if (order != 0) {
    THROW_ERROR("Grid2dXYZ::load_g2zbin: unsupported order={}. path={}",
                order, path_in.string());
  }

  // reserved 4 bytes -- read and discard
  io_binary::read_binary<uint32_t>(ifs);

  // -- Grid2d header --
  const std::string grid_name = io_binary::read_string(ifs);

  // Read canonical axes (pass-through via Grid1d::load()).
  Grid1d x_axis_new;
  x_axis_new.load(ifs, tolerance_ratio);
  Grid1d y_axis_new;
  y_axis_new.load(ifs, tolerance_ratio);

  set_x_axis(x_axis_new);
  set_y_axis(y_axis_new);

  // Set Grid2d and Grid2dXYZ names
  set_name(grid_name);

  // -- z raster --
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  vec_vec_memory_allocate();

  if (precision == 4) {
    std::vector<float> buf(nbinx);
    for (int iy = 0; iy < nbiny; ++iy) {
      ifs.read(reinterpret_cast<char*>(buf.data()),
               nbinx * static_cast<std::streamsize>(sizeof(float)));
      if (ifs.fail()) {
        THROW_ERROR("Grid2dXYZ::load_g2zbin: failed to read z row iy={}. path={}",
                    iy, path_in.string());
      }
      for (int ix = 0; ix < nbinx; ++ix) {
        vec_vec_z.at(iy).at(ix) = static_cast<double>(buf[ix]);
      }
    }
  } else {
    for (int iy = 0; iy < nbiny; ++iy) {
      ifs.read(reinterpret_cast<char*>(vec_vec_z.at(iy).data()),
               nbinx * static_cast<std::streamsize>(sizeof(double)));
      if (ifs.fail()) {
        THROW_ERROR("Grid2dXYZ::load_g2zbin: failed to read z row iy={}. path={}",
                    iy, path_in.string());
      }
    }
  }

  LOG_INFO("Grid2dXYZ::load_g2zbin: loaded {}. nbinx={}, nbiny={}, precision={}",
           path_in.string(), nbinx, nbiny, precision);
  get_x_axis().out_info(spdlog::level::debug);
  get_y_axis().out_info(spdlog::level::debug);
}

