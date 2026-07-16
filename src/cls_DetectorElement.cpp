/// @file cls_DetectorElement.cpp
/// @brief Implementation of DetectorElement class methods
/// @details
/// Implements calculation routines, I/O operations, and validation logic for DetectorElement.
/// Key features:
/// - Solid angle calculation with multiple angle unit support
/// - Effective area computation accounting for angular divergence
/// - Muon flux interpolation and signal calculations
/// - Projected density reconstruction (DL/PL ratio)
/// - Binary serialization for save/load operations
/// - Diagnostic inequality operator with optional logging
#include <cassert>
#include <iostream>

#include "ns_myapp.hpp"
#include "ns_angle_util.hpp"
#include "cls_DetectorElement.hpp"
#include "ns_mymacro.hpp"
#include "ns_mylogger.hpp"
#include "spdlog_pch.hpp"
#include "ns_stats_util.hpp"
//##################################################
//##################################################
// class DetectorElement
//##################################################
//##################################################

//==================================================
// operators
//==================================================

/// @brief Inequality operator
bool DetectorElement::operator!=(const DetectorElement& other) const
{
  #ifdef NODEBUG
    if (unique_index != other.unique_index) return true;
    if (detid != other.detid) return true;
    if (id_in_this_detector != other.id_in_this_detector) return true;
    if (ray3d != other.ray3d) return true;
    if (angle_unit != other.angle_unit) return true;
    if (txmin != other.txmin) return true;
    if( txmax != other.txmax) return true;
    if (tymin != other.tymin) return true;
    if (tymax != other.tymax) return true;
    if (effective_area_m2 != other.effective_area_m2) return true;
    if (solid_angle != other.solid_angle) return true;
    if (exposure_time_sec != other.exposure_time_sec) return true;
    if (PL != other.PL) return true;
    if (DL != other.DL) return true;
    if (penetrating_muon_flux != other.penetrating_muon_flux) return true;
    if (signal != other.signal) return true;
    if (noise_det != other.noise_det) return true;
    if (noise_poi != other.noise_poi) return true;
    if (proj_density != other.proj_density) return true;
    if (proj_density_upper != other.proj_density_upper) return true;
    if (proj_density_lower != other.proj_density_lower) return true;
    if (eff_low != other.eff_low) return true;
    if (eff_cnt != other.eff_cnt) return true;
    if (eff_upp != other.eff_upp) return true;
    if (vec_tf_in_PL != other.vec_tf_in_PL) return true;
  #else
    if (unique_index != other.unique_index) { LOG_WARN("DetectorElement: unique_index differs"); return true; }
    if (detid != other.detid) { LOG_WARN("DetectorElement: detid differs"); return true; }
    if (id_in_this_detector != other.id_in_this_detector) { LOG_WARN("DetectorElement: id_in_this_detector differs"); return true; }
    if (ray3d != other.ray3d) { LOG_WARN("DetectorElement: ray3d differs"); return true; }
    if (angle_unit != other.angle_unit) { LOG_WARN("DetectorElement: angle_unit differs"); return true; }
    if (txmin != other.txmin) { LOG_WARN("DetectorElement: txmin differs"); return true; }
    if (txmax != other.txmax) { LOG_WARN("DetectorElement: txmax differs"); return true; }
    if (tymin != other.tymin) { LOG_WARN("DetectorElement: tymin differs"); return true; }
    if (tymax != other.tymax) { LOG_WARN("DetectorElement: tymax differs"); return true; }
    if (effective_area_m2 != other.effective_area_m2) { LOG_WARN("DetectorElement: effective_area_m2 differs"); return true; }
    if (solid_angle != other.solid_angle) { LOG_WARN("DetectorElement: solid_angle differs"); return true; }
    if (exposure_time_sec != other.exposure_time_sec) { LOG_WARN("DetectorElement: exposure_time_sec differs"); return true; }
    if (PL != other.PL) { LOG_WARN("DetectorElement: PL differs"); return true; }
    if (DL != other.DL) { LOG_WARN("DetectorElement: DL differs"); return true; }
    if (penetrating_muon_flux != other.penetrating_muon_flux) { LOG_WARN("DetectorElement: penetrating_muon_flux differs"); return true; }
    if (signal != other.signal) { LOG_WARN("DetectorElement: signal differs"); return true; }
    if (noise_det != other.noise_det) { LOG_WARN("DetectorElement: noise_det differs"); return true; }
    if (noise_poi != other.noise_poi) { LOG_WARN("DetectorElement: noise_poi differs"); return true; }
    if (proj_density != other.proj_density) { LOG_WARN("DetectorElement: proj_density differs"); return true; }
    if (proj_density_upper != other.proj_density_upper) { LOG_WARN("DetectorElement: proj_density_upper differs"); return true; }
    if (proj_density_lower != other.proj_density_lower) { LOG_WARN("DetectorElement: proj_density_lower differs"); return true; }
    if (eff_low != other.eff_low) { LOG_WARN("DetectorElement: eff_low differs"); return true; }
    if (eff_cnt != other.eff_cnt) { LOG_WARN("DetectorElement: eff_cnt differs"); return true; }
    if (eff_upp != other.eff_upp) { LOG_WARN("DetectorElement: eff_upp differs"); return true; }
    if (vec_tf_in_PL != other.vec_tf_in_PL) { LOG_WARN("DetectorElement: vec_tf_in_PL differs"); return true; }
  #endif
  return false;
}

//----------------------------------------------------------

/// @brief Write ASCII output
void DetectorElement::out( FILE *fout ) const
{
  // detector index
  fprintf(fout,"%d "
  , unique_index);

  // v3_direction
  fprintf(fout,"%7.4lf %7.4lf %7.4lf "
  , ray3d.vx(), ray3d.vy(), ray3d.vz());

  // v3_position
  fprintf(fout,"%E %E %E "
  , ray3d.x(), ray3d.y(), ray3d.z());

  // v4_tangent (unit based output)
  switch(angle_unit) {
    case AngleUnit::Tangent:
      fprintf(fout, "Tangent: %7.4lf %7.4lf %7.4lf %7.4lf "
      , txmin, txmax, tymin, tymax);
      break;
    case AngleUnit::Radian:
      fprintf(fout, "Radian: %7.4lf %7.4lf %7.4lf %7.4lf "
      , txmin, txmax, tymin, tymax);
      break;
    case AngleUnit::Degree:
      fprintf(fout, "Degree: %7.4lf %7.4lf %7.4lf %7.4lf "
      , txmin, txmax, tymin, tymax);
      break;
  }
  
  // effective_area_m2, solid_angle, exposure time
  fprintf(fout,"%E %E %E "
  , effective_area_m2, solid_angle, exposure_time_sec);

  // PL, DL, penetrating_muon_flux, signal, noise (det floor + poi bucket = total)
  fprintf(fout,"%E %E %E %E %E "
  ,PL, DL, penetrating_muon_flux, signal, noise_det + noise_poi);
  // proj_density, proj_density_lower, proj_density_upper
  fprintf(fout,"%E %E %E "
  , proj_density, proj_density_lower, proj_density_upper
  );

  // eff_low, eff_cnt, eff_upp
  fprintf(fout,"%E %E %E "
  , eff_low, eff_cnt, eff_upp
  );

  fprintf(fout,"\n");
}

/// @brief Construct diagnostic name string
std::string DetectorElement::get_name() const
{
  std::string name =
      "DetectorElement_uqid=" + std::to_string(unique_index)
      + "det_id=" + std::to_string(detid)
      + "id_in_this_detector=" + std::to_string(id_in_this_detector);
  return name;
}

/// @brief Calculate solid angle based on angle_unit
double DetectorElement::calc_solid_angle() const
{
  if( angle_unit == AngleUnit::Degree ) {
    // call angle_util::calc_omega_degree
    const double txmin_rad = txmin * M_PI / 180.0;
    const double txmax_rad = txmax * M_PI / 180.0;
    const double tymin_rad = tymin * M_PI / 180.0;
    const double tymax_rad = tymax * M_PI / 180.0;
    return angle_util::calc_omega_radian(txmin_rad, txmax_rad, tymin_rad, tymax_rad);
  }
  if( angle_unit == AngleUnit::Radian ) {
    // call angle_util::calc_omega_radian
    return angle_util::calc_omega_radian(txmin, txmax, tymin, tymax);
  }
  if( angle_unit == AngleUnit::Tangent ) {
    // call angle_util::calc_omega_tangent
    return angle_util::calc_omega_tangent(txmin, txmax, tymin, tymax);
  }
  THROW_ERROR_NAME("DetectorElement::angle_unit is not defined");
}

/// @brief Calculate solid angle using alternative formula
double DetectorElement::calc_solid_angle_alternative() const
{
  if( angle_unit == AngleUnit::Degree ) {
    // call angle_util::calc_omega_degree
    const double txmin_rad = txmin * M_PI / 180.0;
    const double txmax_rad = txmax * M_PI / 180.0;
    const double tymin_rad = tymin * M_PI / 180.0;
    const double tymax_rad = tymax * M_PI / 180.0;
    return angle_util::calc_omega_radian(txmin_rad, txmax_rad, tymin_rad, tymax_rad);
  }
  if( angle_unit == AngleUnit::Radian ) {
    // call angle_util::calc_omega_radian
    return angle_util::calc_omega_radian(txmin, txmax, tymin, tymax);
  }
  if( angle_unit == AngleUnit::Tangent ) {
    return angle_util::calc_omega_tangent_alternative(txmin, txmax, tymin, tymax);
  } else {
    THROW_ERROR_NAME("angle_unit is not tangent");
  }
}

/// @brief Calculate effective area with angular corrections
double DetectorElement::calc_effective_area(
    const Eigen::Vector3d v3_ele_direction
  , const Eigen::Vector3d v3_det_length
  , const double n_unit ) const
{
  // check v3_ele_direction is not zero vector
  if( v3_ele_direction.norm() < 1.0e-6 ) THROW_ERROR_NAME("v3_ele_direction.norm() < 1.0e-6");

  // dot product of direction vectors
  const double dot_product = v3_ele_direction.dot(ray3d.dir());
  
  double tanx = get_tx();
  double tany = get_ty();
  // If angle_unit is Degree, convert to tangent: tanx=tan(tx), tany=tan(ty)
  if( angle_unit==AngleUnit::Degree ){
    tanx = tan(tanx*M_PI/180.0);
    tany = tan(tany*M_PI/180.0);
  }
  // If angle_unit is Radian, convert to tangent
  else if (angle_unit==AngleUnit::Radian) {
    tanx = tan(tanx);
    tany = tan(tany);
  }
  // If angle_unit is Tangent, use tanx, tany as-is
  else if (angle_unit==AngleUnit::Tangent) {
    // do nothing
  }
  else {
    THROW_ERROR_NAME("angle_unit is not defined");
  }

  // length of x and y dimension
  const double xlen = v3_det_length.x() - v3_det_length.z()*fabs(tanx);
  const double ylen = v3_det_length.y() - v3_det_length.z()*fabs(tany);
  return dot_product * xlen * ylen * n_unit; // m^2
}

// Check if position is inside bounding box
bool DetectorElement::is_inside(
   const double xmin, const double xmax
  ,const double ymin, const double ymax
  ,const double zmin, const double zmax ) const
{
  const double x = get_x();
  if( x <  xmin ) return false;
  if( x >= xmax ) return false;
  const double y = get_y();
  if( y <  ymin ) return false;
  if( y >= ymax ) return false;
  const double z = get_z();
  if( z <  zmin ) return false;
  if( z >= zmax ) return false;
  return true;
}

// Check if position is inside voxel
bool DetectorElement::is_inside( const Grid3dVoxel &g3vox ) const
{
  return is_inside(
      g3vox.get_xmin(), g3vox.get_xmax() 
    , g3vox.get_ymin(), g3vox.get_ymax() 
    , g3vox.get_zmin(), g3vox.get_zmax() 
  );
}

// Validate vec_tf_in_PL alternates and ends with false
void DetectorElement::check_alternating_tf_in() const
{
  bool tf_in_tmp_prev = false;
  for(const auto &[tf_in_tmp,dist_tmp] : this->get_vec_tf_in_PL() ){
    if( tf_in_tmp_prev==tf_in_tmp ){
      LOG_ERROR("tf_in should be flipped from false to true, or true to false.");
      LOG_DEBUG("det_id={:02d}, tx={:7.4f}, ty={:7.4f}",get_detid(),get_tx(),get_ty());
      THROW_ERROR("check_alternating_tf_in tf_in_prev==tf_in");
    }
    tf_in_tmp_prev = tf_in_tmp;
  }
  // the last tf_in should be false
  if( tf_in_tmp_prev ){
    LOG_ERROR("check_alternating_tf_in : the last tf_in should be false.");
    LOG_DEBUG("det_id={:02d}, tx={:7.4f}, ty={:7.4f}",get_detid(),get_tx(),get_ty());
    THROW_ERROR("calc_DL_with_vec_tf_in_PL the last tf_in==true");
  }
}

// Convert angular range to target angle unit
std::array<double,4>
DetectorElement::get_txmin_txmax_tymin_tymax(
  const AngleUnit& angle_unit_in) const
{
  double txmin_tmp = txmin;
  double txmax_tmp = txmax;
  double tymin_tmp = tymin;
  double tymax_tmp = tymax;
  if( angle_unit == AngleUnit::Degree ) {
    if( angle_unit_in == AngleUnit::Radian ) {
      txmin_tmp = txmin * M_PI / 180.0;
      txmax_tmp = txmax * M_PI / 180.0;
      tymin_tmp = tymin * M_PI / 180.0;
      tymax_tmp = tymax * M_PI / 180.0;
    }
    else if( angle_unit_in == AngleUnit::Tangent ) {
      txmin_tmp = tan(txmin * M_PI / 180.0);
      txmax_tmp = tan(txmax * M_PI / 180.0);
      tymin_tmp = tan(tymin * M_PI / 180.0);
      tymax_tmp = tan(tymax * M_PI / 180.0);
    }
  }
  else if( angle_unit == AngleUnit::Radian ) {
    if( angle_unit_in == AngleUnit::Degree ) {
      txmin_tmp = txmin * 180.0 / M_PI;
      txmax_tmp = txmax * 180.0 / M_PI;
      tymin_tmp = tymin * 180.0 / M_PI;
      tymax_tmp = tymax * 180.0 / M_PI;
    }
    else if( angle_unit_in == AngleUnit::Tangent ) {
      txmin_tmp = tan(txmin);
      txmax_tmp = tan(txmax);
      tymin_tmp = tan(tymin);
      tymax_tmp = tan(tymax);
    }
  }
  else if( angle_unit == AngleUnit::Tangent ) {
    if( angle_unit_in == AngleUnit::Degree ) {
      txmin_tmp = atan(txmin) * 180.0 / M_PI;
      txmax_tmp = atan(txmax) * 180.0 / M_PI;
      tymin_tmp = atan(tymin) * 180.0 / M_PI;
      tymax_tmp = atan(tymax) * 180.0 / M_PI;
    }
    else if( angle_unit_in == AngleUnit::Radian ) {
      txmin_tmp = atan(txmin);
      txmax_tmp = atan(txmax);
      tymin_tmp = atan(tymin);
      tymax_tmp = atan(tymax);
    }
  }
  else{
    THROW_ERROR("DetectorElement::get_txmin_txmax_tymin_tymax, angle_unit is not defined");
  }
  return {txmin_tmp, txmax_tmp, tymin_tmp, tymax_tmp};
}

double DetectorElement::calc_set_proj_density(
  const double PL_thres, const double DL_thres )
{
  if (PL < PL_thres){
    set_proj_density(invalid_proj_dens);
    return invalid_proj_dens;
  }
  if (DL < DL_thres){
    set_proj_density(invalid_proj_dens);
    return invalid_proj_dens;
  }
  const double proj_dens = DL / PL;
  set_proj_density(proj_dens);
  return proj_dens;
}

double DetectorElement::calc_peneflux( const Grid2dXYZ &g2flux) const
{
  const double costhz_tmp = this->get_ray3d().vz();
  const double DL_tmp = this->get_DL();
  const double peneflux_tmp
     = g2flux.get_bilinear_interpolated_peneflux(costhz_tmp,DL_tmp);
#ifdef DBGPRINT
mylogger::g_logger->trace("costhz_tmp={:.3f},DL_tmp={:E},peneflux_tmp={:E}",costhz_tmp,DL_tmp,peneflux_tmp);
#endif
  return peneflux_tmp;
}

double DetectorElement::calc_signal() const
{
  const double peneflux_tmp = this->get_peneflux();
  const double effective_area_m2_tmp = this->get_effective_area_m2();
  const double omega_tmp = this->get_solid_angle();
  const double time_sec_tmp = this->get_exposure_time_sec();
  const double signal_tmp
    = peneflux_tmp * effective_area_m2_tmp * omega_tmp * time_sec_tmp;
  return signal_tmp;
}

// calc signal with efficiency
double DetectorElement::calc_signal_eff() const
{
  const double signal = this->calc_signal();
  const double eff_sample = this->calc_eff_sample();
  const double signal_eff = signal * eff_sample;
  return signal_eff;
}

// calc approx. volume from vec_tf_in_PL as a cone volume
double DetectorElement::calc_approx_volume() const
{
  double volume = 0.0;
  const auto [txmin_tmp, txmax_tmp, tymin_tmp, tymax_tmp] 
    = get_txmin_txmax_tymin_tymax(AngleUnit::Radian);
  const double omega = this->calc_solid_angle();
  bool tf_in_prev = false;
  for( const auto [tf_in, length] : vec_tf_in_PL ){
    if( tf_in ){ // tf_in==true, means sky to rock
      volume -= 0.333333*length*length*length*omega;
    } else{      // tf_in==false, means rock to sky
      volume += 0.333333*length*length*length*omega;
    }
    tf_in_prev = tf_in;
  }
  if( tf_in_prev==true ){
    // the last tf_in should be false, means rock to sky
    THROW_ERROR("tf_in_last is true, it should be false");
  }
  // check volume is not NaN or negative
  if( std::isnan(volume) ){
    THROW_ERROR("calc_approx_volume: volume is {}",volume);
  }
  if( volume < 0.0 ){
    THROW_ERROR("calc_approx_volume: volume is {}, it should be non-negative",volume);
  }
  return volume;
}

// Sample efficiency from asymmetric Gaussian
double DetectorElement::calc_eff_sample() const
{
  return stats_util::sample_asymmetric_gaussian(eff_low,eff_cnt,eff_upp);
}


// Save all members to binary stream
void DetectorElement::save( std::ofstream &ofs ) const
{
  // unique_index
  io_binary::write_binary(ofs,unique_index);
  // detector id
  io_binary::write_binary(ofs,detid);
  // id_in_this_detector
  io_binary::write_binary(ofs,id_in_this_detector);
  // ray3d
  ray3d.save(ofs);
  // txmin, txmax, tymin, tymax
  io_binary::write_binary(ofs,txmin);
  io_binary::write_binary(ofs,txmax);
  io_binary::write_binary(ofs,tymin);
  io_binary::write_binary(ofs,tymax);
  
  // angle_unit, 0=Tangent, 1=Radian, 2=Degree
  io_binary::write_binary(ofs,static_cast<int>(angle_unit));

  // effective_area_m2, solid_angle, exposure_time_sec
  io_binary::write_binary(ofs,effective_area_m2);
  io_binary::write_binary(ofs,solid_angle);
  io_binary::write_binary(ofs,exposure_time_sec);
  // PL, DL, penetrating_muon_flux, signal, noise (det floor + poi bucket, serialized separately)
  io_binary::write_binary(ofs,PL);
  io_binary::write_binary(ofs,DL);
  io_binary::write_binary(ofs,penetrating_muon_flux);
  io_binary::write_binary(ofs,signal);
  io_binary::write_binary(ofs,noise_det);
  io_binary::write_binary(ofs,noise_poi);
  // proj_density, proj_density_lower, proj_density_upper
  io_binary::write_binary(ofs,proj_density);
  io_binary::write_binary(ofs,proj_density_lower);
  io_binary::write_binary(ofs,proj_density_upper);
  // eff_low, eff_cnt, eff_upp (added in PIPELINE_VERSION 7)
  io_binary::write_binary(ofs,eff_low);
  io_binary::write_binary(ofs,eff_cnt);
  io_binary::write_binary(ofs,eff_upp);

  // vec_tf_in_PL
  io_binary::write_vec_tp_bool_double(ofs,vec_tf_in_PL);

}

// Load all members from binary stream
void DetectorElement::load( std::ifstream &ifs )
{
  // unique_index
  unique_index = io_binary::read_binary<int>(ifs);
  // detector id
  detid = io_binary::read_binary<int>(ifs);
  // id_in_this_detector
  id_in_this_detector = io_binary::read_binary<int>(ifs);
  // ray3d
  ray3d.load(ifs);
  // txmin, txmax, tymin, tymax
  txmin = io_binary::read_binary<double>(ifs);
  txmax = io_binary::read_binary<double>(ifs);
  tymin = io_binary::read_binary<double>(ifs);
  tymax = io_binary::read_binary<double>(ifs);

  // angle_unit, 0=Tangent, 1=Radian, 2=Degree
  angle_unit = static_cast<AngleUnit>(io_binary::read_binary<int>(ifs));

  // effective_area_m2, solid_angle, exposure_time_sec
  effective_area_m2 = io_binary::read_binary<double>(ifs);
  solid_angle = io_binary::read_binary<double>(ifs);
  exposure_time_sec = io_binary::read_binary<double>(ifs);
  // PL, DL, penetrating_muon_flux, signal, noise (det floor + poi bucket, serialized separately)
  PL = io_binary::read_binary<double>(ifs);
  DL = io_binary::read_binary<double>(ifs);
  penetrating_muon_flux = io_binary::read_binary<double>(ifs);
  signal = io_binary::read_binary<double>(ifs);
  noise_det = io_binary::read_binary<double>(ifs);
  noise_poi = io_binary::read_binary<double>(ifs);
  // proj_density, proj_density_lower, proj_density_upper
  proj_density = io_binary::read_binary<double>(ifs);
  proj_density_lower = io_binary::read_binary<double>(ifs);
  proj_density_upper = io_binary::read_binary<double>(ifs);
  // eff_low, eff_cnt, eff_upp (added in PIPELINE_VERSION 7)
  eff_low = io_binary::read_binary<double>(ifs);
  eff_cnt = io_binary::read_binary<double>(ifs);
  eff_upp = io_binary::read_binary<double>(ifs);

  // vec_tf_in_PL
  vec_tf_in_PL = io_binary::read_vec_tp_bool_double(ifs);
}