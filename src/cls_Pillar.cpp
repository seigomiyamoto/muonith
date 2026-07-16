/// @file cls_Pillar.cpp
/// @brief Implementation of the Pillar class.

#include "cls_Pillar.hpp"
#include "ns_io_binary.hpp"
#include "spdlog_pch.hpp"
#include "ns_mylogger.hpp"


//=====================
// class Pillar
//=====================

void Pillar::set_values(
  const double zmin_in, const double zmax_in, const double density_in )
{
  zmin = zmin_in;
  zmax = zmax_in;
  density = density_in;
}

bool Pillar::is_z_inside( const double z_in ) const
{
  if( z_in <  zmin ) return false;
  if( z_in >= zmax ) return false;
  return true;
}

bool Pillar::is_below_zmax( const double z_in ) const
{
  if( z_in >= zmax ) return false;
  return true;
}

void Pillar::save(std::ofstream& ofs) const
{
  io_binary::write_binary(ofs, density);
  io_binary::write_binary(ofs, zmin);
  io_binary::write_binary(ofs, zmax);
}

void Pillar::load(std::ifstream& ifs)
{
  density = io_binary::read_binary<double>(ifs);
  zmin = io_binary::read_binary<double>(ifs);
  zmax = io_binary::read_binary<double>(ifs);
}
