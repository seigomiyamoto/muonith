/// @file cls_Voxel.cpp
/// @brief Implementation of the Voxel class
#include "cls_Voxel.hpp"

#include <fstream>

#include "ns_myapp.hpp"
#include "ns_mymacro.hpp"
#include "ns_mylogger.hpp"
#include "spdlog_pch.hpp"

//############################################################################
// class Voxel
//############################################################################

//============================================================================
// Operators
//============================================================================

bool Voxel::operator!=(const Voxel& other) const
{
#ifdef NODEBUG
  if (tf_exist != other.tf_exist) return true;
  if (density != other.density) return true;
#else
  if (tf_exist != other.tf_exist) { LOG_WARN("Voxel: tf_exist differs"); return true; }
  if (density != other.density) { LOG_WARN("Voxel: density differs"); return true; }
#endif
  return false;
}

//============================================================================
// Binary I/O
//============================================================================

void Voxel::save(std::ofstream& ofs) const
{
  io_binary::write_bool(ofs, tf_exist);
  io_binary::write_binary(ofs, density);
  if (ofs.fail()) THROW_ERROR("Voxel::save: stream write failed.");
}

void Voxel::load(std::ifstream& ifs)
{
  tf_exist = io_binary::read_bool(ifs);
  density = io_binary::read_binary<double>(ifs);
  if (ifs.fail()) THROW_ERROR("Voxel::load: stream read failed.");
}
