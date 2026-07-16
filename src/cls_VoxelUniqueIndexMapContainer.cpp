/// @file cls_VoxelUniqueIndexMapContainer.cpp
/// @brief Implementation of VoxID class for voxel unique index mapping
#include <cassert>
#include <iostream>
#include <fstream>

#include "cls_VoxelUniqueIndexMapContainer.hpp"
#include "ns_myapp.hpp"
#include "ns_mymacro.hpp"
#include "ns_mylogger.hpp"
#include "spdlog_pch.hpp"

//###################################################################################
// class id_container::VoxID
//###################################################################################

/// @brief Inequality operator
/// @note The name member is not compared
bool id_container::VoxID::operator!=(const id_container::VoxID& other) const
{
#ifdef NODEBUG
  if (uqiv_min != other.uqiv_min) return true;
  if (uqiv_max != other.uqiv_max) return true;
  if (bimap_uqiv_ixiyiz != other.bimap_uqiv_ixiyiz) return true;
#else
  if (uqiv_min != other.uqiv_min) { LOG_WARN("uqiv_min differs"); return true; }
  if (uqiv_max != other.uqiv_max) { LOG_WARN("uqiv_max differs"); return true; }
  if (bimap_uqiv_ixiyiz != other.bimap_uqiv_ixiyiz) { LOG_WARN("bimap_uqiv_ixiyiz differs"); return true; }
#endif
  return false;
}

Grid3d::Ixiyiz id_container::VoxID::get_ixiyiz(const int uqiv_in) const
{
  if (uqiv_in < get_uqiv_min()) {
    THROW_ERROR("VoxID::get_ixiyiz: uqiv_in < uqiv_min. uqiv_in={}, uqiv_min={}", uqiv_in, get_uqiv_min());
  }
  if (uqiv_in > get_uqiv_max()) {
    THROW_ERROR("VoxID::get_ixiyiz: uqiv_in > uqiv_max. uqiv_in={}, uqiv_max={}", uqiv_in, get_uqiv_max());
  }
  return bimap_uqiv_ixiyiz.getABorDefault(uqiv_in);
}

Grid3d::Uqiv id_container::VoxID::get_uqiv(const Grid3d::Ixiyiz& ixiyiz) const
{
  return bimap_uqiv_ixiyiz.getBAorDefault(ixiyiz);
}

Grid3d::Uqiv id_container::VoxID::get_uqiv(
  const int ix_in, const int iy_in, const int iz_in) const
{
  const Grid3d::Ixiyiz ixiyiz = {ix_in, iy_in, iz_in};
  return get_uqiv(ixiyiz);
}

std::vector<Grid3d::Uqiv> id_container::VoxID::get_vec_sorted_uqiv() const
{
  std::vector<Grid3d::Uqiv> keys;
  keys.reserve(bimap_uqiv_ixiyiz.size());

  for (const auto& pair : bimap_uqiv_ixiyiz.getMapAB()) {
    keys.push_back(pair.first);
  }

  const bool is_continuity = myapp::sort_and_check_vec_int_continuity(keys);
  if (!is_continuity) {
    THROW_ERROR("VoxID::get_vec_sorted_uqiv: keys are not continuous");
  }

  return keys;
}

std::set<Grid3d::Uqiv> id_container::VoxID::get_set_uqiv() const
{
  return bimap_uqiv_ixiyiz.get_set_A();
}

std::vector<Grid3d::Uqiv> id_container::VoxID::get_vec_uqiv() const
{
  return bimap_uqiv_ixiyiz.get_vec_A();
}

void id_container::VoxID::reserve_uqiv_umps(const size_t size_in)
{
  bimap_uqiv_ixiyiz.reserve(size_in);
}

void id_container::VoxID::insert_to_uqiv_umps(
  const Grid3d::Uqiv uqiv_in, const Grid3d::Ixiyiz& ixiyiz_in)
{
  bimap_uqiv_ixiyiz.insert(uqiv_in, ixiyiz_in);
}

void id_container::VoxID::clear_uqiv_umps()
{
  bimap_uqiv_ixiyiz.clear();
}

void id_container::VoxID::out_ump_uqiv_ixiyiz(const std::filesystem::path& pathout) const
{
  tuple_int::UmpIntInt3 converted_map;

  for (const auto& [uqiv, ixiyiz] : bimap_uqiv_ixiyiz.getMapAB()) {
    converted_map.emplace(uqiv, std::make_tuple(ixiyiz[0], ixiyiz[1], ixiyiz[2]));
  }

  tuple_int::out_ump_int_int3(pathout, converted_map);
}

void id_container::VoxID::out_ump_ixiyiz_uqiv(const std::filesystem::path& pathout) const
{
  tuple_int::UmpInt3Int converted_map;

  for (const auto& [ixiyiz, uqiv] : bimap_uqiv_ixiyiz.getMapBA()) {
    converted_map.emplace(std::make_tuple(ixiyiz[0], ixiyiz[1], ixiyiz[2]), uqiv);
  }

  tuple_int::out_ump_int3_int(pathout, converted_map);
}

void id_container::VoxID::save(std::ofstream& ofs) const
{
  io_binary::write_string(ofs, name);
  io_binary::write_binary(ofs, uqiv_min);
  io_binary::write_binary(ofs, uqiv_max);

  io_binary::write_binary(ofs, bimap_uqiv_ixiyiz.size());
  for (const auto& [key, value] : bimap_uqiv_ixiyiz.getMapAB()) {
    io_binary::write_binary(ofs, key);
    io_binary::write_binary(ofs, value);
  }

  io_binary::write_binary(ofs, bimap_uqiv_ixiyiz.size());
  for (const auto& [key, value] : bimap_uqiv_ixiyiz.getMapBA()) {
    io_binary::write_binary(ofs, key);
    io_binary::write_binary(ofs, value);
  }

  if (ofs.fail()) {
    THROW_ERROR("VoxID::save: Write to binary stream failed");
  }
}

void id_container::VoxID::load(std::ifstream& ifs)
{
  name = io_binary::read_string(ifs);
  uqiv_min = io_binary::read_binary<int>(ifs);
  uqiv_max = io_binary::read_binary<int>(ifs);

  bimap_uqiv_ixiyiz.clear();

  std::size_t size_ump1 = io_binary::read_binary<std::size_t>(ifs);
  for (std::size_t i = 0; i < size_ump1; ++i) {
    Grid3d::Uqiv key = io_binary::read_binary<int>(ifs);
    Grid3d::Ixiyiz value = io_binary::read_binary<std::array<int, 3>>(ifs);
    bimap_uqiv_ixiyiz.insert(key, value);
  }

  // Read and discard reverse map section for backward compatibility
  std::size_t size_ump2 = io_binary::read_binary<std::size_t>(ifs);
  for (std::size_t i = 0; i < size_ump2; ++i) {
    io_binary::read_binary<std::array<int, 3>>(ifs);
    io_binary::read_binary<int>(ifs);
  }

  if (ifs.fail()) {
    THROW_ERROR("VoxID::load: Read from binary stream failed");
  }
}
