/// @file cls_UqigAvailIndexer.cpp
/// @brief Implementation of UqigAvailIndexer class.
/// @details Implements bidirectional mapping operations, including build, lookup, and serialization.
#include "cls_UqigAvailIndexer.hpp"
#include <algorithm>
#include <stdexcept>
#include "ns_io_binary.hpp"
#include "ns_mymacro.hpp"

bool UqigAvailIndexer::operator!=(const UqigAvailIndexer& other) const
{
#ifdef NODEBUG
  if (bimap_UqigAvail_Index != other.bimap_UqigAvail_Index) return true;
#else
  if (bimap_UqigAvail_Index != other.bimap_UqigAvail_Index) { LOG_WARN("UqigAvailIndexer: bimap_UqigAvail_Index differs"); return true; }
#endif
  return false;
}

void UqigAvailIndexer::build(const std::vector<UqigAvail>& vec_avail)
{
  // Sort in ascending order
  std::vector<UqigAvail> vec_avail_sorted = vec_avail;
  std::sort(vec_avail_sorted.begin(), vec_avail_sorted.end());

  // Reserve capacity for efficiency
  bimap_UqigAvail_Index.reserve(vec_avail_sorted.size());

  // Insert consecutive indices (duplicates are handled by UOBimap::insert)
  for (Index i = 0; i < static_cast<Index>(vec_avail_sorted.size()); ++i) {
    bimap_UqigAvail_Index.insert(vec_avail_sorted[i], i);
  }
}

void UqigAvailIndexer::save(std::ofstream& ofs) const
{
  // Write element count
  int n = static_cast<int>(size());
  ofs.write(reinterpret_cast<const char*>(&n), sizeof(n));

  // Write UqigAvail values in index order
  for (int i = 0; i < n; ++i) {
    UqigAvail a = getUqigAvail(i);
    ofs.write(reinterpret_cast<const char*>(&a), sizeof(a));
  }

  // Check for write errors
  if (ofs.fail()) {
    THROW_ERROR("UqigAvailIndexer::save: write operation failed.");
  }
}

void UqigAvailIndexer::load(std::ifstream& ifs)
{
  int n;
  ifs.read(reinterpret_cast<char*>(&n), sizeof(n));
  if (ifs.fail()) {
    THROW_ERROR("UqigAvailIndexer::load: Failed to read element count.");
  }
  if (n < 0) {
    THROW_ERROR("UqigAvailIndexer::load: Invalid element count (n={}). Must be non-negative.", n);
  }
  std::vector<UqigAvail> tmp(n);
  for (auto& a : tmp) {
    ifs.read(reinterpret_cast<char*>(&a), sizeof(a));
  }
  if (ifs.fail()) {
    THROW_ERROR("UqigAvailIndexer::load: Failed to read UqigAvail data.");
  }
  build(tmp);
}
