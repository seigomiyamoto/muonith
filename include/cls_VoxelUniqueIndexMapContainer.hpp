/// @file cls_VoxelUniqueIndexMapContainer.hpp
/// @brief Voxel unique index mapping container
/// @details Manages bidirectional mappings between voxel grid coordinates (ix, iy, iz)
///          and unique voxel indices (uqiv) for efficient voxel access.
///
/// @par Typical Workflow:
/// 1. Create a VoxID instance
/// 2. Reserve memory for expected voxel count with reserve_uqiv_umps()
/// 3. Insert voxel mappings with insert_to_uqiv_umps()
/// 4. Set min/max unique index bounds with set_uqiv_min()/set_uqiv_max()
/// 5. Query mappings with get_ixiyiz() (uqiv -> coords) or get_uqiv() (coords -> uqiv)
/// 6. Serialize/deserialize with save()/load() for persistence
///
/// @par Terminology:
/// - uqiv: Unique Voxel Index - a unique integer identifier for each voxel
/// - ixiyiz: Voxel grid coordinates as (ix, iy, iz)
/// - ump: Unordered map container
///
/// @par Thread Safety:
/// Not thread-safe. External synchronization required for concurrent access.
#pragma once

#include <map>
#include <fstream>
#include <iostream>
#include <sstream> // istringstream
#include <string>
#include <cstdio>
#include <set>
#include <cmath>
#include <functional>  //for sorting,std::hash
#include <algorithm>//for sorting
#include <vector>
#include <filesystem> // for std::filesystem::path
#include <stdexcept>  // for std::runtime_error
#include <array>
#include <unordered_map>
#include <cstddef>   // std::size_t

#include "ns_tuple_int.hpp"
#include "cls_Grid3d.hpp"
#include "cls_UOBimap.hpp"


//#####################################################################################
/// @namespace id_container
/// @brief Classes for assigning and storing unique IDs
//#####################################################################################
namespace id_container {

  //##################################################################################
  //##################################################################################
  /// @class VoxID
  /// @brief Container for assigning and storing unique IDs to multiple voxels
  /// @details This class manages bidirectional mappings between unique voxel indices
  ///          (uqiv) and voxel grid coordinates (ix, iy, iz). It provides efficient
  ///          lookup in both directions using unordered_map containers.
  ///
  /// @par Typical Workflow:
  /// 1. Create a VoxID instance
  /// 2. Reserve memory using reserve_uqiv_umps() for expected voxel count
  /// 3. Insert voxel mappings using insert_to_uqiv_umps()
  /// 4. Set min/max uqiv values with set_uqiv_min()/set_uqiv_max()
  /// 5. Query mappings using get_ixiyiz() or get_uqiv()
  ///
  /// @par Thread Safety:
  /// Not thread-safe. External synchronization required for concurrent access.
  ///
  /// @ingroup basicTools
  //##################################################################################
  //##################################################################################
  class VoxID {
    public:
    //======================================================================
    /// @name constants
    ///@{
    static constexpr int uqiv_not_assigned = -1;
    ///@} ------------------------------------------------------------------

    private:
    /// @brief name of this class
    std::string name = "VoxID";

    /// @brief min of unique_index for Voxels
    Grid3d::Uqiv uqiv_min = uqiv_not_assigned;

    /// @brief max of unique_index for Voxels
    Grid3d::Uqiv uqiv_max = uqiv_not_assigned;

    /// @brief Bidirectional map between uqiv and (ix,iy,iz).
    /// @note uqiv=the index assigned for all Voxels, uniquely.
    UOBimap<Grid3d::Uqiv, Grid3d::Ixiyiz,
            std::hash<Grid3d::Uqiv>, std::equal_to<Grid3d::Uqiv>,
            Grid3d::IxiyizHash, Grid3d::IxiyizEqual> bimap_uqiv_ixiyiz{
      Grid3d::UqivNotFound, Grid3d::IxiyizNotFound};

  public:

    //======================================================================
    /// @name constructor_destructor
    ///@{

    /// @brief default constructor
    VoxID() = default;

    /// @brief copy constructor
    VoxID(const VoxID &org) = default;

    /// @brief assignment operator
    VoxID& operator=(const VoxID& other) = default;
    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name operator
    ///@{

    /// @brief Inequality operator
    /// @param[in] other The VoxID instance to compare with
    /// @return true if the objects differ, false otherwise
    /// @note The name member is not compared
    bool operator!=(const VoxID& other) const;

    /// @brief Equality operator (defined using inequality operator)
    /// @param[in] other The VoxID instance to compare with
    /// @return true if the objects are equal, false otherwise
    bool operator==(const VoxID& other) const {
      return !(*this != other);
    };
    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name getter
    ///@{
    /// @brief get the name of VoxID
    std::string get_name() const { return name; };

    /// @brief get the min of unique_index for Voxels
    Grid3d::Uqiv get_uqiv_min() const { return uqiv_min; };
    
    /// @brief get the max of unique_index for Voxels
    Grid3d::Uqiv get_uqiv_max() const { return uqiv_max; };
    
    /// @brief get the reference of ump_uqiv_ixiyiz
    const Grid3d::UmpUqivIxiyiz& get_ump_uqiv_ixiyiz_ref() const { return bimap_uqiv_ixiyiz.getMapAB(); };

    /// @brief get the reference of ump_ixiyiz_uqiv
    const Grid3d::UmpIxiyizUqiv& get_ump_ixiyiz_uqiv_ref() const { return bimap_uqiv_ixiyiz.getMapBA(); };

    /// @brief get the copy of ump_uqiv_ixiyiz
    Grid3d::UmpUqivIxiyiz get_ump_uqiv_ixiyiz_copy() const { return bimap_uqiv_ixiyiz.getMapAB(); };

    /// @brief Get a copy of ump_ixiyiz_uqiv
    /// @return Copy of the (ix,iy,iz) to uqiv mapping
    Grid3d::UmpIxiyizUqiv get_ump_ixiyiz_uqiv_copy() const { return bimap_uqiv_ixiyiz.getMapBA(); };

    /// @brief Get const reference to the uqiv-ixiyiz bimap.
    const auto& get_bimap_uqiv_ixiyiz_ref() const { return bimap_uqiv_ixiyiz; };

    /// @brief Get a copy of the uqiv-ixiyiz bimap.
    auto get_bimap_uqiv_ixiyiz_copy() const { return bimap_uqiv_ixiyiz; };

    /// @brief Get (ix,iy,iz) coordinates using uqiv
    /// @param[in] uqiv_in Unique voxel index to look up
    /// @return Voxel coordinates, or Grid3d::IxiyizNotFound if not found
    /// @throws std::runtime_error if uqiv_in is out of [uqiv_min, uqiv_max] range
    Grid3d::Ixiyiz get_ixiyiz( const int uqiv_in ) const;

    /// @brief Get unique_index using ix, iy, iz
    /// @param[in] ix_in X grid index
    /// @param[in] iy_in Y grid index
    /// @param[in] iz_in Z grid index
    /// @return Unique voxel index, or Grid3d::UqivNotFound if not found
    Grid3d::Uqiv get_uqiv( const int ix_in, const int iy_in, const int iz_in ) const;

    /// @brief Get unique_index using (ix, iy, iz) array
    /// @param[in] ixiyiz Voxel grid coordinates
    /// @return Unique voxel index, or Grid3d::UqivNotFound if not found
    Grid3d::Uqiv get_uqiv( const Grid3d::Ixiyiz &ixiyiz ) const;

    /// @brief Get end iterator of ump_uqiv_ixiyiz
    /// @return End iterator of the uqiv to ixiyiz mapping
    auto get_end_ump_uqiv_ixiyiz() const { return bimap_uqiv_ixiyiz.getMapAB().end(); };

    /// @brief Get end iterator of ump_ixiyiz_uqiv
    /// @return End iterator of the ixiyiz to uqiv mapping
    auto get_end_ump_ixiyiz_uqiv() const { return bimap_uqiv_ixiyiz.getMapBA().end(); };

    /// @brief Get the vector of keys from ump_uqiv_ixiyiz
    /// @return Vector of unique voxel indices, sorted in ascending order
    /// @throws std::runtime_error if keys are not continuous (via THROW_ERROR)
    /// @note Keys are sorted in ascending order and checked for continuity.
    ///       Complexity: O(n log n) due to sorting operation.
    std::vector<Grid3d::Uqiv> get_vec_sorted_uqiv() const;

    /// @brief Get the set of all uqiv (unique index of voxel) in ump_uqiv_ixiyiz
    /// @return Set of unique voxel indices, automatically sorted
    /// @note Complexity: O(n log n) for set insertion
    std::set<Grid3d::Uqiv> get_set_uqiv() const;

    /// @brief Get the vector of all uqiv (unique index of voxel) in ump_uqiv_ixiyiz
    /// @return Vector of unique voxel indices (unsorted)
    /// @note Complexity: O(n) for copying map keys
    std::vector<Grid3d::Uqiv> get_vec_uqiv() const;


    ///@} ------------------------------------------------------------------
    
    //======================================================================
    /// @name setter
    ///@{

    /// @brief Set the name of VoxID
    /// @param[in] name_in Name string to set
    void set_name( const std::string &name_in ){ name = name_in; };

    /// @brief Set the minimum unique_index for Voxels
    /// @param[in] value Minimum uqiv value
    void set_uqiv_min( const Grid3d::Uqiv value ){ uqiv_min = value; };

    /// @brief Set the maximum unique_index for Voxels
    /// @param[in] value Maximum uqiv value
    void set_uqiv_max( const Grid3d::Uqiv value ){ uqiv_max = value; };

    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name output
    ///@{
    
    /// @brief Output ump_uqiv_ixiyiz map to file
    /// @param[in] pathout Output file path
    void out_ump_uqiv_ixiyiz(const std::filesystem::path& pathout) const;

    /// @brief Output ump_ixiyiz_uqiv map to file
    /// @param[in] pathout Output file path
    void out_ump_ixiyiz_uqiv(const std::filesystem::path& pathout) const;
    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name reserve,insert,clear,find
    ///@{

    /// @brief Allocate memory for ump_uqiv_ixiyiz and ump_ixiyiz_uqiv
    /// @param[in] size_in Number of elements to reserve
    void reserve_uqiv_umps( const size_t size_in );

    /// @brief Insert {uqiv, (ix,iy,iz)} mapping to both maps
    /// @param[in] uqiv_in Unique voxel index
    /// @param[in] tp_ixiyiz Voxel grid coordinates
    /// @note Inserts into both bidirectional maps. Average complexity: O(1).
    void insert_to_uqiv_umps(
      const Grid3d::Uqiv uqiv_in, const Grid3d::Ixiyiz& tp_ixiyiz);

    /// @brief Find (ix,iy,iz) using uqiv and return iterator
    /// @param[in] uqiv_in Unique voxel index to find
    /// @return Iterator to the found element, or end() if not found
    /// @note Average complexity: O(1) for hash table lookup
    auto find_ixiyiz( const int uqiv_in ) const {
      return bimap_uqiv_ixiyiz.getMapAB().find(uqiv_in);
    };

    /// @brief Clear ump_uqiv_ixiyiz and ump_ixiyiz_uqiv
    void clear_uqiv_umps();
    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name binary_io
    ///@{
    
    /// @brief Save VoxID data to binary output stream
    /// @param[in,out] ofs Output file stream (must be opened in binary mode)
    /// @throws std::runtime_error if write operation fails
    void save( std::ofstream& ofs ) const;

    /// @brief Load VoxID data from binary input stream
    /// @param[in,out] ifs Input file stream (must be opened in binary mode)
    /// @throws std::runtime_error if read operation fails
    void load( std::ifstream &ifs );

    ///@} ------------------------------------------------------------------
}; // End of class VoxID

}; // end of namespace id_container

