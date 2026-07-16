/// @file cls_Voxel.hpp
/// @brief Voxel data structure for 3D grid
/// @details Defines Voxel class storing density and existence flag
///          for muography inversion. Hit-related data (per-detector hit flags,
///          hit-element counters) are managed at Grid3dVoxel level.
///
/// ## Typical workflow
/// 1. Create a Voxel with default constructor
/// 2. Set existence flag via set_tf_exist()
/// 3. Set density (kg/m^3) via set_density()
///
/// ## Units
/// - Density: kg/m^3 (SI units)
///
/// ## Thread safety
/// - Not thread-safe. External synchronization required for concurrent access.
///
/// ## Binary I/O
/// - save() / load() use internal binary format via io_binary utilities.
#pragma once

#include <iosfwd>

#include "ns_type_definitions.hpp"

using namespace index_type_definitions;

/// @brief Type alias for hit count (unsigned short, max 65535)
using Nhit = unsigned short int;

//#########################################################################
/// @class Voxel
/// @brief Data container for a single voxel in a 3D grid
///
/// @details Stores per-voxel attributes for muography inversion:
///          - Material existence flag
///          - Density value (kg/m^3, SI units)
///
/// Hit-related data (per-detector flags, element hit counts) are stored
/// at Grid3dVoxel level.
///
/// ## Usage example
/// @code
///   Voxel v;
///   v.set_tf_exist(true);
///   v.set_density(2700.0);  // kg/m^3
/// @endcode
///
/// @ingroup terrainClasses
/// @ingroup matrixClasses
/// @ingroup geometryClasses
//#########################################################################
class Voxel {
 private:
  /// @brief True if material exists in this voxel
  bool tf_exist = false;

  /// @brief Density value in kg/m^3 (SI units)
  /// @note Historically was g/cm^3, now converted to SI (kg/m^3).
  double density = 0.0;

 public:
  //======================================================================
  /// @name Constructors and Destructor
  ///@{

  /// @brief Default constructor
  Voxel() = default;

  /// @brief Copy constructor
  Voxel(const Voxel& org) = default;

  /// @brief Move constructor
  Voxel(Voxel&& other) noexcept = default;

  /// @brief Destructor
  ~Voxel() = default;
  ///@} ------------------------------------------------------------------

  //======================================================================
  /// @name Operators
  ///@{

  /// @brief Copy all member data from another Voxel
  /// @param[in] org Source Voxel to copy from
  /// @note Prefer copy assignment operator (operator=) for standard semantics
  void copy(const Voxel& org) {
    tf_exist = org.tf_exist;
    density = org.density;
  }

  /// @brief Inequality operator
  /// @param[in] other Voxel to compare against
  /// @return True if any member differs
  bool operator!=(const Voxel& other) const;

  /// @brief Equality operator (defined via inequality)
  /// @param[in] other Voxel to compare against
  /// @return True if all members are equal
  bool operator==(const Voxel& other) const { return !(*this != other); }

  /// @brief Copy assignment operator
  Voxel& operator=(const Voxel& other) = default;

  /// @brief Move assignment operator
  Voxel& operator=(Voxel&& other) noexcept = default;

  ///@} ------------------------------------------------------------------

  //======================================================================
  /// @name Getters
  ///@{

  /// @brief Get material existence flag
  /// @return True if material exists in this voxel
  bool get_tf_exist() const { return tf_exist; }

  /// @brief Get density value
  /// @return Density in kg/m^3 (SI units)
  double get_density() const { return density; }

  ///@} ------------------------------------------------------------------

  //======================================================================
  /// @name Setters
  ///@{

  /// @brief Set material existence flag
  /// @param[in] tf_exist_in True if material exists
  void set_tf_exist(bool tf_exist_in) { tf_exist = tf_exist_in; }

  /// @brief Set density value
  /// @param[in] density_in Density in kg/m^3 (SI units)
  void set_density(double density_in) { density = density_in; }

  /// @brief Add to density value
  /// @param[in] density_in Density increment in kg/m^3 (SI units)
  void add_density(double density_in) { density += density_in; }

  ///@} ------------------------------------------------------------------

  //======================================================================
  /// @name Binary I/O
  ///@{

  /// @brief Save voxel data to binary stream
  /// @param[in,out] ofs Output file stream (must be open in binary mode)
  /// @throws std::runtime_error if write fails
  void save(std::ofstream& ofs) const;

  /// @brief Load voxel data from binary stream
  /// @param[in,out] ifs Input file stream (must be open in binary mode)
  /// @throws std::runtime_error if read fails
  void load(std::ifstream& ifs);
  ///@} ------------------------------------------------------------------
};
