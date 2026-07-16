/// @file ns_dem_operator.hpp
/// @brief Operations for creating shell structures from Grid2dPillar and Grid3dVoxel
/// @details
/// This module provides utilities for geometric operations on grid-based terrain models.
/// The primary workflow is:
/// 1. Subtract a Grid3dVoxel from a Grid2dPillar to create a shell structure
/// 2. Calculate average density of the combined shell + voxel system
///
/// @note Thread-safety: Functions use OpenMP parallel regions internally.
///       Ensure thread-safe usage of input/output objects.
#pragma once

#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>
#include <chrono>
#include <thread>

#define _USE_MATH_DEFINES // Required before <cmath> to use M_PI
#include <cmath>

#include <map>
#include <unordered_map>
#include <algorithm>
#include <cassert>
#include <vector>
#include <filesystem> // for std::filesystem::exists
#include <limits>
#include <tuple>
#include <stdexcept>  // Required for std::runtime_error
#include "cls_Grid2dPillar.hpp"
#include "cls_Grid3dVoxel.hpp"

//##################################################################################
//##################################################################################
/// @namespace dem_operator
/// @brief Functions for creating a shell Grid2dPillar by subtracting Grid3dVoxel from Grid2dPillar
/// @details
/// This namespace contains operations for combining 2D cuboid grids with 3D voxel grids.
/// The main use case is creating a "shell" structure by adjusting the lower z-bounds of
/// cuboids based on the highest occupied voxel in each vertical column.
///
/// @ingroup geometryClasses
/// @ingroup terrainClasses
//##################################################################################
//##################################################################################
namespace dem_operator {

  /// @brief Result structure for triple-shell generation
  /// @details Contains upper, lower, and lateral shell Grid2dPillars and flags
  ///          indicating whether each shell has valid cuboids.
  struct ShellTriples {
    Grid2dPillar upper;   ///< U: Shell above g3vox [z_highest, terrain_surface] for columns with existing voxels.
    Grid2dPillar lower;   ///< D: Shell below g3vox [g2pil.zmin, lower_bound_z] for columns with existing voxels.
    Grid2dPillar lateral; ///< L: Full-height shell [g2pil.zmin, terrain_surface] for columns with all tf_exist=false.
    bool has_upper;       ///< True if any upper shell cuboids have positive height.
    bool has_lower;       ///< True if any lower shell cuboids have positive height.
    bool has_lateral;     ///< True if any lateral shell cuboids have positive height.
  };

  /// @brief Creates upper, lower, and lateral shell Grid2dPillars from g2pil and g3vox
  /// @param[in] g2pil Input 2D cuboid grid (original terrain)
  /// @param[in] g3vox Input 3D voxel grid
  /// @return ShellTriples containing upper, lower, and lateral shells
  /// @details
  /// For each (ix,iy) position:
  /// - If the column has no existing voxels (z_lowest > z_highest), it is a lateral column:
  ///   - Lateral shell covers [terrain_zmin, terrain_zmax]
  ///   - Upper and lower shells are empty
  /// - Otherwise (normal column with existing voxels):
  ///   - Upper shell: If terrain surface > z_highest, creates cuboid [z_highest, terrain_surface]
  ///   - Lower shell: If terrain base < lower_bound_z, creates cuboid [terrain_base, min(terrain_surface, lower_bound_z)]
  ///   - Lateral shell is empty
  /// - Pillars with zero or negative height are marked as empty (zmin==zmax)
  ///
  /// @note Thread-safety: Yes. Uses OpenMP parallel for with collapse(2).
  ///       Each merged (ix,iy) maps to a disjoint set of non-merged (ix2,iy2),
  ///       so different threads write to different cuboids with no data race.
  ///       The has_upper/has_lower/has_lateral flags use OpenMP reduction.
  ShellTriples make_shell_triples_from_g2pil_and_g3vox(
      const Grid2dPillar &g2pil, const Grid3dVoxel &g3vox);

  /// @brief Creates a shell Grid2dPillar by subtracting Grid3dVoxel from Grid2dPillar
  /// @param[in] g2pil Input 2D cuboid grid
  /// @param[in] g3vox Input 3D voxel grid
  /// @return Modified Grid2dPillar where each cuboid's z_min is adjusted to the highest
  ///         occupied voxel z-coordinate in the corresponding (ix,iy) column
  /// @details
  /// For each (ix,iy) position in the voxel grid:
  /// 1. Finds the highest z-coordinate where a voxel exists (tf_exist==true)
  /// 2. Updates all corresponding cuboids in g2pil to have z_min = z_highest
  /// 3. If no voxel exists in a column, uses the voxel grid's z_min
  ///
  /// This effectively "carves out" the voxel volume from the cuboid grid, creating a shell.
  ///
  /// @note Thread-safety: No. This function is not thread-safe due to sequential processing
  ///       and internal state modifications.
  /// @note Progress output: Writes progress to stderr using fprintf (legacy pattern).
  Grid2dPillar make_shell_from_g2pil_and_g3vox(
      const Grid2dPillar &g2pil, const Grid3dVoxel &g3vox);

  /// @brief Calculates the average density of shell + g3vox
  /// @param[in] shell Shell Grid2dPillar (typically output from make_shell_from_g2pil_and_g3vox)
  /// @param[in] g3vox 3D voxel grid
  /// @return Volume-weighted average density [kg/m³] of the combined system
  /// @details
  /// Computes the weighted average density as:
  ///   dens_avr = (mass_shell + mass_vox) / (volume_shell + volume_vox)
  ///
  /// For the shell:
  /// - Volume of each cuboid = dx * dy * (zmax - zmin)
  /// - Mass = volume * density
  ///
  /// For the voxel grid:
  /// - Only voxels with tf_exist==true contribute
  /// - Volume of each voxel = dx * dy * dz
  /// - Mass = volume * density
  ///
  /// @note Units: Assumes consistent units across all grids (typically meters, kg/m³).
  /// @note Thread-safety: Yes. Uses OpenMP parallel reductions with no data races.
  /// @note Complexity: O(nx_shell * ny_shell + nx_vox * ny_vox * nz_vox)
  double calc_avr_dens(
    const Grid2dPillar &shell, const Grid3dVoxel &g3vox);

}; // namespace dem_operator