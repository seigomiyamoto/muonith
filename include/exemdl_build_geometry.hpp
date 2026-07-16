/// @file exemdl_build_geometry.hpp
/// @brief Geometry construction module for terrain and voxel grids.
///
/// @details This file provides the build_geometry namespace which constructs
///          terrain-based Grid2dPillar and voxelized Grid3dVoxel structures
///          from loaded application parameters.
///
///          Typical workflow:
///          1. Load parameters via exemdl::load_parameters.
///          2. Call build_all() with BuildArgs containing app parameters.
///          3. Obtain BuildResult with naive cuboid, input voxel, merged voxel, and shell.
///          4. Use the resulting grids for ray tracing or other computations.
///
///          Coordinate system:
///          - Right-handed, z-axis is vertical (positive upward)
///          - Units are meters
///
///          Thread safety:
///          - build_all() is not thread-safe (modifies shared I/O paths)
///
/// @ingroup ExecModule
#pragma once

#include "exemdl_load_parameters.hpp"
#include "cls_Grid2dPillar.hpp"
#include "cls_Grid3dVoxel.hpp"

/// @namespace exemdl::build_geometry
/// @brief Namespace for geometry construction from loaded parameters.
/// @ingroup ExecModule
namespace exemdl::build_geometry {

  /// @brief Input arguments for geometry construction (Module3).
  struct BuildArgs {
    const exemdl::load_parameters::AppParameters& app_params; ///< Application parameters.
  };

  /// @brief Result of geometry construction (Module3).
  ///
  /// @details Contains geometry objects built from terrain and voxel parameters:
  ///          - g2pil_naive: Grid2dPillar built from raw terrain data.
  ///          - g3vox_input: Grid3dVoxel converted from g2pil_naive with input density.
  ///          - g3vox_merged_input: Lower-resolution merged Grid3dVoxel.
  ///          - g2pil_shell_upper: Upper shell (terrain above g3vox.z_highest).
  ///          - g2pil_shell_lower: Lower shell (terrain below g3vox.zmin).
  ///          - g2pil_shell_lateral: Lateral shell (full-height terrain for columns outside AABB).
  struct BuildResult {
    Grid2dPillar g2pil_naive;        ///< Grid2dPillar built from terrain data.
    Grid3dVoxel  g3vox_input;        ///< Grid3dVoxel with input density structure.
    Grid3dVoxel  g3vox_merged_input; ///< Merged Grid3dVoxel at lower resolution.
    Grid2dPillar g2pil_shell_upper;   ///< Upper shell: terrain above g3vox [z_highest, surface].
    Grid2dPillar g2pil_shell_lower;   ///< Lower shell: terrain below g3vox [base, g3vox.zmin].
    Grid2dPillar g2pil_shell_lateral; ///< Lateral shell: full-height terrain for columns outside AABB.
    bool         has_shell_upper   = false; ///< True if upper shell has valid cuboids.
    bool         has_shell_lower   = false; ///< True if lower shell has valid cuboids.
    bool         has_shell_lateral = false; ///< True if lateral shell has valid cuboids.
  };

  /// @brief Execute all geometry construction steps.
  ///
  /// @details Builds terrain cuboid, converts to voxel grid, merges voxels, and
  ///          extracts outer shell. If prm_g3vox.tf_build_g3vox is false, only
  ///          g2pil_naive is built and cross-section output is generated.
  ///
  /// @param[in] args Input parameters for construction.
  /// @return BuildResult containing all constructed geometry objects.
  ///
  /// @note Thread-safe: No (uses shared I/O directories).
  /// @note Complexity: O(N) where N is the number of voxels.
  BuildResult build_all(const BuildArgs& args);

} // namespace exemdl::build_geometry
