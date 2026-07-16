/// @file exemdl_load_parameters.hpp
/// @brief JSON parameter loading module
/// @details
/// Provides functions for loading and parsing all configuration parameters from JSON files.
/// This module handles the complete deserialization workflow for application parameters,
/// including detector configurations, grid parameters, path calculation settings, noise
/// parameters, flux tables, and inversion parameters.
///
/// Workflow:
/// 1. Call load_all() with a validated JSON object
/// 2. Individual parameter sections are extracted and validated
/// 3. Returns a fully populated AppParameters structure
///
/// Thread-safety: No. JSON parsing and parameter construction are not thread-safe.
/// All loading should be performed during single-threaded initialization.
#pragma once

#include <nlohmann/json.hpp>
#include "cls_DetectorPanelArray.hpp"
#include "cls_Grid2dPillar.hpp"
#include "cls_Grid2dPillarParameters.hpp"
#include "cls_Grid2dBinGroup.hpp"
#include "cls_Grid3dVoxel.hpp"
#include "ns_pathcalc.hpp"
#include "cls_NoiseParameters.hpp"
#include "cls_NagaInvParameters.hpp"

/// @namespace exemdl::load_parameters
/// @brief Namespace for loading parameters from JSON
/// @ingroup ExecModule
namespace exemdl::load_parameters {

  /// @brief Range information for Z_CROSS_SECTION
  /// @details
  /// Defines the z-axis range and step size for cross-sectional output.
  /// Used to specify which z-planes to export during reconstruction.
  ///
  /// Invariants:
  /// - min < max
  /// - step > 0
  ///
  /// @note Units: meters (assumed to match coordinate system units)
  struct ZCrossSection {
    double min;  ///< Minimum z-coordinate
    double max;  ///< Maximum z-coordinate
    double zstep = 0.0; ///< Z-step for cross-section output [meters]. 0 = use z_interval
    bool output_binary = false; ///< If true, output in binary format; if false, text format
  };

  /// @brief Application-wide parameter collection
  /// @details
  /// Aggregates all parameters needed for muon tomography reconstruction.
  /// Contains detector geometry, grid configurations, path calculation settings,
  /// noise models, flux tables, and inversion parameters.
  ///
  /// Typical usage:
  /// @code
  /// nlohmann::json js = /* load from file */;
  /// auto params = exemdl::load_parameters::load_all(js);
  /// // Use params.prm_det, params.prm_g3vox, etc.
  /// @endcode
  ///
  /// @note All nested parameter structures are validated during construction.
  /// Invalid or missing required fields will throw exceptions.
  struct AppParameters {
    ZCrossSection                  zcross; ///< Range information for Z_CROSS_SECTION
    DetectorPanel::ParameterLists  prm_det; ///< Detector panel parameter lists
    Grid2dPillar::Parameters       prm_g2pil; ///< Grid2dPillar parameters
    Grid2dBinGroup::Parameters     prm_bingroup; ///< Grid2dBinGroup parameters
    Grid3dVoxel::Parameters        prm_g3vox; ///< Grid3dVoxel parameters
    Grid3dVoxel::MergeParameters   prm_g3merge; ///< Grid3dVoxel merge parameters
    Grid3dVoxel::ReconstVoxelsParameters prm_reconst_voxels; ///< Grid3dVoxel reconstruction voxels parameters
    pathcalc::Parameters           prm_path; ///< Path calculation parameters
    NoiseParameters                prm_noise; ///< Noise parameters
    FluxTable                      ft_real; ///< Real flux table
    FluxTable                      ft_prior;  ///< Prior flux table
    std::vector<NagaInvParameters> vec_prm_nagainv; ///< List of NagaInv parameters
    bool tf_run_inversion_prior_error = true; ///< If true, compute prior error with lower/upper bounds; if false, center only
  };

  /// @brief Loads all parameters from JSON and populates the structure
  /// @details
  /// Extracts and validates all required parameter sections from the JSON object.
  /// Performs range validation on Z_CROSS_SECTION parameters (min < max, step > 0).
  /// Constructs nested parameter objects by delegating to their respective constructors.
  ///
  /// @param[in] js Loaded and validated JSON object containing all configuration sections
  /// @return Fully populated AppParameters object
  /// @throws std::runtime_error If Z_CROSS_SECTION section is missing
  /// @throws std::runtime_error If z_cross_min >= z_cross_max
  /// @throws std::runtime_error If z_step <= 0
  /// @throws std::runtime_error If any required parameter section is missing or invalid
  /// @note Thread-safety: No. Call during single-threaded initialization only.
  AppParameters load_all(const nlohmann::json& js);

  /// @brief Loads a list of NagaInvParameters from JSON
  /// @details
  /// Extracts an array of inversion parameter sets from the specified JSON section.
  /// Each element in the array is passed to NagaInvParameters constructor for validation.
  /// The list must contain at least one parameter set.
  ///
  /// @param[in] js Loaded and validated JSON object
  /// @param[in] section_name Section name (e.g., "NAGAINV_PARAMETERS")
  /// @return List of NagaInvParameters (at least one element)
  /// @throws std::runtime_error If section_name is empty
  /// @throws std::runtime_error If section is not found in JSON
  /// @throws std::runtime_error If no NagaInvParameters found in section (empty array)
  /// @note Thread-safety: No. Call during single-threaded initialization only.
  std::vector<NagaInvParameters>
    load_nagainv_parameters(const nlohmann::json& js, const std::string& section_name);

}; // namespace exemdl::load_parameters
