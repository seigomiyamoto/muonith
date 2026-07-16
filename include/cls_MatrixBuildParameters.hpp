/// @file cls_MatrixBuildParameters.hpp
/// @brief Matrix construction parameter configuration
/// @details
/// Defines pathcalc::MatrixBuildParameters, a parameter bundle class that aggregates
/// references to all data structures required for matrix construction in muon path length
/// calculations. This class avoids expensive copying of large objects by holding only
/// references.
///
/// ## Typical workflow
/// 1. Create instances of DetectorPanelArray, Grid3dVoxel, FluxTable, and sparse matrices
/// 2. Construct MatrixBuildParameters with references to these objects
/// 3. Pass the MatrixBuildParameters to calc_dNdD functions
/// 4. Optionally save/load the observation matrix from disk using the path specified
///
/// ## Thread safety
/// - **Not thread-safe for writes**: Multiple threads must not modify the same
///   MatrixBuildParameters instance concurrently
/// - **Read-only safe**: If all referenced objects are immutable, multiple threads
///   can read safely
/// - **OpenMP usage**: The referenced objects (e.g., vec_spmat_PL) may be accessed
///   within OpenMP parallel regions in pathcalc functions
///
/// ## Memory management
/// - This class does **not** own any resources; it only holds references
/// - The caller must ensure all referenced objects outlive this parameter bundle
/// - vec_spmat_PL elements are gradually released (resize(0,0)) during processing
///   to reduce memory footprint
#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include "cls_DetectorPanelArray.hpp"
#include "cls_Grid3dVoxel.hpp"
#include "cls_FluxTable.hpp"

namespace pathcalc {
  class MatrixBuildParameters; // forward declaration
};


//####################################
/// @class pathcalc::MatrixBuildParameters
/// @brief Parameter bundle for matrix construction in path calculation
/// @details
/// This class holds references to core data structures required for matrix construction:
/// DetectorPanelArray, Grid3dVoxel, sparse path length matrices (std::vector<SpMatf>),
/// prior density vectors, and FluxTable. It is passed by reference to pathcalc functions
/// to avoid copying large data structures.
///
/// ## Responsibilities
/// - Aggregate references to detector geometry, voxel grid, flux table, and matrices
/// - Control observation matrix load/save behavior via boolean flags
/// - Provide a unified interface for passing multiple large objects to pathcalc functions
///
/// ## Usage example
/// @code
/// // 1. Create core data structures
/// DetectorPanelArray detectors;
/// Grid3dVoxel voxelGrid;
/// FluxTable fluxTable;
/// std::vector<SpMatf> pathLengthMatrices;
/// std::vector<Eigen::VectorXf> priorDensities;
///
/// // 2. Configure and populate these structures
/// // ... (load or compute data)
///
/// // 3. Create parameter bundle
/// MatrixBuildParameters params(
///   "myParams",                          // name
///   false,                               // do not load from file
///   true,                                // save result to file
///   "output/obs_matrix.bin",             // output path
///   detectors, voxelGrid,                // geometry references
///   pathLengthMatrices, priorDensities,  // matrix references
///   fluxTable,                           // flux reference
///   true                                 // apply efficiency
/// );
///
/// // 4. Use in pathcalc functions
/// calc_dNdD::make_grouped_alldet_mat_dNdD(params);
/// @endcode
///
/// @note This class does NOT own the referenced objects. Ensure all referenced data
///       structures remain valid for the lifetime of this parameter bundle.
/// @ingroup parameterClasses
//####################################
class pathcalc::MatrixBuildParameters {
  public:

    /// @brief Name identifier for this instance
    std::string name;

    /// @brief If true, load the observation matrix from binary file instead of computing it
    /// @details
    /// - true: Skip calc_dNdD::make_grouped_alldet_mat_dNdD and load from file
    /// - false: Compute the matrix without loading
    /// @note Used to avoid expensive recomputation when the matrix is already available.
    ///       The file must exist at path_bin_obs_mat_dNdD if this is true.
    bool tf_load_bin_obs_mat_dNdD;

    /// @brief If true, save the observation matrix (large_merged_mat_dNdD) as binary file
    /// @details
    /// - true: Save the computed matrix to path_bin_obs_mat_dNdD
    /// - false: Do not save to disk
    /// @note Useful for caching expensive matrix computations. Both tf_load_bin_obs_mat_dNdD
    ///       and tf_save_bin_obs_mat_dNdD can be true (load first, then save if recomputed).
    bool tf_save_bin_obs_mat_dNdD;

    /// @brief File path for binary observation matrix data (obs_mat_dNdD)
    /// @details
    /// This path is used for both saving and loading the matrix computed by
    /// calc_dNdD::make_grouped_alldet_mat_dNdD. The file is typically large
    /// (size depends on detector count × voxel count).
    /// @note File format is binary (platform-dependent endianness). Ensure sufficient
    ///       disk space before enabling tf_save_bin_obs_mat_dNdD.
    std::filesystem::path path_bin_obs_mat_dNdD;

    /// @brief Reference to DetectorPanelArray data
    /// @details Contains detector panel geometry, positions, orientations, and pixel layouts.
    /// @note Non-const because detector state may be modified during matrix construction.
    DetectorPanelArray& arrdet;

    /// @brief Reference to Grid3dVoxel data
    /// @details Represents the 3D voxelized space through which muon paths are calculated.
    ///          Contains voxel dimensions, positions, and material properties.
    /// @note Non-const because voxel state may be updated during processing.
    Grid3dVoxel& g3vox;

    /// @brief Reference to path length matrix for each detector
    /// @details
    /// Sparse matrices (SpMatf = Eigen::SparseMatrix<float>) storing path lengths
    /// through voxels for each detector element. Each element vec_spmat_PL[i] corresponds
    /// to detector i.
    /// @note **Memory release pattern**: Elements are resized to (0,0) after use to
    ///       progressively reduce memory footprint during matrix construction. This is why
    ///       the reference is non-const.
    /// @note **Thread safety**: Elements may be accessed in OpenMP parallel regions.
    ///       Ensure proper synchronization if multiple threads access the same element.
    std::vector<SpMatf>& vec_spmat_PL;

    /// @brief Reference to prior density vector (Eigen::VectorXf vec_vecxf_DL_prior)
    /// @details
    /// Each element vec_vecxf_DL_prior[i] contains the prior density distribution
    /// for a group of voxels. Used in Bayesian inversion or regularization.
    /// @note Const reference because priors are not modified during matrix construction.
    const std::vector<Eigen::VectorXf>& vec_vecxf_DL_prior;

    /// @brief Reference to FluxTable data
    /// @details
    /// Contains muon flux as a function of energy, zenith angle, and azimuth angle.
    /// Used to weight the observation matrix by expected muon flux.
    /// @note Const reference because flux table is read-only during matrix construction.
    const FluxTable& ft;

    /// @brief If true, apply detector efficiency correction to the signal
    /// @details
    /// When enabled, the observation matrix elements are multiplied by detector-specific
    /// efficiency factors (typically energy- and angle-dependent).
    /// @note Efficiency data is expected to be available in the DetectorPanelArray.
    bool tf_apply_eff;

    /// @brief If true, scale the observation matrix by central efficiency eff_cnt
    /// @details Deterministic (NOT the dice path). Gated by tf_eff_cn_diag so that the
    ///          matrix A and the forward count share the same central efficiency.
    bool tf_apply_eff_cnt;

    /// @brief Constructor from input parameters
    /// @param name_in Instance name identifier (for logging/debugging)
    /// @param tf_load_bin_obs_mat_dNdD_in If true, load matrix from file instead of computing
    /// @param tf_save_bin_obs_mat_dNdD_in If true, save computed matrix to file
    /// @param path_bin_obs_mat_dNdD_in Path for loading/saving binary matrix
    /// @param arrdet_in Reference to DetectorPanelArray (must outlive this object)
    /// @param g3vox_in Reference to Grid3dVoxel (must outlive this object)
    /// @param vec_spmat_PL_in Reference to path length matrix vector (non-const for memory release)
    /// @param vec_vecxf_DL_prior_in Reference to prior density vector (must outlive this object)
    /// @param frdt_in Reference to FluxTable (must outlive this object)
    /// @param tf_apply_eff_in If true, apply detector efficiency correction
    /// @note vec_spmat_PL_in is passed as non-const reference because memory is gradually
    ///       released during processing (resize(0,0) is called on elements).
    /// @note **Lifetime requirement**: All referenced objects (arrdet_in, g3vox_in, etc.)
    ///       must remain valid for the entire lifetime of this MatrixBuildParameters instance.
    /// @note **Thread safety**: This constructor itself is not thread-safe. Do not construct
    ///       multiple instances with overlapping references from different threads without
    ///       external synchronization.
    MatrixBuildParameters(
          const std::string& name_in
        , bool tf_load_bin_obs_mat_dNdD_in
        , bool tf_save_bin_obs_mat_dNdD_in
        , const std::filesystem::path& path_bin_obs_mat_dNdD_in
        , DetectorPanelArray& arrdet_in
        , Grid3dVoxel& g3vox_in
        , std::vector<SpMatf>& vec_spmat_PL_in
        , const std::vector<Eigen::VectorXf>& vec_vecxf_DL_prior_in
        , const FluxTable& frdt_in
        , bool tf_apply_eff_in
        , bool tf_apply_eff_cnt_in = false
      )
      : name(name_in)
      , tf_load_bin_obs_mat_dNdD(tf_load_bin_obs_mat_dNdD_in)
      , tf_save_bin_obs_mat_dNdD(tf_save_bin_obs_mat_dNdD_in)
      , path_bin_obs_mat_dNdD(path_bin_obs_mat_dNdD_in)
      , arrdet(arrdet_in)
      , g3vox(g3vox_in)
      , vec_spmat_PL(vec_spmat_PL_in)
      , vec_vecxf_DL_prior(vec_vecxf_DL_prior_in)
      , ft(frdt_in)
      , tf_apply_eff(tf_apply_eff_in)
      , tf_apply_eff_cnt(tf_apply_eff_cnt_in)
      {}

    /// @brief Copy constructor
    /// @param org Source instance to copy from
    /// @note **Reference semantics**: All reference members are copied (not deep-copied),
    ///       meaning the new instance refers to the same underlying data as the original.
    ///       Modifying referenced objects through either instance affects both.
    /// @note **Use case**: Useful when you need to create a parameter bundle with slightly
    ///       different flags (e.g., different load/save settings) but the same data references.
    MatrixBuildParameters( const MatrixBuildParameters& org )
      : name(org.name)
      , tf_load_bin_obs_mat_dNdD(org.tf_load_bin_obs_mat_dNdD)
      , tf_save_bin_obs_mat_dNdD(org.tf_save_bin_obs_mat_dNdD)
      , path_bin_obs_mat_dNdD(org.path_bin_obs_mat_dNdD)
      , arrdet(org.arrdet)
      , g3vox(org.g3vox)
      , vec_spmat_PL(org.vec_spmat_PL)
      , vec_vecxf_DL_prior(org.vec_vecxf_DL_prior)
      , ft(org.ft)
      , tf_apply_eff(org.tf_apply_eff)
      , tf_apply_eff_cnt(org.tf_apply_eff_cnt)
      {}

    /// @brief Destructor
    /// @note Default destructor is sufficient because this class only holds references
    ///       and does not own any dynamically allocated resources
    ~MatrixBuildParameters() = default;
};

