/// @file cls_NagaInv.hpp
/// @brief Linear inversion solver for muon tomography
/// @details Defines the NagaInv class implementing regularized linear inversion
/// for density reconstruction from muon flux observations.
///
/// @section nagainv_workflow Typical Workflow
/// 1. Construct NagaInv from Grid3dVoxel and NagaInvParameters
/// 2. Build covariance matrices via mp_build_mat_cov_muon() and mp_build_mat_cov_dens_also_diag()
/// 3. Execute density reconstruction via mp_reconst_density_float() or mp_reconst_density_double()
/// 4. Retrieve results from ReconstResult structure
///
/// @section nagainv_formula Mathematical Background
/// The inversion formula is:
/// @code
/// rho' = rho_0 + (^tA * C_N^(-1) * A + C_rho^(-1))^(-1) * ^tA * C_N^(-1) * (N_obs - N_prior)
/// @endcode
/// where A is the projection matrix, C_N is the muon observation covariance,
/// and C_rho is the prior density covariance.
///
/// @section nagainv_units Units and Coordinate System
/// - Density: kg/m^3
/// - Distance: meters
/// - Coordinate system: inherited from Grid3dVoxel (typically z-up, right-handed)
///
/// @section nagainv_threading Thread Safety
/// - mp_build_mat_cov_dens_also_diag() and mp_build_mat_cov_muon() use OpenMP
/// - reconst_density_core() is const and thread-safe for concurrent reads
/// - Setter functions are NOT thread-safe
#pragma once

#include <map>
#include <fstream>
#include <iostream>
#include <sstream> // istringstream
#include <string>
#include <cstdio>
#include <cmath>
#include <functional>  //for sorting
#include <algorithm>//for sorting
#include <vector>
#include <filesystem>
#include <memory>
#include <nlohmann/json.hpp> // for nlohmann::json

#include <Eigen/Dense>
#include "cls_Grid3dVoxel.hpp"
#include "ns_io_binary.hpp"
#include "cls_NagaInvParameters.hpp"
#include "ns_eigen_blas.hpp"

using namespace index_type_definitions;

/// @brief Regularized linear inversion solver for muon tomography density reconstruction
/// @ingroup matrixClasses
/// @details This class performs MAP (Maximum A Posteriori) estimation for density
/// reconstruction from observed muon counts. It inherits voxel grid structure from
/// Grid3dVoxel and adds inversion-specific data and algorithms.
///
/// @par Class Invariants
/// - mat_dNdD has dimensions (num_observations x num_voxels)
/// - mat_cov_muon is square with size num_observations
/// - mat_cov_dens is square with size num_voxels
/// - vecxf_nmuon_obs and vecxf_nmuon_prior have size num_observations
///
/// @par Example Usage
/// @code
/// NagaInv inv(name, voxel_grid, params, dNdD_matrix, obs_muons, prior_muons);
/// auto result = inv.mp_reconst_density_float(params, prior_density);
/// if (result.is_valid) {
///   // Use result.vecxf_dens_rec for reconstructed density
/// }
/// @endcode
///
/// @see Grid3dVoxel, NagaInvParameters, ReconstResult
class NagaInv : public Grid3dVoxel {
  public:
    /// @brief Maximum voxel count for dense covariance matrix construction
    /// @details mat_cov_dens requires num_voxel^2 * 4 bytes (float).
    ///          50000 voxels -> ~9.3 GB, 51811 voxels -> 10 GB.
    static constexpr int MAX_VOXEL_FOR_DENSE_COV = 50000;

  private:
    /// @brief Observed muon count vector [dimensionless]
    /// @details Size equals number of observation bins (bin groups)
    Eigen::VectorXf vecxf_nmuon_obs = Eigen::VectorXf(0);

    /// @brief Expected muon count vector assuming prior density [dimensionless]
    /// @details Computed from initial density assumption including shell regions
    Eigen::VectorXf vecxf_nmuon_prior = Eigen::VectorXf(0);

    /// @brief Observation matrix A: dN/dDensity (muon count derivative w.r.t. density)
    /// @details Dimensions: (num_observations x num_voxels). Converts density changes to muon count changes.
    Eigen::MatrixXf mat_dNdD = Eigen::MatrixXf::Zero(0,0);

    /// @brief Covariance matrix for muon observations C_N
    /// @details Square matrix (num_observations x num_observations). Diagonal for independent bins.
    Eigen::MatrixXf mat_cov_muon = Eigen::MatrixXf::Zero(0,0);

    /// @brief Prior covariance matrix for density C_rho [kg^2/m^6]
    /// @details Square matrix (num_voxels x num_voxels). Encodes spatial correlation.
    Eigen::MatrixXf mat_cov_dens = Eigen::MatrixXf::Zero(0,0);

  public:
    
    //======================================================================
    /// @name constructor_destructor
    ///@{
    
    /// @brief default constructor
    NagaInv() = default;
    
    /// @brief copy constructor
    NagaInv(const NagaInv &org) = default;

    /// @brief move constructor
    NagaInv(NagaInv &&other) noexcept = default;
    
    /// @brief Constructor from parameters
    /// @param[in] name_in Instance name
    /// @param[in] g3vox_in Target voxel grid
    /// @param[in] prm_nagainv Parameter table for inversion
    /// @param[in] mat_dNdD_in Observation matrix converted to muon counts
    /// @param[in] vecxf_nmuon_obs_in Observed muon counts within defined bin groups
    /// @param[in] vecxf_nmuon_prior_in Expected muon counts assuming initial density (including shell)
    /// @note This function delegates to NagaInv::initialize
    NagaInv(const std::string &name_in
          , const Grid3dVoxel &g3vox_in
          , const NagaInvParameters &prm_nagainv
          , const Eigen::MatrixXf &mat_dNdD_in
          , const Eigen::VectorXf &vecxf_nmuon_obs_in
          , const Eigen::VectorXf &vecxf_nmuon_prior_in)
    {
      initialize(
          name_in
        , g3vox_in
        , prm_nagainv
        , mat_dNdD_in
        , vecxf_nmuon_obs_in
        , vecxf_nmuon_prior_in
        );
    }; 
    
    /// @brief Initialize instance from parameters
    /// @param[in] name_in Instance name
    /// @param[in] g3vox_in Target voxel grid
    /// @param[in] prm_nagainv Parameter table for inversion
    /// @param[in] mat_dNdD Observation matrix converted to muon counts
    /// @param[in] vecxf_nmuon_obs Observed muon counts within defined bin groups
    /// @param[in] vecxf_nmuon_prior Expected muon counts assuming initial density (including shell)
    void initialize(
        const std::string &name_in
      , const Grid3dVoxel &g3vox_in
      , const NagaInvParameters &prm_nagainv
      , const Eigen::MatrixXf &mat_dNdD
      , const Eigen::VectorXf &vecxf_nmuon_obs
      , const Eigen::VectorXf &vecxf_nmuon_prior);

    /// @brief destructor
    ~NagaInv() = default;

    /// @brief Tolerance for equality comparison
    static constexpr float eps_eq = 1.0e-6;

    ///@} ------------------------------------------------------------------
    
    //======================================================================
    /// @name operators
    ///@{
    
    /// @brief Inequality operator (does not compare name or NagaInvParameters)
    /// @param[in] other Object to compare against
    /// @return true if this and other are not equal
    bool operator!=(const NagaInv& other ) const;

    /// @brief Equality operator (does not compare name; defined via inequality operator)
    bool operator==(const NagaInv& other) const {
      return !(*this != other);
    };
    
    /// @brief assignment operator
    NagaInv& operator=(const NagaInv& other) = default;
    ///@} ------------------------------------------------------------------
    
    //======================================================================
    /// @name getter functions
    ///@{

    /// @brief return vecxf_nmuon_obs.rows();
    int get_num_obs() const { return vecxf_nmuon_obs.rows(); };
    
    /// @brief return vecxf_nmuon_obs(i)
    double get_nmuon_obs( const int i ) const;

    /// @brief return vecxf_nmuon_prior(i)
    double get_nmuon_prior( const int i ) const;

    /// @brief  get the copy of mat_dNdD
    Eigen::MatrixXf get_mat_dNdD() const { return mat_dNdD; };

    /// @brief get the reference of mat_dNdD
    const Eigen::MatrixXf& get_mat_dNdD_ref() const { return mat_dNdD; }

    /// @brief get the copy of mat_cov_muon
    Eigen::MatrixXf get_mat_cov_muon() const { return mat_cov_muon; }

    /// @brief get the reference of mat_cov_muon
    const Eigen::MatrixXf& get_mat_cov_muon_ref() const { return mat_cov_muon; }

    /// @brief get the copy of vecxf_nmuon_obs
    Eigen::VectorXf get_vecxf_nmuon_obs() const { return vecxf_nmuon_obs; }

    /// @brief get the reference of vecxf_nmuon_obs
    const Eigen::VectorXf& get_vecxf_nmuon_obs_ref() const { return vecxf_nmuon_obs; }

    /// @brief get the copy of vecxf_nmuon_prior
    Eigen::VectorXf get_vecxf_nmuon_prior() const { return vecxf_nmuon_prior; }
    
    /// @brief get the reference of vecxf_nmuon_prior
    const Eigen::VectorXf& get_vecxf_nmuon_prior_ref() const { return vecxf_nmuon_prior; }

    /// @brief get the copy of mat_cov_dens
    Eigen::MatrixXf get_mat_cov_dens() const { return mat_cov_dens; }

    /// @brief get the reference of mat_cov_dens
    const Eigen::MatrixXf& get_mat_cov_dens_ref() const { return mat_cov_dens; };

    /// @brief get the immutable reference of base class Grid3d
    const Grid3dVoxel& getBaseRef() const { return *this; };

    /// @brief get the copy of Grid3dVoxel
    Grid3dVoxel getBaseCopy() const { return *this; };

    ///@} -------------------------------------------------------------------


    //======================================================================
    /// @name clone functions
    ///@{

    /// @brief Clone with optionally modified data
    /// @param[in] new_dNdD New observation matrix (nullptr to keep unchanged)
    /// @param[in] new_matcov_muon New muon covariance matrix (nullptr to keep unchanged)
    /// @param[in] new_obs New observed muon vector (nullptr to keep unchanged)
    /// @param[in] new_prior New prior muon vector (nullptr to keep unchanged)
    /// @return NagaInv clone with the specified modifications applied
    NagaInv clone_with_modified_data(
      const Eigen::MatrixXf *new_dNdD = nullptr
    , const Eigen::MatrixXf *new_matcov_muon = nullptr
    , const Eigen::VectorXf *new_obs = nullptr
    , const Eigen::VectorXf *new_prior = nullptr ) const; 
    ///@} -------------------------------------------------------------------

    //======================================================================
    /// @name setter functions
    ///@{
    
    /// @brief call mutable reference of base class Grid3d
    Grid3dVoxel& callBaseRef() { return *this; };

    /// @brief set the observed muon vector
    void set_vecxf_nmuon_obs(const Eigen::VectorXf &v) { vecxf_nmuon_obs = v; }

    /// @brief set the prior muon vector
    void set_vecxf_nmuon_prior(const Eigen::VectorXf &v) { vecxf_nmuon_prior = v; }

    /// @brief set the kernel matrix mat_dNdD
    void set_mat_dNdD(const Eigen::MatrixXf &m) { mat_dNdD = m; }

    /// @brief set the covariance matrix mat_cov_muon
    void set_mat_cov_muon(const Eigen::MatrixXf &m) { mat_cov_muon = m; }

    /// @brief set the covariance matrix mat_cov_dens
    void set_mat_cov_dens(const Eigen::MatrixXf &m) { mat_cov_dens = m; }

    ///@} -------------------------------------------------------------------

    //======================================================================
    /// @name build functions
    ///@{

    
    /// @brief Build density covariance matrix (mat_cov_dens) with OpenMP parallelization
    /// @param[in,out] prm_nagainv NagaInv parameter structure (including dynamic function settings)
    /// @note Diagonal elements are computed using sigma_rho_diag.
    /// @note Uses OpenMP for parallel computation. Thread-safe for writing to mat_cov_dens.
    void mp_build_mat_cov_dens_also_diag(NagaInvParameters &prm_nagainv);

    /// @brief Build muon observation covariance matrix (mat_cov_muon) with OpenMP parallelization
    /// @param[in] nmuon_thres Threshold below which observed muon count is clipped
    /// @param[in] nmuon_lt_thres Value to use when nmuon_group < nmuon_thres (prevents division by zero)
    /// @note Uses OpenMP for parallel computation of diagonal elements.
    /// @note Off-diagonal elements are zero (independent observation bins assumed).
    /// @note Variance is N (linear-space Poisson).
    void mp_build_mat_cov_muon(
            const float nmuon_thres = 1.0E-10
          , const float nmuon_lt_thres = 1.0E-10 );

    /// @brief Add efficiency-uncertainty variance to the C_N (mat_cov_muon) diagonal.
    /// @param[in] vecxf_var_eff Per-row efficiency variance, row-aligned with mat_cov_muon.
    ///            An empty vector is a no-op (backward compatible).
    /// @note Added post-construction so initialize()/mp_build_mat_cov_muon() signatures
    ///       stay unchanged. Uses OpenMP for the diagonal update.
    /// @throws std::runtime_error If vecxf_var_eff size does not match mat_cov_muon dimension.
    void add_var_eff_to_cov_muon(const Eigen::VectorXf& vecxf_var_eff);
    ///@} -------------------------------------------------------------------
    
    ///======================================================================
    /// @name output_functions
    ///@{

    /// @brief Output mat_cov_dens to file
    /// @param[in] prefix Filename prefix (output: prefix_name.tmp)
    void out_mat_cov_dens(const std::string &prefix = "mat_cov_dens") const;

    /// @brief Output mat_cov_muon to file
    /// @param[in] prefix Filename prefix (output: prefix_name.tmp)
    void out_mat_cov_muon(const std::string &prefix = "mat_cov_muon") const;

    //======================================================================
    /// @name density reconstruction functions
    ///@{

    public:

/// @brief Structure holding density reconstruction results
/// @details Contains output vectors from density reconstruction, difference information, and validity flag.
    struct ReconstResult {
      /// @brief Reconstructed density vector [kg/m^3]
      Eigen::VectorXf vecxf_dens_rec;

      /// @brief Square root of diagonal elements of posterior density covariance [kg/m^3]
      Eigen::VectorXf vecxf_diag_sqrt_cov_dens;

      /// @brief Difference vector from prior density (rec - prior)
      /// @details Stores the difference between reconstruction and prior density for visualization/evaluation.
      Eigen::VectorXf vecxf_delta_dens_prior;

      /// @brief Difference vector from true density (rec - real)
      /// @details Valid only when true density is provided.
      Eigen::VectorXf vecxf_diff_from_real;

      /// @brief Effective number of parameters p_eff = trace(R), R = C' A^T C_N^-1 A [dimensionless]
      /// @details Computed only when tf_calc_chi2ndf is set (0 otherwise), reusing the gain
      ///          matrix G = C' A^T C_N^-1 already formed for the reconstruction.
      double p_eff = 0.0;

      /// @brief Validity flag for the result
      /// @details Set to false if computation fails or NaN is detected.
      bool is_valid = true;
    };

    private:
    /// @brief Core function for density reconstruction (BLAS/LAPACK-accelerated)
    /// @details Performs MAP estimation-based density reconstruction using prior density,
    /// observations, and covariance matrices. Assumes mat_dNdD and covariance matrices
    /// are pre-built externally.
    /// @n Formula: rho' = rho_0 + (^tA * C_N^(-1) * A + C_rho^(-1))^(-1) * ^tA * C_N^(-1) * (N_obs - N_prior)
    /// @param[in] input_mat_dNdD Projection matrix A (muon counts vs voxels)
    /// @param[in] input_mat_cov_muon Covariance matrix of observed muon counts C_N
    /// @param[in] input_mat_cov_dens Prior covariance matrix of density C_rho
    /// @param[in] input_vecxf_nmuon_obs Observed muon counts N_obs
    /// @param[in] input_vecxf_nmuon_prior Predicted muon counts N_prior
    /// @param[in] input_vecxf_dens_prior Prior voxel density vector (rho_prior) [kg/m^3]
    /// @param[in] prm_nagainv Inversion parameters (regularization, thresholds, etc.)
    /// @param[in] fn_inverse_matrix Function returning inverse matrix (typically BLAS/LAPACK-based)
    /// @param[in] input_real_dens Optional true density; if provided, computes rec - real difference
    /// @param[in] tf_save_tmp_data_bin Reserved for future: save intermediate matrices
    /// @param[in] tf_verbose Enable verbose logging of reconstruction progress
    /// @param[in] precomputed_mat_cov_dens_inv Optional precomputed inverse of input_mat_cov_dens.
    ///            If provided, skips the C_rho^{-1} inversion (useful when the same mat_cov_dens
    ///            is reused across multiple calls, e.g. in NagaInvLooper).
    /// @return ReconstResult containing estimated density, covariance, and differences
    ReconstResult reconst_density_core(
        const Eigen::MatrixXf& input_mat_dNdD
      , const Eigen::MatrixXf& input_mat_cov_muon
      , const Eigen::MatrixXf& input_mat_cov_dens
      , const Eigen::VectorXf& input_vecxf_nmuon_obs
      , const Eigen::VectorXf& input_vecxf_nmuon_prior
      , const Eigen::VectorXf& input_vecxf_dens_prior
      , const NagaInvParameters& prm_nagainv
      , std::function<Eigen::MatrixXf(const Eigen::MatrixXf&)> fn_inverse_matrix
      , const std::optional<Eigen::VectorXf>& input_real_dens = std::nullopt
      , const bool tf_save_tmp_data_bin = false // Reserved for future: save intermediate matrices to binary files
      , const bool tf_verbose = true
      , const std::optional<Eigen::MatrixXf>& precomputed_mat_cov_dens_inv = std::nullopt ) const;

    public:

    /// @brief Execute Nagahara's linear inversion to reconstruct density (float precision)
    /// @details Performs MAP estimation using the formula:
    /// @n rho' = rho_0 + (^tA * C_N^(-1) * A + C_rho^(-1))^(-1) * ^tA * C_N^(-1) * (N_obs - N_prior)
    /// @n Uses BLAS/LAPACK for matrix inversion with float precision.
    /// @param[in] prm_nagainv Inversion parameters (regularization, thresholds)
    /// @param[in] input_vecxf_dens_prior Prior density vector [kg/m^3]
    /// @param[in] input_vecxf_real_dens Optional true density for error evaluation
    /// @param[in] tf_save_tmp_data_bin Reserved for future use
    /// @param[in] tf_verbose Enable verbose logging
    /// @param[in] precomputed_mat_cov_dens_inv Optional precomputed inverse of mat_cov_dens.
    ///            Skips the O(n^3) C_rho inversion when provided.
    /// @return ReconstResult containing reconstructed density and diagnostics
    /// @note Complexity: O(n^3) where n = max(num_observations, num_voxels)
    /// @see mp_reconst_density_double for double precision version
    ReconstResult mp_reconst_density_float(
      const NagaInvParameters &prm_nagainv
    , const Eigen::VectorXf& input_vecxf_dens_prior
    , const std::optional<Eigen::VectorXf>& input_vecxf_real_dens
      = std::nullopt
    , const bool tf_save_tmp_data_bin=false
    , const bool tf_verbose=true
    , const std::optional<Eigen::MatrixXf>& precomputed_mat_cov_dens_inv
      = std::nullopt );

    /// @brief Execute Nagahara's linear inversion to reconstruct density (double precision)
    /// @details Performs MAP estimation using the formula:
    /// @n rho' = rho_0 + (^tA * C_N^(-1) * A + C_rho^(-1))^(-1) * ^tA * C_N^(-1) * (N_obs - N_prior)
    /// @n Uses BLAS/LAPACK for matrix inversion with double precision for improved numerical stability.
    /// @param[in] prm_nagainv Inversion parameters (regularization, thresholds)
    /// @param[in] input_vecxf_dens_prior Prior density vector [kg/m^3]
    /// @param[in] input_vecxf_real_dens Optional true density for error evaluation
    /// @param[in] tf_save_tmp_data_bin Reserved for future use
    /// @param[in] tf_verbose Enable verbose logging
    /// @return ReconstResult containing reconstructed density and diagnostics
    /// @note Complexity: O(n^3) where n = max(num_observations, num_voxels)
    /// @note Slower than float version but more numerically stable for ill-conditioned matrices
    /// @see mp_reconst_density_float for float precision version
    ReconstResult mp_reconst_density_double(
      const NagaInvParameters &prm_nagainv
    , const Eigen::VectorXf& input_vecxf_dens_prior
    , const std::optional<Eigen::VectorXf>& input_vecxf_real_dens
      = std::nullopt
    , const bool tf_save_tmp_data_bin=false
    , const bool tf_verbose=true );

    ///@} -------------------------------------------------------------------

    //======================================================================
    /// @name binary_io_functions
    ///@{

    /// @brief Save NagaInv to std::ofstream
    /// @param[out] ofs Output file stream
    void save( std::ofstream &ofs ) const;

    /// @brief Save NagaInv to std::filesystem::path
    /// @param[in] pathout Output file path
    void save( const std::filesystem::path& pathout ) const;

    /// @brief Load NagaInv from std::ifstream
    /// @param[in,out] ifs Input file stream
    void load( std::ifstream &ifs );

    /// @brief Load NagaInv from std::filesystem::path
    /// @param[in] path_in Input file path
    /// @throws std::runtime_error if architecture is incompatible
    void load( const std::filesystem::path &path_in );
    ///@} -------------------------------------------------------------------

};
