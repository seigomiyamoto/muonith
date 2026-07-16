/// @file cls_NagaInv.cpp
/// @brief Implementation of NagaInv class for muon tomography density reconstruction
#include "cls_NagaInv.hpp"
#include "ns_myapp.hpp"
#include "ns_mymacro.hpp"
#include "ns_mylogger.hpp"
#include "spdlog_pch.hpp"


bool NagaInv::operator!=(const NagaInv& other ) const
{
  #ifdef NODEBUG
    if (!vecxf_nmuon_obs.isApprox(other.vecxf_nmuon_obs, eps_eq)) return true;
    if (!vecxf_nmuon_prior.isApprox(other.vecxf_nmuon_prior, eps_eq)) return true;
    if (!mat_dNdD.isApprox(other.mat_dNdD, eps_eq)) return true;
    if (!mat_cov_muon.isApprox(other.mat_cov_muon, eps_eq)) return true;
    if (!mat_cov_dens.isApprox(other.mat_cov_dens, eps_eq)) return true;
  #else
    if (!vecxf_nmuon_obs.isApprox(other.vecxf_nmuon_obs, eps_eq)) { LOG_WARN("NagaInv: vecxf_nmuon_obs differs"); return true;}
    if (!vecxf_nmuon_prior.isApprox(other.vecxf_nmuon_prior, eps_eq)) { LOG_WARN("NagaInv: vecxf_nmuon_prior differs"); return true;}
    if (!mat_dNdD.isApprox(other.mat_dNdD, eps_eq)) { LOG_WARN("NagaInv: mat_dNdD differs"); return true;}
    if (!mat_cov_muon.isApprox(other.mat_cov_muon, eps_eq)) { LOG_WARN("NagaInv: mat_cov_muon differs"); return true;}
    if (!mat_cov_dens.isApprox(other.mat_cov_dens, eps_eq)) { LOG_WARN("NagaInv: mat_cov_dens differs"); return true;}
  #endif
  return false;
}

// Initialize from parameters
void NagaInv::initialize(
  const std::string &name_in // Instance name
, const Grid3dVoxel &g3vox_in // Target voxel grid
, const NagaInvParameters &prm_nagainv // Parameter table
, const Eigen::MatrixXf &mat_dNdD_in // Observation matrix converted to muon counts
, const Eigen::VectorXf &vecxf_nmuon_obs_in // Observed muon counts in defined bin groups
, const Eigen::VectorXf &vecxf_nmuon_prior_in) // Expected muon count for initial density (including shell)
{
  LOG_INFO(
    "...set parameters from prm_nagainv(name={})"
    ,prm_nagainv.get_name());

  // if tf_exec==false, nothing to do
  if(prm_nagainv.get_tf_exec()==false){
    LOG_WARN("prm_nagainv.tf_exec==false, NagaInv is not executed.");
    return;
  }

  // if Grid3dVoxel::get_uqiv_not_assigned(), THORW_ERROR
  if( g3vox_in.get_uqiv_not_assigned()==g3vox_in.get_uqiv_max() ){
    THROW_ERROR_NAME("g3vox_in.get_uqiv_not_assigned()");
  };

  // set g3vox from g3vox_in
  LOG_INFO(
    "...set g3vox from g3vox_in(name={})"
    ,g3vox_in.get_name());
  this->Grid3dVoxel::set(g3vox_in);

  // set the name
  this->Grid3dVoxel::set_name(name_in);

  // set mat_dNdD from mat_dNdD_in
  LOG_INFO(
    "...set mat_dNdD from mat_dNdD_in");
  this->mat_dNdD = mat_dNdD_in;
  if( mat_dNdD.allFinite()==false )
    THROW_ERROR_NAME("mat_dNdD.allFinite()==false");

  // set vecxf_nmuon from vecxf_nmuon_obs
  LOG_INFO(
    "...set vecxf_nmuon_obs from vecxf_nmuon_obs_in");
  this->vecxf_nmuon_obs = vecxf_nmuon_obs_in;
  if( vecxf_nmuon_obs.allFinite()==false )
    THROW_ERROR_NAME("vecxf_nmuon_obs.allFinite()==false");

  // Set vecxf_nmuon_prior from vecxf_nmuon_prior_in
  LOG_INFO(
    "...set vecxf_nmuon_prior from vecxf_nmuon_prior_in");
  this->vecxf_nmuon_prior = vecxf_nmuon_prior_in;
  if( vecxf_nmuon_prior.allFinite()==false )
    THROW_ERROR_NAME("vecxf_nmuon_prior.allFinite()==false");

  // mp_build mat_cov_muon
  LOG_INFO("...mp_build_mat_cov_muon");
  mp_build_mat_cov_muon(
      prm_nagainv.get_nmuon_thres()
    , prm_nagainv.get_nmuon_under_thres());
  
  if(mat_cov_muon.allFinite()==false)
    THROW_ERROR_NAME("mat_cov_muon.allFinite()==false");

  // build mat_cov_dens
  LOG_INFO("...mp_build_mat_cov_dens");
  NagaInvParameters prm_nagainv_tmp = prm_nagainv.clone();
  mp_build_mat_cov_dens_also_diag(prm_nagainv_tmp);

  if(mat_cov_dens.allFinite()==false)
    THROW_ERROR_NAME("mat_cov_dens.allFinite()==false");
}

/// @brief return vecxf_nmuon_obs(i)
double NagaInv::get_nmuon_obs( const int i ) const
{
  if (i < 0 || i >= vecxf_nmuon_obs.size()) {
    THROW_ERROR_NAME("Index out of bounds in get_nmuon_obs");
  }
  return vecxf_nmuon_obs(i);
}

/// @brief return vecxf_nmuon_prior(i)
double NagaInv::get_nmuon_prior( const int i ) const
{
  if (i < 0 || i >= vecxf_nmuon_prior.size()) {
    THROW_ERROR_NAME("Index out of bounds in get_nmuon_prior");
  }
  return vecxf_nmuon_prior(i);
}

/// @brief Clone with optionally modified data
/// @param new_dNdD New observation matrix (nullptr to keep unchanged)
/// @param new_matcov_muon New muon covariance matrix (nullptr to keep unchanged)
/// @param new_obs New observed muon vector (nullptr to keep unchanged)
/// @param new_prior New prior muon vector (nullptr to keep unchanged)
/// @return NagaInv clone with the specified modifications applied
NagaInv NagaInv::clone_with_modified_data(
  const Eigen::MatrixXf *new_dNdD
, const Eigen::MatrixXf *new_matcov_muon
, const Eigen::VectorXf *new_obs
, const Eigen::VectorXf *new_prior) const
{
  NagaInv copy = *this;
  if (new_dNdD)   copy.mat_dNdD = *new_dNdD;
  if (new_matcov_muon) copy.mat_cov_muon = *new_matcov_muon;
  if (new_obs)   copy.vecxf_nmuon_obs = *new_obs;
  if (new_prior) copy.vecxf_nmuon_prior = *new_prior;
  return copy;
}

void NagaInv::mp_build_mat_cov_dens_also_diag(NagaInvParameters &prm_nagainv)
{
  LOG_INFO("start...");
  prm_nagainv.print_json_params(spdlog::level::debug);

  // Initialize matrix and iterate over voxel pairs
  const int num_voxel = mat_dNdD.cols();

  // Check voxel count to prevent excessive memory usage
  if (num_voxel > MAX_VOXEL_FOR_DENSE_COV) {
    const double estimated_gb = static_cast<double>(num_voxel) * num_voxel * 4.0 / (1024.0 * 1024.0 * 1024.0);
    THROW_ERROR("NagaInv::mp_build_mat_cov_dens_also_diag: "
      "num_voxel({}) exceeds MAX_VOXEL_FOR_DENSE_COV({}). "
      "Estimated memory for mat_cov_dens: {:.1f} GB. "
      "Please reduce voxel count by coarsening grid resolution or narrowing target region.",
      num_voxel, MAX_VOXEL_FOR_DENSE_COV, estimated_gb);
  }

  mat_cov_dens = Eigen::MatrixXf::Zero(num_voxel, num_voxel);

  const int n_threads = omp_get_max_threads();
  int ivox2 = 0;

  const double sigma_in = prm_nagainv.get_sigma_rho();
  const double sigma_sq = sigma_in * sigma_in;
  const bool tf_aniso = prm_nagainv.is_anisotropic();

  // Precompute inverse correlation lengths outside the loop
  const double corr_length_inv = 1.0 / prm_nagainv.get_corr_length();
  const double corr_length_xy_inv = 1.0 / prm_nagainv.get_corr_length_xy();
  const double corr_length_z_inv  = 1.0 / prm_nagainv.get_corr_length_z();
  const bool tf_ellipsoidal = (prm_nagainv.get_aniso_cov_type() == "ellipsoidal");

  if (tf_aniso) {
    LOG_INFO("Anisotropic covariance: corr_length_xy={}, corr_length_z={}, type={}",
      prm_nagainv.get_corr_length_xy(), prm_nagainv.get_corr_length_z(),
      prm_nagainv.get_aniso_cov_type());
  } else {
    LOG_INFO("Isotropic covariance: corr_length={}", prm_nagainv.get_corr_length());
  }

  #pragma omp parallel for private(ivox2)
  for (int ivox1 = 0; ivox1 < num_voxel - 1; ivox1++) {
    for (ivox2 = ivox1 + 1; ivox2 < num_voxel; ivox2++) {
      const Ixiyiz pos0 = get_ixiyiz(ivox1);
      const Ixiyiz pos1 = get_ixiyiz(ivox2);

      double result;
      if (tf_aniso) {
        const double dxy = calc_dist_dxy(pos0, pos1);
        const double dz  = calc_dist_dz(pos0, pos1);
        if (tf_ellipsoidal) {
          // C = sigma^2 * exp(-sqrt((dxy/L_xy)^2 + (dz/L_z)^2))
          const double sxy = dxy * corr_length_xy_inv;
          const double sz  = dz  * corr_length_z_inv;
          result = sigma_sq * exp(-sqrt(sxy * sxy + sz * sz));
        } else {
          // Separable: C = sigma^2 * exp(-dxy/L_xy - dz/L_z)
          result = sigma_sq * exp(-dxy * corr_length_xy_inv - dz * corr_length_z_inv);
        }
      } else {
        // Isotropic: C = sigma^2 * exp(-dxyz/L)
        const double dxyz = calc_dist_dxyz(pos0, pos1);
        result = sigma_sq * exp(-dxyz * corr_length_inv);
      }
      mat_cov_dens(ivox1, ivox2) = static_cast<float>(result);

      if (ivox1 % n_threads == 0 && ivox2 % n_threads == 0) {
        fprintf(stderr,
          "mp_build_mat_cov_dens_also_diag : ivox1=%d / %d, ivox2=%d / %d ...\r"
          ,ivox1, num_voxel - 1, ivox2, num_voxel);
      }
    }
  }

  // Reflect symmetric matrix by adding transpose
  mat_cov_dens += mat_cov_dens.transpose().eval();

  // Fill diagonal elements with self-variance
  double sigma_rho_diag = prm_nagainv.get_sigma_rho_diag();

  if (sigma_rho_diag <= 0.0) {
    LOG_WARN("sigma_rho_diag is non-positive ({}), setting to {}(= sigma_rho)"
    , sigma_rho_diag, prm_nagainv.get_sigma_rho());
    // Use default value
    sigma_rho_diag = prm_nagainv.get_sigma_rho();
  }

  #pragma omp parallel for
  for (int ivox = 0; ivox < num_voxel; ivox++) {
    mat_cov_dens(ivox, ivox) = sigma_rho_diag * sigma_rho_diag;
  }
}

// build mat_cov_number
void NagaInv::mp_build_mat_cov_muon(
  const float nmuon_thres
, const float nmuon_lt_thres )
{
  const int n_data = vecxf_nmuon_obs.rows();
  LOG_INFO("=================================================");
  LOG_INFO("start...");
  LOG_INFO("n_data (=vecxf_nmuon_obs.rows()) = {}",n_data);
  LOG_INFO("=================================================");
  // SLEEP_MSEC(500);
  
  // mat_cov_muon.resize(n_data,n_data);
  mat_cov_muon = Eigen::MatrixXf::Zero(n_data,n_data);
  LOG_INFO("make blank mat_cov_muon({} x {})",n_data,n_data);

  // set diagonal elements
  #pragma omp parallel for
  for(int i=0;i<n_data;i++){
    float nmuon_group = vecxf_nmuon_obs(i);
    if(nmuon_group < nmuon_thres ){
      LOG_WARN("nmuon_group({})={:E} < nmuon_thres={:E}, set nmuon_group = {:E}"
        ,i,nmuon_group,nmuon_thres,nmuon_lt_thres);
        nmuon_group = nmuon_lt_thres;
    }
    double sigma_nmuon_sq = fabs(nmuon_group);
    mat_cov_muon(i,i) =  sigma_nmuon_sq;
  }
}

// Add efficiency-uncertainty variance to the C_N diagonal (post-construction).
// Wiring choice B: keeps initialize()/mp_build_mat_cov_muon() signatures unchanged.
void NagaInv::add_var_eff_to_cov_muon(const Eigen::VectorXf& vecxf_var_eff)
{
  // empty vector -> efficiency C_N diagonal disabled (backward compatible)
  if (vecxf_var_eff.size() == 0) {
    return;
  }

  const int n_data = mat_cov_muon.rows();
  if (vecxf_var_eff.size() != n_data) {
    THROW_ERROR("NagaInv::add_var_eff_to_cov_muon: size mismatch. "
      "vecxf_var_eff.size()={}, mat_cov_muon.rows()={}",
      vecxf_var_eff.size(), n_data);
  }

  LOG_INFO("add efficiency variance to mat_cov_muon diagonal (n_data={})", n_data);
  #pragma omp parallel for
  for (int i = 0; i < n_data; i++) {
    mat_cov_muon(i,i) += vecxf_var_eff(i);
  }

  if (mat_cov_muon.allFinite() == false) {
    THROW_ERROR_NAME("mat_cov_muon.allFinite()==false after adding efficiency variance");
  }
}

NagaInv::ReconstResult NagaInv::mp_reconst_density_float(
  const NagaInvParameters& prm_nagainv
, const Eigen::VectorXf& input_vecxf_dens_prior
, const std::optional<Eigen::VectorXf>& input_vecxf_dens_real
, const bool tf_save_tmp_data_bin
, const bool tf_verbose
, const std::optional<Eigen::MatrixXf>& precomputed_mat_cov_dens_inv )
{
  std::chrono::system_clock::time_point start_time_func, end_time_func;
  start_time_func = time_now;

  const ReconstResult result = reconst_density_core(
      mat_dNdD
    , mat_cov_muon
    , mat_cov_dens
    , vecxf_nmuon_obs
    , vecxf_nmuon_prior
    , input_vecxf_dens_prior
    , prm_nagainv
    , [](const Eigen::MatrixXf& mat) {
          return eigen_blas::getInverseMatrixFloat(mat);
        } // Inverse matrix function for float precision
    , input_vecxf_dens_real // real density vector
    , tf_save_tmp_data_bin
    , tf_verbose // verbose
    , precomputed_mat_cov_dens_inv
  );

  end_time_func = time_now;
  if (tf_verbose) {
    LOG_INFO("reconst_density_core (float) finished.");
    myapp::cast_time_msec(spdlog::level::debug,"(float) execution done."
    , start_time_func, end_time_func);
  }
  
  return result;
}

NagaInv::ReconstResult NagaInv::mp_reconst_density_double(
  const NagaInvParameters& prm_nagainv
, const Eigen::VectorXf& input_vecxf_dens_prior
, const std::optional<Eigen::VectorXf>& input_vecxf_dens_real
, const bool tf_save_tmp_data_bin
, const bool tf_verbose )
{
  std::chrono::system_clock::time_point start_time_func, end_time_func;
  start_time_func = time_now;

  const Eigen::VectorXf vecxf_prior_dens = Grid3dVoxel::get_vecxf_density();
  
  const ReconstResult result = reconst_density_core(
      mat_dNdD
    , mat_cov_muon
    , mat_cov_dens
    , vecxf_nmuon_obs
    , vecxf_nmuon_prior
    , input_vecxf_dens_prior
    , prm_nagainv
    , [](const Eigen::MatrixXf& mat) {
          return eigen_blas::getInverseMatrixDouble(mat);
        } // Inverse matrix function for double precision
    , input_vecxf_dens_real // real density vector
    , tf_save_tmp_data_bin // save tmp bin file
    , tf_verbose // verbose
  );
  end_time_func = time_now;
  if (tf_verbose) {
    LOG_INFO("reconst_density_core (double) finished.");
    myapp::cast_time_msec(spdlog::level::debug,"(double) execution done."
    , start_time_func, end_time_func);
  }

  return result;
}


NagaInv::ReconstResult NagaInv::reconst_density_core(
    const Eigen::MatrixXf& input_mat_dNdD
  , const Eigen::MatrixXf& input_mat_cov_muon
  , const Eigen::MatrixXf& input_mat_cov_dens
  , const Eigen::VectorXf& input_vecxf_nmuon_obs
  , const Eigen::VectorXf& input_vecxf_nmuon_prior
  , const Eigen::VectorXf& input_vecxf_dens_prior
  , const NagaInvParameters& prm_nagainv
  , std::function<Eigen::MatrixXf(const Eigen::MatrixXf&)> fn_inverse_matrix
  , const std::optional<Eigen::VectorXf>& input_vecxf_dens_real
  , const bool tf_save_tmp_data_bin // Reserved for future: save intermediate matrices to binary files
  , const bool tf_verbose
  , const std::optional<Eigen::MatrixXf>& precomputed_mat_cov_dens_inv ) const
{
  std::string msg;

  // 0.0 Check that all input vectors and matrices are finite
  if (!input_mat_dNdD.allFinite()) THROW_ERROR_NAME("input_mat_dNdD invalid");
  if (!input_mat_cov_muon.allFinite()) THROW_ERROR_NAME("input_mat_cov_muon invalid");
  if (!input_mat_cov_dens.allFinite()) THROW_ERROR_NAME("input_mat_cov_dens invalid");
  if (!input_vecxf_nmuon_obs.allFinite()) THROW_ERROR_NAME("input_vecxf_nmuon_obs invalid");
  if (!input_vecxf_nmuon_prior.allFinite()) THROW_ERROR_NAME("input_vecxf_nmuon_prior invalid");
  if (!input_vecxf_dens_prior.allFinite()) THROW_ERROR_NAME("input_vecxf_dens_prior invalid");

  // 0.1.1 check size of vector vs vector
  if (input_vecxf_nmuon_obs.size() != input_vecxf_nmuon_prior.size()) {
    msg = fmt::format(
      "input_vecxf_nmuon_obs.size() != input_vecxf_nmuon_prior.size() : {} != {}"
      , input_vecxf_nmuon_obs.size(), input_vecxf_nmuon_prior.size());
    THROW_ERROR_NAME(msg);
  }

  // 0.1.2 check size of matrix itself
  if (input_mat_cov_dens.rows() != input_mat_cov_dens.cols()) {
    msg = fmt::format(
      "input_mat_cov_dens is not square matrix : {} != {}"
      , input_mat_cov_dens.rows(), input_mat_cov_dens.cols());
    THROW_ERROR_NAME(msg);
  }
  if (input_mat_cov_muon.rows() != input_mat_cov_muon.cols()) {
    msg = fmt::format(
      "input_mat_cov_muon is not square matrix : {} != {}"
      , input_mat_cov_muon.rows(), input_mat_cov_muon.cols());
    THROW_ERROR_NAME(msg);
  }

  // 0.1.3 check size of vector vs matrix
  if( input_mat_cov_muon.rows() != input_vecxf_nmuon_obs.size() ){
    msg = fmt::format(
      "input_mat_cov_muon.rows() != input_vecxf_nmuon_obs.size() : {} != {}"
      , input_mat_cov_muon.rows(), input_vecxf_nmuon_obs.size());
    THROW_ERROR_NAME(msg);
  }
  if( input_mat_cov_muon.cols() != input_vecxf_nmuon_obs.size() ){
    msg = fmt::format(
      "input_mat_cov_muon.cols() != input_vecxf_nmuon_obs.size() : {} != {}"
      , input_mat_cov_muon.cols(), input_vecxf_nmuon_obs.size());
    THROW_ERROR_NAME(msg);
  }

  if (input_mat_dNdD.rows() != input_vecxf_nmuon_obs.size()) {
    msg = fmt::format(
      "input_mat_dNdD.rows() != input_vecxf_nmuon_obs.size() : {} != {}"
      , input_mat_dNdD.rows(), input_vecxf_nmuon_obs.size());
    THROW_ERROR_NAME(msg);
  }
  if( input_mat_dNdD.cols() != input_vecxf_dens_prior.size() ){
    msg = fmt::format(
      "input_mat_dNdD.cols() != input_vecxf_dens_prior.size() : {} != {}"
      , input_mat_dNdD.cols(), input_vecxf_dens_prior.size());
    THROW_ERROR_NAME(msg);
  }

  // 1. Compute inverse matrices
  if (tf_verbose) LOG_INFO("Start inverse mat_cov_muon...");
  Eigen::MatrixXf mat_cov_muon_inv = fn_inverse_matrix(input_mat_cov_muon);
  if (!mat_cov_muon_inv.allFinite()) THROW_ERROR_NAME("mat_cov_muon_inv invalid");

  Eigen::MatrixXf mat_cov_dens_inv;
  if (precomputed_mat_cov_dens_inv.has_value()) {
    if (tf_verbose) LOG_INFO("Using precomputed mat_cov_dens_inv (skipping inversion).");
    mat_cov_dens_inv = *precomputed_mat_cov_dens_inv;
  } else {
    if (tf_verbose) LOG_INFO("Start inverse mat_cov_dens...");
    mat_cov_dens_inv = fn_inverse_matrix(input_mat_cov_dens);
  }
  if (!mat_cov_dens_inv.allFinite()) THROW_ERROR_NAME("mat_cov_dens_inv invalid");

  // 2. (^tA * C_N^(-1))
  if (tf_verbose) LOG_INFO("Start ^tA * C_N^(-1) * A ...");
  Eigen::MatrixXf mat_temp1 = eigen_blas::multiplyMatrixColMajor(
    input_mat_dNdD.transpose(), mat_cov_muon_inv);

  // 3. (^tA * C_N^(-1) * A)
  Eigen::MatrixXf mat_cov_dens_dash_inv = eigen_blas::multiplyMatrixColMajor(
    mat_temp1, input_mat_dNdD);

  // 3.1 resize of mat_temp1(0,0)
  mat_temp1.resize(0, 0);

  // 4. + C_rho^(-1)
  if (tf_verbose) LOG_INFO("Add C_rho^(-1) to (^tA * C_N^(-1) * A) ...");
  mat_cov_dens_dash_inv += mat_cov_dens_inv;

  // 4.1 resize of mat_cov_dens_inv(0,0)
  mat_cov_dens_inv.resize(0, 0);

  // 5. Inverse of the above
  if (tf_verbose) LOG_INFO("Start to calc inverse (^tA * C_N^(-1) * A + C_rho^(-1)) ...");
  Eigen::MatrixXf mat_cov_dens_dash = fn_inverse_matrix(mat_cov_dens_dash_inv);
  if (!mat_cov_dens_dash.allFinite()) THROW_ERROR_NAME("mat_cov_dens_dash invalid");

  // 5.1 resize of mat_cov_dens_dash_inv(0,0)
  mat_cov_dens_dash_inv.resize(0, 0);

  // 5.2 build Grid3dVoxel from mat_cov_dens_dash
  LOG_INFO(" 5.2 save diagonal elements of sqrt of diag elements of mat_cov_dens_dash to vecxf_diag_sqrt_cov_dens_dash ...");
  const Eigen::VectorXf vecxf_diag_sqrt_cov_dens_dash 
    = mat_cov_dens_dash.diagonal().array().abs().sqrt().matrix().eval();

  // 6. mat_all_in_one = mat_cov_dens_dash * ^tA * C_N^(-1)
  if (tf_verbose) LOG_INFO("6. Start to calc mat_all_in_one = mat_cov_dens_dash * ^tA * C_N^(-1) ...");
  Eigen::MatrixXf mat_temp2 = eigen_blas::multiplyMatrixColMajor(
    mat_cov_dens_dash, input_mat_dNdD.transpose());

  Eigen::MatrixXf mat_all_in_one = eigen_blas::multiplyMatrixColMajor(
    mat_temp2, mat_cov_muon_inv);
  if (!mat_all_in_one.allFinite()) THROW_ERROR_NAME("mat_all_in_one invalid");

  // 6.05 p_eff = trace(R) = trace(G A), reusing the gain G = mat_all_in_one = C' A^T C_N^-1.
  // trace(G A) = sum_{i,k} G(i,k) A(k,i): elementwise product of G and A^T, summed in double.
  // Disabled for the public release (bd id-ll39vv): the trace(R) computation is retained in
  // source but not executed, so p_eff stays 0 and never reaches any input or output.
  double d_p_eff = 0.0;
  // if (prm_nagainv.get_tf_calc_chi2ndf()) {
  //   d_p_eff = (mat_all_in_one.array().cast<double>()
  //            * input_mat_dNdD.transpose().array().cast<double>()).sum();
  // }

  // 6.1 resize of mat_temp2(0,0)
  mat_temp2.resize(0, 0);

  // 7. ΔN_obs = vecxf_nmuon_obs - vecxf_nmuon_prior
  Eigen::VectorXf vecxf_delta_nmuon = input_vecxf_nmuon_obs - input_vecxf_nmuon_prior;

  // 8. Δρ = mat_all_in_one * ΔN_obs
  Eigen::VectorXf vecxf_delta_rec_dens
   = eigen_blas::multiplyMatrixColMajor(mat_all_in_one, vecxf_delta_nmuon);

  // 9. Final reconstructed density
  Eigen::VectorXf vecxf_dens_rec = input_vecxf_dens_prior + vecxf_delta_rec_dens;

  ReconstResult result = {
    .vecxf_dens_rec = vecxf_dens_rec,
    .vecxf_diag_sqrt_cov_dens  = vecxf_diag_sqrt_cov_dens_dash
  };

  // Store difference from prior
  result.vecxf_delta_dens_prior = vecxf_dens_rec - input_vecxf_dens_prior;

  // Effective number of parameters (0 unless tf_calc_chi2ndf was set above)
  result.p_eff = d_p_eff;

  // Compare with real density if provided
  if (input_vecxf_dens_real.has_value()) {
    if (input_vecxf_dens_real->size() == vecxf_dens_rec.size()) {
      result.vecxf_diff_from_real = vecxf_dens_rec - *input_vecxf_dens_real;
    } else {
      LOG_WARN("input_vecxf_dens_real size mismatch: expected {}, got {}"
        ,vecxf_dens_rec.size(), input_vecxf_dens_real->size());
    }
  }

  if (tf_verbose)
    LOG_INFO("finished.");
  
  // 10. Validate result
  if (!vecxf_dens_rec.allFinite()) {
    LOG_WARN("In this trial, vecxf_dens_rec is not finite.");
    result.is_valid = false;
    return result;
  }

  return result;
}

// save function to std::ofstream
void NagaInv::save( std::ofstream &ofs ) const
{
  // Grid3dVoxel
  this->Grid3dVoxel::save(ofs);

  // vecxf_nmuon_obs
  io_binary::write_vecxf_stream(ofs,vecxf_nmuon_obs);

  // vecxf_nmuon_prior
  io_binary::write_vecxf_stream(ofs,vecxf_nmuon_prior);

  // mat_dNdD
  io_binary::write_matxf_stream(ofs,mat_dNdD);

  // mat_cov_muon
  io_binary::write_matxf_stream(ofs,mat_cov_muon);

  // mat_cov_dens
  io_binary::write_matxf_stream(ofs,mat_cov_dens);
}

// load from std::ifstream
void NagaInv::load( std::ifstream &ifs )
{
  // Grid3dVoxel
  this->Grid3dVoxel::load(ifs);

  // vecxf_nmuon_obs
  vecxf_nmuon_obs = io_binary::read_vecxf_stream(ifs);

  // vecxf_nmuon_prior
  vecxf_nmuon_prior = io_binary::read_vecxf_stream(ifs);

  // mat_dNdD
  mat_dNdD = io_binary::read_matxf_stream(ifs);

  // mat_cov_muon
  mat_cov_muon = io_binary::read_matxf_stream(ifs);

  // mat_cov_dens
  mat_cov_dens = io_binary::read_matxf_stream(ifs);
}

void NagaInv::save(const std::filesystem::path& pathout) const
{
  std::ofstream ofs = io_binary::open_ofstream(pathout);
  // Write architecture info first
  io_binary::ArchitectureInfo currentInfo = io_binary::get_current_architecture_info();
  io_binary::write_architecture_info(ofs, currentInfo);
  // Then write data
  save(ofs);
  ofs.close();
}

void NagaInv::load(const std::filesystem::path& path_in)
{
  std::ifstream ifs = io_binary::open_ifstream(path_in);
  // 1) Read architecture info stored in file
  io_binary::ArchitectureInfo fileInfo = io_binary::read_architecture_info(ifs);
  // 2) Check compatibility with current environment (throws on mismatch)
  io_binary::check_architecture_compatibility_or_throw(fileInfo);
  // 3) Load data
  load(ifs);
  ifs.close();
}

void NagaInv::out_mat_cov_dens(const std::string& prefix) const
{
  std::filesystem::path pathout = prefix + "_" + get_name() + ".tmp";
  eigen_blas::out_matxf(pathout, mat_cov_dens);
}

void NagaInv::out_mat_cov_muon(const std::string& prefix) const
{
  std::filesystem::path pathout = prefix + "_" + get_name() + ".tmp";
  eigen_blas::out_matxf(pathout, mat_cov_muon);
}

