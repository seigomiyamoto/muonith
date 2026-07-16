// calc_dNdD.cpp
#include <cassert>
#include <iostream>

#include "ns_calc_dNdD.hpp"
#include "ns_pathcalc.hpp"
#include "ns_mymacro.hpp"
#include "ns_mylogger.hpp"
#include "ns_myapp.hpp"
#include <chrono>
#include "cls_MatrixBuildParameters.hpp"
#include "spdlog_pch.hpp"
#include "ns_eigen_blas.hpp"
#include "ns_detector_indexing.hpp"

// make mat_dNdD matrix
// from mat_PL and dFdR_R_costhz in a detector unit
// Create muon number per density matrix from path length matrix.
Eigen::MatrixXf calc_dNdD::mp_make_mat_dNdD(
  const DetectorPanel& panel, const Eigen::MatrixXf &mat_PL
, const Eigen::VectorXf &vecxf_DL_prior, const FluxTable &ft )
{
  const Detid detid = panel.get_detid();
  LOG_DEBUG("detid={} / {}",detid,panel.get_name());

  const int n_ele = mat_PL.rows();
  const int n_vox = mat_PL.cols();
  Eigen::MatrixXf mat_dNdD = Eigen::MatrixXf::Zero(n_ele,n_vox);

  // get number of bin in x and y in the det
  const int nbinx = panel.get_nbinx();
  const int nbiny = panel.get_nbiny();

  // loop of DetectorElement
  const int n_threads = omp_get_max_threads();
  // * Do NOT carelessly use collapse(3) here 2025-03-05 15:48:52
  #pragma omp parallel for collapse(2) schedule(static)
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      const DetectorElement& ele = panel.getDetectorElement(ix,iy);
      const Inthis inthis = ele.get_id_in_this_detector();
      const Eigen::VectorXf vecxf_mat_inthis = mat_PL.row(inthis);
      const float DL_prior = vecxf_DL_prior(inthis);
      const double dFdR = ft.calc_dFdR(ele,DL_prior);
      const double SOT = ele.get_SOT();
      // loop of voxels
      for(int ivox=0;ivox<n_vox;ivox++){
        mat_dNdD(inthis,ivox) = vecxf_mat_inthis(ivox) * dFdR * SOT;
      }
    }
  }
  return mat_dNdD;
}

// make mat_dNdD matrix
// from mat_PL and dFdR_R_costhz in a detector unit
// Create muon number per density matrix from path length matrix.
SpMatf calc_dNdD::mp_make_spmat_dNdD(
  const DetectorPanel& panel, const SpMatf &spmat_PL,
  const Eigen::VectorXf &vecxf_DL_prior, const FluxTable &ft )
{
  const Detid detid = panel.get_detid();
  LOG_DEBUG("detid={} / {}", detid, panel.get_name());

  const int n_ele = spmat_PL.rows();
  const int n_vox = spmat_PL.cols();

  // Calculate scaling factor (dFdR * SOT) for each detector element
  Eigen::VectorXf factors(n_ele);
  const int nbinx = panel.get_nbinx();
  const int nbiny = panel.get_nbiny();
  #pragma omp parallel for collapse(2) schedule(static)
  for (int iy = 0; iy < nbiny; iy++) {
    for (int ix = 0; ix < nbinx; ix++) {
      const DetectorElement& ele = panel.getDetectorElement(ix, iy);
      const Inthis inthis = ele.get_id_in_this_detector();
      const float DL_prior = vecxf_DL_prior(inthis);
      const double dFdR = ft.calc_dFdR(ele, DL_prior);
      const double SOT = ele.get_SOT();
      factors(inthis) = static_cast<float>(dFdR * SOT);
    }
  }

  // Directly scale each row of spmat_PL
  SpMatf spmat_dNdD = spmat_PL;
  #pragma omp parallel for
  for (int k = 0; k < spmat_dNdD.outerSize(); ++k) {
    for (SpMatf::InnerIterator it(spmat_dNdD, k); it; ++it) {
      it.valueRef() *= factors(it.row());
    }
  }

  return spmat_dNdD;
}

// Create function to merge mat_dNdD
// from DetectorPanel and mat_dNdD.
// DetectorPanel holds information on how to merge.
// non openmp version
Eigen::MatrixXf calc_dNdD::make_grouped_mat_dNdD(
    const DetectorPanel& panel, const Eigen::MatrixXf &mat_dNdD )
{
  const Detid detid = panel.get_detid();

  // Get number of rows (= number of groups in g2det)
  const int n_group = panel.get_dic().get_n_group(detid);

  // Get number of columns (= number of voxels)
  const int ncol = mat_dNdD.cols();

  // Define return matrix
  Eigen::MatrixXf mat_dNdD_grouped = Eigen::MatrixXf::Zero(n_group,ncol);

  // loop of column (voxel index)
  for(int i_col=0;i_col<ncol;i_col++){
    if(i_col%10==0) fprintf(stderr,"i_col=%d/%d...\r",i_col,ncol);

    // cut mat_dNdD and make a vector for calc speed up
    const Eigen::VectorXf vecxf_icol = mat_dNdD.col(i_col);

    // bin group loop
    for(Igroup igroup=0;igroup<n_group;igroup++){ // rows raster after merged
      // fprintf(stderr,"igroup=%d/%d...\n",igroup,n_group);
      // const int n_pair = panel.get_n_uqid(igroup);
      mat_dNdD_grouped(igroup,i_col) = 0.0;

      const std::vector<std::array<int,2>> vec_ixiy
         = panel.get_g2bg().get_vec_ixiy(igroup);
      for(const auto& [ix,iy] : vec_ixiy){
        const int inthis = panel.get_id_in_this_detector(ix,iy);
        mat_dNdD_grouped(igroup,i_col) += vecxf_icol(inthis);
      }
    }
  }
  return mat_dNdD_grouped;
}

// Create function to merge mat_dNdD
// from DetectorPanel and mat_dNdD.
// DetectorPanel holds information on how to merge.
// openmp version, but not so fast
Eigen::MatrixXf calc_dNdD::mp_make_grouped_mat_dNdD(
  const DetectorPanel& panel,
  const Eigen::MatrixXf &mat_dNdD )
{
  // Get number of rows after merge
  const Detid detid = panel.get_detid();
  const int n_group = panel.get_dic().get_n_group(detid);
  if (n_group == 0) {
    THROW_ERROR("panel.get_dic().get_n_group(panel.get_detid())==0");
  }

  // Number of columns (= number of voxel elements)
  const int ncol = mat_dNdD.cols();

  const std::string det_name = panel.get_name();

  // Define merged return matrix
  Eigen::MatrixXf mat_dNdD_grouped = Eigen::MatrixXf::Zero(n_group, ncol);

  // Get group information
  const auto vec_igroup_vec_ixiy = panel.get_g2bg().get_vec_igroup_vec_ixiy();

  // 1) Pre-build table converting (ix, iy) in BinGroup to ID
  //    This avoids calling get_id_in_this_detector() inside the loop
  std::vector<std::vector<Inthis>> vec_igroup_vec_inthis(n_group);
  for (Igroup igroup = 0; igroup < n_group; igroup++) {
    const auto &vec_ixiy = vec_igroup_vec_ixiy.at(igroup);
    vec_igroup_vec_inthis.at(igroup).resize(vec_ixiy.size());
    #pragma omp parallel for
    for (size_t i = 0; i < vec_ixiy.size(); i++) {
      const auto [ix, iy] = vec_ixiy.at(i);
      const Inthis inthis = panel.get_id_in_this_detector(ix, iy);
      vec_igroup_vec_inthis.at(igroup).at(i) = inthis;
    }
  }

  // 2) Parallelize outer loop (column i_voxel) with OpenMP
  #pragma omp parallel for
  for (int i_voxel = 0; i_voxel < ncol; i_voxel++) {
    // Reduce progress display frequency
    if (i_voxel % 100 == 0) {
      // Be cautious as multiple threads may compete for display
      fprintf(stderr, "det_name=%s, detid=%02d, i_voxel=%d/%d...\r",
        det_name.c_str(), detid, i_voxel, ncol);
    }

    // Get column vector once
    const Eigen::VectorXf vecxf_matNPL_ivoxel = mat_dNdD.col(i_voxel);

    // 3) Calculate sum for each BinGroup
    // igroup is the row index
    for (Igroup igroup = 0; igroup < n_group; igroup++) {
      float sum = 0.0f;
      // Add based on IDs stored in vec_igroup_vec_inthis[igroup]
      const auto &vec_inthis = vec_igroup_vec_inthis.at(igroup);

      // SIMD can be used if needed (syntax varies by environment)
      // #pragma omp simd reduction(+: sum)
      for (const Inthis& inthis : vec_inthis) {
        sum += vecxf_matNPL_ivoxel(inthis);
      }
      mat_dNdD_grouped(igroup, i_voxel) = sum;
    }
  }

  return mat_dNdD_grouped;
}

Eigen::MatrixXf calc_dNdD::mp_make_grouped_mat_dNdD(
  const DetectorPanel& panel, const SpMatf &spmat_dNdD )
{
  const std::string det_name = panel.get_name();
  LOG_DEBUG("det_name={} / {}", det_name, panel.get_detid());

  const Detid detid = panel.get_detid();
  const int n_group = panel.get_dic().get_n_group(detid); // number of merged angle-bin groups
  if (n_group == 0) {
    THROW_ERROR("panel.get_dic().get_n_group(detid)==0");
  }

  const int ncol = spmat_dNdD.cols(); // = number of voxel elements
  const int nrow = spmat_dNdD.rows(); // = number of detector elements (before grouping)

  Eigen::MatrixXf mat_dNdD_grouped = Eigen::MatrixXf::Zero(n_group, ncol); // (n_group × ncol) output

  const auto vec_igroup_vec_ixiy // [igroup] → vector of (ix,iy) pixels in that group
   = panel.get_g2bg().get_vec_igroup_vec_ixiy();

  // reverse lookup: row index (inthis) → group index (igroup), -1 = ungrouped
  std::vector<Igroup> row_to_group(nrow, -1);
  for (Igroup igroup = 0; igroup < n_group; igroup++) {
    const auto &vec_ixiy = vec_igroup_vec_ixiy.at(igroup);
    for (size_t i = 0; i < vec_ixiy.size(); i++) {
      const auto [ix, iy] = vec_ixiy.at(i);
      const Inthis inthis = panel.get_id_in_this_detector(ix, iy); // (ix,iy) → flat row index
      if (inthis >= 0 && inthis < nrow) {
        row_to_group[inthis] = igroup; // register which group this row belongs to
      }
    }
  }

  #pragma omp parallel for // each thread writes to different columns — no data race
  for (int i_voxel = 0; i_voxel < ncol; i_voxel++) {
    if (i_voxel % 100 == 0) {
      fprintf(stderr, "det_name=%s, detid=%02d, i_voxel=%d/%d...\r",
        det_name.c_str(), detid, i_voxel, ncol);
    }

    // InnerIterator: walks only non-zero entries in this column (ColMajor sparse)
    for (SpMatf::InnerIterator it(spmat_dNdD, i_voxel); it; ++it) {
      const Igroup ig = row_to_group[it.row()]; // O(1) flat array lookup
      if (ig >= 0) { // skip rows not belonging to any group
        mat_dNdD_grouped(ig, i_voxel) += it.value(); // accumulate into group
      }
    }
  }

  return mat_dNdD_grouped;
}

// Legacy version: converts sparse columns to dense VectorXf before grouping.
// Retained for regression testing.
Eigen::MatrixXf calc_dNdD::mp_make_grouped_mat_dNdD_legacy(
  const DetectorPanel& panel, const SpMatf &spmat_dNdD )
{
  const std::string det_name = panel.get_name();
  LOG_DEBUG("det_name={} / {}", det_name, panel.get_detid());

  const Detid detid = panel.get_detid();
  const int n_group = panel.get_dic().get_n_group(detid);
  if (n_group == 0) {
    THROW_ERROR("panel.get_dic().get_n_group(detid)==0");
  }

  const int ncol = spmat_dNdD.cols();

  Eigen::MatrixXf mat_dNdD_grouped = Eigen::MatrixXf::Zero(n_group, ncol);

  const auto vec_igroup_vec_ixiy
   = panel.get_g2bg().get_vec_igroup_vec_ixiy();

  std::vector<std::vector<Inthis>> vec_igroup_vec_inthis(n_group);
  for (Igroup igroup = 0; igroup < n_group; igroup++) {
    const auto &vec_ixiy = vec_igroup_vec_ixiy.at(igroup);
    vec_igroup_vec_inthis.at(igroup).resize(vec_ixiy.size());
    #pragma omp parallel for
    for ( size_t i = 0; i < vec_ixiy.size(); i++ ) {
      const auto [ix, iy] = vec_ixiy.at(i);
      const Inthis inthis = panel.get_id_in_this_detector(ix, iy);
      vec_igroup_vec_inthis.at(igroup).at(i) = inthis;
    }
  }

  #pragma omp parallel for
  for (int i_voxel = 0; i_voxel < ncol; i_voxel++) {
    if (i_voxel % 100 == 0) {
      fprintf(stderr, "det_name=%s, detid=%02d, i_voxel=%d/%d...\r",
        det_name.c_str(), detid, i_voxel, ncol);
    }

    const Eigen::VectorXf vecxf_matNPL_ivoxel = spmat_dNdD.col(i_voxel);

    for (Igroup igroup = 0; igroup < n_group; igroup++) {
      float sum = 0.0f;
      const auto &vec_inthis = vec_igroup_vec_inthis.at(igroup);

      for (const Inthis inthis : vec_inthis) {
        sum += vecxf_matNPL_ivoxel(inthis);
      }
      mat_dNdD_grouped(igroup, i_voxel) = sum;
    }
  }

  return mat_dNdD_grouped;
}

// Calculate mat_PL for g3vox for all detectors and compute expected signal from initial density.
// Return value is vector of mat_PL for each detector.
// Workflow for each detector:
// 1. Execute g3vox::mp_make_mat_PL(det)
// 2. Execute panel.mp_calc_set_peneflux_signal_from_DL(ft) to calculate signal amount
// 3. Execute panel.copy_signal_noise_to_g2bg() to initialize BinGroup
// 4. Execute panel.auto_divide_by_signal_noise_group_all(panel.get_prm_bingrp()) to optimize BinGroup
// 5. Execute panel.calc_set_vec_signal_noise_group() to calculate signal/noise_poi_group and set to vec_signal/noise_poi_group
// 6. Store calculated mat_PL in vector
// 7. Repeat above for all detectors
std::vector<Eigen::MatrixXf> calc_dNdD::calc_mat_PL_and_set_signal(
    DetectorPanelArray &arrdet, Grid3dVoxel &g3vox
  , const FluxTable &ft
  , const pathcalc::Parameters &prm
  , const bool tf_apply_eff )
{
  // Eigen openmp
  const int n_threads = omp_get_max_threads();
  // myapp::set_threads_Eigen(n_threads);

  // Get number of columns = number of g3vox(true)
  // const int ncol = g3vox.get_n_vox_exist();

  // get number of detector
  const int n_detector = arrdet.get_n_det();

  // Vector for storage
  std::vector<Eigen::MatrixXf> vec_spmat_PL;
  vec_spmat_PL.reserve(n_detector);

  // initilization of matrix to be returned
  Eigen::MatrixXf mat_dNdD_grouped_alldet(0,0);

  // Time measurement variables
  std::chrono::system_clock::time_point start, end;

  // Start time measurement
  start = time_now;

  // Loop for each detector
  for(Detid detid=0;detid<n_detector;detid++){
    LOG_INFO("---------------------------------------------------------------");
    LOG_INFO(" calc_mat_PL_and_set_signal detid={}",detid);
    LOG_INFO("---------------------------------------------------------------");
    // Call each detector
    DetectorPanel& panel = arrdet.callDetectorPanel(detid);

    // Create path length matrix for each detector
    // Whether to increment n_hit for voxels hit by beam
    // and whether to update path length of detector elements is specified by prm.
    Eigen::MatrixXf mat_PL
      = pathcalc::g3vox::mp_make_mat_PL(panel,g3vox,prm);

    // Calculate signal amount
    panel.mp_calc_set_peneflux_signal_from_DL(ft,tf_apply_eff);

    // Initialize BinGroup and set signal to vec_vec_signal.
    panel.copy_signal_noise_to_g2bg();

    // Optimize BinGroup
    panel.call_g2bg().auto_divide_by_signal_noise_group_all(panel.get_prm_bingrp());

    // Calculate signal/noise_poi_group and set to vec_signal/noise_poi_group.
    panel.call_g2bg().calc_set_vec_signal_noise_group();

    vec_spmat_PL.push_back(mat_PL);
  }
  // Stop time measurement
  end = time_now;

  // Display calculation time
  myapp::cast_time_msec(stderr,"calc_mat_PL_and_set_signal finished",start,end);

  return vec_spmat_PL;
}


// For all detectors, perform make_grouped_mat_dNdD
// and create the merged matrix.
// 1. Call g2det pointer for each detector.
// 2. Get mat_PL for the detector using Eigen::MatrixXf& mat_PL = vec_spmat_PL.at(detid);.
// 3. Execute calc_dNdD::mp_make_mat_dNdD(panel,mat_PL,ft)
//    to create muon number per density matrix.
// 4. Create number per density matrix for that detector,
Eigen::MatrixXf
  calc_dNdD::make_grouped_alldet_mat_dNdD(
      DetectorPanelArray &arrdet, Grid3dVoxel &g3vox
    , const std::vector<SpMatf> &vec_spmat_PL
    , const std::vector<Eigen::VectorXf> &vec_vecxf_DL_prior
    , const FluxTable &ft )
{
  // get number of detector
  const int n_detector = arrdet.get_n_det();
  LOG_DEBUG(
    "arrdet.get_n_det()={}"
    , n_detector);

  // initilization of matrix to be returned
  Eigen::MatrixXf mat_dNdD_grouped_alldet(0,0);

  std::chrono::system_clock::time_point start, end;
  start = time_now;

  char fnameout_tmp[512];

  // Process for each detector
  for(Detid detid=0;detid<n_detector;detid++){
    LOG_INFO("---------------------------------------------------------------");
    LOG_INFO(" make_grouped_alldet_mat_dNdD detid={}",detid);
    LOG_INFO("---------------------------------------------------------------");

    // Call each detector
    DetectorPanel& panel = arrdet.callDetectorPanel(detid);

    Eigen::MatrixXf mat_PL = eigen_blas::convertSparseToDense(vec_spmat_PL.at(detid));

    // get the reference of vecxf_DL_prior in detid
    const Eigen::VectorXf vecxf_DL_prior = vec_vecxf_DL_prior.at(detid);

    // make mat_dNdD matrix
    // from mat_PL and logdFdR_table in a detector unit
    start = time_now;
    Eigen::MatrixXf mat_dNdD
      = mp_make_mat_dNdD(panel,mat_PL,vecxf_DL_prior,ft);
    end = time_now;
    std::string msg ="calc_dNdD::mp_make_mat_dNdD finished in detid="
                     + std::to_string(detid);
    myapp::cast_time_msec(spdlog::level::debug,msg,start,end);

    // memory release of mat_PL
    mat_PL.resize(0,0);

    // make_merged mat_dNdD
    start = time_now;
    Eigen::MatrixXf mat_dNdD_grouped
        = mp_make_grouped_mat_dNdD(panel,mat_dNdD);
    end = time_now;
    msg = "mp_make_grouped_mat_dNdD finished in detid="
          + std::to_string(detid);
    myapp::cast_time_msec(spdlog::level::info,msg,start,end);
    LOG_INFO("detid={} mat_dNdD_grouped.rows()={}, cols()={}",detid,mat_dNdD_grouped.rows(),mat_dNdD_grouped.cols());
sprintf(fnameout_tmp,"mat_dNdD_grouped_det%02d.tmp",detid);
eigen_blas::out_matrix_non_zero(fnameout_tmp,mat_dNdD_grouped);

    // memory release of mat_dNdD
    mat_dNdD.resize(0,0);

    //======================================================================
    // combine the merged matrix
    //======================================================================
    eigen_blas::append_matrix_rows(mat_dNdD_grouped_alldet,mat_dNdD_grouped);
    LOG_DEBUG("detid={}, mat_dNdD_grouped_alldet.rows()={}, cols()={}",
    detid,mat_dNdD_grouped_alldet.rows(),mat_dNdD_grouped_alldet.cols());

    // memory release of mat_dNdD_grouped
    mat_dNdD_grouped.resize(0,0);

    // info about mat_dNdD_grouped_alldet
    LOG_INFO("detid={}, mat_dNdD_grouped_alldet.rows()={}, cols()={}",detid,mat_dNdD_grouped_alldet.rows(),mat_dNdD_grouped_alldet.cols());

    // SLEEP_MSEC(500);
  }

  // From irow=uqig, remove those that do not exist in map_uqig_avail_uqig.
  mat_dNdD_grouped_alldet = arrdet.get_matrix_uqig_avail_from_matrix_uqig(mat_dNdD_grouped_alldet);
  LOG_INFO("after remove unavailable uqig, mat_dNdD_grouped_alldet.rows()={}, cols()={}",mat_dNdD_grouped_alldet.rows(),mat_dNdD_grouped_alldet.cols());
  SLEEP_MSEC(500);

sprintf(fnameout_tmp,"mat_dNdD_grouped_alldet.tmp");
eigen_blas::out_matrix_non_zero(fnameout_tmp,mat_dNdD_grouped_alldet);
  return mat_dNdD_grouped_alldet;
}

// For all detectors, perform make_grouped_mat_dNdD
// and create the merged matrix.
// 1. Call g2det pointer for each detector.
// 2. Get mat_PL for the detector using Eigen::MatrixXf& mat_PL = vec_spmat_PL.at(detid);.
// 3. Execute calc_dNdD::mp_make_mat_dNdD(panel,mat_PL,ft)
//    to create muon number per density matrix.
// 4. Create number per density matrix for that detector,
Eigen::MatrixXf
  calc_dNdD::make_grouped_alldet_mat_dNdD_sprs(
      DetectorPanelArray &arrdet, Grid3dVoxel &g3vox
    , const std::vector<SpMatf> &vec_spmat_PL
    , const std::vector<Eigen::VectorXf> &vec_vecxf_DL_prior
    , const FluxTable &ft
    , const bool tf_apply_eff
    , const bool tf_apply_eff_cnt )
{
  arrdet.display_status();

  // get number of detector
  const int n_detector = arrdet.get_n_det();
  LOG_DEBUG(
    "arrdet.get_n_det()={}"
    , n_detector);

LOG_DEBUG("vec_spmat_PL.size()={}", vec_spmat_PL.size());
LOG_DEBUG("vec_vecxf_DL_prior.size()={}", vec_vecxf_DL_prior.size());

if (static_cast<size_t>(n_detector) != vec_spmat_PL.size())
  THROW_ERROR("size mismatch: n_detector={} vs vec_spmat_PL.size()={}"
  , n_detector, vec_spmat_PL.size());

if (static_cast<size_t>(n_detector) != vec_vecxf_DL_prior.size())
  THROW_ERROR("size mismatch: n_detector={} vs vec_vecxf_DL_prior.size()={}"
  , n_detector, vec_vecxf_DL_prior.size());

  // initilization of matrix to be returned
  Eigen::MatrixXf mat_dNdD_grouped_alldet(0,0);

  std::chrono::system_clock::time_point start, end;
  start = time_now;

  // Process for each detector
  for(Detid detid=0;detid<n_detector;detid++){
    LOG_INFO("---------------------------------------------------------------");
    LOG_INFO(" make_grouped_alldet_mat_dNdD_sprs detid={}",detid);
    LOG_INFO("---------------------------------------------------------------");

    // Call each detector
    DetectorPanel& panel = arrdet.callDetectorPanel(detid);

    const SpMatf& spmat_PL = vec_spmat_PL.at(detid);

    // get the reference of vecxf_DL_prior in detid
    const Eigen::VectorXf& vecxf_DL_prior = vec_vecxf_DL_prior.at(detid);

    // make spmat_dNdD matrix
    // from mat_PL and logdFdR_table in a detector unit
    start = time_now;
    SpMatf spmat_dNdD
      = mp_make_spmat_dNdD(panel,spmat_PL,vecxf_DL_prior,ft);
    end = time_now;
    std::string msg ="calc_dNdD::mp_make_spmat_dNdD(sprs version) finished in detid="
                     + std::to_string(detid);
    myapp::cast_time_msec(spdlog::level::debug,msg,start,end);

    start = time_now;
    Eigen::MatrixXf mat_dNdD_grouped;
    if (tf_apply_eff_cnt) {
      // central deterministic efficiency: scale A by eff_cnt before grouping (not the dice path)
      const SpMatf spmat_eff_cnt = panel.get_spmat_eff_cnt();
      const SpMatf spmat_dNdD_eff = spmat_eff_cnt * spmat_dNdD;
      mat_dNdD_grouped = mp_make_grouped_mat_dNdD(panel, spmat_dNdD_eff);
    } else if (tf_apply_eff) {
      // apply efficiency table in the detector
      const SpMatf mat_eff_sample = panel.get_spmat_eff_sample();
      const SpMatf spmat_dNdD_eff = mat_eff_sample * spmat_dNdD;
      mat_dNdD_grouped = mp_make_grouped_mat_dNdD(panel, spmat_dNdD_eff);
    } else {
      // effciency = 100%
      mat_dNdD_grouped = mp_make_grouped_mat_dNdD(panel, spmat_dNdD);
    }

    end = time_now;
    msg = "mp_make_grouped_mat_dNdD(sprs version) finished in detid="
          + std::to_string(detid);
    myapp::cast_time_msec(spdlog::level::info,msg,start,end);
    LOG_INFO("detid={} mat_dNdD_grouped.rows()={}, cols()={}",detid,mat_dNdD_grouped.rows(),mat_dNdD_grouped.cols());

    //======================================================================
    // combine the merged matrix
    //======================================================================
    eigen_blas::append_matrix_rows(mat_dNdD_grouped_alldet,mat_dNdD_grouped);
    LOG_DEBUG("detid={}, mat_dNdD_grouped_alldet.rows()={}, cols()={}",
    detid,mat_dNdD_grouped_alldet.rows(),mat_dNdD_grouped_alldet.cols());

    // memory release of mat_dNdD_grouped
    mat_dNdD_grouped.resize(0,0);

    // info about mat_dNdD_grouped_alldet
    LOG_INFO("detid={}, mat_dNdD_grouped_alldet.rows()={}, cols()={}",detid,mat_dNdD_grouped_alldet.rows(),mat_dNdD_grouped_alldet.cols());

    // SLEEP_MSEC(500);
  }

  LOG_INFO("before remove unavailable uqig, mat_dNdD_grouped_alldet.rows()={}, cols()={}"
    , mat_dNdD_grouped_alldet.rows(), mat_dNdD_grouped_alldet.cols());

  // From irow=uqig, remove those that do not exist in map_uqig_avail_uqig.
  mat_dNdD_grouped_alldet
    = arrdet.get_matrix_uqig_avail_from_matrix_uqig(mat_dNdD_grouped_alldet);

  LOG_INFO("after remove unavailable uqig, mat_dNdD_grouped_alldet.rows()={}, cols()={}"
    , mat_dNdD_grouped_alldet.rows(), mat_dNdD_grouped_alldet.cols());
  SLEEP_MSEC(500);

  return mat_dNdD_grouped_alldet;
}

// Factory function to create numbered path length matrix containing information for all detectors and voxels
Eigen::MatrixXf calc_dNdD::create_grouped_mat_dNdD_alldet(
  const pathcalc::MatrixBuildParameters& prm_mat )
{
  std::string msg;

  // final matrix to be returned
  Eigen::MatrixXf mat_grouped_alldet_dNdD = Eigen::MatrixXf::Zero(0, 0);

  // Check if file exists
  const bool tf_load_file_found = std::filesystem::exists(prm_mat.path_bin_obs_mat_dNdD);
  LOG_DEBUG("tf_load_file_found={}", tf_load_file_found);

  //======================================================================
  // Load without calculation.
  //======================================================================
  if (prm_mat.tf_load_bin_obs_mat_dNdD) {
    if (tf_load_file_found) {
      LOG_INFO("mat_grouped_alldet_dNdD is loaded from {}.",
        prm_mat.path_bin_obs_mat_dNdD.string());

      mat_grouped_alldet_dNdD = io_binary::read_matxf_bin(prm_mat.path_bin_obs_mat_dNdD);

      LOG_DEBUG("mat_grouped_alldet_dNdD.rows()={}", mat_grouped_alldet_dNdD.rows());
      LOG_DEBUG("mat_grouped_alldet_dNdD.cols()={}", mat_grouped_alldet_dNdD.cols());

      // Check size of mat_grouped_alldet_dNdD
      if (mat_grouped_alldet_dNdD.rows() == 0)
        THROW_ERROR("mat_grouped_alldet_dNdD.rows()==0");
      if (mat_grouped_alldet_dNdD.cols() == 0)
        THROW_ERROR("mat_grouped_alldet_dNdD.cols()==0");
      if (mat_grouped_alldet_dNdD.data() == nullptr)
        THROW_ERROR("mat_grouped_alldet_dNdD.data()==nullptr");

      return mat_grouped_alldet_dNdD;
    } else {
      // tf_load_bin_obs_mat_dNdD == true but file not found
      LOG_WARN("tf_load_bin_obs_mat_dNdD = true but file {} not found. Proceeding with calculation.",
        prm_mat.path_bin_obs_mat_dNdD.string());
    }
  }

  //======================================================================
  // Calculate (tf_load_bin_obs_mat_dNdD == false or file not found)
  //======================================================================
  LOG_INFO("create_grouped_mat_dNdD_alldet make_grouped_alldet_mat_dNdD start.");
  mat_grouped_alldet_dNdD
    = calc_dNdD::make_grouped_alldet_mat_dNdD(
        prm_mat.arrdet,
        prm_mat.g3vox,
        prm_mat.vec_spmat_PL,
        prm_mat.vec_vecxf_DL_prior,
        prm_mat.ft
      );

  LOG_DEBUG("mat_grouped_alldet_dNdD.rows()={}", mat_grouped_alldet_dNdD.rows());
  LOG_DEBUG("mat_grouped_alldet_dNdD.cols()={}", mat_grouped_alldet_dNdD.cols());
  msg = "make_grouped_alldet_mat_dNdD done.";
  LOG_INFO(msg);

  if (prm_mat.tf_save_bin_obs_mat_dNdD) {
    LOG_INFO("mat_grouped_alldet_dNdD is saved to {}.",
      prm_mat.path_bin_obs_mat_dNdD.string());
    io_binary::out_matxf_bin(prm_mat.path_bin_obs_mat_dNdD, mat_grouped_alldet_dNdD);
  }

  return mat_grouped_alldet_dNdD;
}

// Factory function to create numbered path length matrix containing information for all detectors and voxels
Eigen::MatrixXf
  calc_dNdD::create_grouped_mat_dNdD_alldet_sprs(
    const pathcalc::MatrixBuildParameters& prm_mat )
{
  std::string msg;

  // final matrix to be returned
  Eigen::MatrixXf mat_grouped_alldet_dNdD = Eigen::MatrixXf::Zero(0, 0);

  // Check if file exists
  const bool tf_load_file_found = std::filesystem::exists(prm_mat.path_bin_obs_mat_dNdD);
  LOG_DEBUG("tf_load_file_found={}", tf_load_file_found);

  //======================================================================
  // Load without calculation.
  //======================================================================
  if (prm_mat.tf_load_bin_obs_mat_dNdD) {
    if (tf_load_file_found) {
      LOG_INFO("mat_grouped_alldet_dNdD is loaded from {}.",
        prm_mat.path_bin_obs_mat_dNdD.string());

      mat_grouped_alldet_dNdD = io_binary::read_matxf_bin(prm_mat.path_bin_obs_mat_dNdD);

      LOG_DEBUG("mat_grouped_alldet_dNdD.rows()={}", mat_grouped_alldet_dNdD.rows());
      LOG_DEBUG("mat_grouped_alldet_dNdD.cols()={}", mat_grouped_alldet_dNdD.cols());

      // Check size of mat_grouped_alldet_dNdD
      if (mat_grouped_alldet_dNdD.rows() == 0)
        THROW_ERROR("mat_grouped_alldet_dNdD.rows()==0");
      if (mat_grouped_alldet_dNdD.cols() == 0)
        THROW_ERROR("mat_grouped_alldet_dNdD.cols()==0");
      if (mat_grouped_alldet_dNdD.data() == nullptr)
        THROW_ERROR("mat_grouped_alldet_dNdD.data()==nullptr");

      return mat_grouped_alldet_dNdD;
    } else {
      // tf_load_bin_obs_mat_dNdD == true but file not found
      LOG_WARN("tf_load_bin_obs_mat_dNdD = true but file {} not found. Proceeding with calculation.",
        prm_mat.path_bin_obs_mat_dNdD.string());
    }
  }

  //======================================================================
  // Calculate (tf_load_bin_obs_mat_dNdD == false or file not found)
  //======================================================================
  LOG_INFO("create_grouped_mat_dNdD_alldet make_grouped_alldet_mat_dNdD start.");
  mat_grouped_alldet_dNdD
    = calc_dNdD::make_grouped_alldet_mat_dNdD_sprs(
        prm_mat.arrdet
      , prm_mat.g3vox
      , prm_mat.vec_spmat_PL
      , prm_mat.vec_vecxf_DL_prior
      , prm_mat.ft
      , prm_mat.tf_apply_eff
      , prm_mat.tf_apply_eff_cnt
      );

  LOG_DEBUG("mat_grouped_alldet_dNdD.rows()={}", mat_grouped_alldet_dNdD.rows());
  LOG_DEBUG("mat_grouped_alldet_dNdD.cols()={}", mat_grouped_alldet_dNdD.cols());
  msg = "make_grouped_alldet_mat_dNdD done.";
  LOG_INFO(msg);

  if (prm_mat.tf_save_bin_obs_mat_dNdD) {
    LOG_INFO("mat_grouped_alldet_dNdD is saved to {}.",
      prm_mat.path_bin_obs_mat_dNdD.string());
    io_binary::out_matxf_bin(prm_mat.path_bin_obs_mat_dNdD, mat_grouped_alldet_dNdD);
  }

  return mat_grouped_alldet_dNdD;
}
