// pathcalc.cpp
#include <cassert>
#include <iostream>

// 2022-09-07 15:03:20
#include "ns_pathcalc.hpp"
#include "ns_mymacro.hpp"
#include "ns_mylogger.hpp" 
#include "ns_myapp.hpp" 
#include <chrono>
#include "cls_MatrixBuildParameters.hpp"
#include "spdlog_pch.hpp"
#include "ns_mylogger.hpp"
#include "ns_eigen_blas.hpp"
#include "ns_detector_indexing.hpp"

// calc path length and density length of DetectorElement&
void pathcalc::g3vox::add_PLDL(
  DetectorElement &ele, const Grid3dVoxel &g3vox, const double eps )
{
  const Ray3d ray3d = ele.get_ray3d();

  double sum_PL = 0.0;
  double sum_dens_length = 0.0;

  // single-pass DDA with inline path-length (opt #1+#3)
  g3vox.traverse_ray_with_pathlength(ray3d, eps,
    [&](const Grid3d::Ixiyiz& ixiyiz_hit_vox, double delta_path) {
      const Voxel& vox = g3vox.getVoxel(ixiyiz_hit_vox);
      if (!vox.get_tf_exist()) return;

      const double dens = vox.get_density();
      if (!std::isfinite(delta_path)) {
        THROW_ERROR("g3vox::add_PLDL: delta_path is non-finite. delta_path={}", delta_path);
      }
      if (!std::isfinite(dens)) {
        THROW_ERROR("g3vox::add_PLDL: density is non-finite. dens={}", dens);
      }

      sum_PL += delta_path;
      sum_dens_length += delta_path * dens;
    });

  ele.add_PL(sum_PL);
  ele.add_DL(sum_dens_length);
}

// calc path length and density length of DetectorPanel&
void pathcalc::g3vox::mp_add_PLDL(
  DetectorPanel& panel, const Grid3dVoxel &g3vox, const double eps )
{
  // get the number of DetectorElement in DetectorPanel
  const int nbinx = panel.get_nbinx();
  const int nbiny = panel.get_nbiny();
  #pragma omp parallel for collapse(2) schedule(static)
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      DetectorElement& ele = panel.callDetectorElement(ix,iy);
      add_PLDL(ele,g3vox, eps);
    }
  }
}

// calc path length and density length of DetectorPanelArray&
void pathcalc::g3vox::mp_add_PLDL(
  DetectorPanelArray &arrdet, const Grid3dVoxel &g3vox, const double eps )
{
  // get the number of DetectorPanel in DetectorPanelArray
  const int n_detector = arrdet.get_n_det();
  for(Detid detid=0;detid<n_detector;detid++){
    LOG_INFO("detid={}",detid);
    DetectorPanel& panel = arrdet.callDetectorPanel(detid);
    mp_add_PLDL(panel,g3vox, eps);
  }
}

// BL_max overload: add_PLDL for DetectorElement
// BL_max check requires a separate ray-AABB intersection, so we use
// the optimized get_hit_voxel_index (which does single intersection for BL_max)
// combined with traverse_ray_with_pathlength for the path-length computation.
void pathcalc::g3vox::add_PLDL(
  DetectorElement &ele, const Grid3dVoxel &g3vox, const double BL_max, const double eps )
{
  const Ray3d ray3d = ele.get_ray3d();

  // BL_max warning check via the optimized single-intersection path (opt #4)
  if (BL_max > 0.0) {
    auto [tf_hit, tmin_bl, tmax_bl] = ray3d.is_intersect(g3vox.get_cached_aabb3d());
    if (tf_hit && tmax_bl > BL_max) {
      LOG_WARN_ND("Ray beam length ({:.1f}m) exceeds BL_max ({:.1f}m). Ray3d={}",
                  tmax_bl, BL_max, ray3d.to_string());
    }
  }

  double sum_PL = 0.0;
  double sum_dens_length = 0.0;

  // single-pass DDA with inline path-length (opt #1+#3)
  g3vox.traverse_ray_with_pathlength(ray3d, eps,
    [&](const Grid3d::Ixiyiz& ixiyiz_hit_vox, double delta_path) {
      const Voxel& vox = g3vox.getVoxel(ixiyiz_hit_vox);
      if (!vox.get_tf_exist()) return;

      const double dens = vox.get_density();
      if (!std::isfinite(delta_path)) {
        THROW_ERROR("g3vox::add_PLDL: delta_path is non-finite. delta_path={}", delta_path);
      }
      if (!std::isfinite(dens)) {
        THROW_ERROR("g3vox::add_PLDL: density is non-finite. dens={}", dens);
      }

      sum_PL += delta_path;
      sum_dens_length += delta_path * dens;
    });

  ele.add_PL(sum_PL);
  ele.add_DL(sum_dens_length);
}

// BL_max overload: mp_add_PLDL for DetectorPanel
void pathcalc::g3vox::mp_add_PLDL(
  DetectorPanel& panel, const Grid3dVoxel &g3vox, const double BL_max, const double eps )
{
  const int nbinx = panel.get_nbinx();
  const int nbiny = panel.get_nbiny();
  #pragma omp parallel for collapse(2) schedule(static)
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      DetectorElement& ele = panel.callDetectorElement(ix,iy);
      add_PLDL(ele, g3vox, BL_max, eps);
    }
  }
}

// BL_max overload: mp_add_PLDL for DetectorPanelArray
void pathcalc::g3vox::mp_add_PLDL(
  DetectorPanelArray &arrdet, const Grid3dVoxel &g3vox, const double BL_max, const double eps )
{
  const int n_detector = arrdet.get_n_det();
  for(Detid detid=0; detid<n_detector; detid++){
    LOG_INFO("detid={}",detid);
    DetectorPanel& panel = arrdet.callDetectorPanel(detid);
    mp_add_PLDL(panel, g3vox, BL_max, eps);
  }
}

// it is called in mp_make_mat_PL
// Create path length matrix for each DetectorElement.
// Increments n_hit_ele and n_hit_det for voxels hit by the beam.
// Stores path length per voxel in mat_PL and returns it.
void pathcalc::g3vox::calc_mat_PL(
  DetectorElement &ele, Grid3dVoxel &g3vox, Eigen::MatrixXf &mat_PL
  , const pathcalc::Parameters &prm, const double eps )
{
  // Get ID within g2det
  const int id_this_det = ele.get_id_in_this_detector();
  
  // get detid
  const Detid detid = ele.get_detid();

  // No inside-grid guard: a detector may legitimately sit outside g3vox
  // (consistent with mp_make_spmat_PL); ray traversal handles outside origins.


  const Ray3d ray3d = ele.get_ray3d();

  double sum_PL = 0.0;
  double sum_dens_length = 0.0;

  // single-pass DDA with inline path-length (opt #1+#3) and flat uqiv (opt #2)
  g3vox.traverse_ray_with_pathlength(ray3d, eps,
    [&](const Grid3d::Ixiyiz& ixiyiz_hit_vox, double delta_path) {
      Voxel& vox = g3vox.callVoxel(ixiyiz_hit_vox);
      if (!vox.get_tf_exist()) return;

      const int uqiv = g3vox.get_uqiv_fast(ixiyiz_hit_vox);
      if (uqiv == Grid3d::UqivNotFound) return;

      mat_PL(id_this_det, uqiv) = delta_path;

      sum_PL += delta_path;
      sum_dens_length += delta_path * vox.get_density();

      if (prm.tf_incr_nhit_ele) g3vox.record_hit_ele(uqiv);
      if (prm.tf_incr_nhit_det) g3vox.record_hit_det(uqiv, detid);
    });
  
  // add path length and density length in det
  if( prm.tf_add_PLDL==true ){
    ele.add_PL(sum_PL);
    ele.add_DL(sum_dens_length);
  }
}

// called in mp_make_mat_PL,
// make_grouped_alldet_mat_dNdD_non_mp_merged,
// make_grouped_alldet_mat_dNdD,
// mp_make_grouped_alldet_mat_dNdD
// Create path length matrix for each detector.
// Increments n_hit_ele and n_hit_det for voxels hit by the beam.
// Stores path length per voxel in the matrix and returns it.
// openmp version
// Internally calculates path length matrix,
// which simplifies the calculation of g3vox.incr_n_hit_det.
Eigen::MatrixXf pathcalc::g3vox::mp_make_mat_PL(
  DetectorPanel& panel, Grid3dVoxel &g3vox
  , const pathcalc::Parameters &prm, const double eps )
{
  // Set number of threads used by Eigen library
  const int n_threads = omp_get_max_threads();
  Eigen::setNbThreads(n_threads);
  
  // Get total number of detector elements
  const int n_ele = panel.get_n_element();
  LOG_DEBUG("n_ele={}",n_ele);

  // Number of voxels where rock exists
  const int n_voxel_exist = g3vox.get_n_vox_exist();
  LOG_DEBUG("n_voxel_exist={}",n_voxel_exist);

  // for debug and memory check
  long int n_maxrix_elements = (long int)n_ele*(long int)n_voxel_exist;
  LOG_DEBUG("n_maxrix_elements={:.1E}",(double)n_maxrix_elements);

  // Initialize matrix mat_PL to store calculation results
  // memory allocation of Eigen::MatrixXf mat_PL
  Eigen::MatrixXf mat_PL = Eigen::MatrixXf::Zero( n_ele, n_voxel_exist );
  LOG_DEBUG("mat_PL.rows()={}, mat_PL.cols()={}",mat_PL.rows(),mat_PL.cols());

  // enable multi-thread
  #pragma omp parallel for collapse(2) schedule(static)
  // all detector element loop start
  // Process each detector element in 2D detector grid
  for(int iy=0;iy<panel.get_nbiny(); iy++ ){
    for(int ix=0;ix<panel.get_nbinx(); ix++ ){
      // call the pointer of detector element
      DetectorElement& ele = panel.callDetectorElement(ix,iy);

      // calc path length matrix for a detector element
      calc_mat_PL(ele,g3vox,mat_PL,prm, eps);

    } // detector xloop end
  } //detector yloop end

  // detector element index
  const int nrow = mat_PL.rows();
  const int nele = panel.get_n_element();
  if(nrow!=nele) THROW_ERROR3("nrow!=nele",nrow,nele);

  return mat_PL;
}

// Sparse matrix version of pathcalc::g3vox::mp_make_mat_PL
SpMatf pathcalc::g3vox::mp_make_spmat_PL(
  DetectorPanel& panel, Grid3dVoxel &g3vox,
  const pathcalc::Parameters &prm, const double eps )
{

  const int n_threads = omp_get_max_threads();
  Eigen::setNbThreads(n_threads);

  const int n_ele = panel.get_n_element();
  const int n_voxel_exist = g3vox.get_n_vox_exist();

  // Triplet list for each thread
  std::vector<std::vector<Eigen::Triplet<float>>> triplets_per_thread(n_threads);

  // Parallelization
  #pragma omp parallel for collapse(2) schedule(static)
  for (int iy = 0; iy < panel.get_nbiny(); iy++) {
    for (int ix = 0; ix < panel.get_nbinx(); ix++) {
      int tid = omp_get_thread_num();
      auto &local_triplets = triplets_per_thread[tid];

      DetectorElement& ele = panel.callDetectorElement(ix, iy);
      const int id_this_det = ele.get_id_in_this_detector();
      const Ray3d ray3d = ele.get_ray3d();

      double sum_PL = 0.0;
      double sum_dens_length = 0.0;

      // Single-pass DDA with inline path-length (opt #1+#3) and flat uqiv (opt #2)
      const Detid detid_this = ele.get_detid();
      g3vox.traverse_ray_with_pathlength(ray3d, eps,
        [&](const Grid3d::Ixiyiz& ixiyiz_hit_vox, double delta_path) {
          Voxel& vox = g3vox.callVoxel(ixiyiz_hit_vox);
          if (!vox.get_tf_exist()) return;

          // flat-array O(1) lookup replaces unordered_map hash search (opt #2)
          const int uqiv = g3vox.get_uqiv_fast(ixiyiz_hit_vox);
          if (uqiv == Grid3d::UqivNotFound) return;

          local_triplets.push_back(Eigen::Triplet<float>(
              id_this_det, uqiv, static_cast<float>(delta_path)));

          sum_PL += delta_path;
          sum_dens_length += delta_path * vox.get_density();

          if (prm.tf_incr_nhit_ele)
            g3vox.record_hit_ele(uqiv);
          if (prm.tf_incr_nhit_det)
            g3vox.record_hit_det(uqiv, detid_this);
        });

      if (prm.tf_add_PLDL) {
        ele.add_PL(sum_PL);
        ele.add_DL(sum_dens_length);
      }
    }
  }

  // Merge triplets
  std::vector<Eigen::Triplet<float>> all_triplets;
  for (const auto &local_triplets : triplets_per_thread) {
    all_triplets.insert(all_triplets.end(), local_triplets.begin(), local_triplets.end());
  }

  // Generate SpMatf
  SpMatf spmat_PL(n_ele, n_voxel_exist);
  spmat_PL.setFromTriplets(all_triplets.begin(), all_triplets.end());

  return spmat_PL;
}

// mp_make_mat_PL for DetectorPanelArray
std::vector<SpMatf> pathcalc::g3vox::mp_make_mat_PL(
    DetectorPanelArray &arrdet, Grid3dVoxel &g3vox
  , const pathcalc::Parameters &prm, const double eps )
{
  // Eigen OpenMP
  const int n_threads = omp_get_max_threads();
  
  // check zmin and z_det
  const bool tf_above = is_all_detector_above_zmin(arrdet,g3vox);
  if( tf_above==false ) THROW_ERROR("tf_above==false");
  

  // get num of dectecor
  const int n_det = arrdet.get_n_det();

  // Vector for storage
  std::vector<SpMatf> vec_spmat_PL;

  // Time measurement variables
  std::chrono::system_clock::time_point start, end;

  // Start time measurement
  start = time_now;

  Eigen::MatrixXf mat_PL = Eigen::MatrixXf::Zero(0,0);

  // DetectorPanel loop

  // ! DO NOT use "#pragma omp parallel for" here. Memory deallocation seems to cause issues.
  for( Detid detid=0; detid < n_det; detid++ ){
    LOG_INFO("pathcalc::g3vox::mp_make_mat_PL detid={} / {}....",detid,n_det);
    // call the pointer of DetectorPanel
    DetectorPanel& panel = arrdet.callDetectorPanel(detid);
    LOG_DEBUG("detid={} / {}",detid,panel.get_name());

    LOG_DEBUG("mp_make_mat_PL for DetectorPanel detid={} started",detid);
    mat_PL = mp_make_mat_PL(panel,g3vox,prm, eps);
    LOG_DEBUG("mp_make_mat_PL for DetectorPanel detid={} finished",detid);


    // Convert Eigen::MatrixXf to SpMatf
    LOG_DEBUG("convert Eigen::MatrixXf to SpMatf");
    SpMatf spmat_PL = eigen_blas::convertDenseToSparse(mat_PL,prm.reference_matPL_sparse,prm.epsilon_matPL_sparse);
    LOG_DEBUG("conversion finished");

    // Add to vec_spmat_PL
    vec_spmat_PL.push_back(spmat_PL);

    // Free memory of mat_PL
    mat_PL.resize(0,0);

  } // end of detector element loop
  // Stop time measurement
  end = time_now;

  // Display calculation time
  std::string msg = "pathcalc::g3vox::mp_make_mat_PL finished";
  myapp::cast_time_msec(spdlog::level::debug,msg,start,end);

  return vec_spmat_PL;
}

// mp_make_vec_spmat_PL for DetectorPanelArray
std::vector<SpMatf> pathcalc::g3vox::mp_make_vec_spmat_PL(
    DetectorPanelArray &arrdet, Grid3dVoxel &g3vox
  , const pathcalc::Parameters &prm, const double eps ) 
{
  // Eigen openmp
  const int n_threads = omp_get_max_threads();
  // myapp::set_threads_Eigen(n_threads);
  
  // check zmin and z_det
  // ? is it necessary to check here again ?
  // const bool tf_above = is_all_detector_above_zmin(arrdet,g3vox);
  // if( tf_above==false ) THROW_ERROR("tf_above==false");

  // get num of dectecor
  const int n_det = arrdet.get_n_det();

  // Vector for storage
  std::vector<SpMatf> vec_spmat_PL;

  // Time measurement variables
  std::chrono::system_clock::time_point start, end;

  // Start time measurement
  start = time_now;

  // DetectorPanel loop

  // ! DO NOT use "#pragma omp parallel for" here. Memory deallocation seems to cause issues.
  for( Detid detid=0; detid < n_det; detid++ ){
    LOG_INFO("pathcalc::g3vox::mp_make_mat_PL detid={} / {}....",detid,n_det);
    // call the pointer of DetectorPanel
    DetectorPanel& panel = arrdet.callDetectorPanel(detid);
    LOG_DEBUG("detid={} / {}",detid,panel.get_name());

    LOG_DEBUG("mp_make_mat_PL for DetectorPanel detid={} started",detid);
    SpMatf spmat_PL = mp_make_spmat_PL(panel,g3vox,prm, eps);
    LOG_DEBUG("mp_make_mat_PL for DetectorPanel detid={} finished",detid);

    // Add to vec_spmat_PL
    vec_spmat_PL.push_back(spmat_PL);

  } // end of detector element loop
  // Stop time measurement
  end = time_now;

  // Display calculation time
  std::string msg = "pathcalc::g3vox::mp_make_vec_spmat_PL finished";
  myapp::cast_time_msec(spdlog::level::debug,msg,start,end);

  arrdet.set(FlgProg::tf_calc_PL_DL, true);

  return vec_spmat_PL;
}

/// @brief calc and add DL for VerticalEllipticCylinderCapped to DetectorElement&
void pathcalc::vcyl::add_DL( DetectorElement& ele
  , const VerticalEllipticCylinderCapped &vcyl )
{
  const double PL = geom_util::path_length_inside( vcyl, ele.get_ray3d() );
  const double DL = PL * vcyl.density;
  ele.add_DL(DL);
}

/// @brief calc and add DL for VerticalEllipticCylinderCapped to DetectorPanel&
void pathcalc::vcyl::mp_add_DL( DetectorPanel& panel
  , const VerticalEllipticCylinderCapped &vcyl )
{
  const int nbinx = panel.get_nbinx();
  const int nbiny = panel.get_nbiny();
  #pragma omp parallel for collapse(2) schedule(static)
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      DetectorElement& ele = panel.callDetectorElement(ix,iy);
      add_DL(ele,vcyl);
    }
  }
}

/// @brief calc and add DL for VerticalEllipticCylinderCapped to DetectorPanelArray&
void pathcalc::vcyl::mp_add_DL( DetectorPanelArray &arrdet
, const VerticalEllipticCylinderCapped &vcyl )
{
  const int n_det = arrdet.get_n_det();
  for(Detid detid=0;detid<n_det;detid++){
    LOG_INFO("detid={}",detid);
    DetectorPanel& panel = arrdet.callDetectorPanel(detid);
    mp_add_DL(panel,vcyl);
  }
}



/// @brief calc and add PL/DL for VerticalEllipticCylinderCapped to DetectorElement&
void pathcalc::vcyl::add_PLDL( DetectorElement& ele
  , const VerticalEllipticCylinderCapped &vcyl )
{
  const double PL = geom_util::path_length_inside( vcyl, ele.get_ray3d() );
  const double DL = PL * vcyl.density;
  ele.add_PL(PL);
  ele.add_DL(DL);
}

/// @brief calc and add PL/DL for VerticalEllipticCylinderCapped to DetectorPanel&
void pathcalc::vcyl::mp_add_PLDL( DetectorPanel& panel
  , const VerticalEllipticCylinderCapped &vcyl )
{
  const int nbinx = panel.get_nbinx();
  const int nbiny = panel.get_nbiny();
  #pragma omp parallel for collapse(2) schedule(static)
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      DetectorElement& ele = panel.callDetectorElement(ix,iy);
      add_PLDL(ele,vcyl);
    }
  }
}

/// @brief calc and add PL/DL for VerticalEllipticCylinderCapped to DetectorPanelArray&
void pathcalc::vcyl::mp_add_PLDL( DetectorPanelArray &arrdet
  , const VerticalEllipticCylinderCapped &vcyl )
{
  const int n_det = arrdet.get_n_det();
  for(Detid detid=0;detid<n_det;detid++){
    LOG_INFO("detid={}",detid);
    DetectorPanel& panel = arrdet.callDetectorPanel(detid);
    mp_add_PLDL(panel,vcyl);
  }
}

// NOTE: pathcalc::matrix functions have been moved to the calc_dNdD namespace.
// See ns_calc_dNdD.cpp for implementations.

// add_PLDL
// return the size of the cube
int pathcalc::g2pil::add_PLDL(
  DetectorElement& ele, const Grid2dPillar &g2pil, const double eps )
{
  // Initialize density length and path length
  double sum_PL = 0.0; // meters
  double sum_DL = 0.0; // kg/m2

  // Get detector position and direction in x,y space
  const Ray2d ray2d = ele.get_ray2d();

  // Get detector position and direction in x,y,z space
  const Ray3d ray3d = ele.get_ray3d();

  // Single-pass DDA traversal (optimization #3: no per-ray vector allocation)
  int n_hit = 0;
  g2pil.traverse_ray_2d(ray2d, eps,
    [&](const Grid2d::Ixiy& ixiy) {
      // calc path length via 3D ray-AABB intersection
      const double delta_path = g2pil.get_delta_path(ixiy, ray3d);

      // DEBUG: check for NaN sources
      const double dens = g2pil.get_density(ixiy);
      if (!std::isfinite(delta_path)) {
        THROW_ERROR("g2pil::add_PLDL: delta_path is non-finite. delta_path={}", delta_path);
      }
      if (!std::isfinite(dens)) {
        THROW_ERROR("g2pil::add_PLDL: density is non-finite. dens={}", dens);
      }

      // sum path length and density length
      sum_PL += delta_path;
      sum_DL += delta_path * dens;
      ++n_hit;
    });

  // Set cumulative path length and density length to detector
  ele.add_PL(sum_PL);
  ele.add_DL(sum_DL);
  return n_hit;
}

/// @brief add_PLDL with beam length warning (optimization #4/#5: single AABB intersection)
/// @details If BL_max > 0 and beam length exceeds BL_max, a warning is logged
int pathcalc::g2pil::add_PLDL(
  DetectorElement& ele, const Grid2dPillar &g2pil, const double BL_max, const double eps )
{
  // Initialize density length and path length
  double sum_PL = 0.0; // meters
  double sum_DL = 0.0; // kg/m2

  // Get detector position and direction in x,y space
  const Ray2d ray2d = ele.get_ray2d();

  // Get detector position and direction in x,y,z space
  const Ray3d ray3d = ele.get_ray3d();

  // BL_max check via single intersection with cached bounds (optimization #4)
  if (BL_max > 0.0) {
    auto [tf_hit, tmin_bl, tmax_bl] = ray2d.is_intersect(g2pil.get_cached_aabb2d());
    if (tf_hit && tmax_bl > BL_max) {
      LOG_WARN_ND("Ray beam length ({:.1f}m) exceeds BL_max ({:.1f}m). Ray2d={}",
                  tmax_bl, BL_max, ray2d.to_string());
    }
  }

  // Single-pass DDA traversal (optimization #3: no per-ray vector allocation)
  int n_hit = 0;
  g2pil.traverse_ray_2d(ray2d, eps,
    [&](const Grid2d::Ixiy& ixiy) {
      // calc path length via 3D ray-AABB intersection
      const double delta_path = g2pil.get_delta_path(ixiy, ray3d);

      // DEBUG: check for NaN sources
      const double dens = g2pil.get_density(ixiy);
      if (!std::isfinite(delta_path)) {
        THROW_ERROR("g2pil::add_PLDL: delta_path is non-finite. delta_path={}", delta_path);
      }
      if (!std::isfinite(dens)) {
        THROW_ERROR("g2pil::add_PLDL: density is non-finite. dens={}", dens);
      }

      // sum path length and density length
      sum_PL += delta_path;
      sum_DL += delta_path * dens;
      ++n_hit;
    });

  // Set cumulative path length and density length to detector
  ele.add_PL(sum_PL);
  ele.add_DL(sum_DL);
  return n_hit;
}

/// @brief add_PLDL and change path length and density length of *det with vec_tf_in_PL
/// @details Uses traverse_ray_2d visitor (optimization #3: no per-ray vector allocation)
///          and single ray-AABB intersection per cuboid (optimization #1 from 2D-1).
int pathcalc::g2pil::add_PLDL_with_vec_tf_in_PL(
  DetectorElement& ele, const Grid2dPillar &g2pil, const int ixiy_dist_thres, const double eps )
{
  // Initialize density length and path length
  double sum_PL = 0.0; // meters
  double sum_DL = 0.0; // kg/m2

  // Get detector position and direction in x,y space
  const Ray2d ray2d = ele.get_ray2d();

  // Get detector position and direction in x,y,z space
  const Ray3d ray3d = ele.get_ray3d();

  // State variables for tf_in tracking (same logic as before, but inside visitor)
  bool tf_in = false;
  bool tf_in_prev = false;
  bool tf_hit_prev = false;
  double delta_path_prev = -1.0;
  Ixiy ixiy_prev = Ixiy{-1,-1};
  double dist_min_prev = -1.0;
  double dist_max_prev = -1.0;
  int n_hit = 0;

  // Single-pass DDA traversal (optimization #3: no per-ray vector allocation)
  g2pil.traverse_ray_2d(ray2d, eps,
    [&](const Grid2d::Ixiy& ixiy) {
      // Single ray-AABB intersection (optimization 2D-1)
      const AABB3d aabb3d = g2pil.get_AABB3d(ixiy);
      auto [tf_hit,dist_min,dist_max] = ray3d.is_intersect(aabb3d);

      // Compute delta_path from intersection result (replaces get_delta_path call)
      const double delta_path = tf_hit ? (dist_max - std::max(0.0, dist_min)) : 0.0;

      // sum path length and density length
      sum_PL += delta_path;
      sum_DL += delta_path * g2pil.get_density(ixiy);

      // make dist_min/max positive
      dist_min = std::max(0.0,dist_min);
      dist_max = std::max(0.0,dist_max);

      // Initialize prev variables only on first iteration
      if (n_hit == 0) {
        dist_min_prev = dist_min;
        dist_max_prev = dist_max;
        delta_path_prev = delta_path;
      }

      const bool tf_adjacent = g2pil.is_adjacent(ixiy,ixiy_prev,ixiy_dist_thres);

      // If tf_in was false (outside material) in previous loop and tf_hit==true (delta_path>0), entered material
      if (!tf_in_prev && tf_hit) {
        tf_in = true; // Set to inside material
        ele.insert_tf_in_PL(tf_in,dist_min_prev);
      }
      // If tf_in was true (inside material) in previous loop and ixiy and ixiy_prev are not adjacent, exited material
      if (tf_in_prev && ((!tf_hit && tf_adjacent) || !tf_adjacent)) {
        tf_in = false; // Set to outside material
        ele.insert_tf_in_PL(tf_in, dist_max_prev);
      }

      // update prev values
      tf_in_prev = tf_in;
      tf_hit_prev = tf_hit;
      ixiy_prev = ixiy;
      delta_path_prev = delta_path;
      dist_min_prev = dist_min;
      dist_max_prev = dist_max;
      ++n_hit;
    });

  // If inside material at last loop, set to outside material
  if (tf_in_prev) {
    tf_in = false; // Set to outside material
    ele.insert_tf_in_PL(tf_in, dist_max_prev);
  }

  // check alternating tf_in and the last tf_in==false
  ele.check_alternating_tf_in();

  // Set cumulative path length and density length to detector
  ele.add_PL(sum_PL);
  ele.add_DL(sum_DL);
  return n_hit;
}


// 2022-09-30 17:54:01
// add_PLDL
// return the inedx pair of the cube
std::vector<std::tuple<int,int,double>>
pathcalc::g2pil::add_PLDL_with_vec_tp( 
  DetectorElement& ele, const Grid2dPillar &g2pil, const double eps )
{
  // Initialize density length and path length
  double sum_PL = 0.0; // meters
  double sum_DL = 0.0; // kg/m2

  // Get detector position and direction in x,y space
  const Ray2d ray2d = ele.get_ray2d();

  // Get detector position and direction in x,y,z space
  const Ray3d ray3d = ele.get_ray3d();
  
  std::vector<Ixiy> vec_ixiy_hitcub = g2pil.get_hit_boxes_index(ray2d, eps);

  std::vector<std::tuple<int,int,double>> vec_ix_iy_PL;
  vec_ix_iy_PL.reserve(vec_ixiy_hitcub.size());

  // hit voxel index loop
  for(const auto& ixiy : vec_ixiy_hitcub ){
    // calc path length
    const double delta_path = g2pil.get_delta_path(ixiy,ray3d);

    // sum path length and density length
    sum_PL += delta_path;
    sum_DL += delta_path * g2pil.get_density(ixiy);

    const auto [ix,iy] = ixiy;
    vec_ix_iy_PL.push_back(std::make_tuple(ix,iy,delta_path));
  }
  // ele.set_PL(sum_PL);
  // ele.set_DL(sum_DL);
  ele.add_PL(sum_PL);
  ele.add_DL(sum_DL);
  return vec_ix_iy_PL;
}

// it calls int 
//  pathcalc::debug::add_PLDL(
// DetectorElement& ele, const Grid2dPillar &g2pil );
void pathcalc::g2pil::mp_add_PLDL(
  DetectorPanel& panel, const Grid2dPillar &g2pil, const double eps )
{
  // for debug
  panel.get_x_axis().out_info(spdlog::level::debug);
  panel.get_y_axis().out_info(spdlog::level::debug);

  // Get detector ID
  const Detid detid = panel.get_detid();
  LOG_DEBUG("detid={}",detid);


  // Get number of bins (nbinx, nbiny) in 2D detector grid
  const int nbinx = panel.get_nbinx();
  const int nbiny = panel.get_nbiny();

  // Parallelize for loop
  #pragma omp parallel for collapse(2) schedule(static)

  // Process each detector element in 2D detector grid
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      // Get detector element information
      DetectorElement& ele = panel.callDetectorElement(ix,iy);
      Nhit n_hit_cub = add_PLDL(ele,g2pil, eps);

      // Display information
      if( iy%20==0 && ix%20==0 ) fprintf(stderr
        , "detid=%d, n_cub_hit=%4d: %s: ax=%7.4lf ay=%7.4lf, PL=%8.1lf, DL=%3.1E\r"
        , detid, n_hit_cub, ele.get_angle_unit().c_str()
        , ele.get_tx(),ele.get_ty(),ele.get_PL(), ele.get_DL());
    }
  }
}

// mp_add_PLDL for DetectorPanelArray
void pathcalc::g2pil::mp_add_PLDL(
  DetectorPanelArray &arrdet, const Grid2dPillar &g2pil, const double eps )
{
  std::chrono::system_clock::time_point start, end;
  start = time_now;
  if( arrdet.get(FlgProg::tf_calc_PL_DL) ){
    LOG_WARN("arrdet.get(FlgProg::tf_calc_PL_DL)==true");
    SLEEP_MSEC(500);
  }
  for( Detid detid =0; detid < arrdet.get_n_det(); detid++ ){
    DetectorPanel& panel = arrdet.callDetectorPanel(detid);
    mp_add_PLDL(panel,g2pil, eps);
  }
  arrdet.set(FlgProg::tf_calc_PL_DL, true);
  end = time_now;
  myapp::cast_time_msec(spdlog::level::info,"pathcalc::g2pil::mp_add_PLDL(arrdet) finished",start,end);
}

/// @brief mp_add_PLDL with beam length warning for DetectorPanel
void pathcalc::g2pil::mp_add_PLDL(
  DetectorPanel& panel, const Grid2dPillar &g2pil, const double BL_max, const double eps )
{
  // for debug
  panel.get_x_axis().out_info(spdlog::level::debug);
  panel.get_y_axis().out_info(spdlog::level::debug);

  // Get detector ID
  const Detid detid = panel.get_detid();
  LOG_DEBUG("detid={}, BL_max={}", detid, BL_max);


  // Get number of bins (nbinx, nbiny) in 2D detector grid
  const int nbinx = panel.get_nbinx();
  const int nbiny = panel.get_nbiny();

  // Parallelize for loop
  #pragma omp parallel for collapse(2) schedule(static)

  // Process each detector element in 2D detector grid
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      // Get detector element information
      DetectorElement& ele = panel.callDetectorElement(ix,iy);
      Nhit n_hit_cub = add_PLDL(ele, g2pil, BL_max, eps);

      // Display information
      if( iy%20==0 && ix%20==0 ) fprintf(stderr
        , "detid=%d, n_cub_hit=%4d: %s: ax=%7.4lf ay=%7.4lf, PL=%8.1lf, DL=%3.1E\r"
        , detid, n_hit_cub, ele.get_angle_unit().c_str()
        , ele.get_tx(),ele.get_ty(),ele.get_PL(), ele.get_DL());
    }
  }
}

/// @brief mp_add_PLDL with beam length warning for DetectorPanelArray
void pathcalc::g2pil::mp_add_PLDL(
  DetectorPanelArray &arrdet, const Grid2dPillar &g2pil, const double BL_max, const double eps )
{
  std::chrono::system_clock::time_point start, end;
  start = time_now;
  if( arrdet.get(FlgProg::tf_calc_PL_DL) ){
    LOG_WARN("arrdet.get(FlgProg::tf_calc_PL_DL)==true");
    SLEEP_MSEC(500);
  }
  for( Detid detid =0; detid < arrdet.get_n_det(); detid++ ){
    DetectorPanel& panel = arrdet.callDetectorPanel(detid);
    mp_add_PLDL(panel, g2pil, BL_max, eps);
  }
  arrdet.set(FlgProg::tf_calc_PL_DL, true);
  end = time_now;
  myapp::cast_time_msec(spdlog::level::info,"pathcalc::g2pil::mp_add_PLDL(arrdet,BL_max) finished",start,end);
}

// ----- PL-only functions (no side effects on DetectorElement) -----

/// @brief Calculate path length (PL) only for a single DetectorElement (no side effects).
double pathcalc::g2pil::calc_PL(
  const DetectorElement& ele, const Grid2dPillar &g2pil, const double eps )
{
  double sum_PL = 0.0; // meters

  const Ray2d ray2d = ele.get_ray2d();
  const Ray3d ray3d = ele.get_ray3d();

  g2pil.traverse_ray_2d(ray2d, eps,
    [&](const Grid2d::Ixiy& ixiy) {
      const double delta_path = g2pil.get_delta_path(ixiy, ray3d);
      if (!std::isfinite(delta_path)) {
        THROW_ERROR("g2pil::calc_PL: delta_path is non-finite. delta_path={}", delta_path);
      }
      sum_PL += delta_path;
    });

  return sum_PL;
}

/// @brief Calculate path length (PL) only with beam length warning (no side effects).
double pathcalc::g2pil::calc_PL(
  const DetectorElement& ele, const Grid2dPillar &g2pil, const double BL_max, const double eps )
{
  double sum_PL = 0.0; // meters

  const Ray2d ray2d = ele.get_ray2d();
  const Ray3d ray3d = ele.get_ray3d();

  // BL_max check
  if (BL_max > 0.0) {
    auto [tf_hit, tmin_bl, tmax_bl] = ray2d.is_intersect(g2pil.get_cached_aabb2d());
    if (tf_hit && tmax_bl > BL_max) {
      LOG_WARN_ND("Ray beam length ({:.1f}m) exceeds BL_max ({:.1f}m). Ray2d={}",
                  tmax_bl, BL_max, ray2d.to_string());
    }
  }

  g2pil.traverse_ray_2d(ray2d, eps,
    [&](const Grid2d::Ixiy& ixiy) {
      const double delta_path = g2pil.get_delta_path(ixiy, ray3d);
      if (!std::isfinite(delta_path)) {
        THROW_ERROR("g2pil::calc_PL: delta_path is non-finite. delta_path={}", delta_path);
      }
      sum_PL += delta_path;
    });

  return sum_PL;
}

/// @brief Calculate PL vector for all elements in a DetectorPanel (OpenMP).
Eigen::VectorXf pathcalc::g2pil::mp_calc_PL(
  const DetectorPanel& panel, const Grid2dPillar &g2pil,
  const double BL_max, const double eps )
{
  const int nbinx = panel.get_nbinx();
  const int nbiny = panel.get_nbiny();
  const int n_ele = nbinx * nbiny;
  Eigen::VectorXf result = Eigen::VectorXf::Zero(n_ele);

  #pragma omp parallel for collapse(2) schedule(static)
  for (int iy = 0; iy < nbiny; iy++) {
    for (int ix = 0; ix < nbinx; ix++) {
      const DetectorElement& ele = panel.getDetectorElement(ix, iy);
      const double pl = calc_PL(ele, g2pil, BL_max, eps);
      result(iy * nbinx + ix) = static_cast<float>(pl);
    }
  }

  LOG_DEBUG("g2pil::mp_calc_PL: detid={}, n_ele={}", panel.get_detid(), n_ele);
  return result;
}

/// @brief Calculate PL for all elements in a DetectorPanelArray (flat VectorXf).
Eigen::VectorXf pathcalc::g2pil::mp_calc_PL(
  const DetectorPanelArray& arrdet, const Grid2dPillar &g2pil,
  const double BL_max, const double eps )
{
  std::chrono::system_clock::time_point start, end;
  start = time_now;

  const int n_det = arrdet.get_n_det();
  int total_ele = 0;
  for (int detid = 0; detid < n_det; detid++) {
    total_ele += arrdet.getDetectorPanel(detid).get_n_element();
  }

  Eigen::VectorXf result = Eigen::VectorXf::Zero(total_ele);
  int offset = 0;
  for (int detid = 0; detid < n_det; detid++) {
    const DetectorPanel& panel = arrdet.getDetectorPanel(detid);
    const int n_ele = panel.get_n_element();
    result.segment(offset, n_ele) = mp_calc_PL(panel, g2pil, BL_max, eps);
    offset += n_ele;
  }

  end = time_now;
  myapp::cast_time_msec(spdlog::level::info, "pathcalc::g2pil::mp_calc_PL(arrdet) finished", start, end);

  return result;
}

/// @brief add_PLDL and change path length and density length of DetectorPanel*
/// @details openmp version. \n it calls g2pil::add_PLDL_with_vec_tf_in_PL(det,g2pil)
/// @param panel DetectorPanel*
/// @param g2pil Grid2dPillar
/// @param ixiy_dist_thres the distance threshold for adjacent cubes
/// @param eps small value to avoid missing boxes on the edge
void pathcalc::g2pil::mp_add_PLDL_with_vec_tf_in_PL(
  DetectorPanel& panel, const Grid2dPillar &g2pil
, const int ixiy_dist_thres, const double eps )
{
  // Get detector ID
  const Detid detid = panel.get_detid();
  LOG_DEBUG("detid={}",detid);


  // Get number of bins (nbinx, nbiny) in 2D detector grid
  const int nbinx = panel.get_nbinx();
  const int nbiny = panel.get_nbiny();

  // Parallelize for loop
  #pragma omp parallel for collapse(2) schedule(static)

  // Process each detector element in 2D detector grid
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      // Get detector element information
      DetectorElement& ele = panel.callDetectorElement(ix,iy);
      Nhit n_hit_cub = add_PLDL_with_vec_tf_in_PL(ele,g2pil,ixiy_dist_thres, eps);

      // Display information
      if( iy%20==0 && ix%20==0 ) fprintf(stderr
        , "detid=%d, n_cub_hit=%4d: %s: ax=%7.4lf ay=%7.4lf, PL=%8.1lf, DL=%3.1E\r"
        , detid, n_hit_cub, ele.get_angle_unit().c_str()
        , ele.get_tx(),ele.get_ty(),ele.get_PL(), ele.get_DL());
    }
  }
}

/// @brief add_PLDL and change path length and density length of DetectorPanelArray&
/// @details openmp version. \n it calls mp_add_PLDL_with_vec_tf_in_PL(det,g2pil);
void pathcalc::g2pil::mp_add_PLDL_with_vec_tf_in_PL(
  DetectorPanelArray &arrdet , const Grid2dPillar &g2pil
, const int ixiy_dist_thres, const double eps )
{
  std::chrono::system_clock::time_point start, end;
  start = time_now;
  if( arrdet.get(FlgProg::tf_calc_PL_DL) ){
    LOG_WARN("arrdet.get(FlgProg::tf_calc_PL_DL)==true");
    SLEEP_MSEC(500);
  }
  for( Detid detid =0; detid < arrdet.get_n_det(); detid++ ){
    DetectorPanel& panel = arrdet.callDetectorPanel(detid);
    mp_add_PLDL_with_vec_tf_in_PL(panel,g2pil,ixiy_dist_thres, eps);
  }
  arrdet.set(FlgProg::tf_calc_PL_DL, true);
  end = time_now;
  myapp::cast_time_msec(spdlog::level::info,"pathcalc::g2pil::mp_add_PLDL_with_vec_tf_in_PL(arrdet) finished",start,end);
}



//=============================
// naive version
//=============================

// naive means calculate path length by incrementing the path length
// like, path += delta_path;

// calc_PL for a DetectorElement &ele,
// change the value of PL,
// and change the value of DL
// if out of x_axis or y_axis, PL = -9999; and end the calculation
// return number of the loop
// too naive, not fast
int pathcalc::naive::calc_PL(
  DetectorElement &ele, const Grid2dPillar &g2pil, const pathcalc::Parameters &prm )
{
  double path = prm.PL_min; // path length variable
  int nloop = 0; // number of loop
  Eigen::Vector3d v3_pos; // path position variable
  int ix,iy;
  const Eigen::Vector3d v3_pos0 = ele.get_v3_pos(); // initial position of the detector element
  const Eigen::Vector3d v3_dir0 = ele.get_v3_dir();// direction of the detector element
  double density;
  // path incremental loop start
  // path(=prm.PL_min) += prm.PL_pit, until path <= prm.PL_max
  while( path <= prm.PL_max ){
    v3_pos = v3_pos0 + path * v3_dir0;
    // check v3_pos in the range of x_axis.
    ix = g2pil.get_ix( v3_pos.x() );
    if( ix == Grid1d::OUT_OF_RANGE_LOWER || ix == Grid1d::OUT_OF_RANGE_UPPER ){
      nloop++;
      continue;
    }

    // check v3_pos in the range of y_axis.
    iy = g2pil.get_iy( v3_pos.y() );
    if( iy == Grid1d::OUT_OF_RANGE_LOWER || iy == Grid1d::OUT_OF_RANGE_UPPER ){
      nloop++;
      continue;
    }
// printf("ix=%d, iy=%d\n",ix,iy);
    const Pillar& cub = g2pil.getPillar(ix,iy);

    path += prm.PL_pit;
    nloop++;

    // if v3_pos of muon is outside the cube, continue
    if( false == cub.is_z_inside(v3_pos.z()) ) continue;

    // if v3_pos of muon is inside the cube,
    ele.add_PL(prm.PL_pit); // PL += prm.PL_pit

    // get density
    density = cub.get_density();

    // sum DL
    ele.add_DL(density*prm.PL_pit); // DL += density*prm.PL_pit
  } // path incremental loop end

  return nloop;
}

// calc Path-length for a DetectorPanel,
// which calls naive::calc_PL(ele,g2pil,prm)
// too naive, not fast
void pathcalc::naive::mp_calc_PL(
    DetectorPanel &det, const Grid2dPillar &g2pil
  , const pathcalc::Parameters &prm )
{
  // for debug
  det.get_x_axis().out_info(spdlog::level::debug);
  det.get_y_axis().out_info(spdlog::level::debug);

  int nloop=-1;
  int ix,iy;
  #pragma omp parallel for collapse(2) schedule(static)
  for(iy=0;iy<det.get_nbiny();iy++){
    for(ix=0;ix<det.get_nbinx();ix++){
      DetectorElement& ele = det.callDetectorElement(ix,iy);
      nloop = naive::calc_PL(ele,g2pil,prm); //call the function written above.
      LOG_DEBUG("nloop={}: {}: ax={:.4f} ay={:.4f}, PL={:.1f}, vx={:.4f}, vy={:.4f}, vz={:.4f}"
      , nloop, ele.get_angle_unit()
      , ele.get_tx(), ele.get_ty(), ele.get_PL()
      , ele.get_vx(), ele.get_vy(), ele.get_vz());
    }
  }
}


/// @brief add_PLDL and change path length and density length of DetectorPanel&
/// @details openmp version.
/// @ingroup pathCalculation
void pathcalc::naive::mp_add_DL( DetectorPanel& panel, const double DL)
{
  // Get detector ID
  const Detid detid = panel.get_detid();
  LOG_DEBUG("detid={}",detid);

  const int n_threads = omp_get_max_threads();

  // Get number of bins (nbinx, nbiny) in 2D detector grid
  const int nbinx = panel.get_nbinx();
  const int nbiny = panel.get_nbiny();

  // Parallelize for loop
  #pragma omp parallel for collapse(2) schedule(static)

  // Process each detector element in 2D detector grid
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      // Get detector element information
      DetectorElement& ele = panel.callDetectorElement(ix,iy);
      // add DL
      ele.add_DL(DL);
    }
  }
}


// true = all detector is above zmin of g3vox
// false at least one of the detector is below the zmin
bool pathcalc::g3vox::is_all_detector_above_zmin(
  const DetectorPanelArray &arrdet, const Grid3dVoxel &g3vox )
{
  const double zmin = g3vox.get_zmin();
  for(Detid detid=0;detid<arrdet.get_n_det();detid++){
    const DetectorPanel& panel = arrdet.getDetectorPanel(detid);
    const double z_det = panel.get_v3_pos().z();
    if( z_det < zmin ) return false;
  }
  return true;
}
