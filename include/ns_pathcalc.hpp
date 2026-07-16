/// @file ns_pathcalc.hpp
/// @brief Path length calculation namespace
/// @details Namespace for muon path length calculations through voxelized terrain and matrix construction.
#pragma once

// for openmp
#include <omp.h> 

#include <map>
#include <fstream>
#include <iostream>
#include <sstream> // istringstream
#include <string.h>
#include <string>
#include <cstdio>
#include <cmath>
#include <functional>  //for sorting
#include <algorithm>//for sorting
#include <vector>

#include <Eigen/Dense>

#include "cls_DetectorElement.hpp"
#include "cls_DetectorPanel.hpp"
#include "cls_DetectorPanelArray.hpp"
#include "cls_Grid2dPillar.hpp"
#include "cls_Voxel.hpp"
#include "cls_Grid3dVoxel.hpp"
#include "cls_Grid2dXYZ.hpp"
#include "cls_Grid1dXZ.hpp"
#include "cls_FluxTable.hpp"
#include "cls_PathLengthParameters.hpp"
#include "ns_geom_util.hpp"
// #include "cls_MatrixBuildParameters.hpp"

//#####################################################################################
//#####################################################################################
/// @namespace pathcalc
/// @brief Namespace containing functions for calculating path length, path length matrix (including number version), density length, etc.
/// @details This namespace provides comprehensive tools for muon path length calculations through voxelized terrain,
/// including both dense and sparse matrix representations. Most functions support OpenMP parallelization.
/// @ingroup matrixClasses
/// @ingroup pathCalculation
//#####################################################################################
//#####################################################################################
namespace pathcalc {
  class Parameters; // forward declaration
  class MatrixBuildParameters; // forward declaration

  //===================================================
  /// @namespace pathcalc::g2pil
  /// @brief Functions for Grid2dPillar
  /// @details Provides path length and density length calculation functions
  /// for 2D cuboid grids. Includes OpenMP parallelized versions for panel and array operations.
  //===================================================
  namespace g2pil {
    /// @brief Calculate and add path length (PL) and density length (DL) for a detector element
    /// @param ele Reference to DetectorElement
    /// @param g2pil Reference to Grid2dPillar
    /// @param eps Small value to avoid missing boxes on the edge (default: 1.0e-6)
    /// @return The number of hit cubes
    /// @ingroup pathCalculation
    int add_PLDL(
      DetectorElement& ele, const Grid2dPillar &g2pil, const double eps=1.0e-6 );

    /// @brief Calculate and add PL/DL with beam length warning
    /// @param ele Reference to DetectorElement
    /// @param g2pil Reference to Grid2dPillar
    /// @param BL_max Maximum beam length [meters]. If BL_max > 0 and beam length exceeds this,
    ///               a warning is logged but processing continues normally. If BL_max <= 0, no check.
    /// @param eps Small value to avoid missing boxes on the edge
    /// @return The number of hit cubes
    /// @ingroup pathCalculation
    int add_PLDL(
      DetectorElement& ele, const Grid2dPillar &g2pil, const double BL_max, const double eps );

    /// @brief Calculate and add PL/DL and return vector of hit tuples
    /// @param ele Reference to DetectorElement
    /// @param g2pil Reference to Grid2dPillar
    /// @param eps Small value to avoid missing boxes on the edge (default: 1.0e-6)
    /// @return Vector of hit (ix of g2pil, iy of g2pil, path_length in a cuboid)
    /// @ingroup pathCalculation
    std::vector<std::tuple<int,int,double>> add_PLDL_with_vec_tp(
      DetectorElement& ele, const Grid2dPillar &g2pil, const double eps=1.0e-6 );

    /// @brief Calculate and add PL/DL for a detector panel
    /// @details OpenMP version. Calls g2pil::calc_PLDL(panel,g2pil)
    /// @param panel Reference to DetectorPanel
    /// @param g2pil Reference to Grid2dPillar
    /// @param eps Small value to avoid missing boxes on the edge (default: 1.0e-6)
    /// @note Uses OpenMP
    /// @ingroup pathCalculation
    void mp_add_PLDL(
      DetectorPanel& panel, const Grid2dPillar &g2pil, const double eps=1.0e-6 );

    /// @brief Calculate and add PL/DL for a detector panel with beam length warning
    /// @details OpenMP version with BL_max warning check
    /// @param panel Reference to DetectorPanel
    /// @param g2pil Reference to Grid2dPillar
    /// @param BL_max Maximum beam length [meters]. If BL_max > 0, rays exceeding this log a warning.
    /// @param eps Small value to avoid missing boxes on the edge
    /// @note Uses OpenMP
    /// @ingroup pathCalculation
    void mp_add_PLDL(
      DetectorPanel& panel, const Grid2dPillar &g2pil, const double BL_max, const double eps );

    /// @brief Calculate and add PL/DL for a detector panel array
    /// @details OpenMP version. Calls mp_add_PLDL(det,g2pil)
    /// @param arrdet Reference to DetectorPanelArray
    /// @param g2pil Reference to Grid2dPillar
    /// @param eps Small value to avoid missing boxes on the edge (default: 1.0e-6)
    /// @note Uses OpenMP
    /// @ingroup pathCalculation
    void mp_add_PLDL(
      DetectorPanelArray &arrdet , const Grid2dPillar &g2pil, const double eps=1.0e-6);

    /// @brief Calculate and add PL/DL for a detector panel array with beam length warning
    /// @details OpenMP version with BL_max warning check
    /// @param arrdet Reference to DetectorPanelArray
    /// @param g2pil Reference to Grid2dPillar
    /// @param BL_max Maximum beam length [meters]. If BL_max > 0, rays exceeding this log a warning.
    /// @param eps Small value to avoid missing boxes on the edge
    /// @note Uses OpenMP
    /// @ingroup pathCalculation
    void mp_add_PLDL(
      DetectorPanelArray &arrdet , const Grid2dPillar &g2pil, const double BL_max, const double eps);

    /// @brief Calculate and add PL/DL with adjacency tracking in vec_tf_in_PL
    /// @param ele Reference to DetectorElement
    /// @param g2pil Reference to Grid2dPillar
    /// @param ixiy_dist_thres Threshold for index distance used in adjacency determination. If the index difference between current and previous cuboid is less than or equal to this value, they are considered adjacent (default: 1)
    /// @param eps Small value to avoid missing boxes on the edge (default: 1.0e-6)
    /// @return The number of hit cubes
    /// @ingroup pathCalculation
    int add_PLDL_with_vec_tf_in_PL(
        DetectorElement& ele, const Grid2dPillar &g2pil
      , const int ixiy_dist_thres=1, const double eps=1.0e-6 );

    /// @brief Calculate and add PL/DL with adjacency tracking for a detector panel
    /// @details OpenMP version. Calls g2pil::add_PLDL_with_vec_tf_in_PL(det,g2pil)
    /// @param panel Reference to DetectorPanel
    /// @param g2pil Reference to Grid2dPillar
    /// @param ixiy_dist_thres Threshold for index distance used in adjacency determination. If the index difference between current and previous cuboid is less than or equal to this value, they are considered adjacent (default: 1)
    /// @param eps Small value to avoid missing boxes on the edge (default: 1.0e-6)
    /// @note Uses OpenMP
    /// @ingroup pathCalculation
    void mp_add_PLDL_with_vec_tf_in_PL(
        DetectorPanel& panel, const Grid2dPillar &g2pil
      , const int ixiy_dist_thres=1, const double eps=1.0e-6 );

    /// @brief Calculate and add PL/DL with adjacency tracking for a detector panel array
    /// @details OpenMP version. Calls mp_add_PLDL_with_vec_tf_in_PL(det,g2pil)
    /// @param arrdet Reference to DetectorPanelArray
    /// @param g2pil Reference to Grid2dPillar
    /// @param ixiy_dist_thres Threshold for index distance used in adjacency determination. If the index difference between current and previous cuboid is less than or equal to this value, they are considered adjacent (default: 1)
    /// @param eps Small value to avoid missing boxes on the edge (default: 1.0e-6)
    /// @note Uses OpenMP
    /// @ingroup pathCalculation
    void mp_add_PLDL_with_vec_tf_in_PL(
        DetectorPanelArray &arrdet , const Grid2dPillar &g2pil
      , const int ixiy_dist_thres=1, const double eps=1.0e-6 );

    // ----- PL-only functions (no side effects on DetectorElement) -----

    /// @brief Calculate path length (PL) only for a single DetectorElement (no side effects).
    /// @param[in] ele DetectorElement (read-only: position/direction used, PL/DL not modified).
    /// @param[in] g2pil Grid2dPillar to trace through.
    /// @param[in] eps Small value to avoid missing boxes on the edge (default: 1.0e-6).
    /// @return Total path length through the grid (meters).
    double calc_PL(
      const DetectorElement& ele, const Grid2dPillar& g2pil, const double eps = 1.0e-6);

    /// @brief Calculate path length (PL) only for a single DetectorElement with BL_max warning.
    /// @param[in] ele DetectorElement (read-only).
    /// @param[in] g2pil Grid2dPillar to trace through.
    /// @param[in] BL_max Maximum beam length [meters]. If > 0 and exceeded, a warning is logged.
    /// @param[in] eps Small value to avoid missing boxes on the edge.
    /// @return Total path length through the grid (meters).
    double calc_PL(
      const DetectorElement& ele, const Grid2dPillar& g2pil, const double BL_max, const double eps);

    /// @brief Calculate PL vector for all elements in a DetectorPanel.
    /// @param[in] panel DetectorPanel (read-only).
    /// @param[in] g2pil Grid2dPillar to trace through.
    /// @param[in] BL_max Maximum beam length [meters].
    /// @param[in] eps Small value to avoid missing boxes on the edge.
    /// @return VectorXf containing PL per element (size = nbinx * nbiny).
    /// @note Uses OpenMP for parallel computation.
    Eigen::VectorXf mp_calc_PL(
      const DetectorPanel& panel, const Grid2dPillar& g2pil,
      const double BL_max, const double eps = 1.0e-6);

    /// @brief Calculate PL for all elements in a DetectorPanelArray (flat).
    /// @param[in] arrdet DetectorPanelArray (read-only).
    /// @param[in] g2pil Grid2dPillar to trace through.
    /// @param[in] BL_max Maximum beam length [meters].
    /// @param[in] eps Small value to avoid missing boxes on the edge.
    /// @return Single flat VectorXf containing PL for all panels concatenated.
    ///         Panel segments can be accessed via segment(offset, n_ele).
    /// @note Uses OpenMP for parallel computation.
    Eigen::VectorXf mp_calc_PL(
      const DetectorPanelArray& arrdet, const Grid2dPillar& g2pil,
      const double BL_max, const double eps = 1.0e-6);
  };

  //===================================================
  /// @namespace pathcalc::g3vox
  /// @brief Functions for Grid3dVoxel
  /// @details Provides path length and density length calculation functions
  /// for 3D voxel grids. Includes matrix generation capabilities and OpenMP parallelized versions.
  //===================================================
  namespace g3vox {
    /// @brief Determine whether all detectors are above zmin
    /// @param arrdet Instance of DetectorPanelArray
    /// @param g3vox Instance of Grid3dVoxel
    /// @return Returns true if all detectors are above zmin, false otherwise
    bool is_all_detector_above_zmin(
      const DetectorPanelArray &arrdet, const Grid3dVoxel &g3vox );

    /// @brief Calculate path length and density length for DetectorElement
    /// @param ele Reference to DetectorElement
    /// @param g3vox Reference to Grid3dVoxel
    /// @param eps Small value to avoid missing boxes on the edge (default: 1.0e-6)
    void add_PLDL( DetectorElement& ele, const Grid3dVoxel &g3vox, const double eps=1.0e-6 );

    /// @brief Calculate path length and density length for DetectorPanel
    /// @param panel Reference to DetectorPanel
    /// @param g3vox Reference to Grid3dVoxel
    /// @param eps Small value to avoid missing boxes on the edge (default: 1.0e-6)
    /// @note Uses OpenMP
    void mp_add_PLDL( DetectorPanel& panel, const Grid3dVoxel &g3vox, const double eps=1.0e-6 );

    /// @brief Calculate path length and density length for DetectorPanelArray
    /// @param arrdet Reference to DetectorPanelArray
    /// @param g3vox Reference to Grid3dVoxel
    /// @param eps Small value to avoid missing boxes on the edge (default: 1.0e-6)
    /// @note Uses OpenMP
    void mp_add_PLDL( DetectorPanelArray &arrdet, const Grid3dVoxel &g3vox, const double eps=1.0e-6 );

    /// @brief Calculate path/density length for DetectorElement with beam length warning.
    /// @param ele DetectorElement to update (PL/DL are added).
    /// @param g3vox Grid3dVoxel containing voxel data.
    /// @param BL_max Maximum beam length threshold [meters]. If <= 0, warning is disabled.
    /// @param eps Small scalar to avoid missing boundary voxels.
    void add_PLDL(
      DetectorElement &ele, const Grid3dVoxel &g3vox, const double BL_max, const double eps );

    /// @brief Calculate path/density length for DetectorPanel with beam length warning.
    /// @param panel DetectorPanel to update.
    /// @param g3vox Grid3dVoxel containing voxel data.
    /// @param BL_max Maximum beam length threshold [meters].
    /// @param eps Small scalar to avoid missing boundary voxels.
    /// @note Uses OpenMP for parallel element processing.
    void mp_add_PLDL(
      DetectorPanel& panel, const Grid3dVoxel &g3vox, const double BL_max, const double eps );

    /// @brief Calculate path/density length for DetectorPanelArray with beam length warning.
    /// @param arrdet DetectorPanelArray to update.
    /// @param g3vox Grid3dVoxel containing voxel data.
    /// @param BL_max Maximum beam length threshold [meters].
    /// @param eps Small scalar to avoid missing boundary voxels.
    void mp_add_PLDL(
      DetectorPanelArray &arrdet, const Grid3dVoxel &g3vox, const double BL_max, const double eps );

    /// @brief Calculate path length matrix for each DetectorElement
    /// @param ele Reference to DetectorElement
    /// @param g3vox Reference to Grid3dVoxel
    /// @param mat_PL Reference to Eigen::MatrixXf storing the calculated path length matrix
    /// @param prm pathcalc::Parameters used in the calculation
    /// @param eps Small value to avoid missing boxes on the edge (default: 1.0e-6)
    /// @note Called by mp_make_mat_PL function
    void calc_mat_PL(
      DetectorElement& ele, Grid3dVoxel &g3vox, Eigen::MatrixXf &mat_PL
      , const pathcalc::Parameters &prm, const double eps=1.0e-6 );

    /// @brief Calculate path length matrix for each DetectorPanel and return it
    /// @param panel Reference to DetectorPanel
    /// @param g3vox Reference to Grid3dVoxel
    /// @param prm pathcalc::Parameters used in the calculation
    /// @param eps Small value to avoid missing boxes on the edge (default: 1.0e-6)
    /// @return Returns the calculated path length matrix stored in Eigen::MatrixXf
    /// @note Uses OpenMP. Internally calculates path length matrix and facilitates n_hit_det calculation in g3vox
    Eigen::MatrixXf mp_make_mat_PL(
        DetectorPanel& panel, Grid3dVoxel &g3vox
      , const pathcalc::Parameters &prm, const double eps=1.0e-6 );

    /// @brief Calculate path length matrix for each DetectorPanel and return as SpMatf
    /// @param panel Reference to DetectorPanel
    /// @param g3vox Reference to Grid3dVoxel
    /// @param prm pathcalc::Parameters used in the calculation
    /// @param eps Small value to avoid missing boxes on the edge (default: 1.0e-6)
    /// @return Returns the calculated path length matrix as SpMatf
    /// @note Uses OpenMP and sparse matrix. Internally calculates path length matrix and facilitates n_hit_det calculation in g3vox
    SpMatf mp_make_spmat_PL(
        DetectorPanel& panel, Grid3dVoxel &g3vox
      , const pathcalc::Parameters &prm, const double eps=1.0e-6 );

    /// @brief Calculate path length matrix for all DetectorPanels and store in vector
    /// @param arrdet Reference to DetectorPanelArray
    /// @param g3vox Reference to Grid3dVoxel
    /// @param prm pathcalc::Parameters used in the calculation
    /// @param eps Small value to avoid missing boxes on the edge (default: 1.0e-6)
    /// @return Returns path length matrices for each detector stored in vector
    /// @note Calls mp_make_mat_PL(&ele)
    /// @note Uses OpenMP
    std::vector<SpMatf> mp_make_mat_PL(
        DetectorPanelArray &arrdet, Grid3dVoxel &g3vox
        , const pathcalc::Parameters &prm, const double eps=1.0e-6 );

    /// @brief Calculate path length matrix for all DetectorPanels and store in vector
    /// @param arrdet Reference to DetectorPanelArray
    /// @param g3vox Reference to Grid3dVoxel
    /// @param prm pathcalc::Parameters used in the calculation
    /// @param eps Small value to avoid missing boxes on the edge (default: 1.0e-6)
    /// @return Returns path length matrices for each detector stored in vector
    /// @note Calls mp_make_spmat_PL(&ele). Does not use Eigen::MatrixXf.
    /// @note Uses OpenMP
    std::vector<SpMatf> mp_make_vec_spmat_PL(
        DetectorPanelArray &arrdet, Grid3dVoxel &g3vox
        , const pathcalc::Parameters &prm, const double eps=1.0e-6 );
  };

  /// @brief calc path length for cylinder
  namespace vcyl {
    /// @brief calc and add DL for VerticalEllipticCylinderCapped to DetectorElement&
    void add_DL( DetectorElement& ele, const VerticalEllipticCylinderCapped &vcyl );

    /// @brief calc and add DL for VerticalEllipticCylinderCapped to DetectorPanel&
    /// @note Uses OpenMP
    void mp_add_DL( DetectorPanel& panel, const VerticalEllipticCylinderCapped &vcyl );

    /// @brief calc and add DL for VerticalEllipticCylinderCapped to DetectorPanelArray&
    void mp_add_DL( DetectorPanelArray &arrdet, const VerticalEllipticCylinderCapped &vcyl );

    /// @brief calc and add PL/DL for VerticalEllipticCylinderCapped to DetectorElement&
    void add_PLDL( DetectorElement& ele, const VerticalEllipticCylinderCapped &vcyl );

    /// @brief calc and add PL/DL for VerticalEllipticCylinderCapped to DetectorPanel&
    /// @note Uses OpenMP
    void mp_add_PLDL( DetectorPanel& panel, const VerticalEllipticCylinderCapped &vcyl );

    /// @brief calc and add PL/DL for VerticalEllipticCylinderCapped to DetectorPanelArray&
    void mp_add_PLDL( DetectorPanelArray &arrdet, const VerticalEllipticCylinderCapped &vcyl );
  };

  // NOTE: pathcalc::matrix namespace has been moved to the independent calc_dNdD namespace.
  // See ns_calc_dNdD.hpp for declarations.


  ///=============================
  /// @namespace pathcalc::naive
  /// @brief Naive version
  /// @details Naive means calculate path length by incrementing the path length like:
  /// @code
  /// path = PL_min;
  /// while(path < PL_max) {
  ///   path += delta_path;
  ///   // ...
  /// }
  /// @endcode
  ///=============================
  namespace naive {
    /// @brief Calculate path length for DetectorElement
    /// @param ele Reference to DetectorElement
    /// @param g2pil Reference to Grid2dPillar
    /// @param prm Instance of pathcalc::Parameters
    /// @return Number of calculation loop iterations
    /// @note Modifies PL value and terminates calculation with PL = -9999 if it goes outside x or y axis. Not fast.
    int calc_PL(
      DetectorElement& ele, const Grid2dPillar &g2pil, const pathcalc::Parameters &prm);

    /// @brief Calculate path length and DL, and modify PL value
    /// @param panel Reference to DetectorPanel
    /// @param g2pil Reference to Grid2dPillar
    /// @param prm Instance of pathcalc::Parameters
    /// @note Uses OpenMP
    void mp_calc_PL(
          DetectorPanel &ele, const Grid2dPillar &g2pil
        , const pathcalc::Parameters &prm );

    /// @brief Calculate and add DL to DetectorPanel
    /// @details OpenMP version. Calls all_DL(ele,DL)
    /// @param panel Reference to DetectorPanel
    /// @param DL Density length value
    /// @note Uses OpenMP
    /// @ingroup pathCalculation
    void mp_add_DL( DetectorPanel& panel, const double DL);
  };

};

