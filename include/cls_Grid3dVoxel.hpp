/// @file cls_Grid3dVoxel.hpp
/// @brief 3D voxel grid with density and path length tracking.
///
/// @details Defines Grid3dVoxel, a voxelized Grid3d container that stores a
///          Voxel payload per cell and maintains mappings between unique
///          indices and (ix, iy, iz) coordinates.
///
///          Typical workflow:
///          1. Build from Grid1d axes or convert from Grid2dPillar.
///          2. Allocate voxel storage and initialize densities/existence flags.
///          3. Build unique index maps with build_uqiv_umps().
///          4. Optionally add density structures (checkerboard/ellipsoid/cylinder).
///          5. Export cross sections or trace ray paths for path length queries.
///          6. Save/load binary snapshots as needed.
///
///          Coordinate system and units:
///          - Right-handed, z-axis is vertical (positive upward).
///          - Axis units follow Grid1d values (typically meters).
///          - Density values are in kg/m^3 unless noted otherwise.
///
///          Memory layout:
///          - Voxels are stored as vec_vec_vec_Voxel[iz][iy][ix] (z-major).
///
///          Thread safety:
///          - const methods are thread-safe for read-only access.
///          - mutating methods are not thread-safe.
///          - OpenMP is used in some methods; concurrent mutation of the same
///            instance is not safe.
///
///          I/O formats:
///          - ASCII output via out_* functions (column data).
///          - Binary save/load for full snapshots.
#pragma once

#include <array>
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
#include <atomic>

#include <Eigen/Dense>
#include "cls_Ray.hpp"
#include "cls_AABB.hpp"
#include "ns_tuple_int.hpp"
#include "cls_Voxel.hpp"
#include "cls_Grid3d.hpp"
#include "cls_Grid2dPillar.hpp"
#include "cls_Grid2dVoxel.hpp"
#include "cls_VoxelUniqueIndexMapContainer.hpp"
namespace fs = std::filesystem;

//########################################################################
//########################################################################
/// @class Grid3dVoxel
/// @brief Voxelized 3D grid storing per-cell Voxel payloads and unique indices.
/// @details Grid3dVoxel derives from Grid3d and stores Voxel data for each
///          grid cell along with mappings between unique indices and
///          (ix, iy, iz). It supports conversions from 2D cuboid DEM data,
///          cross-section exports, and ray-path length queries.
///
///          Example:
///          @code
///          Grid1d gx("x", 10, 0.0, 10.0, 1.0);
///          Grid1d gy("y", 10, 0.0, 10.0, 1.0);
///          Grid1d gz("z", 5, 0.0, 5.0, 1.0);
///          Grid3dVoxel g3vox(gx, gy, gz, 1);
///          g3vox.build_uqiv_umps(true, 0);
///          const Voxel& vox = g3vox.getVoxel(0, 0, 0);
///          (void)vox;
///          @endcode
/// @ingroup terrainClasses
//########################################################################
//########################################################################
class Grid3dVoxel : public Grid3d {
  public:
    //======================================================================
    /// @name forward declaration of parameter class for build instance
    ///@{
    class Parameters;
    class MergeParameters;
    class ReconstVoxelsParameters;
    class CheckerBoard3dParameters;
    class EllipsoidParameters;
    class CylinderParameters;
    class CuboidParameters;

    ///@} ------------------------------------------------------------------
  private:
    /// @brief name of this instance
    std::string name = "g3vox";

    /// @brief n_detector
    int n_detector = 0;

    /// @brief num of all Voxel with tf_exist==true
    /// @note Cached value for performance. Could be computed on-demand if memory is a concern.
    int n_vox_exist = 0;

    /// @brief scaler variable, vector of vector
    /// @details be careful of index order, .at(iz).at(iy).at(ix)
    std::vector<std::vector<std::vector<Voxel>>> vec_vec_vec_Voxel = {};

    /// @brief it stores the the pair of unique_index and ix,iy,iz of Voxel
    id_container::VoxID uqiv_container;

    /// @brief Flat array for O(1) uqiv lookup by (ix, iy, iz).
    /// @details Replaces unordered_map hash lookup in the hot path.
    ///          Index formula: ix + iy * nbinx + iz * nbinx * nbiny.
    ///          Non-existing voxels store UqivNotFound (-1).
    /// @note Built by build_flat_uqiv(), called at end of build_uqiv_umps().
    std::vector<Uqiv> flat_uqiv_;

    /// @brief Per-detector hit bitmask, indexed by (uqiv - uqiv_min).
    /// @details bit i = 1 means detector i has a hit in this voxel.
    ///          Size = uqiv_max - uqiv_min + 1 (allocated after build_uqiv_umps()).
    std::vector<uint64_t> vec_hit_det_;

    /// @brief Per-voxel element-hit counter, indexed by (uqiv - uqiv_min).
    std::vector<Nhit> vec_n_hit_ele_;

  public:
    //======================================================================
    /// @name constructor_destructor
    ///@{

    /// @brief default constructor using c++11 default, also initialize base class
    Grid3dVoxel() = default;

    /// @brief constructor from Grid3dVoxelParameters, Grid2dPillar DEM data, and num of detector
    Grid3dVoxel( const Grid2dPillar &g2pil
    , const Grid3dVoxel::Parameters &grdprm
    , const int n_detector );

    /// @brief copy constructor using c++11 default
    Grid3dVoxel( const Grid3dVoxel &org ) = default;

    /// @brief move constructor
    Grid3dVoxel( Grid3dVoxel &&org ) noexcept = default;

    /// @brief destructor
    virtual ~Grid3dVoxel() = default;  // Virtual destructor for base class

    /// @brief build from Grid1d x, y, z
    Grid3dVoxel(
      const Grid1d &x_axis_in, const Grid1d &y_axis_in , const Grid1d &z_axis_in
    , const int n_detector);
    
    /// @brief build from binary file
    Grid3dVoxel(const fs::path &path_in){ load(path_in); };
    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name operator
    ///@{

    /// @brief Inequality operator
    /// @note Name matching is not checked.
    bool operator!=(const Grid3dVoxel& other) const;

    /// @brief Equality operator (defined using inequality operator)
    /// @note Name matching is not checked.
    bool operator==(const Grid3dVoxel& other) const {
      return !(*this != other);
    }

    /// @brief assignment operator
    /// @details copy assignment operator using copy & swap idiom
    Grid3dVoxel& operator=(const Grid3dVoxel& other) = default;
    ///@} ------------------------------------------------------------------


    //======================================================================
    /// @name build functions
    ///@{

    /// @brief allocate memory for std::vector<std::vector<Voxel>> vec_vec_Voxel;
    /// @param n_det is number of detector
    void vec_vec_vec_memory_allocate(const int n_det);
    
    /// @brief allocate memory for std::vector<std::vector<Voxel>> vec_vec_Voxel;
    /// @details openmp version
    /// @param n_det is number of detector
    void mp_vec_vec_vec_memory_allocate(const int n_det);

    /// @brief build vec_vec_vec_Voxel elements
    /// @details openmp
    void build_element( const Grid3dVoxel::Parameters &prm, const int n_det  );

    /// @brief convert from Grid2dPillar
    /// @param g2pil is Grid2dPillar
    /// @param grdprm is Grid3dVoxel::Parameters
    /// @param n_detector is number of detector
    void convert_from_Grid2dPillar( const Grid2dPillar &g2pil
    , const Grid3dVoxel::Parameters &grdprm, const int n_detector );

    /// @brief build sqrt(mat_cov_dens_dash) from Eigen::VectorXf
    /// @param vecxf_diag_mat_cov_dens_dash is Eigen::VectorXf which has diagonal elements of mat_cov_dens_dash
    // ! 2025-06-18 14:11:35 Not used
    // Grid3dVoxel build_dens_err_g3vox(
    //   const Eigen::VectorXf &vecxf_diag_mat_cov_dens_dash ) const;

    ///@} ------------------------------------------------------------------



    //======================================================================
    /// @name getter_functions
    ///@{

    /// @brief get name
    std::string get_name() const { return name; };

    /// @brief get the number of detector
    int get_n_det() const { return n_detector; };

    /// @brief get uqiv_container.uqiv_not_assigned
    int get_uqiv_not_assigned() const {
      return id_container::VoxID::uqiv_not_assigned;
    };

    /// @brief get the pointer of uqiv_container.ump_uqiv_ixiyiz
    // const std::unordered_map<int, Ixiyiz>*
    //   get_ump_uqiv_ixiyiz_ref_pointer() const {
    //     return uqiv_container.get_ump_uqiv_ixiyiz_ref_pointer();
    //   };
    
    /// @brief get the pointer of uqiv_container.ump_ixiyiz_uqiv
    // const tuple_int::UmpInt3Int*
    //   get_ump_ixiyiz_uqiv_ref_pointer() const {
    //     return uqiv_container.get_ump_ixiyiz_uqiv_ref_pointer();
    //   };
    
    /// @brief get_function can't change member variable from unique_index
    const Voxel& getVoxel( const Uqiv unique_index_in ) const;

    /// @brief get_function can't change member variable from ix,iy,iz
    const Voxel& getVoxel( const int ix, const int iy, const int iz ) const;

    /// @brief get_function can't change member variable from (ix,iy,iz)
    const Voxel& getVoxel( const Ixiyiz &ixiyiz ) const{
      const auto [ix,iy,iz] = ixiyiz;
      return getVoxel(ix,iy,iz);
    };

    /// @brief get_function can't change member variable from x,y,z
    const Voxel& getVoxel( const double x, const double y, const double z ) const;

    /// @brief return n_vox_exist
    int get_n_vox_exist() const { return n_vox_exist; };

    /// @brief return VctorXf of density of all exist voxel
    /// @note This only applies to Voxels that have been assigned Uqiv.
    Eigen::VectorXf get_vecxf_density() const;

    /// @brief Get AABB3d(ix,iy,iz)
    /// @param ix, iy, iz
    AABB3d get_AABB3d(const int ix, const int iy, const int iz) const;
    
    /// @brief Get Eigen::Vector3d v3_AABB_min and v3_AABB_max
    /// @param tpl_vox(ix,iy,iz)
    AABB3d get_AABB3d( const Ixiyiz &ixiyiz ) const{
        const auto [ix,iy,iz] = ixiyiz;
        return get_AABB3d(ix,iy,iz);
      };
    
    /// @brief get voxel density from voxel ix,iy,iz
    double get_density( const int ix, const int iy, const int iz ) const {
      return getVoxel(ix,iy,iz).get_density();
    };

    /// @brief get voxel density from tpl_vox(ix,iy,iz)
    double get_density( const Ixiyiz &ixiyiz ) const {
      return getVoxel(ixiyiz).get_density();
    };

    /// @brief Get (density, AABB3d)
    /// @param ix, iy, iz
    std::tuple< double, AABB3d >
      get_density_AABB3d(const int ix, const int iy, const int iz) const;

    /// @brief Get (density, AABB3d)
    /// @param tpl_vox(ix,iy,iz)
    std::tuple< double, AABB3d >
      get_density_AABB3d( const Ixiyiz &ixiyiz ) const {
        const auto [ix,iy,iz] = ixiyiz;
        return get_density_AABB3d(ix,iy,iz);
      };
    
    /// @brief get the set of uqiv(unique index of voxel) all in uqiv_container.
    std::set<int> get_set_uqiv() const {
      return uqiv_container.get_set_uqiv();
    };

    /// @brief get the vector of uqiv(unique index of voxel) all in uqiv_container.
    std::vector<int> get_vec_uqiv() const {
      return uqiv_container.get_vec_uqiv();
    };

    /// @brief Get voxels at a given z height as Grid2dVoxel
    Grid2dVoxel get_Grid2dVoxel_z(const int iz) const;

    /// @brief Get voxels at a given z height as Grid2dVoxel
    Grid2dVoxel get_Grid2dVoxel_z(const double z_cross) const{
      const int iz = get_iz(z_cross);
      return get_Grid2dVoxel_z(iz);
    };

    /// @brief get the difference of density between two Grid3dVoxel
    /// @note Uses OpenMP for parallel iteration over voxels.
    Grid3dVoxel get_delta_dens(const Grid3dVoxel &g3vox_in) const;

    /// @brief get the immutable reference of base class Grid3d
    const Grid3d& getBaseRef() const { return *this; };

    /// @brief get the immutable ref of uqiv_container
    const id_container::VoxID&
      get_uqiv_container() const { return uqiv_container; };

    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name call_functions
    /// @brief call changable pointer of Voxel from ix, iy, iz
    ///@{
    
    /// @brief call changable ref of Voxel from ix, iy, iz
    Voxel& callVoxel( const int ix, const int iy, const int iz );

    /// @brief call changable reference of Voxel from (ix,iy,iz)
    Voxel& callVoxel( const Ixiyiz &ixiyiz ){
      const auto [ix,iy,iz] = ixiyiz;
      return callVoxel(ix,iy,iz);
    };

    /// @brief call changable reference of Voxel from x,y,z
    Voxel& callVoxel( const double x, const double y, const double z );

    /// @brief call changable reference of Voxel from unique_index
    Voxel& callVoxel( const int unique_index_in );

    /// @brief call mutable reference of base class Grid3d
    Grid3dVoxel& callBaseRef() { return *this; };

    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name check_functions
    ///@{

    /// @brief bool check inside or out of range
    /// @return if one of them is out of range, return true, otherwise return false.
    bool is_ixiyiz_out_of_range( const int ix, const int iy, const int iz ) const;
    /// @brief Function to check if vec_vec_vec_Voxel is the same, used in operator==, !=
    bool is_vec_vec_vec_Voxel_same( const Grid3dVoxel &g3vox_in ) const;

    /// @brief bool check the size of two Grid3dVoxel
    bool is_same_size( const Grid3dVoxel &g3vox_in ) const;

    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name read_write_functions
    ///@{

    /// @brief output Voxel(ix,iy,iz) info to FILE*
    void out_voxel( FILE *fout, const int ix, const int iy,  const int iz) const;
    /// @brief output Voxel(ix,iy,iz) info to FILE* with index info
    void out_voxel_with_index( FILE *fout, const int ix, const int iy,  const int iz) const;
    /// @brief output Voxel(uqiv) info to FILE*
    void out_voxel_with_index( FILE *fout, const int unique_index_in ) const;

    /// @brief output all Voxel(tf_exist==true) info to FILE*
    void out_all_exist( FILE *fout ) const;
    /// @brief output all Voxel(tf_exist==true) info to ascii file
    void out_all_exist( const fs::path& pathout ) const;

    /// @brief output all Voxel info to FILE* with index info
    void out_all_with_index( FILE *fout ) const;
    /// @brief output all Voxel info to ascii file with index info
    void out_all_with_index( const fs::path& pathout ) const;

    /// @brief output all Voxel(tf_exist==true) info to FILE* with index info
    void out_all_exist_with_index( FILE *fout ) const;
    /// @brief output all Voxel(tf_exist==true) info to ascii file with index info
    void out_all_exist_with_index( const fs::path& pathout ) const;

    /// @brief output all Voxel(tf_exist==true) info to FILE* with index info, n_hit_det > n_hit_det_thres
    void out_all_exist_with_index( FILE *fout, const Nhit n_hit_det_thres ) const;
    /// @brief output all Voxel(tf_exist==true) info to ascii file with index info, n_hit_det > n_hit_det_thres
    void out_all_exist_with_index( const fs::path& pathout, const Nhit n_hit_det_thres ) const;

    /// @brief output highest iz Voxel(tf_exist==true) info in each ix,iy to FILE* with index info
    void out_highest_exist_with_index( FILE *fout ) const;
    /// @brief output highest iz Voxel(tf_exist==true) info in each ix,iy to ascii file with index info
    void out_highest_exist_with_index( const fs::path& pathout ) const;

    /// @brief output lowest iz Voxel(tf_exist==true) info in each ix,iy to FILE* with index info
    void out_lowest_exist_with_index( FILE *fout, const Nhit n_hit_det_thres ) const;
    
    /// @brief output lowest iz Voxel(tf_exist==true) info in each ix,iy to ascii file with index info
    void out_lowest_exist_with_index( const fs::path& pathout, const Nhit n_hit_det_thres ) const;

    /// @brief output all Voxel(tf_exust==tf_exist_in) info to FILE*
    /// @details in the loop, out_voxel_with_index is called
    void out_all_with_ump_uqiv_ixiyiz( FILE *fout, const bool tf_exist_in ) const;
    /// @brief output all Voxel(tf_exust==tf_exist_in) info to ascii file
    /// @details in the loop, out_voxel_with_index is called
    void out_all_with_ump_uqiv_ixiyiz( const fs::path& pathout, const bool tf_exist_in ) const;

    /// @brief output all Voxel info to FILE*, the voxels are sorted by zxy
    void out_all_sort_zxy( FILE *fout ) const;
    /// @brief output all Voxel info to ascii file, the voxels are sorted by zxy
    void out_all_sort_zxy( const fs::path& pathout ) const;

    /// @brief output all Voxel(exist) info to FILE*, the voxels are sorted by zxy
    void out_all_sort_zxy_exist( FILE *fout ) const;
    /// @brief output all Voxel(exist) info to ascii file, the voxels are sorted by zxy
    void out_all_sort_zxy_exist( const fs::path& pathout ) const;

    /// @brief output all Voxel(exist) info to FILE*, the voxels are sorted by xyz
    void out_all_sort_xyz_exist( FILE *fout ) const;
    /// @brief output all Voxel(exist) info to ascii file, the voxels are sorted by xyz
    void out_all_sort_xyz_exist( const fs::path& pathout ) const;

    /// @brief Output x, y, density data at elevation z_in.
    /// @param pathout Output file path
    /// @param z_in z-coordinate
    /// @param tf_only_exist If true, output only voxels with tf_exist==true.
    /// @param dens_false Density value when tf_exist==false
    /// @param n_det_thres Threshold for n_hit_det
    void out_cross_section_z(
      const fs::path& pathout, const double z_in
      , const bool tf_only_exist=true, const double dens_false=-9999.0
      , const Nhit n_det_thres=0 ) const;

    /// @brief Output x, y, 0(tf_exist==false) or 1(tf_exist==true) data at elevation z_in.
    /// @param pathout Output file path
    /// @param z_in z-coordinate
    void out_cross_section_z_mask(
      const fs::path& pathout, const double z_in) const;

    //----------------------------------------------------------------------
    /// @brief Structure for passing parameters to out_cross_section_z
    struct CrossSectionZParameters {
      double xmin = 999.0; ///< Minimum x value for output
      double xmax =-999.0; ///< Maximum x value for output
      double xstep=-1.0; ///< Step size for x
      double ymin = 999.0; ///< Minimum y value for output
      double ymax =-999.0; ///< Maximum y value for output
      double ystep=-1.0; ///< Step size for y
      double zmin = 999.0; ///< Minimum z value for output
      double zmax =-999.0; ///< Maximum z value for output
      double zstep=-1.0; ///< Step size for z
      int n_detector = 0; ///< Number of detectors
      bool output_binary = false; ///< If true, output in binary format; otherwise ASCII
    };
    //----------------------------------------------------------------------


    /// @brief Output header information for out_cross_section_z.
    void out_cross_section_z_header( FILE *fout, const CrossSectionZParameters& prm_zcross ) const;

    // Output cross-sections from z_cross_min to z_cross_max with z_step increment.
    void out_cross_section_z_all(
      const fs::path& pathout
    , const CrossSectionZParameters &prm_zcross ) const;

    /// @brief Output cross-sections in binary format.
    /// @details Uses the same record layout as Grid2dPillar binary cross-section:
    ///   header (magic + grid info) followed by per-voxel binary records.
    /// @param pathout Output file path
    /// @param prm_zcross Cross-section parameters
    void out_cross_section_z_all_binary(
      const fs::path& pathout
    , const CrossSectionZParameters &prm_zcross ) const;

    /// @brief Output x, y, density data at elevation z_in.
    void out_cross_section_z_7clm(
      const fs::path& pathout, const double z_in) const;

    /// @brief out_cross_section_z_10clm \n
    /// Output x, y, density data at elevation z_in.
    void out_cross_section_z_10clm(
      const fs::path& pathout, const double z_in) const;

    /// @brief Output cross-section image using density vector
    /// @details Duplicate the current Grid3dVoxel object as a base,
    ///          set the specified density vector, and output the Z-direction cross-section.
    ///          Output filename is `<prefix>_zcross_all.tmp` (ASCII) or
    ///          `<prefix>_zcross_all.tmpbin` (when output_binary is true).
    /// @param[in] prefix Prefix for output file and voxel name(e.g., "g3vox_nagainv0"）
    /// @param[in] vec_density Density vector to set（Eigen::VectorXf）
    /// @param[in] prm_zcross Parameter structure for cross-section output
    /// @param[in] tf_zero_init If true, initialize to 0 before setting density vector
    /// @return Output Grid3dVoxel instance
    /// @note This function does not modify the caller object state(const function)
    ///       Internally duplicated Grid3dVoxel instance is used for output.
    Grid3dVoxel write_density_to_cross_section(
        const std::string& prefix
      , const Eigen::VectorXf& vec_density
      , const CrossSectionZParameters& prm_zcross
      , bool tf_zero_init = true
      , const std::array<double,4>& density_quad = {}) const;

    /// @brief Classify tf_exist=false voxels as upper/lower/lateral shell and assign density.
    /// @param[in] density_quad Shell density [prior, upper, lower, lateral].
    ///            Skipped if all shell values (indices 1-3) are zero.
    void apply_shell_density(const std::array<double,4>& density_quad);

    /// @brief   just display voxels in the nearest column from x_in, y_in
    /// @details for debug
    void disp_nearest_column( const double x_in, const double y_in ) const;

    /// @brief Output (x, y, mat_cov_dens(base_uqiv, others)) in ASCII for each z value
    /// @param prefix Prefix for output file
    /// @param uqiv_base Base uqiv to target
    /// @param mat_cov_dens Target matrix (assumes ColMajor)
    /// @param zmin Minimum z value for output
    /// @param zmax Maximum z value for output
    /// @param zstep Step size for z
    void out_xy_sqrt_matcovdens_by_z(
      const std::string &prefix,
      const Uqiv uqiv_base,
      const Eigen::MatrixXf &mat_cov_dens,
      const double zmin, const double zmax, const double zstep) const;

    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name calculation_functions
    ///@{

    /// @brief count num of voxel
    /// @return get_nbinx()*get_nbiny()*get_nbinz();
    int count_n_voxel_all() const { return get_nbinx()*get_nbiny()*get_nbinz(); };

    /// @brief count num of voxel with tf_exist==tf_exist_in
    int count_n_voxel( const bool tf_exist_in ) const;

    ///@} ------------------------------------------------------------------

    // int get_uqiv_min() const { return ump_uqiv_ixiyiz.begin()->first; };
    // int get_uqiv_max() const { return ump_uqiv_ixiyiz.end()->first; };

    //======================================================================
    /// @name setter_functions
    ///@{

    /// @brief copy from const Grid3dVoxel &prg
    void set( const Grid3dVoxel &org );

    /// @brief set name of this instance
    void set_name( const std::string name_in ){ name=name_in; };

    /// @brief set n_detector
    void set_n_detector( const int n_detector_in ){
      n_detector = n_detector_in;
    };

    /// @brief set n_vox_exist = n_vox_exist_in;
    void set_n_vox_exist( const int n_vox_exist_in ){
      n_vox_exist = n_vox_exist_in;
    };

    /// @brief set density = density_in, for all uqiv
    /// @note SI unit. kg/m3
    void set_uniform_density_uqiv( const double density_in);

    /// @brief set density = density_in, for all uqiv
    /// @details openmp version
    /// @note SI unit. kg/m3
    void mp_set_uniform_density_uqiv( const double density_in );
    
    /// @brief set density = density_in, for all exist voxels
    /// @note SI unit. kg/m3
    /// @note Uses OpenMP for parallel iteration over voxels.
    void set_uniform_density_ixiyiz( const double density_in);

    /// @brief set density = density_in, for all exist voxels
    /// @details openmp version
    /// @note SI unit. kg/m3
    void mp_set_uniform_density_ixiyiz( const double density_in );
    
    /// @brief set density of all exist voxel from Eigen::VectorXf vecxf_dens_in
    /// @note SI unit. kg/m3
    void set_density( const Eigen::VectorXf &vecxf_dens_in );

    /// @brief set n_hit_ele for voxel(ix,iy,iz)  tf_exist is true
    void set_n_hit_ele(const int ix,const int iy ,const int iz
      , const Nhit n_hit_ele_in );
    
    /// @brief set n_hit_ele for all voxel tf_exist is true
    void set_n_hit_ele( const Nhit n_hit_ele_in );

    /// @brief set uqiv_container
    void set_uqiv_container( const Grid3dVoxel &g3vox_in ) {
      uqiv_container = g3vox_in.get_uqiv_container();
    };

    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name merge_split_functions
    ///@{

    /// @brief merge voxels
    /// @return merged Grid3dVoxel
    /// @details in the remain instance, tf_exist become false in the merged voxels
    /// @note Not used 2024-12-03 10:52:03
    // Grid3dVoxel merge( 
    //     const int ix_offset, const int ix_factor
    //   , const int iy_offset, const int iy_factor
    //   , const int iz_offset, const int iz_factor );

    /// @brief merge voxels type2
    /// @param xcnt, ycnt, zcnt is the center position of merge
    /// @param merge_factor_x, merge_factor_y, merge_factor_z are merge factors in each axis
    /// @return merged Grid3dVoxel
    /// @details in the remain instance, tf_exist become false in the merged voxels \n
    /// This function merges (combines) 3D grid voxel data (Grid3dVoxel) \n
    /// to create a new Grid3dVoxel object. \n
    /// Merging is performed based on specified offset and multiplier. \n
    /// The number of bins along each axis of the merged voxels decreases, and the voxel size increases.
    /// @note Uses OpenMP for parallel iteration over merged voxel blocks.
    Grid3dVoxel merge(
        const double xcnt, const int merge_factor_x
      , const double ycnt, const int merge_factor_y
      , const double zcnt, const int merge_factor_z );

    /// @brief Merge voxels to lower resolution, then optionally apply region masks.
    /// @param[in] prm_merge Merge configuration (center, factors).
    /// @param[in] prm_reconst Region mask parameters (AABB and/or elliptic cylinder).
    /// @return New Grid3dVoxel with merged (and optionally masked) voxels.
    /// @details When prm_reconst.tf_aabb is true, voxels outside the AABB defined by
    ///          (x/y_aabb_cells * merged_interval) centered at merge_center are set
    ///          to tf_exist=false. When prm_reconst.tf_cylinder is true, voxels outside
    ///          the elliptic cylinder defined by (cylinder_radius_x/y_cells * merged_interval)
    ///          are also set to tf_exist=false. Unique index maps are rebuilt once after
    ///          all masks are applied.
    /// @note The AABB and cylinder extents are specified in post-merge voxel counts.
    Grid3dVoxel merge(
        const Grid3dVoxel::MergeParameters &prm_merge,
        const Grid3dVoxel::ReconstVoxelsParameters &prm_reconst );

    /// @brief Function to execute Grid3dVoxel::merge multiple times
    /// @note 2025-03-19 12:21:08 Not used
    // std::vector<Grid3dVoxel>
    //   get_merged_all( const Grid3dVoxel::Parameters &prm_g3vox );

    /// @brief Overlay merged-grid densities onto this (pre-merge) grid.
    /// @details For each tf_exist=true voxel in g3vox_merged, maps its density
    ///          to the corresponding pre-merge voxels using Grid1d::get_original_index_min_max().
    ///          Each merged voxel covers factor^3 pre-merge voxels, all set to the same density.
    /// @param[in] g3vox_merged Merged grid with reconstruction densities.
    /// @return Copy of this grid with merged densities overlaid.
    /// @note Voxels in this grid not covered by any merged voxel retain their original density.
    /// @note Uses OpenMP for parallel iteration over merged voxels.
    Grid3dVoxel overlay_merged_density(const Grid3dVoxel& g3vox_merged) const;

    /// @brief Overlay merged density with uniform prior background.
    /// @details Pass 1: fills reconst_volume (tf_exist=true in g3vox_merged_all)
    ///          with uniform_prior_density.
    ///          Pass 2: overlays reconstruction results from g3vox_rec.
    /// @param[in] g3vox_merged_all  Full reconst_volume grid (all merged voxels).
    /// @param[in] g3vox_rec         Reconstruction result (subset with tf_exist=true).
    /// @param[in] uniform_prior_density  Uniform prior density value (kg/m^3).
    /// @return Copy of this grid with prior background and reconstruction overlay.
    /// @note Uses OpenMP for parallel iteration over merged voxels.
    Grid3dVoxel overlay_merged_density(
        const Grid3dVoxel& g3vox_merged_all,
        const Grid3dVoxel& g3vox_rec,
        double uniform_prior_density) const;

    /// @brief Increase Grid3dVoxel resolution. n_hit_det and n_hit_ele are reset.
    Grid3dVoxel get_split_g3vox(
      const int x_factor, const int y_factor, const int z_factor ) const;

    ///@} ------------------------------------------------------------------


    //======================================================================
    /// @name make_void_g2pil_functions
    ///@{

    /// @brief return highest voxel in a bottom pixel ix,iy
    double get_highest_exist_z( const int ix, const int iy ) const;

    /// @brief return highest voxel in a bottom x, y
    double get_highest_exist_z( const double x, const double y ) const;

    /// @brief Return lowest z-coordinate where tf_exist==true in column (ix,iy).
    /// @param ix Grid index in x-direction
    /// @param iy Grid index in y-direction
    /// @return Lowest z-coordinate (zlow of the voxel) where voxel exists,
    ///         or zmax if no voxel exists in the column.
    double get_lowest_exist_z( const int ix, const int iy ) const;

    /// @brief Return lowest z-coordinate where tf_exist==true at position (x,y).
    /// @param x x-coordinate
    /// @param y y-coordinate
    /// @return Lowest z-coordinate where voxel exists, or zmax if none.
    double get_lowest_exist_z( const double x, const double y ) const;

    /// @brief Converted from g2pil to g3vox, and then \n
    /// create a cavity in g2pil from merged g3vox_input (this instance).
    /// @return Function that returns the outer shell part as g2pil.
    /// @details zmax_voxel_with_tf_exist_true in g3vox become zmin in g2pil
    /// @note 2023-12-27 17:15:47 Functionality moved to \"dem_operator::make_shell_from_g2pil_and_g3vox\".
    Grid2dPillar make_void_g2pil(const Grid2dPillar &g2pil ) const;
    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name add_structure_functions
    ///@{

    /// @brief add CheckerBoard3dParameters structure
    void add_density_structure( const CheckerBoard3dParameters &prm );

    /// @brief add Ellipsoid structure
    /// @note Uses OpenMP for parallel iteration over voxels.
    void add_density_structure( const Grid3dVoxel::EllipsoidParameters &prm );

    /// @brief add Cylinder structure
    /// @note Uses OpenMP for parallel iteration over voxels.
    void add_density_structure( const Grid3dVoxel::CylinderParameters &prm );

    /// @brief add Cuboid (rectangular box) structure
    /// @note Uses OpenMP for parallel iteration over voxels.
    void add_density_structure( const Grid3dVoxel::CuboidParameters &prm );

    /// @brief add density structure of checker board, ellipsoid, cylinder, and cuboid
    /// @param grdprm Parameter container with structure definitions
    void add_density_structure_all( const Grid3dVoxel::Parameters &grdprm );

    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name pathcalc_constants
    ///@{

    /// @brief constant path length for not intersecting any voxels
    static constexpr double PATH_NO_HIT = 0.0;
    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name path_length_calc_functions
    ///@{

    /// @brief get intersect path length between Ray3d and Voxel(ix,iy,iz)
    /// @details Mainly called by pathcalc. Find intersection between ray3d and Voxel(tpl=ix,iy,iz).
    /// @return If no intersection, return PATH_NO_HIT. If there is an intersection, return the computed delta_path.
    /// @ingroup ray_tracing_functions
    double get_delta_path(
      const Ixiyiz &ixiyiz, const Ray3d &ray3d) const;

    /// @brief If the voxel in g3vox is hit by the ray from g2det, \n
    /// increment n_hit_det of the hit voxel.
    /// @details Information about whether it was hit is held by mat_path_len.
    // void incr_n_hit_det(
    //   const Eigen::MatrixXf &mat_path_len );
    ///@} ------------------------------------------------------------------


    //======================================================================
    /// @name uqiv_container_Grid3dVoxel
    ///@{
    
    /// @brief get reference of uqiv_container
    id_container::VoxID& get_uqiv_container(){
      return uqiv_container;
    };

    /// @brief get reference of ump_uqiv_ixiyiz
    const UmpUqivIxiyiz& get_ump_uqiv_ixiyiz_ref() const {
      return uqiv_container.get_ump_uqiv_ixiyiz_ref();
    };

    /// @brief  get reference of ump_ixiyiz_uqiv
    const UmpIxiyizUqiv& get_ump_ixiyiz_uqiv_ref() const {
      return uqiv_container.get_ump_ixiyiz_uqiv_ref();
    };

    /// @brief get copy of ump_uqiv_ixiyiz
    UmpUqivIxiyiz get_ump_uqiv_ixiyiz_copy() const {
      return uqiv_container.get_ump_uqiv_ixiyiz_copy();
    };

    /// @brief  get reference of ump_ixiyiz_uqiv
    UmpIxiyizUqiv get_ump_ixiyiz_uqiv_copy() const {
      return uqiv_container.get_ump_ixiyiz_uqiv_copy();
    };

    /// @brief Get const reference to the uqiv-ixiyiz bimap.
    const auto& get_bimap_uqiv_ixiyiz_ref() const {
      return uqiv_container.get_bimap_uqiv_ixiyiz_ref();
    };

    /// @brief Get a copy of the uqiv-ixiyiz bimap.
    auto get_bimap_uqiv_ixiyiz_copy() const {
      return uqiv_container.get_bimap_uqiv_ixiyiz_copy();
    };

    /// @brief clear ump_uqiv_ixiyiz and ump_ixiyiz_uqiv
    /// @details it calls VoxID::clear_uqiv_umps()
    void clear_uqiv_umps(){
      uqiv_container.clear_uqiv_umps();
    };

    /// @brief set min unique_index
    void set_uqiv_min( const int min_unique_index_in ){
      uqiv_container.set_uqiv_min(min_unique_index_in);
    };
    
    /// @brief set max unique_index
    void set_uqiv_max( const int max_unique_index_in ){
      uqiv_container.set_uqiv_max(max_unique_index_in);
    };

    /// @brief allocate memory for ump_uqiv_ixiyiz and ump_ixiyiz_uqiv
    void reserve_uqiv_umps( const size_t size_in ){
      uqiv_container.reserve_uqiv_umps(size_in);
    };

    /// @brief insert key and values to ump_uqiv_ixiyiz and ump_ixiyiz_uqiv
    void insert_to_uqiv_umps( const Uqiv uqiv, const Ixiyiz& ixiyiz ){
      uqiv_container.insert_to_uqiv_umps(uqiv, ixiyiz);
    };

    /// @brief build ump_uqiv_ixiyiz and ump_ixiyiz_uqiv \n
    /// for only voxels which tf_exist==tf_exist_in
    /// @return n_voxel_exist
    int build_uqiv_umps(
      const bool tf_exist_in=true, const int uqiv_start = 0 );

    /// @brief Build flat_uqiv_ array for O(1) uqiv lookup.
    /// @details Populates flat_uqiv_ (size = nbinx * nbiny * nbinz) from the
    ///          existing unordered_map. Called automatically at the end of
    ///          build_uqiv_umps(). Non-existing voxels are set to UqivNotFound.
    /// @note Thread-safe: No (mutates internal state).
    void build_flat_uqiv();

    /// @brief Filter voxels by n_hit condition and reassign unique_index.
    /// @return Map of old unique_index and validity bool (used to remove matrix later)
    std::map<Grid3d::Uqiv,Grid3d::Uqiv> re_assign_uqiv_by_nhit_det(
        const bool tf_exist_in, const Grid3dVoxel::Parameters &prm_g3vox
      , const int uqiv_start );

    /// @brief get (ix,iy,iz) from unique_index
    Ixiyiz get_ixiyiz( const Uqiv uqiv_in ) const {
      return uqiv_container.get_ixiyiz(uqiv_in);
    };

    /// @brief get unique index of voxels from ix, iy ,iz
    Uqiv get_uqiv( const int ix_in, const int iy_in, const int iz_in ) const;

    /// @brief get unique index of voxels from (ix, iy ,iz)
    Uqiv get_uqiv( const Ixiyiz& ixiyiz ) const {
      return get_uqiv(ixiyiz[0],ixiyiz[1],ixiyiz[2]);
    };

    /// @brief O(1) uqiv lookup via flat array — no hash computation.
    /// @param ix_in X grid index.
    /// @param iy_in Y grid index.
    /// @param iz_in Z grid index.
    /// @return Unique voxel index, or UqivNotFound if not assigned.
    /// @note Requires build_flat_uqiv() to have been called (done automatically
    ///       at the end of build_uqiv_umps). Falls back to get_uqiv() if the
    ///       flat array has not been built.
    /// @note Thread-safe: Yes (read-only). Time complexity: O(1).
    Uqiv get_uqiv_fast( const int ix_in, const int iy_in, const int iz_in ) const;

    /// @brief O(1) uqiv lookup via flat array — overload for Ixiyiz.
    /// @param ixiyiz Voxel grid coordinates (ix, iy, iz).
    /// @return Unique voxel index, or UqivNotFound if not assigned.
    /// @note Thread-safe: Yes (read-only). Time complexity: O(1).
    Uqiv get_uqiv_fast( const Ixiyiz& ixiyiz ) const {
      return get_uqiv_fast(ixiyiz[0], ixiyiz[1], ixiyiz[2]);
    };

    /// @brief get min unique_index
    Uqiv get_uqiv_min() const {
      return uqiv_container.get_uqiv_min();
    };

    /// @brief get max unique_index
    Uqiv get_uqiv_max() const {
      return uqiv_container.get_uqiv_max();
    };

    /// @brief output ump_uqiv_ixiyiz
    void out_ump_uqiv_ixiyiz( const fs::path& pathout ) const{
      uqiv_container.out_ump_uqiv_ixiyiz(pathout);
    };

    /// @brief output ump_ixiyiz_uqiv
    void out_ump_ixiyiz_uqiv( const fs::path& pathout ) const{
      uqiv_container.out_ump_ixiyiz_uqiv(pathout);
    };

    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name hit_data_management (uqiv-indexed)
    ///@{

    /// @brief Allocate vec_hit_det_ and vec_n_hit_ele_ based on uqiv range.
    /// @details Called automatically at the end of build_uqiv_umps().
    void allocate_hit_data();

    /// @brief Zero-fill vec_hit_det_ and vec_n_hit_ele_.
    void clear_hit_data();

    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name hit_data_write (OpenMP-safe via atomic_ref)
    ///@{

    /// @brief Record a detector hit for the voxel identified by uqiv.
    /// @details Thread-safe: uses std::atomic_ref with relaxed ordering.
    /// @param[in] uqiv   Unique voxel index.
    /// @param[in] detid  Detector index (0-based, must be < 64).
    void record_hit_det(Uqiv uqiv, int detid);

    /// @brief Increment the element-hit counter for the voxel identified by uqiv.
    /// @details Thread-safe: uses std::atomic_ref with relaxed ordering.
    /// @param[in] uqiv  Unique voxel index.
    void record_hit_ele(Uqiv uqiv);

    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name hit_data_read
    ///@{

    /// @brief Return the raw hit bitmask for the voxel identified by uqiv.
    /// @param[in] uqiv  Unique voxel index.
    /// @return Hit bitmask (bit i = 1 means detector i has a hit).
    uint64_t get_hit_det(Uqiv uqiv) const;

    /// @brief Return true if detector det_id has a hit in the voxel identified by uqiv.
    /// @param[in] uqiv    Unique voxel index.
    /// @param[in] det_id  Detector index (0-based, must be < 64).
    /// @return true if detector det_id has a hit.
    bool get_tf_hit_grid(Uqiv uqiv, int det_id) const;

    /// @brief Return the number of detectors that hit the voxel identified by uqiv.
    /// @param[in] uqiv  Unique voxel index.
    /// @return Number of detectors with hits (popcount of bitmask).
    int get_n_hit_det_grid(Uqiv uqiv) const;

    /// @brief Return the element-hit count for the voxel identified by uqiv.
    /// @param[in] uqiv  Unique voxel index.
    /// @return Element-hit count.
    Nhit get_n_hit_ele_grid(Uqiv uqiv) const;

    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name binary_io_Grid3dVoxel
    ///@{

    /// @brief save vec_vec_vec_Voxel to std::ofstream, is used higher class
    void save_vec_vec_vec_Voxel( std::ofstream &ofs ) const;
    
    /// @brief load vec_vec_vec_Voxel to std::ifstream, is used higher class
    void load_vec_vec_vec_Voxel(std::ifstream& ifs);

    /// @brief save Grid3dVoxel to std::ofstream, is used \n 
    /// void save( const fs::path& pathout) const;
    void save( std::ofstream &ofs ) const;

    /// @brief save Grid3dVoxel to fs::path
    /// @details output file is binary
    void save( const fs::path& pathout) const;

    /// @brief load Grid3dVoxel from std::ifstream, \n
    /// is in void load( const fs::path &path_in );
    void load( std::ifstream &ifs );

    /// @brief load Grid3dVoxel from fs::path
    /// @details input file should be binary
    void load( const fs::path &path_in );
    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name other_functions
    ///@{
    
    /// @brief   make matrix_voxel_distance, Probably not used
    /// @note Uses OpenMP for parallel distance computation.
    Eigen::MatrixXf mp_make_voxel_distance_matrix( ) const;

    ///@} ------------------------------------------------------------------
    

};
