// cls_Grid3dVoxel.cpp
#include "cls_Grid3dVoxel.hpp"
#include "ns_myapp.hpp"
#include "ns_mymacro.hpp"
#include "ns_mylogger.hpp"
#include "cls_Grid3dVoxelParameters.hpp"
#include "spdlog_pch.hpp"
#include "ns_iodir.hpp"
#include "ns_geom_util.hpp"
#include <atomic>
#include <set>

namespace {
/// @brief Integer floor division (rounds toward negative infinity).
int floor_div(int a, int b) {
  int q = a / b;
  int r = a % b;
  // Adjust when remainder exists and signs differ.
  if (r != 0 && ((a < 0) != (b < 0))) q -= 1;
  return q;
}
} // anonymous namespace

//####################################################################################
//####################################################################################
// class Grid3dVoxel
// as derived class of the base class Grid3d
//####################################################################################
//####################################################################################

// build from Grid1d x3
Grid3dVoxel::Grid3dVoxel(
  const Grid1d &x_axis_in, const Grid1d &y_axis_in , const Grid1d &z_axis_in
, const int n_detector_in)
: Grid3d(x_axis_in,y_axis_in,z_axis_in)
{
  set_n_detector(n_detector_in);
  mp_vec_vec_vec_memory_allocate(n_detector_in);
}

// constructor from Grid3dVoxel::Parameters, Grid2dPillar DEM data
Grid3dVoxel::Grid3dVoxel(
  const Grid2dPillar &g2pil, const Grid3dVoxel::Parameters &grdprm
, const int n_detector_in )
: Grid3dVoxel()
{
  LOG_INFO("constructor from Parameters and Grid2dPillar DEM data");
  set_name("g3vox_from_"+g2pil.get_name());
  set_n_detector(n_detector_in);
  set_n_vox_exist(0);
  set_uqiv_min(-1);
  set_uqiv_max(-1);

  LOG_INFO("constructor: call convert_from_Grid2dPillar");
  convert_from_Grid2dPillar(g2pil,grdprm,n_detector_in);

  LOG_INFO("Grid3dVoxel constructor : add density structure");
  add_density_structure_all(grdprm);
}

// add density structure of checker board and ellipsoid
void Grid3dVoxel::add_density_structure_all( const Grid3dVoxel::Parameters &grdprm )
{
  // add checker board density structure
  for( auto param : grdprm.vec_checkerboard_prm ){
    if( param.tf_exec==false ) continue;
    LOG_INFO("adding checker board density structure: {}", param.name);

    // Compute actual cell dimensions from grid intervals and multipliers
    const double x_interval = get_x_interval();
    const double y_interval = get_y_interval();
    const double z_interval = get_z_interval();

    param.v3_length_computed.x() = x_interval * static_cast<double>(param.v3_len_interval_mult.x());
    param.v3_length_computed.y() = y_interval * static_cast<double>(param.v3_len_interval_mult.y());
    param.v3_length_computed.z() = z_interval * static_cast<double>(param.v3_len_interval_mult.z());

    LOG_INFO("  interval multipliers: ({}, {}, {})",
             param.v3_len_interval_mult.x(), param.v3_len_interval_mult.y(), param.v3_len_interval_mult.z());
    LOG_INFO("  grid intervals: ({:.6f}, {:.6f}, {:.6f}) m", x_interval, y_interval, z_interval);
    LOG_INFO("  computed cell dimensions: ({:.3f}, {:.3f}, {:.3f}) m",
             param.v3_length_computed.x(), param.v3_length_computed.y(), param.v3_length_computed.z());

    // Snap coordinates to grid if enabled
    if (param.tf_snap_to_grid) {
      // Snap center position to grid center
      Eigen::Vector3d& v3_cnt = param.v3_pos_cnt;
      v3_cnt.x() = get_xcnt(v3_cnt.x());
      v3_cnt.y() = get_ycnt(v3_cnt.y());
      v3_cnt.z() = get_zcnt(v3_cnt.z());
      LOG_INFO("  snapped center: ({:.3f}, {:.3f}, {:.3f}) m", v3_cnt.x(), v3_cnt.y(), v3_cnt.z());

      if (param.v3_len_cells == Eigen::Vector3i(0, 0, 0)) {
        // Snap AABB boundaries to grid boundaries
        AABB3d& aabb = param.aabb3d;
        const int ix_min = std::max(0, get_ix(aabb.xmin()));
        const int ix_max = std::min(get_nbinx() - 1, get_ix(aabb.xmax()));
        const int iy_min = std::max(0, get_iy(aabb.ymin()));
        const int iy_max = std::min(get_nbiny() - 1, get_iy(aabb.ymax()));
        const int iz_min = std::max(0, get_iz(aabb.zmin()));
        const int iz_max = std::min(get_nbinz() - 1, get_iz(aabb.zmax()));
        aabb.set_xmin(get_xlow(ix_min));
        aabb.set_xmax(get_xup(ix_max));
        aabb.set_ymin(get_ylow(iy_min));
        aabb.set_ymax(get_yup(iy_max));
        aabb.set_zmin(get_zlow(iz_min));
        aabb.set_zmax(get_zup(iz_max));
        LOG_INFO("  snapped AABB: x[{:.3f}, {:.3f}], y[{:.3f}, {:.3f}], z[{:.3f}, {:.3f}]",
                 aabb.xmin(), aabb.xmax(), aabb.ymin(), aabb.ymax(), aabb.zmin(), aabb.zmax());
      }
    }

    // Compute AABB from cells (after center snap)
    if (param.v3_len_cells != Eigen::Vector3i(0, 0, 0)) {
      const Eigen::Vector3d& v3_cnt = param.v3_pos_cnt;
      const Eigen::Vector3d& v3_cell_size = param.v3_length_computed;
      param.aabb3d.set_xmin(v3_cnt.x() - 0.5 * param.v3_len_cells.x() * v3_cell_size.x());
      param.aabb3d.set_xmax(v3_cnt.x() + 0.5 * param.v3_len_cells.x() * v3_cell_size.x());
      param.aabb3d.set_ymin(v3_cnt.y() - 0.5 * param.v3_len_cells.y() * v3_cell_size.y());
      param.aabb3d.set_ymax(v3_cnt.y() + 0.5 * param.v3_len_cells.y() * v3_cell_size.y());
      param.aabb3d.set_zmin(v3_cnt.z() - 0.5 * param.v3_len_cells.z() * v3_cell_size.z());
      param.aabb3d.set_zmax(v3_cnt.z() + 0.5 * param.v3_len_cells.z() * v3_cell_size.z());
      LOG_INFO("  AABB from cells ({},{},{}): x[{:.3f}, {:.3f}], y[{:.3f}, {:.3f}], z[{:.3f}, {:.3f}]",
               param.v3_len_cells.x(), param.v3_len_cells.y(), param.v3_len_cells.z(),
               param.aabb3d.xmin(), param.aabb3d.xmax(),
               param.aabb3d.ymin(), param.aabb3d.ymax(),
               param.aabb3d.zmin(), param.aabb3d.zmax());
    }

    // Compute elliptic cylinder radii and override AABB XY bounds if region_type == "cylinder"
    if (param.region_type == "cylinder") {
      const double radius_x = param.radius_x_meters;
      const double radius_y = param.radius_y_meters;
      // Override AABB XY with elliptic bounding box (used as quick rejection)
      param.aabb3d.set_xmin(param.v3_pos_cnt.x() - radius_x);
      param.aabb3d.set_xmax(param.v3_pos_cnt.x() + radius_x);
      param.aabb3d.set_ymin(param.v3_pos_cnt.y() - radius_y);
      param.aabb3d.set_ymax(param.v3_pos_cnt.y() + radius_y);
      LOG_INFO("  cylinder region: radius_x={:.3f} m, radius_y={:.3f} m",
               radius_x, radius_y);
      LOG_INFO("  cylinder bounding AABB: x[{:.3f}, {:.3f}], y[{:.3f}, {:.3f}], z[{:.3f}, {:.3f}]",
               param.aabb3d.xmin(), param.aabb3d.xmax(),
               param.aabb3d.ymin(), param.aabb3d.ymax(),
               param.aabb3d.zmin(), param.aabb3d.zmax());
    }

    add_density_structure(param);
  }

  // add ellipsoid density structure
  for( const auto& param : grdprm.vec_ellipsoid_prm ){
    if( param.tf_exec==false ) continue;
    LOG_INFO("adding ellipsoid density structure: {}", param.name);
    add_density_structure(param);
  }

  // add cylinder density structure
  for( const auto& param : grdprm.vec_cylinder_prm ){
    if( param.tf_exec==false ) continue;
    LOG_INFO("adding cylinder density structure: {}", param.name);
    add_density_structure(param);
  }

  // add cuboid density structure
  for( const auto& param : grdprm.vec_cuboid_prm ){
    if( param.tf_exec==false ) continue;
    LOG_INFO("adding cuboid density structure: {}", param.name);
    add_density_structure(param);
  }
}

// call Voxel function 1
Voxel& Grid3dVoxel::callVoxel( const int ix, const int iy, const int iz )
{
  // if out of range, throw error
  check_ix_inside(ix);
  check_iy_inside(iy);
  check_iz_inside(iz);
  return vec_vec_vec_Voxel.at(iz).at(iy).at(ix);
}

// call Voxel function 2
Voxel& Grid3dVoxel::callVoxel( const double x, const double y, const double z )
{
  auto [ix,iy,iz] = Grid3d::get_ixiyiz(x,y,z);
  return vec_vec_vec_Voxel.at(iz).at(iy).at(ix);
}

// call Voxel function 3
Voxel& Grid3dVoxel::callVoxel( const Uqiv uqiv_in )
{
  auto [ix,iy,iz] = get_ixiyiz(uqiv_in);
  // if out of range, throw error
  check_ix_inside(ix);
  check_iy_inside(iy);
  check_iz_inside(iz);
  return vec_vec_vec_Voxel.at(iz).at(iy).at(ix);
}

//
// get Voxel functions
//
const Voxel& Grid3dVoxel::getVoxel( const Uqiv uqiv_in ) const
{
  auto [ix,iy,iz] = get_ixiyiz(uqiv_in);
  // if out of range, throw error
  check_ix_inside(ix);
  check_iy_inside(iy);
  check_iz_inside(iz);
  return vec_vec_vec_Voxel.at(iz).at(iy).at(ix);
}

// get funtion
// if remove const before Voxel&, build error
const Voxel& Grid3dVoxel::getVoxel( const int ix, const int iy, const int iz ) const
{
  // if out of range, throw error
  check_ix_inside(ix);
  check_iy_inside(iy);
  check_iz_inside(iz);
  return vec_vec_vec_Voxel.at(iz).at(iy).at(ix);
}

// if remove const before Voxel&, build error
const Voxel& Grid3dVoxel::getVoxel( const double x, const double y, const double z ) const
{
  auto [ix,iy,iz] = Grid3d::get_ixiyiz(x,y,z);
  // if out of range, throw error
  check_ix_inside(ix);
  check_iy_inside(iy);
  check_iz_inside(iz);
  return vec_vec_vec_Voxel.at(iz).at(iy).at(ix);
}

bool Grid3dVoxel::is_ixiyiz_out_of_range( const int ix, const int iy, const int iz ) const
{
  if( ix < 0 ) return true;
  if( iy < 0 ) return true;
  if( iz < 0 ) return true;
  if( ix > get_nbinx()-1 ) return true;
  if( iy > get_nbiny()-1 ) return true;
  if( iz > get_nbinz()-1 ) return true;
  return false;
}

/// @brief Get voxels at a given z height as Grid2dVoxel
Grid2dVoxel Grid3dVoxel::get_Grid2dVoxel_z(const int iz) const
{
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  const double z_low = get_zlow(iz);
  const double z_up  = get_zup(iz);
  Grid2dVoxel g2vox(get_x_axis(),get_y_axis(),z_low,z_up,get_n_det());
  g2vox.set_vec_vec_Voxel( vec_vec_vec_Voxel.at(iz) );
  return g2vox;
}

/// @brief get the difference of density between two Grid3dVoxel
Grid3dVoxel Grid3dVoxel::get_delta_dens(const Grid3dVoxel &g3vox_in) const
{
  // check the size of two Grid3dVoxel
  if( !is_same_size(g3vox_in) ){
    LOG_ERROR(", size is different.");
    THROW_ERROR("Grid3dVoxel::get_delta_density, size is different.");
  }
  // make Grid3dVoxel from this
  Grid3dVoxel g3vox_delta(*this);
  // set name
  g3vox_delta.set_name("diff_"+get_name()+"_"+g3vox_in.get_name());
  // set density difference
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  const int nbinz = get_nbinz();
  #pragma omp parallel for collapse(3)
  for(int iz=0;iz<nbinz;iz++){
    for(int iy=0;iy<nbiny;iy++){
      for(int ix=0;ix<nbinx;ix++){
        const Voxel& vox1 =    this->getVoxel(ix,iy,iz);
        const Voxel& vox2 = g3vox_in.getVoxel(ix,iy,iz);
        const double delta_dens = vox1.get_density() - vox2.get_density();
        g3vox_delta.callVoxel(ix,iy,iz).set_density(delta_dens);
      }
    }
  }
  return g3vox_delta;
}



// allocate memory for std::vector<std::vector<Voxel>> vec_vec_vec_Voxel;
void Grid3dVoxel::vec_vec_vec_memory_allocate(const int n_det){
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  const int nbinz = get_nbinz();
  vec_vec_vec_Voxel.resize(nbinz); // allocate z direction
  for(int iz=0;iz<nbinz;iz++){
    if( iz%10 == 0 )fprintf(stderr,"progress %d/%d\r",iz,nbinz);
    vec_vec_vec_Voxel.at(iz).resize(nbiny); // allocate y direction
    for(int iy=0;iy<nbiny;iy++){
      vec_vec_vec_Voxel.at(iz).at(iy).resize(nbinx); // allocate x direction
      for(int ix=0;ix<nbinx;ix++){
        callVoxel(ix,iy,iz).set_tf_exist(false); // set default value "false"
        callVoxel(ix,iy,iz).set_density(0.); // set default value 0.0"
      }
    }
  }
}

// allocate memory for std::vector<std::vector<Voxel>> vec_vec_vec_Voxel;
void Grid3dVoxel::mp_vec_vec_vec_memory_allocate(const int n_det)
{
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  const int nbinz = get_nbinz();
  vec_vec_vec_Voxel.resize(nbinz); // allocate z direction
  #pragma omp parallel for
  for(int iz=0;iz<nbinz;iz++){
    if( iz%10 == 0 )fprintf(stderr,"progress %d/%d\r",iz,nbinz);
    vec_vec_vec_Voxel.at(iz).resize(nbiny); // allocate y direction  
    for(int iy=0;iy<nbiny;iy++){
      vec_vec_vec_Voxel.at(iz).at(iy).resize(nbinx); // allocate x direction
      for(int ix=0;ix<nbinx;ix++){
        callVoxel(ix,iy,iz).set_tf_exist(false); // set default value "false"
        callVoxel(ix,iy,iz).set_density(0.); // set default value 0.0"
      }
    }
  }
}

/// @brief bool check the size of two Grid3dVoxel
bool Grid3dVoxel::is_same_size( const Grid3dVoxel &g3vox_in ) const
{
  // check base class Grid3d
  if(this->getBaseRef() != g3vox_in.getBaseRef()){
    LOG_WARN(", Grid3d is different.");
    return false;
  }

  // check the size of z in vec_vec_vec_Voxel
  const int nbinz = vec_vec_vec_Voxel.size();
  const int nbinz_in = g3vox_in.vec_vec_vec_Voxel.size();
  if( nbinz != nbinz_in ){
    LOG_WARN(", vec_vec_vec_Voxel.size() is different.");
    return false;
  }

  // #pragma omp parallel for collapse(2) schedule(static)
  // ! Cannot use break in parallel for loop
  for(int iz=0;iz<nbinz;iz++){
    // check the size of y in vec_vec_vec_Voxel
    const int nbiny = vec_vec_vec_Voxel.at(iz).size();
    const int nbiny_in = g3vox_in.vec_vec_vec_Voxel.at(iz).size();
    if( nbiny != nbiny_in ){
      LOG_WARN(", vec_vec_vec_Voxel.at(iz={}).size() is different.",iz);
      return false;
    }
    for(int iy=0;iy<nbiny;iy++){
      // check the size of x in vec_vec_vec_Voxel
      const int nbinx = vec_vec_vec_Voxel.at(iz).at(iy).size();
      const int nbinx_in = g3vox_in.vec_vec_vec_Voxel.at(iz).at(iy).size();
      if( nbinx != nbinx_in ){
        LOG_WARN(", vec_vec_vec_Voxel.at(iz={},iy={}).size() is different.",iz,iy);
        return false;
      }
      for(int ix=0;ix<nbinx;ix++){
        // check voxel(ix,iy,iz).tf_exist is same
        const Voxel& vox = getVoxel(ix,iy,iz);
        const Voxel& vox_in = g3vox_in.getVoxel(ix,iy,iz);
        if( vox.get_tf_exist() != vox_in.get_tf_exist() ){
          LOG_WARN(", voxel(ix={},iy={},iz={}) tf_exist is different.",ix,iy,iz);
          return false;
        }
        // n_det is now a grid-level property; per-voxel check removed.
      }
    }
  }
  return true;
}

// output function
void Grid3dVoxel::out_voxel( FILE *fout, const int ix, const int iy,  const int iz) const
{
  const Voxel& vox = getVoxel(ix,iy,iz);
  const bool tf_exist = vox.get_tf_exist();
  const double x_low = get_xlow(ix);
  const double y_low = get_ylow(iy);
  const double z_low = get_zlow(iz);
  const double x_up  = get_xup(ix);
  const double y_up  = get_yup(iy);
  const double z_up  = get_zup(iz);
  const double density = vox.get_density();
  const Uqiv uqiv_ov = get_uqiv_fast(ix,iy,iz);
  const int  n_hit_ele = (uqiv_ov != UqivNotFound) ? static_cast<int>(get_n_hit_ele_grid(uqiv_ov)) : 0;
  const int  n_hit_det = (uqiv_ov != UqivNotFound) ? get_n_hit_det_grid(uqiv_ov) : 0;
  fprintf(fout,"%d %E %E %E %E %E %E %E %d %d\n"
  ,tf_exist,x_low,x_up,y_low,y_up,z_low,z_up,density,n_hit_ele,n_hit_det);
}

void Grid3dVoxel::out_voxel_with_index( FILE *fout, const int ix, const int iy,  const int iz) const
{
  const Voxel& vox = getVoxel(ix,iy,iz);
  fprintf(fout,"%d %d %d ",ix,iy,iz);
  bool tf_exist = vox.get_tf_exist();
  double x_low = get_xlow(ix);
  double y_low = get_ylow(iy);
  double z_low = get_zlow(iz);
  double x_up  = get_xup(ix);
  double y_up  = get_yup(iy);
  double z_up  = get_zup(iz);
  double density = vox.get_density();
  const Uqiv uqiv_ovi = get_uqiv_fast(ix,iy,iz);
  const int  n_hit_ele = (uqiv_ovi != UqivNotFound) ? static_cast<int>(get_n_hit_ele_grid(uqiv_ovi)) : 0;
  const int  n_hit_det = (uqiv_ovi != UqivNotFound) ? get_n_hit_det_grid(uqiv_ovi) : 0;
  fprintf(fout,"%d %E %E %E %E %E %E %E %d %d\n"
  ,tf_exist,x_low,x_up,y_low,y_up,z_low,z_up,density,n_hit_ele,n_hit_det);
}

void Grid3dVoxel::out_voxel_with_index( FILE *fout, const int unique_index_in ) const
{
  auto [ix,iy,iz] = get_ixiyiz(unique_index_in);
  out_voxel_with_index(fout,ix,iy,iz);
}


void Grid3dVoxel::out_all_exist( FILE *fout ) const
{
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  const int nbinz = get_nbinz();
  for(int iz=0;iz<nbinz;iz++){
    for(int iy=0;iy<nbiny;iy++){
      for(int ix=0;ix<nbinx;ix++){
        const Voxel& vox = getVoxel(ix,iy,iz);
        if( vox.get_tf_exist()==false ) continue;
        out_voxel(fout,ix,iy,iz);
      }
    }
  }
}

void Grid3dVoxel::out_all_exist( const fs::path& pathout ) const
{
  FILE *fout = myapp::get_fout(pathout);
  out_all_exist(fout);
  myapp::close(fout,pathout);
}

void Grid3dVoxel::out_all_with_index( FILE *fout ) const
{
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  const int nbinz = get_nbinz();
  for(int iz=0;iz<nbinz;iz++){
    for(int iy=0;iy<nbiny;iy++){
      for(int ix=0;ix<nbinx;ix++){
        out_voxel_with_index(fout,ix,iy,iz);
      }
    }
  }
}

void Grid3dVoxel::out_all_with_index( const fs::path& pathout ) const
{
  FILE *fout = myapp::get_fout(pathout);
  out_all_with_index(fout);
  myapp::close(fout,pathout);
}

void Grid3dVoxel::out_all_exist_with_index( FILE *fout ) const
{
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  const int nbinz = get_nbinz();
  for(int iz=0;iz<nbinz;iz++){
    for(int iy=0;iy<nbiny;iy++){
      for(int ix=0;ix<nbinx;ix++){
        const Voxel& vox = getVoxel(ix,iy,iz);
        if( vox.get_tf_exist()==false ) continue;
        out_voxel_with_index(fout,ix,iy,iz);
      }
    }
  }
}

void Grid3dVoxel::out_all_exist_with_index( const fs::path& pathout ) const
{
  FILE *fout = myapp::get_fout(pathout);
  out_all_exist_with_index(fout);
  myapp::close(fout,pathout);
}

void Grid3dVoxel::out_all_exist_with_index( FILE *fout, const Nhit n_hit_det_thres ) const
{
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  const int nbinz = get_nbinz();
  for(int iz=0;iz<nbinz;iz++){
    for(int iy=0;iy<nbiny;iy++){
      for(int ix=0;ix<nbinx;ix++){
        const Voxel& vox = getVoxel(ix,iy,iz);
        if( vox.get_tf_exist()==false ) continue;
        const Uqiv uqiv_ae = get_uqiv_fast(ix,iy,iz);
        if( uqiv_ae == UqivNotFound ) continue;
        if( get_n_hit_det_grid(uqiv_ae)<n_hit_det_thres ) continue;
        out_voxel_with_index(fout,ix,iy,iz);
      }
    }
  }
}

void Grid3dVoxel::out_all_exist_with_index(
  const fs::path& pathout, const Nhit n_hit_det_thres ) const
{
  FILE *fout = myapp::get_fout(pathout);
  out_all_exist_with_index(fout, n_hit_det_thres);
  myapp::close(fout,pathout);
}

void Grid3dVoxel::out_highest_exist_with_index( FILE *fout ) const
{
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  const int nbinz = get_nbinz();
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      bool flg_output = true;
      for(int iz=nbinz-1;iz>=0;iz--){
        const Voxel& vox = getVoxel(ix,iy,iz);
        if( vox.get_tf_exist()==false ) continue;
        if(flg_output==true){
          out_voxel_with_index(fout,ix,iy,iz);
          flg_output = false;
          break;
        }
      }
    }
  }
}

void Grid3dVoxel::out_highest_exist_with_index(
  const fs::path& pathout ) const
{
  FILE *fout = myapp::get_fout(pathout);
  out_highest_exist_with_index(fout);
  myapp::close(fout,pathout);
}

void Grid3dVoxel::out_lowest_exist_with_index(
  FILE *fout, const Nhit n_hit_det_thres ) const
{
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  const int nbinz = get_nbinz();
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      bool flg_output = true;
      for(int iz=0;iz<nbinz;iz++){
        const Voxel& vox = getVoxel(ix,iy,iz);
        if( vox.get_tf_exist()==false ) continue;
        const Uqiv uqiv_le = get_uqiv_fast(ix,iy,iz);
        if( uqiv_le == UqivNotFound ) continue;
        if( get_n_hit_det_grid(uqiv_le)<n_hit_det_thres ) continue;
        if(flg_output==true){
          out_voxel_with_index(fout,ix,iy,iz);
          flg_output = false;
        }
      }
    }
  }
}

void Grid3dVoxel::out_lowest_exist_with_index( const fs::path& pathout, const Nhit n_hit_det_thres ) const
{
  FILE *fout = myapp::get_fout(pathout);
  out_lowest_exist_with_index(fout,n_hit_det_thres);
  myapp::close(fout,pathout);
}

void Grid3dVoxel::out_all_with_ump_uqiv_ixiyiz( FILE *fout, const bool tf_exist_in ) const
{
  const Uqiv uqiv_min = get_uqiv_min();
  const Uqiv uqiv_max = get_uqiv_max();
  for(Uqiv uqiv=uqiv_min;uqiv<=uqiv_max;uqiv++){
    const Voxel& vox = getVoxel(uqiv);
    if(vox.get_tf_exist()!=tf_exist_in) continue;
    out_voxel_with_index(fout,uqiv);
  }
}

void Grid3dVoxel::out_all_with_ump_uqiv_ixiyiz(
  const fs::path& pathout, const bool tf_exist_in ) const
{
  FILE *fout = myapp::get_fout(pathout);
  out_all_with_ump_uqiv_ixiyiz(fout,tf_exist_in);
  myapp::close(fout,pathout);
}


void Grid3dVoxel::out_all_sort_zxy( FILE *fout ) const
{
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  const int nbinz = get_nbinz();
  // for(int iz=nbinz-1;iz>=0;iz--){
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      for(int iz=0;iz<nbinz;iz++){
        out_voxel_with_index(fout,ix,iy,iz);
      }
    }
  }
}

void Grid3dVoxel::out_all_sort_zxy( const fs::path& pathout ) const
{
  FILE *fout = myapp::get_fout(pathout);
  out_all_sort_zxy(fout);
  myapp::close(fout,pathout);
}

void Grid3dVoxel::out_all_sort_zxy_exist( FILE *fout ) const
{
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  const int nbinz = get_nbinz();
  // for(int iz=nbinz-1;iz>=0;iz--){
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      bool tf_flg_exist = false;
      for(int iz=0;iz<nbinz;iz++){
        const Voxel& vox = getVoxel(ix,iy,iz);
        if( vox.get_tf_exist()==false ) continue;
        out_voxel(fout,ix,iy,iz);
        tf_flg_exist = true;
      }
      if( tf_flg_exist==true ) fprintf(fout,"\n");
    }
  }
}

void Grid3dVoxel::out_all_sort_zxy_exist( const fs::path& pathout ) const
{
  FILE *fout = myapp::get_fout(pathout);
  out_all_sort_zxy_exist(fout);
  myapp::close(fout,pathout);
}

void Grid3dVoxel::out_all_sort_xyz_exist( FILE *fout ) const
{
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  const int nbinz = get_nbinz();
  // for(int iz=nbinz-1;iz>=0;iz--){
  for(int iz=0;iz<nbinz;iz++){
    for(int iy=0;iy<nbiny;iy++){
      bool tf_flg_exist = false;
      for(int ix=0;ix<nbinx;ix++){
        const Voxel& vox = getVoxel(ix,iy,iz);
        if( vox.get_tf_exist()==false ) continue;
        out_voxel(fout,ix,iy,iz);
        tf_flg_exist = true;
      }
      if( tf_flg_exist==true ) fprintf(fout,"\n");
    }
  }
}

void Grid3dVoxel::out_all_sort_xyz_exist( const fs::path& pathout ) const
{
  FILE *fout = myapp::get_fout(pathout);
  out_all_sort_xyz_exist(fout);
  myapp::close(fout,pathout);
}

// build vec_vec_vec_Voxel elements
void Grid3dVoxel::build_element( const Grid3dVoxel::Parameters &prm, const int n_det )
{
  
  // assignment of detector name
  name = prm.name;

  const int nbinx = (int) (( prm.xmax - prm.xmin )/prm.x_pitch);
  const int nbiny = (int) (( prm.ymax - prm.ymin )/prm.y_pitch);
  const int nbinz = (int) (( prm.zmax - prm.zmin )/prm.z_pitch);
  
  // assignment of x and y axis
  const std::string name_x_axis = "x_axis_" + name;
  const std::string name_y_axis = "y_axis_" + name;
  const std::string name_z_axis = "z_axis_" + name;
  
  // setting parameter of Grid1d x_axis and y_axis
  const Grid1d x_axis_tmp(name_x_axis, nbinx, prm.xmin, prm.xmax, prm.x_pitch);
  set_x_axis(x_axis_tmp);
  get_x_axis().out_info(spdlog::level::debug);

  const Grid1d y_axis_tmp(name_y_axis, nbiny, prm.ymin, prm.ymax, prm.y_pitch);
  set_y_axis(y_axis_tmp);
  get_y_axis().out_info(spdlog::level::debug);

  const Grid1d z_axis_tmp(name_z_axis, nbinz, prm.zmin, prm.zmax, prm.z_pitch);
  set_z_axis(z_axis_tmp);
  get_z_axis().out_info(spdlog::level::debug);

  mp_vec_vec_vec_memory_allocate(n_det);
}

// convert from Grid2dPillar
// todo Draw explanatory diagram.
void Grid3dVoxel::convert_from_Grid2dPillar(
    const Grid2dPillar &g2pil
  , const Grid3dVoxel::Parameters &grdprm
  , const int n_detector_in )
{

  const bool tf_use_g2_axis = grdprm.use_grid2d_of_g2pil;

  // setting of x-axis
  const std::string name_x_axis = "x_axis_" + g2pil.get_name();
  if( tf_use_g2_axis ){
    Grid1d g1x(g2pil.get_x_axis());
    g1x.set_name(name_x_axis);
    set_x_axis(g1x);
  }else{
    const int nbinx_tmp = grdprm.get_nbinx();
    const Grid1d g1x(name_x_axis,nbinx_tmp,grdprm.xmin,grdprm.xmax,grdprm.x_pitch);
    set_x_axis(g1x);
  }
  get_x_axis().out_info(spdlog::level::debug);

  // setting of y-axis
  const std::string name_y_axis = "y_axis_" + g2pil.get_name();
  if( tf_use_g2_axis ){
    Grid1d g1y(g2pil.get_y_axis());
    g1y.set_name(name_y_axis);
    set_y_axis(g1y);
  }else{
    const int nbiny_tmp = grdprm.get_nbiny();
    const Grid1d g1y(name_y_axis,nbiny_tmp,grdprm.ymin,grdprm.ymax,grdprm.y_pitch);
    set_y_axis(g1y);
  }
  get_y_axis().out_info(spdlog::level::debug);

  // setting of z-axis
  const std::string name_z_axis = "z_axis_" + g2pil.get_name();
  const int nbinz = grdprm.get_nbinz();
  const Grid1d g1z(name_z_axis
  ,nbinz,grdprm.zmin,grdprm.zmax,grdprm.z_pitch,true);
  set_z_axis(g1z);
  g1z.out_info(spdlog::level::debug);

  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();

  // allocate memory
  mp_vec_vec_vec_memory_allocate(n_detector_in);

  // convert from Grid2dPillar
  for(int iy=0;iy<nbiny;iy++){
    fprintf(stderr,"progress : iy=%d/nbiny=%d\r",iy,nbiny);
    double ycnt = get_ycnt(iy);
    for(int ix=0;ix<nbinx;ix++){
      // center_value = min + ((double)index + 0.5)*interval;
      double xcnt = get_xcnt(ix);
      // const Pillar& cub = g2pil.getPillar(ix,iy);
      // ! The index ix,iy of g3vox is not the index of g2pil, so should be like below
      const Pillar& cub = g2pil.getPillar(xcnt,ycnt);
      double zmin_cub = cub.get_zmin();
      double zmax_cub = cub.get_zmax();
      double density_cub = cub.get_density();
      int izmin = get_iz(zmin_cub);
      if( !is_iz_above_lower_limit(izmin) ) izmin = 0;
      int izmax = get_iz(zmax_cub);
      if( !is_iz_below_upper_limit(izmax) ) izmax = nbinz-1;
      // for(iz=izmin;iz<=izmax;iz++){
      for(int iz=0;iz<nbinz;iz++){
        double z_lower_vox = get_zlow(iz);
        
        // If voxel z-axis index is within cuboid's z min-max range:
        if(izmin<=iz&&iz<izmax){
          callVoxel(ix,iy,iz).set_tf_exist(true); // Material exists
          callVoxel(ix,iy,iz).set_density(density_cub); // Set density to density_cub
        }
        // If cuboid zmax falls between voxel z-axis bin edges:
        else if( izmax == iz ){
          // Calculate what percentage of cuboid fits
          double height_fill_ratio = (zmax_cub - z_lower_vox) / grdprm.z_pitch;
          if( height_fill_ratio < 0.0 ){
            LOG_CRITICAL("height_fill_ratio < 0.00");
            LOG_CRITICAL("ix={}, iy={}, iz={}, zmax_cub={:.1f}, z_lower_vox={:.1f}",ix,iy,iz,zmax_cub,z_lower_vox);
            THROW_ERROR_NAME("Grid3dVoxel::convert_from_Grid2dPillar, height_fill_ratio < 0.00");
          }
          if( height_fill_ratio > 1.0 ){
            LOG_CRITICAL("height_fill_ratio > 1.00");
            LOG_CRITICAL("ix={}, iy={}, iz={}, zmax_cub={:.1f}, z_lower_vox={:.1f}",ix,iy,iz,zmax_cub,z_lower_vox);
            LOG_CRITICAL("please make [GRID3D_VOXEL_PARAMETERS] zmax more larger");
            FILE *fout = fopen("delete_this_file_before_run.tmp","wt");
            if(fout == NULL) THROW_ERROR_NAME("Cannot open file : delete_this_file_before_run.tmp");
            fprintf(fout,"height_fill_ratio > 1.00\n");
            fprintf(fout,"please make [GRID3D_VOXEL_PARAMETERS] zmax more larger\n");
            if(fclose(fout) == EOF) THROW_ERROR_NAME("fclose(fout) == EOF");
            THROW_ERROR_NAME2("Grid3dVoxel::convert_from_Grid2dPillar, height_fill_ratio > 1.00",height_fill_ratio);
          }
          if( height_fill_ratio == 0.0 ){
            callVoxel(ix,iy,iz).set_density(0.);
            callVoxel(ix,iy,iz).set_tf_exist(false);
          }
          // Calculate density for portion that didn't fit in g3vox
          if( height_fill_ratio > 0.0 && height_fill_ratio <= 1.0 ){
            // Set as no material (tf_exist=false).
            // This way it won't be counted when calculating path length.
            // * 2023-09-15 13:56:33 Keep the following comment for now. Skipping due to time constraints.
            // In Grid3dVoxel::add_density_structure,
            // instead of adding/subtracting to give 500, 1500 from uniform_prior_density=1000,
            // if we multiply to give 500, 1500,
            // Voxels with 0.0<height_ratio<1.0 would probably be fine with tf_exist=true
            // ! Leaving it as below is questionable, but keeping it for now.
            callVoxel(ix,iy,iz).set_tf_exist(false); // Set as no material
            callVoxel(ix,iy,iz).set_density(density_cub*height_fill_ratio); // Density is non-zero.
          }
        }
        else{
          callVoxel(ix,iy,iz).set_tf_exist(false);
          callVoxel(ix,iy,iz).set_density(0.);
        }
      } // iz loop end
    } // ix loop end
  } // iy loop end
}


// 2022-11-09 12:27:45
// count num of voxel with tf_exist=true
int Grid3dVoxel::count_n_voxel( const bool tf_exist_in ) const
{
  int n_voxel = 0;
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  const int nbinz = get_nbinz();
  for(int iz=0; iz < nbinz; iz++){
    for(int iy=0; iy < nbiny; iy++){
      for(int ix=0; ix < nbinx; ix++){
        if( getVoxel(ix,iy,iz).get_tf_exist()==tf_exist_in ) n_voxel++;
      }
    }
  }
  return n_voxel;
}

// copy Grid3dVoxel
void Grid3dVoxel::set(const Grid3dVoxel &org)
{
  *this = org;
}

// Inequality operator
bool Grid3dVoxel::operator!=(const Grid3dVoxel& other) const
{
  // Name matching is not checked.
  // if (name != other.name) return true;

  if (n_vox_exist != other.n_vox_exist) return true;

  // Compare vec_vec_vec_Voxel (detailed comparison logic may be needed)
  if (!is_vec_vec_vec_Voxel_same(other)) return true;

  // Compare id_container::VoxID
  if (uqiv_container != other.uqiv_container) return true;

  return false;
}

// 2023-03-17 17:22:15
// dimension of density is SI unit. i.e. kg/m3
void Grid3dVoxel::set_uniform_density_uqiv( const double density_in )
{
  const int nvox = get_n_vox_exist();
  for(Uqiv uqiv=0;uqiv<nvox;uqiv++){
    fprintf(stderr,"set_uniform_density : uqiv=%d / %d...\r",uqiv,nvox);
    Voxel& vox = callVoxel(uqiv);
    if(vox.get_tf_exist()==false) continue;
    vox.set_density(density_in);
  }
}

// dimension of density is SI unit. i.e. kg/m3
void Grid3dVoxel::mp_set_uniform_density_uqiv( const double density_in )
{
  const int nvox = get_n_vox_exist();
  const int n_threads = omp_get_max_threads();
  #pragma omp parallel for
  for(Uqiv uqiv=0;uqiv<nvox;uqiv++){
    if( uqiv%n_threads==0 )
      fprintf(stderr,"set_uniform_density : uqiv=%d / %d...\r",uqiv,nvox);
    Voxel& vox = callVoxel(uqiv);
    if(vox.get_tf_exist()==false) continue;
    vox.set_density(density_in);
  }
}

// dimension of density is SI unit. i.e. kg/m3
void Grid3dVoxel::set_uniform_density_ixiyiz( const double density_in )
{
  const int nbinz = get_nbinz();
  const int nbiny = get_nbiny();
  const int nbinx = get_nbinx();
  #pragma omp parallel for collapse(3)
  for(int iz=0;iz<nbinz;iz++){
    for(int iy=0;iy<nbiny;iy++){
      for(int ix=0;ix<nbinx;ix++){
        Voxel& vox = callVoxel(ix,iy,iz);
        if(vox.get_tf_exist()==false) continue;
        vox.set_density(density_in);
      }
    }
  }
}

// dimension of density is SI unit. i.e. kg/m3
void Grid3dVoxel::mp_set_uniform_density_ixiyiz( const double density_in )
{
  const int nbinz = get_nbinz();
  const int nbiny = get_nbiny();
  const int nbinx = get_nbinx();
  #pragma omp parallel for collapse(3)
  for(int iz=0;iz<nbinz;iz++){
    for(int iy=0;iy<nbiny;iy++){
      for(int ix=0;ix<nbinx;ix++){
        Voxel& vox = callVoxel(ix,iy,iz);
        if(vox.get_tf_exist()==false) continue;
        vox.set_density(density_in);
      }
    }
  }
}

// build ump_uqiv_ixiyiz and ump_ixiyiz_uqiv
// for only voxels which tf_exist==tf_exist_in
// returns n_voxel_exist
int Grid3dVoxel::build_uqiv_umps( const bool tf_exist_in, const Uqiv uqiv_start )
{
  // count nbinx,nbiny,nbinz
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  const int nbinz = get_nbinz();
  
  // Clear existing maps before rebuilding.  Required when called a second
  // time (e.g. after AABB mask in merge()) — without this, stale uqiv
  // values survive because unordered_map::insert ignores duplicate keys.
  clear_uqiv_umps();

  // allocate memory
  reserve_uqiv_umps(nbinx*nbiny*nbinz);
  
  // set uqiv
  Uqiv uqiv = uqiv_start;
  set_uqiv_min( uqiv );

  int n_vox_count = 0;
  for(int iz=0; iz < nbinz; iz++){
    for(int iy=0; iy < nbiny; iy++){
      for(int ix=0; ix < nbinx; ix++){
        if( getVoxel(ix,iy,iz).get_tf_exist()!=tf_exist_in ) continue;

        // below is only tf_exist is true
        const Ixiyiz ixiyiz = {ix,iy,iz};
        
        // insert (key,value) to (ump_uqiv_ixiyiz<ump_ixiyiz_uqiv)
        insert_to_uqiv_umps( uqiv, ixiyiz );

        uqiv++;
        n_vox_count++;
      }
    }
  }
  set_uqiv_max( uqiv - 1 );
  LOG_DEBUG("uqiv_start={}, uqiv_max={}",uqiv_start,uqiv-1); 

  set_n_vox_exist(n_vox_count);
  LOG_DEBUG("n_vox_count={}",n_vox_count);

  // build flat array for O(1) lookup — replaces hash-map access in hot path
  build_flat_uqiv();

  // allocate hit vectors after uqiv system is built
  allocate_hit_data();

  return uqiv;
}

// Build flat_uqiv_ array from the unordered_map for O(1) direct-index lookup.
// Eliminates hash computation + pointer chasing in the ray-tracing hot path.
void Grid3dVoxel::build_flat_uqiv()
{
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  const int nbinz = get_nbinz();
  const std::size_t total = static_cast<std::size_t>(nbinx) * nbiny * nbinz;

  // fill with sentinel — non-existing voxels remain UqivNotFound
  flat_uqiv_.assign(total, UqivNotFound);

  // populate from the authoritative unordered_map
  const auto& ump = get_ump_ixiyiz_uqiv_ref();
  for (const auto& [ixiyiz, uqiv] : ump) {
    const std::size_t idx =
        static_cast<std::size_t>(ixiyiz[0])
      + static_cast<std::size_t>(ixiyiz[1]) * nbinx
      + static_cast<std::size_t>(ixiyiz[2]) * nbinx * nbiny;
    flat_uqiv_[idx] = uqiv;
  }
}

// O(1) uqiv lookup via flat array — avoids unordered_map hash + pointer chase.
// Falls back to the hash-map path if flat_uqiv_ has not been built.
Grid3d::Uqiv Grid3dVoxel::get_uqiv_fast(
    const int ix_in, const int iy_in, const int iz_in ) const
{
  if (flat_uqiv_.empty()) {
    // fallback: flat array not yet built
    return get_uqiv(ix_in, iy_in, iz_in);
  }
  const std::size_t idx =
      static_cast<std::size_t>(ix_in)
    + static_cast<std::size_t>(iy_in) * get_nbinx()
    + static_cast<std::size_t>(iz_in) * get_nbinx() * get_nbiny();
  return flat_uqiv_[idx];
}

// Allocate hit-data vectors based on uqiv range.
void Grid3dVoxel::allocate_hit_data()
{
  const auto uqiv_min = get_uqiv_min();
  const auto uqiv_max = get_uqiv_max();
  if (uqiv_min > uqiv_max) {
    THROW_ERROR("Grid3dVoxel::allocate_hit_data: invalid uqiv range. min={}, max={}", uqiv_min, uqiv_max);
  }
  const size_t n = static_cast<size_t>(uqiv_max - uqiv_min + 1);
  vec_hit_det_.assign(n, 0ULL);
  vec_n_hit_ele_.assign(n, 0);
}

// Zero-fill hit-data vectors.
void Grid3dVoxel::clear_hit_data()
{
  std::fill(vec_hit_det_.begin(), vec_hit_det_.end(), 0ULL);
  std::fill(vec_n_hit_ele_.begin(), vec_n_hit_ele_.end(), 0);
}

// Record a detector hit (thread-safe via atomic_ref).
void Grid3dVoxel::record_hit_det(Uqiv uqiv, int detid)
{
  const auto uqiv_min = get_uqiv_min();
  std::atomic_ref<uint64_t>(vec_hit_det_.at(uqiv - uqiv_min))
    .fetch_or(1ULL << detid, std::memory_order_relaxed);
}

// Increment element-hit counter (thread-safe via atomic_ref).
void Grid3dVoxel::record_hit_ele(Uqiv uqiv)
{
  const auto uqiv_min = get_uqiv_min();
  std::atomic_ref<Nhit>(vec_n_hit_ele_.at(uqiv - uqiv_min))
    .fetch_add(1, std::memory_order_relaxed);
}

// Return raw hit bitmask for a voxel.
uint64_t Grid3dVoxel::get_hit_det(Uqiv uqiv) const
{
  const auto uqiv_min = get_uqiv_min();
  return vec_hit_det_.at(uqiv - uqiv_min);
}

// Return true if a specific detector has a hit.
bool Grid3dVoxel::get_tf_hit_grid(Uqiv uqiv, int det_id) const
{
  return (get_hit_det(uqiv) >> det_id) & 1ULL;
}

// Return the number of detectors that hit the voxel.
int Grid3dVoxel::get_n_hit_det_grid(Uqiv uqiv) const
{
  // __builtin_popcountll: GCC/Clang intrinsic that counts the number of
  // set bits (1-bits) in a 64-bit integer. Maps to hardware POPCNT instruction.
  return __builtin_popcountll(get_hit_det(uqiv));
}

// Return the element-hit count for a voxel.
Nhit Grid3dVoxel::get_n_hit_ele_grid(Uqiv uqiv) const
{
  const auto uqiv_min = get_uqiv_min();
  return vec_n_hit_ele_.at(uqiv - uqiv_min);
}

Grid3d::Uqiv Grid3dVoxel::get_uqiv( const int ix_in, const int iy_in, const int iz_in ) const
{
  // if out of range, throw error
  check_ix_inside(ix_in);
  check_iy_inside(iy_in);
  check_iz_inside(iz_in);
  return uqiv_container.get_uqiv(ix_in,iy_in,iz_in);
}

// merge voxels type2
// returns merged Grid3dVoxel
// in the remain instance, tf_exist become false in the merged voxels
// This function merges (combines) 3D grid voxel data (Grid3dVoxel)
// to create a new Grid3dVoxel object.
// Merging is performed based on specified offset and multiplier.
// The number of bins along each axis of the merged voxels decreases, and the voxel size increases.
Grid3dVoxel Grid3dVoxel::merge( 
    const double xcnt, const int merge_factor_x_in
  , const double ycnt, const int merge_factor_y_in
  , const double zcnt, const int merge_factor_z_in)
{
  int n_merged = 0;
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  const int nbinz = get_nbinz();
  // Check merge_factor_x, y, z. Values are adjusted if too large.
  int merge_factor_x = merge_factor_x_in;
  int merge_factor_y = merge_factor_y_in;
  int merge_factor_z = merge_factor_z_in;
  // Create Grid1d objects for each axis after merging.
  const Grid1d x_axis_merged = get_x_axis().get_merged(xcnt,merge_factor_x);
  const Grid1d y_axis_merged = get_y_axis().get_merged(ycnt,merge_factor_y);
  const Grid1d z_axis_merged = get_z_axis().get_merged(zcnt,merge_factor_z);

  // Minimum x, y, z values after merge
  const double xmin = x_axis_merged.get_min();
  const double ymin = y_axis_merged.get_min();
  const double zmin = z_axis_merged.get_min();
  // ixmin, iymin, izmin represent the pre-merge Grid1d index
  // corresponding to the post-merge minimum (xmin, ymin, zmin).
  const int ixmin = get_x_axis().get_index(xmin);
  const int iymin = get_y_axis().get_index(ymin);
  const int izmin = get_z_axis().get_index(zmin);

// debug
  const auto log_level = spdlog::level::debug;
  LOG_DEBUG(" before_merge :");
  get_x_axis().out_info(log_level);
  get_y_axis().out_info(log_level);
  get_z_axis().out_info(log_level);

  LOG_DEBUG(" after_merge :");
  x_axis_merged.out_info(log_level);
  y_axis_merged.out_info(log_level);
  z_axis_merged.out_info(log_level);
  // object which will be returned
  // Create merged Grid3dVoxel object.
  const int n_det = this->get_n_det();
  Grid3dVoxel g3vox_merged_input(
    x_axis_merged,y_axis_merged,z_axis_merged, n_det
  );

  // Determine endpoints for ix/iy/iz and their sub-index loops.
  // If nbin < merge_factor, set ix_end to finish in one loop iteration.
  // Otherwise, loop from imin to nbin-merge_factor.
  int ix_end, nbinx_sub;
  if( nbinx < merge_factor_x ){
    ix_end = ixmin + 1;
    nbinx_sub = nbinx - ixmin;
  }
  else{ 
    ix_end = nbinx - merge_factor_x;
    nbinx_sub = merge_factor_x;
  }

  int iy_end, nbiny_sub;
  if( nbiny < merge_factor_y ){
    iy_end = iymin + 1;
    nbiny_sub = nbiny - iymin;
  }
  else{ 
    iy_end = nbiny - merge_factor_y;
    nbiny_sub = merge_factor_y;
  }

  int iz_end, nbinz_sub;
  if( nbinz < merge_factor_z ){
    iz_end = izmin + 1;
    nbinz_sub = nbinz - izmin;
  }
  else{ 
    iz_end = nbinz - merge_factor_z;
    nbinz_sub = merge_factor_z;
  }

  // Loop over merged-grid steps; ix/iy/iz refer to pre-merge indices.
  int ix_sub, iy_sub, iz_sub;
  // double n_exist = 0.0, sum_density = 0.0;
  
  // Parallelize
  #pragma omp parallel for \
    private(ix_sub,iy_sub,iz_sub) \
    reduction(+: n_merged)
    // reduction(+: n_exist, sum_density, n_merged)
  
  //
  // Loop over pre-merge indices with post-merge Grid step
  //
  for(int iz=izmin; iz < iz_end; iz += merge_factor_z ){
    // const double zlow_before = get_zlow(iz);
    const double zcnt_before = get_zcnt(iz);
    for(int iy=iymin; iy < iy_end; iy += merge_factor_y ){
      // const double ylow_before = get_ylow(iy);
      const double ycnt_before = get_ycnt(iy);
      for(int ix=ixmin; ix < ix_end; ix += merge_factor_x ){
        //
        // Code within this block is executed by a single thread.
        //

        // const double xlow_before = get_xlow(ix);
        const double xcnt_before = get_xcnt(ix);
        fprintf(stderr,"Grid3dVoxel::merge ix=%d, iy=%d, iz=%d....\r",ix,iy,iz);
        
        // start of local voxel loop
        // Initialize n_exist
        double n_exist = 0.0;
// printf("ix=%d, iy=%d, iz=%d \n",ix,iy,iz);

        // Loop within the set of voxels to be merged
        for(iz_sub=0; iz_sub < nbinz_sub; iz_sub++ ){
          const int izz = iz + iz_sub;
          for(iy_sub=0; iy_sub < nbiny_sub; iy_sub++ ){
            const int iyy = iy + iy_sub;
            for(ix_sub=0; ix_sub < nbinx_sub; ix_sub++ ){
              const int ixx = ix + ix_sub;
              // ixx, iyy, izzis pre-merge index
              const Voxel& vox_tmp = getVoxel(ixx,iyy,izz);
              // Count number of tf_exist==true
              if( vox_tmp.get_tf_exist()) n_exist+=1.0;
            } // end ix_sub loop
          } // end iy_sub loop
        } // end iz_sub loop

        // If all merge targets are true, i.e. 
        // (int)n_exist == nbinx_sub * nbiny_sub * nbinz_sub
        // then change all to false, set g3vox_input tf_exist to true, and average density
        // Do nothing if condition is not met.
        if( (int)n_exist == nbinx_sub * nbiny_sub * nbinz_sub ){
          n_merged++;
          
          // Initialize sum_density
          double sum_density = 0.0;
          
          for(int iz_sub=0; iz_sub < nbinz_sub; iz_sub++ ){
            const int izz = iz + iz_sub;
            for(int iy_sub=0; iy_sub < nbiny_sub; iy_sub++ ){
              const int iyy = iy + iy_sub;
              for(int ix_sub=0; ix_sub < nbinx_sub; ix_sub++ ){
                const int ixx = ix + ix_sub;
                // ixx, iyy, izzis pre-merge index
                
                //Calculate density sum for later averaging
                const double density = getVoxel(ixx,iyy,izz).get_density();
                sum_density += density;

                // Set all g3vox to be merged to tf_exist=false.
                callVoxel(ixx,iyy,izz).set_tf_exist(false);
              } // end ix_sub loop
            } // end iy_sub loop
          } // end iz_sub loop
          // Get post-merge g3vox index from pre-merge index.
          const int ix_merged = x_axis_merged.get_index(xcnt_before);
          const int iy_merged = y_axis_merged.get_index(ycnt_before);
          const int iz_merged = z_axis_merged.get_index(zcnt_before);
          // const int ix_merged_low = x_axis_merged.get_index(xlow_before);
          // const int iy_merged_low = y_axis_merged.get_index(ylow_before);
          // const int iz_merged_low = z_axis_merged.get_index(zlow_before);
          Voxel& vox_merged = g3vox_merged_input.callVoxel(ix_merged,iy_merged,iz_merged);
          // Set post-merge tf_exist to true.
          vox_merged.set_tf_exist(true);
          
          // Average post-merge density.
          const double avr_density = sum_density/n_exist;
          vox_merged.set_density(avr_density);
        }
        // End of if "all merge targets are true"
        // end of local voxel loop

        //
        // Code up to here is executed by a single thread
        //

      } // end ix loop stepping by merged-grid stride (pre-merge index)
    } // end iy loop stepping by merged-grid stride (pre-merge index)
  } // end iz loop stepping by merged-grid stride (pre-merge index)
  if( n_merged == 0 ){
    THROW_ERROR_NAME(
      "Grid3dVoxel::merge, n_merged == 0\n"
      " Please check [GRID3D_VOXEL_MERGE_PARAMETERS_??] again.\n"
      " xcnt = {}\n"
      " merge_factor_x = {}\n"
      " ycnt = {}\n"
      " merge_factor_y = {}\n"
      " zcnt = {}\n"
      " merge_factor_z = {}\n",
      xcnt, merge_factor_x, ycnt, merge_factor_y, zcnt, merge_factor_z
    );
  }

  // build unique index map
  g3vox_merged_input.build_uqiv_umps();

  LOG_INFO(", n_merged = {}",n_merged);
  return g3vox_merged_input;
}

// @brief merge voxels type2 with MergeParameters and ReconstVoxelsParameters
// @param prm_merge is Grid3dVoxel::MergeParameters
// @param prm_reconst is Grid3dVoxel::ReconstVoxelsParameters
// @details calls merge voxels type2, then applies AABB and/or cylinder masks
Grid3dVoxel Grid3dVoxel::merge(
    const Grid3dVoxel::MergeParameters &prm_merge,
    const Grid3dVoxel::ReconstVoxelsParameters &prm_reconst )
{
  Grid3dVoxel result = merge(
      prm_merge.x_merge_center, prm_merge.x_merge_factor
    , prm_merge.y_merge_center, prm_merge.y_merge_factor
    , prm_merge.z_merge_center, prm_merge.z_merge_factor );

  const bool need_mask = prm_reconst.tf_aabb || prm_reconst.tf_cylinder;

  // Apply AABB mask if enabled
  if (prm_reconst.tf_aabb) {
    const double aabb_half_x = prm_reconst.x_aabb_meters * 0.5;
    const double aabb_half_y = prm_reconst.y_aabb_meters * 0.5;

    const double resolved_zmin = (prm_reconst.aabb_zmin_mode == "g3vox_zmin")
      ? result.get_zmin()
      : prm_reconst.aabb_zmin_value;

    const AABB3d aabb(
      Eigen::Vector3d(prm_reconst.x_aabb_cnt - aabb_half_x,
                      prm_reconst.y_aabb_cnt - aabb_half_y,
                      resolved_zmin),
      Eigen::Vector3d(prm_reconst.x_aabb_cnt + aabb_half_x,
                      prm_reconst.y_aabb_cnt + aabb_half_y,
                      prm_reconst.aabb_zmax));

    LOG_INFO("Applying AABB mask after merge: aabb=({},{},{}) - ({},{},{})",
             aabb.xmin(), aabb.ymin(), aabb.zmin(),
             aabb.xmax(), aabb.ymax(), aabb.zmax());

    const int n_before = result.get_n_vox_exist();
    const int nbinx = result.get_nbinx();
    const int nbiny = result.get_nbiny();
    const int nbinz = result.get_nbinz();

    for (int iz = 0; iz < nbinz; iz++) {
      for (int iy = 0; iy < nbiny; iy++) {
        for (int ix = 0; ix < nbinx; ix++) {
          if (!result.getVoxel(ix, iy, iz).get_tf_exist()) continue;
          const Eigen::Vector3d pos(
            result.get_xcnt(ix), result.get_ycnt(iy), result.get_zcnt(iz));
          if (!aabb.is_inside(pos)) {
            result.callVoxel(ix, iy, iz).set_tf_exist(false);
          }
        }
      }
    }

    const int n_after = result.get_n_vox_exist();
    LOG_INFO("AABB mask applied: n_vox_exist {} -> {} (removed {})",
             n_before, n_after, n_before - n_after);
  }

  // Apply elliptic cylinder mask if enabled
  if (prm_reconst.tf_cylinder) {
    const double radius_x = prm_reconst.cylinder_radius_x_meters;
    const double radius_y = prm_reconst.cylinder_radius_y_meters;
    const double inv_rx_sq = (radius_x > 0.0) ? 1.0 / (radius_x * radius_x) : 0.0;
    const double inv_ry_sq = (radius_y > 0.0) ? 1.0 / (radius_y * radius_y) : 0.0;
    const double cx = prm_reconst.x_cyl_cnt;
    const double cy = prm_reconst.y_cyl_cnt;

    LOG_INFO("Applying cylinder mask after merge: center=({},{}), radius=({},{})",
             cx, cy, radius_x, radius_y);

    const int n_before = result.get_n_vox_exist();
    const int nbinx = result.get_nbinx();
    const int nbiny = result.get_nbiny();
    const int nbinz = result.get_nbinz();

    for (int iz = 0; iz < nbinz; iz++) {
      for (int iy = 0; iy < nbiny; iy++) {
        for (int ix = 0; ix < nbinx; ix++) {
          if (!result.getVoxel(ix, iy, iz).get_tf_exist()) continue;
          const double dx = result.get_xcnt(ix) - cx;
          const double dy = result.get_ycnt(iy) - cy;
          if (dx * dx * inv_rx_sq + dy * dy * inv_ry_sq > 1.0) {
            result.callVoxel(ix, iy, iz).set_tf_exist(false);
          }
        }
      }
    }

    const int n_after = result.get_n_vox_exist();
    LOG_INFO("Cylinder mask applied: n_vox_exist {} -> {} (removed {})",
             n_before, n_after, n_before - n_after);
  }

  // Rebuild unique index maps once after all masks
  // (6-arg merge already called build_uqiv_umps(),
  // but tf_exist changes require a fresh rebuild)
  if (need_mask) {
    result.build_uqiv_umps();
  }

  return result;
}

Grid3dVoxel Grid3dVoxel::overlay_merged_density(
    const Grid3dVoxel& g3vox_merged) const
{
  Grid3dVoxel result(*this);

  const Grid1d& x_org = get_x_axis();
  const Grid1d& y_org = get_y_axis();
  const Grid1d& z_org = get_z_axis();
  const Grid1d& x_mrg = g3vox_merged.get_x_axis();
  const Grid1d& y_mrg = g3vox_merged.get_y_axis();
  const Grid1d& z_mrg = g3vox_merged.get_z_axis();

  int n_overlaid = 0;
  const int nbinx_m = g3vox_merged.get_nbinx();
  const int nbiny_m = g3vox_merged.get_nbiny();
  const int nbinz_m = g3vox_merged.get_nbinz();

  #pragma omp parallel for collapse(3) reduction(+: n_overlaid)
  for (int iz_m = 0; iz_m < nbinz_m; iz_m++) {
    for (int iy_m = 0; iy_m < nbiny_m; iy_m++) {
      for (int ix_m = 0; ix_m < nbinx_m; ix_m++) {
        const auto [iz_min, iz_max] = z_org.get_original_index_min_max(z_org, z_mrg, iz_m);
        const auto [iy_min, iy_max] = y_org.get_original_index_min_max(y_org, y_mrg, iy_m);
        const auto& vox_m = g3vox_merged.getVoxel(ix_m, iy_m, iz_m);
        if (!vox_m.get_tf_exist()) continue;

        const double density = vox_m.get_density();
        const auto [ix_min, ix_max] = x_org.get_original_index_min_max(x_org, x_mrg, ix_m);

        for (int iz = iz_min; iz <= iz_max; iz++) {
          for (int iy = iy_min; iy <= iy_max; iy++) {
            for (int ix = ix_min; ix <= ix_max; ix++) {
              auto& vox = result.callVoxel(ix, iy, iz);
              vox.set_density(density);
              vox.set_tf_exist(true);
            }
          }
        }
        n_overlaid++;
      }
    }
  }

  LOG_INFO("overlay_merged_density: overlaid {} merged voxels onto {} pre-merge grid",
           n_overlaid, result.get_name());

  return result;
}

Grid3dVoxel Grid3dVoxel::overlay_merged_density(
    const Grid3dVoxel& g3vox_merged_all,
    const Grid3dVoxel& g3vox_rec,
    double uniform_prior_density) const
{
  Grid3dVoxel result(*this);

  // Replace checker density with uniform_prior_density for all existing voxels
  result.mp_set_uniform_density_ixiyiz(uniform_prior_density);

  // Reset tf_exist for all voxels so that apply_shell_density
  // treats outer mountain voxels as shell candidates.
  // Pass 1/2 below will restore tf_exist=true only for reconst_volume.
  {
    const int nbz = result.get_nbinz();
    const int nby = result.get_nbiny();
    const int nbx = result.get_nbinx();
    #pragma omp parallel for collapse(3)
    for (int iz = 0; iz < nbz; iz++) {
      for (int iy = 0; iy < nby; iy++) {
        for (int ix = 0; ix < nbx; ix++) {
          result.callVoxel(ix, iy, iz).set_tf_exist(false);
        }
      }
    }
  }

  const Grid1d& x_org = get_x_axis();
  const Grid1d& y_org = get_y_axis();
  const Grid1d& z_org = get_z_axis();

  // --- Pass 1: Fill reconst_volume with uniform_prior_density ---
  {
    const Grid1d& x_mrg = g3vox_merged_all.get_x_axis();
    const Grid1d& y_mrg = g3vox_merged_all.get_y_axis();
    const Grid1d& z_mrg = g3vox_merged_all.get_z_axis();

    int n_prior = 0;
    const int nbinx_m = g3vox_merged_all.get_nbinx();
    const int nbiny_m = g3vox_merged_all.get_nbiny();
    const int nbinz_m = g3vox_merged_all.get_nbinz();

    #pragma omp parallel for collapse(3) reduction(+: n_prior)
    for (int iz_m = 0; iz_m < nbinz_m; iz_m++) {
      for (int iy_m = 0; iy_m < nbiny_m; iy_m++) {
        for (int ix_m = 0; ix_m < nbinx_m; ix_m++) {
          if (!g3vox_merged_all.getVoxel(ix_m, iy_m, iz_m).get_tf_exist()) continue;

          const auto [iz_min, iz_max] = z_org.get_original_index_min_max(z_org, z_mrg, iz_m);
          const auto [iy_min, iy_max] = y_org.get_original_index_min_max(y_org, y_mrg, iy_m);
          const auto [ix_min, ix_max] = x_org.get_original_index_min_max(x_org, x_mrg, ix_m);

          for (int iz = iz_min; iz <= iz_max; iz++) {
            for (int iy = iy_min; iy <= iy_max; iy++) {
              for (int ix = ix_min; ix <= ix_max; ix++) {
                auto& vox = result.callVoxel(ix, iy, iz);
                vox.set_density(uniform_prior_density);
                vox.set_tf_exist(true);
              }
            }
          }
          n_prior++;
        }
      }
    }
    LOG_INFO("overlay_merged_density: filled {} reconst voxels with prior density {:.1f}",
             n_prior, uniform_prior_density);
  }

  // --- Pass 2: Overlay reconstruction results ---
  {
    const Grid1d& x_mrg = g3vox_rec.get_x_axis();
    const Grid1d& y_mrg = g3vox_rec.get_y_axis();
    const Grid1d& z_mrg = g3vox_rec.get_z_axis();

    int n_overlaid = 0;
    const int nbinx_m = g3vox_rec.get_nbinx();
    const int nbiny_m = g3vox_rec.get_nbiny();
    const int nbinz_m = g3vox_rec.get_nbinz();

    #pragma omp parallel for collapse(3) reduction(+: n_overlaid)
    for (int iz_m = 0; iz_m < nbinz_m; iz_m++) {
      for (int iy_m = 0; iy_m < nbiny_m; iy_m++) {
        for (int ix_m = 0; ix_m < nbinx_m; ix_m++) {
          const auto& vox_m = g3vox_rec.getVoxel(ix_m, iy_m, iz_m);
          if (!vox_m.get_tf_exist()) continue;

          const double density = vox_m.get_density();
          const auto [iz_min, iz_max] = z_org.get_original_index_min_max(z_org, z_mrg, iz_m);
          const auto [iy_min, iy_max] = y_org.get_original_index_min_max(y_org, y_mrg, iy_m);
          const auto [ix_min, ix_max] = x_org.get_original_index_min_max(x_org, x_mrg, ix_m);

          for (int iz = iz_min; iz <= iz_max; iz++) {
            for (int iy = iy_min; iy <= iy_max; iy++) {
              for (int ix = ix_min; ix <= ix_max; ix++) {
                auto& vox = result.callVoxel(ix, iy, iz);
                vox.set_density(density);
                vox.set_tf_exist(true);
              }
            }
          }
          n_overlaid++;
        }
      }
    }
    LOG_INFO("overlay_merged_density: overlaid {} rec voxels onto {} pre-merge grid",
             n_overlaid, result.get_name());
  }

  return result;
}

double Grid3dVoxel::get_highest_exist_z( const int ix, const int iy ) const
{
  fprintf(stderr,"Grid3dVoxel::get_highest_exist_z ix=%d, iy=%d\r",ix,iy);
  // Flags for determining conditions used later
  bool flg_tf_exit_in  = false; // Flag set when entering a region where material exists
  bool flg_tf_exit_out = false; // Flag set when exiting a region where material exists

  const int nbinz = get_nbinz();

  // Traverse each voxel along the z-axis using a for loop.
  for(int iz=0;iz<nbinz;iz++){
    const Voxel& vox = getVoxel(ix,iy,iz);
// PRINT_INT_E(vox.get_tf_exist());
    if( vox.get_tf_exist()==true ){ // If material exists in the voxel,
      flg_tf_exit_in = true; // in flag ON
    }
    if( flg_tf_exit_in==true && vox.get_tf_exist()==false ){ // if in_flag is ON and no material is found when going to the next upper voxel,
      flg_tf_exit_out = true; // out flag ON
      return get_zlow(iz); // return the z value of the bottom surface of that voxel.
    }
  }
  if( flg_tf_exit_in==true && flg_tf_exit_out==false ){ // If out_flag remains OFF,
    return get_zup(nbinz-1); // Return the upper z coordinate of the topmost voxel along the z axis.
  }
  // THROW_ERROR_NAME("Grid3dVoxel::get_highest_voxel_iz : no candidate.");
  // If no material exists in voxels, return zmin.
  const double zmin = get_zmin();
  return zmin;
}

double Grid3dVoxel::get_highest_exist_z( const double x, const double y ) const
{
  const int ix = get_ix(x);
  const int iy = get_iy(y);
  return get_highest_exist_z(ix,iy);
}

double Grid3dVoxel::get_lowest_exist_z( const int ix, const int iy ) const
{
  fprintf(stderr,"Grid3dVoxel::get_lowest_exist_z ix=%d, iy=%d\r",ix,iy);

  const int nbinz = get_nbinz();

  // Traverse from bottom to top, find first voxel with tf_exist==true
  for(int iz=0; iz<nbinz; iz++){
    const Voxel& vox = getVoxel(ix,iy,iz);
    if( vox.get_tf_exist()==true ){
      return get_zlow(iz); // Return the z value of the bottom surface of that voxel
    }
  }
  // If no voxel exists in the column, return zmax
  const double zmax = get_zmax();
  return zmax;
}

double Grid3dVoxel::get_lowest_exist_z( const double x, const double y ) const
{
  const int ix = get_ix(x);
  const int iy = get_iy(y);
  return get_lowest_exist_z(ix,iy);
}

// 2023-01-16 16:14:57
// for debug
void Grid3dVoxel::disp_nearest_column( const double x_in, const double y_in ) const
{
  const int ix = get_ix(x_in);
  const int iy = get_iy(y_in);
  const double x_cnt = 0.5*( get_xup(x_in) + get_xlow(x_in) );
  const double y_cnt = 0.5*( get_yup(y_in) + get_ylow(y_in) );
  const int nbinz = get_nbinz();
  for(int iz=0;iz<nbinz;iz++){
    const Voxel& vox = getVoxel(ix,iy,iz);
    if( vox.get_tf_exist()==false ) continue;
    const double z_low = get_z_axis().get_lower_value(iz);
    const double z_upp = get_z_axis().get_upper_value(iz);
    LOG_INFO("");
    LOG_INFO("ix={}, iy={}, iz={}",ix,iy,iz);
    LOG_INFO("x_cnt={:.1f}, y_cnt={:.1f}, z_low={:.1f}, z_upp={:.1f}"
    ,x_cnt, y_cnt, z_low, z_upp);
  }
}


void Grid3dVoxel::out_xy_sqrt_matcovdens_by_z(
  const std::string &prefix,
  const Uqiv uqiv_base,
  const Eigen::MatrixXf &mat_cov_dens,
  const double zmin, const double zmax, const double zstep) const
{
  // Number of z samples to generate.
  const int nz = (int)((zmax - zmin) / zstep);
  LOG_DEBUG("nz = {}", nz);

  // Prepare container to store (x, y, value) for each z index value
  std::map<int, std::vector<std::array<double,3>>> iz_to_data;

  // Extract target vector from column
  const Eigen::MatrixXf vecxf_cov_dens = mat_cov_dens.col(uqiv_base);
  const Uqiv n_uqiv = vecxf_cov_dens.size();
  LOG_DEBUG("n_uqiv = {}", n_uqiv);

  // #pragma omp parallel for
  for (Uqiv irow = 0; irow < n_uqiv; irow++) {
    const auto [ix, iy, iz] = get_ixiyiz(irow);
    const double x = get_xcnt(ix);
    const double y = get_ycnt(iy);
    const double z = get_zcnt(iz);
    const double value = sqrt(vecxf_cov_dens(irow));

    // Check if z is within specified range
    if (z < zmin || z >= zmax) continue;

    // Calculate index corresponding to z
    const int iz_idx = (int)((z - zmin) / zstep);
    iz_to_data[iz_idx].push_back({x, y, value});
  }

  // Output to file for each z index
  for (int iz_idx = 0; iz_idx < nz; iz_idx++) {
    // Calculate actual z value
    double z_val = zmin + ((double)iz_idx + 0.5) * zstep;
    // Create filename (e.g., prefix_z_<z_value>.txt)
    char filename[256];
    std::snprintf(filename, sizeof(filename), "%s_uqivbase%d_z%.0lf.tmp", prefix.c_str(), uqiv_base, z_val);

    FILE *fout = std::fopen(filename, "wt");
    if (!fout) {
      LOG_WARN("Failed to open file: {}", filename);
      continue; // Simply skip if file open fails
    }
    LOG_INFO("Writing to file: {}", filename);

    if (iz_to_data.count(iz_idx) > 0) {
      for (const auto &data : iz_to_data[iz_idx]) {
        // data[0]: x, data[1]: y, data[2]: value
        std::fprintf(fout, "%E %E %E\n", data[0], data[1], data[2]);
      }
    }
    std::fclose(fout);
  }
}

Eigen::MatrixXf Grid3dVoxel::mp_make_voxel_distance_matrix() const
{
  const int n_voxel_exist = get_n_vox_exist();
  Eigen::MatrixXf mat_voxel_distance = Eigen::MatrixXf::Zero(n_voxel_exist,n_voxel_exist);
  LOG_INFO("memory allocation of mat_voxel_distance({},{}) {} elements OK"
  ,n_voxel_exist,n_voxel_exist,n_voxel_exist*n_voxel_exist);
  
  int ivox2; // unique index of tf_exist=true;
  // #pragma omp parallel for private(tp1_ixiyiz,tp2_ixiyiz,ivox2)
  #pragma omp parallel for private(ivox2)
  for(int ivox1=0;ivox1<n_voxel_exist-1;ivox1++){
    const Ixiyiz ixiyiz1 = get_ixiyiz(ivox1);
    for(ivox2=ivox1+1;ivox2<n_voxel_exist;ivox2++){
      const Ixiyiz ixiyiz2 = get_ixiyiz(ivox2);
      mat_voxel_distance(ivox1,ivox2) = calc_dist_dxyz(ixiyiz1,ixiyiz2);
    }
  }
  // Include ivox1>ivox2 components that were not calculated.
  mat_voxel_distance += mat_voxel_distance.transpose();
  return mat_voxel_distance;
}

// 2023-03-17 17:35:25
Eigen::VectorXf Grid3dVoxel::get_vecxf_density() const
{
  // loop of uqiv from uqiv_min to uqiv_max using get_ump_uqiv_ixiyiz_ref
  const Uqiv uqiv_min = get_uqiv_min();
  const Uqiv uqiv_max = get_uqiv_max();
  const int n_vox_exist = uqiv_max - uqiv_min + 1;
  Eigen::VectorXf vecxf_dens = Eigen::VectorXf::Zero(n_vox_exist);
  for(Uqiv uqiv=uqiv_min;uqiv<=uqiv_max;uqiv++){
    const Ixiyiz tp_ixiyiz = get_ixiyiz(uqiv);
    const Voxel& vox = getVoxel(tp_ixiyiz);
    vecxf_dens(uqiv) = vox.get_density();
  }
  return vecxf_dens;
}

void Grid3dVoxel::set_density( const Eigen::VectorXf &vecxf_dens_in )
{
  LOG_WARN("is not perfect , should be modified 2024-11-20 12:21:08");
  // loop of uqiv from uqiv_min to uqiv_max using get_ump_uqiv_ixiyiz_ref
  const Uqiv uqiv_min = get_uqiv_min();
  const Uqiv uqiv_max = get_uqiv_max();
  const int n_uqiv_vox = uqiv_max - uqiv_min + 1;
  LOG_DEBUG("uqiv_min={}, uqiv_max={}, n_uqiv_vox={}"
    , uqiv_min, uqiv_max, n_uqiv_vox);

  // check size of vecxf_dens_in and n_vox_exist
  if( vecxf_dens_in.size() != n_uqiv_vox ){
    THROW_ERROR_NAME3(
      "Grid3dVoxel::set_density, size of vecxf_dens_in != n_uqiv_vox"
      , vecxf_dens_in.size(), n_uqiv_vox);
  }

  // assign density to each voxel
  for(Uqiv uqiv=uqiv_min;uqiv<=uqiv_max;uqiv++){
    const Ixiyiz tp_ixiyiz = get_ixiyiz(uqiv);
    Voxel& vox = callVoxel(tp_ixiyiz);
    const double density = vecxf_dens_in(uqiv);
    vox.set_density(density);
  }
}

// Reassign unique indices based on n_hit thresholds for voxels.
// Returns a map from old unique indices to new indices for later cleanup.
std::map<Grid3d::Uqiv,Grid3d::Uqiv> Grid3dVoxel::re_assign_uqiv_by_nhit_det(
    const bool tf_exist_in, const Grid3dVoxel::Parameters &prm_g3vox
  , const Uqiv uqiv_start )
{
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  const int nbinz = get_nbinz();

  // map old uqiv is available or not
  // std::map<int,bool> map_uqiv_available;

  // store available old uqiv and new uqiv before and after re-assign
  std::map<int,int> map_uqiv_old_new_avail;

  const auto bimap_before = get_bimap_uqiv_ixiyiz_copy();

  // Save old uqiv_min before clearing — needed for hit-data vector lookup.
  const Uqiv uqiv_min_old = get_uqiv_min();

  // Initialize uqiv_umps before assigning new indices.
  clear_uqiv_umps();

  // set initial value of uqiv
  Uqiv uqiv_new = uqiv_start;
  set_uqiv_min( uqiv_new );

  // loop for voxels
  int n_vox_count = 0;
  // start to assign new unique index
  // loop for voxels
  for(int iz=0; iz < nbinz; iz++){
    for(int iy=0; iy < nbiny; iy++){
      for(int ix=0; ix < nbinx; ix++){
        const Voxel& vox = getVoxel(ix,iy,iz);
        if( vox.get_tf_exist()!=tf_exist_in ) continue;

        // Look up old uqiv for hit-data access
        const Ixiyiz ixiyiz_key = {ix,iy,iz};
        const auto it_old = bimap_before.getMapBA().find(ixiyiz_key);
        if( it_old == bimap_before.getMapBA().end() ) continue;
        const Uqiv uqiv_old = it_old->second;
        const int nhd = __builtin_popcountll(vec_hit_det_.at(uqiv_old - uqiv_min_old));
        const Nhit nhe = vec_n_hit_ele_.at(uqiv_old - uqiv_min_old);

        if( nhd < prm_g3vox.n_hit_det_min ) continue;
        if( nhd > prm_g3vox.n_hit_det_max ) continue;
        if( nhe < prm_g3vox.n_hit_ele_min ) continue;
        if( nhe > prm_g3vox.n_hit_ele_max ) continue;

        // make tp_ixiyiz
        const Ixiyiz ixiyiz = {ix,iy,iz};

        // insert tp_ixiyiz to ump_uqiv_ixiyiz
        insert_to_uqiv_umps( uqiv_new, ixiyiz );

        // insert old_uqiv and new uqiv to map_uqiv_old_new_avail
        map_uqiv_old_new_avail.insert(std::pair(uqiv_old,uqiv_new));

        uqiv_new++;
        n_vox_count++;
      }
    }
  } // voxel loop end

  if( bimap_before != get_bimap_uqiv_ixiyiz_ref() ){
    LOG_INFO("bimap_uqiv_ixiyiz changed after re-assignment");
  }else{
    LOG_INFO("bimap_uqiv_ixiyiz unchanged after re-assignment");
  }

  set_uqiv_max( uqiv_new - 1 );

  set_n_vox_exist(n_vox_count);
  LOG_INFO(", n_vox_exist = {}",n_vox_count);

  // rebuild flat array after reassigning uqiv
  build_flat_uqiv();

  // re-allocate hit vectors for the new uqiv range
  allocate_hit_data();

  return map_uqiv_old_new_avail;
}

// set n_hit_ele for voxel(ix,iy,iz) tf_exist is true (grid-level storage)
void Grid3dVoxel::set_n_hit_ele(const int ix ,const int iy ,const int iz
  , const Nhit n_hit_ele_in )
{
  const Voxel& vox = getVoxel(ix,iy,iz);
  if(!vox.get_tf_exist()) return;
  const Uqiv uqiv = get_uqiv_fast(ix,iy,iz);
  if(uqiv == UqivNotFound) return;
  const auto uqiv_min = get_uqiv_min();
  vec_n_hit_ele_.at(uqiv - uqiv_min) = n_hit_ele_in;
}

// set n_hit_ele for all voxel tf_exist is true (grid-level storage)
void Grid3dVoxel::set_n_hit_ele( const Nhit n_hit_ele_in )
{
  std::fill(vec_n_hit_ele_.begin(), vec_n_hit_ele_.end(), n_hit_ele_in);
}

// zmax_voxel_with_tf_exist_true in g3vox become zmin in g2pil
// Converted from g2pil to g3vox, and then
// create a cavity in g2pil from merged g3vox_input (this instance).
// Function that returns the outer shell part as g2pil.
Grid2dPillar Grid3dVoxel::make_void_g2pil( const Grid2dPillar &g2pil_org ) const
{
  if( this->get_x_interval() < g2pil_org.get_x_interval() ){
    THROW_ERROR_NAME("this->get_x_interval() < g2pil.get_x_interval()");
  }
  if( this->get_y_interval() < g2pil_org.get_y_interval() ){
    THROW_ERROR_NAME("this->get_y_interval() < g2pil.get_y_interval()");
  }

  Grid2dPillar g2pil_void(g2pil_org);

  // Get number of bins in g3vox
  const int nbinx = this->get_nbinx();
  const int nbiny = this->get_nbiny();
  // ix and iy are the index of g3vox

  // Get g2pil grid size
  const double dx_cub = g2pil_org.get_x_interval();
  const double dy_cub = g2pil_org.get_y_interval();

  // g3vox(Assuming cavity)loop over y_grid of
  for(int iy=0;iy<nbiny;iy++){
    // Get y range
    const double y_lower = this->get_ylow(iy);
    const double y_upper = this->get_yup(iy);
    
    // g3vox(Assuming cavity)loop over x_grid of
    for(int ix=0;ix<nbinx;ix++){
      // Get x range
      const double x_lower = this->get_xlow(ix);
      const double x_upper = this->get_xup(ix);
      
      // Get maximum height of internal cavity,
      // return this->get_zmin() if no cavity
      double zmax_vox = this->get_highest_exist_z(ix,iy);
      
      // Go to next Grid if no cavity exists
      if( zmax_vox == this->get_zmin() ) continue;

      // Further loop inside g3vox
      for(double y=y_lower+0.5*dy_cub; y<y_upper; y+=dy_cub){
        for(double x=x_lower+0.5*dx_cub; x<x_upper; x+=dx_cub){
          Pillar& cub = g2pil_void.callPillar(x,y);
          double zmax_cub = cub.get_zmax();
          if( zmax_vox > zmax_cub ){
            THROW_ERROR_NAME("zmax_vox > zmax_cub"); // Does this error mean an error in this->merge?
          }
          cub.set_zmin(zmax_vox);
          // double zmin_cub = cub.get_zmin();
        }
      }
    } // ix loop in g3vox
  } // iy loop in g3vox
  return g2pil_void;
}

// Get Eigen::Vector3d v3_AABB_min and v3_AABB_max
AABB3d Grid3dVoxel::get_AABB3d(const int ix, const int iy, const int iz) const
{
  // Get coordinates of faces along x and y axes
  const double xmin_AABB = this->get_xlow(ix);
  const double xmax_AABB = this->get_xup(ix);
  const double ymin_AABB = this->get_ylow(iy);
  const double ymax_AABB = this->get_yup(iy);
  const double zmin_AABB = this->get_zlow(iz);
  const double zmax_AABB = this->get_zup(iz);

  // Create AABB for rectangular solid.
  const Eigen::Vector3d v3_AABB_min(xmin_AABB,ymin_AABB,zmin_AABB);
  const Eigen::Vector3d v3_AABB_max(xmax_AABB,ymax_AABB,zmax_AABB);
  const AABB3d aabb3d(v3_AABB_min,v3_AABB_max);

  return aabb3d;
}


// density Get Eigen::Vector3d v3_AABB_min and v3_AABB_max
std::tuple< double, AABB3d >
  Grid3dVoxel::get_density_AABB3d(const int ix, const int iy, const int iz) const
{
  // get voxel information
  const Voxel& vox = this->getVoxel(ix,iy,iz);
  const double density = vox.get_density();
  // return density and AABB3d
  return std::make_tuple(density,get_AABB3d(ix,iy,iz));
}

// add CheckerBoard3dParameters structure
void Grid3dVoxel::add_density_structure( const CheckerBoard3dParameters &prm )
{
  // NaN sources check
  if (!std::isfinite(prm.delta_density)) {
    THROW_ERROR("Grid3dVoxel::add_density_structure(Checkerboard): delta_density is non-finite. name={}, delta_density={}", prm.name, prm.delta_density);
  }
  if (!std::isfinite(prm.delta_density_offset)) {
    THROW_ERROR("Grid3dVoxel::add_density_structure(Checkerboard): delta_density_offset is non-finite. name={}, delta_density_offset={}", prm.name, prm.delta_density_offset);
  }

  if( prm.tf_exec==false ){
    LOG_INFO(", name={}, tf_exec=false",prm.name);
  }
  else{
    // Get number of bins
    const int nbinx = get_nbinx();
    const int nbiny = get_nbiny();
    const int nbinz = get_nbinz();

    // Integer-index-based checkerboard sign determination
    // (avoids floating-point error in floor() with non-integer grid intervals)
    const int ix_cnt = get_ix(prm.v3_pos_cnt.x());
    const int iy_cnt = get_iy(prm.v3_pos_cnt.y());
    const int iz_cnt = get_iz(prm.v3_pos_cnt.z());
    const int mx = prm.v3_len_interval_mult.x();
    const int my = prm.v3_len_interval_mult.y();
    const int mz = prm.v3_len_interval_mult.z();

    // Pre-compute cylinder check parameters (avoid string compare in hot loop)
    const bool tf_cylinder = (prm.region_type == "cylinder");
    const double radius_x = prm.radius_x_meters;
    const double radius_y = prm.radius_y_meters;
    const double inv_rx_sq = (radius_x > 0.0) ? 1.0 / (radius_x * radius_x) : 0.0;
    const double inv_ry_sq = (radius_y > 0.0) ? 1.0 / (radius_y * radius_y) : 0.0;
    const double cyl_cx = prm.v3_pos_cnt.x();
    const double cyl_cy = prm.v3_pos_cnt.y();

    // #pragma omp parallel for
    for(int iz=0;iz<nbinz;iz++){
      double zcnt = get_zcnt(iz);
      if( zcnt <  prm.aabb3d.zmin() ) continue;
      if( zcnt >= prm.aabb3d.zmax() ) continue;
      for(int iy=0;iy<nbiny;iy++){
        double ycnt = get_ycnt(iy);
        if( ycnt <  prm.aabb3d.ymin() ) continue;
        if( ycnt >= prm.aabb3d.ymax() ) continue;
        for(int ix=0;ix<nbinx;ix++){
          double xcnt = get_xcnt(ix);
          if( xcnt <  prm.aabb3d.xmin() ) continue;
          if( xcnt >= prm.aabb3d.xmax() ) continue;

          // Elliptic cylinder region: reject voxels outside XY ellipse
          if (tf_cylinder) {
            const double dx = xcnt - cyl_cx;
            const double dy = ycnt - cyl_cy;
            if (dx * dx * inv_rx_sq + dy * dy * inv_ry_sq > 1.0) continue;
          }

          // access to voxel
          Voxel& vox = callVoxel(ix,iy,iz);

          // If no material, go to next voxel
          if( vox.get_tf_exist()==false ) continue;

          // Checkerboard sign via integer block indices
          const int bx = floor_div(ix - ix_cnt, mx);
          const int by = floor_div(iy - iy_cnt, my);
          const int bz = floor_div(iz - iz_cnt, mz);
          // Checkerboard sign: flip when the sum of block indices is odd.
          const bool odd_block = (bx + by + bz) % 2 != 0;
          const double sign = odd_block ? -1.0 : 1.0;
          
          const double delta_density = prm.delta_density_offset + sign * prm.delta_density;

          // Add delta_density_offset + or - delta_density.
          vox.add_density(delta_density);
        } // ix loop end
      } // iy loop end
    } // iz loop end
  } // if( prm.tf_exec==true) end
}

// add Ellipsoid structure
void Grid3dVoxel::add_density_structure(
  const EllipsoidParameters &prm )
{
  // NaN sources check
  if (!std::isfinite(prm.delta_density)) {
    THROW_ERROR("Grid3dVoxel::add_density_structure(Ellipsoid): delta_density is non-finite. name={}, delta_density={}", prm.name, prm.delta_density);
  }

  if( prm.tf_exec==false ){
    LOG_INFO(", name={}, tf_exec=false",prm.name);
  }
  else{
    // Get number of bins
    const int nbinx = get_nbinx();
    const int nbiny = get_nbiny();
    const int nbinz = get_nbinz();

    // Check if all bins are inside the ellipsoid.
    #pragma omp parallel for
    for(int iz=0;iz<nbinz;iz++){
      double zcnt = get_zcnt(iz);
      for(int iy=0;iy<nbiny;iy++){
        double ycnt = get_ycnt(iy);
        for(int ix=0;ix<nbinx;ix++){
          double xcnt = get_xcnt(ix);

          // access to voxel
          Voxel& vox = callVoxel(ix,iy,iz);

          // If no material, go to next voxel
          if( vox.get_tf_exist()==false ) continue;

          // get 2d position of grid center
          const Eigen::Vector3d v3_pos(xcnt,ycnt,zcnt);

          // judge whether the position is inside the ellipsoid or not
          const bool is_inside
           = geom_util::isInsideEllipsoid(
              prm.v3_pos_cnt, prm.v3_length
            , prm.theta_x, prm.theta_y, prm.theta_z
            , prm.rotation_type, v3_pos);

          // if not, continue
          if(!is_inside ) continue;

          // if inside, add delta_density
          vox.add_density(prm.delta_density);
        } // ix loop end
      } // iy loop end
    } // iz loop end
  } // if( prm.tf_exec==true ) end
}

// add Cylinder structure
void Grid3dVoxel::add_density_structure(
  const CylinderParameters &prm )
{
  // NaN sources check
  if (!std::isfinite(prm.delta_density)) {
    THROW_ERROR("Grid3dVoxel::add_density_structure(Cylinder): delta_density is non-finite. name={}, delta_density={}", prm.name, prm.delta_density);
  }

  if( prm.tf_exec==false ){
    LOG_INFO(", name={}, tf_exec=false",prm.name);
  }
  else{
    // Get number of bins
    const int nbinx = get_nbinx();
    const int nbiny = get_nbiny();
    const int nbinz = get_nbinz();

    // Check if all bins are inside the ellipsoid.
    #pragma omp parallel for
    for(int iz=0;iz<nbinz;iz++){
      double zcnt = get_zcnt(iz);
      for(int iy=0;iy<nbiny;iy++){
        double ycnt = get_ycnt(iy);
        for(int ix=0;ix<nbinx;ix++){
          double xcnt = get_xcnt(ix);

          // access to voxel
          Voxel& vox = callVoxel(ix,iy,iz);

          // If no material, go to next voxel
          if( vox.get_tf_exist()==false ) continue;

          // get 2d position of grid center
          const Eigen::Vector3d v3_pos(xcnt,ycnt,zcnt);

          // judge whether the position is inside the ellipsoid or not
          const bool is_inside
           = geom_util::isInsideCylinder(
              prm.v3_pos_cnt, prm.v3_length
            , prm.theta_x, prm.theta_y, prm.theta_z
            , prm.rotation_type, v3_pos);

          // if not, continue
          if(!is_inside ) continue;

          // if inside, add delta_density
          vox.add_density(prm.delta_density);
        } // ix loop end
      } // iy loop end
    } // iz loop end
  } // if( prm.tf_exec==true ) end
}

void Grid3dVoxel::add_density_structure(
  const CuboidParameters &prm )
{
  // NaN sources check
  if (!std::isfinite(prm.delta_density)) {
    THROW_ERROR("Grid3dVoxel::add_density_structure(Cuboid): delta_density is non-finite. name={}, delta_density={}", prm.name, prm.delta_density);
  }

  if( prm.tf_exec==false ){
    LOG_INFO(", name={}, tf_exec=false",prm.name);
  }
  else{
    // Get number of bins
    const int nbinx = get_nbinx();
    const int nbiny = get_nbiny();
    const int nbinz = get_nbinz();

    // Check if each voxel center is inside the cuboid.
    #pragma omp parallel for
    for(int iz=0;iz<nbinz;iz++){
      double zcnt = get_zcnt(iz);
      for(int iy=0;iy<nbiny;iy++){
        double ycnt = get_ycnt(iy);
        for(int ix=0;ix<nbinx;ix++){
          double xcnt = get_xcnt(ix);

          // access to voxel
          Voxel& vox = callVoxel(ix,iy,iz);

          // If no material, go to next voxel
          if( vox.get_tf_exist()==false ) continue;

          // get 3d position of grid center
          const Eigen::Vector3d v3_pos(xcnt,ycnt,zcnt);

          // judge whether the position is inside the cuboid or not
          const bool is_inside
           = geom_util::isInsideCuboid(
              prm.v3_pos_cnt, prm.v3_length
            , prm.theta_x, prm.theta_y, prm.theta_z
            , prm.rotation_type, v3_pos);

          // if not, continue
          if(!is_inside ) continue;

          // if inside, add delta_density
          vox.add_density(prm.delta_density);
        } // ix loop end
      } // iy loop end
    } // iz loop end
  } // if( prm.tf_exec==true ) end
}

// Output x, y, density data at elevation z_in.
void Grid3dVoxel::out_cross_section_z(
  const fs::path& pathout, const double z_in
  , const bool tf_only_exist, const double dens_false
  , const Nhit n_det_thres ) const
{
  LOG_INFO("output to pathout={}",pathout.string());

  // open file
  FILE *fout = myapp::get_fout(pathout);

  // Get number of bins
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  // const int nbinz = get_nbinz();
  
  // check z_in is in the range of z-axis
  if( z_in <  get_zmin() ){
    LOG_WARN("z_in <  get_zmin(). skip.");
    return ;
  }
  if( z_in >= get_zmax() ){
    LOG_WARN("z_in >= get_zmax(). skip.");
    return ;
  }
  
  // get the index of z_in in z-axis
  const int iz = get_iz(z_in);

  for(int iy=0;iy<nbiny;iy++){
    // const double ylow = get_ylow(iy);
    const double ycnt = get_ycnt(iy);
    for(int ix=0;ix<nbinx;ix++){
      // const double xlow = get_xlow(ix);
      const double xcnt = get_xcnt(ix);
      const Voxel& vox = getVoxel(ix,iy,iz);
      double density = vox.get_density();

      if( vox.get_tf_exist()==false ){
        // if tf_only_exist is true, and voxel does not exist, skip
        if(tf_only_exist==true){ continue;}
        
        // if tf_only_exist is false, and voxel does not exist
        // , set density to dens_false
        else{ density = dens_false; }
      }
      // if n_hit_det < n_det_thres, skip (grid-level hit data)
      const Uqiv uqiv_cs = get_uqiv_fast(ix,iy,iz);
      // nhd_cs: n_hit_det for cross-section output (cs = cross section)
      const int nhd_cs = (uqiv_cs != UqivNotFound) ? get_n_hit_det_grid(uqiv_cs) : 0;
      if( nhd_cs < n_det_thres ) continue;

      fprintf(fout,"%E %E %E\n",xcnt,ycnt,density);
    }
  }
  myapp::close(fout,pathout);
}

// Output x, y, 0(tf_exist==false) or 1(tf_exist==true) data at elevation z_in.
void Grid3dVoxel::out_cross_section_z_mask(
  const fs::path& pathout, const double z_in) const
{
  LOG_INFO("output to pathout={}",pathout.string());

  // open file
  FILE *fout = myapp::get_fout(pathout);

  // Get number of bins
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  // const int nbinz = get_nbinz();
  
  // check z_in is in the range of z-axis
  if( z_in <  get_zmin() ){
    LOG_WARN("z_in <  get_zmin(). skip.");
    return ;
  }
  if( z_in >= get_zmax() ){
    LOG_WARN("z_in >= get_zmax(). skip.");
    return ;
  }
  
  // get the index of z_in in z-axis
  const int iz = get_iz(z_in);

  for(int iy=0;iy<nbiny;iy++){
    // const double ylow = get_ylow(iy);
    const double ycnt = get_ycnt(iy);
    for(int ix=0;ix<nbinx;ix++){
      // const double xlow = get_xlow(ix);
      const double xcnt = get_xcnt(ix);
      const Voxel& vox = getVoxel(ix,iy,iz);
      if( vox.get_tf_exist()==false ){
        fprintf(fout,"%E %E %d\n",xcnt,ycnt,0);
      }else{
        fprintf(fout,"%E %E %d\n",xcnt,ycnt,1);
      }
    }
  }
  myapp::close(fout,pathout);
}

/// @brief Output header information for out_cross_section_z.
void Grid3dVoxel::out_cross_section_z_header( FILE *fout
  , const CrossSectionZParameters& prm_zcross ) const
{
  // Output z information
  fprintf(fout,"# z_info %12.9E %12.9E %12.9E\n"
    ,prm_zcross.zmin,prm_zcross.zmax,prm_zcross.zstep);
  // Output x, y information
  fprintf(fout,"# x_info %12.9E %12.9E %12.9E\n"
    ,prm_zcross.xmin,prm_zcross.xmax,prm_zcross.xstep);
  fprintf(fout,"# y_info %12.9E %12.9E %12.9E\n"
    ,prm_zcross.ymin,prm_zcross.ymax,prm_zcross.ystep);
  // Output detector information
  fprintf(fout,"# n_detector %d\n",prm_zcross.n_detector);
  // Number of z values
  const int nz = (prm_zcross.zmax-prm_zcross.zmin)/prm_zcross.zstep+1;
  fprintf(fout,"# nz %d\n",nz);
}

 // Output cross-sections from z_cross_min to z_cross_max with z_step increment.
void Grid3dVoxel::out_cross_section_z_all(
  const fs::path& pathout, const CrossSectionZParameters &prm_zcross) const
{
  if (prm_zcross.output_binary) {
    // Binary content gets a .tmpbin extension so it is not mistaken for text.
    fs::path pathout_bin = pathout;
    pathout_bin.replace_extension(".tmpbin");
    out_cross_section_z_all_binary(pathout_bin, prm_zcross);
    return;
  }
  LOG_INFO("pathout={} ... ",pathout.string());
  // open file
  FILE *fout = myapp::get_fout(pathout);
  // write header
  out_cross_section_z_header(fout,prm_zcross);

  // write information of each z
  const double zmin = prm_zcross.zmin;
  const double zmax = prm_zcross.zmax;
  const double zstep = prm_zcross.zstep;
  const double g3vox_zmin = get_zmin();
  const double g3vox_zmax = get_zmax();
  const double g3vox_z_interval = get_z_interval();

  // Track processed iz to avoid duplicate output when zstep < z_interval
  std::set<int> processed_iz;

  for(double z=zmin;z<=zmax;z+=zstep){
    // Skip z values outside the g3vox z-axis range
    if(z < g3vox_zmin){
      LOG_WARN_ND("z={} < g3vox_zmin={}, skipping", z, g3vox_zmin);
      continue;
    }
    if(z >= g3vox_zmax){
      LOG_WARN_ND("z={} >= g3vox_zmax={}, skipping", z, g3vox_zmax);
      continue;
    }
    const int iz = get_iz(z);
    if(processed_iz.count(iz) > 0){
      LOG_WARN_ND("iz={} already processed (zstep={} < z_interval={}), skipping z={}",
               iz, zstep, g3vox_z_interval, z);
      continue;
    }
    processed_iz.insert(iz);
    Grid2dVoxel g2vox = get_Grid2dVoxel_z(iz);
    // Override z range so that zcnt = 0.5*(zmin+zmax) equals the query z,
    // ensuring consistent z labels across grids with different pitches.
    const double half_z = 0.5 * g3vox_z_interval;
    g2vox.set_zmin(z - half_z);
    g2vox.set_zmax(z + half_z);
    g2vox.out_voxel_all(fout);
    LOG_INFO("z={} (iz={}) done", z, iz);
  }
}

void Grid3dVoxel::out_cross_section_z_all_binary(
  const fs::path& pathout, const CrossSectionZParameters &prm_zcross) const
{
  LOG_INFO("pathout={} (binary) ... ", pathout.string());
  FILE *fout = myapp::get_fout_binary(pathout);

  auto write_bytes = [fout](const void* ptr, size_t size) {
    if (std::fwrite(ptr, size, 1, fout) != 1) {
      THROW_ERROR("Grid3dVoxel::out_cross_section_z_all_binary: failed to write data");
    }
  };
  auto write_value = [&write_bytes](const auto& value) {
    write_bytes(&value, sizeof(value));
  };

  // Magic
  const char magic[8] = {'G','2','Z','B','I','N','\0','\0'};
  write_bytes(magic, sizeof(magic));

  // Version and reserved
  const std::uint32_t version = 1;
  write_value(version);
  const std::uint32_t reserved = 0;
  write_value(reserved);

  // Header doubles: x/y/z min, max, step
  const double header_doubles[] = {
    prm_zcross.xmin, prm_zcross.xmax, prm_zcross.xstep,
    prm_zcross.ymin, prm_zcross.ymax, prm_zcross.ystep,
    prm_zcross.zmin, prm_zcross.zmax, prm_zcross.zstep
  };
  for (const double value : header_doubles) {
    write_value(value);
  }

  // Collect valid z-values (same duplicate-iz skip logic as text version)
  const double g3vox_zmin = get_zmin();
  const double g3vox_zmax = get_zmax();
  const double g3vox_z_interval = get_z_interval();
  std::set<int> processed_iz;
  std::vector<std::pair<double, int>> z_iz_pairs;

  for (double z = prm_zcross.zmin; z <= prm_zcross.zmax; z += prm_zcross.zstep) {
    if (z < g3vox_zmin || z >= g3vox_zmax) continue;
    const int iz = get_iz(z);
    if (processed_iz.count(iz) > 0) continue;
    processed_iz.insert(iz);
    z_iz_pairs.emplace_back(z, iz);
  }

  const std::int32_t nbinx = get_nbinx();
  const std::int32_t nbiny = get_nbiny();
  const std::int32_t nz = static_cast<std::int32_t>(z_iz_pairs.size());
  const std::int32_t header_ints[] = {
    static_cast<std::int32_t>(prm_zcross.n_detector), nbinx, nbiny, nz
  };
  for (const std::int32_t value : header_ints) {
    write_value(value);
  }

  const std::int64_t n_records =
    static_cast<std::int64_t>(nbinx) *
    static_cast<std::int64_t>(nbiny) *
    static_cast<std::int64_t>(nz);
  write_value(n_records);

  // Write voxel data for each z-level
  for (const auto& [z, iz] : z_iz_pairs) {
    Grid2dVoxel g2vox = get_Grid2dVoxel_z(iz);
    const double half_z = 0.5 * g3vox_z_interval;
    g2vox.set_zmin(z - half_z);
    g2vox.set_zmax(z + half_z);
    g2vox.out_voxel_all_binary(fout);
    LOG_INFO("z={} (iz={}) done", z, iz);
  }

  myapp::close(fout, pathout);
}

// Output x, y, density data at elevation z_in.
void Grid3dVoxel::out_cross_section_z_7clm(
  const fs::path& pathout, const double z_in) const
{
  LOG_INFO("output to pathout={}",pathout.string());

  // open file
  FILE *fout = myapp::get_fout(pathout);

  // Get number of bins
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  // const int nbinz = get_nbinz();
  
  // check z_in is in the range of z-axis
  if( z_in <  get_zmin() ){
    LOG_WARN("z_in <  get_zmin(). skip.");
    return ;
  }
  if( z_in >= get_zmax() ){
    LOG_WARN("z_in >= get_zmax(). skip.");
    return ;
  }
  
  // get the index of z_in in z-axis
  const int iz = get_iz(z_in);

  // z value
  const double zlow = get_zlow(iz);
  const double zupp = get_zup(iz);

  // loop over ix, iy
  for(int iy=0;iy<nbiny;iy++){
    // const double ycnt = get_ycnt(iy);
    const double ylow = get_ylow(iy);
    const double yupp = get_yup(iy);
    for(int ix=0;ix<nbinx;ix++){
      // const double xcnt = get_xcnt(ix);
      const double xlow = get_xlow(ix);
      const double xupp = get_xup(ix);
      // const Voxel& vox = getVoxel(ix,iy,iz);  // error: call to member function 'get_Voxel' is ambiguous
      const Voxel& vox = getVoxel( static_cast<int>(ix),static_cast<int>(iy),static_cast<int>(iz) );
      const double density = vox.get_density();
      // fprintf(fout,"%E %E %E\n",xcnt,ycnt,density);
      // fprintf(fout,"%E %E %E\n",xlow,ylow,density);
      if(vox.get_tf_exist()==false) continue;
      fprintf(fout,"%10.1lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf %E\n"
      ,xlow,xupp,ylow,yupp,zlow,zupp,density);
      // fprintf(fout,"%10.1lf %10.1lf %E\n",xcnt,ycnt,density);
    }
  }
  myapp::close(fout,pathout);
}

// out_cross_section_z_10clm
// Output x, y, density data at elevation z_in.
void Grid3dVoxel::out_cross_section_z_10clm(
  const fs::path& pathout, const double z_in) const
{
  LOG_INFO("output to pathout={}",pathout.string());
  // open file
  FILE *fout = myapp::get_fout(pathout);

  // Get number of bins
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  // const int nbinz = get_nbinz();
  
  // check z_in is in the range of z-axis
  if( z_in <  get_zmin() ){
    LOG_WARN("z_in <  get_zmin(). skip.");
    return ;
  }
  if( z_in >= get_zmax() ){
    LOG_WARN("z_in >= get_zmax(). skip.");
    return ;
  }
  
  // get the index of z_in in z-axis
  const int iz = get_iz(z_in);

  // z value
  const double zlow = get_zlow(iz);
  const double zupp = get_zup(iz);

  // loop over ix, iy
  for(int iy=0;iy<nbiny;iy++){
    // const double ycnt = get_ycnt(iy);
    const double ylow = get_ylow(iy);
    const double yupp = get_yup(iy);
    for(int ix=0;ix<nbinx;ix++){
      // const double xcnt = get_xcnt(ix);
      const double xlow = get_xlow(ix);
      const double xupp = get_xup(ix);
      // const Voxel& vox = getVoxel(ix,iy,iz);  // error: call to member function 'get_Voxel' is ambiguous
      const Voxel& vox = getVoxel( static_cast<int>(ix),static_cast<int>(iy),static_cast<int>(iz) );
      const double density = vox.get_density();
      // fprintf(fout,"%E %E %E\n",xcnt,ycnt,density);
      // fprintf(fout,"%E %E %E\n",xlow,ylow,density);
      if(vox.get_tf_exist()==false) continue;
      fprintf(fout,"%3d %3d %3d %7.1lf %7.1lf %7.1lf %7.1lf %7.1lf %7.1lf %E\n"
      ,ix,iy,iz,xlow,xupp,ylow,yupp,zlow,zupp,density);
      // fprintf(fout,"%10.1lf %10.1lf %E\n",xcnt,ycnt,density);
    }
  }
  myapp::close(fout,pathout);
}

void Grid3dVoxel::apply_shell_density(const std::array<double,4>& density_quad)
{
  // density_quad: [prior, upper, lower, lateral]
  // Skip if all shell values are zero (default-constructed)
  if (density_quad[1] == 0.0 && density_quad[2] == 0.0 && density_quad[3] == 0.0) return;

  const int nbx = get_nbinx();
  const int nby = get_nbiny();
  const int nbz = get_nbinz();
  #pragma omp parallel for collapse(2)
  for (int ix = 0; ix < nbx; ix++) {
    for (int iy = 0; iy < nby; iy++) {
      int iz_min_exist = nbz;
      int iz_max_exist = -1;
      for (int iz = 0; iz < nbz; iz++) {
        if (getVoxel(ix, iy, iz).get_tf_exist()) {
          if (iz < iz_min_exist) iz_min_exist = iz;
          if (iz > iz_max_exist) iz_max_exist = iz;
        }
      }
      for (int iz = 0; iz < nbz; iz++) {
        const auto& vox = getVoxel(ix, iy, iz);
        if (vox.get_tf_exist()) continue;
        if (vox.get_density() == 0.0) continue;
        double shell_dens;
        if (iz_max_exist < 0) {
          shell_dens = density_quad[3]; // lateral: no tf_exist in column
        } else if (iz > iz_max_exist) {
          shell_dens = density_quad[1]; // upper
        } else if (iz < iz_min_exist) {
          shell_dens = density_quad[2]; // lower
        } else {
          shell_dens = density_quad[3]; // lateral: gap between exist voxels
        }
        callVoxel(ix, iy, iz).set_density(shell_dens);
      }
    }
  }
  LOG_INFO("apply_shell_density: applied [upper={}, lower={}, lateral={}]",
           density_quad[1], density_quad[2], density_quad[3]);
}

Grid3dVoxel Grid3dVoxel::write_density_to_cross_section(
  const std::string& prefix
, const Eigen::VectorXf& vec_density
, const CrossSectionZParameters& prm_zcross
, const bool tf_zero_init
, const std::array<double,4>& density_quad ) const
{
  Grid3dVoxel g3vox(*this);
  g3vox.set_name(prefix);
  if (tf_zero_init) {
    g3vox.mp_set_uniform_density_ixiyiz(0.0);
  }
  g3vox.set_density(vec_density);
  g3vox.apply_shell_density(density_quad);

  std::string filename = prefix + "_zcross_all.tmp";
  fs::path pathout = iodir::make_pathout(filename);
  g3vox.out_cross_section_z_all(pathout, prm_zcross);
  return g3vox;
}

// Increase the resolution of Grid3dVoxel.
Grid3dVoxel Grid3dVoxel::get_split_g3vox(
  const int x_factor, const int y_factor, const int z_factor ) const
{

  // split factor should be more than 0
  if( x_factor < 1 ) THROW_ERROR_NAME2("x_factor should be >=1",x_factor);
  if( y_factor < 1 ) THROW_ERROR_NAME2("y_factor should be >=1",y_factor);
  if( z_factor < 1 ) THROW_ERROR_NAME2("z_factor should be >=1",z_factor);

  // make new Grid1d
  const Grid1d x_axis_new = get_x_axis().get_split(x_factor);
  const Grid1d y_axis_new = get_y_axis().get_split(y_factor);
  const Grid1d z_axis_new = get_z_axis().get_split(z_factor);

  // make new Grid3dVoxel. vec_vec_vec_memory_allocate() was done.
  const int n_det = this->get_n_det();
  Grid3dVoxel g3vox_new(x_axis_new,y_axis_new,z_axis_new,n_det);
  g3vox_new.name = this->name + "_splited";

  // transfer the data of vec_vec_vec_Voxel
  const int nbinx_new = g3vox_new.get_nbinx();
  const int nbiny_new = g3vox_new.get_nbiny();
  const int nbinz_new = g3vox_new.get_nbinz();

  // ix_new, iy_new, iz_new is the index of g3vox_new
  for(int iz_new=0;iz_new<nbinz_new;iz_new++){
    const double zlow_new = g3vox_new.get_zlow(iz_new);
    for(int iy_new=0;iy_new<nbiny_new;iy_new++){
      const double ylow_new = g3vox_new.get_ylow(iy_new);
      for(int ix_new=0;ix_new<nbinx_new;ix_new++){
        const double xlow_new = g3vox_new.get_xlow(ix_new);
        
        // get the index of this->vec_vec_vec_Voxel
        const int ix_old = this->get_ix(xlow_new);
        const int iy_old = this->get_iy(ylow_new);
        const int iz_old = this->get_iz(zlow_new);

        // get the pointer of Voxel
        const Voxel vox_old = this->getVoxel(ix_old,iy_old,iz_old);
        Voxel &vox_new = g3vox_new.callVoxel(ix_new,iy_new,iz_new);

        // copy the data (hit data is managed at grid level, not per-voxel)
        vox_new.copy(vox_old);
      }
    }
  }

  return g3vox_new;
}


//===============================
// path length calc function
//===============================

// Primarily used by pathcalc; find the intersection between the ray
// defined by v3_dir/v3_pos and the voxel with indices (ix, iy, iz).
// Returns PATH_NO_HIT when there is no intersection; otherwise returns the
// computed delta_path.
double Grid3dVoxel::get_delta_path(
    const Ixiyiz &ixiyiz_hit, const Ray3d &ray3d ) const
{
  // Get AABB of hit rectangular solid
  AABB3d aabb3d = get_AABB3d(ixiyiz_hit);

  // check whether the detector beam line intersect the voxel or not
  // and return tmin and tmax
  auto [tf_intersect,tmin,tmax] = ray3d.is_intersect(aabb3d);
  
  // if not intersect goto next loop
  if( tf_intersect == false ) return PATH_NO_HIT;

  // if intersect
  if( tmin < 0 ) tmin = 0.0;
  const double delta_path = tmax - tmin;
  
  return delta_path;
}

// If the voxel in g3vox is hit by the ray from g2det,
// increment n_hit_det of the hit voxel.
// Information about whether it was hit is held by mat_path_len.
// void Grid3dVoxel::incr_n_hit_det( const Eigen::MatrixXf &mat_path_len )
// {
//   // Eigen openmp
//   const int n_threads = omp_get_max_threads();
//   // myapp::set_threads_Eigen(n_threads);

//   // number of voxel index
//   const int ncol = mat_path_len.cols(); 

//   // enable multi-thread
//   #pragma omp parallel for
//   for(int icol=0;icol<ncol;icol++){ // voxel index loop
//     const Eigen::VectorXf vecxf_icol = mat_path_len.col(icol);

//     // sum of SUM_i_det{vecxf_icol[i_det]} > 0 means that the voxel is hit by any beam from detector element.
//     if( vecxf_icol.sum() > 0) callVoxel(icol)->incr_n_hit_det();
//   } // voxel index loop end
// }

// Function to check if vec_vec_vec_Voxel is the same
bool Grid3dVoxel::is_vec_vec_vec_Voxel_same( const Grid3dVoxel &g3vox_in ) const
{
  // Get number of bins
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  const int nbinz = get_nbinz();

  // Return false if number of bins differs
  if( nbinx != g3vox_in.get_nbinx() ) return false;
  if( nbiny != g3vox_in.get_nbiny() ) return false;
  if( nbinz != g3vox_in.get_nbinz() ) return false;

  // If number of bins is the same, compare voxel contents.
  for(int iz=0;iz<nbinz;iz++){
    for(int iy=0;iy<nbiny;iy++){
      for(int ix=0;ix<nbinx;ix++){
        const Voxel& vox1 = this->getVoxel(ix,iy,iz);
        const Voxel& vox2 = g3vox_in.getVoxel(ix,iy,iz);
        if( vox1 != vox2 ) return false;
      }
    }
  }
  return true;
}

//=======================================
// binary I/O functions
//=======================================

// save vec_vec_vec_Voxel to std::ofstream
void Grid3dVoxel::save_vec_vec_vec_Voxel( std::ofstream &ofs ) const
{
  // size_z records the depth (z) dimension size.
  const std::size_t size_z = vec_vec_vec_Voxel.size();

  // size_y is derived from the y dimension of vec_vec_vec_Voxel[0].
  // If vec_vec_vec_Voxel is empty (size_z == 0), size_y becomes zero.
  const std::size_t size_y = size_z > 0 ? vec_vec_vec_Voxel.at(0).size() : 0;

  // size_x reflects the x dimension from vec_vec_vec_Voxel[0][0].
  // If the prior dimensions are zero (size_y == 0 or size_z == 0), size_x is also zero.
  const std::size_t size_x = size_y > 0 ? vec_vec_vec_Voxel.at(0).at(0).size() : 0;

  // Save size of each dimension
  io_binary::write_binary(ofs, size_z);
  io_binary::write_binary(ofs, size_y);
  io_binary::write_binary(ofs, size_x);

  // Save Voxel object information
  for (const auto& vec_vec : vec_vec_vec_Voxel) {
    for (const auto& vec : vec_vec) {
      for (const auto& voxel : vec) {
        voxel.save(ofs);
      }
    }
  }
  if( ofs.fail() ) THROW_ERROR("std::ofstream& ofs failed.");
}

// load vec_vec_vec_Voxel to std::ifstream
void Grid3dVoxel::load_vec_vec_vec_Voxel(std::ifstream& ifs)
{
  // Load size of each dimension
  const std::size_t size_z = io_binary::read_binary<std::size_t>(ifs);
  const std::size_t size_y = io_binary::read_binary<std::size_t>(ifs);
  const std::size_t size_x = io_binary::read_binary<std::size_t>(ifs);

  // Reserve vector size
  vec_vec_vec_Voxel.resize(size_z, std::vector<std::vector<Voxel>>(size_y, std::vector<Voxel>(size_x)));

  // Load Voxel object information
  for (auto& vec_vec : vec_vec_vec_Voxel) {
    for (auto& vec : vec_vec) {
      for (auto& voxel : vec) {
        voxel.load(ifs);
      }
    }
  }
  if( ifs.fail() ) THROW_ERROR("std::ifstream& ifs failed.");
}

// save Grid3dVoxel to std::ofstream
void Grid3dVoxel::save( std::ofstream &ofs ) const
{
  // save Grid3d 
  this->Grid3d::save(ofs);

  // save name
  io_binary::write_string(ofs, name);

  // save n_vox_exist
  io_binary::write_binary(ofs, n_vox_exist);

  // save vec_vec_vec_Voxel
  save_vec_vec_vec_Voxel(ofs);

  // save uqiv_container
  uqiv_container.save(ofs);

  if( ofs.fail() ) THROW_ERROR("std::ofstream& ofs failed.");
}

// save Grid3dVoxel to fs::path
void Grid3dVoxel::save( const fs::path& pathout) const
{
  LOG_INFO("to pathout={}",pathout.string());
  std::ofstream ofs = io_binary::open_ofstream(pathout);
  save(ofs); ofs.close();
}

// load Grid3dVoxel from std::ifstream
void Grid3dVoxel::load( std::ifstream &ifs )
{
  // load Grid3d
  this->Grid3d::load(ifs);

  // load name
  name = io_binary::read_string(ifs);

  // load n_vox_exist
  n_vox_exist = io_binary::read_binary<int>(ifs);

  // load vec_vec_vec_Voxel
  load_vec_vec_vec_Voxel(ifs);

  // load uqiv_container
  uqiv_container.load(ifs);

  // rebuild flat array for O(1) uqiv lookup after loading
  build_flat_uqiv();

  // hit counters are not serialized; allocate zero-filled vectors so that
  // record_hit_ele/record_hit_det work after a checkpoint resume
  allocate_hit_data();

  if( ifs.fail() ) THROW_ERROR("std::ifstream& ifs failed.");
}

// load Grid3dVoxel from fs::path
void Grid3dVoxel::load( const fs::path &path_in )
{
  std::ifstream ifs = io_binary::open_ifstream(path_in);
  load(ifs); ifs.close();
}
