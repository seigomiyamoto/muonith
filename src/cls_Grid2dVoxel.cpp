// cls_Grid2dVoxel.cpp
#include "cls_Grid2dVoxel.hpp"
#include "ns_myapp.hpp"
#include "ns_mymacro.hpp"
#include "ns_mylogger.hpp"
#include <cstdio>
#include <cstdint>
#include "spdlog_pch.hpp"

//####################################################################################
//####################################################################################
// class Grid2dVoxel
// as derived class of the base class Grid2d
//####################################################################################
//####################################################################################

// non-equality operator
bool Grid2dVoxel::operator!=(const Grid2dVoxel& other) const
{
  // Name mismatch is intentionally ignored.
  // if (name != other.name) return true;

#ifdef NODEBUG
  if (zmin != other.zmin) return true;
  if (zmax != other.zmax) return true;
  if (!is_vec_vec_Voxel_same(other)) return true;
#else
  if (zmin != other.zmin) { LOG_WARN("!=: zmin mismatch ({} vs {})", zmin, other.zmin); return true; }
  if (zmax != other.zmax) { LOG_WARN("!=: zmax mismatch ({} vs {})", zmax, other.zmax); return true; }
  if (!is_vec_vec_Voxel_same(other)) { LOG_WARN("!=: vec_vec_Voxel content mismatch"); return true; }
#endif

  return false;
}

// mutable voxel access by indices
Voxel& Grid2dVoxel::callVoxel( const int ix, const int iy )
{
  // if out of range, throw error
  check_ix_inside(ix);
  check_iy_inside(iy);
  return vec_vec_Voxel.at(iy).at(ix);
}

// mutable voxel access by coordinates
Voxel& Grid2dVoxel::callVoxel( const double x, const double y)
{
  const int ix = get_ix(x);
  const int iy = get_iy(y);
  // if out of range, throw error
  check_ix_inside(ix);
  check_iy_inside(iy);
  return vec_vec_Voxel.at(iy).at(ix);
}

// const voxel access by indices
const Voxel& Grid2dVoxel::getVoxel( const int ix, const int iy ) const
{
  // if out of range, throw error
  check_ix_inside(ix);
  check_iy_inside(iy);
  return vec_vec_Voxel.at(iy).at(ix);
}

// const voxel access by coordinates
const Voxel& Grid2dVoxel::getVoxel( const double x, const double y) const
{
  const int ix = get_ix(x);
  const int iy = get_iy(y);
  // if out of range, throw error
  check_ix_inside(ix);
  check_iy_inside(iy);
  return vec_vec_Voxel.at(iy).at(ix);
}

// get density vector in (iy, ix) order
Eigen::VectorXf Grid2dVoxel::get_vecxf_density() const
{
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  Eigen::VectorXf vecxf_dens = Eigen::VectorXf::Zero(nbinx * nbiny);
  for (int iy = 0; iy < nbiny; ++iy) {
    for (int ix = 0; ix < nbinx; ++ix) {
      const int index = iy * nbinx + ix;
      vecxf_dens(index) = static_cast<float>(getVoxel(ix, iy).get_density());
    }
  }
  return vecxf_dens;
}

// allocate memory for std::vector<std::vector<Voxel>> vec_vec_Voxel
void Grid2dVoxel::vec_vec_memory_allocate(const int n_det){
  // Algorithm:
  // 1) Resize outer vector for y-direction.
  // 2) Resize inner vectors for x-direction.
  // 3) Initialize voxel flags and densities.
  // Note: n_det is retained for API compatibility but no longer used here;
  //       per-detector hit flags are now stored at Grid3dVoxel level.
  if (n_det < 0) {
    THROW_ERROR("Grid2dVoxel::vec_vec_memory_allocate: n_det must be non-negative");
  }
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  vec_vec_Voxel.resize(nbiny); // allocate y direction
  #pragma omp parallel for
  for(int iy=0;iy<nbiny;iy++){
    vec_vec_Voxel.at(iy).resize(nbinx); // allocate x direction
    for(int ix=0;ix<nbinx;ix++){
      Voxel& vox = callVoxel(ix,iy);
      vox.set_tf_exist(false); // set default value "false"
      vox.set_density(0.); // set default value 0.0
    }
  }
}

// OpenMP-enabled allocation wrapper
void Grid2dVoxel::mp_vec_vec_memory_allocate(const int n_det)
{
  vec_vec_memory_allocate(n_det);
}

// output function
void Grid2dVoxel::out_voxel( FILE *fout, const int ix, const int iy) const
{
  if (fout == nullptr) {
    THROW_ERROR("Grid2dVoxel::out_voxel: fout is null");
  }
  const Voxel& vox = getVoxel(ix,iy);
  const bool tf_exist = vox.get_tf_exist();
  const double xcnt = get_xcnt(ix);
  const double ycnt = get_ycnt(iy);
  const double zcnt = 0.5*(zmin+zmax);
  const double density = vox.get_density();
  // Per-detector hit flags are now stored at Grid3dVoxel level.
  // Output n_det=0 and n_hit_det=0 for format compatibility.
  const int n_det = 0;
  const int n_hit_det = 0;

  fprintf(fout,"%d %12.9E %12.9E %12.9E %12.9E %3d"
  ,(int)tf_exist,xcnt,ycnt,zcnt,density,n_det);

  fprintf(fout," %3d\n",n_hit_det);
}

// apply out_voxel for all ix,iy
void Grid2dVoxel::out_voxel_all( FILE *fout ) const
{
  if (fout == nullptr) {
    THROW_ERROR("Grid2dVoxel::out_voxel_all: fout is null");
  }
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      out_voxel(fout,ix,iy);
    }
  }
}

void Grid2dVoxel::out_voxel_all_binary( FILE *fout ) const
{
  if (fout == nullptr) {
    THROW_ERROR("Grid2dVoxel::out_voxel_all_binary: fout is null");
  }
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();

  auto write_bytes = [fout](const void* ptr, size_t size) {
    if (std::fwrite(ptr, size, 1, fout) != 1) {
      THROW_ERROR("Grid2dVoxel::out_voxel_all_binary failed to write data");
    }
  };

  // Algorithm:
  // 1) Loop over iy/ix in storage order.
  // 2) Write flags and geometry (xcnt, ycnt, zcnt).
  // 3) Write density and detector hit flags.
  // 4) Write total hit count.
  for (int iy = 0; iy < nbiny; ++iy) {
    for (int ix = 0; ix < nbinx; ++ix) {
      const Voxel& vox = getVoxel(ix, iy);
      const std::uint8_t tf_exist = static_cast<std::uint8_t>(vox.get_tf_exist());
      write_bytes(&tf_exist, sizeof(tf_exist));

      const double xcnt = get_xcnt(ix);
      const double ycnt = get_ycnt(iy);
      const double zcnt = 0.5 * (zmin + zmax);
      write_bytes(&xcnt, sizeof(xcnt));
      write_bytes(&ycnt, sizeof(ycnt));
      write_bytes(&zcnt, sizeof(zcnt));

      const double density = vox.get_density();
      write_bytes(&density, sizeof(density));

      // Per-detector hit flags are now stored at Grid3dVoxel level.
      // Output n_det=0 and n_hit_det=0 for format compatibility.
      const std::int32_t n_det = 0;
      write_bytes(&n_det, sizeof(n_det));

      const std::int32_t n_hit_det = 0;
      write_bytes(&n_hit_det, sizeof(n_hit_det));
    }
  }
}

// get AABB2d for (ix, iy)
AABB2d Grid2dVoxel::get_AABB2d(const int ix, const int iy) const
{
  check_ix_inside(ix);
  check_iy_inside(iy);
  // get coordinates aligned to x/y axes
  const double xmin_AABB = this->get_xlow(ix);
  const double xmax_AABB = this->get_xup(ix);
  const double ymin_AABB = this->get_ylow(iy);
  const double ymax_AABB = this->get_yup(iy);

  // build 2D AABB
  const Eigen::Vector2d v2_AABB_min(xmin_AABB,ymin_AABB);
  const Eigen::Vector2d v2_AABB_max(xmax_AABB,ymax_AABB);
  const AABB2d aabb2d(v2_AABB_min,v2_AABB_max);

  return aabb2d;
}

// get density and AABB2d
std::tuple< double, AABB2d >
  Grid2dVoxel::get_density_AABB2d(const int ix, const int iy ) const
{
  // get voxel information
  const Voxel& vox = this->getVoxel(ix,iy);
  const double density = vox.get_density();
  // return density and AABB2d
  return std::make_tuple(density,get_AABB2d(ix,iy));
}




// get AABB3d for (ix, iy)
AABB3d Grid2dVoxel::get_AABB3d(const int ix, const int iy) const
{
  check_ix_inside(ix);
  check_iy_inside(iy);
  // get coordinates aligned to x/y axes
  const double xmin_AABB = this->get_xlow(ix);
  const double xmax_AABB = this->get_xup(ix);
  const double ymin_AABB = this->get_ylow(iy);
  const double ymax_AABB = this->get_yup(iy);
  const double zmin_AABB = this->zmin;
  const double zmax_AABB = this->zmax;

  // build 3D AABB
  const Eigen::Vector3d v3_AABB_min(xmin_AABB,ymin_AABB,zmin_AABB);
  const Eigen::Vector3d v3_AABB_max(xmax_AABB,ymax_AABB,zmax_AABB);
  const AABB3d aabb3d(v3_AABB_min,v3_AABB_max);

  return aabb3d;
}

// get density and AABB3d
std::tuple< double, AABB3d >
  Grid2dVoxel::get_density_AABB3d(const int ix, const int iy ) const
{
  // get voxel information
  const Voxel& vox = this->getVoxel(ix,iy);
  const double density = vox.get_density();
  // return density and AABB3d
  return std::make_tuple(density,get_AABB3d(ix,iy));
}

// check if voxel storage is identical
bool Grid2dVoxel::is_vec_vec_Voxel_same( const Grid2dVoxel &g2vox_in ) const
{
  // get bin counts
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();

  // return false if counts differ
  if( nbinx != g2vox_in.get_nbinx() ) return false;
  if( nbiny != g2vox_in.get_nbiny() ) return false;

  // compare voxel contents
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      const Voxel& vox1 = this->getVoxel(ix,iy);
      const Voxel& vox2 = g2vox_in.getVoxel(ix,iy);
      if( vox1 != vox2 ) return false;
    }
  }
  return true;
}

// deep-copy voxel storage with size checks
void Grid2dVoxel::set_vec_vec_Voxel(
  const std::vector<std::vector<Voxel>> &vec_vec_Voxel_in )
{
  // get bin counts
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();

  if (vec_vec_Voxel_in.empty()) {
    if (nbinx != 0 || nbiny != 0) {
      THROW_ERROR("Grid2dVoxel::set_vec_vec_Voxel: input is empty but grid size is non-zero");
    }
    vec_vec_Voxel.clear();
    return;
  }

  if( nbiny != static_cast<int>(vec_vec_Voxel_in.size()) ) {
    THROW_ERROR("Grid2dVoxel::set_vec_vec_Voxel: nbiny mismatch");
  }

  const int size_x = static_cast<int>(vec_vec_Voxel_in.front().size());
  if( nbinx != size_x ) {
    THROW_ERROR("Grid2dVoxel::set_vec_vec_Voxel: nbinx mismatch");
  }

  for (int iy = 0; iy < nbiny; ++iy) {
    const int row_size = static_cast<int>(vec_vec_Voxel_in.at(iy).size());
    if (row_size != nbinx) {
      THROW_ERROR("Grid2dVoxel::set_vec_vec_Voxel: row size mismatch at iy={}", iy);
    }
  }

  // deep-copy voxel storage
  vec_vec_Voxel = vec_vec_Voxel_in;
}

//=======================================
// binary I/O functions
//=======================================

// save vec_vec_Voxel to std::ofstream
void Grid2dVoxel::save_vec_vec_Voxel( std::ofstream &ofs ) const
{
  // size_y is the y-direction size of vec_vec_Voxel
  const std::size_t size_y = vec_vec_Voxel.size();

  // size_x is the x-direction size of the first row
  const std::size_t size_x = size_y > 0 ? vec_vec_Voxel.at(0).size() : 0;

  // store dimensions
  io_binary::write_binary(ofs, size_y);
  io_binary::write_binary(ofs, size_x);

  // store voxel data in (iy, ix) order
  for (const auto& vec : vec_vec_Voxel) {
    for (const auto& voxel : vec) {
      voxel.save(ofs);
    }
  }
  if( ofs.fail() ) THROW_ERROR("Grid2dVoxel::save_vec_vec_Voxel: stream write failed");
}

// load vec_vec_Voxel to std::ifstream
void Grid2dVoxel::load_vec_vec_Voxel(std::ifstream& ifs)
{
  // read y-direction size
  const std::size_t size_y = io_binary::read_binary<std::size_t>(ifs);

  // read x-direction size
  const std::size_t size_x = io_binary::read_binary<std::size_t>(ifs);

  // read voxel data in (iy, ix) order
  vec_vec_Voxel.resize(size_y);
  for (auto& vec : vec_vec_Voxel) {
    vec.resize(size_x);
    for (auto& voxel : vec) {
      voxel.load(ifs);
    }
  }
  if( ifs.fail() ) THROW_ERROR("Grid2dVoxel::load_vec_vec_Voxel: stream read failed");
}



// save Grid2dVoxel to std::ofstream
void Grid2dVoxel::save( std::ofstream &ofs ) const
{
  // save Grid2d 
  this->Grid2d::save(ofs);

  // save name
  io_binary::write_string(ofs, name);

  // save zmin, zmax
  io_binary::write_binary(ofs, zmin);
  io_binary::write_binary(ofs, zmax);

  // save vec_vec_Voxel
  save_vec_vec_Voxel(ofs);

  if( ofs.fail() ) THROW_ERROR("Grid2dVoxel::save: stream write failed");
}

// save Grid2dVoxel to fs::path
void Grid2dVoxel::save( const fs::path& pathout) const
{
  std::ofstream ofs = io_binary::open_ofstream(pathout);
  save(ofs); ofs.close();
}

// load Grid2dVoxel from std::ifstream
void Grid2dVoxel::load( std::ifstream &ifs )
{
  // load Grid2d
  this->Grid2d::load(ifs);

  // load name
  name = io_binary::read_string(ifs);

  // load zmin, zmax
  zmin = io_binary::read_binary<double>(ifs);
  zmax = io_binary::read_binary<double>(ifs);

  // load vec_vec_Voxel
  load_vec_vec_Voxel(ifs);

  if( ifs.fail() ) THROW_ERROR("Grid2dVoxel::load: stream read failed");
}

// load Grid2dVoxel from fs::path
void Grid2dVoxel::load( const fs::path &path_in )
{
  std::ifstream ifs = io_binary::open_ifstream(path_in);
  load(ifs); ifs.close();
}
