// cls_Grid2dPillar.cpp
#include "ns_myapp.hpp"
#include "ns_mymacro.hpp"
#include "ns_mylogger.hpp"
#include "ns_io_binary.hpp"
#include <cstdint>
#include <cstdio>
#include <vector>
#include "cls_Grid2dPillar.hpp"
#include "cls_Grid2dXYZ.hpp"
#include "cls_Grid2dPillarParameters.hpp"
#include "cls_Grid2dVoxel.hpp"
#include "spdlog_pch.hpp"
#include "ns_geom_util.hpp"
#include "ns_iodir.hpp"

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

//##################################################
//##################################################
// class Grid2dPillar
// as derived class of the base class Grid2d
//##################################################
//##################################################

/// @brief not equal operator
bool Grid2dPillar::operator!=(const Grid2dPillar &other) const
{
  #ifdef NODEBUG
    if (static_cast<const Grid2d&>(*this) != static_cast<const Grid2d&>(other)) return true;
    if (vec_vec_Pillar != other.vec_vec_Pillar) return true;
  #else
    if (static_cast<const Grid2d&>(*this) != static_cast<const Grid2d&>(other)) { LOG_WARN("Grid2dPillar: Base Grid2d part differs"); return true; }
    if (vec_vec_Pillar != other.vec_vec_Pillar) { LOG_WARN("Grid2dPillar: vec_vec_Pillar differs"); return true; }
  #endif
  return false;
};


// constructors
//it calls convert_from_vec_xyz, read_dem
Grid2dPillar::Grid2dPillar( const std::filesystem::path &path_in
  , const double zmin, const double density_in
  , const bool tf_shift_x, const bool tf_shift_y
  , const double tolerance_ratio ) : Grid2dPillar() 
{
  read_dem(path_in,zmin,density_in,tf_shift_x,tf_shift_y,tolerance_ratio);
  set_name(path_in.filename().string());
  LOG_INFO("Grid2dPillar build. name = {}",name);
}

//it calls convert_from_vec_xyz, read_dem, and set_name
Grid2dPillar::Grid2dPillar( const std::string &name_in
  , const std::filesystem::path &path_in
  , const double zmin, const double density_in
  , const bool tf_shift_x, const bool tf_shift_y
  , const double tolerance_ratio ) : Grid2dPillar()
{
  read_dem(path_in, zmin, density_in, tf_shift_x, tf_shift_y,tolerance_ratio);
  set_name(name_in);
  LOG_INFO("Grid2dPillar build. name = {}",name);
}

// build from Grid2dPillar::Parameters
// it calls constructor above,
// it calls convert_from_vec_xyz, read_dem, and set_name
// it also calls add_density_structure for cylinder, dike, checkerboard
Grid2dPillar::Grid2dPillar( const Grid2dPillar::Parameters &prm )
  : Grid2dPillar(
        prm.name
      , prm.path_dem
      , prm.zmin
      , prm.initial_uniform_density
      , prm.tf_shift_x
      , prm.tf_shift_y
      , prm.tolerance_ratio
    )
{
  // add vertical cylinders
  for(const auto& param : prm.vec_vertical_cylinder_parameters) {
    add_density_structure(param);
  }

  // add vertical dikes
  for(const auto& param : prm.vec_vertical_dike_parameters) {
    add_density_structure(param);
  }

  // add vertical checkerboard
  for(const auto& param : prm.vec_vertical_checkerboard_parameters) {
    add_density_structure(param);
  }
}

// set uniform density to all Pillar
void Grid2dPillar::set_uniform_density( const double density_in )
{
  // DEBUG: check for NaN sources
  if (!std::isfinite(density_in)) {
    THROW_ERROR("Grid2dPillar::set_uniform_density: density_in is non-finite. density_in={}", density_in);
  }

  #pragma omp parallel for collapse(2)
  for(int iy=0;iy<get_nbiny();iy++){
    for(int ix=0;ix<get_nbinx();ix++){
      callPillar(ix,iy).set_density(density_in);
    }
  }
}


// allocate memory for std::vector<std::vector<Pillar>> vec_vec_Pillar;
void Grid2dPillar::vec_vec_memory_allocate()
{
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  vec_vec_Pillar.resize(nbiny); // allocate y direction
  for(int iy=0;iy<get_nbiny();iy++){
    vec_vec_Pillar.at(iy).resize(nbinx); // allocate x direction
    for(int ix=0;ix<nbinx;ix++){
      vec_vec_Pillar.at(iy).at(ix).set_values(0.,0.,0.);
    }
  }
}

// get/call functions
// get immutable reference of Pillar
const Pillar& Grid2dPillar::getPillar( const double x, const double y ) const
{
  const int ix = get_ix(x);
  const int iy = get_iy(y);
  // if out of range, throw error
  check_ix_inside(ix);
  check_iy_inside(iy);
  return vec_vec_Pillar.at(iy).at(ix);
}

const Pillar& Grid2dPillar::getPillar( const int ix, const int iy ) const
{
  // if out of range, throw error
  check_ix_inside(ix);
  check_iy_inside(iy);
  return vec_vec_Pillar.at(iy).at(ix);
}

Pillar& Grid2dPillar::callPillar( const double x, const double y )
{
  const int ix = get_ix(x);
  const int iy = get_iy(y);
  // if out of range, throw error
  check_ix_inside(ix);
  check_iy_inside(iy);
  return vec_vec_Pillar.at(iy).at(ix);
}

Pillar& Grid2dPillar::callPillar( const int ix, const int iy )
{
  // if out of range, throw error
  check_ix_inside(ix);
  check_iy_inside(iy);
  return vec_vec_Pillar.at(iy).at(ix);
};

// reconstruct data from std::vector<std::array<double,3>> &vec_xyz
void Grid2dPillar::convert_from_vec_xyz(
    std::vector<std::array<double,3>> &vec_xyz
  , const double zmin, const double density_in
  , const bool tf_shift_x, const bool tf_shift_y
  , const double tolerance_ratio )
{
  // DEBUG: check for NaN sources
  if (!std::isfinite(density_in)) {
    THROW_ERROR("Grid2dPillar::convert_from_vec_xyz: density_in is non-finite. density_in={}", density_in);
  }
  if (!std::isfinite(zmin)) {
    THROW_ERROR("Grid2dPillar::convert_from_vec_xyz: zmin is non-finite. zmin={}", zmin);
  }

  // Algorithm overview:
  // 1. Infer grid structure (x/y axes) from point cloud
  // 2. If irregular, sort points and retry
  // 3. Allocate Pillar array
  // 4. Assign z-values (elevation) to each Pillar

  // Step 1: Try to infer regular grid axes from point positions
  // This assumes vec_xyz contains regularly-spaced terrain points (e.g., from DEM)
  bool tf_interval_check_1st = false;
  bool tf_interval_check_2nd = false;

  tf_interval_check_1st
   = set_xy_axis_from_vec_xyz(
    vec_xyz,tf_shift_x,tf_shift_y, tolerance_ratio);

  // Step 2: If grid inference failed, sort points and retry
  // Points may be unsorted from file input; sorting enables grid detection
  if( tf_interval_check_1st == false ){
    tf_interval_check_1st = true;
    LOG_WARN("... set_xy_axis_from_vec_xyz ... failed.");
    LOG_WARN("need to sort vec_xyz by x and y.");
    LOG_INFO(", applying std::sort ... ");
    std::sort(vec_xyz.begin(), vec_xyz.end(), myapp::compare_array_yxz);

    tf_interval_check_2nd
     = set_xy_axis_from_vec_xyz(
      vec_xyz,tf_shift_x,tf_shift_y,tolerance_ratio);

    if( tf_interval_check_2nd == false ){
      THROW_ERROR_NAME("Grid2dPillar::convert_from_vec_xyz ... set_xy_axis_from_vec_xyz ... failed after sort.\n");
    }
  }
  LOG_INFO("... set_xy_axis_from_vec_xyz ... done. \n");

  // Step 3: Allocate memory for Pillar array (column-major: iy, ix)
  vec_vec_memory_allocate();
  LOG_INFO("... vec_vec_memory_allocate ... done. \n");

  // Display inferred grid information
  get_x_axis().out_info(spdlog::level::info);
  get_y_axis().out_info(spdlog::level::info);

  // Step 4: Assign elevation values to each Pillar from point cloud
  // Note: Some points may lie outside the computed grid range due to tolerance.
  // Such points (e.g., edge artifacts from DEM) are safely skipped.
  for( const auto& [x,y,z] : vec_xyz ){
    const int ix = get_ix(x);
    // Skip points outside grid bounds (e.g., edge points from DEM exceeding computed range)
    if( ix < 0 ) continue;
    if( ix >= get_nbinx() ) continue;
    const int iy = get_iy(y);
    if( iy < 0 ) continue;
    if( iy >= get_nbiny() ) continue;

    // Set Pillar vertical extent: from zmin (base) to z (terrain surface)
    const double z_lower = zmin; // elevation = zmin [meter a.s.l.]
    const double z_upper = z; // elevation,    [meter a.s.l.]
    if( z_lower < z_upper ){
      vec_vec_Pillar.at(iy).at(ix).set_values(
        z_lower, z_upper, density_in );
    }else{
      // If z <= zmin, create zero-density cuboid (handles below-base points)
      vec_vec_Pillar.at(iy).at(ix).set_values(
        z_lower, z_upper, 0.0);
    }
  }
}

// Read terrain from a .g2zbin (binary grid) or an ASCII xyz file.
// g2zbin path uses the Grid2dXYZ constructor; ASCII path uses
// myapp::read_vec_xyz + convert_from_vec_xyz.
void Grid2dPillar::read_dem(
    const std::filesystem::path &path_in
  , const double zmin, const double density_in
  , const bool tf_shift_x, const bool tf_shift_y
  , const double tolerance_ratio )
{
  if(path_in.empty()) THROW_ERROR_NAME("Grid2dPillar::read_dem: Filename not assigned");
  LOG_INFO("reading file... {}\n",path_in.string());

  if (path_in.extension() == ".g2zbin") {
    // g2zbin path: Grid2dXYZ constructor handles header + z raster loading.
    // Then copy axes and convert z values to Pillars.
    Grid2dXYZ g2xyz(path_in, tf_shift_x, tf_shift_y, tolerance_ratio);

    set_x_axis(g2xyz.get_x_axis());
    set_y_axis(g2xyz.get_y_axis());
    vec_vec_memory_allocate();

    const int nbinx = get_nbinx();
    const int nbiny = get_nbiny();

    // One pillar per cell: bottom fixed at zmin, top at the terrain height.
    // A cell with no usable terrain height becomes an empty (zero-density) pillar.
    #pragma omp parallel for collapse(2)
    for (int iy = 0; iy < nbiny; ++iy) {
      for (int ix = 0; ix < nbinx; ++ix) {
        const double z_dem = g2xyz.get_z(ix, iy);  // terrain height of this cell [m]

        const double z_lower = zmin;               // bottom: common base plane
        double z_upper = z_dem;                    // top: terrain height, clamped below

        // Terrain below the base plane (sea, lake bed) is raised up to zmin.
        // NaN (elevation missing in the DEM) fails this test and stays NaN.
        if (z_dem < zmin) z_upper = zmin;

        if (z_lower < z_upper) {
          // Pillar has thickness: fill it with the uniform initial density.
          vec_vec_Pillar.at(iy).at(ix).set_values(z_lower, z_upper, density_in);
        } else {
          // No thickness (terrain at or below zmin, or z_dem is NaN): empty pillar.
          // Force the top to zmin so a NaN z_dem is never stored; a NaN zmax
          // would make Pillar::is_z_inside() return true at every height.
          z_upper = zmin;
          vec_vec_Pillar.at(iy).at(ix).set_values(z_lower, z_upper, 0.0);
        }
      }
    }
  } else {
    // Existing path: read ASCII xyz via myapp::read_vec_xyz, then
    // infer grid structure and convert to Pillars.
    std::vector<std::array<double,3>>
      vec_xyz = myapp::read_vec_xyz(path_in);

    // Clamp terrain heights below zmin up to zmin (Uses OpenMP)
    #pragma omp parallel for
    for (std::size_t i = 0; i < vec_xyz.size(); i++) {
      if (vec_xyz.at(i).at(2) < zmin) vec_xyz.at(i).at(2) = zmin;
    }

    // Store to vec_vec_Pillar in Grid2dPillar
    convert_from_vec_xyz(vec_xyz,zmin,density_in,tf_shift_x,tf_shift_y,tolerance_ratio);
  }

  name = path_in.filename().string();
  std::string fname_str(name);
  std::string xname = "x_axis_" + fname_str;
  std::string yname = "y_axis_" + fname_str;
  set_x_axis_name(xname);
  set_y_axis_name(yname);
}

void Grid2dPillar::out( const std::filesystem::path& pathout ) const
{
  if (pathout.empty()) {
    THROW_ERROR_NAME("Grid2dPillar::out: Output file path is empty");
  }
  LOG_INFO("outputing {} ...",pathout.string());
  FILE *fout = fopen(pathout.c_str(),"wt");
  if(fout == NULL)
    THROW_ERROR_NAME("Grid2dPillar::out: Cannot open file '{}'"
    , pathout.string());
  int ix,iy;
  double xmin,xmax,ymin,ymax,zmin,zmax,density;
  for(iy=0;iy<get_nbiny();iy++){
    for(ix=0;ix<get_nbinx();ix++){
      const Pillar &cub = getPillar(ix,iy);
      xmin = get_xlow(ix);
      xmax = get_xup(ix);
      ymin = get_ylow(iy);
      ymax = get_yup(iy);
      zmin = cub.get_zmin();
      zmax = cub.get_zmax();
      density = cub.get_density();
      fprintf(fout,"%d %d %E %E %E %E %E %E %E\n"
      ,ix,iy,xmin,xmax,ymin,ymax,zmin,zmax,density);
    }
  }
  myapp::close(fout,pathout);
}

// disp nearest point and zmin & zmax value
void Grid2dPillar::disp_nearest_point( const double x_in, const double y_in ) const
{
  const double x_low = get_x_axis().get_lower_value(x_in);
  const double x_up  = get_x_axis().get_upper_value(x_in);
  const double y_low = get_y_axis().get_lower_value(y_in);
  const double y_up  = get_y_axis().get_upper_value(y_in);
  const double x_cub = 0.5*(x_low + x_up);
  const double y_cub = 0.5*(y_low + y_up);
  const Pillar& cub = getPillar(x_cub,y_cub);
  LOG_INFO("x={:.1f}, y={:.1f}, zmin={:.1f}, zmax={:.1f}"
  ,x_cub,y_cub,cub.get_zmin(),cub.get_zmax());
}

// Add vertical cylinder density structure to existing Grid2dPillar.
// For each grid cell whose center lies inside the ellipse, add delta_density.
void Grid2dPillar::add_density_structure( const VerticalCylinderParameters &prm )
{
  // DEBUG: check for NaN sources
  if (!std::isfinite(prm.delta_density)) {
    THROW_ERROR("Grid2dPillar::add_density_structure(Cylinder): delta_density is non-finite. name={}, delta_density={}", prm.name, prm.delta_density);
  }

  if( prm.tf_exec==false ){
    LOG_INFO("name={} tf_exec=false",prm.name);
  }
  else{
    // Get grid dimensions
    const int nbinx = get_nbinx();
    const int nbiny = get_nbiny();

    // Track whether at least one bin is inside the cylinder
    bool flg_1entry_at_least = false;

    // Check all grid cells to see if their centers lie inside the ellipse
    #pragma omp parallel for collapse(2)
    for(int iy=0;iy<nbiny;iy++){
      for(int ix=0;ix<nbinx;ix++){
        double ycnt = get_ycnt(iy);
        double xcnt = get_xcnt(ix);
        const Eigen::Vector2d v2_pos(xcnt,ycnt);
        bool tf_inside = geom_util::isInsideEllipse(prm.v2_pos_cnt,prm.v2_length,prm.angle_rot,v2_pos);

        // If not inside ellipse, skip to next bin
        if(tf_inside==false) continue;

        // If inside ellipse, add delta_density
        Pillar& cub = callPillar(ix,iy);
        cub.add_density(prm.delta_density);
        flg_1entry_at_least = true;
      }
    }

    // Throw error if no bins found inside cylinder
    if( flg_1entry_at_least == false )
      THROW_ERROR_NAME(
        "Grid2dPillar::add_density_structure: No bin entry found in cylinder '{}'"
        , prm.name);
  }
}

// Add vertical dike (rectangular block) density structure
void Grid2dPillar::add_density_structure( const VerticalDikeParameters &prm )
{
  // DEBUG: check for NaN sources
  if (!std::isfinite(prm.delta_density)) {
    THROW_ERROR("Grid2dPillar::add_density_structure(Dike): delta_density is non-finite. name={}, delta_density={}", prm.name, prm.delta_density);
  }

  if( prm.tf_exec==false ){
    LOG_INFO("name={} tf_exec=false",prm.name);
  }
  else{
    // Get grid dimensions
    const int nbinx = get_nbinx();
    const int nbiny = get_nbiny();

    // Track whether at least one bin is inside the dike
    bool flg_1entry_at_least = false;

    // Check all grid cells to see if their centers lie inside the rectangle
    #pragma omp parallel for collapse(2)
    for(int iy=0;iy<nbiny;iy++){
      for(int ix=0;ix<nbinx;ix++){
        double ycnt = get_ycnt(iy);
        double xcnt = get_xcnt(ix);
        const Eigen::Vector2d v2_pos(xcnt,ycnt);
        bool tf_inside = prm.aabb2d.is_inside(v2_pos,prm.angle_rot);

        // If not inside rectangle, skip to next bin
        if(tf_inside==false) continue;

        // If inside rectangle, add delta_density
        callPillar(ix,iy).add_density(prm.delta_density);
        flg_1entry_at_least = true;
      }
    }

    // Throw error if no bins found inside dike
    if( flg_1entry_at_least == false )
      THROW_ERROR_NAME(
        "Grid2dPillar::add_density_structure: No bin entry found in dike '{}'"
        , prm.name);
  }
}

// Add vertical checkerboard density pattern
void Grid2dPillar::add_density_structure( const VerticalCheckerBoardParameters &prm )
{
  // DEBUG: check for NaN sources
  if (!std::isfinite(prm.delta_density)) {
    THROW_ERROR("Grid2dPillar::add_density_structure(Checkerboard): delta_density is non-finite. name={}, delta_density={}", prm.name, prm.delta_density);
  }
  if (!std::isfinite(prm.delta_density_offset)) {
    THROW_ERROR("Grid2dPillar::add_density_structure(Checkerboard): delta_density_offset is non-finite. name={}, delta_density_offset={}", prm.name, prm.delta_density_offset);
  }

  if( prm.tf_exec==false ){
    LOG_INFO("name={} tf_exec=false",prm.name);
  }
  else{
    // Get grid dimensions
    const int nbinx = get_nbinx();
    const int nbiny = get_nbiny();

    // Integer-index-based checkerboard sign determination
    // (avoids floating-point error in floor() with non-integer grid intervals)
    const int ix_cnt = get_ix(prm.v2_pos_cnt.x());
    const int iy_cnt = get_iy(prm.v2_pos_cnt.y());
    const int mx = static_cast<int>(std::lround(prm.v2_length.x() / get_x_interval()));
    const int my = static_cast<int>(std::lround(prm.v2_length.y() / get_y_interval()));

    // Determine which checkerboard cell (positive or negative) each bin belongs to
    #pragma omp parallel for collapse(2)
    for(int iy=0;iy<nbiny;iy++){
      for(int ix=0;ix<nbinx;ix++){
        double ycnt = get_ycnt(iy);
        if( ycnt <  prm.aabb2d.ymin() ) continue;
        if( ycnt >= prm.aabb2d.ymax() ) continue;
        double xcnt = get_xcnt(ix);
        if( xcnt <  prm.aabb2d.xmin() ) continue;
        if( xcnt >= prm.aabb2d.xmax() ) continue;

        // Checkerboard sign via integer block indices
        const int bx = floor_div(ix - ix_cnt, mx);
        const int by = floor_div(iy - iy_cnt, my);
        const bool odd_block = (bx + by) % 2 != 0;
        const double sign = odd_block ? -1.0 : 1.0;

        // Add positive or negative delta_density based on checkerboard pattern
        callPillar(ix,iy).add_density(
          prm.delta_density_offset + sign * prm.delta_density
        );
      }
    }
  }
}

Grid2dVoxel Grid2dPillar::make_cross_section_voxel(
  const double z, const int n_detector) const
{
  Grid2dVoxel g2vox(get_x_axis(), get_y_axis(), z, z, n_detector);
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  for(int iy=0; iy<nbiny; ++iy){
    for(int ix=0; ix<nbinx; ++ix){
      const Pillar& cub = getPillar(ix, iy);
      Voxel& vox = g2vox.callVoxel(ix, iy);
      const bool tf_exist = cub.is_z_inside(z);
      vox.set_tf_exist(tf_exist);
      if (!tf_exist) {
        vox.set_density(0.0);
      } else {
        vox.set_density(cub.get_density());
      }
    }
  }
  return g2vox;
}

void Grid2dPillar::out_cross_section_z_header(
  FILE *fout, const CrossSectionZParameters& prm_zcross) const
{
  fprintf(fout,"# z_info %E %E %E\n",
    prm_zcross.zmin, prm_zcross.zmax, prm_zcross.zstep);
  fprintf(fout,"# x_info %E %E %E\n",
    prm_zcross.xmin, prm_zcross.xmax, prm_zcross.xstep);
  fprintf(fout,"# y_info %E %E %E\n",
    prm_zcross.ymin, prm_zcross.ymax, prm_zcross.ystep);
  fprintf(fout,"# n_detector %d\n", prm_zcross.n_detector);
  const int nz = (prm_zcross.zmax - prm_zcross.zmin) / prm_zcross.zstep + 1;
  fprintf(fout,"# nz %d\n", nz);
}

void Grid2dPillar::out_cross_section_z_all(
  const std::filesystem::path &filepath,
  const CrossSectionZParameters &prm_zcross) const
{
  if (prm_zcross.output_binary) {
    // Binary content gets a .tmpbin extension so it is not mistaken for text.
    std::filesystem::path filepath_bin = filepath;
    filepath_bin.replace_extension(".tmpbin");
    out_cross_section_z_all_binary(filepath_bin, prm_zcross);
    return;
  }
  LOG_INFO("filepath={} ... ", filepath.string());
  FILE *fout = myapp::get_fout(filepath);
  out_cross_section_z_header(fout, prm_zcross);

  const double zmin = prm_zcross.zmin;
  const double zmax = prm_zcross.zmax;
  const double zstep = prm_zcross.zstep;
  for(double z = zmin; z <= zmax; z += zstep){
    const Grid2dVoxel g2vox = make_cross_section_voxel(z, prm_zcross.n_detector);
    g2vox.out_voxel_all(fout);
    LOG_INFO("z={} done", z);
  }
  myapp::close(fout, filepath);
}

void Grid2dPillar::out_cross_section_z_all_binary(
  const std::filesystem::path &filepath,
  const CrossSectionZParameters &prm_zcross) const
{
  LOG_INFO("filepath={} ... ", filepath.string());
  FILE *fout = myapp::get_fout_binary(filepath);

  auto write_bytes = [fout](const void* ptr, size_t size) {
    if (std::fwrite(ptr, size, 1, fout) != 1) {
      THROW_ERROR("Grid2dPillar::out_cross_section_z_all_binary failed to write data");
    }
  };

  auto write_value = [&write_bytes](const auto& value) {
    write_bytes(&value, sizeof(value));
  };

  const char magic[8] = {'G','2','Z','B','I','N','\0','\0'};
  write_bytes(magic, sizeof(magic));

  const std::uint32_t version = 1;
  write_value(version);

  const std::uint32_t reserved = 0;
  write_value(reserved);

  const double header_doubles[] = {
    prm_zcross.xmin, prm_zcross.xmax, prm_zcross.xstep,
    prm_zcross.ymin, prm_zcross.ymax, prm_zcross.ystep,
    prm_zcross.zmin, prm_zcross.zmax, prm_zcross.zstep
  };
  for (const double value : header_doubles) {
    write_value(value);
  }

  const int32_t nbinx = get_nbinx();
  const int32_t nbiny = get_nbiny();
  std::vector<double> z_values;
  for (double z = prm_zcross.zmin; z <= prm_zcross.zmax; z += prm_zcross.zstep) {
    z_values.push_back(z);
  }
  const int32_t nz = static_cast<int32_t>(z_values.size());

  const int32_t header_ints[] = {
    prm_zcross.n_detector, nbinx, nbiny, nz
  };
  for (const int32_t value : header_ints) {
    write_value(value);
  }

  const std::int64_t n_records =
    static_cast<std::int64_t>(nbinx) *
    static_cast<std::int64_t>(nbiny) *
    static_cast<std::int64_t>(nz);
  write_value(n_records);

  for (const double z : z_values) {
    const Grid2dVoxel g2vox = make_cross_section_voxel(z, prm_zcross.n_detector);
    g2vox.out_voxel_all_binary(fout);
    LOG_INFO("z={} done", z);
  }

  myapp::close(fout, filepath);
}

// Output x, y, density data at specified elevation z_in
void Grid2dPillar::out_cross_section_z(
  const std::filesystem::path& filepath, const double z_in) const
{
  LOG_DEBUG(
    "z_in={}, saving filepath={} ..."
    , z_in,filepath.string());
  FILE *fout = fopen(filepath.c_str(),"wt");
  if(fout == NULL)
    THROW_ERROR_NAME(
      "Grid2dPillar::out_cross_section_z: Cannot open file '{}'"
      , filepath.string());
  const int nbinx = get_nbinx();
  const int nbiny = get_nbiny();
  for(int iy=0;iy<nbiny;iy++){
    double ycnt = get_ycnt(iy);
    for(int ix=0;ix<nbinx;ix++){
      double xcnt = get_xcnt(ix);
      const Pillar& cub = getPillar(ix,iy);
      const double zmin = cub.get_zmin();
      const double zmax = cub.get_zmax();
      if( z_in <  zmin ) continue;
      if( z_in >= zmax ) continue;
      const double density = cub.get_density();
      fprintf(fout,"%E %E %E\n",xcnt,ycnt,density);
    }
  }
  myapp::close(fout,filepath);
}

// Get density and AABB3d (axis-aligned bounding box) for specified grid cell
std::tuple< double, AABB3d >
  Grid2dPillar::get_density_AABB3d(const int ix, const int iy) const
{
  // Get face coordinates along x and y axes
  const double xmin_AABB = this->get_xlow(ix);
  const double xmax_AABB = this->get_xup(ix);
  const double ymin_AABB = this->get_ylow(iy);
  const double ymax_AABB = this->get_yup(iy);

  // Get face coordinates along z axis
  const Pillar& cub = this->getPillar(ix,iy);
  const double density = cub.get_density();
  const double zmin_AABB = cub.get_zmin();
  const double zmax_AABB = cub.get_zmax();

  // Construct cuboid AABB
  const Eigen::Vector3d v3_AABB_min(xmin_AABB,ymin_AABB,zmin_AABB);
  const Eigen::Vector3d v3_AABB_max(xmax_AABB,ymax_AABB,zmax_AABB);
  const AABB3d aabb3d(v3_AABB_min,v3_AABB_max);

  return std::make_tuple(density,aabb3d);
}

// return AABB3d from tuple of hit cuboid
AABB3d Grid2dPillar::get_AABB3d(const Ixiy &ixiy_hit_cub) const
{
  // get ix iy from tuple
  const auto [ix,iy] = ixiy_hit_cub;

  // Get face coordinates along x and y axes
  const double xmin_AABB = this->get_xlow(ix);
  const double xmax_AABB = this->get_xup(ix);
  const double ymin_AABB = this->get_ylow(iy);
  const double ymax_AABB = this->get_yup(iy);

  // Get face coordinates along z axis
  const Pillar& cub = this->getPillar(ix,iy);
  const double zmin_AABB = cub.get_zmin();
  const double zmax_AABB = cub.get_zmax();

  // Construct cuboid AABB
  const Eigen::Vector3d v3_AABB_min(xmin_AABB,ymin_AABB,zmin_AABB);
  const Eigen::Vector3d v3_AABB_max(xmax_AABB,ymax_AABB,zmax_AABB);
  const AABB3d aabb3d(v3_AABB_min,v3_AABB_max);

  // return tuple of v3_AABB_min,v3_AABB_max
  return aabb3d;
}

// get dz bewtween zmax of Grid2dPillar from std::array<double,3> &xyz
double Grid2dPillar::get_dz( const std::array<double,3> &xyz ) const
{
  const auto [x,y,z] = xyz;
  // if out of range, return const_dz_no_hit
  if( is_inside(x,y) == false ) return const_dz_no_hit;
  const int ix = get_ix(x);
  const int iy = get_iy(y);

  const Pillar& cub = getPillar(ix,iy);
  const double zcub = cub.get_zmax();
  const double dz = z- zcub;
  return dz;
}

// get dz bewtween zmax of Grid2dPillar from std::tuple<double,double,double,int,int,int> &tp_xyzrgb
double Grid2dPillar::get_dz( const std::tuple<double,double,double,int,int,int> &tp_xyzrgb ) const
{
  const auto [x,y,z,r,g,b] = tp_xyzrgb;
  // if out of range, return const_dz_no_hit
  if( is_inside(x,y) == false ) return const_dz_no_hit;
  const int ix = get_ix(x);
  const int iy = get_iy(y);

  const Pillar& cub = getPillar(ix,iy);
  const double zcub = cub.get_zmax();
  const double dz = z- zcub;
  return dz;
}

/// @brief get minimum zmax of Pillar within radius from (x_center,y_center)
double Grid2dPillar::getMinimumZmaxCircle(
  const double x_center, const double y_center, const double radius ) const
{
  // if x_center,y_center is out of range, THROW_ERROR
  if( is_inside(x_center,y_center) == false )
    THROW_ERROR_NAME("Grid2dPillar::getMinimumZmaxCircle : x_center,y_center is out of range.");
  
  // Collect all grid cells whose center falls inside the circle.
  std::vector< Ixiy > vec_ixiy = Grid2d::get_vec_ixiy_in_circle(x_center, y_center, radius);

  // When radius < grid_interval, the circle may contain no cell center.
  // Fall back to the single cell that contains the center point.
  if (vec_ixiy.empty()) {
    LOG_WARN_ND( "no grid cell found in circle. "
      "center=({:.6E}, {:.6E}), radius={:.2f}, grid_interval=({:.2f}, {:.2f}). "
      "Falling back to nearest cell.",
      x_center, y_center, radius, get_x_interval(), get_y_interval());
    // Use the cell that contains (x_center, y_center) as the best approximation.
    const int ix = get_ix(x_center);
    const int iy = get_iy(y_center);
    return getPillar(ix, iy).get_zmax();
  }

  // Find the minimum zmax among all cells inside the circle.
  double zmax_min = std::numeric_limits<double>::max();
  for( const auto& [ix,iy] : vec_ixiy ){
    const Pillar& cub = getPillar(ix,iy);
    const double zmax = cub.get_zmax();
    if( zmax < zmax_min ) zmax_min = zmax;
  }
  return zmax_min;
}


// get vec_dz from std::vector<std::array<double,3>> &vec_xyz
std::vector<double> Grid2dPillar::get_vec_dz(
  const std::vector<std::array<double,3>> &vec_xyz ) const
{
  std::vector<double> vec_dz;
  for( const auto& tp_xyz : vec_xyz ){
    const double dz = get_dz(tp_xyz);
    // if out of range, skip
    if(dz==const_dz_no_hit) continue;
    vec_dz.push_back(dz);
  }
  return vec_dz;
}

double Grid2dPillar::calc_depth_from_zmax(
  const Eigen::Vector2d& v2_pos_obj,
  const Ray3d& ray_det,
  const double elev_ang_center_rad) const
{
  const Pillar& cub = getPillar(v2_pos_obj.x(), v2_pos_obj.y());
  const double z_surface = cub.get_zmax();

  const double dx = v2_pos_obj.x() - ray_det.x();
  const double dy = v2_pos_obj.y() - ray_det.y();
  const double horizontal_distance = std::hypot(dx, dy);

  const double tan_elev = std::tan(elev_ang_center_rad);
  if( std::isfinite(tan_elev) == false ){
    THROW_ERROR("Grid2dPillar::calc_depth_from_zmax invalid elev_ang_center_rad.");
  }

  const double z_on_ray = ray_det.z() + horizontal_distance * tan_elev;
  return z_surface - z_on_ray;
}

// get vec_xyz,dz from std::vector<std::array<double,3>> &vec_xyz
std::vector<std::array<double,4>>
  Grid2dPillar::get_vec_xyzdz(
    const std::vector<std::array<double,3>> &vec_xyz ) const
{
  std::vector<std::array<double,4>> vec_xyz_dz;
  for( const auto& xyz : vec_xyz ){
    const double dz = get_dz(xyz);
    // if out of range, skip
    if(dz==const_dz_no_hit) continue;
    const auto [x,y,z] = xyz;
    vec_xyz_dz.push_back( {x,y,z,dz} );
  }
  return vec_xyz_dz;
}

// get vec_xyz,dz from std::vector<std::tuple<double,double,double,int,int,int>> &vec_xyzrgb
std::vector<std::tuple<double,double,double,double,int,int,int>>
Grid2dPillar::get_vec_xyzdzrgb(
  const std::vector<std::tuple<double,double,double,int,int,int>> &vec_xyzrgb ) const
{
  std::vector<std::tuple<double,double,double,double,int,int,int>> vec_xyzdzrgb;
  for( const auto& tp_xyzrgb : vec_xyzrgb ){
    const double dz = get_dz(tp_xyzrgb);
    // if out of range, skip
    if(dz==const_dz_no_hit) continue;
    const auto &[x,y,z,r,g,b] = tp_xyzrgb;
    vec_xyzdzrgb.push_back( std::make_tuple(x,y,z,dz,r,g,b) );
  }
  return vec_xyzdzrgb;
}

//===============================
// path length calc function
//===============================

// Compute path length through a specific Pillar intersected by a ray.
// This is primarily called during muon radiography path calculations.
// Algorithm:
// 1. Get AABB (bounding box) for the Pillar at grid indices (ix, iy)
// 2. Compute ray-AABB intersection parameters (tmin, tmax)
// 3. If no intersection, return PATH_NO_HIT (0.0)
// 4. If ray starts inside cuboid (tmin < 0), clamp tmin to 0
// 5. Return path length = tmax - tmin
double Grid2dPillar::get_delta_path(
  const Ixiy &tpl_hit_cub, const Ray3d &ray3d ) const
{
  // Get axis-aligned bounding box for the cuboid at specified grid indices
  const AABB3d aabb3d = get_AABB3d(tpl_hit_cub);

  // Compute ray-box intersection using slab method
  // Returns: [intersection_flag, entry_parameter, exit_parameter]
  auto [tf_intersect,tmin_raw,tmax] = ray3d.is_intersect(aabb3d);

  // No intersection: ray misses the cuboid entirely
  if( tf_intersect == false ) return PATH_NO_HIT;

  // Handle ray starting inside cuboid: clamp entry point to ray origin
  // tmin < 0 means ray origin is inside box; use 0 to measure from origin
  const double tmin = (tmin_raw < 0) ? 0.0 : tmin_raw;
  const double delta_path = tmax - tmin;

  return delta_path;
}

// Get terrain elevation profile (zmax - ray_z) along Ray3d
std::vector<std::tuple<double,double>>
  Grid2dPillar::get_dzmax_profile(
    const Ray3d &ray3d, const double max_distance ) const
{
  // storage of zmax
  std::vector<std::tuple<double,double>> vec_dist_dzmax;

  // get Ray2d from Ray3d
  const Ray2d ray2d = ray3d.toRay2d();

  // get hit boxes index of ray2d from Grid2d
  std::vector<Ixiy> vec_ixiy_hitcub = Grid2d::get_hit_boxes_index(ray2d);

  for( const auto [ix,iy] : vec_ixiy_hitcub ){
    const Pillar& cub = getPillar(ix,iy);
    const double dzmax = cub.get_zmax() - ray3d.z();
    const double x = get_xcnt(ix);
    const double y = get_ycnt(iy);
    const Eigen::Vector2d v2_pos(x,y);
    const Eigen::Vector2d v2_diff = v2_pos - ray2d.pos();
    const double distance = v2_diff.norm();
    if( distance > max_distance ) continue;
    vec_dist_dzmax.push_back( std::make_tuple(distance,dzmax) );
  }
  return vec_dist_dzmax;
}

// save dzmax profile to file
void Grid2dPillar::out_dzmax_profile(
    std::filesystem::path &outputpath
  , const Ray3d &ray3d , const double max_distance ) const
{
  // get zmax profile
  std::vector<std::tuple<double,double>> vec_dist_zmax = get_dzmax_profile(ray3d,max_distance);

  // output to file
  FILE *fout = fopen(outputpath.c_str(),"wt");
  if(fout == NULL)
    THROW_ERROR_NAME(
      "Grid2dPillar::out_dzmax_profile: Cannot open file '{}'"
      , outputpath.string());
  for( const auto& [distance,zmax] : vec_dist_zmax ){
    fprintf(fout,"%E %E\n",distance,zmax);
  }
  myapp::close(fout,outputpath);
}

//----------------------------------------------------------------------
// Binary I/O
//----------------------------------------------------------------------

void Grid2dPillar::save(std::ofstream& ofs) const
{
  Grid2d::save(ofs);
  io_binary::write_string(ofs, name);
  const uint64_t nbiny = vec_vec_Pillar.size();
  io_binary::write_binary(ofs, nbiny);
  for (const auto& vec : vec_vec_Pillar) {
    const uint64_t nbinx = vec.size();
    io_binary::write_binary(ofs, nbinx);
    for (const auto& cub : vec) {
      cub.save(ofs);
    }
  }
}

void Grid2dPillar::load(std::ifstream& ifs)
{
  Grid2d::load(ifs);
  name = io_binary::read_string(ifs);
  const uint64_t nbiny = io_binary::read_binary<uint64_t>(ifs);
  vec_vec_Pillar.resize(nbiny);
  for (auto& vec : vec_vec_Pillar) {
    const uint64_t nbinx = io_binary::read_binary<uint64_t>(ifs);
    vec.resize(nbinx);
    for (auto& cub : vec) {
      cub.load(ifs);
    }
  }
}

void Grid2dPillar::save(const std::filesystem::path& pathout) const
{
  std::ofstream ofs(pathout, std::ios::binary);
  if (!ofs) THROW_ERROR("Grid2dPillar::save: Cannot open file '{}'", pathout.string());
  save(ofs);
}

void Grid2dPillar::load(const std::filesystem::path& path_in)
{
  std::ifstream ifs(path_in, std::ios::binary);
  if (!ifs) THROW_ERROR("Grid2dPillar::load: Cannot open file '{}'", path_in.string());
  load(ifs);
}
