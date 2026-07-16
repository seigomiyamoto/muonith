/// @file exemdl_build_geometry.cpp
/// @brief Implementation of geometry construction module.
#include "exemdl_build_geometry.hpp"
#include "ns_mylogger.hpp"
#include "ns_dem_operator.hpp"
#include "ns_iodir.hpp"

namespace fs = std::filesystem;

namespace {

/// @brief Output a single Grid3dVoxel z-direction cross-section.
void output_g3vox_zcross(
    const Grid3dVoxel& g3vox,
    double zmin, double zmax, double zstep,
    int n_detector, bool output_binary)
{
  Grid3dVoxel::CrossSectionZParameters prm;
  prm.xmin  = g3vox.get_xmin();
  prm.xmax  = g3vox.get_xmax();
  prm.xstep = g3vox.get_x_interval();
  prm.ymin  = g3vox.get_ymin();
  prm.ymax  = g3vox.get_ymax();
  prm.ystep = g3vox.get_y_interval();
  prm.zmin  = zmin;
  prm.zmax  = zmax;
  prm.zstep = zstep;
  prm.n_detector = n_detector;
  prm.output_binary = output_binary;
  LOG_INFO("{} cross-section: zstep={}", g3vox.get_name(), prm.zstep);
  char buf[256];
  snprintf(buf, sizeof(buf), "%s_zcross_all.tmp", g3vox.get_name().c_str());
  const fs::path pathout = iodir::make_pathout(buf);
  g3vox.out_cross_section_z_all(pathout, prm);
}

/// @brief Output a single Grid2dPillar (shell) z-direction cross-section.
void output_shell_zcross(
    const Grid2dPillar& shell, bool has_shell,
    const Grid3dVoxel& g3vox_ref,
    double zmin, double zmax, double zstep,
    int n_detector, bool output_binary)
{
  if (!has_shell) return;
  Grid2dPillar::CrossSectionZParameters prm;
  prm.xmin  = g3vox_ref.get_xmin();
  prm.xmax  = g3vox_ref.get_xmax();
  prm.xstep = g3vox_ref.get_x_interval();
  prm.ymin  = g3vox_ref.get_ymin();
  prm.ymax  = g3vox_ref.get_ymax();
  prm.ystep = g3vox_ref.get_y_interval();
  prm.zmin  = zmin;
  prm.zmax  = zmax;
  prm.zstep = zstep;
  prm.n_detector = n_detector;
  prm.output_binary = output_binary;
  const std::string name = fmt::format("{}_zcross_all.tmp", shell.get_name());
  const fs::path pathout = iodir::make_pathout(name);
  shell.out_cross_section_z_all(pathout, prm);
}

/// @brief Output z-direction cross-sections for g3vox and shell grids.
void output_zcross_sections(
    const exemdl::build_geometry::BuildResult& ret,
    const exemdl::build_geometry::BuildArgs& args)
{
  const int n_detector = args.app_params.prm_det.get_n_det();
  const double zmin = args.app_params.zcross.min;
  const double zmax = args.app_params.zcross.max;
  const double zstep = (args.app_params.zcross.zstep > 0.0)
                      ? args.app_params.zcross.zstep
                      : ret.g3vox_merged_input.get_z_interval();
  const bool output_binary = args.app_params.zcross.output_binary;

  output_g3vox_zcross(ret.g3vox_input, zmin, zmax, zstep, n_detector, output_binary);
  output_g3vox_zcross(ret.g3vox_merged_input, zmin, zmax, zstep, n_detector, output_binary);

  output_shell_zcross(ret.g2pil_shell_upper,   ret.has_shell_upper,   ret.g3vox_input, zmin, zmax, zstep, n_detector, output_binary);
  output_shell_zcross(ret.g2pil_shell_lower,   ret.has_shell_lower,   ret.g3vox_input, zmin, zmax, zstep, n_detector, output_binary);
  output_shell_zcross(ret.g2pil_shell_lateral, ret.has_shell_lateral, ret.g3vox_input, zmin, zmax, zstep, n_detector, output_binary);
}

} // anonymous namespace

exemdl::build_geometry::BuildResult
  exemdl::build_geometry::build_all(const BuildArgs& args)
{
  BuildResult ret;

  // Build naive cuboid from terrain parameters
  ret.g2pil_naive = Grid2dPillar(args.app_params.prm_g2pil);

  if (args.app_params.prm_g3vox.tf_build_g3vox == false) {
    LOG_INFO("prm_g3vox.tf_build_g3vox=false, skipping g3vox construction.");

    LOG_INFO("Outputting cross-section of g2pil_naive");
    Grid2dPillar::CrossSectionZParameters prm_zcross;
    prm_zcross.xmin          = ret.g2pil_naive.get_xmin();
    prm_zcross.xmax          = ret.g2pil_naive.get_xmax();
    prm_zcross.xstep         = ret.g2pil_naive.get_x_interval();
    prm_zcross.ymin          = ret.g2pil_naive.get_ymin();
    prm_zcross.ymax          = ret.g2pil_naive.get_ymax();
    prm_zcross.ystep         = ret.g2pil_naive.get_y_interval();
    prm_zcross.zmin          = args.app_params.zcross.min;
    prm_zcross.zmax          = args.app_params.zcross.max;
    const double z_interval  = args.app_params.prm_g3vox.z_pitch;
    prm_zcross.zstep         = (args.app_params.zcross.zstep > 0.0)
                              ? args.app_params.zcross.zstep : z_interval;
    LOG_INFO("Z cross-section: zstep={}, z_pitch={}", prm_zcross.zstep, z_interval);
    prm_zcross.n_detector    = args.app_params.prm_det.get_n_det();
    prm_zcross.output_binary = args.app_params.zcross.output_binary;

    LOG_DEBUG(
      args.app_params.zcross.output_binary
        ? "Output cross sections as binary files"
        : "Output cross sections as text files");

    ret.g2pil_naive.set_name("g2pil_naive");
    char buf[256];
    snprintf(buf, sizeof(buf), "%s_zcross_all.tmp", ret.g2pil_naive.get_name().c_str());
    fs::path pathout = iodir::make_pathout(buf);
    ret.g2pil_naive.out_cross_section_z_all(pathout, prm_zcross);

    return ret;
  }

  // Build voxel grid from cuboid with input density
  ret.g3vox_input = Grid3dVoxel(
    ret.g2pil_naive, args.app_params.prm_g3vox, args.app_params.prm_det.get_n_det());
  ret.g3vox_input.set_name("g3vox_input");

  // Merge voxels to lower resolution
  ret.g3vox_merged_input = ret.g3vox_input.merge(args.app_params.prm_g3merge, args.app_params.prm_reconst_voxels);
  ret.g3vox_merged_input.set_name("g3vox_merged_input");

  // Extract upper, lower, and lateral shells from cuboid and merged voxel
  dem_operator::ShellTriples shell = dem_operator::make_shell_triples_from_g2pil_and_g3vox(
    ret.g2pil_naive, ret.g3vox_merged_input);

  ret.g2pil_shell_upper   = std::move(shell.upper);
  ret.g2pil_shell_lower   = std::move(shell.lower);
  ret.g2pil_shell_lateral = std::move(shell.lateral);
  ret.has_shell_upper   = shell.has_upper;
  ret.has_shell_lower   = shell.has_lower;
  ret.has_shell_lateral = shell.has_lateral;

  ret.g2pil_shell_upper.set_name("g2pil_shell_upper");
  ret.g2pil_shell_lower.set_name("g2pil_shell_lower");
  ret.g2pil_shell_lateral.set_name("g2pil_shell_lateral");

  LOG_INFO("Shell construction complete: has_shell_upper={}, has_shell_lower={}, has_shell_lateral={}",
    ret.has_shell_upper, ret.has_shell_lower, ret.has_shell_lateral);

  // Set shell densities before cross-section output.
  // trace_path_lengths/compute_prior will independently copy these and set densities again, so this is safe.
  if (!args.app_params.vec_prm_nagainv.empty()) {
    const auto& prm0 = args.app_params.vec_prm_nagainv.front();
    const double upd = prm0.get_uniform_prior_density();
    if (ret.has_shell_upper)   ret.g2pil_shell_upper.set_uniform_density(prm0.get_shell_density_upper(upd));
    if (ret.has_shell_lower)   ret.g2pil_shell_lower.set_uniform_density(prm0.get_shell_density_lower(upd));
    if (ret.has_shell_lateral) ret.g2pil_shell_lateral.set_uniform_density(prm0.get_shell_density_lateral(upd));
  }

  output_zcross_sections(ret, args);

  return ret;
}
