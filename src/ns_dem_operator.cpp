// dem_operator.cpp
#include <stdexcept>
#include <sstream>
#include <thread>
#include "ns_mymacro.hpp"
#include "ns_mylogger.hpp"
#include "ns_dem_operator.hpp"
#include "spdlog_pch.hpp"

// Grid2dPillar dem_operator::make_shell_from_g2pil_and_g3vox
Grid2dPillar dem_operator::make_shell_from_g2pil_and_g3vox(
    const Grid2dPillar &g2pil, const Grid3dVoxel &g3vox)
{
  // Create a copy of g2pil
  Grid2dPillar g2pil_shell(g2pil);

  // Get the number of bins in x and y
  const int nbinx = g3vox.get_nbinx();
  const int nbiny = g3vox.get_nbiny();

  // Loop over ix, iy of the merged g3vox
  for(int iy=0;iy<nbiny;iy++){
    for(int ix=0;ix<nbinx;ix++){
      fprintf(stderr,"dem_operator::make_shell_from_g2pil_and_g3vox ix=%04d, iy=%04d\r",ix,iy);

      // Get the highest z_upper where tf_exist==true.
      // If no voxel with tf_exist==true exists in the column formed by ix,iy, returns zmin.
      const double z_highest = g3vox.get_highest_exist_z(ix,iy);

      // Get Grid2d from g3vox.
      Grid2d g2d_merged(g3vox.get_x_axis(), g3vox.get_y_axis());

      // Get the vector of (ix,iy) indices from g2pil_shell that correspond to this merged cell.
      const std::vector<Grid2d::Ixiy> vec_ixiy
       = g2pil_shell.get_non_merged_vec_ixiy(g2d_merged,ix,iy);

      // Loop over ix,iy in vec_ixiy
      for(const auto [ix2,iy2] : vec_ixiy){
        // Call Pillar from g2pil_shell.
        Pillar& cub = g2pil_shell.callPillar(ix2, iy2);

        // Update the cuboid's z_lower.
        cub.set_zmin(z_highest);
      } // End of loop over ix,iy in vec_ixiy

    } // End of ix loop
  } // End of iy loop
  
  return g2pil_shell;
}

// Creates upper, lower, and lateral shell Grid2dPillars
dem_operator::ShellTriples dem_operator::make_shell_triples_from_g2pil_and_g3vox(
    const Grid2dPillar &g2pil, const Grid3dVoxel &g3vox)
{
  ShellTriples result;
  result.upper = Grid2dPillar(g2pil);    // Copy for upper shell
  result.lower = Grid2dPillar(g2pil);    // Copy for lower shell
  result.lateral = Grid2dPillar(g2pil);  // Copy for lateral shell
  result.has_upper = false;
  result.has_lower = false;
  result.has_lateral = false;

  const int nbinx = g3vox.get_nbinx();
  const int nbiny = g3vox.get_nbiny();
  const double g3vox_zmin = g3vox.get_zmin();

  // Pre-construct Grid2d from g3vox axes (shared read-only across threads)
  const Grid2d g2d_merged(g3vox.get_x_axis(), g3vox.get_y_axis());

  int tf_has_upper = 0;
  int tf_has_lower = 0;
  int tf_has_lateral = 0;

  // Each merged (ix, iy) maps to a disjoint set of non-merged (ix2, iy2),
  // so different threads write to different cuboids — no data race.
  #pragma omp parallel for collapse(2) schedule(static) reduction(|:tf_has_upper, tf_has_lower, tf_has_lateral)
  for(int iy=0; iy<nbiny; iy++){
    for(int ix=0; ix<nbinx; ix++){
      // Get z bounds from voxel grid for this column
      const double z_highest = g3vox.get_highest_exist_z(ix,iy);
      const double z_lowest = g3vox.get_lowest_exist_z(ix,iy);
      // Per-column lower boundary: use z_lowest if column has existing voxels,
      // otherwise fall back to g3vox_zmin to avoid overlap with upper shell.
      const double lower_bound_z = (z_lowest <= z_highest) ? z_lowest : g3vox_zmin;

      const std::vector<Grid2d::Ixiy> vec_ixiy
        = result.upper.get_non_merged_vec_ixiy(g2d_merged, ix, iy);

      for(const auto [ix2, iy2] : vec_ixiy){
        // Get original terrain bounds from input g2pil
        const Pillar& cub_orig = g2pil.getPillar(ix2, iy2);
        const double terrain_zmin = cub_orig.get_zmin();  // Base elevation
        const double terrain_zmax = cub_orig.get_zmax();  // Surface elevation

        Pillar& cub_upper = result.upper.callPillar(ix2, iy2);
        Pillar& cub_lower = result.lower.callPillar(ix2, iy2);
        Pillar& cub_lateral = result.lateral.callPillar(ix2, iy2);

        // Detect lateral column: no existing voxels in this column
        const bool is_lateral_column = (z_lowest > z_highest);

        if(is_lateral_column){
          // === Lateral Shell ===
          // Covers full terrain height for columns with all tf_exist=false
          if(terrain_zmax > terrain_zmin){
            cub_lateral.set_zmin(terrain_zmin);
            cub_lateral.set_zmax(terrain_zmax);
            tf_has_lateral = 1;
          } else {
            cub_lateral.set_zmin(terrain_zmin);
            cub_lateral.set_zmax(terrain_zmin);
          }

          // Upper and lower are empty for lateral columns
          cub_upper.set_zmin(terrain_zmax);
          cub_upper.set_zmax(terrain_zmax);
          cub_lower.set_zmin(terrain_zmin);
          cub_lower.set_zmax(terrain_zmin);
        } else {
          // === Normal column with existing voxels ===
          // Lateral is empty
          cub_lateral.set_zmin(terrain_zmin);
          cub_lateral.set_zmax(terrain_zmin);

          // === Upper Shell Logic ===
          // Upper shell exists if terrain surface > highest existing voxel
          if(terrain_zmax > z_highest && z_highest >= terrain_zmin){
            cub_upper.set_zmin(z_highest);
            // zmax remains as terrain_zmax (already set from copy)
            tf_has_upper = 1;
          } else {
            // Mark as empty (zero height) - set both to same value
            cub_upper.set_zmin(terrain_zmax);
            cub_upper.set_zmax(terrain_zmax);
          }

          // === Lower Shell Logic (using per-column lower_bound_z) ===
          if(terrain_zmin < lower_bound_z){
            // Lower shell covers from terrain base up to either:
            // - lower_bound_z (per-column lowest existing voxel), or
            // - terrain surface (if terrain is entirely below lower_bound_z)
            const double lower_shell_zmax = std::min(terrain_zmax, lower_bound_z);
            if(lower_shell_zmax > terrain_zmin){
              cub_lower.set_zmin(terrain_zmin);
              cub_lower.set_zmax(lower_shell_zmax);
              tf_has_lower = 1;
            } else {
              // Zero height - mark as empty
              cub_lower.set_zmin(terrain_zmin);
              cub_lower.set_zmax(terrain_zmin);
            }
          } else {
            // Terrain base is at or above lower_bound_z - no lower shell needed
            cub_lower.set_zmin(lower_bound_z);
            cub_lower.set_zmax(lower_bound_z);
          }
        }
      } // End of loop over vec_ixiy
    } // End of ix loop
  } // End of iy loop

  result.has_upper = (tf_has_upper != 0);
  result.has_lower = (tf_has_lower != 0);
  result.has_lateral = (tf_has_lateral != 0);

  LOG_INFO("make_shell_triples_from_g2pil_and_g3vox: has_upper={}, has_lower={}, has_lateral={}",
    result.has_upper, result.has_lower, result.has_lateral);

  return result;
}

// Calculate the average density of shell + g3vox.
double dem_operator::calc_avr_dens(
  const Grid2dPillar &shell, const Grid3dVoxel &g3vox)
{
  double volume_sum_shell = 0.0;
  double mass_sum_shell = 0.0;

  // for shell
  const int nx_shell = shell.get_nbinx();
  const int ny_shell = shell.get_nbiny();
  const double dx_shell = shell.get_x_interval();
  const double dy_shell = shell.get_y_interval();
  #pragma omp parallel for collapse(2) schedule(static) reduction(+:volume_sum_shell, mass_sum_shell)
  for(int iy=0;iy<ny_shell;iy++){
    for(int ix=0;ix<nx_shell;ix++){
      const Pillar& cub = shell.getPillar(ix,iy);
      const double zmin = cub.get_zmin();
      const double zmax = cub.get_zmax();
      const double vol = dx_shell*dy_shell*(zmax-zmin);
      const double dens = cub.get_density();
      volume_sum_shell += vol;
      mass_sum_shell += vol*dens;
    }
  }
  const double dens_avr_shell = mass_sum_shell / volume_sum_shell;
  LOG_INFO(
    "volume_sum_shell={:E}, mass_sum_shell={:E}, dens_avr_shell={:E}"
    , volume_sum_shell,mass_sum_shell,dens_avr_shell);

  // for g3vox
  double volume_sum_vox = 0.0;
  double mass_sum_vox = 0.0;

  const int nx_vox = g3vox.get_nbinx();
  const int ny_vox = g3vox.get_nbiny();
  const int nz_vox = g3vox.get_nbinz();
  const double dx_vox = g3vox.get_x_interval();
  const double dy_vox = g3vox.get_y_interval();
  const double dz_vox = g3vox.get_z_interval();
  const double vol = dx_vox*dy_vox*dz_vox;
  #pragma omp parallel for collapse(3) reduction(+:volume_sum_vox, mass_sum_vox)
  for(int iz=0;iz<nz_vox;iz++){
    for(int iy=0;iy<ny_vox;iy++){
      for(int ix=0;ix<nx_vox;ix++){
        const Voxel& vox = g3vox.getVoxel(ix,iy,iz);
        if(!vox.get_tf_exist()) continue;
        const double dens = vox.get_density();
        volume_sum_vox += vol;
        mass_sum_vox += vol * dens;
      }
    }
  }
  const double dens_avr_vox = mass_sum_vox / volume_sum_vox;
  LOG_INFO(
    "volume_sum_vox={:E}, mass_sum_vox={:E}, dens_avr_vox={:E}"
    , volume_sum_vox,mass_sum_vox,dens_avr_vox);

  const double vol_sum = volume_sum_shell + volume_sum_vox;
  const double mass_sum = mass_sum_shell + mass_sum_vox;
  const double dens_avr = mass_sum / vol_sum;
  LOG_INFO(
    "vol_sum={:E}, mass_sum={:E}, dens_avr={:E}"
    , vol_sum, mass_sum, dens_avr);
  return dens_avr;
}