// src/exemdl_load_parameters.cpp
#include "exemdl_load_parameters.hpp"
#include <stdexcept>

exemdl::load_parameters::AppParameters
  exemdl::load_parameters::load_all( const nlohmann::json& js )
{
  exemdl::load_parameters::AppParameters retParams;

  const std::string keyZCross = "Z_CROSS_SECTION";
  if (!js.contains(keyZCross)) {
    LOG_ERROR("Z_CROSS_SECTION section missing in JSON");
    THROW_ERROR("exemdl::load_parameters::load_all: Z_CROSS_SECTION not found in JSON");
  }

  const nlohmann::json& jz = js.at(keyZCross);
  retParams.zcross.min  = jz.at("min").get<double>();
  retParams.zcross.max  = jz.at("max").get<double>();
  retParams.zcross.output_binary = jz.value("output_binary", false);
  if(retParams.zcross.min >= retParams.zcross.max) {
    LOG_ERROR("z_cross_min={} must be less than z_cross_max={}", retParams.zcross.min, retParams.zcross.max);
    THROW_ERROR("exemdl::load_parameters::load_all: z_cross_min >= z_cross_max");
  }

  // Read zstep (optional, default=0 means use z_interval at runtime)
  retParams.zcross.zstep = jz.value("zstep", 0.0);
  if (retParams.zcross.zstep < 0.0) {
    LOG_ERROR("Z_CROSS_SECTION zstep={} must be >= 0", retParams.zcross.zstep);
    THROW_ERROR("exemdl::load_parameters::load_all: zstep < 0");
  }

  retParams.prm_det     = DetectorPanel::ParameterLists(js, "DETECTOR_PARAMETER_LISTS");
  retParams.prm_g2pil   = Grid2dPillar::Parameters(js, "GRID2D_PILLAR_PARAMETERS");
  retParams.prm_bingroup= Grid2dBinGroup::Parameters(js, "BIN_GROUP_PARAMETERS");
  retParams.prm_g3vox   = Grid3dVoxel::Parameters(js, "GRID3D_VOXEL_PARAMETERS");
  retParams.prm_g3merge = retParams.prm_g3vox.prm_merge;
  retParams.prm_reconst_voxels = retParams.prm_g3vox.prm_reconst_voxels;
  retParams.prm_path    = pathcalc::Parameters(js, "PATH_LENGTH_PARAMETERS");
  retParams.prm_noise   = NoiseParameters(js, "NOISE_PARAMETERS");
  retParams.ft_real    = FluxTable(js, "FLUX_RANGE_DATA_TABLE_REAL");
  retParams.ft_prior   = FluxTable(js, "FLUX_RANGE_DATA_TABLE_PRIOR");
  retParams.vec_prm_nagainv = load_nagainv_parameters(js, "NAGAINV_PARAMETERS");

  const std::string keyRunInversion = "RUN_INVERSION";
  if (js.contains(keyRunInversion)) {
    const nlohmann::json& j_inv = js.at(keyRunInversion);
    retParams.tf_run_inversion_prior_error
      = j_inv.value("tf_prior_error", true);
  } else {
    retParams.tf_run_inversion_prior_error = true;
  }

  return retParams;
}

std::vector<NagaInvParameters>
  exemdl::load_parameters::load_nagainv_parameters(
    const nlohmann::json& js, const std::string& section_name)
{
  if(section_name.empty()) {
    LOG_ERROR("Section name is empty");
    THROW_ERROR("exemdl::load_parameters::load_nagainv_parameters: Section name is empty");
  }

  LOG_DEBUG("Loading from section '{}'", section_name);
  if (!js.contains(section_name)) {
    LOG_ERROR("Section '{}' not found in JSON", section_name);
    THROW_ERROR("exemdl::load_parameters::load_nagainv_parameters: Section '" + section_name + "' not found in JSON");
  }

  std::vector<NagaInvParameters> vec_param;
  // NOTE:
  // Each element in the "NAGAINV_PARAMETERS" array is a complete parameter set.
  // Therefore, we pass an empty string to NagaInvParameters(...) to indicate
  // that we are *not* extracting a sub-section from a larger JSON object.
  const std::string sub_section_should_have_no_name = "";
  for (const auto& j : js.at(section_name)) {
    // emplace is a move operation, so it avoids unnecessary copies.
    vec_param.emplace_back(j, sub_section_should_have_no_name);
  }

  if (vec_param.size() == 0) {
    LOG_ERROR("No parameters found in section '{}'", section_name);
    THROW_ERROR("exemdl::load_parameters::load_nagainv_parameters: No parameters found in section '" + section_name + "'");
  }

  return vec_param;
}
