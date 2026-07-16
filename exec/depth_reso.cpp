//
// 2025-04-25 15:57:56
//

// for intel thread building block
// https://stackoverflow.com/questions/50136667/using-intel-tbb-in-debug-mode
// #define TBB_USE_DEBUG 1

#include <chrono>
#include <string>
#include <Eigen/Dense>
#include <nlohmann/json.hpp> // for nlohmann::json
using json = nlohmann::json;

#include "cls_Ray.hpp"
#include "ns_myapp.hpp"
#include "cls_DetectorPanelArray.hpp"
#include "cls_DepthResolutionSweep.hpp"
#include "cls_Grid2dPillar.hpp"
#include "cls_Grid2dPillarParameters.hpp"
#include "cls_Grid3dVoxel.hpp"
#include "ns_pathcalc.hpp"
#include "cls_NagaInv.hpp"
#include "ns_mymacro.hpp"
#include "ns_mylogger.hpp"
#include "ns_dem_operator.hpp"
#include "cls_MatrixBuildParameters.hpp"
#include "cls_FluxTable.hpp"
#include "cls_NagaInvManager.hpp"
#include "cls_NoiseParameters.hpp"
#include "st_NagaInvLooper.hpp"

#include "ns_geom_util.hpp"

#include "spdlog_pch.hpp"
#include "ns_mylogger.hpp"
#include "ns_eigen_blas.hpp"
#include "ns_detector_indexing.hpp"

#include <getopt.h>
#include <cstdlib>
#include "ns_seed.hpp"
#include "ns_param_constants.hpp"

static std::string fname_json;
static bool seed_given = false;
static unsigned seed_value = 42;

int do_depth_reso(const std::string& fname_json)
{
  START_FUNC();

  // Load the JSON.
  const json js = myapp::load_json(fname_json);

  // Decide the seed (CLI > JSON > default).
  if (seed_given) {
    seed::set_global_seed(seed_value);
    LOG_INFO("Seed from CLI: {}", seed_value);
    SLEEP_MSEC(500);
  } else if (js.contains("seed")) {
    unsigned json_seed = js["seed"].get<unsigned>();
    seed::set_global_seed(json_seed);
    LOG_INFO("Seed from JSON: {}", json_seed);
    SLEEP_MSEC(500);
  } else {
    seed::set_global_seed(42);  // fallback
    LOG_WARN("No seed provided. Using default seed: 42");
    SLEEP_MSEC(500);
  }

  // Initialize parameter constants from JSON
  param_constants::init(js);

  START_FUNC();
  std::chrono::system_clock::time_point start_main, end_main, start, end;
  start_main = time_now;

  // * Check the DEBUG mode.
  #ifdef NODEBUG
  LOG_DEBUG("DEBUG mode is not active.");
  #else
  LOG_DEBUG("DEBUG mode is active.");
  #endif
  SLEEP_MSEC(500);

  LOG_INFO("Starting application");

hoge(1200000);

  unsigned int n_max_threads = std::thread::hardware_concurrency();
  LOG_INFO("Maximum number of threads on this PC: {}", n_max_threads);

  const int n_threads = omp_get_max_threads();
  LOG_INFO("Maximum number of threads for this program: {}", n_threads);
  SLEEP_MSEC(500);

  // Set the number of threads used by Eigen.
  const int n_threads_eigen = (int)(n_threads/2);
  LOG_INFO("Number of threads for Eigen: {}", n_threads_eigen);
  Eigen::setNbThreads(n_threads);

hoge(2000000);
  myapp::filecheck(fname_json);
hoge(3000000);
  //======================================================================
  // Load the parameter tables from the JSON. The second argument is the section name.
  //======================================================================

  LOG_INFO("Loading parameters from the runcard...");

  // Read sweep_params first to get angle_bin_override settings
  const auto sweep_params = DepthResolutionSweep::Parameters::from_json(js);

  // Build prm_g2det_list (non-const to allow override)
  DetectorPanel::ParameterLists prm_g2det_list(js,"DETECTOR_PARAMETER_LISTS");

  // Build base distribution using original det params (before override)
  const DetectorPanel::ParameterLists prm_g2det_list_orig = prm_g2det_list;

  // Apply angle_bin_override if enabled
  if (sweep_params.tf_override_angle_bin) {
    prm_g2det_list.set_angle_bin_override(
      sweep_params.angle_bin_override.nbinx,
      sweep_params.angle_bin_override.txmin,
      sweep_params.angle_bin_override.txmax,
      sweep_params.angle_bin_override.nbiny,
      sweep_params.angle_bin_override.tymin,
      sweep_params.angle_bin_override.tymax
    );
  }

  const Grid2dPillar::Parameters prm_g2pil(js,"GRID2D_PILLAR_PARAMETERS");
  const Grid2dBinGroup::Parameters prm_g2bg(js,"BINGROUP_PARAMETERS");
  const NoiseParameters prm_noise(js,"NOISE_PARAMETERS");
  const pathcalc::Parameters prm_path(js,"PATH_LENGTH_PARAMETERS");
  const Grid2dBinGroup::Parameters prm_bingroup(js,"BIN_GROUP_PARAMETERS");

  //======================================================================
  // Load peneflux, range and dFdR.
  //======================================================================
  LOG_INFO("Loading the FluxTable.");
  // FluxTable ft_prior( js, "FLUX_RANGE_DATA_TABLE_PRIOR");
  // ft_prior.set_name("fluxtable_prior");
  // ft_prior.load_tables();
  FluxTable ft_real( js, "FLUX_RANGE_DATA_TABLE_REAL");
  ft_real.set_name("fluxtable_real");
  ft_real.load_tables();

hoge(3100000);
  // The constructor below also calls add_density_structure, cylynder, dike, checkerboard
  LOG_INFO("Loading data from the DEM into g2pil.");
  LOG_INFO("A uniform density is given to Grid2dPillar. No additional structure is added yet.");
  Grid2dPillar g2pil0( prm_g2pil.path_dem
                     , prm_g2pil.zmin
                     , prm_g2pil.initial_uniform_density
                     , prm_g2pil.tf_shift_x
                     , prm_g2pil.tf_shift_y
                     , prm_g2pil.tolerance_ratio );
  g2pil0.set_name("g2pil0");
  
hoge(4000000);
  LOG_INFO("Building the template DetectorPanelArray arrdet.");
  DetectorPanelArray arrdet0(prm_g2det_list);
  arrdet0.set_name("arrdet0");
  arrdet0.display_status();

  // Output base PL/signal distribution using original det params
  if (sweep_params.tf_out_det_PL_signal) {
    DetectorPanelArray arrdet_orig(prm_g2det_list_orig);
    for (Detid d = 0; d < arrdet_orig.get_n_det(); ++d) {
      DetectorPanel panel_work = arrdet_orig.getDetectorPanelCopy(d);
      pathcalc::g2pil::mp_add_PLDL(panel_work, g2pil0);
      panel_work.mp_calc_set_peneflux_signal_from_DL(ft_real, false);
      panel_work.out_txtyPL("fig_");
      panel_work.out_txtySignal("fig_");
    }
    spdlog::info("Base PL/signal distribution written to tmp/");
  }

// hoge(8000000);
  //======================================================================
  LOG_INFO("Running DepthResolutionSweep as a trial.");
  //======================================================================

  DepthResolutionSweep sweep( sweep_params, ft_real, g2pil0, arrdet0);
  sweep.run();

hoge(12000000);
  end_main = time_now;
  const std::string msg = "main function done.";
  myapp::cast_time_msec(spdlog::level::info,msg,start_main,end_main);
  LOG_INFO("Application finished");
  // spdlog::shutdown();.....

  END_FUNC();
  return 0;
}

void parse_args(int argc, char** argv)
{
  const struct option longopts[] = {
    {"json",  required_argument, nullptr, 'j'},
    {"seed",  required_argument, nullptr, 's'},
    {nullptr, 0, nullptr, 0}
  };

  int opt;
  while ((opt = getopt_long(argc, argv, "j:s:", longopts, nullptr)) != -1) {
    switch (opt) {
      case 'j': fname_json = optarg; break;
      case 's': seed_given = true; seed_value = std::stoul(optarg); break;
      default:
        std::cerr << "Usage: " << argv[0]
                  << " --json <file> [--seed <val>]"
                  << std::endl;
        std::exit(1);
    }
  }

  if (fname_json.empty()) {
    std::cerr << "Error: --json must be specified\n";
    std::exit(1);
  }
}

//
// main function (mainly handles try/catch)
//
int main(int argc, char **argv)
{
  START_FUNC();

  // Parse the arguments (only the minimum information is needed at this stage).
  parse_args(argc, argv);

  // Load the JSON (for the initial setup in main).
  const json js = myapp::load_json(fname_json);

  if( !js.contains("LOG_FILE") ) {
    THROW_ERROR("JSON does not contain LOG_FILE section.");
  }

  const auto [stdout_level, stderr_level, file_level] = mylogger::get_log_levels(js, "LOG_FILE");
  std::string basename = std::string("main_") + std::filesystem::path(fname_json).stem().string();

  std::string log_dir_str = js.at("LOG_FILE").at("path_log_dir").get<std::string>();


  bool give_new_num = false;
  give_new_num = js.at("LOG_FILE").at("give_new_number").get<bool>();

  if (log_dir_str.empty()) log_dir_str = "logs"; // default log directory
  std::filesystem::path log_dir(log_dir_str);
  std::filesystem::create_directories(log_dir);
  auto log_filepath = mylogger::generate_unique_log_path(log_dir, basename, ".log", give_new_num);
  bool archive_existing = mylogger::get_log_archive_existing(js, "LOG_FILE");

  auto logger = mylogger::create_logger(stdout_level, stderr_level, file_level, log_filepath, archive_existing);
  spdlog::set_default_logger(logger);
  mylogger::configure_rate_limit(js, "LOG_FILE");

  LOG_INFO("Logging started to {}", log_filepath.string());
  LOG_INFO("Default max count for rate-limited logs: {}", mylogger::get_default_max_count());

  // Confirm that the global logger is not nullptr.
  if (mylogger::g_logger == nullptr) {
    std::cerr << "[FATAL] g_logger is nullptr before entering do_depth_reso()." << std::endl;
    std::exit(1);
  }

  // Seed setting: CLI > JSON > default.
  if (seed_given) {
    seed::set_global_seed(seed_value);
    LOG_INFO("Seed from CLI: {}", seed_value);
    SLEEP_MSEC(500);
  } else if (js.contains("seed")) {
    unsigned json_seed = js["seed"];
    seed::set_global_seed(json_seed);
    LOG_INFO("Seed from JSON: {}", json_seed);
    SLEEP_MSEC(500);
  } else {
    seed::set_global_seed(42);
    LOG_WARN("No seed given. Using default seed 42.");
    SLEEP_MSEC(500);
  }

  // Enter the main body (uses global variables).
  try {
    do_depth_reso(fname_json);
    spdlog::default_logger()->flush();
  } catch (...) {
    spdlog::default_logger()->flush();
    LOG_CRITICAL("Unhandled exception");
    return 1;
  }

  END_FUNC();
  LOG_INFO("Seed was {}", seed::get_global_seed());
  spdlog::shutdown();
  return 0;
} // main
