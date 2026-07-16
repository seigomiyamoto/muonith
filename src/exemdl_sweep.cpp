/// @file exemdl_sweep.cpp
/// @brief Implementation of parameter sweep orchestration.

#include "exemdl_sweep.hpp"
#include "exemdl_pipeline.hpp"
#include "exemdl_load_parameters.hpp"
#include "ns_iodir.hpp"
#include "ns_myapp.hpp"
#include "ns_mymacro.hpp"
#include "spdlog_pch.hpp"

#include <fstream>
#include <algorithm>
#include <cmath>
#include <cctype>
#include <unordered_map>

namespace exemdl::sweep {

//----------------------------------------------------------------------
SweepConfig SweepConfig::from_json(
    const nlohmann::json& root, const std::string& section)
{
  SweepConfig cfg;
  if (!root.contains(section)) return cfg;  // tf_exec = false
  const auto& js = root.at(section);

  cfg.tf_exec = js.value("tf_exec", false);
  if (!cfg.tf_exec) return cfg;

  cfg.base_index = js.value("base_index", 0);

  if (!js.contains("vec_sigma_rho") || !js.contains("vec_corr_length")) {
    THROW_ERROR("exemdl::sweep::from_json: vec_sigma_rho and vec_corr_length are required");
  }
  cfg.vec_sigma_rho = js.at("vec_sigma_rho").get<std::vector<double>>();
  cfg.vec_corr_length = js.at("vec_corr_length").get<std::vector<double>>();

  if (js.contains("vec_uniform_prior_density")) {
    const auto& arr = js.at("vec_uniform_prior_density");
    for (const auto& elem : arr) {
      if (elem.is_number()) {
        // Legacy format: scalar -> expand to quad
        double v = elem.get<double>();
        cfg.vec_density_quad.push_back({v, v, v, v});
      } else if (elem.is_array() && elem.size() == 4) {
        // New format: [prior, upper, lower, lateral]
        cfg.vec_density_quad.push_back({
          elem.at(0).get<double>(),
          elem.at(1).get<double>(),
          elem.at(2).get<double>(),
          elem.at(3).get<double>()
        });
      } else {
        THROW_ERROR("exemdl::sweep::from_json: "
          "vec_uniform_prior_density element must be a number "
          "or array of 4 numbers. got={}", elem.dump());
      }
    }
  }

  cfg.link_diag_to_sigma = js.value("link_diag_to_sigma", true);
  if (!cfg.link_diag_to_sigma) {
    if (!js.contains("vec_sigma_rho_diag")) {
      THROW_ERROR("exemdl::sweep::from_json: vec_sigma_rho_diag required when link_diag_to_sigma=false");
    }
    cfg.vec_sigma_rho_diag = js.at("vec_sigma_rho_diag").get<std::vector<double>>();
  }

  cfg.module8_mode = js.value("module8_mode", std::string("none"));
  cfg.module8_indices = js.value("module8_indices", std::vector<int>{});

  if (cfg.vec_sigma_rho.empty()) {
    THROW_ERROR("exemdl::sweep::from_json: vec_sigma_rho is empty");
  }
  if (cfg.vec_corr_length.empty()) {
    THROW_ERROR("exemdl::sweep::from_json: vec_corr_length is empty");
  }

  return cfg;
}

//----------------------------------------------------------------------
int SweepConfig::total_combinations() const
{
  int n_upd = std::max(static_cast<int>(vec_density_quad.size()), 1);
  int n_sr  = static_cast<int>(vec_sigma_rho.size());
  int n_cl  = static_cast<int>(vec_corr_length.size());
  int n_srd = link_diag_to_sigma
    ? n_sr  // sigma_rho_diag follows sigma_rho, not independent axis
    : std::max(static_cast<int>(vec_sigma_rho_diag.size()), 1);

  if (link_diag_to_sigma) {
    // sigma_rho and sigma_rho_diag are linked: n_sr * n_cl * n_upd
    return n_sr * n_cl * n_upd;
  } else {
    return n_sr * n_cl * n_srd * n_upd;
  }
}

//----------------------------------------------------------------------
std::vector<SweepPoint> SweepConfig::generate_points() const
{
  std::vector<SweepPoint> points;

  // If vec_density_quad is empty, use sentinel {-1,-1,-1,-1} (will use base values)
  auto dq_values = vec_density_quad;
  if (dq_values.empty()) dq_values = {{-1.0, -1.0, -1.0, -1.0}};

  int idx = 0;
  for (const auto& dq : dq_values) {
    for (double sr : vec_sigma_rho) {
      std::vector<double> srd_values;
      if (link_diag_to_sigma) {
        srd_values = {sr};  // sigma_rho_diag = sigma_rho
      } else {
        srd_values = vec_sigma_rho_diag;
      }
      for (double cl : vec_corr_length) {
        for (double srd : srd_values) {
          SweepPoint pt;
          pt.sigma_rho = sr;
          pt.corr_length = cl;
          pt.sigma_rho_diag = srd;
          pt.density_quad = dq;
          pt.index = idx++;
          pt.label = fmt::format("sr{}_cl{}_srd{}_p{}_u{}_l{}_lat{}",
              static_cast<int>(sr), static_cast<int>(cl),
              static_cast<int>(srd),
              static_cast<int>(dq[0]), static_cast<int>(dq[1]),
              static_cast<int>(dq[2]), static_cast<int>(dq[3]));
          points.push_back(pt);
        }
      }
    }
  }
  return points;
}

//----------------------------------------------------------------------
bool SweepConfig::should_run_module8(int sweep_index, int total_count) const
{
  if (module8_mode == "none") return false;
  if (module8_mode == "all") return true;
  if (module8_mode == "last") return (sweep_index == total_count - 1);
  if (module8_mode == "selected") {
    return std::find(module8_indices.begin(), module8_indices.end(), sweep_index)
           != module8_indices.end();
  }
  LOG_WARN("exemdl::sweep: Unknown module8_mode '{}', treating as 'none'", module8_mode);
  return false;
}

//----------------------------------------------------------------------
int resolve_end_stage(const std::string& token)
{
  // Accept a plain non-negative integer (all digits) as-is.
  const bool all_digit = !token.empty() &&
      std::all_of(token.begin(), token.end(),
                  [](unsigned char c) { return std::isdigit(c) != 0; });
  if (all_digit) return std::stoi(token);

  // Otherwise map a feature-name token to its stage number.
  static const std::unordered_map<std::string, int> stage_by_name = {
    {"build_geometry", 3},
    {"trace_path_lengths", 4},
    {"compute_prior", 5},
    {"build_observation_matrix", 6},
    {"invert_density", 7},
    {"analyze_errors", 8},
  };
  const auto it = stage_by_name.find(token);
  if (it == stage_by_name.end()) {
    THROW_ERROR("resolve_end_stage: unknown end_stage token '{}'. "
      "Expected an integer 3-8 or one of: build_geometry, trace_path_lengths, "
      "compute_prior, build_observation_matrix, invert_density, analyze_errors.", token);
  }
  return it->second;
}

//----------------------------------------------------------------------
std::vector<SweepResult> run_sweep(
    const fs::path& json_path, bool seed_given, unsigned seed_value,
    const fs::path& resume_from, int end_stage)
{
  using namespace exemdl::pipeline;

  // --- Init or resume from checkpoint ---
  State state;
  fs::path original_outdir;

  if (!resume_from.empty() && fs::exists(resume_from)) {
    LOG_INFO("exemdl::sweep: Resuming from checkpoint '{}'", resume_from.string());
    state = load(resume_from);

    // Reload app_params from the current JSON (allows config changes on resume)
    state.js = myapp::load_json(json_path);
    state.app_params = exemdl::load_parameters::load_all(state.js);

    original_outdir = iodir::get_default_output_dir();
    LOG_INFO("exemdl::sweep: Resumed (completed_module={}, has_det={}, has_mat={})",
      state.completed_module, state.det.has_value(), state.mat.has_value());
  } else {
    // Initialize pipeline (init/load_parameters/build_geometry)
    state = init(json_path, seed_given, seed_value);
    original_outdir = iodir::get_default_output_dir();

    // Save checkpoint after init/load_parameters/build_geometry
    const fs::path checkpoint_build_geometry = iodir::make_pathout("checkpoint_build_geometry");
    save(state, checkpoint_build_geometry);
    iodir::set_default_output_dir(original_outdir);
    LOG_INFO("exemdl::sweep: Checkpoint saved after init/load_parameters/build_geometry (Module 1-3) to '{}'", checkpoint_build_geometry.string());
  }

  // Resolve end_stage: CLI > JSON5 > default(8)
  if (end_stage < 0) {
    // JSON5 end_stage may be an integer or a feature-name string.
    if (state.js.contains("end_stage") && state.js.at("end_stage").is_string()) {
      end_stage = resolve_end_stage(state.js.at("end_stage").get<std::string>());
    } else {
      end_stage = state.js.value("end_stage", 8);
    }
  }

  // Validate end_stage
  if (end_stage < 3) {
    THROW_ERROR("run_sweep: end_stage={} is invalid. "
      "Modules 1-3 are bundled in init() and cannot be split. "
      "Minimum allowed value is 3.", end_stage);
  }
  if (end_stage != 3 && end_stage != 4 && end_stage != 5 && end_stage != 6 && end_stage != 7 && end_stage != 8) {
    THROW_ERROR("run_sweep: end_stage={} is out of range. "
      "Valid values are 3, 4, 5, 6, 7, or 8.", end_stage);
  }

  // Early return if pipeline stops at build_geometry
  if (end_stage == 3) {
    LOG_INFO("run_sweep: Pipeline completed up to Build Geometry (Module 3) (end_stage=3). Modules 4-8 were not executed.");
    return {};
  }

  // Early return if pipeline stops at trace_path_lengths
  if (end_stage == 4) {
    run_trace_path_lengths(state);
    const fs::path checkpoint_trace_path_lengths = iodir::make_pathout("checkpoint_trace_path_lengths");
    save(state, checkpoint_trace_path_lengths);
    iodir::set_default_output_dir(original_outdir);
    LOG_INFO("run_sweep: Pipeline completed up to Trace Path Lengths (Module 4) (end_stage=4). "
      "Checkpoint saved to '{}'", checkpoint_trace_path_lengths.string());
    return {};
  }

  // Parse sweep config
  SweepConfig cfg = SweepConfig::from_json(state.js);
  if (cfg.tf_exec && end_stage < 7) {
    THROW_ERROR("run_sweep: end_stage={} is incompatible with sweep mode (tf_exec=true). "
      "Sweep requires at least Invert Density (Module 7). "
      "To run modules 1-{} only, set NAGAINV_PARAM_SWEEP.tf_exec=false.",
      end_stage, end_stage);
  }
  if (!cfg.tf_exec) {
    LOG_INFO("exemdl::sweep: NAGAINV_PARAM_SWEEP not active, running single configuration");
    const auto& base_prm = state.app_params.vec_prm_nagainv.at(cfg.base_index);
    double upd = base_prm.get_uniform_prior_density();
    double sr  = base_prm.get_sigma_rho();
    double cl  = base_prm.get_corr_length();
    double srd = base_prm.get_sigma_rho_diag();
    std::array<double,4> base_dq = {upd,
      base_prm.get_shell_density_upper(upd),
      base_prm.get_shell_density_lower(upd),
      base_prm.get_shell_density_lateral(upd)};

    // Run raytrace (density-independent)
    if (!state.det.has_value()) {
      run_trace_path_lengths(state);
      const fs::path checkpoint_trace_path_lengths = iodir::make_pathout("checkpoint_trace_path_lengths");
      save(state, checkpoint_trace_path_lengths);
      iodir::set_default_output_dir(original_outdir);
      LOG_INFO("exemdl::sweep: checkpoint_trace_path_lengths saved to '{}'", checkpoint_trace_path_lengths.string());
    } else {
      LOG_INFO("exemdl::sweep: Raytrace already completed (resumed from checkpoint)");
    }

    run_compute_prior(state, base_dq);
    run_build_observation_matrix(state);

    const fs::path checkpoint_build_observation_matrix = iodir::make_pathout("checkpoint_build_observation_matrix");
    save(state, checkpoint_build_observation_matrix);
    iodir::set_default_output_dir(original_outdir);
    LOG_INFO("exemdl::sweep: checkpoint_build_observation_matrix saved to '{}'", checkpoint_build_observation_matrix.string());

    if (end_stage <= 6) {
      LOG_INFO("run_sweep: Pipeline completed up to Module {} (end_stage={}). "
        "Modules {}-8 were not executed.", end_stage, end_stage, end_stage + 1);
      SweepResult res;
      res.point = {sr, cl, srd, base_dq, 0, "base"};
      res.module8_executed = false;
      return {res};
    }
    auto inv_result = run_invert_density(state, sr, cl, srd);
    save_recon_io(state, inv_result, sr, cl, srd);
    if (end_stage <= 7) {
      LOG_INFO("run_sweep: Pipeline completed up to Module {} (end_stage={}). "
        "Modules {}-8 were not executed.", end_stage, end_stage, end_stage + 1);
      SweepResult res;
      res.point = {sr, cl, srd, base_dq, 0, "base"};
      res.module8_executed = false;
      return {res};
    }

    SweepResult res;
    res.point = {sr, cl, srd, base_dq, 0, "base"};
    res.module8_executed = true;
    run_analyze_errors(state, inv_result);

    return {res};
  }

  // Validate base_index
  if (cfg.base_index < 0 ||
      cfg.base_index >= static_cast<int>(state.app_params.vec_prm_nagainv.size())) {
    THROW_ERROR("exemdl::sweep: base_index {} out of range (size={})",
                cfg.base_index, state.app_params.vec_prm_nagainv.size());
  }

  auto points = cfg.generate_points();
  const int total = static_cast<int>(points.size());
  LOG_INFO("exemdl::sweep: {} combinations to sweep", total);

  // Run raytrace once — unless already done (resume from checkpoint with det)
  if (!state.det.has_value()) {
    run_trace_path_lengths(state);
    const fs::path checkpoint_trace_path_lengths = original_outdir / "checkpoint_trace_path_lengths";
    save(state, checkpoint_trace_path_lengths);
    iodir::set_default_output_dir(original_outdir);
    LOG_INFO("exemdl::sweep: checkpoint_trace_path_lengths saved (raytrace done once)");
  } else {
    LOG_INFO("exemdl::sweep: Raytrace already completed (resumed from checkpoint)");
  }

  std::vector<SweepResult> results;
  results.reserve(total);

  // Track current density_quad to avoid unnecessary prior(M5)+M6 re-runs
  std::array<double,4> current_dq = {
    std::numeric_limits<double>::quiet_NaN(),
    std::numeric_limits<double>::quiet_NaN(),
    std::numeric_limits<double>::quiet_NaN(),
    std::numeric_limits<double>::quiet_NaN()};
  const auto& base_prm = state.app_params.vec_prm_nagainv.at(cfg.base_index);

  // Warn once if end_stage skips analyze_errors but module8_mode is set
  if (end_stage <= 7 && cfg.module8_mode != "none") {
    LOG_WARN("run_sweep: end_stage={} skips Analyze Errors (Module 8), "
      "but module8_mode='{}' is set. Analyze Errors (Module 8) will not run.",
      end_stage, cfg.module8_mode);
  }

  for (const auto& pt : points) {
    auto dq = pt.density_quad;
    // Sentinel -1 means use base parameter values
    if (dq[0] < 0) {
      double upd = base_prm.get_uniform_prior_density();
      dq = {upd,
        base_prm.get_shell_density_upper(upd),
        base_prm.get_shell_density_lower(upd),
        base_prm.get_shell_density_lateral(upd)};
    }

    // Skip sweep points whose output directory already exists
    const std::string sweep_dir = fmt::format("sweep_{:03d}_{}", pt.index, pt.label);
    const fs::path sweep_path = original_outdir / sweep_dir;
    if (fs::exists(sweep_path) && !fs::is_empty(sweep_path)) {
      LOG_INFO("exemdl::sweep: [{}/{}] Skipping (already exists: {})",
               pt.index + 1, total, sweep_dir);
      SweepResult sr;
      sr.point = pt;
      sr.module8_executed = false;  // unknown, mark as false
      results.push_back(sr);
      // Still track dq so checkpoint_m6 reuse works for subsequent points
      current_dq = dq;
      continue;
    }

    // Re-run prior(M5)+M6 only when density_quad changes (raytrace already done)
    if (dq != current_dq) {
      iodir::set_default_output_dir(original_outdir);
      const std::string dq_suffix = fmt::format("p{}_u{}_l{}_lat{}",
        static_cast<int>(dq[0]), static_cast<int>(dq[1]),
        static_cast<int>(dq[2]), static_cast<int>(dq[3]));
      const fs::path checkpoint_build_observation_matrix =
        original_outdir / ("checkpoint_build_observation_matrix_" + dq_suffix);
      if (fs::exists(checkpoint_build_observation_matrix / "state_meta.bin")) {
        // Reuse existing checkpoint
        LOG_INFO("exemdl::sweep: Loading checkpoint {}", dq_suffix);
        auto loaded = load(checkpoint_build_observation_matrix);
        state.det = std::move(loaded.det);
        state.mat = std::move(loaded.mat);
        state.completed_module = loaded.completed_module;
        iodir::set_default_output_dir(original_outdir);
      } else {
        // Run prior(M5) + M6 only (raytrace already done)
        LOG_INFO("exemdl::sweep: density_quad=[{},{},{},{}], running prior(M5)+M6",
          dq[0], dq[1], dq[2], dq[3]);
        run_compute_prior(state, dq);
        run_build_observation_matrix(state);
        save(state, checkpoint_build_observation_matrix);
        iodir::set_default_output_dir(original_outdir);
        LOG_INFO("exemdl::sweep: Checkpoint saved after prior(M5)+M6 to '{}'",
          checkpoint_build_observation_matrix.string());
      }
      current_dq = dq;
    }

    // Set output subdirectory for this sweep point
    iodir::set_default_output_dir(original_outdir / sweep_dir);

    LOG_INFO("exemdl::sweep: [{}/{}] sigma_rho={}, corr_length={}, sigma_rho_diag={}, dq=[{},{},{},{}]",
             pt.index + 1, total, pt.sigma_rho, pt.corr_length, pt.sigma_rho_diag,
             dq[0], dq[1], dq[2], dq[3]);

    auto inv_result = run_invert_density(state, pt.sigma_rho, pt.corr_length, pt.sigma_rho_diag);
    save_recon_io(state, inv_result, pt.sigma_rho, pt.corr_length, pt.sigma_rho_diag);

    SweepResult sr;
    sr.point = pt;

    if (end_stage >= 8 && cfg.should_run_module8(pt.index, total)) {
      const std::string prefix = fmt::format("{:03d}", pt.index);
      run_analyze_errors(state, inv_result, prefix);
      sr.module8_executed = true;
    }

    results.push_back(sr);
  }

  // Restore original output directory
  iodir::set_default_output_dir(original_outdir);

  // Write summary CSV (only when full pipeline completed)
  if (end_stage >= 8) {
    const fs::path csv_path = iodir::make_pathout("sweep_summary.csv");
    write_sweep_summary_csv(results, csv_path);
    LOG_INFO("exemdl::sweep: Summary written to {}", csv_path.string());
  }
  if (end_stage < 8) {
    LOG_INFO("run_sweep: Sweep completed up to Module {} (end_stage={}). "
      "Modules {}-8 were not executed. Sweep summary CSV was not written.",
      end_stage, end_stage, end_stage + 1);
  }

  return results;
}

//----------------------------------------------------------------------
void write_sweep_summary_csv(
    const std::vector<SweepResult>& results, const fs::path& output_path)
{
  std::ofstream ofs(output_path);
  if (!ofs.is_open()) {
    THROW_ERROR("exemdl::sweep::write_sweep_summary_csv: Cannot open {}", output_path.string());
  }

  ofs << "index,sigma_rho,corr_length,sigma_rho_diag,"
      << "prior_density,shell_density_upper,shell_density_lower,shell_density_lateral,"
      << "label,module8_executed\n";

  for (const auto& r : results) {
    ofs << r.point.index << ","
        << r.point.sigma_rho << ","
        << r.point.corr_length << ","
        << r.point.sigma_rho_diag << ","
        << r.point.density_quad[0] << ","
        << r.point.density_quad[1] << ","
        << r.point.density_quad[2] << ","
        << r.point.density_quad[3] << ","
        << r.point.label << ","
        << (r.module8_executed ? "true" : "false") << "\n";
  }
}

} // namespace exemdl::sweep
