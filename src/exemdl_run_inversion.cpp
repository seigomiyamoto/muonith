// src/exemdl_run_inversion.cpp
#include "exemdl_run_inversion.hpp"
#include "ns_myapp.hpp"
#include <spdlog/spdlog.h>
#include <chrono>
#include <fstream>
#include "ns_iodir.hpp"

exemdl::run_inversion::InversionResults
  exemdl::run_inversion::build_inv_res(
    const BuildArgs& args
  , const std::string& name_NagaInv
  , const Eigen::VectorXf& vecxf_nmuon_prior_in
  , const Eigen::MatrixXf& mat_grouped_dNdD_prior_in
  , const Grid3dVoxel& g3vox_prior)
{
  std::chrono::system_clock::time_point start, end;

  // Create NagaInvManager instance within this module
  NagaInvManager manager(args.params.vec_prm_nagainv);
  LOG_INFO("Created NagaInvManager internally");

  // Create NagaInv instance
  std::unique_ptr<NagaInv> pNagainv
    = manager.create_nagainv_instance(
        args.index_run
      , name_NagaInv
      , g3vox_prior
      , mat_grouped_dNdD_prior_in
      , args.vecxf_nmuon_obs
      , vecxf_nmuon_prior_in
      );

  // Wiring B: add efficiency-uncertainty variance to the C_N diagonal post-construction.
  // An empty args.vecxf_var_eff (efficiency C_N diagonal disabled) is a no-op.
  pNagainv->add_var_eff_to_cov_muon(args.vecxf_var_eff);

  start = time_now;

  LOG_INFO("=============================================================");
  LOG_INFO("Density reconstruction started");
  LOG_INFO("=============================================================");
  // These parameters are currently hardcoded. They should be loaded from JSON
  // configuration once a unified parameter loading mechanism for these flags
  // is established in the AppParameters structure.
  const bool tf_save_tmp_data_bin = true;
  const bool tf_verbose = true;

  NagaInvParameters prm_nagainv
    = args.params.vec_prm_nagainv.at(args.index_run).clone();

  const Eigen::VectorXf vecxf_dens_prior
    = g3vox_prior.get_vecxf_density();
  
  std::optional<Eigen::VectorXf> vecxf_dens_input;
  if (args.tf_use_dens_input) {
    vecxf_dens_input = args.g3vox_input.get_vecxf_density();
  } else {
    vecxf_dens_input = std::nullopt; // No input density
  }
  
  NagaInv::ReconstResult reconst_res =
    pNagainv->mp_reconst_density_float(
        prm_nagainv
      , vecxf_dens_prior
      , vecxf_dens_input
      , tf_save_tmp_data_bin
      , tf_verbose
    );
  end = time_now;

  // Log reconstruction completion time
  myapp::cast_time_msec(
      spdlog::level::info
    , "mp_reconst_density_float for " + name_NagaInv + " completed", start, end
  );

  // Define return value
  InversionResults invRes;

  // Set cross-section output parameters
  Grid3dVoxel::CrossSectionZParameters prm_zcross;
  prm_zcross.xmin       = args.g3vox_input.get_xmin();
  prm_zcross.xmax       = args.g3vox_input.get_xmax();
  prm_zcross.xstep      = args.g3vox_input.get_x_interval();
  prm_zcross.ymin       = args.g3vox_input.get_ymin();
  prm_zcross.ymax       = args.g3vox_input.get_ymax();
  prm_zcross.ystep      = args.g3vox_input.get_y_interval();
  prm_zcross.zmin       = args.params.zcross.min;
  prm_zcross.zmax       = args.params.zcross.max;
  const double g3vox_z_interval = args.g3vox_input.get_z_interval();
  prm_zcross.zstep      = (args.params.zcross.zstep > 0.0)
                         ? args.params.zcross.zstep : g3vox_z_interval;
  prm_zcross.n_detector = args.params.prm_det.get_n_det();
  prm_zcross.output_binary = args.params.zcross.output_binary;

  LOG_INFO("Z cross-section: zstep={}, z_interval={}",
           prm_zcross.zstep, g3vox_z_interval);

  invRes.prm_zcross = prm_zcross;

  std::string name;
  fs::path pathout;

  // Output cross-sections for prior model
  LOG_INFO("Outputting prior model cross-sections");
  name = fmt::format("g3vox_prior_{}_{:02d}_zcross_all.tmp", name_NagaInv, args.index_run);
  pathout = iodir::make_pathout(name);
  g3vox_prior.out_cross_section_z_all(pathout, prm_zcross);

  std::string prefix;
  LOG_INFO("Setting and outputting reconstructed density to g3vox_rec");
  prefix = fmt::format("g3vox_rec_{}_{:02d}", name_NagaInv, args.index_run);
  invRes.g3vox_rec = args.g3vox_input.write_density_to_cross_section(
    prefix, reconst_res.vecxf_dens_rec, prm_zcross, true, args.density_quad
  );

  LOG_INFO("Setting and outputting prior difference to g3vox_delta_prior");
  prefix = fmt::format("g3vox_delta_prior_{}_{:02d}", name_NagaInv, args.index_run);
  invRes.g3vox_delta_prior = args.g3vox_input.write_density_to_cross_section(
      prefix, reconst_res.vecxf_delta_dens_prior, prm_zcross
    );

  LOG_INFO("Setting and outputting observation difference to g3vox_diff_real");
  prefix = fmt::format("g3vox_diff_real_{}_{:02d}", name_NagaInv, args.index_run);
  invRes.g3vox_diff_real = args.g3vox_input.write_density_to_cross_section(
      prefix, reconst_res.vecxf_diff_from_real, prm_zcross
  );
  LOG_INFO("Setting and outputting uncertainty to g3vox_dens_error");
  prefix = fmt::format("g3vox_dens_error_{}_{:02d}", name_NagaInv, args.index_run);
  invRes.g3vox_dens_err = args.g3vox_input.write_density_to_cross_section(
      prefix, reconst_res.vecxf_diag_sqrt_cov_dens, prm_zcross
  );

  // Self-evaluation: muon chi-square + effective number of parameters of rho'.
  // Gated by tf_calc_chi2ndf (default off): no extra computation when disabled.
  if (prm_nagainv.get_tf_calc_chi2ndf()) {
    // Warn (physics output unchanged) when the observed counts carry no injected
    // measurement noise: for noiseless forward-projected synthetic data, the residual
    // sits inside the assumed C_N, so chi2/ndf is structurally < 1 and can mislead.
    if (!prm_nagainv.get_tf_signal_poisson()) {
      LOG_WARN("Model self-eval ({}): tf_signal_poisson=false, observed counts have no "
               "injected Poisson noise; chi2/ndf can be structurally < 1 for synthetic data.",
               name_NagaInv);
    }

    // Post-fit residual, consistent with the linearized forward model
    // N(rho) = N_prior + A (rho - rho0): Delta = (N_obs - N_prior) - A (rho' - rho0).
    // A (= dN/drho) and C_N come from the same NagaInv used for the fit
    // (C_N already includes the add_var_eff_to_cov_muon update above).
    const int n_obs = pNagainv->get_num_obs();          // ndf = number of observation bins
    invRes.ndf_muon   = static_cast<double>(n_obs);
    invRes.p_eff_muon = reconst_res.p_eff;              // trace(R), reused gain (0 if unavailable)

    const Eigen::VectorXf& vecxf_nmuon_prior = pNagainv->get_vecxf_nmuon_prior_ref();
    const Eigen::MatrixXf& mat_dNdD_A        = pNagainv->get_mat_dNdD_ref();     // A = dN/drho
    const Eigen::MatrixXf& mat_cov_N         = pNagainv->get_mat_cov_muon_ref(); // C_N (diagonal)

    const bool tf_chi2_inputs_ok =
         reconst_res.is_valid
      && n_obs > 0
      && mat_dNdD_A.rows() == n_obs
      && mat_dNdD_A.cols() == reconst_res.vecxf_delta_dens_prior.size()
      && vecxf_nmuon_prior.size() == n_obs
      && args.vecxf_nmuon_obs.size() == n_obs
      && mat_cov_N.rows() == n_obs;

    if (tf_chi2_inputs_ok) {
      // post-fit residual in muon-count space
      const Eigen::VectorXf vecxf_resid_nmuon
        = (args.vecxf_nmuon_obs - vecxf_nmuon_prior)
        - mat_dNdD_A * reconst_res.vecxf_delta_dens_prior;

      // C_N is diagonal (independent bins): chi2 = sum_i resid_i^2 / C_N(i,i)
      const Eigen::ArrayXd arr_var   = mat_cov_N.diagonal().cast<double>().array();
      const Eigen::ArrayXd arr_resid = vecxf_resid_nmuon.cast<double>().array();

      if (arr_var.minCoeff() > 0.0) {
        invRes.chi2_muon = (arr_resid.square() / arr_var).sum();
        // ndf = N_obs. trace(R)/p_eff reporting is disabled for the public release (bd id-ll39vv).
        LOG_INFO("Model self-eval ({}): chi2_muon={:.3f}, ndf={}, chi2/ndf={:.4f}",
                 name_NagaInv, invRes.chi2_muon, n_obs,
                 invRes.chi2_muon / static_cast<double>(n_obs));
      } else {
        LOG_WARN("Model self-eval ({}): non-positive C_N variance; chi2 skipped", name_NagaInv);
      }
    } else {
      LOG_WARN("Model self-eval ({}): input dimension/validity check failed; chi2 skipped (n_obs={})",
               name_NagaInv, n_obs);
    }
  }

  invRes.reconst_res = std::move(reconst_res);
  invRes.pNagainv = std::move(pNagainv);

  return invRes;
}

namespace {
  // Write a per-run CSV summary of the muon chi-square self-evaluation.
  // One row per independently reconstructed model; lower/upper are omitted when
  // they merely mirror center (they would carry identical chi2).
  void write_chi2_selfeval_csv(
      const exemdl::run_inversion::InversionResultsAll& allRes, const int index_run,
      const bool tf_signal_poisson)
  {
    // Write the chi2 self-eval CSV in the current working directory, i.e. next to
    // run_prg.sh (the run dir), not under the per-run output subdir: bd id-ll39vv.
    const fs::path pathout(fmt::format("chi2_selfeval_{:02d}.csv", index_run));
    std::ofstream ofs(pathout);
    if (!ofs) {
      LOG_WARN("write_chi2_selfeval_csv: cannot open {}", pathout.string());
      return;
    }
    // Definition (Nishiyama et al. 2017): chi2 = (N_obs - N(rho'))^T C_N^-1 (N_obs - N(rho')), ndf = N_obs.
    // trace(R)/p_eff columns are disabled for the public release (bd id-ll39vv).
    ofs << "# muon chi-square self-eval; chi2=(N_obs-N(rho'))^T C_N^-1 (N_obs-N(rho')), ndf=N_obs\n";
    if (!tf_signal_poisson) {
      // Self-documenting caveat so the number is not misread when the CSV is read
      // without the run log (bd id-ll39vv).
      ofs << "# WARNING: tf_signal_poisson=false; observed counts carry no injected Poisson "
             "noise, so chi2/ndf is structurally < 1 and is NOT a goodness-of-fit measure.\n";
    }
    ofs << "model,chi2_muon,ndf_muon,chi2_per_ndf\n";
    const auto write_row = [&ofs](const char* model,
                                  const exemdl::run_inversion::InversionResults& inv) {
      const double ndf     = inv.ndf_muon;                          // = N_obs
      const double ratio   = (ndf > 0.0) ? inv.chi2_muon / ndf : 0.0;
      ofs << fmt::format("{},{:.6g},{:.0f},{:.6g}\n",
                         model, inv.chi2_muon, ndf, ratio);
    };
    // Only independently reconstructed models are written; lower/upper are skipped
    // when they merely mirror center (bd id-ll39vv).
    write_row("center", allRes.center);
    if (allRes.has_lower) write_row("lower", allRes.lower);
    if (allRes.has_upper) write_row("upper", allRes.upper);
    LOG_INFO("Wrote chi2 self-eval CSV: {}", pathout.string());
  }
} // namespace

exemdl::run_inversion::InversionResultsAll
  exemdl::run_inversion::build_all( const BuildArgs& args )
{
  InversionResultsAll allRes;

  const bool tf_prior_error = args.params.tf_run_inversion_prior_error;
  const bool tf_calc_chi2ndf =
    args.params.vec_prm_nagainv.at(args.index_run).get_tf_calc_chi2ndf();
  const bool tf_signal_poisson =
    args.params.vec_prm_nagainv.at(args.index_run).get_tf_signal_poisson();

  if (!tf_prior_error) {
    LOG_INFO("RUN_INVERSION.tf_prior_error=false: center prior only");

    Grid3dVoxel g3vox_prior_center = args.g3vox_input; // copy constructor
    g3vox_prior_center.set(args.prior_info_all.center.g3vox);

    allRes.center = build_inv_res(
        args
      , "center"
      , args.prior_info_all.center.vecxf_nmuon
      , args.res_mat.mat_dNdD_grouped_center
      , g3vox_prior_center
      );

    // lower/upper slots mirror center for downstream consumers
    auto clone_center_to = [&](InversionResults& dst) {
      dst.reconst_res = allRes.center.reconst_res;
      dst.prm_zcross  = allRes.center.prm_zcross;
      dst.g3vox_rec   = allRes.center.g3vox_rec;
      dst.g3vox_delta_prior = allRes.center.g3vox_delta_prior;
      dst.g3vox_diff_real   = allRes.center.g3vox_diff_real;
      dst.g3vox_dens_err    = allRes.center.g3vox_dens_err;
      dst.chi2_muon = allRes.center.chi2_muon;
      dst.ndf_muon  = allRes.center.ndf_muon;
      dst.p_eff_muon = allRes.center.p_eff_muon;
      dst.pNagainv.reset();
    };

    clone_center_to(allRes.lower);
    clone_center_to(allRes.upper);
    allRes.has_lower = false;
    allRes.has_upper = false;

    if (tf_calc_chi2ndf) write_chi2_selfeval_csv(allRes, args.index_run, tf_signal_poisson);
    return allRes;
  } else {
    LOG_INFO("RUN_INVERSION.tf_prior_error=true: lower/upper prior error enabled");
  }

  // Prepare three types of prior models
  Grid3dVoxel g3vox_prior_low = args.g3vox_input; // copy constructor
  g3vox_prior_low.set(args.prior_info_all.lower.g3vox);

  Grid3dVoxel g3vox_prior_center = args.g3vox_input; // copy constructor
  g3vox_prior_center.set(args.prior_info_all.center.g3vox);

  Grid3dVoxel g3vox_prior_upper = args.g3vox_input; // copy constructor
  g3vox_prior_upper.set(args.prior_info_all.upper.g3vox);

  // Create NagaInv instances and run density reconstruction
  allRes.lower = build_inv_res(
      args
    , "lower"
    , args.prior_info_all.lower.vecxf_nmuon
    , args.res_mat.mat_dNdD_grouped_lower
    , g3vox_prior_low
    );
  allRes.has_lower = true;

  allRes.center = build_inv_res(
      args
    , "center"
    , args.prior_info_all.center.vecxf_nmuon
    , args.res_mat.mat_dNdD_grouped_center
    , g3vox_prior_center
    );

  allRes.upper = build_inv_res(
      args
    , "upper"
    , args.prior_info_all.upper.vecxf_nmuon
    , args.res_mat.mat_dNdD_grouped_upper
    , g3vox_prior_upper
    );
  allRes.has_upper = true;

  if (tf_calc_chi2ndf) write_chi2_selfeval_csv(allRes, args.index_run, tf_signal_poisson);
  return allRes;
}
