// ${work_dir_name}/swp001/prm_muonith.json5
// muonith.exe configuration for $site_name_cap
{
  "seed" : 42
  // Maximum module to execute (valid: 3, 4, 5, 6, 7, 8). Default: 8 (all modules).
  // CLI --end-stage takes precedence over this value.
  $end_stage_line
  // trace < debug < info < warn < err < critical < off
  // For example, if stdout_level=debug and stderr_level=error,
  // messages at error level or above are output to both stdout and stderr,
  // so if stdout is not redirected to a file, they appear twice on the terminal.
, "LOG_FILE": { // Log file settings
      "path_log_dir": "logs" // Log file directory
    , "give_new_number": false // Whether to assign a new number to log files
    , "log_level": {
        "stdout_level": "debug" // stdout log level
      , "stderr_level": "error" // stderr log level
      , "file_level":   "trace" // file log level
    },
    "default_max_count": 50,
    "archive_existing": false
  }
, "RUN_INVERSION": {
    "tf_prior_error": $tf_prior_error // If true, compute prior error including lower/upper; if false, center only
  }
, "DETECTOR_PARAMETER_LISTS": {
    "name": "arrdet_${det_list_name}" // Name for the detector parameter list. Set freely.
  , "det_files": [
$det_files_block
    ]
  , "tf_apply_eff": false // Whether to apply efficiency to signal
    // false stops the txty ASCII dumps under tmp/. With the current swp001 settings the
    // per-element PL / density / signal figures are plotted from the saved binary,
    // so no figure is lost.
  , "tf_out_txty_ascii": $tf_out_txty_ascii
    // false stops the g2bg ASCII dumps under tmp/. With the current swp001 settings the
    // efficiency figures are plotted from the saved binary, so no figure is lost.
  , "tf_out_g2bg_ascii": $tf_out_g2bg_ascii
    // true saves det/arrdet_g3vox_prior<suffix>.bin (about 160 MB each, three files
    // when tf_prior_error is true). Needed to plot the prior array from the binary.
  , "tf_save_arrdet_prior": false
  },
  "PATH_LENGTH_PARAMETERS": {
    "name": "range2000_pitch1" // Name for PATH_LENGTH_PARAMETERS section. Set freely.
  , "tf_add_PLDL": true // Whether to add computed path length to detector (default:true)
  , "tf_incr_nhit_det": true // Whether to increment nhit_det of g3vox (default:true)
  , "tf_incr_nhit_ele": true // Whether to increment nhit_ele of g3vox (default:true)
  , "reference_matPL_sparse": 1 // Typical value order for sparse matrix
  , "epsilon_matPL_sparse": 1.0e-4 // Precision threshold to treat path-length sparse matrix entries as zero
  , "BL_max": $BL_max // meters: maximum density length
  , "PL_min": 0 // meters: minimum path length
  , "PL_max": $PL_max // meters: maximum path length
  , "PL_pit": 1.0 // meters: path length pitch (currently unused)
  , "tf_add_shell": true // Whether to add path length / density length of g2pil_shell
  , "tf_load_arrdet_g2pil": false // Whether to load DetectorPanelArray with g2pil path length (default:false)
  , "tf_save_arrdet_g2pil": false // Whether to save g2pil DetectorPanelArray (default:false)
  , "path_arrdet_g2pil_bin": "arrdet_g2pil.tmp.bin" // Binary filepath for arrdet_g2pil
  , "tf_load_arrdet_g3vox": false // Whether to load arrdet_g3vox and vec_spmat_PL
  , "tf_save_arrdet_g3vox": false // Whether to save arrdet_g3vox and vec_spmat_PL
  , "path_arrdet_g3vox_bin": "arrdet_g3vox.tmp.bin" // Binary filepath for arrdet_g3vox
  , "path_vec_spmat_PL_bin": "vec_spmat_PL.tmp.bin" // Binary filepath for vec_spmat_PL
  , "tf_load_bin_obs_mat_dNdD": false // Whether to load instead of running calc_dNdD::make_grouped_alldet_mat_dNdD
  , "tf_save_bin_obs_mat_dNdD": false // Whether to save large_merged_mat_dNdD as a binary file
  , "path_bin_obs_mat_dNdD": "obs_mat_dNdD.tmp.bin" // Filepath for binary matrix data
  },
  "BIN_GROUP_PARAMETERS": {
    "name": "prm_bingroup_g2pil" // Name for BIN_GROUP_PARAMETERS section. Set freely.
  , "signal_init": 0 // Initial value of signal
  , "noise_init": 0 // Initial value of noise
  , "is_avail_init": true // Initial is_avail
  , "PL_thres": 0.0 // PathLength threshold (meters) or DensityLength threshold (kg/m2)
  , "DL_thres": 10000 // DensityLength threshold value
  , "is_avail_under_thres": false // is_avail setting when below PL/DL threshold
  , "signal_under_thres": -1.0 // signal value when below threshold
  , "noise_under_thres": -1.0 // noise value when below threshold
  , "tf_run_1st_grouping": true // Whether to perform initial bin grouping
  , "tf_run_auto_grouping": true // Whether to perform automatic bin grouping
  , "igroup_start": 0 // Starting group ID
  , "nx_div_init": 4 // Initial number of divisions (x direction)
  , "ny_div_init": 2 // Initial number of divisions (y direction)
  , "signal_noise_group_trig": 200 // Groups with signal+noise above this value are subdivided
  , "ixlen_min": 5 // Bins with ixlen below this are not divided
  , "iylen_min": 5 // Bins with iylen below this are not divided
  , "tf_prefer_split_x": false // Whether to split in x direction when xlen==ylen in AutoGrouping
  , "nloop_limit": 10000 // Maximum loop count
  , "n_detector_grouping_manual": 0 // Number of detectors for manual grouping file
  , "vec_tf_read_bin_group_list": [] // Flags for manual grouping per detector
  , "vec_file_path_bin_group_list": [] // File paths for manual grouping lists
  },
  "GRID2D_PILLAR_PARAMETERS": {
    "name": "g2pil_${z_pitch}m_${base_density_int}kgm3" // Name for GRID2D_PILLAR_PARAMETERS. Set freely.
  , "initial_uniform_density": $base_density_int // density: kg/m3 (e.g. water 1.0 g/cm3 = 1000 kg/m3)
  , "path_dem": "$dem_path" // Path to DEM file
  , "zmin": $g2pil_zmin // zmin value (sea level for g2pil)
  , "tf_shift_x": true // Whether to shift in x direction
  , "tf_shift_y": true // Whether to shift in y direction
  , "tolerance_ratio": 1.0e-4 // tolerance ratio
    // Below are optional density structure parameters
  , "vertical_dike_params": [
    ],
    "vertical_cylinder_params": [
    ],
    "vertical_checkerboard_params": [
    ]
  },
  "GRID3D_VOXEL_PARAMETERS": {
    "tf_build_g3vox": true // Whether to build g3vox
  , "name": "g3vox_${z_pitch}m_${base_density_int}kgm3" // Name for GRID3D_VOXEL_PARAMETERS instance. Set freely.
  , "use_grid2d_of_g2pil": true // Use x/y bounds from Grid2dPillar
  , "zmin": $g3vox_zmin // zmin value
  , "zmax": $g3vox_zmax // zmax value
  , "z_pitch": $z_pitch // z pitch
  , "n_hit_det_min": $n_hit_det_min // Minimum nhit_det
  , "n_hit_det_max": 9999 // Maximum nhit_det
  , "n_hit_ele_min": 0 // Minimum nhit_ele
  , "n_hit_ele_max": 999999 // Maximum nhit_ele
  , "tf_end_after_merged": false // Whether to terminate after building g3vox_input
    // Optional density structure parameters
  , "checkerboard_3d_params": [
$checkerboard_3d_entries
    ],
    "cuboid_params": [
$cuboid_params_entries
    ],
    "merge_params": {
        "tf_exec": true // Execution flag
      , "name": "merge_01" // Merge parameter name
      , "x_merge_center": $center_x // Merge center x
      , "y_merge_center": $center_y // Merge center y
      , "z_merge_center": $z_merge_center // Merge center z
      , "x_merge_factor": $merge_factor // Merge factor in x direction
      , "y_merge_factor": $merge_factor // Merge factor in y direction
      , "z_merge_factor": $merge_factor // Merge factor in z direction
    }
  $reconst_voxels_block
  },

  "FLUX_RANGE_DATA_TABLE_PRIOR": {
    "pathin_log_peneflux": "../../../fluxtable/$flux_groom/costhz-kgm2-log10peneflux-allinone.g2zbin" // Path to log penetrating flux
  , "tf_xcnt_peneflux": true // If true, costhz returns xcnt instead of xlow
  , "tf_ycnt_peneflux": false // If false, y returns ylow
  , "pathin_dFdR_R_costhz": "../../../fluxtable/$flux_groom/g2_dFdR_R_costhz.g2zbin" // Path to dF/dR vs R vs costhz
  , "tf_xcnt_dFdR": true // If true, costhz returns xcnt instead of xlow
  , "tf_ycnt_dFdR": false // If false, y returns ylow
  , "tf_check_dFdR_divergence": true // Enable dF/dR divergence (noise) check after load. Do not change under normal use.
  , "dFdR_divergence_threshold": 0.05 // Sign-flip ratio threshold per costhz slice. Do not change under normal use.
  , "tf_dFdR_divergence_fatal": true // true = abort on a divergent dF/dR table (honda table must stay healthy). Do not change under normal use.
  },

  "FLUX_RANGE_DATA_TABLE_REAL": {
    "pathin_log_peneflux": "../../../fluxtable/$flux_groom/costhz-kgm2-log10peneflux-allinone.g2zbin" // Path to log penetrating flux
  , "tf_xcnt_peneflux": true // If true, costhz returns xcnt instead of xlow
  , "tf_ycnt_peneflux": false // If false, y returns ylow
  , "pathin_dFdR_R_costhz": "../../../fluxtable/$flux_groom/g2_dFdR_R_costhz.g2zbin" // Path to dF/dR vs R vs costhz
  , "tf_xcnt_dFdR": true // If true, costhz returns xcnt instead of xlow
  , "tf_ycnt_dFdR": false // If false, y returns ylow
  , "tf_check_dFdR_divergence": true // Enable dF/dR divergence (noise) check after load. Do not change under normal use.
  , "dFdR_divergence_threshold": 0.05 // Sign-flip ratio threshold per costhz slice. Do not change under normal use.
  , "tf_dFdR_divergence_fatal": true // true = abort on a divergent dF/dR table (honda table must stay healthy). Do not change under normal use.
  },

  "NOISE_PARAMETERS": {
    "name": "noise_tmp" // Name for NOISE_PARAMETERS section. Set freely.
  , "tf_exec": $noise_tf_exec // Execution flag
  , "flux_proport_ratio_floor" : $flux_proport_ratio_floor // Flux-proportional noise, deterministic floor (no fluctuation)
  , "flux_proport_ratio_poisson" : $flux_proport_ratio_poisson // Flux-proportional noise, Poisson-fluctuated component
  , "SOT_proport_noise_ratio_floor" : $SOT_proport_noise_ratio_floor // Angular-independent noise, deterministic floor (no fluctuation)
  , "SOT_proport_noise_ratio_poisson" : $SOT_proport_noise_ratio_poisson // Angular-independent noise, Poisson-fluctuated component
  , "user_defined_noise_flux_ratio" : 0.0 // User-defined noise as a fraction of flux
  , "path_user_defined_noise_distribution": [
    //   "../../noise/user_defined_noise.txt"
    ]
  },

  "NAGAINV_PARAMETERS": [
    {
        "name": "nagainv72_0"
      , "tf_exec": true
      , "tf_logN": false
      , "tf_signal_poisson": $nagainv_tf_signal_poisson
      // Add the efficiency-uncertainty variance to the C_N diagonal (default false = unchanged).
      , "tf_eff_cn_diag": $nagainv_tf_eff_cn_diag
      // Treat the efficiency uncertainty as independent within a merged group (default false = common).
      , "tf_eff_cn_diag_independent": $nagainv_tf_eff_cn_diag_independent
      , "tf_calc_chi2ndf": $nagainv_tf_calc_chi2ndf // chi2/ndf + p_eff=trace(R) self-eval; default false = off, no extra work
      , "nmuon_thres": 1.0E-10
      , "nmuon_under_thres": 1.0E-10
      , "uniform_prior_density": $base_density_int
$shell_density_lines
      , "corr_length": $corr_length
      , "sigma_rho": $sigma_rho
      , "sigma_rho_diag": $sigma_rho
      , "tf_aniso": false
      , "corr_length_xy": $corr_length
      , "corr_length_z": $corr_length
      , "aniso_cov_type": "separable"
    }
  ]
, "NAGAINV_PARAM_SWEEP": {
    "tf_exec": false                        // Set true to enable parameter sweep
  , "base_index": 0                         // Index in NAGAINV_PARAMETERS to use as template
  , "vec_sigma_rho": [200, 300]             // Values to sweep [kg/m³] (required)
  , "vec_corr_length": [50, 70]             // Values to sweep [m] (required)
  // Density sweep values [kg/m³]. Each element can be:
  //   Scalar: applies to voxels and all shells equally (expanded to [v, v, v, v])
  //   Quad:   [prior_density, shell_upper, shell_lower, shell_lateral]
  // Scalar and Quad can be mixed in the same array.
  , "vec_uniform_prior_density": [$base_density_int]
  // "vec_uniform_prior_density": [         // Mixed example
  //   2000,                                //   scalar: all shells = 2000
  //   [2000,  500, 2500, 1800]             //   quad: per-shell override
  // ]
  , "link_diag_to_sigma": true              // true = sigma_rho_diag follows sigma_rho
  // , "vec_sigma_rho_diag": [200, 300]     // Required when link_diag_to_sigma=false [kg/m³]
  , "module8_mode": "$module8_mode"                   // "none" | "all" | "last" | "selected"
  , "module8_indices": []                   // Sweep indices for module8_mode="selected" (e.g. [0, 2, 5])
  }
, "PROJ_DENS_EVAL_GROUPED": {
    "tf_exec": $proj_dens_tf_exec // Execution flag
  , "tf_signal_poisson": $proj_dens_tf_signal_poisson // Whether to apply Poisson error to signal counts
  , "dens_min": 0 // Minimum search range
  , "dens_max": 10000 // Maximum search range
  , "dens_steps": [200, 50, 10] // Search step sizes; single value means uniform search
  , "range_factor": 2.0 // Search range: best_dens +/- (previous step) * range_factor
  , "sigma": 1.5 // sigma of sqrt(nmuon)
  }
, "Z_CROSS_SECTION": {
    "min": $zcross_min // z_cross_min value
  , "max": $zcross_max // z_cross_max value
  , "output_binary": true
  , "zstep": $zcross_zstep // Z-step [meters]. 0 = use z_interval
  }
}
