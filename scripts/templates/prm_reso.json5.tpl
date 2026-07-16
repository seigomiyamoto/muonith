// ${work_dir_name}/${depth_work_subdir}/prm_reso.json5
// depth_reso.exe configuration for $site_name_cap (heatmap analysis)
{
  "seed": 42,
  "PARAM_CONSTANTS": {          // ! DO NOT change these values unless necessary
    "DL_min": 5000.0,           // [kg/m^2] Minimum density length
    "DL_max": 4999000.0,        // [kg/m^2] Maximum density length
    "costhz_min": 0.0,          // Minimum cos(zenith angle)
    "costhz_max": 1.0,          // Maximum cos(zenith angle)
    "BL_max_default": 5000.0    // [m] Maximum beam length (default)
  },
  "LOG_FILE": {
    "path_log_dir": "logs",
    "give_new_number": false,
    "log_level": {
      "stdout_level": "debug",
      "stderr_level": "error",
      "file_level": "trace"
    },
    "default_max_count": 20,
    "archive_existing": false
  },
  "DETECTOR_PARAMETER_LISTS": {
    "name": "arrdet_${det_list_name}",
    "det_files": [
$det_files_block
    ],
    "tf_apply_eff": false
  },
  "PATH_LENGTH_PARAMETERS": {
    "name": "range2000_pitch1",
    "tf_add_PLDL": true,
    "tf_incr_nhit_det": true,
    "tf_incr_nhit_ele": true,
    "reference_matPL_sparse": 1,
    "epsilon_matPL_sparse": 0.0001,
    "PL_min": 0,
    "PL_max": 2000,
    "PL_pit": 1.0,
$BL_max_line
    "tf_add_shell": true,
    "tf_load_arrdet_g2pil": false,
    "tf_save_arrdet_g2pil": false,
    "path_arrdet_g2pil_bin": "arrdet_g2pil.tmp.bin",
    "tf_load_arrdet_g3vox": false,
    "tf_save_arrdet_g3vox": false,
    "path_arrdet_g3vox_bin": "arrdet_g3vox.tmp.bin",
    "path_vec_spmat_PL_bin": "vec_spmat_PL.tmp.bin",
    "tf_load_bin_obs_mat_dNdD": false,
    "tf_save_bin_obs_mat_dNdD": false,
    "path_bin_obs_mat_dNdD": "obs_mat_dNdD.tmp.bin"
  },
  "BIN_GROUP_PARAMETERS": {
    "name": "prm_bingroup_g2pil",
    "signal_init": 0,
    "noise_init": 0,
    "is_avail_init": true,
    "PL_thres": 0.0,
    "DL_thres": 10000,
    "is_avail_under_thres": false,
    "signal_under_thres": -1.0,
    "noise_under_thres": -1.0,
    "tf_run_1st_grouping": true,
    "tf_run_auto_grouping": true,
    "igroup_start": 0,
    "nx_div_init": 4,
    "ny_div_init": 2,
    "signal_noise_group_trig": 100,
    "ixlen_min": 5,
    "iylen_min": 5,
    "tf_prefer_split_x": false,
    "nloop_limit": 10000,
    "n_detector_grouping_manual": 0,
    "vec_tf_read_bin_group_list": [],
    "vec_file_path_bin_group_list": []
  },
  "BINGROUP_PARAMETERS": {},
  "NOISE_PARAMETERS": {},
  "GRID2D_PILLAR_PARAMETERS": {
    "name": "$grid2d_pillar_name",
    "initial_uniform_density": $grid2d_initial_uniform_density,
    "path_dem": "$dem_path",
    "zmin": 0.0,
    "tf_shift_x": true,
    "tf_shift_y": true,
    "tolerance_ratio": 0.01,
    "vertical_dike_params": [],
    "vertical_cylinder_params": [],
    "vertical_checkerboard_params": []
  },
  "FLUX_RANGE_DATA_TABLE_REAL": {
    "pathin_log_peneflux": "../../../fluxtable/$flux_groom/costhz-kgm2-log10peneflux-allinone.g2zbin",
    "tf_xcnt_peneflux": true,
    "tf_ycnt_peneflux": false,
    "pathin_dFdR_R_costhz": "../../../fluxtable/$flux_groom/g2_dFdR_R_costhz.g2zbin",
    "tf_xcnt_dFdR": true,
    "tf_ycnt_dFdR": false,
    "tf_dFdR_divergence_fatal": false  // depth-reso path is not in active use: warn only, do not abort on divergence
  },
  "DEPTH_RESOLUTION_SWEEP": {
      "mode": "ddl"  // "cylinder" or "ddl"/"baumkuchen", default is "cylinder"
    , "common": {
        "output_ascii_prefix": "depth_res"
      , "x_cnt_obj": $x_cnt_obj      // $site_name_cap target center x [m] (EPSG:$epsg)
      , "y_cnt_obj": $y_cnt_obj     // $site_name_cap target center y [m] (EPSG:$epsg)
      , "base_density": $base_density   // Base density [kg/m3]
      , "obj_size_upper_limit": $obj_size_upper_limit  // Max anomaly size [m]
      , "obj_size_lower_limit": $obj_size_lower_limit    // Min anomaly size [m]
      , "elev_center_step": $elev_center_step     // Elevation step [rad]
      , "angle_between_cut_factor": 1.0
      , "sweep_range_factor": 1.0 // ! DO NOT TOUCH THIS FACTOR !! KEEP IT 1.0 FOR NOW !!
      , "vec_delta_density": $vec_delta_density
      , "stat_alphas": [0.3, 0.1]
      , "both_side": false           // true = two-sided test, false(default) = one-sided test
$tf_out_det_PL_signal_line
      // Angle bin override: set to true to use finer angle bins without separate detector files
      , "tf_override_angle_bin": $tf_override_angle_bin
      , "angle_bin_override": {
          "nbinx": $abo_nbinx, "txmin": $abo_txmin, "txmax": $abo_txmax
        , "nbiny": $abo_nbiny, "tymin": $abo_tymin,  "tymax": $abo_tymax
      }
      , "signal_noise_amplifiers": [
$signal_noise_amplifiers_block
      ]
    }
    , "cylinder": {}
    , "ddl": {
        "obj_size_step": $obj_size_step  // Size step [m] (DDL mode only)
      , "depth_step": $depth_step     // Depth step [m]
    }
  }
}
