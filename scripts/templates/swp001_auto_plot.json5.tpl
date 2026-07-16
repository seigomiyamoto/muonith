{
  // ==========================================================================
  // $site_name_cap swp001 configuration for auto_plot.py
  // Based on showa/swp001/auto_plot.json5
  // ==========================================================================

  "hist2d_exec": true,
  "g2bg_exec": true,
  "g3vox_exec": true,
  // det: binary-direct txty/g2bg figures via plot_det_arrdet.py, reading
  // det/arrdet_*.bin checkpoints. Enabled because the swp001 prm disables the
  // tmp/ ASCII dumps (tf_out_txty_ascii / tf_out_g2bg_ascii false), so the
  // hist2d/g2bg families above find no .tmp files. Set false if you re-enable the
  // ASCII dumps, to avoid producing the same figures twice.
  "det_exec": true,

  // ===========================================================================
  // hist2d: 2D histogram plots
  // ===========================================================================
  "hist2d": {
    "save_csv": false,
    "template": {
      "log": false,
      "binx": "$nbinx $txmin $txmax",
      "biny": "$nbiny $tymin $tymax",
      "binz": "20 0 $PL_hist_max",
      "dpi": 300,
      "colormap": "jet",
      "color_over": "black",
      "color_under": "white",
    },
    "patterns": {
      "arrdet_*_signal_*": {
        "log": true,
        "colormap": "gnuplot2",
        "color_over": "#e0ffff",
        "color_under": "#2f4f4f",
        "binz": "$signal_binz",
      },
      "arrdet_*_dens_*": {
        "binz": "20 1500.0 2500",
      },
      "arrdet_*_PL_*": {
        "colormap": "gnuplot",
        "color_over": "#ffff99",
        "color_under": "#000033",
        "binz": "20 0.0 $PL_hist_max",
      },
    },
    "exclude_patterns": [
      "arrdet_g3vox_prior_center_txty_*",
      "arrdet_g3vox_prior_upper_txty_*",
      "arrdet_g3vox_prior_lower_txty_*",
    ],
  },

  // ===========================================================================
  // g2bg: Grid2D background plots
  // ===========================================================================
  "g2bg": {
    "auto_detect": true,
    "exclude_patterns": [],
    "template": {
      "input_dir": "tmp",
      "input_suffix": ".tmp",
      "output_dir": "figs",
      "output_suffix": "",
      "png_dpi": 300,
      "make_pdf_per_detid": false,
      "make_gif_per_detid": false,
      "make_pdf_per_field": true,
      "make_gif_per_field": true,
      "pdf_suffix": "",
      "pdf_dpi": 300,
      "gif_suffix": "",
      "gif_delay": 400,
      "gif_loop": 0,
      "gif_dpi": 300,
      "exec_erace_pngs": false,
      "rect_edge_color": "gray",
      "bg_color": "silver",
      "grid_color": "gray",
      "gridline_type": "dotted",
      "gridline_width": 1.0,
      "plot_fields": {
        "sig":      { "exec": true,  "logz": true,  "pallet": "gnuplot2", "color_under": "#2f4f4f", "color_over": "#e0ffff", "ngrad": 20, "vmin": 1.0E+1, "vmax": 1.0E+3, "alpha": 1.0 },
        "noi":      { "exec": true,  "logz": true,  "pallet": "cividis", "color_under": "white", "color_over": "orange", "ngrad": 20, "vmin": 1.0E-1, "vmax": 1.0E+3, "alpha": 1.0 },
        "sig_over_noi": { "exec": true,  "logz": true,  "pallet": "RdYlGn", "color_under": "#4d0000", "color_over": "cyan", "ngrad": 20, "vmin": 1.0E-1, "vmax": 1.0E+3, "alpha": 1.0 },
        "signoi":   { "exec": true,  "logz": true,  "pallet": "gnuplot2", "color_under": "#2f4f4f", "color_over": "#e0ffff", "ngrad": 20, "vmin": 1.0E+1, "vmax": 1.0E+3, "alpha": 1.0 },
        "dens_cnt": { "exec": true,  "logz": false, "pallet": "jet", "color_under": "silver", "color_over": "#4d0000", "ngrad": 20, "vmin": $g2bg_dens_cnt_vmin, "vmax": $g2bg_dens_cnt_vmax, "alpha": 1.0 },
        "dens_err": { "exec": true,  "logz": false, "pallet": "inferno", "color_under": "black", "color_over": "white", "ngrad": 15, "vmin": 0.0, "vmax": 300.0, "alpha": 1.0 },
        "eff_cnt":  { "exec": true,  "logz": false, "pallet": "viridis", "color_under": "silver", "color_over": "#4d0000", "ngrad": 20, "vmin": 0.87, "vmax": 0.99, "alpha": 1.0 },
        "eff_err":  { "exec": true,  "logz": false, "pallet": "inferno", "color_under": "black", "color_over": "white", "ngrad": 20, "vmin": 0.0, "vmax": 0.05, "alpha": 1.0 },
      },
    },
  },

  // ===========================================================================
  // g3vox: Grid3D voxel plots
  // ===========================================================================
  "g3vox": {
    "auto_binning": "data",
    "auto_detect": true,
    "det_info_path": "prm_muonith.json5",
    "exclude_patterns": [
      "g3vox_*_lower_*",
      "g3vox_*_upper_*",
    ],

    "render_patterns": {
      "*_error*": {
        "colormap": "YlOrRd",
        "vmin":    0,
        "vmax": $error_vmax,
        "cbar_tick_step": $error_cbar_step,
        "ngrad":  20,
        "overflow_color": "black",
        "underflow_color": "white",
      },
      "*_diff*": {
        "colormap": "seismic",
        "vmin": -500,
        "vmax": 500,
        "ngrad": 50,
      },
      "*_delta*": {
        "colormap": "seismic",
        "vmin": -500,
        "vmax": 500,
        "ngrad": 50,
      },
      "*_dens*": {
        "colormap": "jet",
        "vmin": 1500,
        "vmax": 2500,
        "ngrad": 60,
      },
      "*_prior*": {
        "colormap": "jet",
        "vmin": 0,
        "vmax": 3000,
        "ngrad": 60,
      },
      "*_rec*": {
        "colormap": "jet",
        "vmin": 0,
        "vmax": 3000,
        "ngrad": 60,
      },
      "*_input*": {
        "colormap": "jet",
        "vmin": 0,
        "vmax": 3000,
        "ngrad": 60,
      },
      "*_shell*": {
        "colormap": "jet",
        "vmin": 0,
        "vmax": 3000,
        "ngrad": 60,
      },
    },

    "template": {
      "run_params": {
        "filename_in": "",
        "tf_out_each_det": false,
        "n_hit_det_thres": 0,
        "cross_section": ["-z"],

        "zmin": $g3vox_plot_zmin,
        "zmax": $g3vox_plot_zmax,
        "zstep": $g3vox_plot_zstep

        , "xmin": $g3vox_plot_xmin
        , "xmax": $g3vox_plot_xmax
        , "xstep": $g3vox_plot_xstep

        , "ymin": $g3vox_plot_ymin
        , "ymax": $g3vox_plot_ymax
        , "ystep": $g3vox_plot_ystep
      },

      "render_params": {
        "output_png_prefix": "",
        "dpi": 300,
        "save_png_dir": "figs",

        "xy_unit": "km",

        "output_formats": ["gif"],

        "colormap": "jet",
        "bg_color": "silver",
        "overflow_color": "black",
        "underflow_color": "white",

        "grid_color": "grey",
        "gridline_type": "dotted",
        "gridline_width": 0.5,

        "ngrad": 60,
        "vmin": 0.0,
        "vmax": 3000,

        "mask_color": "white",
        "mask_alpha": 0.3,

        "main_width": 10,
        "cbar_width": 0.5,

        "title_fontsize": 18,
        "xy_label_fontsize": 16,
        "xy_tick_fontsize": 14,
        "cbar_label_fontsize": 16,
        "cbar_tick_fontsize": 14,

        "det_pos_marker_type": "*",
        "det_pos_makker_size": 16,
        "det_pos_makker_color": "black",

        "det_name_show": true,
        "det_name_fontsize": 10,
        "det_name_offset_y": 0.02,
        "det_name_color": "black",
      },

      "gifski_params": {
        "fps": 0.5,
        "quality": 100,
      },
    },
  },

  // ===========================================================================
  // det: binary-direct figures via plot_det_arrdet.py (gated by det_exec above)
  // ===========================================================================
  "det": {
    "auto_detect": true,
    // Where the figures go; a relative path is anchored at this file's folder.
    "out_dir": "figs",
    // Drawing workers: positive = that many, 1 = one by one, 0 = all CPU
    // cores, negative = leave that many cores free (-2 on 32 cores -> 30).
    "n_jobs": -2,
    // The g2pil_naive checkpoint figures are not needed in swp001.
    "exclude_patterns": ["arrdet_g2pil_*"],
  },
}
