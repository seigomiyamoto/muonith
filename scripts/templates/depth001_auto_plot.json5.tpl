{
  // $site_name_cap depth001: auto_plot config for base PL/signal distribution
  "hist2d_exec": true,
  "g2bg_exec": false,
  "g3vox_exec": false,

  "hist2d": {
    "save_csv": false,
    "template": {
      "log": false,
      "binx": "320 -1.6 1.6",
      "biny": "160 0.00 1.6",
      "binz": "20 0 $PL_hist_max",
      "dpi": 300,
      "colormap": "jet",
      "color_over": "black",
      "color_under": "white",
    },
    "patterns": {
      "fig_*_txty_signal_*": {
        "log": true,
        "colormap": "gnuplot2",
        "color_over": "#e0ffff",
        "color_under": "#2f4f4f",
        "binz": "18 1.0E-3 1.0E+3",
      },
      "fig_*_txty_PL_*": {
        "colormap": "gnuplot",
        "color_over": "#ffff99",
        "color_under": "#000033",
        "binz": "20 0.0 $PL_hist_max",
      },
    },
    "exclude_patterns": [],
  },
}
