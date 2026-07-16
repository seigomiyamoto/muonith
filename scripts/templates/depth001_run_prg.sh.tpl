#!/usr/bin/env bash
# ${work_dir_name}/${depth_work_subdir}/run_prg.sh
# Depth vs resolution analysis for $site_name_cap (heatmap version)

function usage() {
  echo "Usage: $$0 <runexe: true|false> [seed]"
  exit 1
}

if [ $$# -lt 1 ]; then
  usage
fi

# Argument parsing
RUN_EXE="$$1"
SEED_INPUT="$$2"

epsg_code=$epsg  # $site_name_cap

echo "RUN_EXE: $$RUN_EXE"
echo "SEED_INPUT: $$SEED_INPUT"

data_dir="tmp"
fig_dir="figs"
log_dir="logs"

mkdir -p "$$data_dir" "$$fig_dir" "$$log_dir"

# Execution parameters
EXE="../../../build-release/exec/depth_reso.exe"
PRMRC="prm_reso.json5"

# Plot detector positions and terrain
uv run python ../../../scripts/plot_dem_with_detectors.py \
 --config $$PRMRC --output "fig_${site_name}_detectors.png" \
 --contour --cint 10 --dpi 300 --show_distance \
 --shade --shade-alpha 0.2 \
 --zmin $plot_zmin --zmax $plot_zmax --xyunit "km" \
 --xcnt $center_x --ycnt $center_y  --epsg $$epsg_code

if [ "$$(uname -s)" = "Darwin" ]; then
  NUM_THREADS=$$(( $$(sysctl -n hw.logicalcpu) - 3 ))
else
  NUM_THREADS=$$(( $$(nproc) - 3 ))
fi
export OMP_NUM_THREADS=$$NUM_THREADS
export OPENBLAS_NUM_THREADS=$$NUM_THREADS
echo "NUM_THREADS: $$NUM_THREADS (logical CPUs - 3)"

set -e
###################################################################
if [ "$$RUN_EXE" = "true" ]; then
  rm -rf "$$data_dir"
  rm -f "$$log_dir"/*.*
  mkdir -p "$$data_dir"
  if [ $$# -ge 2 ]; then
    SEED_INPUT="$$2"
    echo "Seed specified: $$SEED_INPUT"
    /usr/bin/time -v -o "$$log_dir/time_report.log" $$EXE --json $$PRMRC --seed $$SEED_INPUT
  else
    echo "No seed specified. Using seed from JSON."
    /usr/bin/time -v -o "$$log_dir/time_report.log" $$EXE --json $$PRMRC
  fi
else
  echo "Executable run skipped (runexe=$$RUN_EXE)"
fi
#####################################################################

set +e  # Ignore errors for plotting
rm -f "$$fig_dir"/*.*

# Plot heatmaps (detectors and signal_amps are read from prm_reso.json5)
uv run python ../../../scripts/plot_detection_limit_heatmap.py \
  --main_config prm_reso.json5 \
  --plot_config heatmap_config.json5 \
  --surface_elevation $surface_elevation \
  --datadir $$data_dir \
  --outdir $$fig_dir

echo "Done. Heatmap plots saved in $$fig_dir/"

# auto_plot: 2D histogram of base PL/signal distribution
uv run python ../../../scripts/auto_plot.py --config auto_plot.json5 --verbose \
 --error-log "$$log_dir/error_auto_plot.log"
