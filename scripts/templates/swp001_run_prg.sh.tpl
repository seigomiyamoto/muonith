#!/usr/bin/env bash

function usage() {
  echo "Usage:"
  echo "  $$0 <runexe: true|false> [seed] [end_stage]"
  echo "  $$0 <prm.json5> <runexe: true|false> [seed] [end_stage]"
  echo "  end_stage: stop after module N (valid: 3, 4, 5, 6, 7, 8). Default: run all."
  exit 1
}

if [ $$# -lt 1 ]; then
  usage
fi

if [ "$$1" = "true" ] || [ "$$1" = "false" ]; then
  PRMRC="prm_muonith.json5"
  RUN_EXE="$$1"
  SEED_INPUT="$${2:-}"
  END_STAGE="$${3:-}"
else
  if [ $$# -lt 2 ]; then
    usage
  fi
  PRMRC="$$1"
  RUN_EXE="$$2"
  SEED_INPUT="$${3:-}"
  END_STAGE="$${4:-}"
fi

if [ "$$RUN_EXE" != "true" ] && [ "$$RUN_EXE" != "false" ]; then
  usage
fi

echo "PRMRC: $$PRMRC"
echo "RUN_EXE: $$RUN_EXE"
echo "SEED_INPUT: $$SEED_INPUT"
if [ -n "$$END_STAGE" ]; then
  echo "END_STAGE: $$END_STAGE"
fi

data_dir="tmp"
fig_dir="figs"
log_dir="logs"

mkdir -p "$$data_dir" "$$fig_dir" "$$log_dir"

EXE="../../../build-release/exec/muonith.exe"

${dem_gen_block}EPSG=$epsg
uv run python ../../../scripts/plot_dem_with_detectors.py \
 --config $$PRMRC --output "fig_${site_name}_${swp_work_subdir}_detectors.png" \
 --contour --cint 10 --dpi 300 --show_distance \
 --shade --shade-alpha 0.2 \
 --zmin $plot_zmin --zmax $plot_zmax --xyunit "km" \
 --xcnt $center_x --ycnt $center_y --epsg $$EPSG

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
  END_STAGE_FLAG=""
  if [ -n "$$END_STAGE" ]; then
    END_STAGE_FLAG="--end-stage $$END_STAGE"
  fi
  if [ "$$(uname -s)" = "Darwin" ]; then
    TIME_CMD="gtime" # GNU time on macOS (brew install gnu-time); same -v output as Linux /usr/bin/time
  else
    TIME_CMD="/usr/bin/time"
  fi
  if [ -n "$$SEED_INPUT" ]; then
    echo "Seed specified: $$SEED_INPUT"
    $$TIME_CMD -v -o "$$log_dir/time_report.log" $$EXE --json $$PRMRC --seed $$SEED_INPUT $$END_STAGE_FLAG
  else
    echo "No seed specified. Using seed from JSON."
    $$TIME_CMD -v -o "$$log_dir/time_report.log" $$EXE --json $$PRMRC $$END_STAGE_FLAG
  fi
else
  echo "Executable run skipped (runexe=$$RUN_EXE)"
fi
#####################################################################
set +e

rm -f "$$fig_dir"/*.*

uv run python ../../../scripts/auto_plot.py --config auto_plot.json5 --verbose \
 --error-log "$$log_dir/error_auto_plot.log"
