#!/usr/bin/env bash
# ============================================================
# Set up a volcano analysis work site
#
# Usage:
#   bash setup_station.sh <station_name> [--make-detparams] [--force] [--clean]
#   bash setup_station.sh <station_name> [--skip-env-check]
#
# --make-detparams builds the detector runcards (det_XX.json5) from the
# KML/GeoJSON detector positions. It does not start the analysis; that is
# done later by run_prg.sh inside work/<station_name>/.
#
# Examples:
#   bash setup_station.sh tutorial                  # license-free synthetic tutorial (start here)
#   bash setup_station.sh tarumae_base              # skeleton only
#   bash setup_station.sh tarumae_base --make-detparams  # skeleton + detector runcards
#   bash setup_station.sh meakan --make-detparams --force  # regenerate everything
#   bash setup_station.sh omuro --make-detparams    # first-time setup with checks
#   bash setup_station.sh omuro --make-detparams --skip-env-check  # skip apt/venv checks
#
# Prerequisites (checked automatically unless --skip-env-check):
#   /usr/bin/time    -- apt install -y time
#   python3          -- apt install -y python3 python3-venv
#   .venv/           -- created automatically if absent
#   param_sites/<station_name>.json5  -- copied from template if absent
#   data/<station_name>/              -- created if absent; DEM (+KML) required for --make-detparams
#
# For a brand new site:
#   1. bash setup_station.sh <name>          # auto-copies template; follow instructions
#   2. Edit param_sites/<name>.json5
#   3. Copy DEM (.g2zbin) and KML/GeoJSON into data/<name>/
#      (KML/GeoJSON not required if "skip_detparams": true in the json5)
#   4. bash setup_station.sh <name> --make-detparams
# ============================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
if [ $# -lt 1 ]; then
  echo "Usage: bash setup_station.sh <station_name> [--make-detparams] [--force] [--clean] [--skip-env-check]"
  echo ""
  echo "Available stations:"
  for f in param_sites/*.json5; do
    name="$(basename "$f" .json5)"
    [ "$name" = "template" ] && continue
    echo "  $name"
  done
  exit 1
fi

STATION_NAME="$1"
shift

SKIP_ENV_CHECK=0
EXTRA_ARGS=()
for arg in "$@"; do
  case "$arg" in
    --skip-env-check) SKIP_ENV_CHECK=1 ;;
    *) EXTRA_ARGS+=("$arg") ;;
  esac
done

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
STATION_JSON5="$SCRIPT_DIR/param_sites/$STATION_NAME.json5"
VENV_DIR="$SCRIPT_DIR/.venv"
VENV_PYTHON="$VENV_DIR/bin/python3"
VENV_PIP="$VENV_DIR/bin/pip"

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
info()  { echo "[INFO]  $*"; }
warn()  { echo "[WARN]  $*" >&2; }
error() { echo "[ERROR] $*" >&2; exit 1; }

# ---------------------------------------------------------------------------
# Phase 0: Environment check
# ---------------------------------------------------------------------------
if [ "$SKIP_ENV_CHECK" -eq 0 ]; then
  echo "=== Phase 0: Environment check ==="

  # 0-1. /usr/bin/time
  if [ ! -x /usr/bin/time ]; then
    error "/usr/bin/time not found. Install with: sudo apt install -y time"
  fi
  info "/usr/bin/time ... OK"

  # 0-2. python3
  if ! command -v python3 > /dev/null 2>&1; then
    error "python3 not found. Install with: sudo apt install -y python3 python3-venv"
  fi
  info "python3 ... OK ($(python3 --version 2>&1))"

  # 0-3. taskset (warn only; present on Ubuntu 24.04 by default)
  if ! command -v taskset > /dev/null 2>&1; then
    warn "taskset not found (util-linux). run_prg.sh may fail."
  else
    info "taskset ... OK"
  fi

  # 0-4. nproc (warn only)
  if ! command -v nproc > /dev/null 2>&1; then
    warn "nproc not found (coreutils). run_prg.sh may fail."
  else
    info "nproc ... OK ($(nproc) logical CPUs)"
  fi

  # 0-5. Python venv
  # The package list lives in pyproject.toml (repo root); uv.lock pins versions.
  # With uv available, "uv sync" creates/updates .venv from that single source.
  # Without uv, fall back to venv+pip (keep the fallback list in sync with
  # pyproject.toml [project] dependencies).
  echo ""
  echo "=== Phase 0-5: Python venv ==="
  if command -v uv > /dev/null 2>&1; then
    info "uv found. Syncing .venv from pyproject.toml + uv.lock..."
    # --inexact: do not uninstall packages outside the lock (e.g. the docs
    # group installed via "uv sync --group docs").
    (cd "$SCRIPT_DIR" && uv sync --quiet --inexact)
    info ".venv/ synced via uv."
  elif [ ! -f "$VENV_PYTHON" ]; then
    warn "uv not found. Falling back to venv+pip (versions not pinned)."
    info ".venv/ not found. Creating..."
    python3 -m venv "$VENV_DIR"
    info "Installing required Python packages..."
    "$VENV_PIP" install --quiet --upgrade pip
    "$VENV_PIP" install --quiet numpy matplotlib pandas scipy requests json5 pyproj rasterio Pillow
    info ".venv/ created and packages installed."
  else
    warn "uv not found. Falling back to import check (versions not pinned)."
    info ".venv/ exists. Checking required packages..."
    MISSING=()
    for pkg in numpy matplotlib pandas scipy requests json5 pyproj rasterio PIL; do
      if ! "$VENV_PYTHON" -c "import $pkg" 2>/dev/null; then
        MISSING+=("$pkg")
      fi
    done
    if [ ${#MISSING[@]} -gt 0 ]; then
      warn "Missing packages: ${MISSING[*]}. Installing..."
      "$VENV_PIP" install --quiet "${MISSING[@]/#PIL/Pillow}"
    else
      info "Required Python packages ... OK"
    fi
  fi
fi

# ---------------------------------------------------------------------------
# Phase 1: param_sites/<name>.json5
# ---------------------------------------------------------------------------
echo ""
echo "=== Phase 1: Site config ==="

if [ ! -f "$STATION_JSON5" ]; then
  TEMPLATE="$SCRIPT_DIR/param_sites/template.json5"
  if [ ! -f "$TEMPLATE" ]; then
    error "Template not found: $TEMPLATE"
  fi
  cp "$TEMPLATE" "$STATION_JSON5"
  echo ""
  echo "  Created: $STATION_JSON5"
  echo ""
  echo "  *** ACTION REQUIRED ***"
  echo "  Edit the following fields in $STATION_JSON5:"
  echo "    station_name       -- set to \"$STATION_NAME\""
  echo "    epsg               -- projected CRS code (e.g. 6681 for JGD2011 zone 13)"
  echo "                          OR set to 0 with \"skip_detparams\": true for custom CRS"
  echo "    dem_file           -- relative path to .g2zbin under data/$STATION_NAME/"
  echo "    center_x, center_y -- volcano center in projected CRS [meters]"
  echo "    surface_elevation  -- summit elevation [m asl]"
  echo "    detector_input     -- relative path to .kml or .geojson (omit if skip_detparams)"
  echo ""
  echo "  Then re-run:"
  echo "    bash setup_station.sh $STATION_NAME --make-detparams"
  echo ""
  exit 0
fi
info "param_sites/$STATION_NAME.json5 ... OK"

# Detect skip_detparams flag (simple grep — does not parse JSON5 syntax fully)
SKIP_DETPARAMS=0
if grep -Eq '^[[:space:]]*"skip_detparams"[[:space:]]*:[[:space:]]*true' "$STATION_JSON5"; then
  SKIP_DETPARAMS=1
  info "skip_detparams=true detected; KML/GeoJSON requirement relaxed."
fi

# Detect data_source key (multiple stations may share one data dir; e.g.
# meakan01/02/03 all use data/meakan/, tarumae-* all use data/tarumae/).
# Falls back to STATION_NAME if data_source is absent.
DATA_SOURCE=$(grep -E '^[[:space:]]*"data_source"[[:space:]]*:' "$STATION_JSON5" 2>/dev/null \
              | sed -E 's/.*"data_source"[[:space:]]*:[[:space:]]*"([^"]+)".*/\1/' \
              | head -1 || true)
DATA_DIR_NAME="${DATA_SOURCE:-$STATION_NAME}"

# ---------------------------------------------------------------------------
# Phase 2: data/<name>/ and required input files
# ---------------------------------------------------------------------------
echo ""
echo "=== Phase 2: Data directory ==="

DATA_DIR="$SCRIPT_DIR/data/$DATA_DIR_NAME"
WORK_DIR="$SCRIPT_DIR/work/$STATION_NAME"
if [ ! -d "$DATA_DIR" ]; then
  if [ "$SKIP_DETPARAMS" -eq 1 ]; then
    info "data/$DATA_DIR_NAME/ not found; skip_detparams=true allows existing work/ inputs."
  else
    mkdir -p "$DATA_DIR"
    touch "$DATA_DIR/.gitkeep"
    info "Created: data/$DATA_DIR_NAME/"
  fi
fi
if [ "$DATA_DIR_NAME" != "$STATION_NAME" ]; then
  info "data_source=$DATA_DIR_NAME (shared with other stations)"
fi

# Check DEM and detector input exist when --make-detparams is requested
RUN_REQUESTED=0
for arg in "${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"}"; do
  [ "$arg" = "--make-detparams" ] && RUN_REQUESTED=1
done

if [ "$RUN_REQUESTED" -eq 1 ]; then
  # Count actual data files (resolve symlinks via -L; ignore .gitkeep)
  N_G2ZBIN=0
  N_KML=0
  if [ -d "$DATA_DIR" ]; then
    N_G2ZBIN=$(find -L "$DATA_DIR" -maxdepth 2 -name "*.g2zbin" 2>/dev/null | wc -l)
    N_KML=$(find -L "$DATA_DIR" -maxdepth 2 \( -name "*.kml" -o -name "*.geojson" \) 2>/dev/null | wc -l)
  fi
  N_WORK_G2ZBIN=0
  if [ "$SKIP_DETPARAMS" -eq 1 ] && [ -d "$WORK_DIR/dem" ]; then
    N_WORK_G2ZBIN=$(find -L "$WORK_DIR/dem" -maxdepth 1 -name "*.g2zbin" 2>/dev/null | wc -l)
  fi

  # Stations with a "dem_generator" key build their DEM inside run_prg.sh, so a
  # committed .g2zbin is never required (synthetic demo stations, e.g. eg-*).
  HAS_DEM_GENERATOR=0
  if grep -Eq '^[[:space:]]*"dem_generator"[[:space:]]*:' "$STATION_JSON5"; then
    HAS_DEM_GENERATOR=1
    info "dem_generator detected; DEM is generated by run_prg.sh."
  fi

  MISSING_INPUT=0
  if [ "$N_G2ZBIN" -eq 0 ] && [ "$N_WORK_G2ZBIN" -eq 0 ] && [ "$HAS_DEM_GENERATOR" -eq 0 ]; then
    MISSING_INPUT=1
  fi
  if [ "$SKIP_DETPARAMS" -eq 0 ] && [ "$N_KML" -eq 0 ]; then
    MISSING_INPUT=1
  fi

  if [ "$MISSING_INPUT" -eq 1 ]; then
    echo ""
    echo "  *** ACTION REQUIRED ***"
    echo "  data/$DATA_DIR_NAME/ is missing required input files."
    echo "  Please copy:"
    if [ "$N_G2ZBIN" -eq 0 ] && [ "$N_WORK_G2ZBIN" -eq 0 ]; then
      echo "    - DEM file (.g2zbin)  -->  data/$DATA_DIR_NAME/"
    fi
    if [ "$SKIP_DETPARAMS" -eq 0 ] && [ "$N_KML" -eq 0 ]; then
      echo "    - Detector file (.kml or .geojson)  -->  data/$DATA_DIR_NAME/"
    fi
    echo ""
    echo "  After copying, re-run:"
    echo "    bash setup_station.sh $STATION_NAME --make-detparams"
    echo ""
    exit 1
  fi
fi

info "data/$DATA_DIR_NAME/ ... OK"
if [ "$RUN_REQUESTED" -eq 1 ] && [ "$N_G2ZBIN" -eq 0 ] && [ "$N_WORK_G2ZBIN" -gt 0 ]; then
  info "Using existing work/$STATION_NAME/dem/ DEM because skip_detparams=true."
fi

# ---------------------------------------------------------------------------
# Phase 3: Run init_work_site.py
# ---------------------------------------------------------------------------
echo ""
echo "=== Phase 3: Generating work site: $STATION_NAME ==="

# Use venv python if available; fall back to system python3
if [ -f "$VENV_PYTHON" ]; then
  PYTHON="$VENV_PYTHON"
else
  PYTHON="python3"
fi

"$PYTHON" "$SCRIPT_DIR/scripts/init_work_site.py" --verbose "${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"}" "$STATION_JSON5"

# ---------------------------------------------------------------------------
# Done
# ---------------------------------------------------------------------------
echo ""
echo "=== Done: work/$STATION_NAME/ ==="
echo ""
echo "Next steps:"
RUN_SUBDIRS=$("$PYTHON" - "$STATION_JSON5" <<'PY'
import sys
import json5

cfg = json5.load(open(sys.argv[1], encoding="utf-8"))
for run in cfg.get("swp001_runs") or []:
  print(run["work_subdir"])
if not cfg.get("swp001_runs"):
  print(cfg.get("swp001", {}).get("work_subdir", "swp001"))
if "depth_sweep" in cfg:
  print(cfg.get("depth_sweep", {}).get("work_subdir", "depth001"))
PY
)

while IFS= read -r subdir; do
  [ -n "$subdir" ] || continue
  if [ -f "$SCRIPT_DIR/work/$STATION_NAME/$subdir/run_prg.sh" ]; then
    echo "  cd work/$STATION_NAME/$subdir && bash run_prg.sh true 42"
  else
    echo "  after generation: cd work/$STATION_NAME/$subdir && bash run_prg.sh true 42"
  fi
done <<< "$RUN_SUBDIRS"
