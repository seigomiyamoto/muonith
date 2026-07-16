#!/usr/bin/env python3
"""Scaffolding tool that generates a volcano analysis work site from site.json5.

Generates detparams/, depth001/, and swp001/ directory trees with all
configuration files, shell scripts, and templates needed to run the
muonith analysis pipeline.

Usage:
  # New style: config in param_sites/, data in data/, output to work/
  python3 scripts/init_work_site.py param_sites/tarumae_base.json5
  python3 scripts/init_work_site.py param_sites/tarumae_base.json5 --run --verbose

  # Legacy style: config inside work dir (backward compatible)
  python3 scripts/init_work_site.py work/mysite/site.json5
  python3 scripts/init_work_site.py work/mysite/site.json5 --dry-run
"""

import argparse
import logging
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from string import Template

log = logging.getLogger("init_work_site")

# ---------------------------------------------------------------------------
# Template helpers
# ---------------------------------------------------------------------------

_TEMPLATE_DIR = Path(__file__).parent / "templates"


def _load_template(filename: str) -> Template:
  """Load a string.Template from the templates directory."""
  return Template((_TEMPLATE_DIR / filename).read_text(encoding="utf-8"))


# ---------------------------------------------------------------------------
# Template strings
# ---------------------------------------------------------------------------

TPL_MK_COORDINATE_LIST = _load_template("mk_coordinate_list.sh.tpl")

TPL_TEMPLATE_DET = _load_template("template-det.json5.tpl")

TPL_MKJSON5 = _load_template("mkjson5.sh.tpl")

TPL_DEPTH001_RUN = _load_template("depth001_run_prg.sh.tpl")

TPL_DEPTH001_AUTO_PLOT = _load_template("depth001_auto_plot.json5.tpl")

TPL_PRM_RESO = _load_template("prm_reso.json5.tpl")

TPL_REC001_RUN = _load_template("swp001_run_prg.sh.tpl")

TPL_PRM_MUONITH = _load_template("prm_muonith.json5.tpl")

TPL_AUTO_PLOT = _load_template("swp001_auto_plot.json5.tpl")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_file(path: Path, content: str, dry_run: bool, force: bool) -> bool:
  """Write *content* to *path*.  Return True if file was written."""
  if path.exists() and not force:
    log.warning("SKIP (exists): %s  -- use --force to overwrite", path)
    return False
  if dry_run:
    log.info("DRY-RUN would create: %s", path)
    return False
  path.parent.mkdir(parents=True, exist_ok=True)
  path.write_text(content, encoding="utf-8")
  log.info("Created: %s", path)
  return True


def _copy_file(src: Path, dst: Path, dry_run: bool, force: bool) -> bool:
  """Copy binary file *src* -> *dst*.  Return True if file was copied."""
  if dst.exists() and not force:
    log.warning("SKIP (exists): %s  -- use --force to overwrite", dst)
    return False
  if dry_run:
    log.info("DRY-RUN would copy: %s -> %s", src, dst)
    return False
  dst.parent.mkdir(parents=True, exist_ok=True)
  shutil.copy2(src, dst)
  log.info("Copied: %s -> %s", src, dst)
  return True


ANGLE_BIN_KEYS = ("nbinx", "txmin", "txmax", "nbiny", "tymin", "tymax")
ANGLE_BIN_FLOAT_TOL = 5e-5


def _angle_bins_from_detector_params(prm: dict) -> dict:
  return {key: prm[key] for key in ANGLE_BIN_KEYS}


def _build_uniform_eff_table(angle_bins: dict) -> str:
  """Build a unit-efficiency table on the detector angle grid."""
  nbinx = int(angle_bins["nbinx"])
  nbiny = int(angle_bins["nbiny"])
  txmin = float(angle_bins["txmin"])
  txmax = float(angle_bins["txmax"])
  tymin = float(angle_bins["tymin"])
  tymax = float(angle_bins["tymax"])
  if nbinx <= 0 or nbiny <= 0 or txmax <= txmin or tymax <= tymin:
    log.error("Invalid detector angle grid: %s", angle_bins)
    sys.exit(1)
  dx = (txmax - txmin) / nbinx
  dy = (tymax - tymin) / nbiny
  lines = []
  for iy in range(nbiny):
    y0 = tymin + iy * dy
    y1 = y0 + dy
    for ix in range(nbinx):
      x0 = txmin + ix * dx
      x1 = x0 + dx
      lines.append(
        f"{x0: .4f} {x1: .4f} {y0: .4f} {y1: .4f}  "
        "1.0000  1.0000  1.0000"
      )
  return "\n".join(lines) + "\n"


def _write_uniform_eff_table(path: Path, angle_bins: dict, dry_run: bool, force: bool) -> bool:
  """Write the generated unit-efficiency table, replacing stale grids."""
  content = _build_uniform_eff_table(angle_bins)
  if path.exists() and not force:
    current_bins = _read_eff_table_angle_bins(path)
    if _angle_bins_match(angle_bins, current_bins):
      log.warning("SKIP (exists): %s  -- use --force to overwrite", path)
      return False
    if dry_run:
      log.info("DRY-RUN would update stale eff_table grid: %s", path)
      return False
    path.write_text(content, encoding="utf-8")
    log.info("Updated stale eff_table grid: %s", path)
    return True
  return _write_file(path, content, dry_run, force)


def _read_eff_table_angle_bins(path: Path) -> dict:
  """Read the rectangular angle grid from a whitespace eff_table."""
  try:
    xbins = set()
    ybins = set()
    n_rows = 0
    for line in path.read_text(encoding="utf-8").splitlines():
      s = line.strip()
      if not s or s.startswith("#"):
        continue
      cols = s.split()
      if len(cols) < 4:
        continue
      x0, x1, y0, y1 = (float(cols[i]) for i in range(4))
      xbins.add((x0, x1))
      ybins.add((y0, y1))
      n_rows += 1
  except (OSError, ValueError) as exc:
    log.error("Cannot read eff_table angle grid: %s: %s", path, exc)
    sys.exit(1)
  if not xbins or not ybins:
    log.error("Cannot read eff_table angle grid: %s has no bins", path)
    sys.exit(1)
  if n_rows != len(xbins) * len(ybins):
    log.error(
      "eff_table is not a complete rectangular grid: %s has %d rows, "
      "expected %d x %d = %d",
      path, n_rows, len(xbins), len(ybins), len(xbins) * len(ybins),
    )
    sys.exit(1)
  return {
    "nbinx": len(xbins),
    "txmin": min(x0 for x0, _ in xbins),
    "txmax": max(x1 for _, x1 in xbins),
    "nbiny": len(ybins),
    "tymin": min(y0 for y0, _ in ybins),
    "tymax": max(y1 for _, y1 in ybins),
  }


def _angle_bins_match(det_bins: dict, eff_bins: dict) -> bool:
  for key in ("nbinx", "nbiny"):
    if int(det_bins[key]) != int(eff_bins[key]):
      return False
  for key in ("txmin", "txmax", "tymin", "tymax"):
    if abs(float(det_bins[key]) - float(eff_bins[key])) > ANGLE_BIN_FLOAT_TOL:
      return False
  return True


def _resolve_eff_table_path(work_dir: Path, run_subdir: str, path_eff_table: str) -> Path:
  return (work_dir / run_subdir / path_eff_table).resolve()


def _validate_detparams_eff_table_grid(work_dir: Path, cfg: dict, spec_label: str) -> None:
  """Ensure every runcard's path_eff_table has the same angle grid."""
  import json5  # already validated by load_site_config()
  det_dir = work_dir / "detparams" / spec_label
  candidates = sorted(det_dir.glob("det_*.json5")) + sorted(det_dir.glob("runcard_det_*.json5"))
  template = det_dir / "template-det.json5"
  if template.exists():
    candidates.append(template)
  if not candidates:
    log.warning("eff_table grid check skipped: no detector json5 in %s", det_dir)
    return
  run_subdir = _swp_work_subdir(cfg)
  checked = 0
  for det_path in candidates:
    try:
      with det_path.open(encoding="utf-8") as f:
        prm = json5.load(f)["DETECTOR_PARAMETERS"]
    except (OSError, ValueError, KeyError, TypeError) as exc:
      log.error("Cannot read detector runcard for eff_table grid check: %s: %s", det_path, exc)
      sys.exit(1)
    path_eff_table = str(prm.get("path_eff_table", "none"))
    if path_eff_table.lower() == "none":
      continue
    try:
      det_bins = _angle_bins_from_detector_params(prm)
    except KeyError as exc:
      log.error("Detector runcard is missing angle-bin key %s: %s", exc, det_path)
      sys.exit(1)
    eff_path = _resolve_eff_table_path(work_dir, run_subdir, path_eff_table)
    eff_bins = _read_eff_table_angle_bins(eff_path)
    if not _angle_bins_match(det_bins, eff_bins):
      log.error(
        "Detector angle grid and eff_table grid differ: %s vs %s\n"
        "  detector: %s\n  eff_table: %s",
        det_path, eff_path, det_bins, eff_bins,
      )
      sys.exit(1)
    checked += 1
  log.info("eff_table grid check: %d detector file(s) matched %s", checked, spec_label)


def _fmt_num(v):
  """Format a number for JSON5 output: int as int, float as float."""
  if isinstance(v, float):
    if v == int(v) and abs(v) < 1e15:
      # Show as e.g. 2000.0 not 2000
      return f"{v}"
    return f"{v}"
  return str(v)


def _fmt_num_one_decimal(v):
  """Format a numeric plotting limit with one decimal place."""
  return f"{float(v):.1f}"


def _swp_work_subdir(cfg: dict) -> str:
  """Return the swp work subdirectory name from cfg["swp001"]."""
  return _validate_work_subdir(
    str(cfg.get("swp001", {}).get("work_subdir", "swp001")),
    "swp001.work_subdir",
  )


def _validate_work_subdir(subdir: str, field_name: str) -> str:
  """Validate a generated work subdirectory name."""
  p = Path(subdir)
  if p.is_absolute() or p.name != subdir or subdir in ("", ".", ".."):
    log.error("Invalid %s: %s", field_name, subdir)
    sys.exit(1)
  return subdir


def _iter_swp001_run_cfgs(cfg: dict):
  """Yield one cfg per swp001 run, preserving the legacy single-run behavior."""
  runs = cfg.get("swp001_runs")
  if not runs:
    yield cfg
    return
  if not isinstance(runs, list):
    log.error("swp001_runs must be an array")
    sys.exit(1)
  if len(runs) == 0:
    log.error("swp001_runs must contain at least one run")
    sys.exit(1)

  base_rec = cfg["swp001"]
  base_rv = cfg.get("reconst_voxels", {})
  for idx, run in enumerate(runs):
    if not isinstance(run, dict):
      log.error("swp001_runs[%d] must be an object", idx)
      sys.exit(1)
    if "work_subdir" not in run:
      log.error("swp001_runs[%d] is missing work_subdir", idx)
      sys.exit(1)

    rec_overrides = {
      k: v for k, v in run.items()
      if k not in ("name", "description", "reconst_voxels")
    }
    merged_rec = {**base_rec, **rec_overrides}
    merged_rec["work_subdir"] = _validate_work_subdir(
      str(merged_rec["work_subdir"]),
      f"swp001_runs[{idx}].work_subdir",
    )
    rv_overrides = run.get("reconst_voxels", {})
    if not isinstance(rv_overrides, dict):
      log.error("swp001_runs[%d].reconst_voxels must be an object", idx)
      sys.exit(1)
    merged_rv = {**base_rv, **rv_overrides}
    yield {**cfg, "swp001": merged_rec, "reconst_voxels": merged_rv}


def _swp_work_subdirs(cfg: dict) -> list[str]:
  """Return all swp work subdirectories that will be generated."""
  subdirs: list[str] = []
  for run_cfg in _iter_swp001_run_cfgs(cfg):
    subdir = _swp_work_subdir(run_cfg)
    if subdir in subdirs:
      log.error("Duplicate swp001 work_subdir: %s", subdir)
      sys.exit(1)
    subdirs.append(subdir)
  return subdirs


def _depth_work_subdir(cfg: dict) -> str:
  """Return the depth work subdirectory name from cfg["depth_sweep"]."""
  subdir = str(cfg.get("depth_sweep", {}).get("work_subdir", "depth001"))
  p = Path(subdir)
  if p.is_absolute() or p.name != subdir or subdir in ("", ".", ".."):
    log.error("Invalid depth_sweep.work_subdir: %s", subdir)
    sys.exit(1)
  return subdir


# ---------------------------------------------------------------------------
# Config loading and validation
# ---------------------------------------------------------------------------

def load_site_config(path: str) -> dict:
  """Load and return site.json5 as a dict."""
  try:
    import json5
  except ImportError:
    log.error("json5 module not found.  Install it:  pip install json5")
    sys.exit(1)
  p = Path(path)
  if not p.exists():
    log.error("site.json5 not found: %s", p)
    sys.exit(1)
  with open(p, encoding="utf-8") as f:
    return json5.load(f)


def validate_config(cfg: dict) -> None:
  """Validate required keys, types, and value ranges.

  When cfg["skip_detparams"] is True, the EPSG / detector_input /
  detector_spec checks are relaxed because the station uses a non-EPSG
  custom orthogonal CRS and reuses pre-existing work/<name>/detparams/
  contents (KML/CSV/det conversion is bypassed).

  depth_sweep is required only when generate_depth001 will run; its
  presence becomes optional so stations can opt out of depth analysis.
  """
  skip_dp = bool(cfg.get("skip_detparams", False))

  required_top = [
    "station_name", "dem_file",
    "center_x", "center_y", "surface_elevation",
    "plot_zmin", "plot_zmax",
    "swp001",
  ]
  if not skip_dp:
    required_top += ["epsg", "detector_input", "z_pos_offset", "detector_spec"]

  missing = [k for k in required_top if k not in cfg]
  if missing:
    log.error("Missing required top-level keys: %s", ", ".join(missing))
    sys.exit(1)

  if not skip_dp:
    spec = cfg["detector_spec"]
    spec_keys = [
      "label", "length_hori", "length_vert", "length_dept",
      "n_unit", "days", "nbinx", "txmin", "txmax",
      "nbiny", "tymin", "tymax",
    ]
    missing_spec = [k for k in spec_keys if k not in spec]
    if missing_spec:
      log.error("Missing detector_spec keys: %s", ", ".join(missing_spec))
      sys.exit(1)

  if "depth_sweep" in cfg:
    ds = cfg["depth_sweep"]
    ds_keys = [
      "base_density", "vec_delta_density", "signal_noise_amplifiers",
      "depth_step", "obj_size_step",
    ]
    missing_ds = [k for k in ds_keys if k not in ds]
    if missing_ds:
      log.error("Missing depth_sweep keys: %s", ", ".join(missing_ds))
      sys.exit(1)

  rec = cfg["swp001"]
  rec_keys = [
    "BL_max", "PL_max", "g3vox_zmin", "g3vox_zmax", "z_pitch",
    "merge_factor", "z_merge_center", "corr_length", "sigma_rho",
    "PL_hist_max",
  ]
  missing_rec = [k for k in rec_keys if k not in rec]
  if missing_rec:
    log.error("Missing swp001 keys: %s", ", ".join(missing_rec))
    sys.exit(1)

  # Value range checks (skipped when using custom CRS)
  if not skip_dp:
    if cfg["epsg"] < 1000:
      log.error("epsg code looks invalid: %d", cfg["epsg"])
      sys.exit(1)


def convert_center_if_degree(cfg: dict) -> None:
  """Convert center_x/y from lon/lat to projected CRS if center_unit is 'degree'."""
  raw_unit = cfg.get("center_unit", "meter")
  unit = raw_unit.strip().lower()
  if unit in ("meter", "meters", "m"):
    return
  if unit in ("km", "kilometer", "kilometers"):
    cfg["center_x"] = cfg["center_x"] * 1000.0
    cfg["center_y"] = cfg["center_y"] * 1000.0
    log.info("Converted center: km -> m (%s, %s)", cfg["center_x"], cfg["center_y"])
    return
  if unit not in ("degree", "degrees", "deg"):
    log.error("center_unit must be 'degree'/'meter'/'km', got '%s'", raw_unit)
    sys.exit(1)
  if cfg.get("skip_detparams", False):
    log.error(
      "center_unit='degree' is incompatible with skip_detparams=true "
      "(no EPSG to transform to). Use 'meter' and supply local coordinates."
    )
    sys.exit(1)
  try:
    from pyproj import Transformer
  except ImportError:
    log.error("pyproj is required for center_unit='degree'.  pip install pyproj")
    sys.exit(1)
  epsg = cfg["epsg"]
  transformer = Transformer.from_crs("EPSG:4326", f"EPSG:{epsg}", always_xy=True)
  lon, lat = cfg["center_x"], cfg["center_y"]
  x, y = transformer.transform(lon, lat)
  log.info("Converted center: (%s, %s) degree -> (%s, %s) EPSG:%d", lon, lat, x, y, epsg)
  cfg["center_x"] = x
  cfg["center_y"] = y


# ---------------------------------------------------------------------------
# Path resolution (param_sites/ + data/ support)
# ---------------------------------------------------------------------------

def _find_repo_root(start: Path) -> Path:
  """Return the repository root, i.e. the parent directory of scripts/.

  This script always lives in <repo root>/scripts/, so the root is derived from
  its own location rather than from a .git directory: a release downloaded as a
  ZIP or tarball has no .git, and the lookup would then fail.  *start* is kept
  for call-site compatibility and is not used.
  """
  del start
  return Path(__file__).resolve().parent.parent


def resolve_work_dir(site_json5_path: Path) -> Path:
  """Determine the work directory from the site.json5 path.

  - If site.json5 is under param_sites/ (e.g. param_sites/tarumae_base.json5),
    return work/{stem}/ under the repo root.
  - Otherwise (legacy), return the parent directory of site.json5.
  """
  resolved = site_json5_path.resolve()
  parent_name = resolved.parent.name
  if parent_name == "param_sites":
    repo_root = resolved.parent.parent
    site_name = resolved.stem
    return repo_root / "work" / site_name
  return resolved.parent


def resolve_data_dir(site_name: str, repo_root: Path) -> Path | None:
  """Return data/{site_name}/ if it exists, else None."""
  data_dir = repo_root / "data" / site_name
  if data_dir.is_dir():
    return data_dir
  return None


def count_detectors_in_kml(kml_path: Path) -> int:
  """Count detector Placemarks in a KML file.

  Supports both GSI (det_00, det_01, ...) and Google Earth (det00, det01, ...)
  naming conventions.  Searches all descendants so Folder-nested layouts work too.
  """
  import re
  import xml.etree.ElementTree as ET
  ns = {"kml": "http://www.opengis.net/kml/2.2"}
  det_pattern = re.compile(r"^det_?\d{2,}$")
  tree = ET.parse(kml_path)
  count = 0
  for pm in tree.getroot().findall(".//kml:Placemark", ns):
    name_el = pm.find("kml:name", ns)
    if name_el is not None and name_el.text and det_pattern.match(name_el.text.strip()):
      count += 1
  if count == 0:
    log.error("No detector Placemarks (det_XX / detXX) found in %s", kml_path)
    sys.exit(1)
  return count


def _symlink_file(src: Path, dst: Path, dry_run: bool, force: bool) -> bool:
  """Create a symlink dst -> src.  Return True if link was created.

  The link target is stored as a path relative to dst's parent so that the
  link remains valid when the entire repository tree is moved (e.g. another
  user clones it under a different absolute path).
  """
  if dst.exists() or dst.is_symlink():
    if not force:
      log.warning("SKIP (exists): %s  -- use --force to overwrite", dst)
      return False
    dst.unlink()
  if dry_run:
    log.info("DRY-RUN would symlink: %s -> %s", dst, src)
    return False
  dst.parent.mkdir(parents=True, exist_ok=True)
  import os
  rel = os.path.relpath(src.resolve(), start=dst.parent.resolve())
  dst.symlink_to(rel)
  log.info("Symlinked: %s -> %s", dst, rel)
  return True


def link_input_files(
  data_dir: Path, work_dir: Path, cfg: dict,
  dry_run: bool, force: bool,
) -> None:
  """Symlink DEM and detector input files from data/ to work/.

  When skip_detparams=true and detector_input is absent, the KML/GeoJSON
  symlink step is skipped (existing detparams/<label>/ is reused).
  """
  dem_basename = Path(cfg["dem_file"]).name

  # DEM -> work_dir/dem/
  dem_src = data_dir / dem_basename
  dem_dst = work_dir / "dem" / dem_basename
  if dem_src.exists():
    _symlink_file(dem_src, dem_dst, dry_run, force)
  else:
    log.warning("DEM not found in data dir: %s", dem_src)

  # KML/GeoJSON -> work_dir/detparams/  (skipped when skip_detparams=true)
  det_input = cfg.get("detector_input", "")
  if not det_input:
    return
  det_basename = Path(det_input).name
  det_src = data_dir / det_basename
  det_dst = work_dir / "detparams" / det_basename
  if det_src.exists():
    _symlink_file(det_src, det_dst, dry_run, force)
  else:
    log.warning("Detector input not found in data dir: %s", det_src)


# ---------------------------------------------------------------------------
# Directory creation
# ---------------------------------------------------------------------------

def create_directories(work_dir: Path, cfg: dict) -> None:
  """Create the directory skeleton under work_dir.

  When skip_detparams=true the per-label directory is NOT pre-created here;
  it is materialised later as a symlink to data/<site>/detparams/<label>/ by
  link_detparams_for_skip().  Only the parent detparams/ container is made.
  """
  skip_dp = bool(cfg.get("skip_detparams", False))
  depth_work_subdir = _depth_work_subdir(cfg)
  dirs = [
    work_dir / depth_work_subdir,
    work_dir / "dem",
    work_dir / "detparams",
  ]
  dirs.extend(work_dir / subdir for subdir in _swp_work_subdirs(cfg))
  spec = cfg.get("detector_spec")
  if spec and "label" in spec and not skip_dp:
    dirs.insert(0, work_dir / "detparams" / spec["label"])
  for d in dirs:
    d.mkdir(parents=True, exist_ok=True)
    log.info("Directory: %s", d)


def link_detparams_for_skip(
  repo_root: Path, src_detparams_root: Path, work_dir: Path, spec_label: str,
  dry_run: bool, force: bool,
  renumber_runcards: list[str] | None = None,
  detector_overrides: dict | None = None,
) -> None:
  """Create relative symlinks under work_dir/detparams/ for skip_detparams flow.

  Default mode (renumber_runcards is None, detector_overrides empty):
    - work/<site>/detparams/<label>/   -> data/<site>/detparams/<label>/  (dir symlink)
    - work/<site>/detparams/eff_table/ -> data/_shared/eff_table/         (dir symlink)

  Renumber mode (renumber_runcards is a list of original runcard filenames
  sorted alphanumerically, detector_overrides empty):
    - work/<site>/detparams/<label>/  is created as a real directory
      containing per-file symlinks det_00.json5 -> ../../../../data/...,
      det_01.json5 -> ..., one per entry in renumber_runcards.  This lets
      the generated prm reference legacy det_NN.json5 names while keeping
      the canonical runcard_det_*.json5 storage in data/.
    - work/<site>/detparams/eff_table/ -> data/_shared/eff_table/         (dir symlink)

  Override mode (detector_overrides is a non-empty dict, e.g.
  {"tymin": 0.0, "tymax": 1.6}):
    - work/<site>/detparams/<label>/  is created as a real directory
      containing per-file *copies* with the given top-level
      DETECTOR_PARAMETERS keys overridden via in-place regex replacement
      (preserving comments and formatting).  data/ is never modified.
    - When combined with renumber_runcards, copies are written as
      det_00.json5..det_NN.json5; otherwise the original filenames are kept.
    - work/<site>/detparams/eff_table/ -> data/_shared/eff_table/         (dir symlink)
  """
  shared_eff = repo_root / "data" / "_shared" / "eff_table"
  src_label = src_detparams_root / spec_label
  dst_label = work_dir / "detparams" / spec_label
  dst_eff = work_dir / "detparams" / "eff_table"

  if not src_label.is_dir():
    log.error(
      "skip_detparams: source dir not found: %s "
      "(expected data/<site>/detparams/<label>/ or work/<site>/detparams/<label>/)",
      src_label,
    )
    sys.exit(1)
  if not shared_eff.is_dir():
    log.error(
      "skip_detparams: shared eff_table not found: %s "
      "(expected data/_shared/eff_table/)", shared_eff,
    )
    sys.exit(1)

  # Helper to clear a destination path consistently across modes.
  def _clear(dst: Path) -> bool:
    if dst.is_symlink() or dst.exists():
      if force:
        if dry_run:
          log.info("DRY-RUN: would remove %s", dst)
          return True
        if dst.is_symlink() or dst.is_file():
          dst.unlink()
        else:
          import shutil
          shutil.rmtree(dst)
        return True
      log.warning("SKIP (exists): %s  -- use --force to overwrite", dst)
      return False
    return True

  overrides_active = bool(detector_overrides)
  same_label_dir = src_label.resolve() == dst_label.resolve()

  def _apply_overrides(src_text: str, ov: dict) -> str:
    """Replace top-level DETECTOR_PARAMETERS values via line-anchored regex.

    Matches lines like  , "tymin": 0.0  // optional comment
    and rewrites only the value, preserving leading whitespace, comma,
    quotes and trailing comment.  Numbers are formatted via json.dumps
    (so 0 -> '0', 1.6 -> '1.6', true -> 'true').
    """
    import json as _json
    out = src_text
    for k, v in ov.items():
      pattern = re.compile(
        rf'^(\s*,?\s*"{re.escape(k)}"\s*:\s*)([^,/\n]+?)(\s*(?://.*)?)$',
        flags=re.MULTILINE,
      )
      def _sub(m, _v=v):
        return f"{m.group(1)}{_json.dumps(_v)}{m.group(3)}"
      out_new, n = pattern.subn(_sub, out)
      if n == 0:
        log.warning("detector_overrides: key %r not found in runcard", k)
      out = out_new
    return out

  renumber_active = renumber_runcards is not None

  # When overrides_active without explicit renumber, enumerate runcards
  # in src_label and keep their original names.
  if not renumber_active and overrides_active:
    runcards = sorted(p.name for p in src_label.glob("runcard_det_*.json5"))
    if not runcards:
      runcards = sorted(p.name for p in src_label.glob("det_*.json5"))
  else:
    runcards = renumber_runcards or []

  # ---- label dir: dir-symlink OR real dir + per-file symlinks/copies ----
  if same_label_dir:
    if renumber_active:
      log.error(
        "skip_detparams: cannot renumber detector files in-place: %s",
        dst_label,
      )
      sys.exit(1)
    if overrides_active:
      override_targets = [dst_label / rc_name for rc_name in runcards]
      template_path = dst_label / "template-det.json5"
      if template_path.exists():
        override_targets.append(template_path)
      for rc_path in override_targets:
        if dry_run:
          log.info(
            "DRY-RUN: would update %s in-place with overrides %s",
            rc_path, detector_overrides,
          )
        else:
          new_text = _apply_overrides(
            rc_path.read_text(encoding="utf-8"), detector_overrides,
          )
          rc_path.write_text(new_text, encoding="utf-8")
          log.info("Updated in-place: %s (overrides %s)", rc_path, detector_overrides)
    else:
      log.info("skip_detparams: reusing existing detparams directory: %s", dst_label)
  else:
    if not renumber_active and not overrides_active:
      # mode 1: full dir-symlink (no per-file work)
      if _clear(dst_label):
        rel = os.path.relpath(src_label.resolve(), start=dst_label.parent.resolve())
        if dry_run:
          log.info("DRY-RUN: would symlink %s -> %s", dst_label, rel)
        else:
          dst_label.symlink_to(rel, target_is_directory=True)
          log.info("Symlinked: %s -> %s", dst_label, rel)
    else:
      # mode 2/3: real dir + per-file entries
      if _clear(dst_label):
        if dry_run:
          log.info("DRY-RUN: would mkdir %s", dst_label)
        else:
          dst_label.mkdir(parents=True, exist_ok=True)
          log.info("Created dir: %s", dst_label)

        for i, rc_name in enumerate(runcards):
          # Filename in dst: det_NN.json5 if renumber_active, else original rc_name.
          det_name = f"det_{i:02d}.json5" if renumber_active else rc_name
          link_path = dst_label / det_name
          src_path = src_label / rc_name

          if overrides_active:
            # Copy with in-place key overrides; keeps comments intact.
            if dry_run:
              log.info(
                "DRY-RUN: would copy %s -> %s with overrides %s",
                src_path, link_path, detector_overrides,
              )
            else:
              new_text = _apply_overrides(
                src_path.read_text(encoding="utf-8"), detector_overrides,
              )
              link_path.write_text(new_text, encoding="utf-8")
              log.info(
                "Wrote copy: %s (overrides %s, source %s)",
                link_path, detector_overrides, rc_name,
              )
          else:
            # Per-file symlink (renumber-only, no edits)
            target_abs = src_path.resolve()
            rel = os.path.relpath(target_abs, start=link_path.parent.resolve())
            if dry_run:
              log.info("DRY-RUN: would symlink %s -> %s (was %s)", link_path, rel, rc_name)
            else:
              if link_path.is_symlink() or link_path.exists():
                link_path.unlink()
              link_path.symlink_to(rel)
              log.info("Symlinked: %s -> %s (was %s)", link_path, rel, rc_name)

  # ---- eff_table: always dir-symlink ----
  if _clear(dst_eff):
    rel = os.path.relpath(shared_eff.resolve(), start=dst_eff.parent.resolve())
    if dry_run:
      log.info("DRY-RUN: would symlink %s -> %s", dst_eff, rel)
    else:
      dst_eff.symlink_to(rel, target_is_directory=True)
      log.info("Symlinked: %s -> %s", dst_eff, rel)

  # ---- eff_model: dir-symlink when the site provides one (optional) ----
  # Sites that set path_eff_model in their runcards keep the analytic
  # efficiency model json5 under data/<site>/detparams/eff_model/; mirror it
  # into the work dir the same way as eff_table so the relative path
  # ../detparams/eff_model/<name>.json5 resolves from the swp run directory.
  src_eff_model = src_detparams_root / "eff_model"
  dst_eff_model = work_dir / "detparams" / "eff_model"
  if src_eff_model.is_dir():
    if src_eff_model.resolve() == dst_eff_model.resolve():
      log.info("skip_detparams: reusing existing eff_model directory: %s", dst_eff_model)
    elif _clear(dst_eff_model):
      rel = os.path.relpath(src_eff_model.resolve(), start=dst_eff_model.parent.resolve())
      if dry_run:
        log.info("DRY-RUN: would symlink %s -> %s", dst_eff_model, rel)
      else:
        dst_eff_model.symlink_to(rel, target_is_directory=True)
        log.info("Symlinked: %s -> %s", dst_eff_model, rel)


# ---------------------------------------------------------------------------
# Detector files list builder
# ---------------------------------------------------------------------------

def _build_det_files_block(
  n_detectors: int,
  prefix: str,
  filenames: list[str] | None = None,
) -> str:
  """Build the det_files JSON array content (standard JSON5, 2-space indent).

  Inserted at nesting level 3 (root > section > det_files), so 6-space indent.
  prefix is the path prefix, e.g. '../detparams/0.2msq120days/'.
  If filenames is provided (skip_detparams=True path), entries use those
  literal names; otherwise det_{NN}.json5 is generated for legacy pipeline.
  """
  indent = " " * 6
  if filenames is None:
    names = [f"det_{i:02d}.json5" for i in range(n_detectors)]
  else:
    names = list(filenames)
  lines = []
  for i, fn in enumerate(names):
    trailing = "," if i < len(names) - 1 else ""
    lines.append(f'{indent}"{prefix}{fn}"{trailing}')
  return "\n".join(lines)


def _build_signal_noise_block(amplifiers: list) -> str:
  """Build signal_noise_amplifiers array lines (standard JSON5, 2-space indent).

  Inserted at nesting level 4, so 8-space indent.
  """
  n = len(amplifiers)
  indent = " " * 8
  lines = []
  for i, pair in enumerate(amplifiers):
    a, b = pair
    trailing = "," if i < n - 1 else ""
    lines.append(f"{indent}[{a}, {b}]{trailing}")
  return "\n".join(lines)


def _det_list_name(spec_label: str) -> str:
  """Convert spec label like '0.2msq120days' to list name like '02msq120days'."""
  return spec_label.replace(".", "")


# ---------------------------------------------------------------------------
# Generator: detparams
# ---------------------------------------------------------------------------

def generate_detparams(work_dir: Path, cfg: dict, dry_run: bool, force: bool) -> None:
  """Generate detparams/ files."""
  site_name = cfg["station_name"]
  epsg = cfg["epsg"]
  spec = cfg["detector_spec"]
  spec_label = spec["label"]
  det_input = cfg["detector_input"]
  det_input_basename = Path(det_input).name
  dem_file = cfg["dem_file"]
  dem_basename = Path(dem_file).name

  # Determine script type based on extension
  ext = Path(det_input).suffix.lower()
  if ext == ".kml":
    kml_or_geojson_script = "kml_to_csv.py"
    input_flag = "kml_in"
  elif ext == ".geojson":
    kml_or_geojson_script = "geojson_to_csv.py"
    input_flag = "geojson_in"
  else:
    log.warning("Unknown detector_input extension '%s', defaulting to kml_to_csv.py", ext)
    kml_or_geojson_script = "kml_to_csv.py"
    input_flag = "kml_in"

  # DEM relative path from detparams/ to dem/
  dem_relative_path = f"../dem/{dem_basename}"

  # 1) mk_coordinate_list.sh
  content = TPL_MK_COORDINATE_LIST.substitute(
    epsg=epsg,
    kml_or_geojson_script=kml_or_geojson_script,
    input_flag=input_flag,
    detector_input_basename=det_input_basename,
    site_name=site_name,
    dem_relative_path=dem_relative_path,
  )
  _write_file(work_dir / "detparams" / "mk_coordinate_list.sh", content, dry_run, force)

  # 2) template-det.json5
  content = TPL_TEMPLATE_DET.substitute(
    nbinx=spec["nbinx"],
    txmin=spec["txmin"],
    txmax=spec["txmax"],
    nbiny=spec["nbiny"],
    tymin=spec["tymin"],
    tymax=spec["tymax"],
    length_hori=spec["length_hori"],
    length_vert=spec["length_vert"],
    length_dept=spec["length_dept"],
    n_unit=spec["n_unit"],
    days=_fmt_num(float(spec["days"])),
    spec_label=spec_label,
  )
  _write_file(
    work_dir / "detparams" / spec_label / "template-det.json5",
    content, dry_run, force,
  )

  # 3) mkjson5.sh
  csv_basename = f"{site_name}_table_epsg{epsg}_gsi.csv"
  content = TPL_MKJSON5.substitute(
    work_dir_name=work_dir.name,
    spec_label=spec_label,
    csv_basename=csv_basename,
    z_pos_offset=_fmt_num(cfg["z_pos_offset"]),
  )
  _write_file(
    work_dir / "detparams" / spec_label / "mkjson5.sh",
    content, dry_run, force,
  )

  # 4) eff_table_sample00.tmp -- generated on the detector angle grid
  eff_dst = work_dir / "detparams" / spec_label / "eff_table_sample00.tmp"
  _write_uniform_eff_table(eff_dst, spec, dry_run, force)


# ---------------------------------------------------------------------------
# Generator: depth001
# ---------------------------------------------------------------------------

def generate_depth001(work_dir: Path, cfg: dict, n_detectors: int, dry_run: bool, force: bool) -> None:
  """Generate depth analysis files."""
  site_name = cfg["station_name"]
  site_name_cap = site_name.replace("_", " ").title()
  epsg = cfg["epsg"]
  spec = cfg["detector_spec"]
  spec_label = spec["label"]
  ds = cfg["depth_sweep"]
  depth_work_subdir = _depth_work_subdir(cfg)

  det_files_prefix = f"../detparams/{spec_label}/"
  det_files_block = _build_det_files_block(
    n_detectors, det_files_prefix, cfg.get("_det_filenames_override"),
  )
  det_list_name = _det_list_name(spec_label)

  base_density = ds["base_density"]
  base_density_int = int(base_density)

  # Angle bin override
  tf_override = ds.get("tf_override_angle_bin", False)
  abo = ds.get("angle_bin_override", {})

  sna_block = _build_signal_noise_block(ds["signal_noise_amplifiers"])

  # Format vec_delta_density as JSON array string
  vdd = ds["vec_delta_density"]
  vdd_str = "[" + ", ".join(str(v) for v in vdd) + "]"

  # 1) run_prg.sh
  content = TPL_DEPTH001_RUN.substitute(
    work_dir_name=work_dir.name,
    depth_work_subdir=depth_work_subdir,
    site_name=site_name,
    site_name_cap=site_name_cap,
    epsg=epsg,
    center_x=cfg["center_x"],
    center_y=cfg["center_y"],
    plot_zmin=cfg["plot_zmin"],
    plot_zmax=cfg["plot_zmax"],
    surface_elevation=cfg["surface_elevation"],
  )
  _write_file(work_dir / depth_work_subdir / "run_prg.sh", content, dry_run, force)

  # 2) prm_reso.json5
  dem_path = f"../dem/{Path(cfg['dem_file']).name}"
  path_length_bl_max = ds.get("BL_max")
  bl_max_line = ""
  if path_length_bl_max is not None:
    bl_max_line = f'    "BL_max": {_fmt_num(path_length_bl_max)},'
  tf_out_det_pl_signal_line = '      , "tf_out_det_PL_signal": true  // true = output base PL/signal distribution'
  if "tf_out_det_PL_signal" in ds:
    tf_out_det_pl_signal = ds["tf_out_det_PL_signal"]
    if tf_out_det_pl_signal is None:
      tf_out_det_pl_signal_line = ""
    else:
      tf_out_det_pl_signal_line = (
        "      , \"tf_out_det_PL_signal\": "
        f"{str(bool(tf_out_det_pl_signal)).lower()}  // true = output base PL/signal distribution"
      )
  content = TPL_PRM_RESO.substitute(
    work_dir_name=work_dir.name,
    depth_work_subdir=depth_work_subdir,
    flux_groom=cfg.get("flux_groom", "daemon_groom"),
    site_name=site_name,
    site_name_cap=site_name_cap,
    epsg=epsg,
    det_list_name=det_list_name,
    det_files_block=det_files_block,
    base_density_int=base_density_int,
    dem_path=dem_path,
    grid2d_pillar_name=ds.get("grid2d_pillar_name", f"{site_name}01"),
    grid2d_initial_uniform_density=_fmt_num(
      ds.get("grid2d_initial_uniform_density", base_density_int),
    ),
    x_cnt_obj=_fmt_num(ds.get("x_cnt_obj", cfg["center_x"])),
    y_cnt_obj=_fmt_num(ds.get("y_cnt_obj", cfg["center_y"])),
    base_density=_fmt_num(base_density),
    obj_size_upper_limit=_fmt_num(ds["obj_size_upper_limit"]),
    obj_size_lower_limit=_fmt_num(ds["obj_size_lower_limit"]),
    elev_center_step=_fmt_num(ds.get("elev_center_step", 0.0025)),
    BL_max_line=bl_max_line,
    tf_out_det_PL_signal_line=tf_out_det_pl_signal_line,
    vec_delta_density=vdd_str,
    tf_override_angle_bin=str(tf_override).lower(),
    abo_nbinx=abo.get("nbinx", 640),
    abo_txmin=abo.get("txmin", -0.4),
    abo_txmax=abo.get("txmax", 0.4),
    abo_nbiny=abo.get("nbiny", 640),
    abo_tymin=abo.get("tymin", 0.0),
    abo_tymax=abo.get("tymax", 0.8),
    signal_noise_amplifiers_block=sna_block,
    obj_size_step=_fmt_num(ds["obj_size_step"]),
    depth_step=_fmt_num(ds["depth_step"]),
  )
  _write_file(work_dir / depth_work_subdir / "prm_reso.json5", content, dry_run, force)

  # 3) heatmap_config.json5 -- copy from scripts/templates/, then override
  script_dir = Path(__file__).resolve().parent
  hm_src = script_dir / "templates" / "heatmap_config.json5"
  hm_dst = work_dir / depth_work_subdir / "heatmap_config.json5"
  if hm_src.exists():
    _copy_file(hm_src, hm_dst, dry_run, force)
    # Override size_max / depth_max from station config if present
    if not dry_run and hm_dst.exists():
      size_max = ds.get("size_max")
      depth_max = ds.get("depth_max")
      contour_nan_margin = ds.get("contour_nan_margin")
      if size_max is not None or depth_max is not None or contour_nan_margin is not None:
        text = hm_dst.read_text()
        if size_max is not None:
          text = re.sub(r'"size_max"\s*:\s*\d+', f'"size_max": {size_max}', text)
        if depth_max is not None:
          text = re.sub(r'"depth_max"\s*:\s*\d+', f'"depth_max": {depth_max}', text)
        if contour_nan_margin is not None:
          text = re.sub(
            r'"contour_nan_margin"\s*:\s*\d+',
            f'"contour_nan_margin": {contour_nan_margin}',
            text,
          )
        hm_dst.write_text(text)
        log.info(
          "heatmap_config.json5: size_max=%s, depth_max=%s, contour_nan_margin=%s",
          size_max, depth_max, contour_nan_margin,
        )
  else:
    log.warning("Template not found: %s  -- skipping heatmap_config copy", hm_src)

  # 4) auto_plot.json5
  pl_hist_max = cfg.get("swp001", {}).get("PL_hist_max", 2000)
  content = TPL_DEPTH001_AUTO_PLOT.substitute(
    site_name_cap=site_name_cap,
    PL_hist_max=pl_hist_max,
  )
  _write_file(work_dir / depth_work_subdir / "auto_plot.json5", content, dry_run, force)


# ---------------------------------------------------------------------------
# Generator: swp001
# ---------------------------------------------------------------------------

def _build_shell_density_lines(rec: dict, base_density_int: int) -> str:
  """Build optional shell-density override lines for NAGAINV_PARAMETERS."""
  keys = [
    ("shell_density_upper", "upper"),
    ("shell_density_lower", "lower"),
    ("shell_density_lateral", "lateral"),
  ]
  if any(key in rec for key, _ in keys):
    return "\n".join(
      f'      , "{key}": {_fmt_num(rec.get(key, base_density_int))}   '
      f'// {label} shell prior density [kg/m³]'
      for key, label in keys
    )
  return "\n".join(
    f'      // , "{key}": {base_density_int}   '
    f'// Optional: {label} shell density [kg/m³]'
    for key, label in keys
  )


def _make_eff_case_detparams(
  work_dir: Path, cfg: dict, case_label: str, eff_table_basename: str,
  dry_run: bool, force: bool,
) -> None:
  """Create a per-case detparams dir for the skip_detparams flow.

  Copies the runcards into work/<site>/detparams/<case_label>/ while swapping
  path_eff_table and disabling path_eff_model.  The shared eff_table symlink
  (created by link_detparams_for_skip) is reused: path_eff_table is resolved
  relative to the swp run directory, so the variant dir holds the det_NN.json5
  copies only.  No source file under data/ is modified.
  """
  src_root = cfg["_det_source_root"]
  src_label = cfg["detector_spec"]["label"]
  src_dir = src_root / src_label
  runcard_names = cfg.get("_det_runcard_names") or []
  renumber = bool(cfg.get("detparams_renumber_to_det_nn", False))
  dst_label = work_dir / "detparams" / case_label
  new_path = f"../detparams/eff_table/{eff_table_basename}"
  # path_eff_table value contains slashes, so the generic detector_overrides
  # regex (which excludes '/') cannot match it; use a string-value-aware pattern.
  pat = re.compile(r'^(\s*,?\s*"path_eff_table"\s*:\s*")[^"\n]*(")', flags=re.MULTILINE)
  # An eff-case run sweeps eff_table files, so the analytic efficiency model
  # (which takes precedence over path_eff_table in C++) must be disabled in
  # the per-case copies.
  pat_eff_model = re.compile(r'^(\s*,?\s*"path_eff_model"\s*:\s*")[^"\n]*(")', flags=re.MULTILINE)

  if dst_label.exists() and not force:
    log.warning("SKIP (exists): %s  -- use --force to overwrite", dst_label)
    return
  if dry_run:
    log.info(
      "DRY-RUN: would create eff-case detparams %s (path_eff_table -> %s, %d runcards)",
      dst_label, new_path, len(runcard_names),
    )
    return
  if dst_label.exists():
    import shutil
    shutil.rmtree(dst_label)
  dst_label.mkdir(parents=True, exist_ok=True)
  for i, rc_name in enumerate(runcard_names):
    det_name = f"det_{i:02d}.json5" if renumber else rc_name
    src_text = (src_dir / rc_name).read_text(encoding="utf-8")
    new_text, n = pat.subn(rf'\g<1>{new_path}\g<2>', src_text)
    if n == 0:
      log.warning("eff-case detparams: path_eff_table not found in %s", rc_name)
    new_text = pat_eff_model.sub(r'\g<1>none\g<2>', new_text)
    (dst_label / det_name).write_text(new_text, encoding="utf-8")
  log.info(
    "Created eff-case detparams: %s (%d runcards, eff_table=%s)",
    dst_label, len(runcard_names), eff_table_basename,
  )
  _validate_detparams_eff_table_grid(work_dir, cfg, case_label)


def _make_angle_case_detparams(
  work_dir: Path, cfg: dict, case_label: str,
  dry_run: bool, force: bool,
) -> None:
  """Create a degree-unit detparams dir for the skip_detparams flow.

  Copies the materialized runcards from work/<site>/detparams/<label>/
  (detector_overrides already applied by link_detparams_for_skip) into
  work/<site>/detparams/<case_label>/, rewriting angle_unit to "degree"
  and converting txmin/txmax/tymin/tymax from tangent to degree
  (atan(v) * 180 / pi).  nbinx/nbiny are kept as-is.  No source file
  under data/ is modified.
  """
  src_label = cfg["detector_spec"]["label"]
  src_dir = work_dir / "detparams" / src_label
  det_names = cfg.get("_det_filenames_override") or []
  dst_label = work_dir / "detparams" / case_label

  pat_unit = re.compile(r'^(\s*,?\s*"angle_unit"\s*:\s*")[^"\n]*(")', flags=re.MULTILINE)
  pat_range = {
    key: re.compile(rf'^(\s*,?\s*"{key}"\s*:\s*)([^,/\n]+?)(\s*(?://.*)?)$', flags=re.MULTILINE)
    for key in ("txmin", "txmax", "tymin", "tymax")
  }

  if dst_label.exists() and not force:
    log.warning("SKIP (exists): %s  -- use --force to overwrite", dst_label)
    return
  if dry_run:
    log.info(
      "DRY-RUN: would create degree-case detparams %s (%d runcards)",
      dst_label, len(det_names),
    )
    return
  if dst_label.exists():
    shutil.rmtree(dst_label)
  dst_label.mkdir(parents=True, exist_ok=True)
  for det_name in det_names:
    src_text = (src_dir / det_name).read_text(encoding="utf-8")
    new_text, n = pat_unit.subn(r'\g<1>degree\g<2>', src_text)
    if n == 0:
      log.warning("degree-case detparams: angle_unit not found in %s", det_name)
    for key, pat in pat_range.items():
      m = pat.search(new_text)
      if m is None:
        log.warning("degree-case detparams: %s not found in %s", key, det_name)
        continue
      deg_value = math.degrees(math.atan(float(m.group(2))))
      new_text = pat.sub(rf'\g<1>{deg_value:.6g}\g<3>', new_text, count=1)
    (dst_label / det_name).write_text(new_text, encoding="utf-8")
  log.info(
    "Created degree-case detparams: %s (%d runcards, tangent -> degree)",
    dst_label, len(det_names),
  )


def _read_runcard_angle_bins(det_dir: Path) -> dict | None:
  """Read the angle-axis frame from the first det_*.json5 in det_dir.

  Returns {nbinx, txmin, txmax, nbiny, tymin, tymax} taken from the
  DETECTOR_PARAMETERS block, or None when no runcard is readable (e.g.
  --dry-run, or the non-skip_detparams flow where runcards are generated
  later by run_pipeline) so the caller can fall back to detector_spec.
  """
  import json5  # already validated by load_site_config()
  candidates = sorted(det_dir.glob("det_*.json5")) if det_dir.is_dir() else []
  if not candidates:
    return None
  path = candidates[0]
  try:
    with path.open(encoding="utf-8") as f:
      prm = json5.load(f)["DETECTOR_PARAMETERS"]
    return {key: prm[key] for key in ("nbinx", "txmin", "txmax", "nbiny", "tymin", "tymax")}
  except (OSError, ValueError, KeyError, TypeError) as exc:
    log.warning("auto_plot angle bins: cannot read %s: %s", path, exc)
    return None


def _build_dem_gen_block(cfg: dict) -> str:
  """Build the run_prg.sh DEM-generation snippet for synthetic demo stations.

  Opt-in via the ``dem_generator`` site-config key.  Returns an empty string for
  normal stations so their run_prg.sh stays byte-identical.  When present, the
  snippet regenerates the (git-ignored) synthetic DEM into ../dem/ before the
  detector plot and the executable run, calling
  param_sites/tutorial/make_primitive_dem.py so the binary never needs to be
  committed.
  """
  gen = cfg.get("dem_generator")
  if not gen:
    return ""
  dem_basename = Path(cfg["dem_file"]).name
  script_rel = "../../../param_sites/tutorial/make_primitive_dem.py"
  groups = [f"--shape {gen['shape']}"]
  groups.append(
    f"--x_min {_fmt_num(gen['x_min'])} --x_max {_fmt_num(gen['x_max'])}"
    f" --y_min {_fmt_num(gen['y_min'])} --y_max {_fmt_num(gen['y_max'])}"
  )
  groups.append(
    f"--x_interval {_fmt_num(gen['x_interval'])}"
    f" --y_interval {_fmt_num(gen['y_interval'])}"
  )
  groups.append(
    f"--xcnt {_fmt_num(gen['xcnt'])} --ycnt {_fmt_num(gen['ycnt'])}"
    f" --base {_fmt_num(gen['base'])} --height {_fmt_num(gen['height'])}"
  )
  # Size basis: --width (unified footprint) plus any explicit per-shape override.
  size_tokens = []
  if "width" in gen:
    size_tokens.append(f"--width {_fmt_num(gen['width'])}")
  for opt in ("radius", "half_width", "sigma_x", "sigma_y"):
    if opt in gen:
      size_tokens.append(f"--{opt} {_fmt_num(gen[opt])}")
  if size_tokens:
    groups.append(" ".join(size_tokens))
  if gen.get("float64"):
    groups.append("--float64")
  groups.append(f"--outbin ../dem/{dem_basename}")
  cmd = "uv run python " + script_rel + " \\\n  " + " \\\n  ".join(groups)
  header = (
    "# Regenerate the synthetic demo DEM before the run (not committed to the\n"
    "# repo; driven by the dem_generator key in the site config).\n"
    "mkdir -p ../dem"
  )
  return f"{header}\n{cmd}\n\n"


def generate_swp001(work_dir: Path, cfg: dict, n_detectors: int, dry_run: bool, force: bool) -> None:
  """Generate swp001/ files."""
  site_name = cfg["station_name"]
  site_name_cap = site_name.replace("_", " ").title()
  epsg = cfg["epsg"]
  spec = cfg["detector_spec"]
  spec_label = spec["label"]
  rec = cfg["swp001"]
  swp_work_subdir = _swp_work_subdir(cfg)

  # Per-run efficiency case: when eff_table_basename is set (skip_detparams flow),
  # build a per-case detparams variant with path_eff_table swapped and point this
  # run at it; otherwise reuse the shared detparams dir (backward compatible).
  eff_table_basename = rec.get("eff_table_basename")
  # Per-run angle-unit case: detparams_angle_unit="degree" builds a degree-unit
  # detparams variant (see _make_angle_case_detparams).  skip_detparams only;
  # not combinable with eff_table_basename.
  angle_unit_case = rec.get("detparams_angle_unit")
  if angle_unit_case is not None:
    if str(angle_unit_case).lower() != "degree":
      log.error('detparams_angle_unit supports only "degree": %s', angle_unit_case)
      sys.exit(1)
    if eff_table_basename:
      log.error("detparams_angle_unit cannot be combined with eff_table_basename")
      sys.exit(1)
    if not cfg.get("skip_detparams", False):
      log.error("detparams_angle_unit requires skip_detparams=true")
      sys.exit(1)
  if eff_table_basename and bool(cfg.get("skip_detparams", False)):
    tag = eff_table_basename.removeprefix("eff_table_").removesuffix(".txt")
    case_label = f"{spec_label}__{tag}"
    _make_eff_case_detparams(
      work_dir, cfg, case_label, eff_table_basename, dry_run, force,
    )
    det_files_prefix = f"../detparams/{case_label}/"
  elif angle_unit_case:
    case_label = f"{spec_label}__degree"
    _make_angle_case_detparams(work_dir, cfg, case_label, dry_run, force)
    det_files_prefix = f"../detparams/{case_label}/"
  else:
    det_files_prefix = f"../detparams/{spec_label}/"
  det_files_block = _build_det_files_block(
    n_detectors, det_files_prefix, cfg.get("_det_filenames_override"),
  )
  det_list_name = _det_list_name(spec_label)

  base_density_int = int(rec.get("initial_uniform_density", 2000))
  g2bg_dens_cnt_delta = float(rec.get("g2bg_dens_cnt_delta", 500.0))
  g2bg_dens_cnt_vmin = base_density_int - g2bg_dens_cnt_delta
  g2bg_dens_cnt_vmax = base_density_int + g2bg_dens_cnt_delta

  # Default to station center; allow swp001 override for cases where the
  # reconstruction grid is aligned to a different origin (e.g. blas_swp003
  # uses (-2.5, -2.5) regardless of summit position).
  center_x = rec.get("x_merge_center", cfg["center_x"])
  center_y = rec.get("y_merge_center", cfg["center_y"])

  z_pitch = rec["z_pitch"]

  # 1) run_prg.sh
  content = TPL_REC001_RUN.substitute(
    site_name=site_name,
    swp_work_subdir=swp_work_subdir,
    epsg=epsg,
    center_x=cfg["center_x"],
    center_y=cfg["center_y"],
    plot_zmin=cfg["plot_zmin"],
    plot_zmax=cfg["plot_zmax"],
    dem_gen_block=_build_dem_gen_block(cfg),
  )
  _write_file(work_dir / swp_work_subdir / "run_prg.sh", content, dry_run, force)

  # 2) prm_muonith.json5
  dem_path = f"../dem/{Path(cfg['dem_file']).name}"

  # reconst_voxels block (region mask after merge)
  rv = cfg.get("reconst_voxels", {})
  tf_aabb = rv.get("tf_aabb", rec.get("tf_aabb", False))
  tf_cylinder = rv.get("tf_cylinder", False)

  # Center coordinates: fall back to merge center (already in meters)
  x_aabb_cnt = _fmt_num(rv.get("x_aabb_cnt", center_x))
  y_aabb_cnt = _fmt_num(rv.get("y_aabb_cnt", center_y))
  x_cyl_cnt  = _fmt_num(rv.get("x_cyl_cnt",  center_x))
  y_cyl_cnt  = _fmt_num(rv.get("y_cyl_cnt",  center_y))

  reconst_voxels_block = (
    ', "reconst_voxels": {\n'
    f'        "tf_aabb": {"true" if tf_aabb else "false"}\n'
    f'      , "x_aabb_cnt": {x_aabb_cnt}\n'
    f'      , "y_aabb_cnt": {y_aabb_cnt}\n'
    f'      , "x_aabb_meters": {_fmt_num(rv.get("x_aabb_meters", 0.0))}\n'
    f'      , "y_aabb_meters": {_fmt_num(rv.get("y_aabb_meters", 0.0))}\n'
    f'      , "aabb_zmin_mode": "{rv.get("aabb_zmin_mode", "g3vox_zmin")}"  // "g3vox_zmin" or "manual" (manual requires "aabb_zmin_value")\n'
    f'      , "aabb_zmax": {_fmt_num(float(rv.get("aabb_zmax", 9999)))}\n'
    f'      , "tf_cylinder": {"true" if tf_cylinder else "false"}\n'
    f'      , "x_cyl_cnt": {x_cyl_cnt}\n'
    f'      , "y_cyl_cnt": {y_cyl_cnt}\n'
    f'      , "cylinder_radius_x_meters": {_fmt_num(rv.get("cylinder_radius_x_meters", 0.0))}\n'
    f'      , "cylinder_radius_y_meters": {_fmt_num(rv.get("cylinder_radius_y_meters", 0.0))}\n'
    '    }\n'
  )

  # Migration warning: detect legacy AABB fields in swp001
  _legacy_aabb_keys = {"tf_aabb", "x_aabb_cells", "y_aabb_cells", "aabb_zmin_mode", "aabb_zmax"}
  found_legacy = _legacy_aabb_keys & set(rec.keys())
  if found_legacy:
    log.warning(
      "[%s] swp001 contains legacy AABB fields: %s. "
      "Move them to 'reconst_voxels' section.",
      site_name, sorted(found_legacy),
    )

  # checkerboard_3d_params block (one or more entries from cfg)
  cb_list = cfg.get("checkerboard_3d_params") or [{
    "tf_exec":              True,
    "name":                 "cb_3d_500",
    "delta_density_offset": 0.0,
    "delta_density":        500,
    "xlen_cells":           6,
    "ylen_cells":           6,
    "zlen_cells":           4,
    "xcnt":                 center_x,
    "ycnt":                 center_y,
    "zcnt":                 rec["z_merge_center"],
    "xlen_interval_mult":   40,
    "ylen_interval_mult":   40,
    "zlen_interval_mult":   40,
    "tf_snap_to_grid":      True,
  }]
  cb_entries = []
  for i, cb in enumerate(cb_list):
    trailing = "," if i < len(cb_list) - 1 else ""
    cb_entries.append(
      "      {\n"
      f'        "tf_exec": {"true" if cb.get("tf_exec", True) else "false"} // Execution flag\n'
      f'      , "name": "{cb.get("name", "cb_3d_500")}" // Checkerboard name\n'
      f'      , "delta_density_offset": {_fmt_num(float(cb.get("delta_density_offset", 0.0)))} // delta_density_offset value\n'
      f'      , "delta_density": {cb.get("delta_density", 500)} // delta_density value\n'
      f'      , "xlen_cells": {cb.get("xlen_cells", 6)}\n'
      f'      , "ylen_cells": {cb.get("ylen_cells", 6)}\n'
      f'      , "zlen_cells": {cb.get("zlen_cells", 4)}\n'
      f'      , "xcnt": {_fmt_num(float(cb.get("xcnt", center_x)))} // Center x coordinate\n'
      f'      , "ycnt": {_fmt_num(float(cb.get("ycnt", center_y)))} // Center y coordinate\n'
      f'      , "zcnt": {_fmt_num(float(cb.get("zcnt", rec["z_merge_center"])))} // Center z coordinate\n'
      f'      , "xlen_interval_mult": {cb.get("xlen_interval_mult", 40)}\n'
      f'      , "ylen_interval_mult": {cb.get("ylen_interval_mult", 40)}\n'
      f'      , "zlen_interval_mult": {cb.get("zlen_interval_mult", 40)}\n'
      f'      , "tf_snap_to_grid": {"true" if cb.get("tf_snap_to_grid", True) else "false"}\n'
      "      }" + trailing
    )
  checkerboard_3d_entries = "\n".join(cb_entries)

  # cuboid_params block (opt-in per site; empty by default -> no anomaly).
  # Fields mirror Grid3dVoxel::CuboidParameters (cnt / len / theta / rotation_type).
  cub_list = cfg.get("cuboid_params") or []
  cub_entries = []
  for i, cub in enumerate(cub_list):
    trailing = "," if i < len(cub_list) - 1 else ""
    cub_entries.append(
      "      {\n"
      f'        "tf_exec": {"true" if cub.get("tf_exec", True) else "false"} // Execution flag\n'
      f'      , "name": "{cub.get("name", "cuboid_01")}" // Cuboid name\n'
      f'      , "delta_density": {_fmt_num(float(cub.get("delta_density", 0.0)))} // delta_density value (kg/m3)\n'
      f'      , "xcnt": {_fmt_num(float(cub.get("xcnt", center_x)))} // Center x coordinate\n'
      f'      , "ycnt": {_fmt_num(float(cub.get("ycnt", center_y)))} // Center y coordinate\n'
      f'      , "zcnt": {_fmt_num(float(cub.get("zcnt", rec["z_merge_center"])))} // Center z coordinate\n'
      f'      , "xlen": {_fmt_num(float(cub.get("xlen", 100.0)))} // Full side length x (m)\n'
      f'      , "ylen": {_fmt_num(float(cub.get("ylen", 100.0)))} // Full side length y (m)\n'
      f'      , "zlen": {_fmt_num(float(cub.get("zlen", 100.0)))} // Full side length z (m)\n'
      f'      , "theta_x_deg": {_fmt_num(float(cub.get("theta_x_deg", 0.0)))} // Rotation around x (deg)\n'
      f'      , "theta_y_deg": {_fmt_num(float(cub.get("theta_y_deg", 0.0)))} // Rotation around y (deg)\n'
      f'      , "theta_z_deg": {_fmt_num(float(cub.get("theta_z_deg", 0.0)))} // Rotation around z (deg)\n'
      f'      , "rotation_type": "{cub.get("rotation_type", "LOCAL")}" // LOCAL or GLOBAL\n'
      "      }" + trailing
    )
  cuboid_params_entries = "\n".join(cub_entries)

  # end_stage line: comment if absent, active line otherwise.
  if "end_stage" in rec:
    end_stage_line = f', "end_stage": {int(rec["end_stage"])}'
  else:
    end_stage_line = '// , "end_stage": 8'

  # module8_mode: default "all" preserves prior template behavior.
  module8_mode = str(rec.get("module8_mode", "all"))

  # RUN_INVERSION.tf_prior_error: default false keeps generated swp001 runs
  # centered unless a site explicitly opts into lower/upper prior-error passes.
  tf_prior_error = "true" if rec.get("tf_prior_error", False) else "false"

  # NAGAINV_PARAMETERS.tf_signal_poisson: default false preserves prior template behavior.
  nagainv_tf_signal_poisson = (
    "true" if rec.get("nagainv_tf_signal_poisson", False) else "false"
  )
  noise_tf_exec = "true" if rec.get("noise_tf_exec", False) else "false"
  # PROJ_DENS_EVAL_GROUPED.tf_exec: default true preserves prior behavior for
  # all sites; tutorial opts out so the first run matches the paper's <1 min.
  proj_dens_tf_exec = "true" if rec.get("proj_dens_tf_exec", True) else "false"
  # DETECTOR_PARAMETER_LISTS.tf_out_txty_ascii / tf_out_g2bg_ascii: default true
  # preserves prior template behavior; a site opts out explicitly (e.g. omuro).
  tf_out_txty_ascii = "true" if rec.get("tf_out_txty_ascii", True) else "false"
  tf_out_g2bg_ascii = "true" if rec.get("tf_out_g2bg_ascii", True) else "false"
  # PROJ_DENS_EVAL_GROUPED.tf_signal_poisson: default true preserves prior
  # template behavior; a site opts out explicitly (e.g. omuro).
  proj_dens_tf_signal_poisson = (
    "true" if rec.get("proj_dens_tf_signal_poisson", True) else "false"
  )
  # NAGAINV_PARAMETERS.tf_eff_cn_diag: default false preserves prior template behavior.
  nagainv_tf_eff_cn_diag = (
    "true" if rec.get("nagainv_tf_eff_cn_diag", False) else "false"
  )
  nagainv_tf_eff_cn_diag_independent = (
    "true" if rec.get("nagainv_tf_eff_cn_diag_independent", False) else "false"
  )
  nagainv_tf_calc_chi2ndf = "true" if rec.get("nagainv_tf_calc_chi2ndf", False) else "false"  # default false: p_eff/chi2 self-eval, opt-in per site
  # NOISE_PARAMETERS: noise split into floor (deterministic) and poisson (fluctuated)
  # for both flux-proportional and SOT-proportional components. Default 0.0 = no noise.
  flux_proport_ratio_floor = rec.get("flux_proport_ratio_floor", "0.0")
  flux_proport_ratio_poisson = rec.get("flux_proport_ratio_poisson", "0.0")
  sot_proport_noise_ratio_floor = rec.get("SOT_proport_noise_ratio_floor", "0.0")
  sot_proport_noise_ratio_poisson = rec.get("SOT_proport_noise_ratio_poisson", "0.0")

  content = TPL_PRM_MUONITH.substitute(
    work_dir_name=work_dir.name,
    site_name=site_name,
    site_name_cap=site_name_cap,
    flux_groom=cfg.get("flux_groom", "daemon_groom"),
    det_list_name=det_list_name,
    det_files_block=det_files_block,
    base_density_int=base_density_int,
    dem_path=dem_path,
    BL_max=rec["BL_max"],
    PL_max=rec["PL_max"],
    g3vox_zmin=rec["g3vox_zmin"],
    g3vox_zmax=rec["g3vox_zmax"],
    z_pitch=z_pitch,
    center_x=center_x,
    center_y=center_y,
    z_merge_center=_fmt_num(rec["z_merge_center"]),
    merge_factor=rec["merge_factor"],
    g2pil_zmin=_fmt_num(float(rec.get("g2pil_zmin", 0.0))),
    n_hit_det_min=rec.get("n_hit_det_min", 3),
    corr_length=rec["corr_length"],
    sigma_rho=rec["sigma_rho"],
    zcross_min=_fmt_num(float(rec.get("zcross_min", rec["g3vox_zmin"]))),
    zcross_max=_fmt_num(float(rec.get("zcross_max", rec["g3vox_zmax"]))),
    zcross_zstep=_fmt_num(float(rec.get("zcross_zstep", rec["z_pitch"] * rec["merge_factor"]))),
    reconst_voxels_block=reconst_voxels_block,
    checkerboard_3d_entries=checkerboard_3d_entries,
    cuboid_params_entries=cuboid_params_entries,
    end_stage_line=end_stage_line,
    module8_mode=module8_mode,
    tf_prior_error=tf_prior_error,
    nagainv_tf_signal_poisson=nagainv_tf_signal_poisson,
    nagainv_tf_eff_cn_diag=nagainv_tf_eff_cn_diag,
    nagainv_tf_eff_cn_diag_independent=nagainv_tf_eff_cn_diag_independent,
    nagainv_tf_calc_chi2ndf=nagainv_tf_calc_chi2ndf,
    noise_tf_exec=noise_tf_exec,
    proj_dens_tf_exec=proj_dens_tf_exec,
    proj_dens_tf_signal_poisson=proj_dens_tf_signal_poisson,
    tf_out_txty_ascii=tf_out_txty_ascii,
    tf_out_g2bg_ascii=tf_out_g2bg_ascii,
    flux_proport_ratio_floor=flux_proport_ratio_floor,
    flux_proport_ratio_poisson=flux_proport_ratio_poisson,
    SOT_proport_noise_ratio_floor=sot_proport_noise_ratio_floor,
    SOT_proport_noise_ratio_poisson=sot_proport_noise_ratio_poisson,
    shell_density_lines=_build_shell_density_lines(rec, base_density_int),
  )
  _write_file(work_dir / swp_work_subdir / "prm_muonith.json5", content, dry_run, force)

  # 3) auto_plot.json5
  # Angle-axis frame follows the runcards this run actually references
  # (det_files_prefix target), so unit variants (e.g. degree) get a matching
  # hist2d frame.  Fall back to detector_spec when no runcard is readable.
  det_dir = work_dir / "detparams" / det_files_prefix.rstrip("/").rsplit("/", 1)[-1]
  angle_bins = _read_runcard_angle_bins(det_dir)
  if angle_bins is None:
    log.info(
      "auto_plot angle bins: no readable runcard in %s; falling back to detector_spec",
      det_dir,
    )
    angle_bins = {key: spec[key] for key in ("nbinx", "txmin", "txmax", "nbiny", "tymin", "tymax")}
  content = TPL_AUTO_PLOT.substitute(
    site_name_cap=site_name_cap,
    nbinx=angle_bins["nbinx"],
    txmin=angle_bins["txmin"],
    txmax=angle_bins["txmax"],
    nbiny=angle_bins["nbiny"],
    tymin=angle_bins["tymin"],
    tymax=angle_bins["tymax"],
    PL_hist_max=rec["PL_hist_max"],
    # hist2d binz for the per-element signal figures ("<nbin> <min> <max>",
    # log-scaled).  Sites whose muon counts span a different decade range
    # override it in their swp001 block; the default preserves the historical
    # value so existing stations regenerate byte-identical auto_plot.json5.
    signal_binz=rec.get("signal_binz", "18 1.0E-3 1.0E+3"),
    error_vmax=int(rec.get("error_vmax", 1000)),
    error_cbar_step=int(rec.get("error_cbar_step", 100)),
    g2bg_dens_cnt_vmin=_fmt_num_one_decimal(g2bg_dens_cnt_vmin),
    g2bg_dens_cnt_vmax=_fmt_num_one_decimal(g2bg_dens_cnt_vmax),
    g3vox_plot_zmin=rec.get("plot_zmin", rec["g3vox_zmin"] - 110),
    g3vox_plot_zmax=rec.get("plot_zmax", rec["g3vox_zmax"] + 10),
    g3vox_plot_zstep=rec.get("plot_zstep", 20),
    g3vox_plot_xmin=rec.get("plot_xmin", center_x - 1500),
    g3vox_plot_xmax=rec.get("plot_xmax", center_x + 1500),
    g3vox_plot_xstep=rec.get("plot_xstep", 20),
    g3vox_plot_ymin=rec.get("plot_ymin", center_y - 2000),
    g3vox_plot_ymax=rec.get("plot_ymax", center_y + 2000),
    g3vox_plot_ystep=rec.get("plot_ystep", 20),
  )
  _write_file(work_dir / swp_work_subdir / "auto_plot.json5", content, dry_run, force)


# ---------------------------------------------------------------------------
# Pipeline execution (--run)
# ---------------------------------------------------------------------------

def run_pipeline(work_dir: Path, cfg: dict) -> None:
  """Execute the KML/GeoJSON -> CSV -> det_XX.json5 pipeline."""
  detparams_dir = work_dir / "detparams"
  spec_label = cfg["detector_spec"]["label"]
  spec_dir = detparams_dir / spec_label

  # Step 1: Run mk_coordinate_list.sh
  mk_script = detparams_dir / "mk_coordinate_list.sh"
  log.info("Running: %s", mk_script)
  result = subprocess.run(
    ["bash", str(mk_script)],
    cwd=str(detparams_dir),
    capture_output=True, text=True,
  )
  if result.returncode != 0:
    log.error("mk_coordinate_list.sh failed (exit %d):\n%s",
              result.returncode, result.stderr)
    sys.exit(1)
  if result.stdout:
    log.info("mk_coordinate_list.sh stdout:\n%s", result.stdout)

  # Step 2: Run mkjson5.sh
  mkjson5_script = spec_dir / "mkjson5.sh"
  log.info("Running: %s", mkjson5_script)
  result = subprocess.run(
    ["bash", str(mkjson5_script)],
    cwd=str(spec_dir),
    capture_output=True, text=True,
  )
  if result.returncode != 0:
    log.error("mkjson5.sh failed (exit %d):\n%s",
              result.returncode, result.stderr)
    sys.exit(1)
  if result.stdout:
    log.info("mkjson5.sh stdout:\n%s", result.stdout)

  log.info("Pipeline complete. Detector JSON5 files created in %s", spec_dir)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
  parser = argparse.ArgumentParser(
    description="Generate a volcano analysis work site from site.json5",
  )
  parser.add_argument("site_json5", help="Path to site.json5")
  parser.add_argument("--run", action="store_true",
                      help="Execute KML->CSV->det_XX.json5 pipeline")
  parser.add_argument("--dry-run", action="store_true",
                      help="Show what would be created without writing")
  parser.add_argument("--force", action="store_true",
                      help="Overwrite existing files")
  parser.add_argument("--clean", action="store_true",
                      help="Remove tmp/, logs/, figs/ under swp001 before generation")
  parser.add_argument("--verbose", action="store_true",
                      help="Print detailed progress")
  args = parser.parse_args()

  # Logging setup
  level = logging.DEBUG if args.verbose else logging.INFO
  logging.basicConfig(
    format="%(levelname)s: %(message)s",
    level=level,
  )

  # Load and validate
  cfg = load_site_config(args.site_json5)
  # Backward compatibility: accept "site_name" as alias for "station_name"
  if "station_name" not in cfg and "site_name" in cfg:
    cfg["station_name"] = cfg.pop("site_name")
  validate_config(cfg)
  convert_center_if_degree(cfg)

  # Determine work directory
  site_json5_path = Path(args.site_json5).resolve()
  work_dir = resolve_work_dir(site_json5_path)
  repo_root = _find_repo_root(site_json5_path)

  log.info("Site: %s", cfg["station_name"])
  log.info("Work dir: %s", work_dir)

  # Create directory skeleton
  if not args.dry_run:
    create_directories(work_dir, cfg)

  # Copy input files from data/ if available
  data_key = cfg.get("data_source", cfg["station_name"])
  data_dir = resolve_data_dir(data_key, repo_root)
  if data_dir is not None:
    log.info("Data dir: %s", data_dir)
    link_input_files(data_dir, work_dir, cfg, args.dry_run, args.force)
  else:
    log.debug("No data dir found for '%s' (data/%s/)", data_key, data_key)

  skip_dp = bool(cfg.get("skip_detparams", False))

  # Detect number of detectors:
  #   - skip_detparams=True: count existing runcard_det_*.json5 files
  #     and auto-populate cfg["detector_spec"] / cfg["epsg"] from the first
  #     runcard so downstream generators (generate_swp001) do not need to
  #     special-case skip_detparams.
  #   - else: parse the KML/GeoJSON
  if skip_dp:
    import json5  # already validated by load_site_config()
    # In skip_detparams flow, prefer data/<site>/detparams/<label>/.
    # Some legacy stations keep the canonical detector files only under
    # work/<site>/detparams/<label>/; support that as a read-only source too.
    src_detparams_roots = []
    if data_dir is not None:
      src_detparams_roots.append(data_dir / "detparams")
    src_detparams_roots.append(work_dir / "detparams")
    n_detectors = 0
    spec_label_dir = None
    first_runcard = None
    src_detparams = None
    # Prefer the label specified in cfg.detector_spec.label if present;
    # otherwise pick the first directory containing runcard_det_*.json5.
    preferred_label = cfg.get("detector_spec", {}).get("label", "") if isinstance(cfg.get("detector_spec"), dict) else ""
    candidates = []
    for root in src_detparams_roots:
      if root.is_dir():
        if preferred_label:
          pref_dir = root / preferred_label
          if pref_dir.is_dir():
            candidates.append(pref_dir)
        candidates.extend(
          d for d in sorted(root.iterdir())
          if d.is_dir() and d.name != "eff_table" and d not in candidates
        )
    for label_dir in candidates:
      runcards = sorted(label_dir.glob("runcard_det_*.json5"))
      if not runcards:
        runcards = sorted(label_dir.glob("det_*.json5"))
      if runcards:
        n_detectors = len(runcards)
        spec_label_dir = label_dir.name
        first_runcard = runcards[0]
        src_detparams = label_dir.parent
        break
    if n_detectors == 0 or first_runcard is None:
      log.error(
        "skip_detparams=true but no runcard_det_*.json5 or det_*.json5 found "
        "under data/<site>/detparams/<label>/ or work/<site>/detparams/<label>/",
      )
      sys.exit(1)
    log.info(
      "skip_detparams=true: reusing %d detector json5 files from %s/%s/",
      n_detectors, src_detparams, spec_label_dir,
    )
    # Optional: rename references to legacy det_NN.json5 in the generated prm
    # (per-file symlinks are created later by link_detparams_for_skip()).
    renumber = bool(cfg.get("detparams_renumber_to_det_nn", False))
    runcard_names = [p.name for p in runcards]
    cfg["_det_runcard_names"] = runcard_names  # original names, for symlink targets
    cfg["_det_source_root"] = src_detparams
    if renumber:
      cfg["_det_filenames_override"] = [f"det_{i:02d}.json5" for i in range(n_detectors)]
      log.info(
        "skip_detparams=true: detparams_renumber_to_det_nn=true -> "
        "prm will reference det_00..det_%02d.json5 (per-file symlinks)",
        n_detectors - 1,
      )
    else:
      cfg["_det_filenames_override"] = runcard_names

    # Auto-populate detector_spec from the first runcard, then mirror
    # detector_overrides because those values are written to generated runcards.
    with open(first_runcard, encoding="utf-8") as fh:
      rc = json5.load(fh)
    dp = rc.get("DETECTOR_PARAMETERS", rc)
    auto_spec = {
      "label":       spec_label_dir,
      "length_hori": dp.get("length_hori", 0.333),
      "length_vert": dp.get("length_vert", 0.333),
      "length_dept": dp.get("length_dept", 0.03),
      "n_unit":      dp.get("n_unit", 2),
      "days":        dp.get("days", 120),
      "nbinx":       dp.get("nbinx", 320),
      "txmin":       dp.get("txmin", -1.6),
      "txmax":       dp.get("txmax", 1.6),
      "nbiny":       dp.get("nbiny", 160),
      "tymin":       dp.get("tymin", 0.0),
      "tymax":       dp.get("tymax", 1.6),
    }
    user_spec = cfg.get("detector_spec") or {}
    detector_spec_overrides = cfg.get("detector_overrides") or {}
    merged = {**auto_spec, **user_spec}  # user values override auto values
    for key in [
      "length_hori", "length_vert", "length_dept", "n_unit", "days",
      "nbinx", "txmin", "txmax", "nbiny", "tymin", "tymax",
    ]:
      if key in detector_spec_overrides:
        merged[key] = detector_spec_overrides[key]
    cfg["detector_spec"] = merged
    log.info(
      "skip_detparams=true: auto-populated detector_spec from %s (label=%s)",
      first_runcard.name, merged["label"],
    )
    if "epsg" not in cfg:
      cfg["epsg"] = 0  # custom CRS marker; only used as a string token in templates
  else:
    det_input_name = Path(cfg["detector_input"]).name
    kml_path = None
    if data_dir is not None:
      kml_path = data_dir / det_input_name
    if kml_path is None or not kml_path.exists():
      kml_path = work_dir / cfg["detector_input"]
    if kml_path.exists():
      n_detectors = count_detectors_in_kml(kml_path)
      log.info("Detected %d detectors from %s", n_detectors, kml_path.name)
    else:
      log.error("KML not found at %s", kml_path)
      sys.exit(1)

  # Clean swp output directories if requested
  if args.clean:
    for swp_subdir in _swp_work_subdirs(cfg):
      swp001_dir = work_dir / swp_subdir
      for subdir in ("tmp", "logs", "figs"):
        target = swp001_dir / subdir
        if target.is_dir():
          if args.dry_run:
            log.info("DRY-RUN: would remove %s", target)
          else:
            import shutil
            shutil.rmtree(target)
            log.info("Removed: %s", target)

  # Generate files
  dry_run = args.dry_run
  force = args.force

  if skip_dp:
    log.info(
      "skip_detparams=true: skipping generate_detparams() — reusing existing detparams/"
    )
    repo_root = Path(__file__).resolve().parent.parent
    renumber_list = (
      cfg.get("_det_runcard_names")
      if cfg.get("detparams_renumber_to_det_nn", False)
      else None
    )
    overrides = cfg.get("detector_overrides") or None
    link_detparams_for_skip(
      repo_root, cfg["_det_source_root"], work_dir, cfg["detector_spec"]["label"],
      dry_run, force, renumber_runcards=renumber_list,
      detector_overrides=overrides,
    )
  else:
    generate_detparams(work_dir, cfg, dry_run, force)
  if not dry_run:
    _validate_detparams_eff_table_grid(work_dir, cfg, cfg["detector_spec"]["label"])

  if "depth_sweep" in cfg:
    generate_depth001(work_dir, cfg, n_detectors, dry_run, force)
  else:
    log.info("depth_sweep section absent: skipping generate_depth001()")

  for run_cfg in _iter_swp001_run_cfgs(cfg):
    generate_swp001(work_dir, run_cfg, n_detectors, dry_run, force)

  # Optionally run the pipeline (KML -> CSV -> det.json5)
  if args.run:
    if skip_dp:
      log.info(
        "skip_detparams=true: skipping pipeline (KML/CSV/det conversion bypassed)"
      )
    elif dry_run:
      log.info("DRY-RUN: would execute pipeline (mk_coordinate_list.sh + mkjson5.sh)")
    else:
      run_pipeline(work_dir, cfg)

  log.info("Done.")


if __name__ == "__main__":
  main()
