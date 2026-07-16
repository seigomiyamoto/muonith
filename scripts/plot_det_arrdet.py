#!/usr/bin/env python3
"""Plot per-element quantities straight from det/ arrdet_*.bin (direct read).

This reads a serialized DetectorPanelArray (checkpoint det/arrdet_g3vox_input.bin
or det/arrdet_g2pil_naive.bin, ~160 MB each) field-by-field, mirroring the C++
write order, then reproduces the standard txty figures by calling the functions
of scripts/hist2d.py directly on the in-memory (tx, ty, value) arrays, with
the hist2d parameters taken from the same auto_plot.json5 that the normal
auto_plot.py flow uses. No intermediate txty text file is written: the whole
point of this tool is to go binary -> figure directly. Output naming and style
are identical to the established figures (one PNG per detector per field, e.g.
arrdet_g3vox_input_txty_PL_det00.png).

Unlike scripts/demo_recon_*.py (which read the named-column mat/ element table)
and scripts/plot_g2bg.py (which reads the grouped ASCII g2bg text), this tool
reads the raw binary directly, so it also works on checkpoints written before
the mat/ element-table export existed. For per-element scalars the values are
identical to mat/mat_element_table_det{NN}.bin; pass --check-mat to verify that
row-for-row.

Binary format (little-endian, LP64), mirrored from the C++ save methods:
  - file header          io_binary::write_architecture_info  ns_io_binary.cpp:605
  - DetectorPanelArray   src/cls_DetectorPanelArray.cpp:2553
  - DetectorPanel        src/cls_DetectorPanel.cpp:2228
  - Grid2dBinGroup       src/cls_Grid2dBinGroup.cpp:1993 (+ Grid2d :803, Grid1d :405)
  - Grid2dBinGroup::Parameters  src/cls_Grid2dBinGroupParameters.cpp:201
  - Ray3d                src/cls_Ray.cpp:493
  - DetectorElement      src/cls_DetectorElement.cpp:427 (192-byte fixed record)
The vec<bool> count is a uint32 (ns_io_binary.cpp:97); every other container /
string count is a size_t (8 bytes).

Usage:
  python3 scripts/plot_det_arrdet.py <arrdet.bin | checkpoint_dir | det_dir> \
      [--field PL,dens,signal] [--det 0,1,2] [--out DIR] [--jobs N] \
      [--plot-config auto_plot.json5] [--check-mat <mat_dir>] [--list-fields]

Also invoked automatically by auto_plot.py (the "det" task family) when the
auto_plot config sets det_exec: true.
"""

import argparse
import os
import struct
import sys
from pathlib import Path

import numpy as np

# txty field name (same as the C++ out_txty* suffix) -> per-element value.
# The extractors mirror cls_DetectorPanelArray.hpp out_txtyPL/DL/Signal/Noise/
# SignalPlusNoise/Dens (:819-858). Rendering parameters (binz, colormap, log,
# ...) are NOT fixed here: they come from the same auto_plot.json5 the normal
# auto_plot.py flow uses, so the figures match the established ones exactly.
TXTY_FIELDS = {
  "PL": lambda ele: np.asarray(ele["PL"], np.float64),
  "DL": lambda ele: np.asarray(ele["DL"], np.float64),
  "signal": lambda ele: np.asarray(ele["signal"], np.float64),
  # get_noise() = noise_det + noise_poi
  "noise": lambda ele: np.asarray(ele["noise_det"], np.float64)
                       + np.asarray(ele["noise_poi"], np.float64),
  "signal_plus_noise": lambda ele: np.asarray(ele["signal"], np.float64)
                                   + np.asarray(ele["noise_det"], np.float64)
                                   + np.asarray(ele["noise_poi"], np.float64),
  # out_txtyDens uses get_proj_density()
  "dens": lambda ele: np.asarray(ele["proj_density"], np.float64),
}

# One DetectorElement fixed record = 216 bytes, packed with no padding in the
# exact write order of DetectorElement::save (src/cls_DetectorElement.cpp:427).
# eff_low/eff_cnt/eff_upp were appended in PIPELINE_VERSION 7.
# The trailing variable-length vec_tf_in_PL is read and skipped separately.
ELE_DTYPE = np.dtype({
  "names": [
    "unique_index", "detid", "id_in_this",
    "pos_x", "pos_y", "pos_z", "dir_x", "dir_y", "dir_z",
    "txmin", "txmax", "tymin", "tymax", "angle_unit",
    "effective_area_m2", "solid_angle", "exposure_time_sec",
    "PL", "DL", "penetrating_muon_flux", "signal", "noise_det", "noise_poi",
    "proj_density", "proj_density_lower", "proj_density_upper",
    "eff_low", "eff_cnt", "eff_upp",
  ],
  "formats": [
    "<i4", "<i4", "<i4",
    "<f8", "<f8", "<f8", "<f8", "<f8", "<f8",
    "<f8", "<f8", "<f8", "<f8", "<i4",
    "<f8", "<f8", "<f8",
    "<f8", "<f8", "<f8", "<f8", "<f8", "<f8",
    "<f8", "<f8", "<f8",
    "<f8", "<f8", "<f8",
  ],
  "offsets": [0, 4, 8, 12, 20, 28, 36, 44, 52, 60, 68, 76, 84, 92, 96, 104,
              112, 120, 128, 136, 144, 152, 160, 168, 176, 184, 192, 200, 208],
  "itemsize": 216,
})
ELE_FIXED_BYTES = 216


class Reader:
  """Sequential little-endian reader mirroring src/ns_io_binary.cpp helpers."""

  def __init__(self, f):
    self.f = f

  def raw(self, n):
    b = self.f.read(n)
    if len(b) != n:
      raise EOFError("unexpected end of file: wanted %d bytes, got %d" % (n, len(b)))
    return b

  def i32(self):
    return struct.unpack("<i", self.raw(4))[0]

  def u16(self):
    return struct.unpack("<H", self.raw(2))[0]

  def u32(self):
    return struct.unpack("<I", self.raw(4))[0]

  def u64(self):
    return struct.unpack("<Q", self.raw(8))[0]

  def f64(self):
    return struct.unpack("<d", self.raw(8))[0]

  def boolean(self):
    return self.raw(1)[0] != 0

  def string(self):  # write_string: size_t length + raw bytes (no NUL)
    n = self.u64()
    return self.raw(n).decode("utf-8", "replace")

  def vec_f64(self):  # write_vec<double>: size_t count + count doubles
    n = self.u64()
    return np.frombuffer(self.raw(n * 8), dtype="<f8") if n else np.empty(0, "<f8")

  def vec_bool(self):  # write_vec_bool: uint32 count (NOT size_t) + count bytes
    n = self.u32()
    return np.frombuffer(self.raw(n), dtype=np.uint8) if n else np.empty(0, np.uint8)

  def vec_vec_f64(self):  # size_t outer count + inner write_vec<double>
    n = self.u64()
    return [self.vec_f64() for _ in range(n)]

  def vec_vec_bool(self):  # size_t outer count + inner write_vec_bool
    n = self.u64()
    return [self.vec_bool() for _ in range(n)]

  def multimap_int_int2(self):  # size_t count + count * (int key, int, int)
    n = self.u64()
    return np.empty((0, 3), "<i4") if n == 0 else \
        np.frombuffer(self.raw(n * 12), dtype="<i4").reshape(n, 3)

  def skip_vec_tp_bool_double(self):  # size_t count + count * (1B bool + 8B double)
    n = self.u64()
    self.f.seek(n * 9, 1)


def read_grid1d(r):  # src/cls_Grid1d.cpp:400
  r.string()               # name
  nbin = r.i32()           # nbin
  raw_min = r.f64()        # raw_min
  raw_max = r.f64()        # raw_max (meaning depends on PIPELINE_VERSION, see below)
  interval = r.f64()       # interval
  # The Grid1d layout is identical across versions, but the stored max changed:
  #   PIPELINE_VERSION <= 5: save stores max - interval (non-half-shift,
  #     cls_Grid1d.cpp before 31b8dda4), so max = raw_max + interval.
  #   PIPELINE_VERSION >= 6: save is a pass-through (canonical min/max as-is,
  #     31b8dda4), so max = raw_max.
  # Both are distinguished unambiguously from the Grid1d invariant
  # max - min == nbin * interval, so no external version input is needed.
  span = raw_max - raw_min
  tol = 0.01 * interval if interval > 0 else 1e-12
  if abs(span - nbin * interval) <= tol:
    axis_max = raw_max                      # v6 pass-through
  elif abs(span - (nbin - 1) * interval) <= tol:
    axis_max = raw_max + interval           # v5 (max - interval was stored)
  else:
    raise ValueError(
      "Grid1d axis is inconsistent: nbin=%d, min=%g, max=%g, interval=%g "
      "(neither the v5 nor the v6 Grid1d::save layout matches)"
      % (nbin, raw_min, raw_max, interval))
  # "min" then feeds the cell edges get_lower_value(i) = min + i*interval,
  # and "min"/"max" give the g2bg header.
  return {"nbin": nbin, "min": raw_min, "max": axis_max, "interval": interval}


def read_g2bg(r):  # src/cls_Grid2dBinGroup.cpp:1993 (+ Grid2d :803)
  r.string()                 # Grid2d name
  x_axis = read_grid1d(r)    # x_axis
  y_axis = read_grid1d(r)    # y_axis
  r.string()                 # BinGroup name
  # bimap_ maps igroup -> list of (ix, iy) cells; each serialized triple is
  # (igroup, ix, iy) (save_bimap -> write_multimap_int_int2, cls_Grid2dBinGroup.cpp:2082).
  bimap_triples = r.multimap_int_int2()
  r.i32()                    # detid_
  r.vec_bool()               # vec_is_avail_group (per-group; writer recomputes, so skip)
  for _ in range(5):         # vec_signal_group / noise_poi / noise_det / signal_poisson / noise_poi_poisson
    r.vec_f64()
  vec_vec_is_avail = r.vec_vec_bool()     # [iy][ix] availability
  vec_vec_signal = r.vec_vec_f64()        # [iy][ix] signal
  vec_vec_noise_poi = r.vec_vec_f64()     # [iy][ix] Poisson noise
  r.vec_vec_f64()                         # vec_vec_noise_det (not needed)
  r.boolean()                # done_grouping
  groups = {}
  groups["dens_lower"] = r.vec_f64()
  groups["dens_center"] = r.vec_f64()
  groups["dens_upper"] = r.vec_f64()
  groups["delta_nmuon_lower"] = r.vec_f64()
  groups["delta_nmuon_center"] = r.vec_f64()
  groups["delta_nmuon_upper"] = r.vec_f64()
  groups["volume"] = r.vec_f64()
  # group efficiencies, appended in PIPELINE_VERSION 7
  groups["eff_low"] = r.vec_f64()
  groups["eff_cnt"] = r.vec_f64()
  groups["eff_upp"] = r.vec_f64()
  # Group igroup -> list of (ix, iy), the membership used to rebuild the g2bg rows.
  bimap = {}
  for igroup, ix, iy in bimap_triples:
    bimap.setdefault(int(igroup), []).append((int(ix), int(iy)))
  return {
    "nbinx": x_axis["nbin"], "nbiny": y_axis["nbin"],
    "x_axis": x_axis, "y_axis": y_axis, "groups": groups, "bimap": bimap,
    "vec_vec_is_avail": vec_vec_is_avail, "vec_vec_signal": vec_vec_signal,
    "vec_vec_noise_poi": vec_vec_noise_poi,
  }


def read_prm_bingrp(r):  # src/cls_Grid2dBinGroupParameters.cpp:201 — advance only
  r.string()                 # name
  r.f64(); r.f64()           # signal_init, noise_init
  r.boolean()                # is_avail_init
  for _ in range(4):         # PL_thres, DL_thres, signal_under_thres, noise_under_thres
    r.f64()
  r.boolean(); r.boolean(); r.boolean()  # is_avail_under_thres, tf_run_1st_grouping, tf_run_auto_grouping
  r.i32(); r.i32(); r.i32()  # igroup_start, nx_div_init, ny_div_init
  r.f64()                    # signal_noise_group_trig
  r.i32(); r.i32()           # ixlen_min, iylen_min
  r.boolean()                # tf_prefer_split_x
  r.i32(); r.i32()           # nloop_limit, n_detector_grouping_manual
  r.vec_bool()               # vec_tf_read_bin_group_list
  n_path = r.u64()           # vec_file_path_bin_group_list count
  if n_path != 0:
    raise ValueError(
      "vec_file_path_bin_group_list has %d entries; std::vector<fs::path> is "
      "serialized as raw non-portable bytes and cannot be parsed here." % n_path)


def read_panel(r):  # src/cls_DetectorPanel.cpp:2228
  g2bg = read_g2bg(r)                    # g2bg_
  nbinx, nbiny = g2bg["nbinx"], g2bg["nbiny"]
  detid = r.i32()                       # detid_
  name = r.string()                     # name
  r.raw(48)                             # ray3d: v3_pos + v3_dir (6 doubles)
  r.raw(24)                             # v3_det_length (3 doubles)
  r.f64()                               # n_unit
  r.f64()                               # days
  r.i32()                               # angle_unit
  n_element = r.i32()                   # n_element
  # No count precedes the elements; the shape is nbiny x nbinx (row-major iy,ix).
  n = nbiny * nbinx
  buf = bytearray()
  for _ in range(n):
    buf += r.raw(ELE_FIXED_BYTES)       # fixed 192-byte record
    r.skip_vec_tp_bool_double()         # trailing vec_tf_in_PL
  ele = np.frombuffer(bytes(buf), dtype=ELE_DTYPE)
  read_prm_bingrp(r)                    # prm_bingrp (advance to next panel)
  return {
    "detid": detid, "name": name, "nbinx": nbinx, "nbiny": nbiny,
    "n_element": n_element, "ele": ele, "groups": g2bg["groups"], "g2bg": g2bg,
  }


def read_pipeline_version(arrdet_path):
  """Read PIPELINE_VERSION from the checkpoint's state_meta.bin, or None.

  A checkpoint bundle is laid out as <bundle>/state_meta.bin + <bundle>/det/
  (exemdl_pipeline.cpp save()), so state_meta.bin sits next to the det/ dir
  holding the arrdet binary. Standalone det/ dirs (tf_save_arrdet_* exports)
  have no state_meta.bin; return None for those.
  """
  meta = Path(arrdet_path).resolve().parent.parent / "state_meta.bin"
  if not meta.is_file():
    return None
  with open(meta, "rb") as f:
    r = Reader(f)
    r.string()               # architectureName
    r.boolean()              # isLittleEndian
    r.u16()                  # pointerSize
    return r.u32()           # PIPELINE_VERSION


MIN_PIPELINE_VERSION = 7  # v7 added eff to DetectorElement and Grid2dBinGroup::save


def read_arrdet(path):
  """Read arrdet_*.bin; return (array_name, list of detector-panel dicts)."""
  path = Path(path)
  version = read_pipeline_version(path)
  if version is not None and version < MIN_PIPELINE_VERSION:
    raise SystemExit(
      "ERROR: checkpoint %s has PIPELINE_VERSION %d, but this reader supports "
      "version %d or newer (older files lack the eff_low/eff_cnt/eff_upp "
      "fields added to DetectorElement::save and Grid2dBinGroup::save in v7 "
      "and cannot be parsed). Regenerate the checkpoint with the current code."
      % (path, version, MIN_PIPELINE_VERSION))
  with open(path, "rb") as f:
    r = Reader(f)
    # File header: write_architecture_info (name + isLittleEndian + pointerSize).
    r.string()
    little_endian = r.boolean()
    pointer_size = r.u16()
    if not little_endian or pointer_size != 8:
      raise ValueError("unsupported architecture in %s: little_endian=%s, "
                       "pointer_size=%d (this reader assumes LE + LP64)"
                       % (path.name, little_endian, pointer_size))
    # DetectorPanelArray::save
    arr_name = r.string()     # name (e.g. "arrdet_g3vox_input")
    r.i32()                   # n_all_element
    n_path = r.u64()          # vec_parameter_file_path count
    for _ in range(n_path):
      r.string()              # write_path == write_string
    n_det = r.u64()           # vec_panel count
    panels = [read_panel(r) for _ in range(n_det)]
    # dic_ / flags_ / set_uqiv follow but are not needed for per-element plots.
  return arr_name, panels


def resolve_arrdet_path(arg):
  """Accept an arrdet .bin, a checkpoint dir, or a det/ dir; return the .bin."""
  p = Path(arg)
  if p.is_file():
    return p
  candidates = ["det/arrdet_g3vox_input.bin", "det/arrdet_g2pil_naive.bin",
                "arrdet_g3vox_input.bin", "arrdet_g2pil_naive.bin"]
  for rel in candidates:
    cand = p / rel
    if cand.is_file():
      return cand
  raise SystemExit("ERROR: no arrdet_*.bin found at or under %s" % p)


def resolve_plot_config(arrdet_path, arg):
  """Pick the auto_plot json5: --plot-config, nearest auto_plot.json5, default.

  Walks up from the arrdet .bin (checkpoint tmp dir -> work dir) looking for
  auto_plot.json5, the same file the normal auto_plot.py flow would use there;
  falls back to scripts/auto_plot_default.json5.
  """
  if arg is not None:
    p = Path(arg)
    if not p.is_file():
      raise SystemExit("ERROR: --plot-config not found: %s" % p)
    return p
  for parent in Path(arrdet_path).resolve().parents:
    cand = parent / "auto_plot.json5"
    if cand.is_file():
      return cand
  return Path(__file__).resolve().parent / "auto_plot_default.json5"


def load_hist2d_config(config_path):
  """Load the whole auto_plot config plus the hist2d template/patterns.

  Returns (config, template, patterns). config is the full dict so the caller
  can honor the same switches as auto_plot.py: the top-level hist2d_exec /
  g2bg_exec gates and the exclude_patterns lists (global + per section).
  """
  sys.path.insert(0, str(Path(__file__).resolve().parent))
  from auto_plot import load_json_file
  config = load_json_file(str(config_path))
  hist2d = config.get("hist2d", {})
  return config, hist2d.get("template", {}), hist2d.get("patterns", {})


def render_with_hist2d(out_dir, arr_name, panel, field, template, patterns):
  """Render one detector/field via hist2d.py functions, no data file.

  Follows hist2d.main() histogram mode step for step (compute_histogram ->
  compute_norm -> setup_colormap -> render_plot) on the in-memory element
  arrays, so the PNG is identical to the auto_plot.py output for the same
  data. The pattern merge is the same as auto_plot.collect (template +
  first matching pattern).
  """
  import matplotlib
  matplotlib.use("Agg")
  from auto_plot import match_pattern
  import hist2d as h2

  # Same source name auto_plot would see; drives pattern match, title, output.
  src_name = "%s_txty_%s_det%02d.tmp" % (arr_name, field, panel["detid"])
  matched = match_pattern(src_name, patterns)
  merged = {**template, **matched} if matched else dict(template)
  use_log = merged.get("log", False)
  out_name = src_name.replace(".tmp", "") + ("_log" if use_log else "") + ".png"

  # Build the same argparse namespace the CLI would produce.
  ngrad, vmin, vmax = merged.get("binz", "20 0 2000").split()
  nbinx, xmin, xmax = merged.get("binx", "320 -1.6 1.6").split()
  nbiny, ymin, ymax = merged.get("biny", "160 0.0 1.6").split()
  argv = [
    src_name,
    "--nbinx", nbinx, "--xmin", xmin, "--xmax", xmax,
    "--nbiny", nbiny, "--ymin", ymin, "--ymax", ymax,
    "--ngrad", ngrad, "--vmin", vmin, "--vmax", vmax,
    "--dpi", str(merged.get("dpi", 300)),
    "--colormap", merged.get("colormap", "jet"),
    "--color-over", merged.get("color_over", "black"),
    "--color-under", merged.get("color_under", "white"),
    # hist2d always prefixes "figs/"; "../" lands the PNG in out_dir.
    "--output", "../" + out_name,
    "--no-csv",
  ]
  if use_log:
    argv.append("--log")
  args = h2.build_parser().parse_args(argv)
  h2.resolve_bins(args)
  args.title = src_name

  # In-memory (tx, ty, value) triplets instead of load_data(file).
  ele = panel["ele"]
  tx = 0.5 * (np.asarray(ele["txmin"], np.float64)
              + np.asarray(ele["txmax"], np.float64))
  ty = 0.5 * (np.asarray(ele["tymin"], np.float64)
              + np.asarray(ele["tymax"], np.float64))
  data = np.column_stack([tx, ty, TXTY_FIELDS[field](ele)])

  # histogram mode of hist2d.main(), unchanged.
  H, xedges, yedges = h2.compute_histogram(data, args)
  H_masked = np.ma.masked_where(H == 0, H)
  norm = h2.compute_norm(H.ravel(), args)
  X, Y = np.meshgrid(0.5 * (xedges[:-1] + xedges[1:]),
                     0.5 * (yedges[:-1] + yedges[1:]))
  cmap = h2.setup_colormap(args)
  contour_levels = h2.compute_contour_levels(args)

  cwd = os.getcwd()
  os.chdir(out_dir)  # render_plot writes figs/{output} relative to cwd
  try:
    h2.render_plot(X, Y, H_masked.T, xedges, yedges, norm, cmap,
                   contour_levels, args, mask_mesh=(X, Y))
  finally:
    os.chdir(cwd)
  # The empty figs/ directory hist2d creates inside out_dir is removed once in
  # main() after all jobs finish: removing it here would race parallel workers
  # that are about to write through the figs/../ relative path.
  out_path = Path(out_dir) / out_name
  print("  figure saved:", out_path)
  return out_path


def build_g2bg_dataframe(panel):
  """Rebuild the plot_g2bg per-group table from the binary g2bg data.

  Mirrors DetectorPanelArray::out_g2bg_all (src/cls_DetectorPanelArray.cpp:2196):
  one row per group, the rectangle from the group's member cell edges, and
  sig/noi summed over those cells. Returns (header, DataFrame). eff_* come
  from the vec_eff_*_group vectors serialized since PIPELINE_VERSION 7
  (-1.0 means the group efficiency was never computed).
  """
  import pandas as pd
  g2bg = panel["g2bg"]
  xa, ya = g2bg["x_axis"], g2bg["y_axis"]
  groups = g2bg["groups"]
  sig_cells = g2bg["vec_vec_signal"]
  noi_cells = g2bg["vec_vec_noise_poi"]
  avail_cells = g2bg["vec_vec_is_avail"]
  n_group = len(groups["dens_center"])
  rows = []
  for igroup in sorted(g2bg["bimap"]):
    cells = g2bg["bimap"][igroup]
    if not cells or igroup >= n_group:
      continue
    # rectangle = min/max of member cell edges (Grid2dBinGroup::get_xmin_xmax_ymin_ymax,
    # cell edge get_lower/upper_value(i) = min + i*interval / min + (i+1)*interval).
    xlow = min(xa["min"] + ix * xa["interval"] for ix, iy in cells)
    xup = max(xa["min"] + (ix + 1) * xa["interval"] for ix, iy in cells)
    ylow = min(ya["min"] + iy * ya["interval"] for ix, iy in cells)
    yup = max(ya["min"] + (iy + 1) * ya["interval"] for ix, iy in cells)
    # sig/noi summed over the group cells (calc_signal_group / calc_noise_poi_group);
    # get_signal(ix,iy) = vec_vec_signal[iy][ix] (src/cls_Grid2dBinGroup.cpp:828).
    sig = float(sum(sig_cells[iy][ix] for ix, iy in cells))
    noi = float(sum(noi_cells[iy][ix] for ix, iy in cells))
    ok = all(avail_cells[iy][ix] != 0 for ix, iy in cells)  # is_avail_group = AND
    rows.append({
      "ig": igroup, "xlow": xlow, "xup": xup, "ylow": ylow, "yup": yup,
      "sig": sig, "noi": noi, "signoi": sig + noi,
      "dens_low": float(groups["dens_lower"][igroup]),
      "dens_cnt": float(groups["dens_center"][igroup]),
      "dens_upp": float(groups["dens_upper"][igroup]),
      "delta_nmuon_low": float(groups["delta_nmuon_lower"][igroup]),
      "delta_nmuon_cnt": float(groups["delta_nmuon_center"][igroup]),
      "delta_nmuon_upp": float(groups["delta_nmuon_upper"][igroup]),
      "volume": float(groups["volume"][igroup]),
      "eff_low": float(groups["eff_low"][igroup]) if igroup < len(groups["eff_low"]) else np.nan,
      "eff_cnt": float(groups["eff_cnt"][igroup]) if igroup < len(groups["eff_cnt"]) else np.nan,
      "eff_upp": float(groups["eff_upp"][igroup]) if igroup < len(groups["eff_upp"]) else np.nan,
      "is_avail": 1 if ok else 0,
    })
  data = pd.DataFrame(rows)
  if not data.empty:
    data["dens_err"] = (data["dens_upp"] - data["dens_low"]) * 0.5
    data["eff_err"] = (data["eff_upp"] - data["eff_low"]) * 0.5
    data["sig_over_noi"] = data["sig"] / data["noi"].where(data["noi"] > 0.0)
  header = {"xmin": xa["min"], "xmax": xa["max"],
            "ymin": ya["min"], "ymax": ya["max"]}
  return header, data


def render_g2bg(out_dir, arr_name, panel, config):
  """Draw the g2bg group figures for one detector straight from the binary.

  Reuses plot_g2bg.plot_field so the PNGs match the normal g2bg flow. Field
  set and color axes come from the g2bg.template.plot_fields block of the same
  auto_plot.json5.
  """
  import matplotlib
  matplotlib.use("Agg")
  sys.path.insert(0, str(Path(__file__).resolve().parent))
  import plot_g2bg as g2

  template = config.get("g2bg", {}).get("template", {})
  plot_fields = template.get("plot_fields", {})
  png_dpi = template.get("png_dpi", 300)
  header, data = build_g2bg_dataframe(panel)
  if data.empty:
    print("  det%02d: no g2bg groups; skipped" % panel["detid"])
    return 0
  str_det_id = "%02d" % panel["detid"]
  output_prefix = "fig_%s_" % arr_name
  datafile = "%s_det%s (binary direct)" % (arr_name, str_det_id)
  n = 0
  for field, params in plot_fields.items():
    if not params.get("exec", True):
      continue
    g2.plot_field(field, output_prefix, str_det_id, params, header, data,
                  png_dpi, str(out_dir), datafile)
    n += 1
  return n


# Worker context for the drawing pool. Populated in main() before the pool is
# created; the fork start method lets workers inherit it, so the ~160 MB panel
# data is shared instead of being pickled per job.
_JOB_CTX = {}


def _run_render_job(job):
  """Render one figure job: ("txty", ipanel, field) or ("g2bg", ipanel)."""
  ctx = _JOB_CTX
  if job[0] == "txty":
    _, ipanel, field = job
    done = render_with_hist2d(ctx["out_dir"], ctx["arr_name"],
                              ctx["panels"][ipanel], field,
                              ctx["template"], ctx["patterns"])
    return 0 if done is None else 1
  _, ipanel = job
  return render_g2bg(ctx["out_dir"], ctx["arr_name"], ctx["panels"][ipanel],
                     ctx["config"])


def resolve_n_jobs(n_jobs):
  """Map the n_jobs convention to a worker count (>= 1).

  Positive: that many workers. 0: all CPU cores. Negative: leave that many
  cores free (e.g. -2 on a 32-core machine -> 30 workers).
  """
  n_cpu = os.cpu_count() or 1
  if n_jobs == 0:
    return n_cpu
  if n_jobs < 0:
    return max(1, n_cpu + n_jobs)
  return n_jobs


def check_against_mat(panels, mat_dir):
  """Compare det/ per-element scalars against mat/ element table, row for row."""
  sys.path.insert(0, str(Path(__file__).resolve().parent))
  from demo_recon_nnls import load_manifest, load_matrix

  mat_dir = Path(mat_dir)
  manifest = load_manifest(mat_dir)
  files = manifest["files"]
  by_detid = {p["detid"]: p for p in panels}
  # det/ record field -> mat/ element-table column name.
  pairs = [("PL", "path_length_m"), ("DL", "density_length_kg_m2"),
           ("signal", "signal"), ("noise_det", "noise_det"),
           ("noise_poi", "noise_poi"), ("proj_density", "proj_density_kg_m3")]
  print("[check-mat] comparing det/ direct read vs mat/ element table")
  worst = 0.0
  for js_det in manifest["detectors"]:
    detid = js_det["detid"]
    if detid not in by_detid or "file_element_table" not in js_det:
      continue
    name = js_det["file_element_table"]
    mm, n_ele, n_col = load_matrix(mat_dir / name)
    arr = np.asarray(mm, dtype=np.float64)
    columns = files[name]["columns"]
    col_idx = {c: j for j, c in enumerate(columns)}
    ele = by_detid[detid]["ele"]
    if len(ele) != n_ele:
      print("  det%02d SIZE MISMATCH: det=%d rows, mat=%d rows"
            % (detid, len(ele), n_ele))
      worst = float("inf")
      continue
    msgs = []
    for det_field, mat_col in pairs:
      a = np.asarray(ele[det_field], np.float64)
      b = arr[:, col_idx[mat_col]]
      denom = np.maximum(np.abs(b), 1e-30)
      rel = float(np.max(np.abs(a - b) / denom))
      worst = max(worst, rel)
      msgs.append("%s %.2e" % (det_field, rel))
    print("  det%02d max rel diff: %s" % (detid, ", ".join(msgs)))
  tol = 1e-6  # float32 stores ~7 significant digits
  verdict = "OK" if worst <= tol else "MISMATCH"
  print("[check-mat] worst rel diff = %.3e (tol %.0e) -> %s" % (worst, tol, verdict))
  return worst <= tol


def main(argv=None):
  parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
  parser.add_argument("source",
                      help="arrdet_*.bin, a checkpoint dir, or a det/ dir")
  parser.add_argument("--field", default="PL,dens,signal",
                      help="comma-separated txty fields to plot "
                           "(default: PL,dens,signal — the same trio the C++ "
                           "exports; see --list-fields)")
  parser.add_argument("--det", default=None,
                      help="comma-separated detids to keep (default: all)")
  parser.add_argument("--out", default=None,
                      help="output directory (default: det.out_dir from the "
                           "auto_plot config, resolved against the config's "
                           "folder; else <source_dir>/plot_det_arrdet)")
  parser.add_argument("--jobs", type=int, default=None,
                      help="figure-drawing workers: positive = that many, "
                           "1 = one by one, 0 = all CPU cores, negative = "
                           "leave that many cores free (default: det.n_jobs "
                           "from the auto_plot config, else -2)")
  parser.add_argument("--check-mat", default=None,
                      help="mat/ directory; validate det/ values against the "
                           "element table row for row")
  parser.add_argument("--plot-config", default=None,
                      help="auto_plot.json5 for the hist2d parameters "
                           "(default: nearest auto_plot.json5 above the "
                           "arrdet .bin, else scripts/auto_plot_default.json5)")
  parser.add_argument("--g2bg", action="store_true",
                      help="also draw the g2bg group figures (sig/noi/signoi/"
                           "dens/dens_err/... straight from the binary). By "
                           "default g2bg follows the config's g2bg_exec switch")
  parser.add_argument("--no-g2bg", action="store_true",
                      help="never draw g2bg figures, even if g2bg_exec is true")
  parser.add_argument("--no-hist2d", action="store_true",
                      help="never draw the per-element txty (hist2d) figures, "
                           "even if hist2d_exec is true")
  parser.add_argument("--list-fields", action="store_true",
                      help="print the plottable field names and exit")
  args = parser.parse_args(argv)

  if args.list_fields:
    print("plottable txty fields (same names as the C++ out_txty* suffixes):")
    for k in TXTY_FIELDS:
      print("  %s" % k)
    return 0

  fields = [x.strip() for x in args.field.split(",") if x.strip()]
  unknown = [x for x in fields if x not in TXTY_FIELDS]
  if unknown:
    print("ERROR: unknown field(s) %s; choose from: %s"
          % (unknown, ", ".join(TXTY_FIELDS)))
    return 1

  arrdet_path = resolve_arrdet_path(args.source)
  print("[1] Read", arrdet_path)
  arr_name, panels = read_arrdet(arrdet_path)
  for p in panels:
    print("  det%02d: %d x %d elements (n_element=%d)"
          % (p["detid"], p["nbiny"], p["nbinx"], p["n_element"]))
  print("  %d detectors loaded" % len(panels))

  if args.det is not None:
    keep = {int(x) for x in args.det.split(",")}
    panels = [p for p in panels if p["detid"] in keep]
    if not panels:
      print("ERROR: no detectors match --det %s" % args.det)
      return 1

  ok = True
  if args.check_mat:
    ok = check_against_mat(panels, args.check_mat)

  config_path = resolve_plot_config(arrdet_path, args.plot_config)
  config, template, patterns = load_hist2d_config(config_path)
  print("[2] plot parameters from", config_path)
  sys.path.insert(0, str(Path(__file__).resolve().parent))
  from auto_plot import is_excluded

  # Output directory: --out beats det.out_dir from the config; a relative
  # out_dir is anchored at the config's folder so the figures land in the work
  # dir's figs/ no matter where the .bin sits.
  det_cfg = config.get("det", {})
  if args.out:
    out_dir = Path(args.out)
  elif det_cfg.get("out_dir"):
    out_dir = Path(det_cfg["out_dir"])
    if not out_dir.is_absolute():
      out_dir = Path(config_path).resolve().parent / out_dir
  else:
    out_dir = arrdet_path.parent / "plot_det_arrdet"
  out_dir.mkdir(parents=True, exist_ok=True)

  jobs = []

  # Per-element txty (hist2d) figures. Honor the same switches as auto_plot.py:
  # the top-level hist2d_exec gate and the hist2d exclude_patterns list.
  do_hist2d = config.get("hist2d_exec", True) and not args.no_hist2d
  if not do_hist2d:
    print("[3] hist2d figures disabled (hist2d_exec=false or --no-hist2d)")
  else:
    for field in fields:  # one PNG per detector per field, standard txty style
      for ipanel, panel in enumerate(panels):
        src_name = "%s_txty_%s_det%02d.tmp" % (arr_name, field, panel["detid"])
        if is_excluded(src_name, config, "hist2d", []):
          print("  [exclude] %s" % src_name)
          continue
        jobs.append(("txty", ipanel, field))
    print("[3] txty figure jobs: %d (one per detector per field)" % len(jobs))

  # g2bg group figures (binary-direct). Default follows the config g2bg_exec
  # switch; --g2bg forces on, --no-g2bg forces off. Same exclude_patterns apply.
  do_g2bg = (not args.no_g2bg) and (args.g2bg or config.get("g2bg_exec", True))
  if do_g2bg:
    n_txty_jobs = len(jobs)
    for ipanel, panel in enumerate(panels):
      gname = "g2bg_%s_det%02d.tmp" % (arr_name, panel["detid"])
      if is_excluded(gname, config, "g2bg", []):
        print("  [exclude] %s" % gname)
        continue
      jobs.append(("g2bg", ipanel))
    print("[4] g2bg figure jobs: %d (one per detector)"
          % (len(jobs) - n_txty_jobs))
  else:
    print("[4] g2bg figures disabled (g2bg_exec=false or --no-g2bg)")

  # --jobs beats det.n_jobs from the config; without either, -2 leaves two
  # CPU cores free. Workers inherit _JOB_CTX through fork (no per-job pickling
  # of the panel data).
  n_workers = resolve_n_jobs(args.jobs if args.jobs is not None
                             else det_cfg.get("n_jobs", -2))
  _JOB_CTX.update(out_dir=out_dir, arr_name=arr_name, panels=panels,
                  template=template, patterns=patterns, config=config)
  if n_workers > 1 and len(jobs) > 1:
    import multiprocessing as mp
    try:
      pool_ctx = mp.get_context("fork")
    except ValueError:  # platform without fork: fall back to one-by-one
      pool_ctx = None
  else:
    pool_ctx = None
  if pool_ctx is None:
    results = [_run_render_job(job) for job in jobs]
  else:
    n_workers = min(n_workers, len(jobs))
    print("[5] drawing %d figures with %d workers -> %s"
          % (len(jobs), n_workers, out_dir))
    with pool_ctx.Pool(n_workers) as pool:
      results = pool.map(_run_render_job, jobs)
  n_fig = sum(results)
  try:  # drop the empty figs/ directory hist2d creates inside out_dir
    (out_dir / "figs").rmdir()
  except OSError:
    pass

  print("DONE (%d figures)" % n_fig)
  return 0 if ok else 2


if __name__ == "__main__":
  sys.exit(main())
