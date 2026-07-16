try:
  import pandas as pd
except ImportError:
  import sys
  sys.exit(
    "ERROR: pandas is required.\n"
    "  pip install pandas  (or activate venv: source .venv/bin/activate)"
  )
import re
import json
import argparse
import os
import sys
import struct
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import copy
import subprocess
from concurrent.futures import ProcessPoolExecutor
from matplotlib.ticker import MaxNLocator
import shutil
try:
  import json5
except ImportError:
  json5 = None
  import warnings
  warnings.warn(
    "json5 module not found. JSON5 config files will fail to parse.\n"
    "  Hint: activate the venv — source .venv/bin/activate",
    stacklevel=2,
  )

def _sanitize_run_tag(text, fallback="run"):
  if not text:
    return fallback
  slug = re.sub(r'[^A-Za-z0-9._-]+', '_', str(text))
  slug = slug.strip("._-")
  return slug or fallback

def _resolve_output_directory(config_path, run_params, render_params):
  user_tag = run_params.get("run_tag") or render_params.get("run_tag")
  if not user_tag and config_path:
    config_base = os.path.splitext(os.path.basename(config_path))[0]
    user_tag = config_base
  run_tag = _sanitize_run_tag(user_tag)
  base_dir = render_params.get('save_png_dir')
  if base_dir:
    base_dir = os.path.abspath(base_dir)
  else:
    base_dir = os.path.abspath(render_params.get('output_png_prefix', run_tag))
  os.makedirs(base_dir, exist_ok=True)
  candidate_dir = os.path.join(base_dir, run_tag)
  if os.path.exists(candidate_dir) and not os.path.isdir(candidate_dir):
    raise NotADirectoryError(f"[ERROR] {candidate_dir} exists but is not a directory.")
  unique_dir = candidate_dir
  os.makedirs(unique_dir, exist_ok=True)
  render_params['run_tag'] = run_tag
  render_params['save_png_dir_base'] = base_dir
  render_params['save_png_dir'] = unique_dir
  render_params['run_output_dir'] = unique_dir
  print(f"[INFO] Output directory resolved: {unique_dir}")
  return unique_dir

# Load config file
def load_config(config_file):
  with open(config_file, 'r') as f:
    content = f.read()
  if json5 is not None:
    return json5.loads(content)
  if config_file.endswith(".json5"):
    raise ImportError(
      f"Cannot parse JSON5 file without json5 module: {config_file}\n"
      "  Install: pip install json5\n"
      "  Or activate venv: source .venv/bin/activate"
    )
  return json.loads(content)

# Magic number for identifying binary files
MAGIC_BINARY = b"G2ZBIN\x00\x00"

def _read_exact(handle, size):
  data = handle.read(size)
  if len(data) != size:
    raise ValueError(f"Unexpected EOF while reading {size} bytes")
  return data

def parse_cross_section_spec(cross_section):
  if not isinstance(cross_section, str) or not cross_section.strip():
    raise ValueError("cross_section must be a non-empty string")
  spec = cross_section.strip()
  descending = False
  if spec[0] in "+-":
    descending = spec[0] == '-'
    spec = spec[1:]
  axis = spec.lower()
  if axis not in {"x", "y", "z"}:
    raise ValueError(f"Invalid cross_section specification: {cross_section}")
  return axis, descending

# Read header and data lines from a text file
def read_text_file(filename):
  # Extract header lines (start with '#')
  header_info = {}
  with open(filename, 'r') as f:
    for line in f:
      line = line.strip()
      if not line:
        continue
      if line.startswith("#"):
        parts = re.split(r'\s+', line[1:].strip())
        key = parts[0]
        values = []
        for v in parts[1:]:
          try:
            values.append(float(v))
          except ValueError:
            values.append(v)
        header_info[key] = values
      else:
        break  # Stop at first data line
  header_info.setdefault("format", ["text"])

  # Bulk parse data with pandas C engine
  df_raw = pd.read_csv(filename, comment='#', sep=r'\s+',
                        header=None, dtype=np.float64, engine='c')
  data_lines = df_raw.values.tolist()
  return header_info, data_lines

def read_binary_file(filename):
  header_info = {}
  with open(filename, 'rb') as f:
    magic = _read_exact(f, len(MAGIC_BINARY))
    if magic != MAGIC_BINARY:
      raise ValueError("Invalid magic number for binary cross section file")

    version = struct.unpack('<I', _read_exact(f, 4))[0]
    _reserved = struct.unpack('<I', _read_exact(f, 4))[0]

    double_vals = struct.unpack('<9d', _read_exact(f, 9 * 8))
    xmin, xmax, xstep, ymin, ymax, ystep, zmin, zmax, zstep = double_vals

    int_vals = struct.unpack('<4i', _read_exact(f, 4 * 4))
    n_detector, nbinx, nbiny, nz = int_vals

    n_records = struct.unpack('<q', _read_exact(f, 8))[0]

    header_info["format"] = ["binary"]
    header_info["version"] = [version]
    header_info["x_info"] = [xmin, xmax, xstep]
    header_info["y_info"] = [ymin, ymax, ystep]
    header_info["z_info"] = [zmin, zmax, zstep]
    header_info["n_detector"] = [n_detector]
    header_info["nbinx"] = [nbinx]
    header_info["nbiny"] = [nbiny]
    header_info["nz"] = [nz]
    header_info["n_records"] = [n_records]

    # Determine actual per-record flag count from file size.
    # The header n_detector may differ from per-record n_det when hit flags
    # are stored at grid level (n_det=0 per record).
    data_start = f.tell()
    f.seek(0, 2)
    data_size = f.tell() - data_start
    f.seek(data_start)

    fixed_size = 1 + 8 + 8 + 8 + 8 + 4 + 4  # tf_exist + xyz + density + n_det + n_hit_det
    if n_records > 0:
      actual_record_size = data_size // n_records
      n_flags = actual_record_size - fixed_size
    else:
      n_flags = n_detector

    if n_flags == n_detector:
      record_dtype = np.dtype([
        ('tf_exist',  'u1'),
        ('xcnt',      '<f8'),
        ('ycnt',      '<f8'),
        ('zcnt',      '<f8'),
        ('density',   '<f8'),
        ('n_det_row', '<i4'),
        ('flags',     'u1', (n_detector,)),
        ('n_hit_det', '<i4'),
      ])
    elif n_flags == 0:
      record_dtype = np.dtype([
        ('tf_exist',  'u1'),
        ('xcnt',      '<f8'),
        ('ycnt',      '<f8'),
        ('zcnt',      '<f8'),
        ('density',   '<f8'),
        ('n_det_row', '<i4'),
        ('n_hit_det', '<i4'),
      ])
    else:
      raise ValueError(
        f"Unexpected record size: {actual_record_size} bytes "
        f"(fixed={fixed_size}, n_flags={n_flags}, header n_detector={n_detector})")

    raw = f.read(n_records * record_dtype.itemsize)
    if len(raw) != n_records * record_dtype.itemsize:
      raise ValueError(f"Unexpected EOF: expected {n_records * record_dtype.itemsize} bytes, got {len(raw)}")
    records = np.frombuffer(raw, dtype=record_dtype)

    # Build dict of arrays (consumed by main() to construct DataFrame)
    data = {
      'tf_exist': records['tf_exist'].astype(np.float64),
      'xcnt':     records['xcnt'].copy(),
      'ycnt':     records['ycnt'].copy(),
      'zcnt':     records['zcnt'].copy(),
      'density':  records['density'].copy(),
      'n_det':    records['n_det_row'].astype(np.float64),
    }
    if n_flags > 0:
      for i in range(n_flags):
        data[f'flag_{i}'] = records['flags'][:, i].astype(np.float64)
    else:
      for i in range(n_detector):
        data[f'flag_{i}'] = np.zeros(n_records, dtype=np.float64)
    data['n_hit_det'] = records['n_hit_det'].astype(np.float64)

  return header_info, data

def read_cross_section_file(filename):
  try:
    with open(filename, 'rb') as f:
      prefix = f.read(len(MAGIC_BINARY))
  except FileNotFoundError:
    raise

  if prefix == MAGIC_BINARY:
    print("[INFO] Detected cross-section format: binary (start reading)", flush=True)
    return read_binary_file(filename)
  print("[INFO] Detected cross-section format: text (start reading)", flush=True)
  return read_text_file(filename)

# Filter by n_hit_det and (optionally) detector flags
def filter_plot_data(df, run_params, det_id=None):
  df_filtered = df[df['n_hit_det'] >= run_params['n_hit_det_thres']]
  if det_id is not None:
    flag_col = f'flag_{det_id}'
    df_filtered = df_filtered[df_filtered[flag_col] == 1]
  return df_filtered

# Utilities for numeric sequences and z-slice determination
def _generate_sequence(start, stop, step, tol=1e-6, max_iter=1000000):
  if step is None:
    raise ValueError("step is None")
  step = float(step)
  if step <= 0:
    raise ValueError("step must be positive.")
  start = float(start)
  stop = float(stop)
  values = []
  for idx in range(max_iter):
    current = start + idx * step
    if current > stop + tol:
      break
    values.append(current)
  else:
    raise ValueError("sequence generation exceeded maximum iterations")
  return values

def _match_targets_to_data(targets, data_values, tol=1e-6):
  matched = []
  data_index = 0
  data_len = len(data_values)
  for target in targets:
    found = False
    while data_index < data_len:
      candidate = data_values[data_index]
      diff = candidate - target
      if abs(diff) <= tol:
        matched.append(candidate)
        data_index += 1
        found = True
        break
      if candidate < target - tol:
        data_index += 1
        continue
      break
    if not found:
      return None
  return matched

def determine_slice_values(df, run_params, cross_section, tol=1e-6):
  axis, _descending = parse_cross_section_spec(cross_section)

  # use_header_bins or use_data_bins mode: get slice values directly from data
  # (these modes compute bin edges, but slice values should come from actual data coordinates)
  if run_params.get("use_header_bins", False) or run_params.get("use_data_bins", False):
    if axis == "z":
      return sorted(df['zcnt'].unique())
    elif axis == "x":
      return sorted(df['xcnt'].unique())
    elif axis == "y":
      return sorted(df['ycnt'].unique())
    else:
      raise ValueError(f"Unknown cross_section: {cross_section}")

  # --- auto_binning: false (default) ---
  col = {"z": "zcnt", "x": "xcnt", "y": "ycnt"}[axis]
  all_slices = sorted(df[col].unique())

  smin = run_params.get(f"{axis}min")
  smax = run_params.get(f"{axis}max")
  if smin is not None:
    all_slices = [v for v in all_slices if v >= float(smin) - 1e-6]
  if smax is not None:
    all_slices = [v for v in all_slices if v <= float(smax) + 1e-6]

  if not all_slices:
    print(f"[WARN] No slices found in range {axis}=[{smin}, {smax}]. "
          f"Data range: {df[col].min():.1f} - {df[col].max():.1f}",
          file=sys.stderr)

  return all_slices

# Function for creating cross-section plots
from matplotlib.ticker import MaxNLocator

def plot_data_for_slice(df_slice, run_params, render_params, slice_value, cross_section, det_id=None, det_xy=None):
  # Determine the fixed column, histogram axes, and bin settings based on cross_section
  axis, _descending = parse_cross_section_spec(cross_section)
  # Axis mapping: fixed column, histogram axes
  axis_map = {
    "z": ("zcnt", "xcnt", "ycnt"),
    "x": ("xcnt", "ycnt", "zcnt"),
    "y": ("ycnt", "xcnt", "zcnt"),
  }
  if axis not in axis_map:
    raise ValueError("Unknown cross_section value: must be one of 'x', 'y', 'z'")
  fixed_col, hist_col1, hist_col2 = axis_map[axis]

  # Get step from run_params or header fallback (may be None if no step info)
  def _get_step(a):
    return run_params.get(f"{a}step", run_params.get(f"header_{a}step"))
  a1, a2 = hist_col1[0], hist_col2[0]
  bin1_min = run_params.get(f"{a1}min", 0)
  bin1_max = run_params.get(f"{a1}max", 0)
  bin1_step = _get_step(a1)
  bin2_min = run_params.get(f"{a2}min", 0)
  bin2_max = run_params.get(f"{a2}max", 0)
  bin2_step = _get_step(a2)

  if bin1_step is not None:
    bin1_step = abs(float(bin1_step))
  if bin2_step is not None:
    bin2_step = abs(float(bin2_step))

  # Get xy_unit setting (default: "km")
  xy_unit = render_params.get('xy_unit', 'km').lower()
  unit_scale = 0.001 if xy_unit == 'km' else 1.0  # Convert meters to km if needed
  unit_label = 'km' if xy_unit == 'km' else 'meters'

  # df_slice is already pre-filtered by slice_value (done in caller)
  df_slice_filtered = filter_plot_data(df_slice, run_params, det_id=det_id)
  
  # Bin edges: use header-derived edges if available, otherwise np.arange
  # hist_col1 is "xcnt", "ycnt", etc. hist_col1[0] gives "x", "y", etc.
  bin1_key = f"{hist_col1[0]}_bin_edges"
  bin2_key = f"{hist_col2[0]}_bin_edges"
  if bin1_key in run_params:
    bin1_edges = run_params[bin1_key]
    if bin1_step is None:
      bin1_step = float(bin1_edges[1] - bin1_edges[0]) if len(bin1_edges) > 1 else 1.0
  else:
    bin1_edges = np.arange(bin1_min, bin1_max + bin1_step, bin1_step)
  if bin2_key in run_params:
    bin2_edges = run_params[bin2_key]
    if bin2_step is None:
      bin2_step = float(bin2_edges[1] - bin2_edges[0]) if len(bin2_edges) > 1 else 1.0
  else:
    bin2_edges = np.arange(bin2_min, bin2_max + bin2_step, bin2_step)

  # Create histogram (fill with zeros if there is no data)
  if df_slice_filtered.empty:
    H = np.zeros((len(bin1_edges)-1, len(bin2_edges)-1))
  else:
    x_vals = df_slice_filtered[hist_col1].values
    y_vals = df_slice_filtered[hist_col2].values
    weights = df_slice_filtered['density'].values
    H, _, _ = np.histogram2d(x_vals, y_vals, bins=[bin1_edges, bin2_edges], weights=weights)
  H = np.ma.masked_where(H == 0, H)

  # Apply unit scaling for display (bin edges for plotting)
  bin1_edges_display = bin1_edges * unit_scale
  bin2_edges_display = bin2_edges * unit_scale
  # Display range: use JSON min/max if present, otherwise fall back to bin edges
  bin1_min_display = run_params.get(f"{hist_col1[0]}min", bin1_edges[0])
  if isinstance(bin1_min_display, (int, float)):
    bin1_min_display = float(bin1_min_display) * unit_scale
  bin1_max_display = run_params.get(f"{hist_col1[0]}max", bin1_edges[-1])
  if isinstance(bin1_max_display, (int, float)):
    bin1_max_display = float(bin1_max_display) * unit_scale
  bin2_min_display = run_params.get(f"{hist_col2[0]}min", bin2_edges[0])
  if isinstance(bin2_min_display, (int, float)):
    bin2_min_display = float(bin2_min_display) * unit_scale
  bin2_max_display = run_params.get(f"{hist_col2[0]}max", bin2_edges[-1])
  if isinstance(bin2_max_display, (int, float)):
    bin2_max_display = float(bin2_max_display) * unit_scale
  bin1_step_display = bin1_step * unit_scale
  bin2_step_display = bin2_step * unit_scale

  # Use cell edges as-is for the mesh grid (display units)
  X, Y = np.meshgrid(bin1_edges_display, bin2_edges_display)
  
  # Figure size: compute the vertical size from the data ratio relative to main_width (horizontal size)
  main_width = render_params.get('main_width', 10)
  # Compute the aspect ratio from the bin edge range
  data_aspect = (bin2_edges[-1] - bin2_edges[0]) / (bin1_edges[-1] - bin1_edges[0])
  main_height = main_width * data_aspect
  cbar_width = render_params.get('cbar_width', 0.5)
  total_width = main_width + cbar_width
  
  title_fontsize = render_params.get('title_fontsize', 18)
  label_fontsize = render_params.get('xy_label_fontsize', 16)
  tick_fontsize  = render_params.get('xy_tick_fontsize', 14)
  
  fig = plt.figure(figsize=(total_width, main_height))
  gs = gridspec.GridSpec(1, 2, width_ratios=[main_width, cbar_width], wspace=0.05)
  ax = fig.add_subplot(gs[0])
  cbar_ax = fig.add_subplot(gs[1])
  
  ax.set_xlim(bin1_min_display, bin1_max_display)
  ax.set_ylim(bin2_min_display, bin2_max_display)
  ax.set_facecolor(render_params.get('bg_color', 'silver'))
  
  cmap = copy.copy(plt.get_cmap(render_params['colormap']))
  cmap.set_under(render_params['underflow_color'])
  cmap.set_over(render_params['overflow_color'])
  bounds = np.linspace(render_params['vmin'], render_params['vmax'], render_params['ngrad'] + 1)
  norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='both')
  
  ax.grid(False)
  mesh = ax.pcolormesh(X, Y, H.T, cmap=cmap, norm=norm, shading='auto')
  cbar = fig.colorbar(mesh, cax=cbar_ax, extend='both')
  cbar.set_label('density', fontsize=render_params.get('cbar_label_fontsize', 16))
  cbar.ax.tick_params(labelsize=render_params.get('cbar_tick_fontsize', 14))
  # Optional fixed colorbar tick spacing (e.g. 100 -> 0,100,...,vmax).
  cbar_tick_step = render_params.get('cbar_tick_step')
  if cbar_tick_step:
    cbar.set_ticks(np.arange(render_params['vmin'],
                             render_params['vmax'] + cbar_tick_step * 0.5,
                             cbar_tick_step))

  # Ticks use the edge values as-is, but adjust with MaxNLocator if there are too many ticks
  ax.xaxis.set_major_locator(MaxNLocator(nbins=10))
  ax.yaxis.set_major_locator(MaxNLocator(nbins=10))
  
  # Or, when using bin1_edges, bin2_edges directly (commented out below)
  # ax.set_xticks(bin1_edges)
  # ax.set_yticks(bin2_edges)
  
  def build_axis_label(column_name, base_label, fallback_step_display, display_unit):
    prefix = column_name.lower()[0] if column_name else ""
    binsize_configs = {
      'x': (('binsize_x', 'xstep', 'x_pitch'), 'dx'),
      'y': (('binsize_y', 'ystep', 'y_pitch'), 'dy'),
    }
    if prefix not in binsize_configs:
      return f"{base_label} [{display_unit}]"
    candidate_keys, suffix = binsize_configs[prefix]
    for key in candidate_keys:
      if key in run_params:
        try:
          value = float(run_params[key]) * unit_scale
          return f"{base_label} [{display_unit}] ({suffix}={value:.3f} {display_unit})"
        except (TypeError, ValueError):
          continue
    if fallback_step_display is not None:
      try:
        value = float(fallback_step_display)
        return f"{base_label} [{display_unit}] ({suffix}={value:.3f} {display_unit})"
      except (TypeError, ValueError):
        pass
    return f"{base_label} [{display_unit}]"

  xlabel = build_axis_label(hist_col1, hist_col1, bin1_step_display, unit_label)
  ylabel = build_axis_label(hist_col2, hist_col2, bin2_step_display, unit_label)
  ax.set_xlabel(xlabel, fontsize=label_fontsize)
  ax.set_ylabel(ylabel, fontsize=label_fontsize)
  title_suffix = f" det {det_id:03d}" if det_id is not None else ""
  ax.set_title(f"{render_params['output_png_prefix']} {cross_section}={slice_value:05.0f}{title_suffix}", fontsize=title_fontsize)
  ax.tick_params(axis='both', labelsize=tick_fontsize)
  
  # Mask overlay (df_slice is unfiltered by n_hit_det, which is correct for mask)
  df_mask = df_slice
  if not df_mask.empty:
    H_mask, _, _ = np.histogram2d(df_mask[hist_col1].values, df_mask[hist_col2].values,
                                  bins=[bin1_edges, bin2_edges],
                                  weights=df_mask['tf_exist'].values)
    overlay_mask = np.where(H_mask == 0, 0, np.nan)
    white_cmap = mpl.colors.ListedColormap([render_params.get('mask_color', 'white')])
    ax.pcolormesh(X, Y, np.ma.masked_invalid(overlay_mask).T,
                  cmap=white_cmap, shading='auto',
                  alpha=render_params.get('mask_alpha', 0.3))
  

  # Draw grid lines (using render_params' grid_color, gridline_width)
  ax.grid(True, color=render_params.get('grid_color', 'grey'), linewidth=render_params.get('gridline_width', 1.5), linestyle=render_params.get('gridline_style', '--'))

  # Star marker for detector position (only when cross=="z", as an example)
  if det_xy is not None and axis == "z":
    show_name = render_params.get('det_name_show', True)
    name_fontsize = render_params.get('det_name_fontsize', 10)
    name_offset_y = render_params.get('det_name_offset_y', 0.02)
    name_color = render_params.get('det_name_color', 'black')
    if det_id is None:
      for det_info in det_xy:
        dx, dy = det_info[0], det_info[1]
        dname = det_info[2] if len(det_info) > 2 else ""
        if dx is not None and dy is not None:
          ax.plot(dx * unit_scale, dy * unit_scale,
                  marker=render_params.get('det_pos_marker_type', '*'),
                  markersize=render_params.get('det_pos_marker_size', 16),
                  color=render_params.get('det_pos_marker_color', 'black'))
          if show_name and dname:
            ax.text(dx * unit_scale, dy * unit_scale + name_offset_y, dname,
                    fontsize=name_fontsize, color=name_color,
                    ha='center', va='bottom')
    else:
      if det_id < len(det_xy):
        det_info = det_xy[det_id]
        dx, dy = det_info[0], det_info[1]
        dname = det_info[2] if len(det_info) > 2 else ""
        if dx is not None and dy is not None:
          ax.plot(dx * unit_scale, dy * unit_scale,
                  marker=render_params.get('det_pos_marker_type', '*'),
                  markersize=render_params.get('det_pos_marker_size', 16),
                  color=render_params.get('det_pos_marker_color', 'black'))
          if show_name and dname:
            ax.text(dx * unit_scale, dy * unit_scale + name_offset_y, dname,
                    fontsize=name_fontsize, color=name_color,
                    ha='center', va='bottom')
  
  # set_box_aspect: use display range (may differ from full bin edges)
  box_aspect = (bin2_max_display - bin2_min_display) / (bin1_max_display - bin1_min_display)
  ax.set_box_aspect(box_aspect)
  
  return fig, ax

# Function for the cross-section creation task (for parallel execution)
def plot_save_task_cross(args_tuple):
  slice_value, cross_section, det_id, run_params, render_params, df_slice, det_xy = args_tuple
  fig, ax = plot_data_for_slice(df_slice, run_params, render_params, slice_value, cross_section, det_id, det_xy)
  save_plot_cross(fig, slice_value, cross_section, run_params, render_params, det_id)
  return (slice_value, cross_section, det_id)

# Parallel processing function for cross_section
def plot_and_save_all_cross_sections(df, run_params, render_params, n_detector, det_xy, cross_section):
  tasks = []
  unique_vals = determine_slice_values(df, run_params, cross_section)

  # * Align the generation order with the descending-order setting (for consistency in visible logs and sequential saves)
  _axis, desc = parse_cross_section_spec(cross_section)
  ordered_vals = list(reversed(unique_vals)) if desc else list(unique_vals)

  # Pre-split DataFrame by slice value to avoid pickling the full df per worker
  fixed_col = {"z": "zcnt", "x": "xcnt", "y": "ycnt"}[_axis]

  for val in ordered_vals:
    df_slice = df[df[fixed_col] == val]
    tasks.append((val, cross_section, None, run_params, render_params, df_slice, det_xy))
  if run_params.get('tf_out_each_det', False):
    for val in ordered_vals:
      df_slice = df[df[fixed_col] == val]
      for det_id in range(n_detector):
        tasks.append((val, cross_section, det_id, run_params, render_params, df_slice, det_xy))

  with ProcessPoolExecutor() as executor:
    futures = [executor.submit(plot_save_task_cross, task) for task in tasks]
    for future in futures:
      try:
        result = future.result()
        print(f"Completed plot for {result[1]}={result[0]}, det_id={result[2]}")
      except Exception as e:
        print(f"Error in plotting: {e}")
  return unique_vals

def fmt(val,offset=10000):
  # Add offset 1000, convert to a 6-digit zero-padded string
  return f"{val + offset:06.0f}"

def create_filename_list_cross(df, run_params, render_params, n_detector, cross_section, unique_vals=None):
  if unique_vals is None:
    unique_vals = determine_slice_values(df, run_params, cross_section)
  _axis, descending = parse_cross_section_spec(cross_section)
  ordered_vals = list(reversed(unique_vals)) if descending else list(unique_vals)
  out_dir = render_params.get('save_png_dir', render_params['output_png_prefix'])
  file_lists = {}
  # The PNG prefix is always "fig_"
  png_prefix = "fig_"
  file_list_all = [f"{out_dir}/{png_prefix}{fmt(val)}_{cross_section}={val:.0f}.png" for val in ordered_vals]
  file_lists['all'] = file_list_all
  if run_params.get('tf_out_each_det', False):
    for det_id in range(n_detector):
      file_list = [f"{out_dir}/{png_prefix}{fmt(val)}_{cross_section}={val:.0f}_det{det_id:03d}.png" for val in ordered_vals]
      file_lists[det_id] = file_list
  return file_lists

def save_plot_cross(fig, slice_value, cross_section, run_params, render_params, det_id=None):
  out_dir = render_params.get('save_png_dir', render_params['output_png_prefix'])
  if not os.path.exists(out_dir):
    os.makedirs(out_dir)
  png_prefix = "fig_"
  if det_id is None:
    filename = f"{out_dir}/{png_prefix}{fmt(slice_value)}_{cross_section}={slice_value:.0f}.png"
  else:
    filename = f"{out_dir}/{png_prefix}{fmt(slice_value)}_{cross_section}={slice_value:.0f}_det{det_id:03d}.png"
  plt.savefig(filename, dpi=render_params.get('dpi', 300))
  plt.close(fig)
  print(f"Saved plot to {filename}")

# Save the filename list for cross_section
def save_filename_list_to_file_cross(file_lists, render_params, cross_section):
  out_dir = render_params.get('save_png_dir', render_params['output_png_prefix'])
  filepath = os.path.join(out_dir, f"file_list_{cross_section}.tmp")
  with open(filepath, 'w') as f:
    for fname in file_lists['all']:
      f.write(f"{fname}\n")
  print(f"Saved file list for {cross_section} to {filepath}")
  for key in file_lists:
    if key == 'all':
      continue
    filepath = os.path.join(out_dir, f"file_list_{cross_section}_det_{key:03d}.tmp")
    with open(filepath, 'w') as f:
      for fname in file_lists[key]:
        f.write(f"{fname}\n")
    print(f"Saved file list for {cross_section} det_id {key} to {filepath}")

def _should_reverse_from_steps(run_params, cross_section):
  axis, _desc = parse_cross_section_spec(cross_section)
  key = {"x": "xstep", "y": "ystep", "z": "zstep"}.get(axis)
  if key is None:
    return False
  try:
    return float(run_params.get(key, 0)) < 0
  except (TypeError, ValueError):
    return False

# Run gifski for cross_section
def run_gifski_cross(file_lists, render_params, run_params, cross_section):
  out_dir = render_params.get('save_png_dir', render_params['output_png_prefix'])
  gifski_executable = render_params.get('gifski_executable', 'gifski')
  gifski_params = render_params.get("gifski_params", {})
  # Debug output: make explicit the final order passed to GIF (to rule out suspected internal sorting by the tool)
  def _debug_dump(tag, files):
    if files:
      print(f"[DEBUG] {tag}: first={os.path.basename(files[0])}, last={os.path.basename(files[-1])}, n={len(files)}")
  def build_command(output_gif, input_files):
    cmd = [gifski_executable]
    cmd.extend([
      "--fps", str(gifski_params.get("fps", 0.5)),
      "--quality", str(gifski_params.get("quality", 100))
    ])
    # Preserve input PNG resolution (gifski defaults to ~800x600)
    if input_files:
      try:
        from PIL import Image as _PILImage
        with _PILImage.open(input_files[0]) as _img:
          cmd.extend(["--width", str(_img.width), "--height", str(_img.height)])
      except Exception:
        pass
    cmd.extend(["-o", output_gif])
    # Explicitly mark the end of options with "--". Avoids internal misinterpretation and future compatibility issues
    cmd.append("--")
    cmd.extend(input_files)
    return cmd
  # Fix the playback order: create sequentially numbered links (or copies) like frame_000000.png ... in a temporary directory
  # Returns: (realized_files, tmp_dir) - tmp_dir is for cleanup after gifski
  def materialize_sequenced_frames(tag, src_files):
    if not src_files:
      return [], None
    tmp_dir = os.path.join(out_dir, f"_gifski_{cross_section}_{tag}_frames")
    if os.path.exists(tmp_dir):
      shutil.rmtree(tmp_dir)
    os.makedirs(tmp_dir, exist_ok=True)
    realized = []
    for i, src in enumerate(src_files):
      dst = os.path.join(tmp_dir, f"frame_{i:06d}.png")
      try:
        # A hard link is fast and space-efficient if on the same filesystem
        os.link(src, dst)
      except Exception:
        try:
          # Use a symbolic link if a hard link is not possible
          os.symlink(os.path.abspath(src), dst)
        except Exception:
          # Copy in environments where that is also not possible
          shutil.copyfile(src, dst)
      realized.append(dst)
    # Returns the "sequentially numbered files sorted in playback order" and the temp directory path
    return realized, tmp_dir

  if "all" in file_lists:
    frames = list(file_lists["all"])  # order was already fixed when file_lists was built (handles descending order)
    if _should_reverse_from_steps(run_params, cross_section):
      frames = list(reversed(frames))
    _debug_dump("gifski frames (all, original)", frames)
    seq_frames, tmp_dir = materialize_sequenced_frames("all", frames)
    _debug_dump("gifski frames (all, sequenced)", seq_frames)
    output_gif = os.path.join(out_dir, f"{render_params['output_png_prefix']}_{cross_section}_anim.gif")
    command = build_command(output_gif, seq_frames)
    print("Running:", " ".join(command))
    subprocess.run(command, check=True)
    print(f"Created animated gif: {output_gif}")
    # Cleanup temporary directory
    if tmp_dir and os.path.exists(tmp_dir):
      shutil.rmtree(tmp_dir)
  if run_params.get("tf_out_each_det", False):
    for key, files in file_lists.items():
      if key == "all":
        continue
      frames = list(files)
      if _should_reverse_from_steps(run_params, cross_section):
        frames = list(reversed(frames))
      _debug_dump(f"gifski frames (det {key}, original)", frames)
      seq_frames, tmp_dir = materialize_sequenced_frames(f"det{key:03d}", frames)
      _debug_dump(f"gifski frames (det {key}, sequenced)", seq_frames)
      output_gif = os.path.join(out_dir, f"{render_params['output_png_prefix']}_{cross_section}_det{key:03d}_anim.gif")
      command = build_command(output_gif, seq_frames)
      print("Running:", " ".join(command))
      subprocess.run(command, check=True)
      print(f"Created animated gif: {output_gif}")
      # Cleanup temporary directory
      if tmp_dir and os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)

# Run convert (PDF creation) for cross_section
def run_convert_pdf_cross(file_lists, render_params, run_params, cross_section):
  out_dir = render_params.get('save_png_dir', render_params['output_png_prefix'])
  convert_options = render_params.get('convert_options', [])
  pdf_processes = []
  reverse_order = _should_reverse_from_steps(run_params, cross_section)
  for key, files in file_lists.items():
    files_use = list(reversed(files)) if reverse_order else files
    if key == "all":
      output_pdf = os.path.join(out_dir, f"{render_params['output_png_prefix']}_{cross_section}_anim.pdf")
    else:
      output_pdf = os.path.join(out_dir, f"{render_params['output_png_prefix']}_{cross_section}_det{key:03d}_anim.pdf")
    command = ["convert"] + convert_options + files_use + [output_pdf]
    print("Running convert command:", " ".join(command))
    proc = subprocess.Popen(command)
    pdf_processes.append(proc)
  for proc in pdf_processes:
    proc.wait()
  print(f"All PDF convert processes finished for {cross_section}.")

def build_bin_edges_from_header(run_params, header_info):
  """Build bin edges from .tmp file header (authoritative grid info from C++).

  Header min/max are CELL EDGES (not centers).
  Grid1d: center(i) = min + (i + 0.5) * interval
  So bin edges are simply: min, min+step, min+2*step, ..., max

  Stores {x,y,z}_bin_edges in run_params without overwriting
  user-specified display range (xmin/xmax etc.).
  """
  for axis in ("x", "y", "z"):
    info = header_info.get(f"{axis}_info")
    if not info or len(info) < 3:
      continue
    hmin, hmax, hstep = map(float, info[:3])
    n_bins = int(round((hmax - hmin) / hstep))
    edges = np.linspace(hmin, hmax, n_bins + 1)
    run_params[f"{axis}_bin_edges"] = edges
    run_params[f"header_{axis}min"] = hmin
    run_params[f"header_{axis}max"] = hmax
    run_params[f"header_{axis}step"] = hstep
  return run_params

def apply_header_bins(run_params, header_info):
  """Apply bin settings from file header (x_info, y_info, z_info)."""
  if not run_params.get("use_header_bins", False):
    return run_params
  updated = dict(run_params)
  axes = {
    "x": header_info.get("x_info"),
    "y": header_info.get("y_info"),
    "z": header_info.get("z_info"),
  }
  for axis, info in axes.items():
    if not info or len(info) < 3:
      continue
    try:
      vmin, vmax, vstep = map(float, info[:3])
    except (TypeError, ValueError):
      continue
    updated[f"{axis}min"] = vmin
    updated[f"{axis}max"] = vmax
    updated[f"{axis}step"] = vstep
    print(f"[INFO] use_header_bins: {axis}min/max/step -> {vmin}, {vmax}, {vstep}", file=sys.stderr)
  return updated


def apply_data_bins(run_params, df):
  """Apply bin settings computed from actual data coordinates.

  For non-uniform grids, this computes bin edges as midpoints between
  consecutive unique values, ensuring each data point falls into exactly one bin.
  """
  if not run_params.get("use_data_bins", False):
    return run_params
  updated = dict(run_params)

  for axis, col in [("x", "xcnt"), ("y", "ycnt"), ("z", "zcnt")]:
    unique_vals = np.array(sorted(df[col].unique()))
    if len(unique_vals) < 2:
      continue

    # Compute bin edges as midpoints between consecutive unique values
    # This handles non-uniform grids correctly
    midpoints = (unique_vals[:-1] + unique_vals[1:]) / 2
    # Add edges at the extremes (half the nearest gap)
    first_gap = unique_vals[1] - unique_vals[0]
    last_gap = unique_vals[-1] - unique_vals[-2]
    bin_edges = np.concatenate([
      [unique_vals[0] - first_gap / 2],
      midpoints,
      [unique_vals[-1] + last_gap / 2]
    ])

    # Store bin edges directly for use in plotting
    updated[f"{axis}_bin_edges"] = bin_edges

    # Also compute approximate min/max/step for compatibility
    diffs = np.diff(unique_vals)
    vstep = float(np.median(diffs))
    vmin = float(bin_edges[0])
    vmax = float(bin_edges[-1])
    updated[f"{axis}min"] = vmin
    updated[f"{axis}max"] = vmax
    updated[f"{axis}step"] = vstep
    print(f"[INFO] use_data_bins: {axis} n_bins={len(bin_edges)-1}, range=[{vmin:.2f}, {vmax:.2f}]", file=sys.stderr)

  return updated

# Extract detector x,y coordinates and name from det_info_path
def get_detector_xy(det_info_paths):
  if json5 is None:
    print("[ERROR] json5 module not found. Please install it with `pip install json5` or similar.", file=sys.stderr)
    sys.exit(1)

  det_xy = []
  for path in det_info_paths:
    try:
      with open(path, 'r') as f:
        detector_data = json5.load(f)
    except Exception as e:
      print(f"[ERROR] Failed to read {path}: {e}", file=sys.stderr)
      sys.exit(1)

    params = detector_data.get("DETECTOR_PARAMETERS")
    if not isinstance(params, dict):
      print(f"[ERROR] {path} does not have a DETECTOR_PARAMETERS object.", file=sys.stderr)
      sys.exit(1)

    def extract_axis(axis):
      if axis not in params:
        print(f"[ERROR] {path}'s DETECTOR_PARAMETERS does not have \"{axis}\".", file=sys.stderr)
        sys.exit(1)
      value = params[axis]
      if isinstance(value, (int, float)):
        return float(value)
      print(f"[ERROR] {path}'s \"{axis}\" must be specified as a number.", file=sys.stderr)
      sys.exit(1)

    x_val = extract_axis("x")
    y_val = extract_axis("y")
    name_val = params.get("name", "")
    det_xy.append((x_val, y_val, name_val))

  return det_xy

# In main(), create each cross-section plot according to run_params' "cross_section" value
def _cleanup_directories(render_params, keep_png=False):
  """Cleanup intermediate files and move outputs to parent directory.

  This function:
  1. Moves *.gif and *.pdf from subdirectories to base_dir (figs/)
  2. Removes *.png intermediate files in subdirectories (unless keep_png=True)
  3. Removes file_list_*.tmp temporary files
  4. Removes subdirectories after cleanup (unless non-empty due to kept PNGs)
  """
  out_dir = render_params.get('save_png_dir', '')
  base_dir = render_params.get('save_png_dir_base', '')

  # Remove file_list_*.tmp in output directory
  if out_dir and os.path.isdir(out_dir):
    for tmp_file in glob.glob(os.path.join(out_dir, "file_list_*.tmp")):
      try:
        os.remove(tmp_file)
      except Exception:
        pass

  # Process subdirectories in base_dir: move outputs, delete intermediates
  if base_dir and os.path.isdir(base_dir):
    # Remove file_list_*.tmp in base_dir
    for tmp_file in glob.glob(os.path.join(base_dir, "file_list_*.tmp")):
      try:
        os.remove(tmp_file)
      except Exception:
        pass

    # Process each subdirectory
    for subdir in glob.glob(os.path.join(base_dir, "*")):
      if not os.path.isdir(subdir):
        continue

      # Remove file_list_*.tmp in subdirectory
      for tmp_file in glob.glob(os.path.join(subdir, "file_list_*.tmp")):
        try:
          os.remove(tmp_file)
        except Exception:
          pass

      # Move *.gif and *.pdf to base_dir
      for ext in ['*.gif', '*.pdf']:
        for src_file in glob.glob(os.path.join(subdir, ext)):
          dst_file = os.path.join(base_dir, os.path.basename(src_file))
          try:
            shutil.move(src_file, dst_file)
            print(f"[INFO] Moved: {os.path.basename(src_file)} -> {base_dir}")
          except Exception as e:
            print(f"[WARN] Failed to move {src_file}: {e}")

      # Remove *.png intermediate files (unless --keep-png)
      if not keep_png:
        for png_file in glob.glob(os.path.join(subdir, "*.png")):
          try:
            os.remove(png_file)
          except Exception:
            pass

      # Remove subdirectory if empty
      if os.path.isdir(subdir) and not os.listdir(subdir):
        try:
          os.rmdir(subdir)
          print(f"[INFO] Removed subdirectory: {subdir}")
        except Exception:
          pass

  # Remove output directory if empty (only if different from base_dir)
  if out_dir and out_dir != base_dir and os.path.isdir(out_dir) and not os.listdir(out_dir):
    try:
      os.rmdir(out_dir)
      print(f"[INFO] Removed empty output directory: {out_dir}")
    except Exception:
      pass


def main():
  parser = argparse.ArgumentParser(description='Load configuration and create cross-sectional plots')
  parser.add_argument('--config', type=str, default='config.json', help='Path to the JSON configuration file')
  bin_group = parser.add_mutually_exclusive_group()
  bin_group.add_argument('--use-header-bins', action='store_true',
                         help='Use xmin/xmax/xstep etc. from the input file header (x_info, y_info, z_info)')
  bin_group.add_argument('--use-data-bins', action='store_true',
                         help='Compute bin edges from actual data coordinates (recommended for accurate binning)')
  # Keep --bin_as_is as alias for --use-header-bins for backward compatibility
  bin_group.add_argument('--bin_as_is', action='store_true',
                         help=argparse.SUPPRESS)  # Hidden, deprecated alias
  parser.add_argument('--keep-png', '--keep_png', action='store_true',
                      default=True,
                      help='Keep intermediate PNG files (default: True)')
  parser.add_argument('--delete-png', '--delete_png', action='store_true',
                      help='Delete intermediate PNG files after GIF/PDF creation')
  args = parser.parse_args()

  config = load_config(args.config)
  run_params = config.get("run_params", {})
  render_params = config.get("render_params", {})

  # Handle bin mode flags
  if args.use_header_bins or args.bin_as_is:
    run_params["use_header_bins"] = True
    print("[INFO] --use-header-bins: will use xmin/xmax/xstep from input file header")
  elif args.use_data_bins:
    run_params["use_data_bins"] = True
    print("[INFO] --use-data-bins: will compute bin edges from actual data coordinates")

  print("Run Parameters:")
  for key, value in run_params.items():
    print(f"  {key}: {value}")
  print("\nRender Parameters:")
  for key, value in render_params.items():
    print(f"  {key}: {value}")

  _resolve_output_directory(args.config, run_params, render_params)

  try:
    # Get cross_section (defaults to "z" if not specified)
    cross_sections = run_params.get("cross_section", "z")
    if not isinstance(cross_sections, list):
      cross_sections = [cross_sections]

    # * Input leniency: also normalize the case where cross_section is mistakenly passed as the string '["-z"]'
    norm = []
    for cs in cross_sections:
      if isinstance(cs, str):
        s = cs.strip()
        if s.startswith('[') and s.endswith(']'):
          try:
            parsed = json.loads(s)
            norm.extend(parsed if isinstance(parsed, list) else [parsed])
            continue
          except Exception:
            pass
      norm.append(cs)
    cross_sections = norm

    # If det_info_path is present, get the detector positions
    det_info_paths = config.get("det_info_path", [])
    det_xy = None
    if det_info_paths:
      det_xy = get_detector_xy(det_info_paths)
      print("Detector positions:")
      for i, det_info in enumerate(det_xy):
        x, y = det_info[0], det_info[1]
        name = det_info[2] if len(det_info) > 2 else ""
        print(f"  det_id {i}: x={x}, y={y}, name={name}")

    filename_in = run_params.get("filename_in")
    if not filename_in:
      print("[ERROR] run_params.filename_in is not specified in config", file=sys.stderr)
      sys.exit(1)
    if not os.path.exists(filename_in):
      print(f"[ERROR] Input file not found: {filename_in}", file=sys.stderr)
      sys.exit(1)

    header_info, data = read_cross_section_file(filename_in)
    detected_format = header_info.get("format", ["text"])[0]
    print(f"[INFO] Detected cross-section format (verified): {detected_format}")
    n_detector = int(header_info.get("n_detector", [0])[0])
    # if len(det_info_paths) != n_detector:
    #   raise ValueError(f"Number of detector info files ({len(det_info_paths)}) does not match n_detector ({n_detector})")

    if isinstance(data, dict):
      # Binary reader returns dict of arrays
      df = pd.DataFrame(data)
    else:
      # Text reader returns list-of-lists
      columns = ['tf_exist', 'xcnt', 'ycnt', 'zcnt', 'density', 'n_det'] + \
                [f'flag_{i}' for i in range(n_detector)] + ['n_hit_det']
      df = pd.DataFrame(data, columns=columns)

    print("Header Information:")
    for key, value in header_info.items():
      print(f"{key}: {value}")
    print("\nData:")
    print(df)

    # Always build bin edges from header (authoritative grid info)
    run_params = build_bin_edges_from_header(run_params, header_info)

    # "header" / "data" flags: legacy behavior (also overwrites display range)
    if run_params.get("use_header_bins", False):
      run_params = apply_header_bins(run_params, header_info)
    elif run_params.get("use_data_bins", False):
      run_params = apply_data_bins(run_params, df)
    df_filtered = filter_plot_data(df, run_params)
    print(f"\nFiltered data has {len(df_filtered)} rows.")

    # Get output formats (default: png and gif; pdf is optional)
    # Supported formats: "png", "gif", "pdf"
    output_formats = render_params.get("output_formats", ["png", "gif"])
    if isinstance(output_formats, str):
      output_formats = [output_formats]
    output_formats = [fmt.lower() for fmt in output_formats]
    print(f"\nOutput formats: {output_formats}")

    # Run processing for each cross_section
    for cs in cross_sections:
      print(f"\nProcessing cross_section: {cs}")
      # PNG is always generated (required for gif/pdf)
      unique_vals = plot_and_save_all_cross_sections(df, run_params, render_params, n_detector, det_xy, cs)
      file_lists_cs = create_filename_list_cross(df, run_params, render_params, n_detector, cs, unique_vals=unique_vals)
      save_filename_list_to_file_cross(file_lists_cs, render_params, cs)
      # Generate GIF if requested
      if "gif" in output_formats:
        run_gifski_cross(file_lists_cs, render_params, run_params, cs)
      # Generate PDF if requested
      if "pdf" in output_formats:
        run_convert_pdf_cross(file_lists_cs, render_params, run_params, cs)

  finally:
    # Always cleanup, even on error
    keep_png = not args.delete_png
    _cleanup_directories(render_params, keep_png=keep_png)

if __name__ == "__main__":
  main()
