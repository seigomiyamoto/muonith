#!/usr/bin/env python
import json
import re
import sys
import os
import datetime
import copy
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib as mpl
import numpy as np
try:
  import pandas as pd
except ImportError:
  sys.exit(
    "ERROR: pandas is required.\n"
    "  pip install pandas  (or activate venv: source .venv/bin/activate)"
  )
import subprocess
import glob
from concurrent.futures import ProcessPoolExecutor

def load_data(file_name):
  """
  Function that reads data from a file and returns metadata and the data body.
  """
  with open(file_name, 'r') as f:
    lines = f.readlines()

  meta = {}
  tokens = lines[0].split()
  if len(tokens) < 3:
    raise ValueError("Line 1 does not contain enough information")
  meta["xmin"] = float(tokens[1])
  meta["xmax"] = float(tokens[2])

  tokens = lines[1].split()
  if len(tokens) < 3:
    raise ValueError("Line 2 does not contain enough information")
  meta["ymin"] = float(tokens[1])
  meta["ymax"] = float(tokens[2])

  data = pd.read_csv(
    file_name,
    sep=r'\s+',
    skiprows=2,
    header=None,
    names=["ig", "xlow", "xup", "ylow", "yup", "sig", "noi", "signoi"
          , "dens_low", "dens_cnt", "dens_upp"
          , "delta_nmuon_low", "delta_nmuon_cnt", "delta_nmuon_upp"
          , "volume", "eff_low", "eff_cnt", "eff_upp", "is_avail"]
  )

  data["dens_err"] = (data["dens_upp"] - data["dens_low"]) * 0.5
  data["eff_err"] = (data["eff_upp"] - data["eff_low"]) * 0.5

  # signal-to-noise ratio (sig / noi); NaN where noise is 0 so the cell
  # falls back to the background color (no noise -> ratio undefined)
  data["sig_over_noi"] = data["sig"] / data["noi"].where(data["noi"] > 0.0)

  return meta, data

def plot_field(field, output_prefix, str_det_id, params, meta, data, png_dpi
            , output_dir, datafile, data0=None, rect_line_width=0.2, rect_edge_color="black"
            , bg_color="white", grid_color="grey", gridline_type="dotted", gridline_width=1.0):
  """
  Function that generates a plot for a single field.
  """
  n_grad  = params.get("ngrad")
  vmin_in = params.get("vmin")
  vmax_in = params.get("vmax")
  alpha_in = params.get("alpha")

  if n_grad is None or vmin_in is None or vmax_in is None or alpha_in is None:
    print(f"Error: Field '{field}' is missing required parameters.")
    return

  xmin = meta["xmin"]
  xmax = meta["xmax"]
  ymin = meta["ymin"]
  ymax = meta["ymax"]

  fig, ax = plt.subplots()
  ax.set_xlim(xmin, xmax)
  ax.set_ylim(ymin, ymax)
  ax.set_aspect('equal')
  ax.set_facecolor(bg_color)

  hori_inch = 10
  aspect_ratio = (xmax - xmin) / (ymax - ymin)
  vert_inch = hori_inch / aspect_ratio
  fig.set_size_inches(hori_inch, vert_inch)

  # Get the palette specified in the JSON and generate the colormap
  pallet = params.get("pallet", "jet")
  cmap = copy.copy(plt.get_cmap(pallet))
  color_under = params.get("color_under", 'silver')
  color_over = params.get("color_over", 'black')
  cmap.set_under(color_under)
  cmap.set_over(color_over)

  if params.get("logz", False):
    if vmin_in <= 0:
      raise ValueError("Log scale requires vmin > 0")
    norm = mpl.colors.LogNorm(vmin=vmin_in, vmax=vmax_in)
  else:
    bounds = np.linspace(vmin_in, vmax_in, n_grad + 1)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='both')

  n_data = 0
  for index, row in data.iterrows():
    n_data += 1
    # print(f"[DEBUG] ig={row['ig']}, is_avail={row['is_avail']}, field={field}, value={row.get(field, 'N/A')}")
    if row["is_avail"] != 1:
      print(f"ig={row['ig']} is not available.")
      continue

    if field == "signoi":
      value = row["signoi"]
    elif field == "sig":
      value = row["sig"]
    elif field == "noi":
      value = row["noi"]
    elif field == "sig_over_noi":
      value = row["sig_over_noi"]
    elif field == "dens_low":
      value = row["dens_low"]
    elif field == "dens_cnt":
      value = row["dens_cnt"]
    elif field == "dens_upp":
      value = row["dens_upp"]
    elif field == "dens_err":
      value = row["dens_err"]
    elif field == "delta_nmuon_low":
      value = row["delta_nmuon_low"]
    elif field == "delta_nmuon_cnt":
      value = row["delta_nmuon_cnt"]
    elif field == "delta_nmuon_upp":
      value = row["delta_nmuon_upp"]
    elif field == "dens_cnt_diff":
      if data0 is None:
        print("Warning: 'dens_cnt_diff' requested but no datafile0 provided.")
        return
      try:
        row0 = data0.iloc[index]
      except IndexError:
        print(f"Error: data index {index} out of range in datafile0.")
        continue
      value = row["dens_cnt"] - row0["dens_cnt"]
    elif field == "volume":
      value = row["volume"]
    elif field == "eff_low":
      value = row["eff_low"]
    elif field == "eff_cnt":
      value = row["eff_cnt"]
    elif field == "eff_upp":
      value = row["eff_upp"]
    elif field == "eff_err":
      value = row["eff_err"]
    else:
      if field in row:
        value = row[field]
      else:
        print(f"Unknown field: {field}. Skipping row ig={row['ig']}.")
        continue

    # for non-finite values (e.g. sig_over_noi where noise == 0) keep the
    # cell outline but fill it with the background color (no data), so the
    # rectangle grid stays visible instead of vanishing into the background
    is_no_data = not np.isfinite(value)

    # For log scale, replace the value with 1.0e-99 if it is 0 or below
    if not is_no_data and params.get("logz", False) and value <= 0.0:
      value = 1.0e-99
    
    xlow = row["xlow"]
    xup  = row["xup"]
    ylow = row["ylow"]
    yup  = row["yup"]
    xlen = xup - xlow
    ylen = yup - ylow

    rect = patches.Rectangle(
      (xlow, ylow),
      xlen,
      ylen,
      linewidth=rect_line_width,
      edgecolor=rect_edge_color,
      facecolor=bg_color if is_no_data else cmap(norm(value)),
      alpha=alpha_in
    )
    ax.add_patch(rect)

    if n_data % 100 == 0:
      sys.stdout.write(f"Processing n_data={n_data}\r")
      sys.stdout.flush()

  # print(f"[DEBUG] Total data rows passed to plot_field: {len(data)}")

  x_ticks = ax.get_xticks()
  y_ticks = ax.get_yticks()
  ax.vlines(x_ticks, ymin=ymin, ymax=ymax, colors=grid_color, linestyles=gridline_type, linewidth=gridline_width)
  ax.hlines(y_ticks, xmin=xmin, xmax=xmax, colors=grid_color, linestyles=gridline_type, linewidth=gridline_width)

  sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
  sm.set_array([])
  cbar = plt.colorbar(sm, ax=ax, extend='both')
  cbar.set_label(field)

  now = datetime.datetime.now()
  timestamp_str = now.strftime("%Y-%m-%d %H:%M:%S")
  label_str = f"{datafile} ({timestamp_str})"
  ax.text(0.5, 1.02, label_str, transform=ax.transAxes, fontsize=12, ha='center', va='bottom')

  output_file_name = os.path.join(output_dir, f'{output_prefix}det{str_det_id}_{field}.png')
  plt.savefig(output_file_name, dpi=png_dpi)
  plt.close(fig)
  print(f"Saved {output_file_name}")

  # If dump_txt is True, dump the field values to a TXT file
  if params.get("dump_txt", False):
    dump_field_values(data, field, str_det_id, output_dir, output_prefix)

def extract_det_number(basename):
  """
  Function that extracts and returns the number between '_det' and the next '_' from basename.
  """
  m = re.search(r'_det(\d+)_', basename)
  if m:
    return m.group(1)
  return None

def parse_arguments():
  """
  Function that parses command-line arguments in --key=value format.
  """
  args = {}
  for arg in sys.argv[1:]:
    if arg.startswith("--"):
      if "=" in arg:
        key, value = arg[2:].split("=", 1)
        args[key] = value
      else:
        args[arg[2:]] = None
  return args

def dump_field_values(data, field, str_det_id, output_dir, output_prefix):
  output_txt_name = os.path.join(output_dir, f"{output_prefix}det{str_det_id}_{field}-tmp.txt")
  with open(output_txt_name, "w") as f:
    f.write(f"# xlow\txup\tylow\tyup\t{field}\n")
    for index, row in data.iterrows():
      if field not in row or pd.isnull(row[field]):
        continue
      f.write(f"{row['xlow']:.6E}\t{row['xup']:.6E}\t{row['ylow']:.6E}\t{row['yup']:.6E}\t{row[field]:.6E}\n")
  print(f"Saved TXT: {output_txt_name}")


def process_task(args):
  """
  Individual task processing for each (det_id, field)
  """
  (str_det_id, field, params,
   input_dir, input_prefix, input_suffix,
   ref_input_dir, png_dpi, output_prefix, output_dir,
   rect_line_width, rect_edge_color,
   bg_color, grid_color, gridline_type, gridline_width) = args
  
  datafile = os.path.join(input_dir, f"{input_prefix}det{str_det_id}{input_suffix}")
  if not os.path.exists(datafile):
    print(f"Error: Datafile '{datafile}' does not exist.")
    return

  print(f"Processing {datafile} for field '{field}'...", flush=True)

  # Virtual field processing (e.g., signoi_diff_ratio)
  if params.get("is_virtual", False) and field == "signoi_diff_ratio":
    input_dir = params.get("input_dir", ".")
    fileA_prefix = params.get("fileA_prefix")
    fileA_suffix = params.get("fileA_suffix")
    fileB_prefix = params.get("fileB_prefix")
    fileB_suffix = params.get("fileB_suffix")

    fileA = os.path.join(input_dir, f"{fileA_prefix}det{str_det_id}{fileA_suffix}")
    fileB = os.path.join(input_dir, f"{fileB_prefix}det{str_det_id}{fileB_suffix}")

    if not os.path.exists(fileA) or not os.path.exists(fileB):
      print(f"Missing virtual field inputs: {fileA} or {fileB}")
      return

    metaA, dataA = load_data(fileA)
    metaB, dataB = load_data(fileB)

    if metaA != metaB:
      print(f"Meta mismatch between {fileA} and {fileB}")
      return

    if "signoi" not in dataA.columns or "signoi" not in dataB.columns:
      print(f"'signoi' not found in either data.")
      return

    data = dataA.copy()
    epsilon = 1.0e-12
    data["signoi_diff_ratio"] = (dataA["signoi"] - dataB["signoi"]) / (dataA["signoi"] + dataB["signoi"] + epsilon)

    plot_field("signoi_diff_ratio", output_prefix, str_det_id, params
               , metaA, data, png_dpi, output_dir, f"VIRTUAL_diff_det{str_det_id}"
               , rect_line_width=rect_line_width, rect_edge_color=rect_edge_color
               , bg_color=bg_color, grid_color=grid_color
               , gridline_type=gridline_type, gridline_width=gridline_width)
    return

  if ref_input_dir:
    datafile0 = os.path.join(ref_input_dir, f"{input_prefix}det{str_det_id}{input_suffix}")
    if os.path.exists(datafile0):
      meta0, data0 = load_data(datafile0)
      meta, data = load_data(datafile)
      if meta0 != meta:
        print("Error: Meta information differs between datafile0 and datafile.", flush=True)
        print("Falling back to normal plot without reference.", flush=True)
        data0 = None  # diff not possible
      else:
        print(f"xmin={meta['xmin']}, xmax={meta['xmax']}, ymin={meta['ymin']}, ymax={meta['ymax']}", flush=True)
        plot_field(field, output_prefix, str_det_id, params, meta, data, png_dpi
                   , output_dir, datafile, data0=data0, rect_line_width=rect_line_width, rect_edge_color=rect_edge_color
                   , bg_color=bg_color, grid_color=grid_color
                   , gridline_type=gridline_type, gridline_width=gridline_width)
        return  # exit because the reference plot is complete
    else:
      print(f"Warning: ref file '{datafile0}' not found. Falling back to single data plot.", flush=True)

  # Normal processing (common processing when ref is unavailable)
  meta, data = load_data(datafile)
  print(f"xmin={meta['xmin']}, xmax={meta['xmax']}, ymin={meta['ymin']}, ymax={meta['ymax']}", flush=True)
  plot_field(field, output_prefix, str_det_id, params, meta, data, png_dpi, output_dir, datafile
             , rect_line_width=rect_line_width, rect_edge_color=rect_edge_color
             , bg_color=bg_color, grid_color=grid_color
             , gridline_type=gridline_type, gridline_width=gridline_width)

def process_pdf_gif_per_detid(det_id, output_dir
  , output_prefix, pdf_prefix, pdf_suffix, pdf_dpi
  , gif_prefix, gif_suffix, gif_dpi, gif_delay, gif_loop
  , exec_erace_pngs, make_pdf, make_gif
  ):
  """
  For each det_id, creates a PDF and GIF from the output PNGs, and
  (if exec_erace_pngs is True) deletes the PNG files after creation.
  """
  print(f"Processing PDF and GIF for det{det_id}...", flush=True)

  # Example PNG file pattern: files starting with output_prefix + "det" + det_id + "_"
  pattern = os.path.join(output_dir, f"{output_prefix}det{det_id}_*.png")
  print(f"[DEBUG] glob pattern: {pattern}", flush=True)
  png_files = sorted(glob.glob(pattern))
  print(f"[DEBUG] found PNG files: {png_files}", flush=True)
  if not png_files:
    print(f"det_id {det_id}: No PNG files found.", flush=True)
    return

  # Create PDF (using ImageMagick convert)
  if make_pdf:
    pdf_filename = os.path.join(output_dir, f"{pdf_prefix}det{det_id}{pdf_suffix}.pdf")
    cmd_pdf = ["convert", "-density", str(pdf_dpi)] + png_files + [pdf_filename]
    print(f"Creating PDF for det{det_id}: {' '.join(cmd_pdf)}", flush=True)
    subprocess.run(cmd_pdf, check=True)

  # Create animated GIF (using ImageMagick convert)
  if make_gif:
    gif_filename = os.path.join(output_dir, f"{gif_prefix}det{det_id}{gif_suffix}.gif")
    # -delay, -loop must be specified before the input images
    cmd_gif = ["convert", "-density", str(gif_dpi), "-delay", str(gif_delay), "-loop", str(gif_loop)] + png_files + [gif_filename]
    print(f"Creating GIF for det{det_id}: {' '.join(cmd_gif)}", flush=True)
    subprocess.run(cmd_gif, check=True)

  # If exec_erace_pngs is True, delete PNG files after creating PDF and GIF
  if exec_erace_pngs:
    for png in png_files:
      try:
        os.remove(png)
      except Exception as e:
        print(f"Error removing {png}: {e}", flush=True)
    print(f"Removed {len(png_files)} PNG files for det{det_id}.", flush=True)
  else:
    print(f"PNG file removal skipped for det{det_id}.", flush=True)

def process_pdf_gif_per_field(field, det_ids, output_dir, output_prefix
  , pdf_prefix, pdf_suffix, pdf_dpi
  , gif_prefix, gif_suffix, gif_dpi
  , gif_delay, gif_loop, exec_erace_pngs
  , make_pdf, make_gif):
  """
  Function that gathers PNGs for the same field across all det_id (list)
  and generates a single PDF/GIF.
  """
  # First, loop over the list of det_ids to collect PNG files
  png_files = []
  for det_id in det_ids:
    fn = os.path.join(output_dir, f"{output_prefix}det{det_id}_{field}.png")
    if os.path.exists(fn):
      png_files.append(fn)
    else:
      print(f"[GroupByField] Warning: {fn} not found. Skipping.", flush=True)

  if not png_files:
    print(f"[GroupByField] No PNGs found for field='{field}'.", flush=True)
    return

  png_files = sorted(png_files)
  print(f"[GroupByField] Target PNGs: {png_files}", flush=True)

  # Create PDF
  if make_pdf:
    pdf_filename = os.path.join(output_dir, f"{pdf_prefix}{field}{pdf_suffix}.pdf")
    cmd_pdf = ["convert", "-density", str(pdf_dpi)] + png_files + [pdf_filename]
    print(f"[GroupByField] Creating PDF for field '{field}': {' '.join(cmd_pdf)}", flush=True)
    subprocess.run(cmd_pdf, check=True)

  # Create GIF
  if make_gif:
    gif_filename = os.path.join(output_dir, f"{gif_prefix}{field}{gif_suffix}.gif")
    cmd_gif = ["convert", "-density", str(gif_dpi), "-delay", str(gif_delay), "-loop", str(gif_loop)] + png_files + [gif_filename]
    print(f"[GroupByField] Creating GIF for field '{field}': {' '.join(cmd_gif)}", flush=True)
    subprocess.run(cmd_gif, check=True)

  # Delete PNGs
  if exec_erace_pngs:
    for f in png_files:
      try:
        os.remove(f)
      except Exception as e:
        print(f"[GroupByField] Error removing {f}: {e}", flush=True)
    print(f"[GroupByField] Removed {len(png_files)} PNG files for field '{field}'", flush=True)


# Move run_pdf_gif_per_detid to top level
def run_pdf_gif_per_detid(args):
  process_pdf_gif_per_detid(*args)

# Move run_pdf_gif_per_field to top level
def run_pdf_gif_per_field(args):
  process_pdf_gif_per_field(*args)


def main():
  args = parse_arguments()
  if "jsonfile" not in args:
    print("Error: Missing required parameters. --jsonfile must be specified.")
    sys.exit(1)

  jsonfile = args["jsonfile"]
  with open(jsonfile, 'r') as f:
    config = json.load(f)

  input_dir      = config.get("input_dir", "tmp")
  input_prefix   = config.get("input_prefix", "")
  input_suffix  = config.get("input_suffix", ".tmp")
  output_dir     = config.get("output_dir", "figs")
  output_prefix  = config.get("output_prefix", "")
  set_str_det_id = config.get("det_ids", [])

  png_dpi         = config.get("png_dpi")

  make_pdf_per_detid = config.get("make_pdf_per_detid", False)
  make_pdf_per_field = config.get("make_pdf_per_field", False)
  pdf_prefix      = config.get("pdf_prefix", "")
  pdf_suffix     = config.get("pdf_suffix", "")
  pdf_dpi         = config.get("pdf_dpi", 100)

  make_gif_per_detid = config.get("make_gif_per_detid", False)
  make_gif_per_field = config.get("make_gif_per_field", False)
  gif_prefix      = config.get("gif_prefix", "")
  gif_suffix     = config.get("gif_suffix", "")
  gif_delay       = config.get("gif_delay", 200)
  gif_loop        = config.get("gif_loop", 0)
  gif_dpi         = config.get("gif_dpi", 100)

  ref_input_dir   = config.get("ref_input_dir", "")
  plot_fields     = config.get("plot_fields", {})
  exec_erace_pngs = config.get("exec_erace_pngs", False)
  rect_line_width = config.get("rect_line_width", 0.2)
  rect_edge_color = config.get("rect_edge_color", "black")
  bg_color        = config.get("bg_color", "white")
  grid_color      = config.get("grid_color", "grey")
  gridline_type   = config.get("gridline_type", "dotted")
  gridline_width  = config.get("gridline_width", 1.0)

  if png_dpi is None:
    print("Error: Missing dpi required parameters in JSON configuration.")
    sys.exit(1)

  if not plot_fields:
    print("Error: Missing plot fields in JSON")
    sys.exit(1)

  if not os.path.exists(input_dir):
    print(f"Error: Input directory '{input_dir}' does not exist.")
    sys.exit(1)

  if not os.path.exists(output_dir):
    os.makedirs(output_dir)

  # 1. Parallel execution of PNG creation (det_id, field combinations)
  tasks = []
  for str_det_id in set_str_det_id:
    for field, params in plot_fields.items():
      if not params.get("exec", False):
        print(f"[SKIP] '{field}' is disabled (exec=False).")
        continue
      tasks.append((str_det_id, field, params,
                    input_dir, input_prefix, input_suffix,
                    ref_input_dir, png_dpi, output_prefix, output_dir
                    , rect_line_width, rect_edge_color
                    , bg_color, grid_color, gridline_type, gridline_width))

  with ProcessPoolExecutor() as executor:
    # Consume the iterator: a discarded map() never retrieves the futures, so a
    # worker exception stays buried and the run exits 0 with no figures written.
    list(executor.map(process_task, tasks))

  # A. PDF/GIF per det_id
  pdf_gif_tasks = []
  for det_id in set_str_det_id:
    pdf_gif_tasks.append((det_id, output_dir, output_prefix, pdf_prefix, pdf_suffix, pdf_dpi
      , gif_prefix, gif_suffix, gif_dpi, gif_delay, gif_loop
      , exec_erace_pngs, make_pdf_per_detid, make_gif_per_detid))

  with ProcessPoolExecutor() as executor:
    executor.map(run_pdf_gif_per_detid, pdf_gif_tasks)

  # B. PDF/GIF per field
  pdf_gif_field_tasks = []
  if make_pdf_per_field or make_gif_per_field:
    for field, params in plot_fields.items():
      if not params.get("exec", False):
        continue
      pdf_gif_field_tasks.append((
        field
        , set_str_det_id   # <- pass the list of all det_ids
        , output_dir, output_prefix
        , pdf_prefix, pdf_suffix, pdf_dpi
        , gif_prefix, gif_suffix, gif_dpi, gif_delay, gif_loop
        , exec_erace_pngs
        , make_pdf_per_field, make_gif_per_field
      ))

  with ProcessPoolExecutor() as executor:
    executor.map(run_pdf_gif_per_field, pdf_gif_field_tasks)

if __name__ == "__main__":
  main()
