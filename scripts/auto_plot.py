#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Auto-plot script for muonith outputs.
Detects output files in current directory and runs appropriate plotting scripts.

Usage:
  python3 auto_plot.py [--config config.json5] [--dry-run] [--exclude PATTERN]

Examples:
  # Basic usage (after exe execution)
  python3 ../../scripts/auto_plot.py

  # Dry-run to see what would be executed
  python3 ../../scripts/auto_plot.py --dry-run

  # Skip certain plot types
  python3 ../../scripts/auto_plot.py --skip-hist2d

  # Exclude specific file patterns
  python3 ../../scripts/auto_plot.py --exclude "arrdet_g3vox_prior_*"
  python3 ../../scripts/auto_plot.py -e "*_debug_*" -e "*_test_*"

Configuration:
  Uses JSON5 format (supports comments). See auto_plot_default.json5 for defaults.
  Exclude patterns can be set globally or per-section in the config file.
"""

import argparse
import glob
import json
import os
import shlex
import subprocess
import sys
from fnmatch import fnmatch
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
  import json5
  HAS_JSON5 = True
except ImportError:
  HAS_JSON5 = False
  import warnings
  warnings.warn(
    "json5 module not found. JSON5 config files will fail to parse.\n"
    "  Hint: activate the venv — source .venv/bin/activate",
    stacklevel=2,
  )

# Default config file path (relative to this script)
DEFAULT_CONFIG_FILE = Path(__file__).parent / "auto_plot_default.json5"


def load_json_file(filepath: str) -> Dict:
  """Load JSON or JSON5 file."""
  with open(filepath, "r") as f:
    content = f.read()
  if HAS_JSON5:
    return json5.loads(content)
  if filepath.endswith(".json5"):
    raise ImportError(
      f"Cannot parse JSON5 file without json5 module: {filepath}\n"
      "  Install: pip install json5\n"
      "  Or activate venv: source .venv/bin/activate"
    )
  return json.loads(content)


def load_default_config() -> Dict:
  """Load default configuration from auto_plot_default.json5."""
  if not DEFAULT_CONFIG_FILE.exists():
    print(f"ERROR: Default config file not found: {DEFAULT_CONFIG_FILE}", file=sys.stderr)
    sys.exit(1)
  return load_json_file(str(DEFAULT_CONFIG_FILE))


def is_excluded(filename: str, config: Dict, section: str, cli_excludes: List[str]) -> bool:
  """Check if filename matches any exclude pattern.

  Args:
    filename: File path or basename to check
    config: Configuration dictionary
    section: Section name (hist2d, g3vox, g2bg)
    cli_excludes: List of exclude patterns from CLI --exclude option

  Returns:
    True if file should be excluded, False otherwise
  """
  basename = os.path.basename(filename)
  # Check global exclude_patterns
  for pattern in config.get("exclude_patterns", []):
    if fnmatch(basename, pattern):
      return True
  # Check section-specific exclude_patterns
  for pattern in config.get(section, {}).get("exclude_patterns", []):
    if fnmatch(basename, pattern):
      return True
  # Check CLI --exclude patterns
  for pattern in cli_excludes:
    if fnmatch(basename, pattern):
      return True
  return False


def find_script_dir() -> Path:
  """Find the script directory (where this script and other plot scripts live)."""
  return Path(__file__).parent.resolve()


def match_pattern(filename: str, patterns: Dict) -> Optional[Dict]:
  """Match filename against patterns and return config if matched."""
  basename = os.path.basename(filename)
  for pattern, config in patterns.items():
    if fnmatch(basename, pattern + ".tmp"):
      return config
  return None


def parse_triplet(triplet_str: str):
  """Parse triplet string like '320 -1.6 1.6' into (n, min, max)."""
  parts = triplet_str.replace(",", " ").split()
  if len(parts) != 3:
    raise ValueError(f"Invalid triplet format: {triplet_str}")
  return int(parts[0]), float(parts[1]), float(parts[2])


def build_hist2d_command(
  script_dir: Path,
  filepath: str,
  config: Dict,
  use_log: bool = False,
  save_csv: bool = True
) -> str:
  """Build command string for hist2d.py."""
  script_path = script_dir / "hist2d.py"

  # Parse triplet format from config (for backward compatibility with existing json5)
  binz = config.get("binz", "20 0 2000")
  binx = config.get("binx", "320 -1.6 1.6")
  biny = config.get("biny", "160 0.0 1.6")
  dpi = config.get("dpi", 300)
  colormap = config.get("colormap", "jet")
  color_over = config.get("color_over", "black")
  color_under = config.get("color_under", "white")

  ngrad, vmin, vmax = parse_triplet(binz)
  nbinx, xmin, xmax = parse_triplet(binx)
  nbiny, ymin, ymax = parse_triplet(biny)

  # Build output filename from input filename
  basename = os.path.basename(filepath)
  output_name = basename.replace(".tmp", "") + ("_log" if use_log else "") + ".png"

  # Build command string for shell execution
  cmd_parts = [
    "python3", str(script_path), filepath,
    "--nbinx", nbinx, "--xmin", xmin, "--xmax", xmax,
    "--nbiny", nbiny, "--ymin", ymin, "--ymax", ymax,
    "--ngrad", ngrad, "--vmin", vmin, "--vmax", vmax,
    "--dpi", dpi,
    "--colormap", colormap,
    "--color-over", color_over,
    "--color-under", color_under,
    "--output", output_name,
  ]
  if use_log:
    cmd_parts.append("--log")
  if not save_csv:
    cmd_parts.append("--no-csv")
  return " ".join(shlex.quote(str(part)) for part in cmd_parts)


def build_g2bg_command(script_dir: Path, json_file: str) -> List[str]:
  """Build command for plot_g2bg.py."""
  script_path = script_dir / "plot_g2bg.py"
  return ["python3", str(script_path), f"--jsonfile={json_file}"]


def build_g3vox_command(script_dir: Path, json_file: str, bin_mode: str = None, keep_png: bool = True) -> List[str]:
  """Build command for plot_Grid3dVoxel.py.

  Args:
    script_dir: Path to script directory
    json_file: Path to JSON config file
    bin_mode: Binning mode - "data" (use data coordinates), "header" (use file header), or None
    keep_png: If True, keep intermediate PNG files (default: True)
  """
  script_path = script_dir / "plot_Grid3dVoxel.py"
  cmd = ["python3", str(script_path), f"--config={json_file}"]
  if bin_mode == "data":
    cmd.append("--use-data-bins")
  elif bin_mode == "header":
    cmd.append("--use-header-bins")
  if not keep_png:
    cmd.append("--delete-png")
  return cmd


def load_g3vox_template(config: Dict) -> Dict:
  """Load g3vox template from config."""
  import copy
  template = config.get("g3vox", {}).get("template", {})
  if template:
    return copy.deepcopy(template)
  # Fallback to default config template
  default_config = load_default_config()
  return copy.deepcopy(default_config.get("g3vox", {}).get("template", {}))


def match_g3vox_render_pattern(filename: str, render_patterns: Dict) -> Optional[Dict]:
  """Match filename against g3vox render_patterns and return config if matched."""
  for pattern, settings in render_patterns.items():
    if fnmatch(filename, pattern):
      return settings
  return None


def resolve_det_info_path(det_info_path) -> Optional[List[str]]:
  """Resolve det_info_path to a list of detector file paths.

  Args:
    det_info_path: Either a list of detector file paths, or a path to a JSON/JSON5
                   file containing a 'det_files' array.

  Returns:
    List of detector file paths, or None if not specified or invalid.
  """
  if det_info_path is None:
    return None

  # If already a list, return as-is
  if isinstance(det_info_path, list):
    return det_info_path

  # If a string ending with .json or .json5, load and extract det_files
  if isinstance(det_info_path, str):
    if det_info_path.endswith((".json", ".json5")):
      if os.path.exists(det_info_path):
        try:
          data = load_json_file(det_info_path)
          # Check for det_files at root level
          det_files = data.get("det_files")
          if isinstance(det_files, list):
            return det_files
          # Check for det_files under DETECTOR_PARAMETER_LISTS
          det_params = data.get("DETECTOR_PARAMETER_LISTS", {})
          det_files = det_params.get("det_files")
          if isinstance(det_files, list):
            return det_files
          print(f"WARNING: {det_info_path} does not contain 'det_files' array", file=sys.stderr)
          return None
        except Exception as e:
          print(f"WARNING: Failed to load {det_info_path}: {e}", file=sys.stderr)
          return None
      else:
        print(f"ERROR: det_info_path file not found: {det_info_path}", file=sys.stderr)
        sys.exit(1)

  return None


def override_save_png_dir(config: Dict, save_png_dir: str) -> Dict:
  """Create a config copy with overridden save_png_dir for g3vox template."""
  import copy
  new_config = copy.deepcopy(config)
  g3vox = new_config.get("g3vox", {})
  template = g3vox.get("template", {})
  render = template.get("render_params", {})
  render["save_png_dir"] = save_png_dir
  return new_config


def generate_g3vox_json(tmp_file: str, config: Dict) -> str:
  """Generate a temporary JSON config for g3vox .tmp file."""
  basename = os.path.basename(tmp_file)
  name_without_ext = (basename.replace("_zcross_all.tmpbin", "")
                      .replace("_zcross_all.tmp", "")
                      .replace(".tmpbin", "").replace(".tmp", ""))

  # Create config from template
  g3vox_config = load_g3vox_template(config)
  g3vox_config["run_params"]["filename_in"] = tmp_file
  g3vox_config["render_params"]["output_png_prefix"] = f"fig_zcross_{name_without_ext}"

  # Apply det_info_path if specified (resolve from JSON file if needed)
  # Note: det_info_path goes at root level of config, not under run_params
  det_info_path = config.get("g3vox", {}).get("det_info_path")
  resolved_det_files = resolve_det_info_path(det_info_path)
  if resolved_det_files:
    g3vox_config["det_info_path"] = resolved_det_files

  # Apply render_patterns if specified
  render_patterns = config.get("g3vox", {}).get("render_patterns", {})
  matched_settings = match_g3vox_render_pattern(basename, render_patterns)
  if matched_settings:
    for key, value in matched_settings.items():
      g3vox_config["render_params"][key] = value

  # Write temporary JSON (include sweep dir name to avoid collisions)
  parent_dir = os.path.basename(os.path.dirname(tmp_file))
  if parent_dir.startswith("sweep_"):
    json_filename = f"_auto_{parent_dir}_{name_without_ext}.json"
  else:
    json_filename = f"_auto_{name_without_ext}.json"
  with open(json_filename, "w") as f:
    json.dump(g3vox_config, f, indent=2)

  return json_filename


def load_g2bg_template(config: Dict) -> Dict:
  """Load g2bg template from config."""
  import copy
  template = config.get("g2bg", {}).get("template", {})
  if template:
    return copy.deepcopy(template)
  # Fallback to default config template
  default_config = load_default_config()
  return copy.deepcopy(default_config.get("g2bg", {}).get("template", {}))


def generate_g2bg_json(prefix: str, det_ids: List[str], config: Dict) -> str:
  """Generate a temporary JSON config for g2bg files.

  Args:
    prefix: Input file prefix (e.g., "g2bg_arrdet_g3vox_input")
    det_ids: List of detector IDs (e.g., ["00", "01", "02"])
    config: Configuration dictionary

  Returns:
    Path to generated JSON file
  """
  g2bg_config = load_g2bg_template(config)

  # Set input/output parameters
  g2bg_config["input_prefix"] = prefix + "_"
  g2bg_config["det_ids"] = sorted(det_ids)
  g2bg_config["output_prefix"] = f"fig_{prefix}_"
  g2bg_config["pdf_prefix"] = f"fig_{prefix}_"
  g2bg_config["gif_prefix"] = f"fig_{prefix}_"

  # Write temporary JSON
  json_filename = f"_auto_{prefix}.json"
  with open(json_filename, "w") as f:
    json.dump(g2bg_config, f, indent=2)

  return json_filename


def run_command(cmd: List[str], dry_run: bool = False, timeout: int = 300) -> Tuple[List[str], int, str]:
  """Run a command and return (cmd, returncode, output)."""
  if dry_run:
    return cmd, 0, "[dry-run]"

  try:
    result = subprocess.run(
      cmd,
      capture_output=True,
      text=True,
      timeout=timeout  # seconds (default 300)
    )
    output = result.stdout + result.stderr
    return cmd, result.returncode, output
  except subprocess.TimeoutExpired:
    return cmd, -1, "TIMEOUT"
  except Exception as e:
    return cmd, -1, str(e)


def collect_hist2d_tasks(
  script_dir: Path,
  config: Dict,
  tmp_dir: str = "tmp",
  cli_excludes: List[str] = None,
  verbose: bool = False
) -> Tuple[List[Tuple[str, str]], int]:
  """Collect hist2d plotting tasks (returns command strings for shell execution).

  Returns:
    Tuple of (tasks list, number of excluded files)
  """
  if cli_excludes is None:
    cli_excludes = []
  tasks = []
  n_excluded = 0

  # Check exec flag at top level (default: True for backward compatibility)
  if not config.get("hist2d_exec", True):
    return tasks, n_excluded

  hist2d_config = config.get("hist2d", {})
  template = hist2d_config.get("template", {})
  patterns = hist2d_config.get("patterns", {})
  save_csv = hist2d_config.get("save_csv", True)  # Default: save CSV

  # Find all .tmp files in tmp directory
  tmp_files = glob.glob(os.path.join(tmp_dir, "*.tmp"))

  for filepath in sorted(tmp_files):
    basename = os.path.basename(filepath)

    # Check exclude patterns first
    if is_excluded(basename, config, "hist2d", cli_excludes):
      n_excluded += 1
      if verbose:
        print(f"  [hist2d] Excluded: {basename}")
      continue

    matched_config = match_pattern(basename, patterns)

    if matched_config:
      merged = {**template, **matched_config}
      use_log = merged.get("log", False)
      cmd = build_hist2d_command(script_dir, filepath, merged, use_log, save_csv)
      desc = f"hist2d {'(log)' if use_log else '(linear)'}: {basename}"
      tasks.append((cmd, desc))

  return tasks, n_excluded


def collect_g2bg_tasks(
  script_dir: Path,
  config: Dict,
  tmp_dir: str = "tmp",
  cli_excludes: List[str] = None,
  verbose: bool = False
) -> Tuple[List[Tuple[List[str], str]], List[str], int]:
  """Collect g2bg plotting tasks by scanning tmp/g2bg_*.tmp files.

  Returns:
    Tuple of (tasks list, generated JSON files, number of excluded files)
  """
  import re
  if cli_excludes is None:
    cli_excludes = []
  tasks = []
  generated_jsons = []
  n_excluded = 0
  g2bg_config = config.get("g2bg", {})

  # Check exec flag at top level (default: True for backward compatibility)
  if not config.get("g2bg_exec", True):
    return tasks, generated_jsons, n_excluded

  if not g2bg_config.get("auto_detect", True):
    return tasks, generated_jsons, n_excluded

  # Find all g2bg_*.tmp files and group by prefix
  tmp_pattern = os.path.join(tmp_dir, "g2bg_*.tmp")
  tmp_files = glob.glob(tmp_pattern)

  # Group files by prefix (e.g., g2bg_arrdet_g3vox_input_det00.tmp -> g2bg_arrdet_g3vox_input)
  prefix_groups = {}  # prefix -> list of det_ids
  det_pattern = re.compile(r"^(.+)_det(\d+)\.tmp$")

  for tmp_file in sorted(tmp_files):
    basename = os.path.basename(tmp_file)

    # Check exclude patterns
    if is_excluded(basename, config, "g2bg", cli_excludes):
      n_excluded += 1
      if verbose:
        print(f"  [g2bg] Excluded: {basename}")
      continue

    match = det_pattern.match(basename)
    if match:
      prefix = match.group(1)
      det_id = match.group(2)
      if prefix not in prefix_groups:
        prefix_groups[prefix] = []
      prefix_groups[prefix].append(det_id)

  # Generate JSON config for each prefix group
  for prefix in sorted(prefix_groups.keys()):
    det_ids = prefix_groups[prefix]
    json_file = generate_g2bg_json(prefix, det_ids, config)
    generated_jsons.append(json_file)
    cmd = build_g2bg_command(script_dir, json_file)
    desc = f"g2bg: {prefix} ({len(det_ids)} detectors)"
    tasks.append((cmd, desc))

  return tasks, generated_jsons, n_excluded


def collect_g3vox_tasks(
  script_dir: Path,
  config: Dict,
  tmp_dir: str = "tmp",
  cli_excludes: List[str] = None,
  verbose: bool = False,
  keep_png: bool = True
) -> Tuple[List[Tuple[List[str], str]], List[str], int]:
  """Collect g3vox plotting tasks from JSON configs and .tmp files.

  Returns:
    Tuple of (tasks list, generated JSON files, number of excluded files)
  """
  if cli_excludes is None:
    cli_excludes = []
  tasks = []
  generated_jsons = []
  n_excluded = 0
  g3vox_config = config.get("g3vox", {})

  # Check exec flag at top level (default: True for backward compatibility)
  if not config.get("g3vox_exec", True):
    return tasks, generated_jsons, n_excluded

  if not g3vox_config.get("auto_detect", True):
    return tasks, generated_jsons, n_excluded

  # Find all g3vox_*.tmp files and generate configs from template
  tmp_files = sorted(
      glob.glob(os.path.join(tmp_dir, "g3vox_*_zcross_all.tmp")) +
      glob.glob(os.path.join(tmp_dir, "g2pil_*_zcross_all.tmp")) +
      glob.glob(os.path.join(tmp_dir, "g3vox_*_zcross_all.tmpbin")) +
      glob.glob(os.path.join(tmp_dir, "g2pil_*_zcross_all.tmpbin"))
  )

  # Check auto_binning mode: "data", "header", or false/None
  auto_binning = g3vox_config.get("auto_binning", False)
  if auto_binning == "data":
    bin_mode = "data"
    print("[g3vox] auto_binning: 'data' - will use --use-data-bins")
  elif auto_binning == "header":
    bin_mode = "header"
    print("[g3vox] auto_binning: 'header' - will use --use-header-bins")
  else:
    bin_mode = None

  for tmp_file in sorted(tmp_files):
    basename = os.path.basename(tmp_file)

    # Check exclude patterns
    if is_excluded(basename, config, "g3vox", cli_excludes):
      n_excluded += 1
      if verbose:
        print(f"  [g3vox] Excluded: {basename}")
      continue

    # Generate temporary JSON config from template
    json_file = generate_g3vox_json(tmp_file, config)
    generated_jsons.append(json_file)
    cmd = build_g3vox_command(script_dir, json_file, bin_mode=bin_mode, keep_png=keep_png)
    desc = f"g3vox: {basename}"
    tasks.append((cmd, desc))

  return tasks, generated_jsons, n_excluded


def collect_det_tasks(
  script_dir: Path,
  config: Dict,
  config_path: Optional[str] = None,
  tmp_dir: str = "tmp",
  cli_excludes: List[str] = None,
  verbose: bool = False
) -> Tuple[List[Tuple[List[str], str]], int]:
  """Collect det tasks: binary-direct figures via plot_det_arrdet.py.

  Scans det/arrdet_*.bin and checkpoint_*/det/arrdet_*.bin in the working
  directory and under tmp_dir (the sweep driver writes its checkpoint
  bundles inside the tmp directory, e.g. tmp/checkpoint_m4/det/).
  Disabled by default: runs only when the config sets det_exec: true.
  Unlike the other families this reads the ~160 MB binary checkpoints
  directly, so no tmp/*.tmp export is needed.

  Returns:
    Tuple of (tasks list, number of excluded files)
  """
  if cli_excludes is None:
    cli_excludes = []
  tasks = []
  n_excluded = 0

  # det_exec defaults to False: with tmp/*.tmp exports present the same
  # figures would be produced twice (see docs/reference/scripts.md).
  if not config.get("det_exec", False):
    return tasks, n_excluded

  det_config = config.get("det", {})
  if not det_config.get("auto_detect", True):
    return tasks, n_excluded

  bin_files = sorted(
      glob.glob(os.path.join("det", "arrdet_*.bin")) +
      glob.glob(os.path.join("checkpoint_*", "det", "arrdet_*.bin")) +
      glob.glob(os.path.join(tmp_dir, "det", "arrdet_*.bin")) +
      glob.glob(os.path.join(tmp_dir, "checkpoint_*", "det", "arrdet_*.bin"))
  )

  script_path = script_dir / "plot_det_arrdet.py"
  for bin_file in bin_files:
    basename = os.path.basename(bin_file)

    # Check exclude patterns
    if is_excluded(basename, config, "det", cli_excludes):
      n_excluded += 1
      if verbose:
        print(f"  [det] Excluded: {basename}")
      continue

    cmd = ["python3", str(script_path), bin_file]
    if config_path:
      cmd.extend(["--plot-config", config_path])
    desc = f"det: {bin_file}"
    tasks.append((cmd, desc))

  return tasks, n_excluded


def load_config(config_path: Optional[str]) -> Dict:
  """Load configuration from file or use default."""
  import copy
  config = copy.deepcopy(load_default_config())

  if config_path:
    if not os.path.exists(config_path):
      print(f"ERROR: Config file not found: {config_path}", file=sys.stderr)
      sys.exit(1)
    user_config = load_json_file(config_path)
    # Merge user config into default (user overrides default)
    for key, value in user_config.items():
      if isinstance(value, dict) and key in config:
        config[key].update(value)
      else:
        config[key] = value

  return config


def main():
  parser = argparse.ArgumentParser(
    description="Auto-plot muonith outputs",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=__doc__
  )
  parser.add_argument(
    "--config", "-c",
    help="Custom config file (JSON5/JSON) to override defaults"
  )
  parser.add_argument(
    "--dry-run", "-n",
    action="store_true",
    help="Print commands without executing"
  )
  parser.add_argument(
    "--skip-hist2d",
    action="store_true",
    help="Skip hist2d plots"
  )
  parser.add_argument(
    "--skip-g2bg",
    action="store_true",
    help="Skip g2bg plots"
  )
  parser.add_argument(
    "--skip-g3vox",
    action="store_true",
    help="Skip g3vox plots"
  )
  parser.add_argument(
    "--skip-det",
    action="store_true",
    help="Skip det plots (binary-direct plot_det_arrdet.py tasks)"
  )
  parser.add_argument(
    "--tmp-dir",
    default="tmp",
    help="Directory containing .tmp files (default: tmp)"
  )
  parser.add_argument(
    "--verbose", "-v",
    action="store_true",
    help="Verbose output"
  )
  parser.add_argument(
    "--exclude", "-e",
    action="append",
    default=[],
    metavar="PATTERN",
    help="Exclude files matching pattern (can be specified multiple times)"
  )
  parser.add_argument(
    "--keep-png", "--keep_png",
    action="store_true",
    default=True,
    help="Keep intermediate PNG files in figs/ subdirectories (default: True)"
  )
  parser.add_argument(
    "--delete-png", "--delete_png",
    action="store_true",
    help="Delete intermediate PNG files after GIF/PDF creation"
  )
  parser.add_argument(
    "--error-log",
    default="error_auto_plot-tmp.log",
    metavar="FILE",
    help="File to write error logs (default: error_auto_plot-tmp.log). Use empty string to disable."
  )

  args = parser.parse_args()

  # Load configuration
  config = load_config(args.config)
  script_dir = find_script_dir()

  # Collect tasks
  all_tasks = []
  generated_jsons = []
  total_excluded = 0

  if not args.skip_g3vox:
    g3vox_tasks, gen_jsons, n_excl = collect_g3vox_tasks(
      script_dir, config, args.tmp_dir, args.exclude, args.verbose, keep_png=not args.delete_png)
    all_tasks.extend(g3vox_tasks)
    generated_jsons.extend(gen_jsons)
    total_excluded += n_excl
    if g3vox_tasks or n_excl:
      n_auto = len(gen_jsons)
      n_existing = len(g3vox_tasks) - n_auto
      excl_msg = f", {n_excl} excluded" if n_excl else ""
      print(f"[g3vox] Found {len(g3vox_tasks)} task(s) ({n_existing} from JSON, {n_auto} auto-generated{excl_msg})")

    # Sweep subdirectory support for g3vox
    sweep_dirs = sorted(glob.glob(os.path.join(args.tmp_dir, "sweep_*")))
    sweep_dirs = [d for d in sweep_dirs if os.path.isdir(d)]
    for sweep_dir in sweep_dirs:
      sweep_name = os.path.basename(sweep_dir)
      sweep_figs_dir = os.path.join("figs", sweep_name)
      sweep_config = override_save_png_dir(config, sweep_figs_dir)
      sw_tasks, sw_jsons, sw_excl = collect_g3vox_tasks(
        script_dir, sweep_config, sweep_dir, args.exclude, args.verbose, keep_png=not args.delete_png)
      all_tasks.extend(sw_tasks)
      generated_jsons.extend(sw_jsons)
      total_excluded += sw_excl
      if sw_tasks or sw_excl:
        excl_msg = f", {sw_excl} excluded" if sw_excl else ""
        print(f"[g3vox] {sweep_name}: {len(sw_tasks)} task(s){excl_msg}")

  if not args.skip_g2bg:
    g2bg_tasks, g2bg_gen_jsons, n_excl = collect_g2bg_tasks(
      script_dir, config, args.tmp_dir, args.exclude, args.verbose)
    all_tasks.extend(g2bg_tasks)
    generated_jsons.extend(g2bg_gen_jsons)
    total_excluded += n_excl
    if g2bg_tasks or n_excl:
      excl_msg = f", {n_excl} excluded" if n_excl else ""
      print(f"[g2bg] Found {len(g2bg_tasks)} task(s) ({len(g2bg_gen_jsons)} auto-generated{excl_msg})")

  if not args.skip_hist2d:
    hist2d_tasks, n_excl = collect_hist2d_tasks(
      script_dir, config, args.tmp_dir, args.exclude, args.verbose)
    all_tasks.extend(hist2d_tasks)
    total_excluded += n_excl
    if hist2d_tasks or n_excl:
      excl_msg = f", {n_excl} excluded" if n_excl else ""
      print(f"[hist2d] Found {len(hist2d_tasks)} .tmp file(s){excl_msg}")

  if not args.skip_det:
    det_tasks, n_excl = collect_det_tasks(
      script_dir, config, args.config, args.tmp_dir, args.exclude, args.verbose)
    all_tasks.extend(det_tasks)
    total_excluded += n_excl
    if det_tasks or n_excl:
      excl_msg = f", {n_excl} excluded" if n_excl else ""
      print(f"[det] Found {len(det_tasks)} arrdet binary file(s){excl_msg}")

  if not all_tasks:
    print("No files found to plot.")
    return 0

  print(f"\nTotal tasks: {len(all_tasks)}")
  print("-" * 60)

  # Separate sequential and parallel tasks
  # Check prefix of desc to avoid matching filenames containing "g3vox"
  sequential_tasks = [(cmd, desc) for cmd, desc in all_tasks
                      if desc.startswith("g3vox") or desc.startswith("g2bg")
                      or desc.startswith("det")]
  parallel_tasks = [(cmd, desc) for cmd, desc in all_tasks
                    if desc.startswith("hist2d")]

  # Dry-run: just print commands
  if args.dry_run:
    for cmd, desc in sequential_tasks:
      print(f"[DRY-RUN] {desc}")
      cmd_str = ' '.join(cmd) if isinstance(cmd, list) else cmd
      print(f"  {cmd_str}")
    for cmd, desc in parallel_tasks:
      print(f"[DRY-RUN] {desc}")
      cmd_str = cmd if isinstance(cmd, str) else ' '.join(cmd)
      print(f"  {cmd_str} &")
    if parallel_tasks:
      print("[DRY-RUN] wait")
    # Cleanup generated JSON files even in dry-run
    for json_file in generated_jsons:
      try:
        os.remove(json_file)
      except Exception:
        pass
    return 0

  # Execute tasks
  n_success = 0
  n_failed = 0
  error_log_file = args.error_log if args.error_log else None

  # Initialize error log file if enabled
  if error_log_file:
    try:
      with open(error_log_file, "w") as f:
        from datetime import datetime
        f.write(f"# auto_plot.py error log - {datetime.now().isoformat()}\n")
        f.write(f"# Working directory: {os.getcwd()}\n\n")
    except Exception as e:
      print(f"WARNING: Could not create error log file: {e}", file=sys.stderr)
      error_log_file = None

  def write_error_log(desc: str, cmd, retcode: int, output: str):
    """Write error details to log file."""
    if not error_log_file:
      return
    try:
      with open(error_log_file, "a") as f:
        f.write("=" * 70 + "\n")
        f.write(f"TASK: {desc}\n")
        cmd_str = ' '.join(cmd) if isinstance(cmd, list) else cmd
        f.write(f"COMMAND: {cmd_str}\n")
        f.write(f"RETURN CODE: {retcode}\n")
        f.write("-" * 70 + "\n")
        f.write("OUTPUT:\n")
        f.write(output.strip() if output else "(no output)")
        f.write("\n\n")
    except Exception:
      pass

  def write_warning_log(desc: str, cmd, output: str):
    """Write warning details to log file (for successful tasks with warnings)."""
    if not error_log_file:
      return
    # Check if output contains warnings
    if not output or "[WARN]" not in output:
      return
    try:
      with open(error_log_file, "a") as f:
        f.write("=" * 70 + "\n")
        f.write(f"TASK: {desc} (SUCCESS with warnings)\n")
        cmd_str = ' '.join(cmd) if isinstance(cmd, list) else cmd
        f.write(f"COMMAND: {cmd_str}\n")
        f.write("-" * 70 + "\n")
        f.write("WARNINGS:\n")
        # Extract only warning lines
        warning_lines = [line for line in output.split('\n') if '[WARN]' in line]
        f.write('\n'.join(warning_lines))
        f.write("\n\n")
    except Exception:
      pass

  # Sequential execution for g3vox/g2bg/det
  for cmd, desc in sequential_tasks:
    print(f"Running: {desc}")
    # det tasks read ~160 MB binaries per detector array; allow 600 s.
    task_timeout = 600 if desc.startswith("det") else 300
    _, retcode, output = run_command(cmd, timeout=task_timeout)
    if retcode == 0:
      n_success += 1
      # Record warnings even on success
      write_warning_log(desc, cmd, output)
    else:
      n_failed += 1
      print(f"  FAILED (retcode={retcode})")
      # Always show error output on failure (truncated), full output with --verbose
      if output:
        max_len = 2000 if args.verbose else 500
        trimmed = output.strip()
        if len(trimmed) > max_len:
          print(f"  Error output (truncated):\n{trimmed[:max_len]}...")
        else:
          print(f"  Error output:\n{trimmed}")
      # Write full error to log file
      write_error_log(desc, cmd, retcode, output)

  # Parallel execution for hist2d using shell & and wait
  if parallel_tasks:
    print(f"\nRunning {len(parallel_tasks)} hist2d tasks in parallel...")
    # Track each background job because plain "wait" can hide failures.
    shell_commands = ["set +e", "pids=()"]
    for cmd, desc in parallel_tasks:
      shell_commands.append(f"({cmd}) &")
      shell_commands.append('pids+=("$!")')
    shell_commands.extend([
      "status=0",
      'for pid in "${pids[@]}"; do',
      '  wait "$pid" || status=1',
      "done",
      'exit "$status"',
    ])
    shell_script = "\n".join(shell_commands)

    result = subprocess.run(
      shell_script,
      shell=True,
      executable='/bin/bash',
      capture_output=True,
      text=True
    )
    if result.returncode == 0:
      n_success += len(parallel_tasks)
      print(f"  All {len(parallel_tasks)} hist2d tasks completed")
      # Record warnings even on success
      combined_output = result.stdout + result.stderr
      write_warning_log("hist2d (parallel)", shell_script, combined_output)
    else:
      # Can't determine individual failures with shell &, assume all failed
      n_failed += len(parallel_tasks)
      print(f"  Some hist2d tasks failed (returncode={result.returncode})")
      if args.verbose:
        print(f"  stderr: {result.stderr[:500]}")
      # Write error to log file
      combined_output = result.stdout + result.stderr
      write_error_log("hist2d (parallel)", shell_script, result.returncode, combined_output)

  # Cleanup generated JSON files
  if generated_jsons:
    for json_file in generated_jsons:
      try:
        os.remove(json_file)
      except Exception:
        pass

  # Final cleanup: remove intermediate PNG files in subdirectories only
  # NOTE: Do NOT delete PNG files in figs/ directly - they are outputs from hist2d/g2bg
  if args.delete_png:
    figs_dir = "figs"
    if os.path.isdir(figs_dir):
      png_count = 0
      # Only remove PNG files in subdirectories (g3vox intermediate files)
      for png_file in glob.glob(os.path.join(figs_dir, "*", "*.png")):
        try:
          os.remove(png_file)
          png_count += 1
        except Exception:
          pass
      # Remove empty subdirectories
      for subdir in glob.glob(os.path.join(figs_dir, "*")):
        if os.path.isdir(subdir) and not os.listdir(subdir):
          try:
            os.rmdir(subdir)
          except Exception:
            pass
      if png_count > 0:
        print(f"[INFO] Removed {png_count} intermediate PNG files from {figs_dir}/ subdirectories")
  else:
    print("[INFO] Keeping intermediate PNG files (--keep-png)")

  # Summary
  print("-" * 60)
  print(f"Completed: {n_success} success, {n_failed} failed")

  # Handle error log file
  if error_log_file:
    try:
      # Check if file has content (beyond header lines)
      with open(error_log_file, "r") as f:
        content = f.read()
      # File has only header if no task entries (no "=" separator lines)
      has_entries = "=" * 70 in content
      if has_entries:
        if n_failed > 0:
          print(f"Error details written to: {error_log_file}")
        else:
          print(f"Warning details written to: {error_log_file}")
      else:
        # Remove log file if no entries
        os.remove(error_log_file)
    except Exception:
      pass

  return 0 if n_failed == 0 else 1


if __name__ == "__main__":
  sys.exit(main())
