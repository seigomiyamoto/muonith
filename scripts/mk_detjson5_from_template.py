#!/usr/bin/env python3
# scripts/mk_detjson5_from_template.py

"""
Reflects each row of the CSV into the template JSON5, generating <id_det>.json5 files while preserving comments.

Supported:
- X/Y column auto-detection works regardless of EPSG notation (regex + scoring). Examples: "X (m, EPSG:6680)", "Y(m,EPSG:6676)", "Easting [m]", "Northing (m)", "X"
- Output file name comes from the id_det column (or det_name if missing, or derived from name via "candidateXX" -> det_XX if that's missing too)
- Rows with missing values (any of x, y, z, yaw_deg) are safely skipped
- Can be overridden with explicit options: --xcol, --ycol, --epsg, --z_pos_offset

Usage:
  python mk_detjson5_from_template.py --template template-det.json5 --csv asama_table_epsg6676.csv
  python mk_detjson5_from_template.py --template template-det.json5 --csv data.csv --z_pos_offset 10.5
"""

import argparse
import csv
import io
import os
import re
import sys

def parse_args():
  p = argparse.ArgumentParser(description="CSV -> JSON5 (comment-preserving)")
  p.add_argument("--template", required=True, help="Path to the template JSON5 file")
  p.add_argument("--csv", required=True, help="Path to the input CSV file")
  # For overriding flexible detection: if explicitly specified, that takes priority
  p.add_argument("--xcol", help="Explicitly specify the X (Easting) column name (optional)")
  p.add_argument("--ycol", help="Explicitly specify the Y (Northing) column name (optional)")
  p.add_argument("--epsg", type=int, help="Prefer this EPSG code (e.g., 6676, 6680)")
  p.add_argument("--z_pos_offset", type=float, default=0.0, help="Offset value added to the Z coordinate (default: 0.0)")
  return p.parse_args()

def read_text(path):
  with io.open(path, "r", encoding="utf-8") as f:
    return f.read()

def write_text(path, text):
  with io.open(path, "w", encoding="utf-8", newline="\n") as f:
    f.write(text)

def to_det_name_from_candidate(candidate_name):
  """"candidateXX" -> det_XX / None otherwise"""
  m = re.match(r"^\s*candidate\s*(\d{1,3})\s*$", candidate_name or "", flags=re.I)
  if not m:
    return None
  n = int(m.group(1))
  return f"det_{n:02d}"

def sub_field_str(src, key, value):
  """Replace "key": "..." formatted values without breaking comments, etc."""
  pattern = rf'("{re.escape(key)}"\s*:\s*)"[^"]*"'
  def repl(m):
    return m.group(1) + f'"{value}"'
  out, n = re.subn(pattern, repl, src, count=1)
  if n == 0:
    raise ValueError(f'String key "{key}" not found in the template.')
  return out

def sub_field_num(src, key, value):
  """Replace "key": 123.45 formatted numeric values without breaking comments, etc."""
  pattern = rf'("{re.escape(key)}"\s*:\s*)[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
  if value is None or (isinstance(value, float) and (value != value)):
    raise ValueError(f'Value for numeric key "{key}" is missing.')
  if isinstance(value, int):
    val_str = f"{value:d}"
  else:
    val_str = f"{float(value):.10f}".rstrip("0").rstrip(".")
    if val_str == "-0":
      val_str = "0"
    if val_str == "":
      val_str = "0"
  def repl(m):
    return m.group(1) + val_str
  out, n = re.subn(pattern, repl, src, count=1)
  if n == 0:
    raise ValueError(f'Numeric key "{key}" not found in the template.')
  return out

def transform_template(tpl_text, det_id, x, y, z, yaw_deg):
  """Copy the template string and replace only the specified fields (preserving comments)"""
  t = tpl_text
  t = sub_field_str(t, "name", det_id)
  t = sub_field_num(t, "x", x)
  t = sub_field_num(t, "y", y)
  t = sub_field_num(t, "z", z)
  t = sub_field_num(t, "yaw_deg", yaw_deg)
  return t

def read_csv_rows(path):
  with io.open(path, "r", encoding="utf-8-sig", newline="") as f:
    reader = csv.DictReader(f)
    rows = list(reader)
  return rows

def parse_float(cell, field_name, row_name):
  try:
    if cell is None or str(cell).strip() == "":
      raise ValueError("empty")
    return float(cell)
  except Exception as e:
    raise ValueError(f'Failed to parse numeric field "{field_name}" in row "{row_name}": {cell!r}') from e

def _extract_epsg(name: str):
  m = re.search(r"epsg\s*[:\s]\s*(\d{3,5})", name, flags=re.I)
  return int(m.group(1)) if m else None

def _unit_is_m(name: str):
  return bool(re.search(r"\b(m|meter|metre)\b", name, flags=re.I))

def _axis_kind(name: str):
  n = name.lower()
  if re.search(r"\bx(?![a-z])", n) or "east" in n or "easting" in n or "e(?!psg)" in n:
    return "x"
  if re.search(r"\by(?![a-z])", n) or "north" in n or "northing" in n or "n(?!a)" in n:
    return "y"
  return None

def _score_header(name: str, prefer_epsg: int | None, role: str | None):
  score = 0
  epsg = _extract_epsg(name)
  unit_m = _unit_is_m(name)
  axis = _axis_kind(name)
  if role and axis == role: score += 3
  elif axis: score += 2
  if epsg is not None: score += 3
  if unit_m: score += 1
  if re.match(r"^[XY]|^(E(ast)?|N(orth)?)", name, flags=re.I): score += 1
  if prefer_epsg is not None and epsg == prefer_epsg: score += 4
  return score, epsg, axis

def detect_xy_cols(headers, prefer_epsg: int | None = None, override_x: str | None = None, override_y: str | None = None):
  # If explicitly specified, that takes priority (existence check only)
  if override_x:
    if override_x not in headers:
      raise SystemExit(f"--xcol='{override_x}' not found in the CSV header.")
    x_col = override_x
  else:
    x_scores = [(*_score_header(h, prefer_epsg, "x"), h) for h in headers]
    x_scores.sort(reverse=True)
    x_col = x_scores[0][3] if x_scores else None

  if override_y:
    if override_y not in headers:
      raise SystemExit(f"--ycol='{override_y}' not found in the CSV header.")
    y_col = override_y
  else:
    y_scores = [(*_score_header(h, prefer_epsg, "y"), h) for h in headers]
    y_scores.sort(reverse=True)
    y_col = y_scores[0][3] if y_scores else None

  if not override_x and not override_y:
    # Align EPSG codes if possible (search for a pair among the top candidates)
    x_top = x_scores[:5]
    y_top = y_scores[:5]
    best, best_pair = -1, (x_col, y_col, None)
    for sx, ex, ax, hx in x_top:
      for sy, ey, ay, hy in y_top:
        if hx == hy:  # the same column is not allowed
          continue
        pair = sx + sy + (3 if (ex is not None and ey is not None and ex == ey) else 0)
        if pair > best:
          best = pair
          best_pair = (hx, hy, ex if (ex is not None and ex == ey) else ex or ey)
    x_col, y_col, epsg_pair = best_pair

  if not x_col or not y_col:
    raise SystemExit("Failed to auto-detect X/Y columns. Please specify them explicitly with --xcol/--ycol.")
  return x_col, y_col

def derive_id_det(row):
  # Priority: id_det -> det_name -> derived from name
  id_det = (row.get("id_det") or "").strip()
  if id_det:
    return id_det
  det_name = (row.get("det_name") or "").strip()
  if det_name:
    return det_name
  name = (row.get("name") or "").strip()
  dn = to_det_name_from_candidate(name)
  return dn  # if it remains None, skip

def main():
  args = parse_args()
  tpl_text = read_text(args.template)
  rows = read_csv_rows(args.csv)

  if not rows:
    raise SystemExit("CSV is empty.")

  header = list(rows[0].keys())
  # Auto-detect column names (broadly handles EPSG / language / notation variants)
  x_col, y_col = detect_xy_cols(header, prefer_epsg=args.epsg, override_x=args.xcol, override_y=args.ycol)
  z_col = "elevation(m)" if "elevation(m)" in header else None
  yaw_col = "azimuth(deg)" if "azimuth(deg)" in header else None

  print(f"[detect] X: {x_col} | Y: {y_col}")
  epsg_x = _extract_epsg(x_col); epsg_y = _extract_epsg(y_col)
  if epsg_x and epsg_y and epsg_x != epsg_y:
    print(f"[warn] X and Y may have different EPSG codes: X={epsg_x}, Y={epsg_y}")

  missing = [c for c in [x_col or "X(...)", y_col or "Y(...)", z_col or "elevation(m)", yaw_col or "azimuth(deg)"] if not c or c not in header]
  if missing:
    raise SystemExit("Required column(s) not found in CSV: " + ", ".join(missing))

  made = 0
  skipped = 0
  skipped_reason_counts = {"no_id":0, "missing_numeric":0}

  for row in rows:
    det_id = derive_id_det(row)
    if not det_id or re.fullmatch(r"-+", det_id):
      skipped += 1
      skipped_reason_counts["no_id"] += 1
      continue

    try:
      x = parse_float(row.get(x_col), x_col, det_id)
      y = parse_float(row.get(y_col), y_col, det_id)
      z = parse_float(row.get(z_col), z_col, det_id)
      z += args.z_pos_offset
      yaw_deg = parse_float(row.get(yaw_col), yaw_col, det_id)
    except ValueError:
      skipped += 1
      skipped_reason_counts["missing_numeric"] += 1
      continue

    out_text = transform_template(tpl_text, det_id, x, y, z, yaw_deg)
    out_name = f"{det_id}.json5"
    write_text(out_name, out_text)
    made += 1
    print(f"wrote: {out_name}")

  print(f"Generated: {made}, Skipped: {skipped} (no_id: {skipped_reason_counts['no_id']}, missing_numeric: {skipped_reason_counts['missing_numeric']})")

if __name__ == "__main__":
  try:
    main()
  except Exception as e:
    print(f"[ERROR] {e}", file=sys.stderr)
    sys.exit(1)
