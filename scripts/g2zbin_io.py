#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Read/write g2zbin format — a compact binary format for uniform 2D grids.

The g2zbin format stores Grid2d metadata (axis definitions) in a header
and only z values in the data section, achieving roughly one third the
size of an explicit per-point format for uniform grids.

File layout:
  Offset  Size  Type       Description
  ------  ----  ---------  ----------------------------------
  0       8     char[8]    Magic: "G2ZBIN\\0\\0"
  8       2     uint16_t   Version = 2
  10      1     uint8_t    z precision: 8 = float64, 4 = float32
  11      1     uint8_t    storage order: 0 = row-major (iy outer, ix inner)
  12      4     bytes      Reserved (zero)
  16      var              Grid2d header (name + x_axis + y_axis)
  16+H    N*P              z values flat array (row-major)

Version history:
  v1: axis min/max stored raw (pre-tf_shift); the reader had to re-apply
      the half-interval shift from external JSON5 flags.
  v2: axis min/max stored canonical (already shifted, true bin edges);
      pass-through read/write, no external shift information needed.

Grid1d binary format (matches C++ Grid1d::save):
  name:     uint64 (8 bytes, LE) length + chars
  nbin:     int32  (4 bytes, LE)
  min:      float64 (8 bytes, LE)
  max:      float64 (8 bytes, LE)
  interval: float64 (8 bytes, LE)
"""

import struct
from pathlib import Path

import numpy as np

# g2zbin magic bytes (8 bytes)
MAGIC = b"G2ZBIN\x00\x00"
VERSION = 2


def axis_centers(axis_info):
  """Return bin-center coordinates of an axis (v2 semantics).

  v2 axis min/max are the outer bin edges, so the center of bin i is
  min + (i + 0.5) * interval.

  Args:
    axis_info: dict with keys nbin, min, interval (as returned by
               read_grid1d / read_g2zbin).

  Returns:
    numpy 1D array of length nbin with bin-center coordinates.
  """
  nbin = int(axis_info["nbin"])
  interval = float(axis_info["interval"])
  return float(axis_info["min"]) + (np.arange(nbin) + 0.5) * interval


def read_grid1d(f):
  """Read a Grid1d from a binary stream.

  The binary layout matches C++ Grid1d::save():
    - name: uint64 length + chars
    - nbin: int32
    - min, max, interval: float64 each

  Args:
    f: Binary file object opened for reading.

  Returns:
    dict with keys: name (str), nbin (int), min (float), max (float),
    interval (float).
  """
  # name: size_t (uint64 LE) + chars
  (name_len,) = struct.unpack("<Q", f.read(8))
  name = f.read(name_len).decode("utf-8") if name_len > 0 else ""
  # nbin: int32
  (nbin,) = struct.unpack("<i", f.read(4))
  # min, max, interval: float64
  (vmin, vmax, interval) = struct.unpack("<ddd", f.read(24))
  return {
    "name": name,
    "nbin": nbin,
    "min": vmin,
    "max": vmax,
    "interval": interval,
  }


def read_grid2d(f):
  """Read a Grid2d header (name + x_axis + y_axis) from a binary stream.

  Args:
    f: Binary file object opened for reading.

  Returns:
    dict with keys: name (str), x_axis (dict), y_axis (dict).
    Each axis dict has keys: name, nbin, min, max, interval.
  """
  # Grid2d name
  (name_len,) = struct.unpack("<Q", f.read(8))
  name = f.read(name_len).decode("utf-8") if name_len > 0 else ""
  x_axis = read_grid1d(f)
  y_axis = read_grid1d(f)
  return {
    "name": name,
    "x_axis": x_axis,
    "y_axis": y_axis,
  }


def write_grid1d(f, info):
  """Write a Grid1d to a binary stream.

  The binary layout matches C++ Grid1d::save():
    - name: uint64 length + chars
    - nbin: int32
    - min, max, interval: float64 each

  Args:
    f: Binary file object opened for writing.
    info: dict with keys: name (str), nbin (int), min (float),
          max (float), interval (float).
  """
  name_bytes = info["name"].encode("utf-8")
  f.write(struct.pack("<Q", len(name_bytes)))
  f.write(name_bytes)
  f.write(struct.pack("<i", info["nbin"]))
  f.write(struct.pack("<ddd", info["min"], info["max"], info["interval"]))


def write_grid2d(f, info):
  """Write a Grid2d header (name + x_axis + y_axis) to a binary stream.

  Args:
    f: Binary file object opened for writing.
    info: dict with keys: name (str), x_axis (dict), y_axis (dict).
          Each axis dict has keys: name, nbin, min, max, interval.
  """
  name_bytes = info["name"].encode("utf-8")
  f.write(struct.pack("<Q", len(name_bytes)))
  f.write(name_bytes)
  write_grid1d(f, info["x_axis"])
  write_grid1d(f, info["y_axis"])


def read_g2zbin(path, expected_version=VERSION):
  """Read a g2zbin file and return grid metadata and z values.

  Args:
    path: Path to the .g2zbin file (str or Path).
    expected_version: File format version to accept (default: current VERSION).
                      Pass 1 only for migration tooling that converts legacy
                      files (v1 axis min/max are raw, pre-shift values).

  Returns:
    (grid2d_info, Z) where:
      grid2d_info: dict with keys name, x_axis, y_axis (v2: canonical
                   bin-edge coordinates).
      Z: numpy 2D array of shape (nbiny, nbinx) containing z values.

  Raises:
    ValueError: If magic bytes or version are invalid.
  """
  path = Path(path)
  with open(path, "rb") as f:
    # Fixed header (16 bytes)
    magic = f.read(8)
    if magic != MAGIC:
      raise ValueError(f"Invalid magic: expected {MAGIC!r}, got {magic!r}")
    (version,) = struct.unpack("<H", f.read(2))
    if version != expected_version:
      raise ValueError(
        f"Unsupported version: {version} (expected {expected_version}). "
        "v1 files store pre-shift coordinates; regenerate the file from "
        "its source data with the generating script"
      )
    (z_precision,) = struct.unpack("<B", f.read(1))
    (storage_order,) = struct.unpack("<B", f.read(1))
    _reserved = f.read(4)  # skip reserved bytes

    if z_precision not in (4, 8):
      raise ValueError(f"Invalid z precision: {z_precision} (expected 4 or 8)")
    if storage_order != 0:
      raise ValueError(f"Unsupported storage order: {storage_order} (only 0=row-major supported)")

    # Grid2d header
    grid2d_info = read_grid2d(f)

    # z values
    nbinx = grid2d_info["x_axis"]["nbin"]
    nbiny = grid2d_info["y_axis"]["nbin"]
    n_total = nbinx * nbiny
    dtype = np.float32 if z_precision == 4 else np.float64
    Z = np.fromfile(f, dtype=dtype, count=n_total)
    if Z.size != n_total:
      raise ValueError(
        f"Unexpected z data size: expected {n_total}, got {Z.size}"
      )
    Z = Z.reshape(nbiny, nbinx)
    # Always return float64 for consistency
    if z_precision == 4:
      Z = Z.astype(np.float64)

  return grid2d_info, Z


def write_g2zbin(path, grid2d_info, Z, use_float32=True):
  """Write a g2zbin file.

  Args:
    path: Output path (str or Path).
    grid2d_info: dict with keys name, x_axis, y_axis.
                 Each axis dict has keys: name, nbin, min, max, interval.
                 min/max must be canonical bin-edge coordinates (v2 semantics:
                 any half-interval shift already applied).
    Z: numpy 2D array of shape (nbiny, nbinx) containing z values.
    use_float32: If True, store z values as float32 (default). If False, float64.

  Raises:
    ValueError: If Z shape does not match axis nbin values.
  """
  path = Path(path)
  nbinx = grid2d_info["x_axis"]["nbin"]
  nbiny = grid2d_info["y_axis"]["nbin"]
  if Z.shape != (nbiny, nbinx):
    raise ValueError(
      f"Z shape {Z.shape} does not match (nbiny={nbiny}, nbinx={nbinx})"
    )

  z_precision = 4 if use_float32 else 8

  with open(path, "wb") as f:
    # Fixed header (16 bytes)
    f.write(MAGIC)
    f.write(struct.pack("<H", VERSION))
    f.write(struct.pack("<B", z_precision))
    f.write(struct.pack("<B", 0))  # storage order: row-major
    f.write(b"\x00" * 4)  # reserved

    # Grid2d header
    write_grid2d(f, grid2d_info)

    # z values
    if use_float32:
      Z.astype(np.float32).tofile(f)
    else:
      Z.astype(np.float64).tofile(f)
