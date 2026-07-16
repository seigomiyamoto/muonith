#!/usr/bin/env bash
# Regenerate the daemonflux dF/dE table (m^-2) from the source momentum table.
#
# Source: org/daemonflux.txt.zip  (columns: log10(P/(GeV/c)), cos(zenith), dF/dP)
#   NOTE: the source header is mislabeled "(m^2 ...)" but the values are per cm^2.
#   dFdE_interp.py applies CM2_TO_M2 = 1e4 so the output is per m^2 (Honda convention).
# Output: daemonflux_ke_m2.tmp on a (cos theta, log10KE) grid,
#   plus the cross-section slice figures under figs/ listed in slice_tmp.json.

set -euo pipefail
cd "$(dirname "$0")"

unzip -o org/daemonflux.txt.zip daemonflux.txt

python3 dFdE_interp.py \
  --infile=daemonflux.txt \
  --outfile=daemonflux_ke_m2.tmp \
  --json=slice_tmp.json
