#!/bin/bash
# Fetch the PDG muon energy-loss table for standard rock (Groom-Mokhov-Striganov).
# The raw file is not committed to this repository (no explicit license on the PDG
# page); this script pins the 2024 release by md5.
#
# Usage: bash fetch_pdg_muE.sh [material]
#   material defaults to "standard_rock". For other materials (about 350 available,
#   e.g. water_liquid, iron_Fe) the md5 check is skipped.

set -eu

MATERIAL="${1:-standard_rock}"
URL="https://pdg.lbl.gov/2024/AtomicNuclearProperties/MUE/muE_${MATERIAL}.txt"
OUT="muE_${MATERIAL}.txt"
MD5_STANDARD_ROCK="f04593f0f5115073cc45fb67d5ce2cd0"

curl -fsSL -o "${OUT}" "${URL}"
echo "downloaded: ${OUT}"

if [ "${MATERIAL}" = "standard_rock" ]; then
  echo "${MD5_STANDARD_ROCK}  ${OUT}" | md5sum -c -
else
  echo "note: md5 pin exists only for standard_rock; verify ${OUT} manually."
fi
