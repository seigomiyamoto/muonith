#!/usr/bin/env bash

# ${work_dir_name}/detparams/$spec_label/mkjson5.sh

script_path="../../../../scripts/mk_detjson5_from_template.py"
template_path="template-det.json5"
csv_path="../${csv_basename}"

uv run python $$script_path \
 --template $$template_path \
 --csv $$csv_path \
 --z_pos_offset $z_pos_offset

