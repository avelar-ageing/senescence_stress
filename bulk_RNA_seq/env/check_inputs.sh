#!/usr/bin/env bash
# Checks the inputs that are not produced by the pipeline and not stored in git.
# Run from bulk_RNA_seq/. Exit code 1 if anything required is missing or differs.
cd "$(dirname "$0")/.."
EXT=${BULK_EXTERNAL_DIR:-external_inputs}
FINAL=../Final   # R/config.R: DATA_DIR = <repo root>/Final
md5() { if command -v md5sum >/dev/null; then md5sum "$1" | cut -d' ' -f1; else md5 -q "$1"; fi; }
bad=0
while IFS=$'\t' read -r path sum note; do
  case "$path" in \#*|"") continue;; esac
  f=${path/#external_inputs/$EXT}; f=${f/#..\/Final/$FINAL}
  if [ ! -f "$f" ]; then echo "MISSING  $f  ($note)"; bad=1
  elif [ "$(md5 "$f")" != "$sum" ]; then echo "DIFFERS  $f  (md5 $(md5 "$f"), expected $sum)"; bad=1
  else echo "ok       $f"; fi
done < env/external_inputs.tsv
exit $bad
