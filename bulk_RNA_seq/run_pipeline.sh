#!/usr/bin/env bash
# run_pipeline.sh - regenerate every output the manuscript cites, in the order of pipeline.tsv.
#
# Run from bulk_RNA_seq/ (every script sources R/config.R relative to it):
#   bash run_pipeline.sh --dry-run          list what would run and what would be skipped
#   bash run_pipeline.sh                    run every script whose listed outputs are missing
#   bash run_pipeline.sh --force            run every script, overwriting outputs
#   bash run_pipeline.sh --from meta_analysis/18_mortality_clock.py
#                                           force-run that script and everything after it
#   bash run_pipeline.sh --only exploratory/42_coverage_threshold_free.py
#   bash run_pipeline.sh --stage C          only stage A, B or C
#
# A script is skipped when all of its "produces" files already exist. Skipping does not
# check that those files are up to date with their inputs: after changing an upstream
# script, use --from <that script>.
#
# Logs: rerun_outputs/logs/<row>_<script>.log. The driver stops at the first failing script, or
# at the first script that finishes without writing every file listed for it.
#
# Environment (defaults in R/config.R; export to override):
#   TAGE_PYTHON      Python used for .py scripts and by reticulate  [.venv/bin/python]
#   TAGE_MODELS_DIR  the three tAge .pkl files                     [external_inputs/models]
#   FORCE_DOWNLOAD=1 re-fetch recount3 studies instead of the cached raw downloads
#   BULK_ALLOW_VERSION_DRIFT=1  run although packages differ from env/renv.lock / requirements.txt
set -euo pipefail
cd "$(dirname "$0")"
abs() { case "$1" in /*) echo "$1";; *) echo "$PWD/$1";; esac; }

RERUN=rerun_outputs
EXT=${BULK_EXTERNAL_DIR:-external_inputs}
MODELS=${TAGE_MODELS_DIR:-$EXT/models}
PY=${TAGE_PYTHON:-.venv/bin/python}
GEO_CACHE=${GEO_CACHE:-.geo_cache}
XLSX=$EXT/geo/GSE175533_TPM.xlsx
MODELS=$(abs "$MODELS"); XLSX=$(abs "$XLSX")
# R packages: the project library .Rlib (env/setup_macos.sh) comes first, so another
# project's installs into the shared user library cannot change the versions used here
[ -d .Rlib ] && export R_LIBS="$PWD/.Rlib${R_LIBS:+:$R_LIBS}"
export TAGE_PYTHON=$(abs "$PY") TAGE_MODELS_DIR=$MODELS

DRY=0; FORCE=0; FROM=""; ONLY=""; STAGE=""
while [ $# -gt 0 ]; do
  case "$1" in
    --dry-run) DRY=1 ;;
    --force)   FORCE=1 ;;
    --from)    FROM=$2; shift ;;
    --only)    ONLY=$2; shift ;;
    --stage)   STAGE=$2; shift ;;
    -h|--help) sed -n '2,24p' "$0"; exit 0 ;;
    *) echo "unknown option $1"; exit 2 ;;
  esac
  shift
done

if [ $DRY -eq 0 ]; then
  [ -x "$PY" ] || { echo "Python not found at $PY - run env/setup_macos.sh or set TAGE_PYTHON"; exit 1; }
  bash env/check_inputs.sh >/dev/null || { bash env/check_inputs.sh; echo "inputs missing - see REPRODUCE.md"; exit 1; }
  bash env/check_versions.sh || exit 1
fi
mkdir -p "$RERUN/logs"

reached_from=0
[ -z "$FROM" ] && reached_from=1
n_run=0; n_skip=0; row=0
while IFS=$'\t' read -r stage script args produces network; do
  case "$stage" in \#*|stage|"") continue;; esac
  row=$((row+1))
  [ -n "$STAGE" ] && [ "$stage" != "$STAGE" ] && continue
  [ -n "$ONLY" ] && [ "$script" != "$ONLY" ] && continue
  force_this=$FORCE
  if [ -n "$FROM" ]; then
    [ "$script" = "$FROM" ] && reached_from=1
    [ $reached_from -eq 0 ] && continue
    force_this=1
  fi
  [ -n "$ONLY" ] && force_this=1

  missing=0
  for o in $produces; do [ -e "$RERUN/$o" ] || missing=1; done
  if [ $missing -eq 0 ] && [ $force_this -eq 0 ]; then
    printf '  skip  %s %-58s (outputs present)\n' "$stage" "$script"; n_skip=$((n_skip+1)); continue
  fi

  a=${args//\{RERUN\}/$RERUN}; a=${a//\{MODELS\}/$MODELS}
  a=${a//\{GEO_CACHE\}/$GEO_CACHE}; a=${a//\{GSE175533_XLSX\}/$XLSX}
  [ "$a" = "-" ] && a=""
  case "$script" in
    *.R)  cmd=(Rscript "$script") ;;
    *.py) cmd=("$PY" "$script") ;;
  esac
  # shellcheck disable=SC2206
  [ -n "$a" ] && cmd+=($a)

  net=""; [ "$network" != "-" ] && net="  [network: $network]"
  if [ $DRY -eq 1 ]; then
    printf '  RUN   %s %s%s\n' "$stage" "${cmd[*]}" "$net"; n_run=$((n_run+1)); continue
  fi
  log="$RERUN/logs/$(printf %02d $row)_$(echo "$script" | tr '/' '_').log"
  printf '  run   %s %-58s ' "$stage" "$script"
  t0=$(date +%s)
  if ! "${cmd[@]}" >"$log" 2>&1; then
    echo "FAILED after $(( $(date +%s) - t0 ))s - see $log"; tail -20 "$log"; exit 1
  fi
  for o in $produces; do
    [ -e "$RERUN/$o" ] || { echo "finished but did not write $RERUN/$o - see $log"; exit 1; }
  done
  echo "ok ($(( $(date +%s) - t0 ))s)"; n_run=$((n_run+1))
done < pipeline.tsv
[ -n "$FROM" ] && [ $reached_from -eq 0 ] && { echo "--from $FROM: not in pipeline.tsv"; exit 2; }
echo "== $n_run run, $n_skip skipped"
