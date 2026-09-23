#!/usr/bin/env bash
# One-time environment setup for the bulk tAge pipeline on macOS (Apple silicon).
# Run from bulk_RNA_seq/:   bash env/setup_macos.sh
#
# Prerequisites, installed by hand once:
#   R 4.6.x arm64          https://cloud.r-project.org/bin/macosx/   (lockfile records 4.6.0)
#   Xcode command-line tools   xcode-select --install
#   gfortran               https://mac.r-project.org/tools/   (for packages built from source)
#   Python 3.12            brew install python@3.12   (or: uv python install 3.12)
#
# What it does:
#   1. .venv/ with the exact Python packages in env/requirements.txt
#   2. the exact R package versions in env/renv.lock (CRAN, Bioconductor 3.23, and
#      Gladyshev-Lab/tAge at the recorded commit) into bulk_RNA_seq/.Rlib, a library of
#      this project only (run_pipeline.sh puts it first), so installs for other projects
#      cannot change these versions
#   3. checks that the three model files are in external_inputs/models/ with the
#      recorded md5 sums
set -euo pipefail
cd "$(dirname "$0")/.."

PY=${PYTHON:-$(command -v python3.12 || true)}
[ -n "$PY" ] || { echo "python3.12 not found (set PYTHON=/path/to/python3.12)"; exit 1; }
echo "== Python: $PY"
[ -d .venv ] || "$PY" -m venv .venv
.venv/bin/pip install --upgrade pip >/dev/null
.venv/bin/pip install -r env/requirements.txt

echo "== R: $(R --version | head -1)"
R --version | head -1 | grep -q "R version 4\.6\." || echo "WARNING: lockfile was made with R 4.6.0"
mkdir -p .Rlib
Rscript -e '
  lib <- normalizePath(".Rlib"); .libPaths(c(lib, .libPaths()))
  if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv", lib = lib, repos = "https://cloud.r-project.org")
  options(renv.config.pak.enabled = FALSE)
  renv::restore(lockfile = "env/renv.lock", library = lib, prompt = FALSE, clean = FALSE)
'
bash env/check_versions.sh

echo "== models"
bash env/check_inputs.sh
echo "== done. Next: bash run_pipeline.sh --dry-run"
