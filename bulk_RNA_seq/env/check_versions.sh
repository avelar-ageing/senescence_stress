#!/usr/bin/env bash
# check_versions.sh - compare the installed R and Python packages with env/renv.lock and
# env/requirements.txt. Run from bulk_RNA_seq/. Exit 1 on any difference, unless
# BULK_ALLOW_VERSION_DRIFT=1 (then report only).
#
# Why: on 2026-09-23 another project reinstalled msigdbr 25.1.1 over 26.1.0 in the shared
# R user library while this pipeline was running; the Hallmark sets gained 3 rows and
# every set-level output changed. run_pipeline.sh calls this before running anything.
cd "$(dirname "$0")/.."
PY=${TAGE_PYTHON:-.venv/bin/python}
[ -d .Rlib ] && export R_LIBS="$PWD/.Rlib${R_LIBS:+:$R_LIBS}"
bad=0

Rscript --vanilla -e '
  lock <- jsonlite::fromJSON("env/renv.lock")$Packages
  inst <- installed.packages(lib.loc = .libPaths())[, c("Package", "Version", "LibPath")]
  inst <- inst[!duplicated(inst[, "Package"]), , drop = FALSE]   # first on .libPaths() wins
  rownames(inst) <- inst[, "Package"]
  d <- do.call(rbind, lapply(names(lock), function(p) {
    have <- if (p %in% rownames(inst)) inst[p, "Version"] else "MISSING"
    if (have != lock[[p]]$Version) data.frame(package = p, lockfile = lock[[p]]$Version, installed = have,
                                              library = if (p %in% rownames(inst)) inst[p, "LibPath"] else "")
  }))
  if (is.null(d)) cat(sprintf("R: all %d lockfile packages at the locked version\n", length(lock))) else {
    cat(sprintf("R: %d of %d lockfile packages differ:\n", nrow(d), length(lock))); print(d, row.names = FALSE); quit(status = 1) }
' || bad=1

"$PY" - <<'EOF' || bad=1
import sys
from importlib.metadata import version, PackageNotFoundError
want = [l.split("==") for l in open("env/requirements.txt") if "==" in l and not l.startswith("#")]
diff = []
for name, v in want:
    v = v.strip()
    try:
        have = version(name)
    except PackageNotFoundError:
        have = "MISSING"
    if have != v:
        diff.append((name, v, have))
if diff:
    print(f"Python: {len(diff)} of {len(want)} packages differ:", *[f"  {n}: requirements {v}, installed {h}" for n, v, h in diff], sep="\n")
    sys.exit(1)
print(f"Python: all {len(want)} requirements at the pinned version")
EOF

if [ $bad -ne 0 ]; then
  if [ "${BULK_ALLOW_VERSION_DRIFT:-0}" = "1" ]; then echo "(BULK_ALLOW_VERSION_DRIFT=1: continuing)"; exit 0; fi
  echo "Versions differ from the lockfiles. Fix with env/setup_macos.sh (installs into .Rlib),"
  echo "or set BULK_ALLOW_VERSION_DRIFT=1 to run anyway (outputs may then differ)."
  exit 1
fi
