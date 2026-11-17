#!/bin/sh
# Run the empirical ITRDB Stage B workflow on macOS or Linux.
# Usage: sh scripts/run_itrdb_empirical_stage_b_v01.sh PROJECT_ROOT [--skip-plots]

set -eu

if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
    echo "Usage: sh scripts/run_itrdb_empirical_stage_b_v01.sh PROJECT_ROOT [--skip-plots]" >&2
    exit 2
fi

project_input=$1
skip_plots=0
if [ "$#" -eq 2 ]; then
    if [ "$2" != "--skip-plots" ]; then
        echo "Unknown option: $2" >&2
        exit 2
    fi
    skip_plots=1
fi

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd -P)
project_root=$(CDPATH= cd -- "$project_input" && pwd -P)
python_cmd=${PYTHON:-python3}

if ! command -v "$python_cmd" >/dev/null 2>&1; then
    echo "Python executable not found: $python_cmd (set PYTHON to override)." >&2
    exit 1
fi
if ! command -v Rscript >/dev/null 2>&1; then
    echo "Rscript was not found on PATH." >&2
    exit 1
fi

PYTHONPATH="$repo_root/python/src${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONPATH

printf 'Repository:   %s\nProject root: %s\nPython:       %s\n\n' \
    "$repo_root" "$project_root" "$python_cmd"

echo '[1/6] Audit value-blind complete cohorts'
"$python_cmd" "$repo_root/scripts/audit_itrdb_empirical_cohorts_v01.py" \
    --project-root "$project_root" --repo-root "$repo_root"

echo '[2/6] Detrend cores and freeze tree-level data'
Rscript "$repo_root/scripts/prepare_itrdb_empirical_stage_b_v01.R" \
    --project-root "$project_root" --repo-root "$repo_root"

echo '[3/6] Generate explicit empirical memberships'
"$python_cmd" "$repo_root/scripts/generate_itrdb_empirical_memberships_v01.py" \
    --repo-root "$repo_root"

echo '[4/6] Run composition-aware empirical diagnostics'
"$python_cmd" "$repo_root/scripts/run_itrdb_empirical_stage_b_v01.py" \
    --repo-root "$repo_root"

echo '[5/6] Validate deterministic empirical outputs'
"$python_cmd" "$repo_root/scripts/validate_itrdb_empirical_stage_b_v01.py" \
    --repo-root "$repo_root"

echo '[6/6] Build empirical figures'
if [ "$skip_plots" -eq 1 ]; then
    echo 'Figure generation skipped by request.'
elif "$python_cmd" -c 'import matplotlib' >/dev/null 2>&1; then
    "$python_cmd" "$repo_root/scripts/plot_itrdb_empirical_stage_b_v01.py" \
        --repo-root "$repo_root"
else
    echo "matplotlib is not installed; tables are complete but figures were skipped." >&2
fi

echo 'ITRDB empirical Stage B workflow completed.'
printf 'Results: %s\n' \
    "$repo_root/analysis/itrdb_empirical_stage_b_v01/results/empirical_key_results_v01.csv"
printf 'Figures: %s\n' \
    "$repo_root/analysis/itrdb_empirical_stage_b_v01/results/figures"

