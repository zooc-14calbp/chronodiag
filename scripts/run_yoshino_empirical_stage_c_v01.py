"""Run the Yoshino empirical Stage C composition-sensitivity analysis."""
from __future__ import annotations

import argparse
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import platform
import sys

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
PYTHON_SRC = REPO_ROOT / "python" / "src"
if str(PYTHON_SRC) not in sys.path:
    sys.path.insert(0, str(PYTHON_SRC))

from chronodiag import (  # noqa: E402
    leave_one_tree_out,
    split_half_stability,
    subset_stability,
    tree_bootstrap,
)


def load_stage_b_runner(repo_root: Path):
    path = repo_root / "scripts" / "run_itrdb_empirical_stage_b_v01.py"
    spec = importlib.util.spec_from_file_location("chronodiag_stage_b_runner", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def stable_manifest_path(path: Path, repo_root: Path) -> str:
    """Return a repository-relative path, or an absolute path outside the repository."""
    resolved = path.resolve()
    try:
        return resolved.relative_to(repo_root.resolve()).as_posix()
    except ValueError:
        return resolved.as_posix()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    parser.add_argument("--preprocessed-dir", type=Path, default=None)
    parser.add_argument("--membership-dir", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()

    repo_root = args.repo_root.resolve()
    analysis_root = repo_root / "analysis" / "yoshino_empirical_stage_c_v01"
    preprocessed_dir = (
        args.preprocessed_dir.resolve()
        if args.preprocessed_dir is not None
        else analysis_root / "preprocessed"
    )
    membership_dir = (
        args.membership_dir.resolve()
        if args.membership_dir is not None
        else analysis_root / "memberships"
    )
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else analysis_root / "results"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    helper = load_stage_b_runner(repo_root)
    contract_path = repo_root / "data" / "yoshino" / "yoshino_dataset_contract_v01.csv"
    keymap_path = repo_root / "data" / "yoshino" / "yoshino_keymap_v01.csv"
    contract = pd.read_csv(contract_path).iloc[0]
    keymap = pd.read_csv(
        keymap_path,
        dtype={
            "series_id": str,
            "tree_id": str,
            "specimen_id": str,
            "source_group": str,
            "site_id": str,
        },
    )
    tree_path = preprocessed_dir / "tree_level_all_methods_v01.csv"
    tree_data = pd.read_csv(tree_path, dtype={"tree_id": str})
    dataset_data = tree_data.loc[tree_data["dataset"] == "yoshino"].copy()
    if dataset_data.empty:
        raise ValueError("no preprocessed Yoshino data")

    structure_rows = [
        helper.structure_row("yoshino", str(method), block)
        for method, block in dataset_data.groupby(
            "detrending_method", sort=True, observed=True
        )
    ]
    sensitivity = helper.chronology_sensitivity(
        "yoshino", dataset_data, helper.PRIMARY_METHOD
    )

    primary = dataset_data.loc[
        dataset_data["detrending_method"] == helper.PRIMARY_METHOD,
        ["year", "tree_id", "value_z"],
    ].rename(columns={"value_z": "value"})
    min_overlap = int(contract["analysis_years"])

    subset_membership = pd.read_csv(
        membership_dir / "yoshino_subset_memberships_v01.csv.gz",
        dtype={"replicate_id": str, "unit_id": str},
    )
    split_membership = pd.read_csv(
        membership_dir / "yoshino_split_half_memberships_v01.csv.gz",
        dtype={"replicate_id": str, "group": str, "unit_id": str},
    )
    bootstrap_membership = pd.read_csv(
        membership_dir / "yoshino_bootstrap_multiplicities_v01.csv.gz",
        dtype={"replicate_id": str, "unit_id": str},
    )

    subset_result = subset_stability(
        primary,
        subset_membership,
        chronology_method="mean",
        min_depth_fraction=1.0,
        min_overlap=min_overlap,
    )
    split_result = split_half_stability(
        primary,
        split_membership,
        chronology_method="mean",
        min_depth_fraction=1.0,
        min_overlap=min_overlap,
        require_exhaustive=True,
        require_balanced=True,
    )
    bootstrap_result = tree_bootstrap(
        primary,
        bootstrap_membership,
        chronology_method="mean",
        min_depth_fraction=1.0,
        min_overlap=min_overlap,
        require_full_size=True,
    )
    loto_result = leave_one_tree_out(
        primary,
        chronology_method="mean",
        min_depth_fraction=1.0,
        min_overlap=min_overlap,
    )

    raw_tables = {
        "subset": subset_result,
        "split_half": split_result,
        "bootstrap": bootstrap_result,
        "loto": loto_result,
    }
    output_paths: list[Path] = []
    for label, frame in raw_tables.items():
        path = output_dir / f"yoshino_{label}_results_v01.csv.gz"
        helper.write_csv_gzip_reproducible(frame, path)
        output_paths.append(path)

    subset_summary = helper.summarize_subset("yoshino", subset_result)
    split_summary = helper.summarize_split("yoshino", split_result)
    bootstrap_summary = helper.summarize_bootstrap("yoshino", bootstrap_result)
    loto_summary = helper.summarize_loto("yoshino", loto_result)
    structure = next(
        row
        for row in structure_rows
        if row["detrending_method"] == helper.PRIMARY_METHOD
    )

    n_trees = int(primary["tree_id"].nunique())
    half_size = n_trees // 2
    half_row = subset_summary.loc[
        subset_summary["subset_size"] == half_size
    ].iloc[0]
    split_row = split_summary.iloc[0]
    boot_row = bootstrap_summary.iloc[0]
    loto_row = loto_summary.iloc[0]
    key_results = pd.DataFrame(
        [
            {
                "dataset": "yoshino",
                "n_trees": n_trees,
                "analysis_start": int(contract["analysis_start"]),
                "analysis_end": int(contract["analysis_end"]),
                "analysis_years": int(contract["analysis_years"]),
                "rbar": structure["rbar"],
                "eps": structure["eps"],
                "effective_dimension": structure["effective_dimension"],
                "half_subset_size": half_size,
                "half_subset_correlation_q05": half_row["correlation_q05"],
                "half_subset_correlation_median": half_row["correlation_median"],
                "half_subset_correlation_min": half_row["correlation_min"],
                "split_correlation_q05": split_row["correlation_q05"],
                "split_correlation_median": split_row["correlation_median"],
                "split_correlation_min": split_row["correlation_min"],
                "bootstrap_srms_q95": boot_row[
                    "standardized_rms_difference_q95"
                ],
                "bootstrap_variance_ratio_q05": boot_row["variance_ratio_q05"],
                "bootstrap_variance_ratio_q95": boot_row["variance_ratio_q95"],
                "loto_max_srms": loto_row[
                    "maximum_standardized_rms_difference"
                ],
                "loto_max_abs_delta_eps": loto_row[
                    "maximum_absolute_delta_eps"
                ],
                "loto_max_abs_delta_effective_dimension": loto_row[
                    "maximum_absolute_delta_effective_dimension"
                ],
            }
        ]
    )

    selected_map = keymap.loc[
        keymap["include_primary_cohort"].astype(str).str.casefold().eq("true")
    ]
    source_groups = (
        selected_map.groupby(["source_group", "site_id"], sort=True)
        .size()
        .rename("n_specimens")
        .reset_index()
    )

    compact = {
        "structure_summary_v01.csv": pd.DataFrame(structure_rows),
        "detrending_sensitivity_summary_v01.csv": sensitivity,
        "subset_summary_v01.csv": subset_summary,
        "split_half_summary_v01.csv": split_summary,
        "bootstrap_summary_v01.csv": bootstrap_summary,
        "loto_summary_v01.csv": loto_summary,
        "yoshino_key_results_v01.csv": key_results,
        "primary_cohort_source_groups_v01.csv": source_groups,
    }
    for name, frame in compact.items():
        path = output_dir / name
        helper.write_csv(frame, path)
        output_paths.append(path)

    manifest = {
        "analysis": "yoshino_empirical_stage_c",
        "version": "v01",
        "input_level": "specimen_mean",
        "primary_detrending_method": helper.PRIMARY_METHOD,
        "chronology_method": "arithmetic mean",
        "min_depth_fraction": 1.0,
        "input_files": {
            "contract": {
                "path": stable_manifest_path(contract_path, repo_root),
                "sha256": file_sha256(contract_path),
            },
            "keymap": {
                "path": stable_manifest_path(keymap_path, repo_root),
                "sha256": file_sha256(keymap_path),
            },
            "tree_level": {
                "path": stable_manifest_path(tree_path, repo_root),
                "sha256": file_sha256(tree_path),
            },
            "membership_manifest": {
                "path": stable_manifest_path(
                    membership_dir / "membership_manifest_v01.json",
                    repo_root,
                ),
                "sha256": file_sha256(
                    membership_dir / "membership_manifest_v01.json"
                ),
            },
        },
        "outputs": {
            path.name: {
                "path": stable_manifest_path(path, repo_root),
                "sha256": file_sha256(path),
            }
            for path in sorted(output_paths)
        },
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
    }
    (output_dir / "analysis_manifest_v01.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print("Yoshino empirical Stage C analysis completed.")
    print(key_results.to_string(index=False))
    print(source_groups.to_string(index=False))


if __name__ == "__main__":
    main()
