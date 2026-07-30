"""Validate frozen ITRDB empirical Stage B preprocessing and outputs."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
TOLERANCE = 1e-10


def compare_frames(left_path: Path, right_path: Path, name: str) -> None:
    left = pd.read_csv(left_path)
    right = pd.read_csv(right_path)
    if list(left.columns) != list(right.columns):
        raise AssertionError(f"{name}: column mismatch")
    if len(left) != len(right):
        raise AssertionError(f"{name}: row-count mismatch")
    pd.testing.assert_frame_equal(
        left,
        right,
        check_dtype=False,
        check_exact=False,
        atol=TOLERANCE,
        rtol=0.0,
        check_like=False,
        obj=name,
    )


def compare_membership_frames(left_path: Path, right_path: Path, name: str) -> None:
    """Compare decoded membership tables rather than compressed byte streams.

    Gzip streams can differ across Python/zlib platforms even when their decoded
    CSV content is identical.  Scientific reproducibility is therefore defined
    by the ordered explicit memberships, identifiers, and multiplicities.
    Reading every column as text also preserves identifiers with leading zeros.
    """
    left = pd.read_csv(left_path, dtype=str, keep_default_na=False)
    right = pd.read_csv(right_path, dtype=str, keep_default_na=False)
    if list(left.columns) != list(right.columns):
        raise AssertionError(f"{name}: column mismatch")
    if len(left) != len(right):
        raise AssertionError(f"{name}: row-count mismatch")
    pd.testing.assert_frame_equal(
        left,
        right,
        check_dtype=False,
        check_exact=True,
        check_like=False,
        obj=name,
    )


def run_command(arguments: list[str]) -> None:
    completed = subprocess.run(arguments, text=True)
    if completed.returncode != 0:
        raise RuntimeError(
            f"validation command failed with exit code {completed.returncode}: "
            + " ".join(arguments)
        )


def validate_preprocessed(repo_root: Path) -> None:
    analysis_root = repo_root / "analysis" / "itrdb_empirical_stage_b_v01"
    preprocessed = analysis_root / "preprocessed"
    tree_path = preprocessed / "tree_level_all_methods_v01.csv"
    summary_path = preprocessed / "preprocessing_summary_v01.csv"
    contract_path = repo_root / "data" / "itrdb" / "empirical_dataset_contract_v01.csv"
    for path in (tree_path, summary_path, contract_path):
        if not path.is_file():
            raise FileNotFoundError(f"required empirical file is absent: {path}")

    tree = pd.read_csv(tree_path, dtype={"tree_id": str})
    summary = pd.read_csv(summary_path)
    contracts = pd.read_csv(contract_path)
    primary_summary = summary.loc[
        (summary["detrending_method"] == "AgeDepSpline")
        & (summary["status"] == "ok")
    ]
    if set(primary_summary["dataset"]) != set(contracts["dataset"]):
        raise AssertionError("not every dataset has a successful primary preprocessing")

    for contract in contracts.itertuples(index=False):
        block = tree.loc[
            (tree["dataset"] == contract.dataset)
            & (tree["detrending_method"] == "AgeDepSpline")
        ]
        expected_rows = int(contract.n_trees_expected) * int(contract.analysis_years)
        if len(block) != expected_rows:
            raise AssertionError(
                f"{contract.dataset}: expected {expected_rows} primary tree rows, "
                f"observed {len(block)}"
            )
        if block[["year", "tree_id"]].duplicated().any():
            raise AssertionError(f"{contract.dataset}: duplicate year-tree rows")
        if block["value_z"].isna().any() or not np.isfinite(block["value_z"]).all():
            raise AssertionError(f"{contract.dataset}: invalid primary z values")
        if block["year"].min() != int(contract.analysis_start) or block[
            "year"
        ].max() != int(contract.analysis_end):
            raise AssertionError(f"{contract.dataset}: analysis interval mismatch")
        grouped = block.groupby("tree_id", sort=True)["value_z"]
        means = grouped.mean().to_numpy(dtype=float)
        sds = grouped.std(ddof=1).to_numpy(dtype=float)
        if not np.allclose(means, 0.0, atol=1e-12, rtol=0.0):
            raise AssertionError(f"{contract.dataset}: tree z means are not zero")
        if not np.allclose(sds, 1.0, atol=1e-12, rtol=0.0):
            raise AssertionError(f"{contract.dataset}: tree z SDs are not one")


def validate_outputs(repo_root: Path) -> None:
    analysis_root = repo_root / "analysis" / "itrdb_empirical_stage_b_v01"
    original_memberships = analysis_root / "memberships"
    original_results = analysis_root / "results"
    for path in (original_memberships, original_results):
        if not path.is_dir():
            raise FileNotFoundError(f"required empirical directory is absent: {path}")

    with tempfile.TemporaryDirectory(prefix="chronodiag-itrdb-stage-b-") as temp_name:
        temp = Path(temp_name)
        regenerated_memberships = temp / "memberships"
        regenerated_results = temp / "results"

        run_command(
            [
                sys.executable,
                str(repo_root / "scripts" / "generate_itrdb_empirical_memberships_v01.py"),
                "--repo-root",
                str(repo_root),
                "--preprocessed-dir",
                str(analysis_root / "preprocessed"),
                "--output-dir",
                str(regenerated_memberships),
            ]
        )
        run_command(
            [
                sys.executable,
                str(repo_root / "scripts" / "run_itrdb_empirical_stage_b_v01.py"),
                "--repo-root",
                str(repo_root),
                "--preprocessed-dir",
                str(analysis_root / "preprocessed"),
                "--membership-dir",
                str(regenerated_memberships),
                "--output-dir",
                str(regenerated_results),
            ]
        )

        membership_names = sorted(
            path.name
            for path in original_memberships.iterdir()
            if path.is_file() and path.name != "membership_manifest_v01.json"
        )
        for name in membership_names:
            compare_membership_frames(
                original_memberships / name,
                regenerated_memberships / name,
                f"ITRDB empirical Stage B membership: {name}",
            )

        result_names = sorted(
            path.name
            for path in original_results.iterdir()
            if path.is_file() and path.name != "analysis_manifest_v01.json"
        )
        for name in result_names:
            compare_frames(
                original_results / name,
                regenerated_results / name,
                f"ITRDB empirical Stage B: {name}",
            )

    key = pd.read_csv(original_results / "empirical_key_results_v01.csv")
    if not ((key["eps"] > 0) & (key["eps"] <= 1)).all():
        raise AssertionError("empirical EPS values fall outside (0, 1]")
    if not (
        (key["effective_dimension"] >= 1)
        & (key["effective_dimension"] <= key["n_trees"])
    ).all():
        raise AssertionError("empirical effective dimension falls outside [1, N]")
    if not (
        key["split_correlation_min"] <= key["split_correlation_median"]
    ).all():
        raise AssertionError("split-half minimum exceeds median")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    args = parser.parse_args()
    repo_root = args.repo_root.resolve()

    validate_preprocessed(repo_root)
    validate_outputs(repo_root)
    print("ITRDB empirical Stage B validation passed.")


if __name__ == "__main__":
    main()
