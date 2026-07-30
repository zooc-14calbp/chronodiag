"""Validate frozen Yoshino empirical Stage C outputs."""
from __future__ import annotations

import argparse
from hashlib import sha256
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
TOLERANCE = 1e-10


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def run_command(arguments: list[str]) -> None:
    completed = subprocess.run(arguments, text=True)
    if completed.returncode != 0:
        raise RuntimeError(
            f"validation command failed with exit code {completed.returncode}: "
            + " ".join(arguments)
        )


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


def compare_memberships(left_path: Path, right_path: Path, name: str) -> None:
    left = pd.read_csv(left_path, dtype=str, keep_default_na=False)
    right = pd.read_csv(right_path, dtype=str, keep_default_na=False)
    pd.testing.assert_frame_equal(
        left,
        right,
        check_dtype=False,
        check_exact=True,
        check_like=False,
        obj=name,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    args = parser.parse_args()
    repo_root = args.repo_root.resolve()
    analysis_root = repo_root / "analysis" / "yoshino_empirical_stage_c_v01"
    preprocessed = analysis_root / "preprocessed"
    memberships = analysis_root / "memberships"
    results = analysis_root / "results"

    raw_path = repo_root / "data" / "yoshino" / "raw" / "nr_all.rwl"
    raw_manifest = repo_root / "data" / "yoshino" / "yoshino_raw_manifest_v01.csv"
    contract = pd.read_csv(
        repo_root / "data" / "yoshino" / "yoshino_dataset_contract_v01.csv"
    ).iloc[0]
    keymap = pd.read_csv(
        repo_root / "data" / "yoshino" / "yoshino_keymap_v01.csv",
        dtype={"series_id": str, "tree_id": str},
    )
    manifest = pd.read_csv(raw_manifest)
    if file_sha256(raw_path) != str(manifest["sha256"].iloc[0]):
        raise AssertionError("Yoshino raw specimen-mean RWL hash changed")
    if len(keymap) != 60 or keymap["tree_id"].nunique() != 60:
        raise AssertionError("Yoshino key map must contain 60 specimen means")
    selected = keymap[
        keymap["include_primary_cohort"].astype(str).str.casefold().eq("true")
    ]
    if selected["tree_id"].nunique() != int(contract["n_trees_expected"]):
        raise AssertionError("Yoshino primary cohort count differs from contract")

    tree_path = preprocessed / "tree_level_all_methods_v01.csv"
    tree = pd.read_csv(tree_path, dtype={"tree_id": str})
    primary = tree.loc[
        (tree["dataset"] == "yoshino")
        & (tree["detrending_method"] == "AgeDepSpline")
    ]
    expected_rows = int(contract["n_trees_expected"]) * int(contract["analysis_years"])
    if len(primary) != expected_rows:
        raise AssertionError("Yoshino primary tree-level row count differs")
    if primary[["year", "tree_id"]].duplicated().any():
        raise AssertionError("duplicate Yoshino year-tree rows")
    grouped = primary.groupby("tree_id", sort=True)["value_z"]
    if not np.allclose(grouped.mean(), 0.0, atol=1e-12, rtol=0.0):
        raise AssertionError("Yoshino z-scored specimen means are not centered")
    if not np.allclose(grouped.std(ddof=1), 1.0, atol=1e-12, rtol=0.0):
        raise AssertionError("Yoshino z-scored specimen means lack unit SD")

    with tempfile.TemporaryDirectory(prefix="chronodiag-yoshino-stage-c-") as name:
        temp = Path(name)
        temp_memberships = temp / "memberships"
        temp_results = temp / "results"
        run_command(
            [
                sys.executable,
                str(repo_root / "scripts" / "generate_yoshino_empirical_memberships_v01.py"),
                "--repo-root",
                str(repo_root),
                "--preprocessed-dir",
                str(preprocessed),
                "--output-dir",
                str(temp_memberships),
            ]
        )
        run_command(
            [
                sys.executable,
                str(repo_root / "scripts" / "run_yoshino_empirical_stage_c_v01.py"),
                "--repo-root",
                str(repo_root),
                "--preprocessed-dir",
                str(preprocessed),
                "--membership-dir",
                str(temp_memberships),
                "--output-dir",
                str(temp_results),
            ]
        )

        for path in sorted(memberships.glob("*.csv.gz")):
            compare_memberships(
                path,
                temp_memberships / path.name,
                f"Yoshino membership {path.name}",
            )
        for path in sorted(results.glob("*.csv*")):
            if path.name.endswith(".json"):
                continue
            compare_frames(
                path,
                temp_results / path.name,
                f"Yoshino Stage C {path.name}",
            )

    key = pd.read_csv(results / "yoshino_key_results_v01.csv").iloc[0]
    if int(key["n_trees"]) != 50:
        raise AssertionError("Yoshino key result does not use 50 specimens")
    if int(key["analysis_start"]) != 1850 or int(key["analysis_end"]) != 2004:
        raise AssertionError("Yoshino key result interval differs from contract")
    if not (0 < float(key["eps"]) <= 1):
        raise AssertionError("Yoshino EPS falls outside (0,1]")
    if not (1 <= float(key["effective_dimension"]) <= 50):
        raise AssertionError("Yoshino effective dimension falls outside [1,N]")
    if float(key["split_correlation_min"]) > float(key["split_correlation_median"]):
        raise AssertionError("Yoshino split minimum exceeds median")

    print("Yoshino empirical Stage C validation passed.")


if __name__ == "__main__":
    main()
