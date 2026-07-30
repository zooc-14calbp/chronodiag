"""Reproduce and validate the frozen same-EPS Stage A analysis.

This is the end-to-end scientific validation layer. It reruns Stage A in a
temporary directory, compares all covariance-level summaries and bootstrap/LOTO
outputs with the frozen repository tables, verifies the exact bootstrap
membership stream, and confirms that regenerated finite samples reproduce the
four target correlation matrices.

The exact long-form values are not required to be byte-identical across BLAS
implementations; their sample correlation matrices are the scientific contract.
"""
from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path
import sys
import tempfile

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
PYTHON_SRC = REPO_ROOT / "python" / "src"
if str(PYTHON_SRC) not in sys.path:
    sys.path.insert(0, str(PYTHON_SRC))

STAGE_A_PATH = REPO_ROOT / "scripts" / "run_same_eps_stage_a_v01.py"
SPEC = importlib.util.spec_from_file_location("same_eps_stage_a_v01", STAGE_A_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"could not load Stage A module: {STAGE_A_PATH}")
STAGE_A = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = STAGE_A
SPEC.loader.exec_module(STAGE_A)

FROZEN_DIR = REPO_ROOT / "analysis" / "same_eps_stage_a_v01"
FIXTURE_DIR = REPO_ROOT / "shared_fixtures"

COMPARISON_FILES = {
    "scenario_summary_v01.csv": None,
    "subset_summary_v01.csv": None,
    "split_half_summary_v01.csv": None,
    "bootstrap_multiplicities_v01.csv.gz": None,
    "bootstrap_results_v01.csv.gz": None,
    "bootstrap_summary_v01.csv": None,
    "bootstrap_driver_summary_v01.csv": None,
    "loto_results_v01.csv": None,
    "loto_summary_v01.csv": None,
    "key_results_v01.csv": None,
}


def _read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path)


def _compare_frames(
    regenerated: pd.DataFrame,
    frozen: pd.DataFrame,
    *,
    name: str,
    tolerance: float,
) -> None:
    pd.testing.assert_frame_equal(
        regenerated,
        frozen,
        check_dtype=False,
        check_exact=False,
        rtol=0.0,
        atol=tolerance,
        check_like=False,
        obj=name,
    )


def validate(
    *,
    tolerance: float = 1e-11,
    bootstrap_replicates: int = STAGE_A.DEFAULT_BOOTSTRAP_REPLICATES,
) -> None:
    if bootstrap_replicates != STAGE_A.DEFAULT_BOOTSTRAP_REPLICATES:
        raise ValueError(
            "the frozen Stage A contract uses exactly "
            f"{STAGE_A.DEFAULT_BOOTSTRAP_REPLICATES} bootstrap replicates"
        )

    with tempfile.TemporaryDirectory(prefix="chronodiag_stage_a_") as directory:
        regenerated_dir = Path(directory)
        STAGE_A.run_stage_a(
            fixture_dir=FIXTURE_DIR,
            output_dir=regenerated_dir,
            n_years=STAGE_A.DEFAULT_N_YEARS,
            start_year=STAGE_A.DEFAULT_START_YEAR,
            bootstrap_replicates=bootstrap_replicates,
            bootstrap_seed=STAGE_A.DEFAULT_BOOTSTRAP_SEED,
        )

        for filename in COMPARISON_FILES:
            regenerated = _read(regenerated_dir / filename)
            frozen = _read(FROZEN_DIR / filename)
            _compare_frames(
                regenerated,
                frozen,
                name=f"same-EPS Stage A: {filename}",
                tolerance=tolerance,
            )

        exact = _read(regenerated_dir / "exact_same_eps_long_v01.csv.gz")
        matrices = STAGE_A.load_scenario_matrices(FIXTURE_DIR)
        expected_rows = (
            len(STAGE_A.SCENARIOS)
            * STAGE_A.DEFAULT_N_YEARS
            * next(iter(matrices.values())).shape[0]
        )
        if len(exact) != expected_rows:
            raise AssertionError(
                f"regenerated exact sample has {len(exact)} rows; "
                f"expected {expected_rows}"
            )

        for scenario, target in matrices.items():
            block = exact.loc[exact["scenario"] == scenario]
            wide = block.pivot(index="year", columns="tree_id", values="value")
            observed = wide.corr().to_numpy(dtype=float)
            maximum = float(np.max(np.abs(observed - target)))
            if maximum > 5e-11:
                raise AssertionError(
                    f"{scenario} regenerated correlation differs from the "
                    f"target by {maximum}"
                )

    print("same-EPS Stage A end-to-end reproduction passed.")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tolerance", type=float, default=1e-11)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    validate(tolerance=args.tolerance)


if __name__ == "__main__":
    main()
