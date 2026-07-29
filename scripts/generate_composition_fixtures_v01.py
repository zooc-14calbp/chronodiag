"""Regenerate the frozen composition-diagnostic fixture set.

Run from the repository root after installing the Python package in editable mode:

    python scripts/generate_composition_fixtures_v01.py

The script is deterministic and contains no random-number generator. Expected
outputs are frozen reference tables used for R/Python parity; changes to those
outputs require scientific review rather than automatic acceptance.
"""
from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
PYTHON_SRC = REPO_ROOT / "python" / "src"
if str(PYTHON_SRC) not in sys.path:
    sys.path.insert(0, str(PYTHON_SRC))

from chronodiag import split_half_stability, subset_stability  # noqa: E402

FIXTURE_DIR = REPO_ROOT / "shared_fixtures"


def build_example_data() -> pd.DataFrame:
    years = np.arange(2000, 2024)
    t = np.arange(len(years), dtype=float)
    signal = (
        np.sin(2 * np.pi * t / 7.0)
        + 0.35 * np.cos(2 * np.pi * t / 13.0)
        + 0.02 * t
    )
    secondary_signal = np.cos(2 * np.pi * t / 5.0)
    loadings = [1.00, 0.95, 1.05, 0.90, 1.00, 0.85, 1.15, 0.75]
    secondary_loadings = [0.00, 0.05, -0.05, 0.15, -0.15, 0.25, -0.25, 0.35]

    rows: list[dict[str, object]] = []
    for index, (loading, secondary_loading) in enumerate(
        zip(loadings, secondary_loadings, strict=True), start=1
    ):
        residual = (
            0.12 * np.sin(2 * np.pi * (t + index - 1) / 3.0 + (index - 1) * 0.2)
            + 0.05 * np.cos(2 * np.pi * t * index / 11.0)
        )
        values = loading * signal + secondary_loading * secondary_signal + residual
        unit = f"T{index:02d}"
        rows.extend(
            {
                "year": int(year),
                "series_id": unit,
                "tree_id": unit,
                "value": float(value),
            }
            for year, value in zip(years, values, strict=True)
        )

    data = pd.DataFrame(rows)
    for unit, year in [("T03", 2005), ("T07", 2011), ("T08", 2018)]:
        data.loc[(data["tree_id"] == unit) & (data["year"] == year), "value"] = np.nan
    return data


def build_exchangeable_exact_data() -> pd.DataFrame:
    """Return a fixed N=6, rho=0.3 dataset with exact sample correlation."""
    n_years = 20
    n_units = 6
    rho = 0.3
    x = np.arange(1, n_years + 1, dtype=float)
    seed_matrix = np.column_stack(
        [
            np.sin(x * (index + 1) * 0.37)
            + 0.1 * np.cos(x * (index + 2) * 0.19)
            + 0.01 * x ** ((index % 3) + 1)
            for index in range(n_units)
        ]
    )
    seed_matrix -= seed_matrix.mean(axis=0)
    orthonormal, _ = np.linalg.qr(seed_matrix)
    for column in range(n_units):
        pivot = int(np.argmax(np.abs(orthonormal[:, column])))
        if orthonormal[pivot, column] < 0:
            orthonormal[:, column] *= -1

    correlation = (1 - rho) * np.eye(n_units) + rho * np.ones((n_units, n_units))
    values = np.sqrt(n_years - 1) * orthonormal @ np.linalg.cholesky(correlation).T
    observed = np.corrcoef(values, rowvar=False)
    if not np.allclose(observed, correlation, atol=1e-14, rtol=0.0):
        raise RuntimeError("exchangeable fixture construction failed")

    rows: list[dict[str, object]] = []
    for index in range(n_units):
        unit = f"E{index + 1:02d}"
        rows.extend(
            {
                "year": year,
                "series_id": unit,
                "tree_id": unit,
                "value": float(value),
            }
            for year, value in zip(
                range(1900, 1900 + n_years), values[:, index], strict=True
            )
        )
    return pd.DataFrame(rows)


def build_subsets() -> pd.DataFrame:
    definitions = {
        "S02_01": ["T01", "T02"],
        "S02_02": ["T03", "T04"],
        "S02_03": ["T07", "T08"],
        "S04_01": ["T01", "T03", "T05", "T07"],
        "S04_02": ["T02", "T04", "T06", "T08"],
        "S06_01": ["T01", "T02", "T03", "T04", "T05", "T06"],
        "S06_02": ["T03", "T04", "T05", "T06", "T07", "T08"],
    }
    return pd.DataFrame(
        [
            {"replicate_id": replicate_id, "unit_id": unit}
            for replicate_id, units in definitions.items()
            for unit in units
        ]
    )


def build_split_halves() -> pd.DataFrame:
    definitions = {
        "H01": (
            ["T01", "T02", "T03", "T04"],
            ["T05", "T06", "T07", "T08"],
        ),
        "H02": (
            ["T01", "T03", "T05", "T07"],
            ["T02", "T04", "T06", "T08"],
        ),
        "H03": (
            ["T01", "T04", "T06", "T07"],
            ["T02", "T03", "T05", "T08"],
        ),
        "H04": (
            ["T01", "T02", "T07", "T08"],
            ["T03", "T04", "T05", "T06"],
        ),
    }
    rows: list[dict[str, str]] = []
    for replicate_id, (group_a, group_b) in definitions.items():
        rows.extend(
            {"replicate_id": replicate_id, "group": "A", "unit_id": unit}
            for unit in group_a
        )
        rows.extend(
            {"replicate_id": replicate_id, "group": "B", "unit_id": unit}
            for unit in group_b
        )
    return pd.DataFrame(rows)


def main() -> None:
    FIXTURE_DIR.mkdir(parents=True, exist_ok=True)
    data = build_example_data()
    exchangeable = build_exchangeable_exact_data()
    subsets = build_subsets()
    split_halves = build_split_halves()

    data.to_csv(
        FIXTURE_DIR / "composition_example_long_v01.csv",
        index=False,
        float_format="%.15g",
    )
    exchangeable.to_csv(
        FIXTURE_DIR / "composition_exchangeable_exact_v01.csv",
        index=False,
        float_format="%.17g",
    )
    subsets.to_csv(FIXTURE_DIR / "subset_partitions_v01.csv", index=False)
    split_halves.to_csv(FIXTURE_DIR / "split_half_partitions_v01.csv", index=False)

    subset_expected = subset_stability(
        data, subsets, min_overlap=10, min_depth_fraction=1.0
    )
    split_expected = split_half_stability(
        data, split_halves, min_overlap=10, min_depth_fraction=1.0
    )
    # Freeze logical values in R-compatible uppercase form. pandas reads these
    # back as bool, while base R reads them as logical rather than character.
    for table, filename in [
        (subset_expected, "expected_subset_stability_v01.csv"),
        (split_expected, "expected_split_half_stability_v01.csv"),
    ]:
        frozen = table.copy()
        frozen["comparison_exhaustive"] = frozen["comparison_exhaustive"].map(
            {True: "TRUE", False: "FALSE"}
        )
        frozen.to_csv(
            FIXTURE_DIR / filename,
            index=False,
            float_format="%.15g",
        )

    print(f"Wrote composition fixtures to {FIXTURE_DIR}")
    print(subset_expected.to_string(index=False))
    print(split_expected.to_string(index=False))


if __name__ == "__main__":
    main()
