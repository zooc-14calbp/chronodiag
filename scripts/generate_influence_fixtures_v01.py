"""Regenerate the frozen tree-bootstrap and LOTO fixture set.

Run from the repository root after installing the Python package in editable mode:

    python scripts/generate_influence_fixtures_v01.py

The script is deterministic and contains no random-number generator. Expected
outputs are frozen R/Python parity tables; changes require scientific review rather
than automatic acceptance.
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

from chronodiag import leave_one_tree_out, tree_bootstrap  # noqa: E402

FIXTURE_DIR = REPO_ROOT / "shared_fixtures"


def build_example_data() -> pd.DataFrame:
    """Return six coherent series with one deliberately influential unit."""
    years = np.arange(2000, 2024)
    t = np.arange(len(years), dtype=float)
    base = np.sin(2 * np.pi * t / 7.0) + 0.03 * t
    values = {
        "T01": base + 0.10 * np.cos(2 * np.pi * t / 5.0),
        "T02": 0.90 * base - 0.08 * np.sin(2 * np.pi * t / 4.0),
        "T03": 1.05 * base + 0.06 * np.cos(2 * np.pi * t / 3.0),
        "T04": base + 0.12 * np.sin(2 * np.pi * t / 6.0),
        "T05": 0.85 * base + 0.15 * np.cos(2 * np.pi * t / 8.0),
        "T06": 1.25 * base + 0.45 * np.sin(2 * np.pi * t / 4.0) + 0.10,
    }
    values["T03"][5] = np.nan
    values["T05"][12] = np.nan
    values["T06"][18] = np.nan

    rows: list[dict[str, object]] = []
    for unit in sorted(values):
        rows.extend(
            {
                "year": int(year),
                "series_id": unit,
                "tree_id": unit,
                "value": float(value) if np.isfinite(value) else np.nan,
            }
            for year, value in zip(years, values[unit], strict=True)
        )
    return pd.DataFrame(rows)


def build_bootstrap_multiplicities() -> pd.DataFrame:
    definitions = {
        "B01": {"T01": 2, "T02": 1, "T03": 1, "T04": 1, "T05": 1},
        "B02": {unit: 1 for unit in ["T01", "T02", "T03", "T04", "T05", "T06"]},
        "B03": {"T01": 1, "T03": 1, "T05": 1, "T06": 3},
        "B04": {"T02": 2, "T04": 2, "T05": 2},
    }
    return pd.DataFrame(
        [
            {
                "replicate_id": replicate_id,
                "unit_id": unit_id,
                "multiplicity": multiplicity,
            }
            for replicate_id, units in definitions.items()
            for unit_id, multiplicity in units.items()
        ]
    )


def main() -> None:
    FIXTURE_DIR.mkdir(parents=True, exist_ok=True)
    data = build_example_data()
    multiplicities = build_bootstrap_multiplicities()

    data.to_csv(
        FIXTURE_DIR / "influence_example_long_v01.csv",
        index=False,
        float_format="%.15g",
    )
    multiplicities.to_csv(
        FIXTURE_DIR / "bootstrap_multiplicities_v01.csv", index=False
    )

    bootstrap_expected = tree_bootstrap(data, multiplicities)
    loto_expected = leave_one_tree_out(data)

    bootstrap_frozen = bootstrap_expected.copy()
    bootstrap_frozen["all_units_represented"] = bootstrap_frozen[
        "all_units_represented"
    ].map({True: "TRUE", False: "FALSE"})
    bootstrap_frozen.to_csv(
        FIXTURE_DIR / "expected_tree_bootstrap_v01.csv",
        index=False,
        float_format="%.15g",
    )
    loto_expected.to_csv(
        FIXTURE_DIR / "expected_leave_one_tree_out_v01.csv",
        index=False,
        float_format="%.15g",
    )

    print(f"Wrote influence fixtures to {FIXTURE_DIR}")
    print(bootstrap_expected.to_string(index=False))
    print(loto_expected.to_string(index=False))


if __name__ == "__main__":
    main()
