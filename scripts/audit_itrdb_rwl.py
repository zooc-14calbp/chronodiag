"""Audit Tucson RWL files with the public chronodiag input contract.

This script is intentionally separate from the package test suite because the
full ITRDB files are acquired from their dataset records rather than duplicated
as small unit-test fixtures.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from chronodiag import read_tucson_rwl


def audit_one(path: Path) -> dict[str, object]:
    data = read_tucson_rwl(path)
    spans = data.groupby("series_id")["year"].agg(["min", "max"])
    lengths = data.groupby("series_id").size()
    depth = data.groupby("year")["series_id"].nunique()
    full_start = int(spans["min"].max())
    full_end = int(spans["max"].min())
    common_length = max(0, full_end - full_start + 1)
    return {
        "dataset": path.stem.lower(),
        "filename": path.name,
        "n_series": int(data["series_id"].nunique()),
        "first_year": int(data["year"].min()),
        "last_year": int(data["year"].max()),
        "n_observations": int(len(data)),
        "n_zero": int((data["value"] == 0).sum()),
        "median_series_length": float(lengths.median()),
        "min_series_length": int(lengths.min()),
        "max_series_length": int(lengths.max()),
        "max_depth": int(depth.max()),
        "full_common_start": full_start if common_length else None,
        "full_common_end": full_end if common_length else None,
        "full_common_length": common_length,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="+", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    table = pd.DataFrame(audit_one(path) for path in args.files)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output, index=False)
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
