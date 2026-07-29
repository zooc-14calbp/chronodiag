"""Validate the chronodiag Python baseline and frozen ITRDB inputs.

Run from any directory after setting PYTHONPATH to ``python/src`` or installing
the package. The script verifies raw-file hashes, recomputes the six-file audit,
compares it with the frozen audit table, and writes machine-readable results under
``validation_output``.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from chronodiag import read_tucson_rwl

REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = REPO_ROOT / "validation_output"

FILES = [
    "ca533.rwl",
    "germ011.rwl",
    "germ012.rwl",
    "germ6.rwl",
    "fran003.rwl",
    "brit042.rwl",
]

AUDIT_COLUMNS = [
    "dataset",
    "n_series",
    "first_year",
    "last_year",
    "n_observations",
    "n_zero",
    "median_series_length",
    "min_series_length",
    "max_series_length",
    "max_depth",
    "full_common_start",
    "full_common_end",
    "full_common_length",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


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


def normalize_audit(table: pd.DataFrame) -> pd.DataFrame:
    result = table.loc[:, AUDIT_COLUMNS].copy()
    result = result.sort_values("dataset", kind="stable").reset_index(drop=True)
    for column in ("full_common_start", "full_common_end"):
        result[column] = pd.to_numeric(result[column], errors="coerce")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project-root",
        type=Path,
        default=REPO_ROOT.parents[1],
        help="EPS project root containing data/ and software/ (default: inferred from repository location)",
    )
    args = parser.parse_args()
    project_root = args.project_root.resolve()
    raw_dir = project_root / "data" / "itrdb_raw"
    metadata_dir = project_root / "data" / "itrdb_metadata"
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    manifest_path = metadata_dir / "itrdb_raw_file_manifest_v01.csv"
    frozen_path = metadata_dir / "itrdb_raw_audit_v01.csv"
    if not manifest_path.exists() or not frozen_path.exists():
        raise FileNotFoundError("ITRDB manifest or frozen audit table is missing")

    manifest = pd.read_csv(manifest_path).set_index("filename")
    hash_results: dict[str, str] = {}
    for name in FILES:
        path = raw_dir / name
        if not path.exists():
            raise FileNotFoundError(f"missing raw RWL file: {path}")
        observed = sha256(path)
        expected = str(manifest.loc[name, "sha256"]).lower()
        if observed != expected:
            raise ValueError(
                f"SHA-256 mismatch for {name}: expected {expected}, observed {observed}"
            )
        hash_results[name] = observed

    observed_table = pd.DataFrame(audit_one(raw_dir / name) for name in FILES)
    frozen_table = pd.read_csv(frozen_path)
    observed_normalized = normalize_audit(observed_table)
    frozen_normalized = normalize_audit(frozen_table)

    try:
        pd.testing.assert_frame_equal(
            observed_normalized,
            frozen_normalized,
            check_exact=False,
            atol=1e-12,
            rtol=0.0,
            check_dtype=False,
        )
    except AssertionError as exc:
        raise AssertionError("recomputed ITRDB audit differs from frozen v01") from exc

    audit_output = OUTPUT_DIR / "itrdb_raw_audit_local_v02.csv"
    summary_output = OUTPUT_DIR / "pre_github_validation_summary_v02.json"
    observed_table.to_csv(audit_output, index=False)
    summary = {
        "status": "passed",
        "raw_files_verified": len(hash_results),
        "audit_matches_frozen_v01": True,
        "audit_output": str(audit_output),
        "hashes": hash_results,
    }
    summary_output.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("All six raw RWL hashes are verified.")
    print("The recomputed audit matches itrdb_raw_audit_v01.csv.")
    print(observed_table.to_string(index=False))
    print(f"Summary: {summary_output}")


if __name__ == "__main__":
    main()
