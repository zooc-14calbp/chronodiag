"""Generate auditable empirical subset, split-half, and bootstrap memberships."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import combinations
import gzip
import io
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
MASTER_SEED = 20260730
MAX_SUBSETS_PER_SIZE = 250
SPLIT_REPLICATES = 1000
BOOTSTRAP_REPLICATES = 2000
VERSION = "v01"


def stable_seed(*parts: object) -> int:
    payload = "|".join(str(part) for part in (MASTER_SEED, *parts)).encode("utf-8")
    return int.from_bytes(sha256(payload).digest()[:8], "big", signed=False)


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def stable_manifest_path(path: Path, repo_root: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(repo_root.resolve()).as_posix()
    except ValueError:
        return resolved.as_posix()


def write_csv_gzip_reproducible(frame: pd.DataFrame, path: Path) -> None:
    text = frame.to_csv(index=False, lineterminator="\n")
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            zipped.write(text.encode("utf-8"))


def sample_unique_combinations(
    units: tuple[str, ...],
    size: int,
    maximum: int,
    *,
    design_key: str,
) -> list[tuple[str, ...]]:
    total = math.comb(len(units), size)
    if total <= maximum:
        return list(combinations(units, size))

    rng = np.random.default_rng(stable_seed(design_key, size))
    seen: set[tuple[str, ...]] = set()
    indices = np.arange(len(units))
    while len(seen) < maximum:
        selected = tuple(
            sorted(
                (units[int(index)] for index in rng.choice(indices, size=size, replace=False))
            )
        )
        seen.add(selected)
    return sorted(seen)


def generate_subset_table(dataset: str, units: tuple[str, ...]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for subset_size in range(2, len(units) - 1):
        selected = sample_unique_combinations(
            units,
            subset_size,
            MAX_SUBSETS_PER_SIZE,
            design_key=f"{dataset}|subset",
        )
        for replicate_number, subset in enumerate(selected, start=1):
            replicate_id = f"{dataset}_S{subset_size:02d}_{replicate_number:04d}"
            for unit_id in subset:
                rows.append(
                    {
                        "replicate_id": replicate_id,
                        "unit_id": unit_id,
                        "subset_size": subset_size,
                    }
                )
    return pd.DataFrame(rows)


def generate_split_table(dataset: str, units: tuple[str, ...]) -> pd.DataFrame:
    size_a = len(units) // 2
    group_as = sample_unique_combinations(
        units,
        size_a,
        SPLIT_REPLICATES,
        design_key=f"{dataset}|split",
    )
    full = set(units)
    rows: list[dict[str, object]] = []
    for replicate_number, group_a in enumerate(group_as, start=1):
        set_a = set(group_a)
        group_b = tuple(sorted(full.difference(set_a)))
        replicate_id = f"{dataset}_H_{replicate_number:04d}"
        rows.extend(
            {
                "replicate_id": replicate_id,
                "group": "A",
                "unit_id": unit_id,
                "group_size": len(group_a),
            }
            for unit_id in group_a
        )
        rows.extend(
            {
                "replicate_id": replicate_id,
                "group": "B",
                "unit_id": unit_id,
                "group_size": len(group_b),
            }
            for unit_id in group_b
        )
    return pd.DataFrame(rows)


def generate_bootstrap_table(dataset: str, units: tuple[str, ...]) -> pd.DataFrame:
    rng = np.random.default_rng(stable_seed(dataset, "bootstrap"))
    rows: list[dict[str, object]] = []
    unit_array = np.asarray(units, dtype=object)
    for replicate_number in range(1, BOOTSTRAP_REPLICATES + 1):
        draws = rng.choice(unit_array, size=len(units), replace=True)
        unique, counts = np.unique(draws, return_counts=True)
        replicate_id = f"{dataset}_B_{replicate_number:04d}"
        for unit_id, multiplicity in zip(unique.tolist(), counts.tolist()):
            rows.append(
                {
                    "replicate_id": replicate_id,
                    "unit_id": str(unit_id),
                    "multiplicity": int(multiplicity),
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    parser.add_argument("--preprocessed-dir", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()

    repo_root = args.repo_root.resolve()
    analysis_root = repo_root / "analysis" / "itrdb_empirical_stage_b_v01"
    preprocessed_dir = (
        args.preprocessed_dir.resolve()
        if args.preprocessed_dir is not None
        else analysis_root / "preprocessed"
    )
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else analysis_root / "memberships"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    tree_path = preprocessed_dir / "tree_level_all_methods_v01.csv"
    if not tree_path.is_file():
        raise FileNotFoundError(
            f"preprocessed tree-level table not found: {tree_path}; run the R "
            "preprocessing script first"
        )
    tree_data = pd.read_csv(tree_path, dtype={"tree_id": str})
    required = {"dataset", "detrending_method", "year", "tree_id", "value_z"}
    missing = required.difference(tree_data.columns)
    if missing:
        raise ValueError(f"preprocessed table is missing columns {sorted(missing)}")

    primary = tree_data.loc[tree_data["detrending_method"] == "AgeDepSpline"].copy()
    if primary.empty:
        raise ValueError("no AgeDepSpline primary data were found")

    manifest: dict[str, object] = {
        "analysis": "itrdb_empirical_stage_b_memberships",
        "version": VERSION,
        "master_seed": MASTER_SEED,
        "max_subsets_per_size": MAX_SUBSETS_PER_SIZE,
        "split_replicates": SPLIT_REPLICATES,
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "datasets": {},
        "files": {},
    }
    design_rows: list[dict[str, object]] = []

    for dataset, block in primary.groupby("dataset", sort=True, observed=True):
        units = tuple(sorted(block["tree_id"].astype(str).unique()))
        n_trees = len(units)
        if n_trees < 5:
            raise ValueError(f"{dataset}: at least five trees are required")

        expected_rows = len(units) * block["year"].nunique()
        if len(block) != expected_rows:
            raise ValueError(
                f"{dataset}: primary tree-level data are not a complete year-tree grid"
            )

        subset = generate_subset_table(str(dataset), units)
        split = generate_split_table(str(dataset), units)
        bootstrap = generate_bootstrap_table(str(dataset), units)

        paths = {
            "subset": output_dir / f"{dataset}_subset_memberships_v01.csv.gz",
            "split": output_dir / f"{dataset}_split_half_memberships_v01.csv.gz",
            "bootstrap": output_dir / f"{dataset}_bootstrap_multiplicities_v01.csv.gz",
        }
        write_csv_gzip_reproducible(subset, paths["subset"])
        write_csv_gzip_reproducible(split, paths["split"])
        write_csv_gzip_reproducible(bootstrap, paths["bootstrap"])

        subset_counts = (
            subset[["replicate_id", "subset_size"]]
            .drop_duplicates()
            .groupby("subset_size", sort=True)
            .size()
        )
        design_rows.append(
            {
                "dataset": dataset,
                "n_trees": n_trees,
                "subset_min_size": int(subset_counts.index.min()),
                "subset_max_size": int(subset_counts.index.max()),
                "subset_replicates_total": int(subset["replicate_id"].nunique()),
                "subset_replicates_min_per_size": int(subset_counts.min()),
                "subset_replicates_max_per_size": int(subset_counts.max()),
                "split_group_a_size": int(n_trees // 2),
                "split_group_b_size": int(n_trees - n_trees // 2),
                "split_replicates": int(split["replicate_id"].nunique()),
                "bootstrap_replicates": int(bootstrap["replicate_id"].nunique()),
            }
        )

        manifest["datasets"][str(dataset)] = {
            "n_trees": n_trees,
            "tree_ids": list(units),
            "subset_sizes": [int(value) for value in subset_counts.index],
            "dataset_seeds": {
                "bootstrap": stable_seed(dataset, "bootstrap"),
                "split": stable_seed(f"{dataset}|split", n_trees // 2),
            },
        }
        for label, path in paths.items():
            relative = stable_manifest_path(path, repo_root)
            manifest["files"][f"{dataset}_{label}"] = {
                "path": relative,
                "sha256": file_sha256(path),
                "rows": int(
                    {"subset": len(subset), "split": len(split), "bootstrap": len(bootstrap)}[
                        label
                    ]
                ),
            }

    design_table = pd.DataFrame(design_rows)
    design_path = output_dir / "membership_design_summary_v01.csv"
    design_table.to_csv(design_path, index=False, lineterminator="\n")
    manifest["files"]["design_summary"] = {
        "path": stable_manifest_path(design_path, repo_root),
        "sha256": file_sha256(design_path),
        "rows": int(len(design_table)),
    }
    manifest_path = output_dir / "membership_manifest_v01.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    print("ITRDB empirical membership generation completed.")
    print(design_table.to_string(index=False))
    print(f"Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
