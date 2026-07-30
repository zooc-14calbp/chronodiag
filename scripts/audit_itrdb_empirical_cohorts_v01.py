"""Audit and reproduce the frozen empirical ITRDB cohort-selection rule.

The audit uses raw Tucson RWL spans and explicit key maps.  A tree's strict
interval is the intersection of all core/radius spans assigned to that tree.
The selected cohort maximizes tree count, then interval length, subject to:

* interval length >= 50 years; and
* interval length >= 3 * number of trees.

These are analysis-design constraints for the present covariance experiment,
not universal chronology-quality thresholds.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
from pathlib import Path
import sys

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
PYTHON_SRC = REPO_ROOT / "python" / "src"
if str(PYTHON_SRC) not in sys.path:
    sys.path.insert(0, str(PYTHON_SRC))

from chronodiag import read_tucson_rwl  # noqa: E402


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def select_cohort(
    keymap: pd.DataFrame,
    *,
    minimum_years: int,
    minimum_years_per_tree_ratio: float,
) -> tuple[dict[str, object], pd.DataFrame]:
    tree_spans = (
        keymap.groupby("tree_id", sort=True, observed=True)
        .agg(
            strict_start=("series_start", "max"),
            strict_end=("series_end", "min"),
            n_series=("series_id", "nunique"),
        )
        .reset_index()
    )
    tree_spans["strict_length"] = (
        tree_spans["strict_end"] - tree_spans["strict_start"] + 1
    )

    candidates: list[dict[str, object]] = []
    for start in sorted(tree_spans["strict_start"].unique()):
        for end in sorted(tree_spans["strict_end"].unique()):
            length = int(end - start + 1)
            if length < minimum_years:
                continue
            eligible = tree_spans.loc[
                (tree_spans["strict_start"] <= start)
                & (tree_spans["strict_end"] >= end)
            ]
            n_trees = int(len(eligible))
            if n_trees < 2:
                continue
            if length < minimum_years_per_tree_ratio * n_trees:
                continue
            candidates.append(
                {
                    "analysis_start": int(start),
                    "analysis_end": int(end),
                    "analysis_years": length,
                    "n_trees": n_trees,
                    "tree_ids": "|".join(sorted(eligible["tree_id"].astype(str))),
                }
            )

    if not candidates:
        raise ValueError("no cohort satisfies the declared selection constraints")

    table = pd.DataFrame(candidates)
    table = table.sort_values(
        ["n_trees", "analysis_years", "analysis_start", "analysis_end"],
        ascending=[False, False, True, True],
        kind="stable",
    ).reset_index(drop=True)
    best = table.iloc[0].to_dict()
    return best, table


def audit_dataset(
    dataset: str,
    contract: pd.Series,
    *,
    project_root: Path,
    repo_root: Path,
) -> tuple[dict[str, object], pd.DataFrame]:
    raw_path = project_root / "data" / "itrdb_raw" / str(contract["raw_filename"])
    keymap_path = (
        repo_root / "data" / "itrdb" / "keymaps" / f"{dataset}_keymap_v01.csv"
    )
    if not raw_path.is_file():
        raise FileNotFoundError(f"raw RWL file not found: {raw_path}")
    if not keymap_path.is_file():
        raise FileNotFoundError(f"key map not found: {keymap_path}")

    observed_hash = file_sha256(raw_path)
    expected_hash = str(contract["raw_sha256"]).lower()
    if observed_hash != expected_hash:
        raise ValueError(
            f"{dataset}: raw SHA-256 mismatch; expected {expected_hash}, "
            f"observed {observed_hash}"
        )

    data = read_tucson_rwl(raw_path, preserve_zero=True)
    observed_spans = (
        data.groupby("series_id", sort=True, observed=True)["year"]
        .agg(series_start="min", series_end="max")
        .reset_index()
    )
    keymap = pd.read_csv(keymap_path, dtype={"series_id": str, "tree_id": str})
    required = {
        "series_id",
        "tree_id",
        "series_start",
        "series_end",
        "include_primary_cohort",
    }
    missing = required.difference(keymap.columns)
    if missing:
        raise ValueError(f"{dataset}: key map is missing columns {sorted(missing)}")

    observed_ids = set(observed_spans["series_id"].astype(str))
    mapped_ids = set(keymap["series_id"].astype(str))
    if observed_ids != mapped_ids:
        raise ValueError(
            f"{dataset}: raw/key-map ID mismatch; "
            f"raw_only={sorted(observed_ids - mapped_ids)}, "
            f"map_only={sorted(mapped_ids - observed_ids)}"
        )

    merged = keymap.drop(columns=["series_start", "series_end"]).merge(
        observed_spans,
        on="series_id",
        how="inner",
        validate="one_to_one",
    )
    best, candidates = select_cohort(
        merged,
        minimum_years=int(contract["minimum_interval_years"]),
        minimum_years_per_tree_ratio=float(
            contract["minimum_years_per_tree_ratio"]
        ),
    )

    selected_ids = set(str(best["tree_ids"]).split("|"))
    declared_ids = set(
        keymap.loc[
            keymap["include_primary_cohort"].astype(str).str.casefold().eq("true"),
            "tree_id",
        ].astype(str)
    )
    if selected_ids != declared_ids:
        raise ValueError(
            f"{dataset}: selected/declaration tree mismatch; "
            f"selected_only={sorted(selected_ids - declared_ids)}, "
            f"declared_only={sorted(declared_ids - selected_ids)}"
        )

    expected = {
        "analysis_start": int(contract["analysis_start"]),
        "analysis_end": int(contract["analysis_end"]),
        "analysis_years": int(contract["analysis_years"]),
        "n_trees": int(contract["n_trees_expected"]),
    }
    for key, value in expected.items():
        if int(best[key]) != value:
            raise ValueError(
                f"{dataset}: cohort mismatch for {key}; "
                f"expected {value}, observed {best[key]}"
            )

    selected_series = int(
        keymap.loc[
            keymap["tree_id"].astype(str).isin(selected_ids), "series_id"
        ].nunique()
    )
    if selected_series != int(contract["n_series_expected"]):
        raise ValueError(
            f"{dataset}: selected series count mismatch; expected "
            f"{int(contract['n_series_expected'])}, observed {selected_series}"
        )

    summary = {
        "dataset": dataset,
        "raw_filename": raw_path.name,
        "raw_sha256": observed_hash,
        "n_raw_series": int(data["series_id"].nunique()),
        "n_raw_observations": int(len(data)),
        "n_zero_width": int((data["value"] == 0).sum()),
        "selected_start": int(best["analysis_start"]),
        "selected_end": int(best["analysis_end"]),
        "selected_years": int(best["analysis_years"]),
        "selected_trees": int(best["n_trees"]),
        "selected_series": selected_series,
        "selected_tree_ids": str(best["tree_ids"]),
        "status": "ok",
    }
    candidates.insert(0, "dataset", dataset)
    return summary, candidates.head(20)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project-root",
        type=Path,
        required=True,
        help="Project root containing data/itrdb_raw.",
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=REPO_ROOT,
        help="chronodiag repository root.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory; defaults to analysis/itrdb_empirical_stage_b_v01.",
    )
    args = parser.parse_args()

    project_root = args.project_root.resolve()
    repo_root = args.repo_root.resolve()
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else repo_root / "analysis" / "itrdb_empirical_stage_b_v01"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    contract_path = repo_root / "data" / "itrdb" / "empirical_dataset_contract_v01.csv"
    contracts = pd.read_csv(contract_path)
    summaries: list[dict[str, object]] = []
    candidate_tables: list[pd.DataFrame] = []

    for row in contracts.itertuples(index=False):
        contract = pd.Series(row._asdict())
        summary, candidates = audit_dataset(
            str(contract["dataset"]),
            contract,
            project_root=project_root,
            repo_root=repo_root,
        )
        summaries.append(summary)
        candidate_tables.append(candidates)

    summary_table = pd.DataFrame(summaries)
    candidate_table = pd.concat(candidate_tables, ignore_index=True)
    summary_path = output_dir / "cohort_selection_audit_v01.csv"
    candidates_path = output_dir / "cohort_candidates_top20_v01.csv"
    summary_table.to_csv(summary_path, index=False)
    candidate_table.to_csv(candidates_path, index=False)

    print("ITRDB empirical cohort audit passed.")
    print(summary_table.to_string(index=False))
    print(f"Summary: {summary_path}")
    print(f"Top candidates: {candidates_path}")


if __name__ == "__main__":
    main()
