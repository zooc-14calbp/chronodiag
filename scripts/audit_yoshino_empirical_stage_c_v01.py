"""Audit the frozen Yoshino specimen-mean RWL and value-blind cohort."""
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
    candidates: list[dict[str, object]] = []
    starts = sorted(keymap["series_start"].astype(int).unique())
    ends = sorted(keymap["series_end"].astype(int).unique())
    for start in starts:
        for end in ends:
            if end < start:
                continue
            years = int(end - start + 1)
            covered = keymap.loc[
                (keymap["series_start"].astype(int) <= start)
                & (keymap["series_end"].astype(int) >= end)
            ].copy()
            n_trees = int(covered["tree_id"].nunique())
            if n_trees < 2:
                continue
            if years < minimum_years:
                continue
            if years < minimum_years_per_tree_ratio * n_trees:
                continue
            candidates.append(
                {
                    "analysis_start": int(start),
                    "analysis_end": int(end),
                    "analysis_years": years,
                    "n_trees": n_trees,
                    "tree_ids": "|".join(sorted(covered["tree_id"].astype(str).unique())),
                    "source_groups": "|".join(
                        sorted(covered["source_group"].astype(str).unique())
                    ),
                }
            )
    if not candidates:
        raise ValueError("no cohort satisfies the frozen design")
    table = pd.DataFrame(candidates).sort_values(
        ["n_trees", "analysis_years", "analysis_start", "analysis_end"],
        ascending=[False, False, True, True],
        kind="stable",
    ).reset_index(drop=True)
    return table.iloc[0].to_dict(), table


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    parser.add_argument("--raw-path", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()

    repo_root = args.repo_root.resolve()
    raw_path = (
        args.raw_path.resolve()
        if args.raw_path is not None
        else repo_root / "data" / "yoshino" / "raw" / "nr_all.rwl"
    )
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else repo_root / "analysis" / "yoshino_empirical_stage_c_v01"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

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

    if not raw_path.is_file():
        raise FileNotFoundError(f"Yoshino specimen-mean RWL not found: {raw_path}")

    raw = read_tucson_rwl(
        raw_path,
        precision="auto",
        drop_missing=False,
        preserve_zero=True,
        strict=True,
    )
    raw_ids = sorted(raw["series_id"].astype(str).unique())
    map_ids = sorted(keymap["series_id"].astype(str).unique())
    if raw_ids != map_ids:
        raise AssertionError("Yoshino raw RWL and key map contain different series IDs")
    if len(raw_ids) != int(contract["n_raw_series_expected"]):
        raise AssertionError("unexpected number of Yoshino specimen-mean series")
    if set(keymap["input_level"]) != {"specimen_mean"}:
        raise AssertionError("Yoshino key map must declare specimen-mean input")

    spans = (
        raw.loc[raw["value"].notna()]
        .groupby("series_id", sort=True)["year"]
        .agg(series_start="min", series_end="max", n_observations="size")
        .reset_index()
    )
    observed = keymap.merge(spans, on="series_id", suffixes=("_contract", "_observed"))
    if not (
        observed["series_start_contract"].astype(int)
        == observed["series_start_observed"].astype(int)
    ).all():
        raise AssertionError("one or more Yoshino start years differ from the key map")
    if not (
        observed["series_end_contract"].astype(int)
        == observed["series_end_observed"].astype(int)
    ).all():
        raise AssertionError("one or more Yoshino end years differ from the key map")

    best, candidates = select_cohort(
        keymap,
        minimum_years=int(contract["minimum_years"]),
        minimum_years_per_tree_ratio=float(
            contract["minimum_years_per_tree_ratio"]
        ),
    )
    selected = sorted(
        keymap.loc[
            keymap["include_primary_cohort"].astype(str).str.casefold().eq("true"),
            "tree_id",
        ].astype(str)
    )
    if int(best["n_trees"]) != int(contract["n_trees_expected"]):
        raise AssertionError("value-blind cohort tree count differs from contract")
    if int(best["analysis_start"]) != int(contract["analysis_start"]):
        raise AssertionError("value-blind cohort start differs from contract")
    if int(best["analysis_end"]) != int(contract["analysis_end"]):
        raise AssertionError("value-blind cohort end differs from contract")
    if str(best["tree_ids"]).split("|") != selected:
        raise AssertionError("key-map cohort flags do not match value-blind selection")

    source_audit = pd.DataFrame(
        [
            {
                "dataset": "yoshino",
                "raw_path": raw_path.as_posix(),
                "raw_sha256": file_sha256(raw_path),
                "raw_size_bytes": raw_path.stat().st_size,
                "n_specimen_mean_series": len(raw_ids),
                "first_year": int(raw["year"].min()),
                "last_year": int(raw["year"].max()),
                "n_observations_nonmissing": int(raw["value"].notna().sum()),
                "n_zero_width": int((raw["value"] == 0).sum()),
                "input_level": "specimen_mean",
                "status": "ok",
            }
        ]
    )
    cohort_audit = pd.DataFrame([best])
    cohort_audit.insert(0, "dataset", "yoshino")
    group_counts = (
        keymap.loc[
            keymap["include_primary_cohort"].astype(str).str.casefold().eq("true")
        ]
        .groupby(["source_group", "site_id"], sort=True)
        .size()
        .rename("n_specimens")
        .reset_index()
    )

    source_audit.to_csv(output_dir / "source_audit_v01.csv", index=False)
    observed.to_csv(output_dir / "series_span_audit_v01.csv", index=False)
    cohort_audit.to_csv(output_dir / "cohort_selection_audit_v01.csv", index=False)
    candidates.head(20).to_csv(
        output_dir / "cohort_candidates_top20_v01.csv", index=False
    )
    group_counts.to_csv(output_dir / "primary_cohort_source_groups_v01.csv", index=False)

    print("Yoshino Stage C source and cohort audit passed.")
    print(source_audit.to_string(index=False))
    print(cohort_audit.to_string(index=False))
    print(group_counts.to_string(index=False))


if __name__ == "__main__":
    main()
