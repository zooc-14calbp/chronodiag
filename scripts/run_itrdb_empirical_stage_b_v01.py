"""Run the empirical ITRDB Stage B composition-sensitivity analysis."""
from __future__ import annotations

import argparse
from hashlib import sha256
import gzip
import json
import math
from pathlib import Path
import platform
import sys
from typing import Iterable

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
PYTHON_SRC = REPO_ROOT / "python" / "src"
if str(PYTHON_SRC) not in sys.path:
    sys.path.insert(0, str(PYTHON_SRC))

from chronodiag import (  # noqa: E402
    effective_dimension,
    eps_value,
    leave_one_tree_out,
    split_half_stability,
    subset_stability,
    tree_bootstrap,
)

PRIMARY_METHOD = "AgeDepSpline"
TIE_TOLERANCE = 1e-12
VERSION = "v01"


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def stable_manifest_path(path: Path, repo_root: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(repo_root.resolve()).as_posix()
    except ValueError:
        return resolved.as_posix()


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, lineterminator="\n")


def write_csv_gzip_reproducible(frame: pd.DataFrame, path: Path) -> None:
    text = frame.to_csv(index=False, lineterminator="\n")
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            zipped.write(text.encode("utf-8"))


def mean_off_diagonal(correlation: np.ndarray) -> float:
    n = correlation.shape[0]
    return float((correlation.sum() - np.trace(correlation)) / (n * (n - 1)))


def distribution_summary(
    values: Iterable[float],
    *,
    prefix: str = "",
) -> dict[str, float | int]:
    array = np.asarray(list(values), dtype=float)
    array = array[np.isfinite(array)]
    if array.size == 0:
        return {
            f"{prefix}n": 0,
            f"{prefix}mean": math.nan,
            f"{prefix}sd": math.nan,
            f"{prefix}q01": math.nan,
            f"{prefix}q05": math.nan,
            f"{prefix}q25": math.nan,
            f"{prefix}median": math.nan,
            f"{prefix}q75": math.nan,
            f"{prefix}q95": math.nan,
            f"{prefix}q99": math.nan,
            f"{prefix}min": math.nan,
            f"{prefix}max": math.nan,
        }
    quantiles = np.quantile(array, [0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99])
    return {
        f"{prefix}n": int(array.size),
        f"{prefix}mean": float(array.mean()),
        f"{prefix}sd": float(array.std(ddof=1)) if array.size > 1 else 0.0,
        f"{prefix}q01": float(quantiles[0]),
        f"{prefix}q05": float(quantiles[1]),
        f"{prefix}q25": float(quantiles[2]),
        f"{prefix}median": float(quantiles[3]),
        f"{prefix}q75": float(quantiles[4]),
        f"{prefix}q95": float(quantiles[5]),
        f"{prefix}q99": float(quantiles[6]),
        f"{prefix}min": float(array.min()),
        f"{prefix}max": float(array.max()),
    }


def complete_matrix(block: pd.DataFrame, value_col: str = "value_z") -> pd.DataFrame:
    matrix = block.pivot(index="year", columns="tree_id", values=value_col)
    matrix = matrix.sort_index().reindex(sorted(matrix.columns), axis=1)
    if matrix.isna().any().any():
        raise ValueError("tree-level input is not a complete year-tree matrix")
    if (matrix.std(axis=0, ddof=1) <= 0).any():
        raise ValueError("one or more tree-level series are constant")
    return matrix


def structure_row(dataset: str, method: str, block: pd.DataFrame) -> dict[str, object]:
    matrix = complete_matrix(block)
    values = matrix.to_numpy(dtype=float)
    correlation = np.corrcoef(values, rowvar=False)
    eigenvalues = np.linalg.eigvalsh(correlation)[::-1]
    n_trees = correlation.shape[0]
    off_diagonal = correlation[np.triu_indices(n_trees, k=1)]
    rbar = mean_off_diagonal(correlation)
    eps = float(eps_value(rbar, n_trees))
    dimension = float(effective_dimension(correlation))
    chronology = values.mean(axis=1)
    return {
        "dataset": dataset,
        "detrending_method": method,
        "analysis_start": int(matrix.index.min()),
        "analysis_end": int(matrix.index.max()),
        "analysis_years": int(matrix.shape[0]),
        "n_trees": int(n_trees),
        "rbar": rbar,
        "eps": eps,
        "effective_dimension": dimension,
        "effective_dimension_fraction": dimension / n_trees,
        "leading_eigenvalue": float(eigenvalues[0]),
        "second_eigenvalue": float(eigenvalues[1]) if n_trees > 1 else math.nan,
        "leading_eigenvalue_fraction": float(eigenvalues[0] / n_trees),
        "second_eigenvalue_fraction": (
            float(eigenvalues[1] / n_trees) if n_trees > 1 else math.nan
        ),
        "off_diagonal_mean": float(off_diagonal.mean()),
        "off_diagonal_sd": (
            float(off_diagonal.std(ddof=1)) if off_diagonal.size > 1 else 0.0
        ),
        "off_diagonal_min": float(off_diagonal.min()),
        "off_diagonal_max": float(off_diagonal.max()),
        "chronology_variance": float(np.var(chronology, ddof=1)),
    }


def summarize_subset(dataset: str, result: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for subset_size, block in result.groupby("group_a_size", sort=True):
        valid = block.loc[block["status"] == "ok"]
        row: dict[str, object] = {
            "dataset": dataset,
            "subset_size": int(subset_size),
            "n_requested": int(len(block)),
            "n_valid": int(len(valid)),
            "n_invalid": int(len(block) - len(valid)),
        }
        row.update(distribution_summary(valid["correlation"], prefix="correlation_"))
        row.update(
            distribution_summary(
                valid["squared_correlation"], prefix="squared_correlation_"
            )
        )
        rows.append(row)
    return pd.DataFrame(rows)


def summarize_split(dataset: str, result: pd.DataFrame) -> pd.DataFrame:
    valid = result.loc[result["status"] == "ok"]
    row: dict[str, object] = {
        "dataset": dataset,
        "group_a_size": int(result["group_a_size"].iloc[0]),
        "group_b_size": int(result["group_b_size"].iloc[0]),
        "n_requested": int(len(result)),
        "n_valid": int(len(valid)),
        "n_invalid": int(len(result) - len(valid)),
    }
    row.update(distribution_summary(valid["correlation"], prefix="correlation_"))
    row.update(
        distribution_summary(valid["squared_correlation"], prefix="squared_correlation_")
    )
    return pd.DataFrame([row])


def summarize_bootstrap(dataset: str, result: pd.DataFrame) -> pd.DataFrame:
    valid = result.loc[result["status"].isin(["ok", "constant_chronology"])]
    row: dict[str, object] = {
        "dataset": dataset,
        "n_requested": int(len(result)),
        "n_valid_direct_metrics": int(len(valid)),
        "n_correlation_ok": int((result["status"] == "ok").sum()),
    }
    for metric in (
        "correlation",
        "standardized_rms_difference",
        "mean_absolute_difference",
        "max_absolute_difference",
        "variance_ratio",
        "unique_unit_count",
        "kish_effective_unit_count",
        "maximum_multiplicity",
    ):
        row.update(distribution_summary(valid[metric], prefix=f"{metric}_"))
    return pd.DataFrame([row])


def tied_ids(
    frame: pd.DataFrame,
    metric: str,
    *,
    absolute: bool = False,
) -> tuple[str, str, int, float]:
    values = frame[metric].to_numpy(dtype=float)
    transformed = np.abs(values) if absolute else values
    finite = np.isfinite(transformed)
    if not finite.any():
        return "", "", 0, math.nan
    extreme = float(np.max(transformed[finite]))
    tied = frame.loc[
        finite & (np.abs(transformed - extreme) <= TIE_TOLERANCE),
        "omitted_unit_id",
    ].astype(str)
    identifiers = sorted(tied.tolist())
    return identifiers[0], "|".join(identifiers), len(identifiers), extreme


def summarize_loto(dataset: str, result: pd.DataFrame) -> pd.DataFrame:
    direct = result.loc[result["status"].isin(["ok", "constant_chronology"])].copy()
    structural = result.loc[result["structure_status"] == "ok"].copy()
    rms_id, rms_ids, rms_count, rms_value = tied_ids(
        direct, "standardized_rms_difference"
    )
    eps_id, eps_ids, eps_count, eps_value_abs = tied_ids(
        structural, "delta_eps", absolute=True
    )
    dim_id, dim_ids, dim_count, dim_value_abs = tied_ids(
        structural, "delta_effective_dimension", absolute=True
    )
    row: dict[str, object] = {
        "dataset": dataset,
        "n_omissions": int(len(result)),
        "n_direct_valid": int(len(direct)),
        "n_structure_valid": int(len(structural)),
        "minimum_loto_correlation": (
            float(result.loc[result["status"] == "ok", "correlation"].min())
            if (result["status"] == "ok").any()
            else math.nan
        ),
        "maximum_standardized_rms_difference": rms_value,
        "max_rms_unit_id": rms_id,
        "max_rms_unit_ids": rms_ids,
        "max_rms_tie_count": rms_count,
        "maximum_absolute_delta_eps": eps_value_abs,
        "max_abs_delta_eps_unit_id": eps_id,
        "max_abs_delta_eps_unit_ids": eps_ids,
        "max_abs_delta_eps_tie_count": eps_count,
        "maximum_absolute_delta_effective_dimension": dim_value_abs,
        "max_abs_delta_effective_dimension_unit_id": dim_id,
        "max_abs_delta_effective_dimension_unit_ids": dim_ids,
        "max_abs_delta_effective_dimension_tie_count": dim_count,
    }
    for metric in (
        "standardized_rms_difference",
        "delta_rbar",
        "delta_eps",
        "delta_effective_dimension",
    ):
        source = direct if metric == "standardized_rms_difference" else structural
        row.update(distribution_summary(source[metric], prefix=f"{metric}_"))
    return pd.DataFrame([row])


def chronology_sensitivity(
    dataset: str, tree_data: pd.DataFrame, primary_method: str
) -> pd.DataFrame:
    methods = sorted(tree_data["detrending_method"].unique())
    primary_block = tree_data.loc[tree_data["detrending_method"] == primary_method]
    primary_matrix = complete_matrix(primary_block)
    primary_chronology = primary_matrix.mean(axis=1)
    primary_sd = float(primary_chronology.std(ddof=1))
    rows: list[dict[str, object]] = []
    for method in methods:
        block = tree_data.loc[tree_data["detrending_method"] == method]
        matrix = complete_matrix(block)
        if not matrix.index.equals(primary_matrix.index) or list(matrix.columns) != list(
            primary_matrix.columns
        ):
            raise ValueError(f"{dataset}: detrending methods use different grids")
        chronology = matrix.mean(axis=1)
        difference = chronology - primary_chronology
        rows.append(
            {
                "dataset": dataset,
                "primary_method": primary_method,
                "comparison_method": method,
                "chronology_correlation": float(
                    np.corrcoef(primary_chronology, chronology)[0, 1]
                ),
                "rms_difference": float(
                    np.sqrt(np.mean(np.square(difference.to_numpy(dtype=float))))
                ),
                "standardized_rms_difference": (
                    float(
                        np.sqrt(np.mean(np.square(difference.to_numpy(dtype=float))))
                        / primary_sd
                    )
                    if primary_sd > 0
                    else math.nan
                ),
                "variance_ratio": float(
                    chronology.var(ddof=1) / primary_chronology.var(ddof=1)
                ),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    parser.add_argument("--preprocessed-dir", type=Path, default=None)
    parser.add_argument("--membership-dir", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()

    repo_root = args.repo_root.resolve()
    analysis_root = repo_root / "analysis" / "itrdb_empirical_stage_b_v01"
    preprocessed_dir = (
        args.preprocessed_dir.resolve()
        if args.preprocessed_dir is not None
        else analysis_root / "preprocessed"
    )
    membership_dir = (
        args.membership_dir.resolve()
        if args.membership_dir is not None
        else analysis_root / "memberships"
    )
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else analysis_root / "results"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    contract_path = repo_root / "data" / "itrdb" / "empirical_dataset_contract_v01.csv"
    contracts = pd.read_csv(contract_path)
    tree_path = preprocessed_dir / "tree_level_all_methods_v01.csv"
    tree_data = pd.read_csv(tree_path, dtype={"tree_id": str})
    required = {
        "dataset",
        "detrending_method",
        "year",
        "tree_id",
        "value_z",
    }
    missing = required.difference(tree_data.columns)
    if missing:
        raise ValueError(f"tree-level table is missing columns {sorted(missing)}")

    structure_rows: list[dict[str, object]] = []
    sensitivity_tables: list[pd.DataFrame] = []
    subset_summaries: list[pd.DataFrame] = []
    split_summaries: list[pd.DataFrame] = []
    bootstrap_summaries: list[pd.DataFrame] = []
    loto_summaries: list[pd.DataFrame] = []
    key_rows: list[dict[str, object]] = []
    output_paths: list[Path] = []

    for contract in contracts.itertuples(index=False):
        dataset = str(contract.dataset)
        dataset_data = tree_data.loc[tree_data["dataset"] == dataset].copy()
        if dataset_data.empty:
            raise ValueError(f"{dataset}: no preprocessed tree-level data")

        for method, method_block in dataset_data.groupby(
            "detrending_method", sort=True, observed=True
        ):
            structure_rows.append(structure_row(dataset, str(method), method_block))
        sensitivity_tables.append(
            chronology_sensitivity(dataset, dataset_data, PRIMARY_METHOD)
        )

        primary = dataset_data.loc[
            dataset_data["detrending_method"] == PRIMARY_METHOD,
            ["year", "tree_id", "value_z"],
        ].copy()
        primary = primary.rename(columns={"value_z": "value"})
        if primary.empty:
            raise ValueError(f"{dataset}: primary method {PRIMARY_METHOD} is absent")
        min_overlap = int(contract.analysis_years)

        subset_membership = pd.read_csv(
            membership_dir / f"{dataset}_subset_memberships_v01.csv.gz",
            dtype={"replicate_id": str, "unit_id": str},
        )
        split_membership = pd.read_csv(
            membership_dir / f"{dataset}_split_half_memberships_v01.csv.gz",
            dtype={"replicate_id": str, "group": str, "unit_id": str},
        )
        bootstrap_membership = pd.read_csv(
            membership_dir / f"{dataset}_bootstrap_multiplicities_v01.csv.gz",
            dtype={"replicate_id": str, "unit_id": str},
        )

        subset_result = subset_stability(
            primary,
            subset_membership,
            year_col="year",
            unit_col="tree_id",
            value_col="value",
            chronology_method="mean",
            min_depth_fraction=1.0,
            min_overlap=min_overlap,
        )
        split_result = split_half_stability(
            primary,
            split_membership,
            year_col="year",
            unit_col="tree_id",
            value_col="value",
            chronology_method="mean",
            min_depth_fraction=1.0,
            min_overlap=min_overlap,
            require_exhaustive=True,
            require_balanced=True,
        )
        bootstrap_result = tree_bootstrap(
            primary,
            bootstrap_membership,
            year_col="year",
            unit_col="tree_id",
            value_col="value",
            chronology_method="mean",
            min_depth_fraction=1.0,
            min_overlap=min_overlap,
            require_full_size=True,
        )
        loto_result = leave_one_tree_out(
            primary,
            year_col="year",
            unit_col="tree_id",
            value_col="value",
            chronology_method="mean",
            min_depth_fraction=1.0,
            min_overlap=min_overlap,
        )

        raw_tables = {
            "subset": subset_result,
            "split_half": split_result,
            "bootstrap": bootstrap_result,
            "loto": loto_result,
        }
        for label, frame in raw_tables.items():
            path = output_dir / f"{dataset}_{label}_results_v01.csv.gz"
            write_csv_gzip_reproducible(frame, path)
            output_paths.append(path)

        subset_summary = summarize_subset(dataset, subset_result)
        split_summary = summarize_split(dataset, split_result)
        bootstrap_summary = summarize_bootstrap(dataset, bootstrap_result)
        loto_summary = summarize_loto(dataset, loto_result)
        subset_summaries.append(subset_summary)
        split_summaries.append(split_summary)
        bootstrap_summaries.append(bootstrap_summary)
        loto_summaries.append(loto_summary)

        n_trees = int(primary["tree_id"].nunique())
        half_size = n_trees // 2
        half_row = subset_summary.loc[
            subset_summary["subset_size"] == half_size
        ].iloc[0]
        split_row = split_summary.iloc[0]
        boot_row = bootstrap_summary.iloc[0]
        loto_row = loto_summary.iloc[0]
        structure = next(
            row
            for row in structure_rows
            if row["dataset"] == dataset
            and row["detrending_method"] == PRIMARY_METHOD
        )
        key_rows.append(
            {
                "dataset": dataset,
                "n_trees": n_trees,
                "analysis_start": int(contract.analysis_start),
                "analysis_end": int(contract.analysis_end),
                "analysis_years": int(contract.analysis_years),
                "rbar": structure["rbar"],
                "eps": structure["eps"],
                "effective_dimension": structure["effective_dimension"],
                "half_subset_size": half_size,
                "half_subset_correlation_q05": half_row["correlation_q05"],
                "half_subset_correlation_median": half_row["correlation_median"],
                "half_subset_correlation_min": half_row["correlation_min"],
                "split_correlation_q05": split_row["correlation_q05"],
                "split_correlation_median": split_row["correlation_median"],
                "split_correlation_min": split_row["correlation_min"],
                "bootstrap_srms_q95": boot_row[
                    "standardized_rms_difference_q95"
                ],
                "bootstrap_variance_ratio_q05": boot_row["variance_ratio_q05"],
                "bootstrap_variance_ratio_q95": boot_row["variance_ratio_q95"],
                "loto_max_srms": loto_row[
                    "maximum_standardized_rms_difference"
                ],
                "loto_max_abs_delta_eps": loto_row["maximum_absolute_delta_eps"],
                "loto_max_abs_delta_effective_dimension": loto_row[
                    "maximum_absolute_delta_effective_dimension"
                ],
            }
        )

    compact_tables = {
        "structure_summary_v01.csv": pd.DataFrame(structure_rows),
        "detrending_sensitivity_summary_v01.csv": pd.concat(
            sensitivity_tables, ignore_index=True
        ),
        "subset_summary_v01.csv": pd.concat(subset_summaries, ignore_index=True),
        "split_half_summary_v01.csv": pd.concat(split_summaries, ignore_index=True),
        "bootstrap_summary_v01.csv": pd.concat(
            bootstrap_summaries, ignore_index=True
        ),
        "loto_summary_v01.csv": pd.concat(loto_summaries, ignore_index=True),
        "empirical_key_results_v01.csv": pd.DataFrame(key_rows),
    }
    for name, frame in compact_tables.items():
        path = output_dir / name
        write_csv(frame, path)
        output_paths.append(path)

    manifest = {
        "analysis": "itrdb_empirical_stage_b",
        "version": VERSION,
        "primary_detrending_method": PRIMARY_METHOD,
        "chronology_method": "arithmetic mean",
        "min_depth_fraction": 1.0,
        "tie_tolerance": TIE_TOLERANCE,
        "input_files": {
            "dataset_contract": {
                "path": stable_manifest_path(contract_path, repo_root),
                "sha256": file_sha256(contract_path),
            },
            "tree_level": {
                "path": stable_manifest_path(tree_path, repo_root),
                "sha256": file_sha256(tree_path),
            },
            "membership_manifest": {
                "path": stable_manifest_path(
                    membership_dir / "membership_manifest_v01.json", repo_root
                ),
                "sha256": file_sha256(
                    membership_dir / "membership_manifest_v01.json"
                ),
            },
        },
        "outputs": {
            path.name: {
                "path": stable_manifest_path(path, repo_root),
                "sha256": file_sha256(path),
            }
            for path in sorted(output_paths)
        },
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
    }
    manifest_path = output_dir / "analysis_manifest_v01.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    print("ITRDB empirical Stage B analysis completed.")
    print(compact_tables["empirical_key_results_v01.csv"].to_string(index=False))
    print(f"Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
