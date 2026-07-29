from __future__ import annotations

from collections.abc import Iterable, Mapping
import math
from typing import Any, Literal

import numpy as np
import pandas as pd

ChronologyMethod = Literal["mean", "median"]
_STATUS_OK = "ok"
_STATUS_INSUFFICIENT = "insufficient_overlap"
_STATUS_CONSTANT = "constant_chronology"
_CORRELATION_METHOD = "pearson"

_RESULT_COLUMNS = [
    "diagnostic",
    "replicate_id",
    "chronology_method",
    "correlation_method",
    "min_depth_fraction",
    "min_overlap",
    "full_pool_size",
    "group_a_size",
    "group_b_size",
    "shared_unit_count",
    "comparison_exhaustive",
    "correlation",
    "squared_correlation",
    "n_years",
    "start_year",
    "end_year",
    "mean_depth_a",
    "min_depth_a",
    "mean_depth_b",
    "min_depth_b",
    "group_a_ids",
    "group_b_ids",
    "status",
]


def _positive_integer(value: Any, name: str) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(f"{name} must be a positive integer")
    try:
        number = int(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise TypeError(f"{name} must be a positive integer") from exc
    if number < 1 or number != value:
        raise ValueError(f"{name} must be a positive integer")
    return number


def _depth_fraction(value: Any) -> float:
    try:
        fraction = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise TypeError("min_depth_fraction must be a finite number in (0, 1]") from exc
    if not math.isfinite(fraction) or fraction <= 0.0 or fraction > 1.0:
        raise ValueError("min_depth_fraction must lie in (0, 1]")
    return fraction


def _chronology_method(value: str) -> ChronologyMethod:
    method = str(value).casefold()
    if method not in {"mean", "median"}:
        raise ValueError("chronology_method must be 'mean' or 'median'")
    return method  # type: ignore[return-value]


def _nonempty_identifier(value: object, name: str) -> str:
    text = str(value)
    if not text.strip():
        raise ValueError(f"{name} contains empty identifiers")
    return text


def _prepare_matrix(
    data: pd.DataFrame,
    *,
    year_col: str,
    unit_col: str,
    value_col: str,
) -> pd.DataFrame:
    if not isinstance(data, pd.DataFrame):
        raise TypeError("data must be a pandas DataFrame")
    required = {year_col, unit_col, value_col}
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"data is missing required columns: {sorted(missing)}")

    source = data.loc[:, [year_col, unit_col, value_col]].copy()
    if source[unit_col].isna().any():
        raise ValueError(f"{unit_col} contains missing identifiers")
    source[unit_col] = source[unit_col].map(
        lambda value: _nonempty_identifier(value, unit_col)
    )

    years = pd.to_numeric(source[year_col], errors="raise")
    year_values = years.to_numpy(dtype=float)
    if not np.isfinite(year_values).all() or not np.equal(
        year_values, np.floor(year_values)
    ).all():
        raise ValueError(f"{year_col} must contain finite integer years")
    source[year_col] = years.astype(int)

    values = pd.to_numeric(source[value_col], errors="coerce")
    bad = source[value_col].notna() & values.isna()
    if bad.any():
        examples = source.loc[bad, value_col].astype(str).head(3).tolist()
        raise ValueError(f"{value_col} contains non-numeric values: {examples}")
    numeric = values.to_numpy(dtype=float)
    finite_or_missing = np.isfinite(numeric) | values.isna().to_numpy()
    if not finite_or_missing.all():
        raise ValueError(f"{value_col} must contain finite numbers or missing values")
    source[value_col] = values

    duplicate = source.duplicated([year_col, unit_col], keep=False)
    if duplicate.any():
        row = source.loc[duplicate, [year_col, unit_col]].iloc[0]
        raise ValueError(
            "data must contain one value per year and biological unit; "
            f"duplicate found for {row[unit_col]!r} in {int(row[year_col])}"
        )

    matrix = source.pivot(index=year_col, columns=unit_col, values=value_col)
    matrix = matrix.sort_index().reindex(sorted(matrix.columns), axis=1)
    matrix.index = matrix.index.astype(int)
    matrix.columns = matrix.columns.astype(str)
    empty_units = matrix.columns[matrix.notna().sum(axis=0) == 0].tolist()
    if empty_units:
        raise ValueError(
            "biological units contain no finite observations: "
            f"{sorted(empty_units)}"
        )
    return matrix


def _normalize_subsets(
    subsets: pd.DataFrame | Mapping[Any, Iterable[Any] | Any],
) -> pd.DataFrame:
    if isinstance(subsets, Mapping):
        rows: list[dict[str, str]] = []
        if not subsets:
            raise ValueError("subsets must contain at least one replicate")
        for raw_replicate_id, raw_units in subsets.items():
            replicate_id = _nonempty_identifier(raw_replicate_id, "subsets")
            if isinstance(raw_units, (str, bytes)):
                units = [raw_units]
            else:
                try:
                    units = list(raw_units)
                except TypeError:
                    units = [raw_units]
            if not units:
                raise ValueError(
                    f"subset replicate {replicate_id!r} contains no biological units"
                )
            for unit in units:
                rows.append(
                    {
                        "replicate_id": replicate_id,
                        "unit_id": _nonempty_identifier(unit, "subsets"),
                    }
                )
        result = pd.DataFrame(rows, columns=["replicate_id", "unit_id"])
    elif isinstance(subsets, pd.DataFrame):
        required = {"replicate_id", "unit_id"}
        missing = required.difference(subsets.columns)
        if missing:
            raise ValueError(
                f"subsets is missing required columns: {sorted(missing)}"
            )
        result = subsets.loc[:, ["replicate_id", "unit_id"]].copy()
        if result.isna().any().any():
            raise ValueError("subsets contains missing replicate_id or unit_id values")
        result["replicate_id"] = result["replicate_id"].map(
            lambda value: _nonempty_identifier(value, "subsets")
        )
        result["unit_id"] = result["unit_id"].map(
            lambda value: _nonempty_identifier(value, "subsets")
        )
    else:
        raise TypeError("subsets must be a DataFrame or mapping")

    if result.empty:
        raise ValueError("subsets must contain at least one replicate")
    if result.duplicated(["replicate_id", "unit_id"]).any():
        raise ValueError("subsets contains duplicate unit assignments")
    return result


def _normalize_split_partitions(partitions: pd.DataFrame) -> pd.DataFrame:
    if not isinstance(partitions, pd.DataFrame):
        raise TypeError("partitions must be a pandas DataFrame")
    required = {"replicate_id", "group", "unit_id"}
    missing = required.difference(partitions.columns)
    if missing:
        raise ValueError(
            f"partitions is missing required columns: {sorted(missing)}"
        )
    result = partitions.loc[:, ["replicate_id", "group", "unit_id"]].copy()
    if result.isna().any().any():
        raise ValueError("partitions contains missing identifiers")
    result["replicate_id"] = result["replicate_id"].map(
        lambda value: _nonempty_identifier(value, "partitions")
    )
    result["unit_id"] = result["unit_id"].map(
        lambda value: _nonempty_identifier(value, "partitions")
    )
    result["group"] = result["group"].map(
        lambda value: _nonempty_identifier(value, "partitions").upper()
    )
    if result.empty:
        raise ValueError("partitions must contain at least one replicate")
    if not set(result["group"]).issubset({"A", "B"}):
        raise ValueError("partition group labels must be A or B")
    if result.duplicated(["replicate_id", "group", "unit_id"]).any():
        raise ValueError("partitions contains duplicate unit assignments")
    return result


def _validate_units(requested: Iterable[str], available: set[str], name: str) -> list[str]:
    units = sorted(set(str(unit) for unit in requested))
    if not units:
        raise ValueError(f"{name} contains no biological units")
    missing = sorted(set(units).difference(available))
    if missing:
        raise ValueError(f"{name} refers to units absent from data: {missing}")
    return units


def _build_chronology(
    matrix: pd.DataFrame,
    units: list[str],
    *,
    method: ChronologyMethod,
    min_depth_fraction: float,
) -> tuple[pd.Series, pd.Series]:
    selected = matrix.loc[:, units]
    depth = selected.notna().sum(axis=1).astype(int)
    required_depth = int(math.ceil(min_depth_fraction * len(units)))
    if method == "mean":
        chronology = selected.mean(axis=1, skipna=True)
    else:
        chronology = selected.median(axis=1, skipna=True)
    chronology = chronology.where(depth >= required_depth)
    return chronology.astype(float), depth


def _bounded_pearson(a: np.ndarray, b: np.ndarray) -> float:
    correlation = float(np.corrcoef(a, b)[0, 1])
    tolerance = 64.0 * np.finfo(float).eps
    if correlation > 1.0 and correlation <= 1.0 + tolerance:
        return 1.0
    if correlation < -1.0 and correlation >= -1.0 - tolerance:
        return -1.0
    if not math.isfinite(correlation) or correlation < -1.0 or correlation > 1.0:
        raise RuntimeError("Pearson correlation fell outside [-1, 1]")
    return correlation


def _comparison_row(
    *,
    diagnostic: str,
    replicate_id: str,
    matrix: pd.DataFrame,
    group_a: list[str],
    group_b: list[str],
    chronology_method: ChronologyMethod,
    min_depth_fraction: float,
    min_overlap: int,
    full_pool_size: int,
) -> dict[str, object]:
    chronology_a, depth_a = _build_chronology(
        matrix,
        group_a,
        method=chronology_method,
        min_depth_fraction=min_depth_fraction,
    )
    chronology_b, depth_b = _build_chronology(
        matrix,
        group_b,
        method=chronology_method,
        min_depth_fraction=min_depth_fraction,
    )
    valid = chronology_a.notna() & chronology_b.notna()
    n_years = int(valid.sum())
    years = matrix.index[valid]

    correlation = math.nan
    status = _STATUS_INSUFFICIENT
    if n_years >= min_overlap:
        a = chronology_a.loc[valid].to_numpy(dtype=float)
        b = chronology_b.loc[valid].to_numpy(dtype=float)
        if np.ptp(a) == 0.0 or np.ptp(b) == 0.0:
            status = _STATUS_CONSTANT
        else:
            correlation = _bounded_pearson(a, b)
            status = _STATUS_OK

    if n_years:
        mean_depth_a = float(depth_a.loc[valid].mean())
        min_depth_a = int(depth_a.loc[valid].min())
        mean_depth_b = float(depth_b.loc[valid].mean())
        min_depth_b = int(depth_b.loc[valid].min())
        start_year: int | None = int(years.min())
        end_year: int | None = int(years.max())
    else:
        mean_depth_a = math.nan
        min_depth_a = 0
        mean_depth_b = math.nan
        min_depth_b = 0
        start_year = None
        end_year = None

    shared = len(set(group_a).intersection(group_b))
    return {
        "diagnostic": diagnostic,
        "replicate_id": replicate_id,
        "chronology_method": chronology_method,
        "correlation_method": _CORRELATION_METHOD,
        "min_depth_fraction": float(min_depth_fraction),
        "min_overlap": int(min_overlap),
        "full_pool_size": int(full_pool_size),
        "group_a_size": len(group_a),
        "group_b_size": len(group_b),
        "shared_unit_count": shared,
        "comparison_exhaustive": len(set(group_a).union(group_b)) == full_pool_size,
        "correlation": correlation,
        "squared_correlation": correlation * correlation,
        "n_years": n_years,
        "start_year": start_year,
        "end_year": end_year,
        "mean_depth_a": mean_depth_a,
        "min_depth_a": min_depth_a,
        "mean_depth_b": mean_depth_b,
        "min_depth_b": min_depth_b,
        "group_a_ids": "|".join(group_a),
        "group_b_ids": "|".join(group_b),
        "status": status,
    }


def subset_stability(
    data: pd.DataFrame,
    subsets: pd.DataFrame | Mapping[Any, Iterable[Any] | Any],
    *,
    year_col: str = "year",
    unit_col: str = "tree_id",
    value_col: str = "value",
    chronology_method: ChronologyMethod = "mean",
    min_depth_fraction: float = 1.0,
    min_overlap: int = 3,
) -> pd.DataFrame:
    """Compare explicit biological-unit subsets with the full chronology.

    The returned ``squared_correlation`` is an empirical counterpart of SSS only
    under the conditions that define SSS: arithmetic-mean chronologies,
    exchangeable series, a fixed analysis period, and a subset contained in the
    full set. The function therefore reports both Pearson correlation and its
    square, as well as the number of biological units shared by the chronologies.

    Random subsets are deliberately not generated inside this function. Pass an
    explicit membership table so that the same compositions can be audited and
    reproduced in R and Python.
    """
    method = _chronology_method(chronology_method)
    fraction = _depth_fraction(min_depth_fraction)
    overlap = _positive_integer(min_overlap, "min_overlap")
    matrix = _prepare_matrix(
        data, year_col=year_col, unit_col=unit_col, value_col=value_col
    )
    partition_table = _normalize_subsets(subsets)
    available = set(matrix.columns)
    full_units = sorted(available)

    rows: list[dict[str, object]] = []
    for replicate_id, block in partition_table.groupby(
        "replicate_id", sort=True, observed=True
    ):
        subset_units = _validate_units(
            block["unit_id"], available, f"subset replicate {replicate_id!r}"
        )
        rows.append(
            _comparison_row(
                diagnostic="subset_to_full",
                replicate_id=str(replicate_id),
                matrix=matrix,
                group_a=subset_units,
                group_b=full_units,
                chronology_method=method,
                min_depth_fraction=fraction,
                min_overlap=overlap,
                full_pool_size=len(full_units),
            )
        )
    return pd.DataFrame(rows, columns=_RESULT_COLUMNS)


def split_half_stability(
    data: pd.DataFrame,
    partitions: pd.DataFrame,
    *,
    year_col: str = "year",
    unit_col: str = "tree_id",
    value_col: str = "value",
    chronology_method: ChronologyMethod = "mean",
    min_depth_fraction: float = 1.0,
    min_overlap: int = 3,
    require_exhaustive: bool = True,
    require_balanced: bool = True,
) -> pd.DataFrame:
    """Compare explicit, non-overlapping A/B biological-unit partitions.

    Each replicate must contain groups ``A`` and ``B``. Units may not occur in
    both groups. By default the two groups must exhaust all biological units in
    ``data``; set ``require_exhaustive=False`` only for a deliberately partial,
    but still disjoint, comparison. By default the two groups may differ by at
    most one unit, as required for a literal split half.
    """
    if not isinstance(require_exhaustive, (bool, np.bool_)):
        raise TypeError("require_exhaustive must be True or False")
    if not isinstance(require_balanced, (bool, np.bool_)):
        raise TypeError("require_balanced must be True or False")
    method = _chronology_method(chronology_method)
    fraction = _depth_fraction(min_depth_fraction)
    overlap = _positive_integer(min_overlap, "min_overlap")
    matrix = _prepare_matrix(
        data, year_col=year_col, unit_col=unit_col, value_col=value_col
    )
    partition_table = _normalize_split_partitions(partitions)
    available = set(matrix.columns)

    rows: list[dict[str, object]] = []
    for replicate_id, block in partition_table.groupby(
        "replicate_id", sort=True, observed=True
    ):
        labels = set(block["group"])
        if labels != {"A", "B"}:
            raise ValueError(
                f"split-half replicate {replicate_id!r} must contain groups A and B"
            )
        group_a = _validate_units(
            block.loc[block["group"] == "A", "unit_id"],
            available,
            f"split-half replicate {replicate_id!r} group A",
        )
        group_b = _validate_units(
            block.loc[block["group"] == "B", "unit_id"],
            available,
            f"split-half replicate {replicate_id!r} group B",
        )
        shared = sorted(set(group_a).intersection(group_b))
        if shared:
            raise ValueError(
                f"split-half replicate {replicate_id!r} is not disjoint: {shared}"
            )
        if require_balanced and abs(len(group_a) - len(group_b)) > 1:
            raise ValueError(
                f"split-half replicate {replicate_id!r} is not balanced; "
                f"group sizes are {len(group_a)} and {len(group_b)}"
            )
        assigned = set(group_a).union(group_b)
        if require_exhaustive and assigned != available:
            omitted = sorted(available.difference(assigned))
            extra = sorted(assigned.difference(available))
            raise ValueError(
                f"split-half replicate {replicate_id!r} is not exhaustive; "
                f"omitted={omitted}, extra={extra}"
            )
        rows.append(
            _comparison_row(
                diagnostic="disjoint_split_half",
                replicate_id=str(replicate_id),
                matrix=matrix,
                group_a=group_a,
                group_b=group_b,
                chronology_method=method,
                min_depth_fraction=fraction,
                min_overlap=overlap,
                full_pool_size=len(available),
            )
        )
    return pd.DataFrame(rows, columns=_RESULT_COLUMNS)
