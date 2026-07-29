from __future__ import annotations

from collections import Counter
from collections.abc import Iterable, Mapping
import math
from typing import Any

import numpy as np
import pandas as pd

from .composition import (
    ChronologyMethod,
    _bounded_pearson,
    _build_chronology,
    _chronology_method,
    _depth_fraction,
    _nonempty_identifier,
    _positive_integer,
    _prepare_matrix,
)
from .core import effective_dimension, eps_value

_STATUS_OK = "ok"
_STATUS_INSUFFICIENT = "insufficient_overlap"
_STATUS_CONSTANT = "constant_chronology"
_STRUCTURE_OK = "ok"
_STRUCTURE_INSUFFICIENT = "insufficient_complete_period"
_STRUCTURE_CONSTANT = "constant_series"
_STRUCTURE_INVALID = "invalid_structure"
_CORRELATION_METHOD = "pearson"

_BOOTSTRAP_COLUMNS = [
    "diagnostic",
    "replicate_id",
    "chronology_method",
    "correlation_method",
    "min_depth_fraction",
    "min_overlap",
    "full_pool_size",
    "nominal_draw_count",
    "unique_unit_count",
    "omitted_unit_count",
    "duplicated_unit_count",
    "maximum_multiplicity",
    "kish_effective_unit_count",
    "shared_unit_count",
    "all_units_represented",
    "correlation",
    "squared_correlation",
    "n_years",
    "start_year",
    "end_year",
    "mean_depth_bootstrap",
    "min_depth_bootstrap",
    "mean_depth_full",
    "min_depth_full",
    "mean_difference",
    "rms_difference",
    "standardized_rms_difference",
    "mean_absolute_difference",
    "max_absolute_difference",
    "max_abs_year",
    "signed_difference_at_max",
    "variance_ratio",
    "bootstrap_ids",
    "status",
]

_LOTO_COLUMNS = [
    "diagnostic",
    "omitted_unit_id",
    "chronology_method",
    "correlation_method",
    "min_depth_fraction",
    "min_overlap",
    "full_pool_size",
    "retained_unit_count",
    "shared_unit_count",
    "correlation",
    "squared_correlation",
    "n_years",
    "start_year",
    "end_year",
    "mean_depth_reduced",
    "min_depth_reduced",
    "mean_depth_full",
    "min_depth_full",
    "mean_difference",
    "rms_difference",
    "standardized_rms_difference",
    "mean_absolute_difference",
    "max_absolute_difference",
    "max_abs_year",
    "signed_difference_at_max",
    "variance_ratio",
    "structure_n_years",
    "full_rbar",
    "reduced_rbar",
    "delta_rbar",
    "full_eps",
    "reduced_eps",
    "delta_eps",
    "full_effective_dimension",
    "reduced_effective_dimension",
    "delta_effective_dimension",
    "structure_status",
    "status",
]


def _positive_multiplicity(value: Any) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise TypeError("multiplicity must contain positive integers")
    try:
        numeric = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise TypeError("multiplicity must contain positive integers") from exc
    if not math.isfinite(numeric) or numeric < 1 or not numeric.is_integer():
        raise ValueError("multiplicity must contain positive integers")
    return int(numeric)


def _normalize_bootstrap_multiplicities(
    multiplicities: pd.DataFrame | Mapping[Any, Iterable[Any] | Any],
) -> pd.DataFrame:
    """Normalize explicit bootstrap draws to replicate/unit/multiplicity rows."""
    if isinstance(multiplicities, Mapping):
        if not multiplicities:
            raise ValueError("multiplicities must contain at least one replicate")
        rows: list[dict[str, object]] = []
        for raw_replicate_id, raw_draws in multiplicities.items():
            replicate_id = _nonempty_identifier(raw_replicate_id, "multiplicities")
            if isinstance(raw_draws, (str, bytes)):
                draws = [raw_draws]
            else:
                try:
                    draws = list(raw_draws)
                except TypeError:
                    draws = [raw_draws]
            if not draws:
                raise ValueError(
                    f"bootstrap replicate {replicate_id!r} contains no draws"
                )
            normalized = [
                _nonempty_identifier(unit, "multiplicities") for unit in draws
            ]
            for unit_id, count in sorted(Counter(normalized).items()):
                rows.append(
                    {
                        "replicate_id": replicate_id,
                        "unit_id": unit_id,
                        "multiplicity": int(count),
                    }
                )
        result = pd.DataFrame(
            rows, columns=["replicate_id", "unit_id", "multiplicity"]
        )
    elif isinstance(multiplicities, pd.DataFrame):
        required = {"replicate_id", "unit_id", "multiplicity"}
        missing = required.difference(multiplicities.columns)
        if missing:
            raise ValueError(
                "multiplicities is missing required columns: "
                f"{sorted(missing)}"
            )
        result = multiplicities.loc[
            :, ["replicate_id", "unit_id", "multiplicity"]
        ].copy()
        if result[["replicate_id", "unit_id"]].isna().any().any():
            raise ValueError("multiplicities contains missing identifiers")
        result["replicate_id"] = result["replicate_id"].map(
            lambda value: _nonempty_identifier(value, "multiplicities")
        )
        result["unit_id"] = result["unit_id"].map(
            lambda value: _nonempty_identifier(value, "multiplicities")
        )
        result["multiplicity"] = result["multiplicity"].map(
            _positive_multiplicity
        )
    else:
        raise TypeError("multiplicities must be a DataFrame or mapping")

    if result.empty:
        raise ValueError("multiplicities must contain at least one replicate")
    if result.duplicated(["replicate_id", "unit_id"]).any():
        raise ValueError(
            "multiplicities contains duplicate replicate/unit rows; combine them "
            "into one positive multiplicity"
        )
    return result.sort_values(
        ["replicate_id", "unit_id"], kind="stable"
    ).reset_index(drop=True)


def _validate_bootstrap_units(
    block: pd.DataFrame, available: set[str], replicate_id: str
) -> tuple[list[str], np.ndarray]:
    units = block["unit_id"].astype(str).tolist()
    missing = sorted(set(units).difference(available))
    if missing:
        raise ValueError(
            f"bootstrap replicate {replicate_id!r} refers to units absent from "
            f"data: {missing}"
        )
    order = np.argsort(np.asarray(units, dtype=object), kind="stable")
    sorted_units = [units[int(index)] for index in order]
    weights = block["multiplicity"].to_numpy(dtype=int)[order]
    return sorted_units, weights


def _build_weighted_chronology(
    matrix: pd.DataFrame,
    units: list[str],
    weights: np.ndarray,
    *,
    method: ChronologyMethod,
    min_depth_fraction: float,
) -> tuple[pd.Series, pd.Series]:
    selected = matrix.loc[:, units]
    values = selected.to_numpy(dtype=float)
    observed = np.isfinite(values)
    weighted_depth = observed @ weights.astype(int)
    nominal_draw_count = int(weights.sum())
    required_depth = int(math.ceil(min_depth_fraction * nominal_draw_count))

    if method == "mean":
        numerator = np.nansum(values * weights.reshape(1, -1), axis=1)
        chronology_values = np.divide(
            numerator,
            weighted_depth,
            out=np.full(values.shape[0], np.nan, dtype=float),
            where=weighted_depth > 0,
        )
    else:
        chronology_values = np.full(values.shape[0], np.nan, dtype=float)
        for row_index in range(values.shape[0]):
            mask = observed[row_index]
            if not mask.any():
                continue
            expanded = np.repeat(values[row_index, mask], weights[mask])
            chronology_values[row_index] = float(np.median(expanded))

    chronology_values[weighted_depth < required_depth] = np.nan
    chronology = pd.Series(
        chronology_values, index=matrix.index, dtype=float, name="chronology"
    )
    depth = pd.Series(
        weighted_depth.astype(int), index=matrix.index, dtype=int, name="depth"
    )
    return chronology, depth


def _comparison_metrics(
    chronology_a: pd.Series,
    depth_a: pd.Series,
    chronology_b: pd.Series,
    depth_b: pd.Series,
    *,
    min_overlap: int,
    depth_a_label: str,
    depth_b_label: str,
) -> dict[str, object]:
    valid = chronology_a.notna() & chronology_b.notna()
    n_years = int(valid.sum())
    years = chronology_a.index[valid]

    result: dict[str, object] = {
        "correlation": math.nan,
        "squared_correlation": math.nan,
        "n_years": n_years,
        "start_year": int(years.min()) if n_years else None,
        "end_year": int(years.max()) if n_years else None,
        f"mean_depth_{depth_a_label}": (
            float(depth_a.loc[valid].mean()) if n_years else math.nan
        ),
        f"min_depth_{depth_a_label}": (
            int(depth_a.loc[valid].min()) if n_years else 0
        ),
        f"mean_depth_{depth_b_label}": (
            float(depth_b.loc[valid].mean()) if n_years else math.nan
        ),
        f"min_depth_{depth_b_label}": (
            int(depth_b.loc[valid].min()) if n_years else 0
        ),
        "mean_difference": math.nan,
        "rms_difference": math.nan,
        "standardized_rms_difference": math.nan,
        "mean_absolute_difference": math.nan,
        "max_absolute_difference": math.nan,
        "max_abs_year": None,
        "signed_difference_at_max": math.nan,
        "variance_ratio": math.nan,
        "status": _STATUS_INSUFFICIENT,
    }
    if n_years < min_overlap:
        return result

    a = chronology_a.loc[valid].to_numpy(dtype=float)
    b = chronology_b.loc[valid].to_numpy(dtype=float)
    difference = a - b
    absolute = np.abs(difference)
    maximum_index = int(np.argmax(absolute))
    variance_a = float(np.var(a, ddof=1))
    variance_b = float(np.var(b, ddof=1))
    sd_b = math.sqrt(variance_b) if variance_b > 0.0 else 0.0
    rms = float(math.sqrt(float(np.mean(np.square(difference)))))

    result.update(
        {
            "mean_difference": float(np.mean(difference)),
            "rms_difference": rms,
            "standardized_rms_difference": (
                rms / sd_b if sd_b > 0.0 else math.nan
            ),
            "mean_absolute_difference": float(np.mean(absolute)),
            "max_absolute_difference": float(absolute[maximum_index]),
            "max_abs_year": int(years[maximum_index]),
            "signed_difference_at_max": float(difference[maximum_index]),
            "variance_ratio": (
                variance_a / variance_b if variance_b > 0.0 else math.nan
            ),
        }
    )

    if np.ptp(a) == 0.0 or np.ptp(b) == 0.0:
        result["status"] = _STATUS_CONSTANT
        return result

    correlation = _bounded_pearson(a, b)
    result["correlation"] = correlation
    result["squared_correlation"] = correlation * correlation
    result["status"] = _STATUS_OK
    return result


def _mean_off_diagonal(correlation: np.ndarray) -> float:
    n = correlation.shape[0]
    if n < 2:
        raise ValueError("at least two units are required")
    return float((correlation.sum() - np.trace(correlation)) / (n * (n - 1)))


def _loto_structure_metrics(
    matrix: pd.DataFrame,
    omitted_unit: str,
    *,
    min_overlap: int,
) -> dict[str, object]:
    complete = matrix.notna().all(axis=1)
    n_years = int(complete.sum())
    empty = {
        "structure_n_years": n_years,
        "full_rbar": math.nan,
        "reduced_rbar": math.nan,
        "delta_rbar": math.nan,
        "full_eps": math.nan,
        "reduced_eps": math.nan,
        "delta_eps": math.nan,
        "full_effective_dimension": math.nan,
        "reduced_effective_dimension": math.nan,
        "delta_effective_dimension": math.nan,
        "structure_status": _STRUCTURE_INSUFFICIENT,
    }
    if n_years < min_overlap or matrix.shape[1] < 3:
        return empty

    full_values = matrix.loc[complete].to_numpy(dtype=float)
    if np.any(np.ptp(full_values, axis=0) == 0.0):
        empty["structure_status"] = _STRUCTURE_CONSTANT
        return empty

    try:
        full_correlation = np.corrcoef(full_values, rowvar=False)
        reduced_columns = [
            column for column in matrix.columns if column != omitted_unit
        ]
        reduced_values = matrix.loc[complete, reduced_columns].to_numpy(dtype=float)
        reduced_correlation = np.corrcoef(reduced_values, rowvar=False)
        full_rbar = _mean_off_diagonal(full_correlation)
        reduced_rbar = _mean_off_diagonal(reduced_correlation)
        full_eps = float(eps_value(full_rbar, matrix.shape[1]))
        reduced_eps = float(eps_value(reduced_rbar, len(reduced_columns)))
        full_dimension = float(effective_dimension(full_correlation))
        reduced_dimension = float(effective_dimension(reduced_correlation))
    except (ValueError, RuntimeError, FloatingPointError):
        empty["structure_status"] = _STRUCTURE_INVALID
        return empty

    return {
        "structure_n_years": n_years,
        "full_rbar": full_rbar,
        "reduced_rbar": reduced_rbar,
        "delta_rbar": reduced_rbar - full_rbar,
        "full_eps": full_eps,
        "reduced_eps": reduced_eps,
        "delta_eps": reduced_eps - full_eps,
        "full_effective_dimension": full_dimension,
        "reduced_effective_dimension": reduced_dimension,
        "delta_effective_dimension": reduced_dimension - full_dimension,
        "structure_status": _STRUCTURE_OK,
    }


def tree_bootstrap(
    data: pd.DataFrame,
    multiplicities: pd.DataFrame | Mapping[Any, Iterable[Any] | Any],
    *,
    year_col: str = "year",
    unit_col: str = "tree_id",
    value_col: str = "value",
    chronology_method: ChronologyMethod = "mean",
    min_depth_fraction: float = 1.0,
    min_overlap: int = 3,
    require_full_size: bool = True,
) -> pd.DataFrame:
    """Evaluate explicit tree-level bootstrap multiplicities.

    The numerical primitive is an auditable table with columns
    ``replicate_id``, ``unit_id``, and positive integer ``multiplicity``.
    A mapping may instead provide the literal sequence of unit draws for each
    replicate. Random draws are deliberately generated outside this function.

    The bootstrap chronology is compared with the full chronology. Because the
    two chronologies share the unique units represented in the bootstrap draw,
    correlation is descriptive; direct difference metrics are the primary
    composition-sensitivity outputs.
    """
    if not isinstance(require_full_size, (bool, np.bool_)):
        raise TypeError("require_full_size must be True or False")
    method = _chronology_method(chronology_method)
    fraction = _depth_fraction(min_depth_fraction)
    overlap = _positive_integer(min_overlap, "min_overlap")
    matrix = _prepare_matrix(
        data, year_col=year_col, unit_col=unit_col, value_col=value_col
    )
    table = _normalize_bootstrap_multiplicities(multiplicities)
    available = set(matrix.columns)
    full_units = sorted(available)
    full_chronology, full_depth = _build_chronology(
        matrix,
        full_units,
        method=method,
        min_depth_fraction=fraction,
    )

    rows: list[dict[str, object]] = []
    for replicate_id, block in table.groupby(
        "replicate_id", sort=True, observed=True
    ):
        units, weights = _validate_bootstrap_units(
            block, available, str(replicate_id)
        )
        nominal_draw_count = int(weights.sum())
        if require_full_size and nominal_draw_count != len(full_units):
            raise ValueError(
                f"bootstrap replicate {replicate_id!r} has "
                f"{nominal_draw_count} draws; expected {len(full_units)}"
            )
        bootstrap_chronology, bootstrap_depth = _build_weighted_chronology(
            matrix,
            units,
            weights,
            method=method,
            min_depth_fraction=fraction,
        )
        metrics = _comparison_metrics(
            bootstrap_chronology,
            bootstrap_depth,
            full_chronology,
            full_depth,
            min_overlap=overlap,
            depth_a_label="bootstrap",
            depth_b_label="full",
        )
        unique_count = len(units)
        weight_sum = float(weights.sum())
        kish = weight_sum * weight_sum / float(np.square(weights).sum())
        row = {
            "diagnostic": "tree_bootstrap_to_full",
            "replicate_id": str(replicate_id),
            "chronology_method": method,
            "correlation_method": _CORRELATION_METHOD,
            "min_depth_fraction": fraction,
            "min_overlap": overlap,
            "full_pool_size": len(full_units),
            "nominal_draw_count": nominal_draw_count,
            "unique_unit_count": unique_count,
            "omitted_unit_count": len(full_units) - unique_count,
            "duplicated_unit_count": int((weights > 1).sum()),
            "maximum_multiplicity": int(weights.max()),
            "kish_effective_unit_count": float(kish),
            "shared_unit_count": unique_count,
            "all_units_represented": set(units) == available,
            "bootstrap_ids": "|".join(
                f"{unit}*{int(weight)}" for unit, weight in zip(units, weights)
            ),
            **metrics,
        }
        rows.append(row)
    return pd.DataFrame(rows, columns=_BOOTSTRAP_COLUMNS)


def leave_one_tree_out(
    data: pd.DataFrame,
    *,
    year_col: str = "year",
    unit_col: str = "tree_id",
    value_col: str = "value",
    chronology_method: ChronologyMethod = "mean",
    min_depth_fraction: float = 1.0,
    min_overlap: int = 3,
) -> pd.DataFrame:
    """Quantify the influence of omitting each biological unit in turn.

    The reduced chronology and full chronology share ``N - 1`` units, so their
    correlation is expected to be very high and is not treated as the primary
    influence measure. Direct chronology differences and common-period changes in
    mean correlation, EPS, and effective dimension are reported instead.
    """
    method = _chronology_method(chronology_method)
    fraction = _depth_fraction(min_depth_fraction)
    overlap = _positive_integer(min_overlap, "min_overlap")
    matrix = _prepare_matrix(
        data, year_col=year_col, unit_col=unit_col, value_col=value_col
    )
    full_units = sorted(matrix.columns)
    if len(full_units) < 2:
        raise ValueError("leave-one-tree-out requires at least two biological units")
    full_chronology, full_depth = _build_chronology(
        matrix,
        full_units,
        method=method,
        min_depth_fraction=fraction,
    )

    rows: list[dict[str, object]] = []
    for omitted_unit in full_units:
        retained = [unit for unit in full_units if unit != omitted_unit]
        reduced_chronology, reduced_depth = _build_chronology(
            matrix,
            retained,
            method=method,
            min_depth_fraction=fraction,
        )
        metrics = _comparison_metrics(
            reduced_chronology,
            reduced_depth,
            full_chronology,
            full_depth,
            min_overlap=overlap,
            depth_a_label="reduced",
            depth_b_label="full",
        )
        structure = _loto_structure_metrics(
            matrix, omitted_unit, min_overlap=overlap
        )
        rows.append(
            {
                "diagnostic": "leave_one_tree_out",
                "omitted_unit_id": omitted_unit,
                "chronology_method": method,
                "correlation_method": _CORRELATION_METHOD,
                "min_depth_fraction": fraction,
                "min_overlap": overlap,
                "full_pool_size": len(full_units),
                "retained_unit_count": len(retained),
                "shared_unit_count": len(retained),
                **metrics,
                **structure,
            }
        )
    return pd.DataFrame(rows, columns=_LOTO_COLUMNS)
