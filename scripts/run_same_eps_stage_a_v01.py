"""Run the exact same-EPS/different-structure Stage A experiment.

The four fixed N=18 correlation matrices have the same mean off-diagonal
correlation and therefore the same plug-in EPS (=0.85), but different covariance
architecture.  Stage A evaluates composition sensitivity without Monte Carlo
uncertainty for subset and balanced split-half diagnostics:

* every non-empty proper subset is enumerated;
* every unique balanced 9-vs-9 split is enumerated;
* tree bootstrap uses one frozen set of explicit N-out-of-N multiplicities;
* leave-one-tree-out evaluates every biological unit.

The script writes machine-readable summaries used by the manuscript.  It does
not assign a composite score or pass/fail label.
"""
from __future__ import annotations

import argparse
from collections.abc import Iterable
import gzip
import io
from dataclasses import dataclass
from hashlib import sha256
from itertools import combinations
import json
import math
from pathlib import Path
import platform
import sys
from typing import Any

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
PYTHON_SRC = REPO_ROOT / "python" / "src"
if str(PYTHON_SRC) not in sys.path:
    sys.path.insert(0, str(PYTHON_SRC))

from chronodiag import effective_dimension, eps_value  # noqa: E402

VERSION = "v01"
SCENARIOS = (
    "exchangeable",
    "two_clusters",
    "unequal_loadings",
    "influential_trees",
)
DEFAULT_N_YEARS = 100
DEFAULT_START_YEAR = 1900
DEFAULT_BOOTSTRAP_REPLICATES = 1000
DEFAULT_BOOTSTRAP_SEED = 20260730
TOLERANCE = 1e-12
SUMMARY_TIE_TOLERANCE = 1e-12


@dataclass(frozen=True)
class ScenarioMetadata:
    scenario: str
    description: str
    unit_class: tuple[str, ...]


def scenario_metadata() -> dict[str, ScenarioMetadata]:
    n = 18
    return {
        "exchangeable": ScenarioMetadata(
            "exchangeable",
            "Equal off-diagonal correlation under a parallel common-signal model.",
            tuple(["exchangeable"] * n),
        ),
        "two_clusters": ScenarioMetadata(
            "two_clusters",
            "Two equal nine-unit clusters with strong within-cluster and weak between-cluster correlation.",
            tuple(["cluster_A"] * 9 + ["cluster_B"] * 9),
        ),
        "unequal_loadings": ScenarioMetadata(
            "unequal_loadings",
            "One-factor model with smoothly increasing unit loadings.",
            tuple(["graded_loading"] * n),
        ),
        "influential_trees": ScenarioMetadata(
            "influential_trees",
            "One-factor model with three high-loading and fifteen lower-loading units.",
            tuple(["high_loading"] * 3 + ["lower_loading"] * 15),
        ),
    }


def _one_factor_loadings(matrix: np.ndarray) -> np.ndarray:
    """Recover positive one-factor loadings from a standardized correlation matrix.

    For an exact one-factor matrix, off-diagonal entries satisfy
    ``R[i, j] = loading[i] * loading[j]``.  Recovering loadings from the frozen
    matrix, rather than duplicating rounded constants in analysis code, keeps the
    composition-driver audit tied to the scientific input itself.
    """
    if matrix.shape[0] < 3:
        raise ValueError("one-factor loading recovery requires at least three units")
    r01 = float(matrix[0, 1])
    r02 = float(matrix[0, 2])
    r12 = float(matrix[1, 2])
    if min(r01, r02, r12) <= 0.0:
        raise ValueError("positive one-factor loading recovery requires positive correlations")
    loading0 = math.sqrt(r01 * r02 / r12)
    loadings = np.empty(matrix.shape[0], dtype=float)
    loadings[0] = loading0
    loadings[1:] = matrix[0, 1:] / loading0

    expected = np.outer(loadings, loadings)
    off_diagonal = ~np.eye(matrix.shape[0], dtype=bool)
    if not np.allclose(
        matrix[off_diagonal], expected[off_diagonal], atol=5e-13, rtol=0.0
    ):
        maximum = float(
            np.max(np.abs(matrix[off_diagonal] - expected[off_diagonal]))
        )
        raise ValueError(
            "matrix is not an exact positive one-factor correlation structure; "
            f"maximum off-diagonal residual={maximum}"
        )
    if np.any(loadings <= 0.0) or np.any(loadings >= 1.0):
        raise ValueError("recovered one-factor loadings must lie in (0, 1)")
    return loadings


def _scenario_loadings(scenario: str, matrix: np.ndarray) -> np.ndarray:
    if scenario == "exchangeable":
        return np.full(matrix.shape[0], math.sqrt(_mean_off_diagonal(matrix)))
    if scenario in {"unequal_loadings", "influential_trees"}:
        return _one_factor_loadings(matrix)
    return np.full(matrix.shape[0], np.nan, dtype=float)


def _mean_off_diagonal(matrix: np.ndarray) -> float:
    n = matrix.shape[0]
    return float((matrix.sum() - np.trace(matrix)) / (n * (n - 1)))


def _matrix_sha256(matrix: np.ndarray) -> str:
    payload = "\n".join(
        ",".join(format(float(value), ".17g") for value in row) for row in matrix
    )
    return sha256((payload + "\n").encode("utf-8")).hexdigest()


def _file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _manifest_path(path: Path) -> str:
    """Return a stable POSIX path relative to the repository when possible."""
    resolved = path.resolve()
    try:
        return resolved.relative_to(REPO_ROOT.resolve()).as_posix()
    except ValueError:
        return resolved.as_posix()


def _validate_correlation(matrix: np.ndarray, name: str) -> None:
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"{name} is not square")
    if not np.isfinite(matrix).all():
        raise ValueError(f"{name} contains non-finite values")
    if not np.allclose(matrix, matrix.T, atol=TOLERANCE, rtol=0.0):
        raise ValueError(f"{name} is not symmetric")
    if not np.allclose(np.diag(matrix), 1.0, atol=TOLERANCE, rtol=0.0):
        raise ValueError(f"{name} diagonal does not equal one")
    minimum = float(np.linalg.eigvalsh(matrix).min())
    if minimum < -1e-10:
        raise ValueError(f"{name} is not positive semidefinite: {minimum}")


def load_scenario_matrices(
    fixture_dir: Path,
) -> dict[str, np.ndarray]:
    matrices: dict[str, np.ndarray] = {}
    for scenario in SCENARIOS:
        path = fixture_dir / f"{scenario}_correlation.csv"
        matrix = np.loadtxt(path, delimiter=",", dtype=float)
        _validate_correlation(matrix, scenario)
        matrices[scenario] = matrix

    sizes = {matrix.shape[0] for matrix in matrices.values()}
    if len(sizes) != 1:
        raise ValueError(f"scenario matrices have inconsistent dimensions: {sizes}")

    rbars = np.array([_mean_off_diagonal(matrix) for matrix in matrices.values()])
    eps_values = np.array(
        [eps_value(rbar, next(iter(sizes))) for rbar in rbars], dtype=float
    )
    if np.ptp(rbars) > 1e-12 or np.ptp(eps_values) > 1e-12:
        raise ValueError("scenario matrices do not have a common rbar and EPS")
    return matrices


def exact_sample_from_correlation(
    correlation: np.ndarray,
    *,
    n_years: int = DEFAULT_N_YEARS,
) -> np.ndarray:
    """Construct centered values whose sample correlation equals ``correlation``.

    One deterministic centered orthonormal temporal basis is used for every
    scenario.  Scenario differences therefore enter only through the correlation
    matrix.
    """
    n_units = correlation.shape[0]
    if n_years <= n_units:
        raise ValueError("n_years must exceed n_units")

    x = np.arange(1, n_years + 1, dtype=float)
    seed = np.column_stack(
        [
            np.sin(x * (index + 1) * 0.173)
            + 0.37 * np.cos(x * (index + 2) * 0.119)
            + 0.01 * ((x - x.mean()) / n_years) ** ((index % 3) + 1)
            for index in range(n_units)
        ]
    )
    seed -= seed.mean(axis=0)
    orthonormal, _ = np.linalg.qr(seed)
    for column in range(n_units):
        pivot = int(np.argmax(np.abs(orthonormal[:, column])))
        if orthonormal[pivot, column] < 0:
            orthonormal[:, column] *= -1

    values = (
        math.sqrt(n_years - 1)
        * orthonormal
        @ np.linalg.cholesky(correlation).T
    )
    values -= values.mean(axis=0)
    observed = np.corrcoef(values, rowvar=False)
    if not np.allclose(observed, correlation, atol=5e-12, rtol=0.0):
        maximum = float(np.max(np.abs(observed - correlation)))
        raise RuntimeError(f"exact sample construction failed; max error={maximum}")
    return values


def _weight_vector(indices: Iterable[int], n_units: int) -> np.ndarray:
    selected = tuple(indices)
    if not selected:
        raise ValueError("weight vector requires at least one unit")
    weights = np.zeros(n_units, dtype=float)
    weights[list(selected)] = 1.0 / len(selected)
    return weights


def _weighted_correlation(
    matrix: np.ndarray,
    weights_a: np.ndarray,
    weights_b: np.ndarray,
) -> float:
    variance_a = float(weights_a @ matrix @ weights_a)
    variance_b = float(weights_b @ matrix @ weights_b)
    covariance = float(weights_a @ matrix @ weights_b)
    correlation = covariance / math.sqrt(variance_a * variance_b)
    if correlation > 1.0 and correlation < 1.0 + 1e-12:
        correlation = 1.0
    if correlation < -1.0 and correlation > -1.0 - 1e-12:
        correlation = -1.0
    if not math.isfinite(correlation) or correlation < -1.0 or correlation > 1.0:
        raise RuntimeError(f"invalid weighted correlation: {correlation}")
    return correlation


def _comparison_metrics(
    matrix: np.ndarray,
    weights_a: np.ndarray,
    weights_b: np.ndarray,
    *,
    n_years: int,
) -> dict[str, float]:
    variance_a = float(weights_a @ matrix @ weights_a)
    variance_b = float(weights_b @ matrix @ weights_b)
    correlation = _weighted_correlation(matrix, weights_a, weights_b)
    difference = weights_a - weights_b
    mean_square_difference = (
        (n_years - 1) / n_years * float(difference @ matrix @ difference)
    )
    rms = math.sqrt(max(0.0, mean_square_difference))
    return {
        "correlation": correlation,
        "squared_correlation": correlation * correlation,
        "standardized_rms_difference": rms / math.sqrt(variance_b),
        "variance_ratio": variance_a / variance_b,
    }


def _distribution_summary(
    values: np.ndarray,
    *,
    prefix: str = "",
) -> dict[str, float | int]:
    if values.ndim != 1 or values.size == 0:
        raise ValueError("distribution summary requires one or more values")
    q01, q05, q25, q50, q75, q95, q99 = np.quantile(
        values, [0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99]
    )
    return {
        f"{prefix}n": int(values.size),
        f"{prefix}mean": float(values.mean()),
        f"{prefix}sd": float(values.std(ddof=1)) if values.size > 1 else 0.0,
        f"{prefix}q01": float(q01),
        f"{prefix}q05": float(q05),
        f"{prefix}q25": float(q25),
        f"{prefix}median": float(q50),
        f"{prefix}q75": float(q75),
        f"{prefix}q95": float(q95),
        f"{prefix}q99": float(q99),
        f"{prefix}min": float(values.min()),
        f"{prefix}max": float(values.max()),
    }


def scenario_summary(
    matrices: dict[str, np.ndarray],
    *,
    n_years: int,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for scenario, matrix in matrices.items():
        n = matrix.shape[0]
        rbar = _mean_off_diagonal(matrix)
        eigenvalues = np.linalg.eigvalsh(matrix)
        off_diagonal = matrix[np.triu_indices(n, k=1)]
        full_weights = np.full(n, 1.0 / n)
        loadings = _scenario_loadings(scenario, matrix)
        if np.isnan(loadings).all():
            latent_model = "not_defined_as_one_common_factor"
            latent_reliability = math.nan
            mean_loading = math.nan
            loading_sd = math.nan
        else:
            latent_model = (
                "parallel_one_factor"
                if scenario == "exchangeable"
                else "congeneric_one_factor"
            )
            signal_variance = float(np.square(full_weights @ loadings))
            residual_variance = float(
                np.sum(np.square(full_weights) * (1.0 - np.square(loadings)))
            )
            latent_reliability = signal_variance / (
                signal_variance + residual_variance
            )
            mean_loading = float(loadings.mean())
            loading_sd = float(loadings.std(ddof=0))
        plug_in_eps = eps_value(rbar, n)
        rows.append(
            {
                "scenario": scenario,
                "n_trees": n,
                "n_years": n_years,
                "mean_interseries_correlation": rbar,
                "eps": plug_in_eps,
                "latent_signal_model": latent_model,
                "latent_mean_reliability": latent_reliability,
                "eps_minus_latent_reliability": (
                    plug_in_eps - latent_reliability
                    if math.isfinite(latent_reliability)
                    else math.nan
                ),
                "mean_loading": mean_loading,
                "loading_sd": loading_sd,
                "effective_dimension": effective_dimension(matrix),
                "leading_eigenvalue_fraction": float(eigenvalues[-1] / n),
                "second_eigenvalue_fraction": float(eigenvalues[-2] / n),
                "minimum_eigenvalue": float(eigenvalues[0]),
                "off_diagonal_sd": float(off_diagonal.std(ddof=0)),
                "off_diagonal_min": float(off_diagonal.min()),
                "off_diagonal_max": float(off_diagonal.max()),
                "mean_chronology_variance": float(
                    full_weights @ matrix @ full_weights
                ),
                "matrix_sha256_canonical": _matrix_sha256(matrix),
            }
        )
    return pd.DataFrame(rows)


def exhaustive_subset_summary(
    matrices: dict[str, np.ndarray],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    n = next(iter(matrices.values())).shape[0]
    full_weights = np.full(n, 1.0 / n)
    for scenario, matrix in matrices.items():
        rbar = _mean_off_diagonal(matrix)
        full_eps = eps_value(rbar, n)
        for subset_size in range(1, n):
            values = np.fromiter(
                (
                    _weighted_correlation(
                        matrix,
                        _weight_vector(indices, n),
                        full_weights,
                    )
                    for indices in combinations(range(n), subset_size)
                ),
                dtype=float,
                count=math.comb(n, subset_size),
            )
            theoretical_sss = eps_value(rbar, subset_size) / full_eps
            row: dict[str, Any] = {
                "scenario": scenario,
                "subset_size": subset_size,
                "theoretical_exchangeable_sss": theoretical_sss,
                "theoretical_exchangeable_correlation": math.sqrt(
                    theoretical_sss
                ),
                "mean_squared_correlation": float(np.mean(values * values)),
            }
            row.update(_distribution_summary(values, prefix="correlation_"))
            rows.append(row)
    return pd.DataFrame(rows)


def exhaustive_balanced_split_summary(
    matrices: dict[str, np.ndarray],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    n = next(iter(matrices.values())).shape[0]
    if n % 2:
        raise ValueError("Stage A requires an even number of units")
    half = n // 2
    all_units = set(range(n))
    # T01 is fixed in A so A/B label reversal is counted only once.
    partitions = list(combinations(range(1, n), half - 1))
    for scenario, matrix in matrices.items():
        values = np.empty(len(partitions), dtype=float)
        for index, rest in enumerate(partitions):
            group_a = (0,) + rest
            group_b = tuple(sorted(all_units.difference(group_a)))
            values[index] = _weighted_correlation(
                matrix,
                _weight_vector(group_a, n),
                _weight_vector(group_b, n),
            )
        rbar = _mean_off_diagonal(matrix)
        row: dict[str, Any] = {
            "scenario": scenario,
            "group_size": half,
            "theoretical_exchangeable_half_eps": eps_value(rbar, half),
        }
        row.update(_distribution_summary(values, prefix="correlation_"))
        rows.append(row)
    return pd.DataFrame(rows)


def generate_bootstrap_multiplicities(
    n_units: int,
    *,
    n_replicates: int = DEFAULT_BOOTSTRAP_REPLICATES,
    seed: int = DEFAULT_BOOTSTRAP_SEED,
) -> tuple[np.ndarray, pd.DataFrame]:
    if n_replicates < 1:
        raise ValueError("n_replicates must be positive")
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, n_units, size=(n_replicates, n_units))
    counts = np.stack(
        [np.bincount(row, minlength=n_units) for row in draws], axis=0
    )
    rows: list[dict[str, Any]] = []
    for replicate_index, count in enumerate(counts, start=1):
        replicate_id = f"B{replicate_index:04d}"
        for unit_index, multiplicity in enumerate(count, start=1):
            if multiplicity:
                rows.append(
                    {
                        "replicate_id": replicate_id,
                        "unit_id": f"T{unit_index:02d}",
                        "multiplicity": int(multiplicity),
                    }
                )
    return counts, pd.DataFrame(rows)


def bootstrap_results(
    matrices: dict[str, np.ndarray],
    counts: np.ndarray,
    *,
    n_years: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    metadata = scenario_metadata()
    raw_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    driver_rows: list[dict[str, Any]] = []
    n = counts.shape[1]
    full_weights = np.full(n, 1.0 / n)

    for scenario, matrix in matrices.items():
        scenario_raw: list[dict[str, Any]] = []
        info = metadata[scenario]
        loadings = _scenario_loadings(scenario, matrix)
        for replicate_index, count in enumerate(counts, start=1):
            weights = count.astype(float) / int(count.sum())
            metrics = _comparison_metrics(
                matrix, weights, full_weights, n_years=n_years
            )
            nonzero = count[count > 0]
            unique_count = int(nonzero.size)
            kish = float(count.sum() ** 2 / np.square(nonzero).sum())
            row: dict[str, Any] = {
                "scenario": scenario,
                "replicate_id": f"B{replicate_index:04d}",
                "nominal_draw_count": int(count.sum()),
                "unique_unit_count": unique_count,
                "omitted_unit_count": n - unique_count,
                "duplicated_unit_count": int(np.sum(count > 1)),
                "maximum_multiplicity": int(count.max()),
                "kish_effective_unit_count": kish,
                **metrics,
                "cluster_A_draw_count": math.nan,
                "absolute_cluster_imbalance_fraction": math.nan,
                "mean_loading_drawn": math.nan,
                "high_loading_draw_count": math.nan,
            }
            if scenario == "two_clusters":
                cluster_a = int(count[:9].sum())
                row["cluster_A_draw_count"] = cluster_a
                row["absolute_cluster_imbalance_fraction"] = abs(
                    2 * cluster_a - n
                ) / n
            elif scenario == "unequal_loadings":
                row["mean_loading_drawn"] = float(count @ loadings / n)
            elif scenario == "influential_trees":
                row["high_loading_draw_count"] = int(count[:3].sum())
                row["mean_loading_drawn"] = float(count @ loadings / n)
            elif scenario == "exchangeable":
                row["mean_loading_drawn"] = float(loadings[0])
            scenario_raw.append(row)
            raw_rows.append(row)

        frame = pd.DataFrame(scenario_raw)
        summary: dict[str, Any] = {
            "scenario": scenario,
            "n_replicates": len(frame),
        }
        for column in (
            "correlation",
            "standardized_rms_difference",
            "variance_ratio",
            "unique_unit_count",
            "kish_effective_unit_count",
        ):
            summary.update(
                _distribution_summary(
                    frame[column].to_numpy(dtype=float),
                    prefix=f"{column}_",
                )
            )
        summary_rows.append(summary)

        if scenario == "two_clusters":
            axis = "absolute_cluster_imbalance_fraction"
            label = "absolute cluster draw imbalance"
        elif scenario == "unequal_loadings":
            axis = "mean_loading_drawn"
            label = "mean loading of bootstrap draws"
        elif scenario == "influential_trees":
            axis = "high_loading_draw_count"
            label = "number of high-loading draws"
        else:
            axis = "kish_effective_unit_count"
            label = "Kish effective unit count"

        x = frame[axis].to_numpy(dtype=float)
        driver_rows.append(
            {
                "scenario": scenario,
                "composition_axis": axis,
                "composition_axis_label": label,
                "pearson_with_standardized_rms_difference": float(
                    np.corrcoef(
                        x,
                        frame["standardized_rms_difference"].to_numpy(
                            dtype=float
                        ),
                    )[0, 1]
                ),
                "pearson_with_variance_ratio": float(
                    np.corrcoef(
                        x, frame["variance_ratio"].to_numpy(dtype=float)
                    )[0, 1]
                ),
            }
        )

    return (
        pd.DataFrame(raw_rows),
        pd.DataFrame(summary_rows),
        pd.DataFrame(driver_rows),
    )


def _tied_extreme_unit_ids(
    frame: pd.DataFrame,
    value_column: str,
    *,
    absolute: bool = False,
    tolerance: float = SUMMARY_TIE_TOLERANCE,
) -> list[str]:
    """Return every unit tied for a replicate-level extremum.

    Exact symmetric scenarios can produce values that are mathematically equal
    but differ by a few floating-point ulps across BLAS/platform combinations.
    Selecting a raw ``idxmax`` would therefore assign an arbitrary biological
    unit to the summary.  Units are treated as tied when their score is within
    the declared absolute, scale-aware numerical tolerance of the extremum.
    """
    if value_column not in frame.columns:
        raise ValueError(f"missing extremum column: {value_column}")
    values = frame[value_column].to_numpy(dtype=float)
    if values.ndim != 1 or values.size == 0 or not np.isfinite(values).all():
        raise ValueError(
            f"{value_column} must contain one or more finite values"
        )
    scores = np.abs(values) if absolute else values
    target = float(scores.max())
    absolute_tolerance = tolerance * max(1.0, abs(target))
    tied = np.isclose(
        scores,
        target,
        rtol=0.0,
        atol=absolute_tolerance,
    )
    unit_ids = sorted(
        frame.loc[tied, "omitted_unit_id"].astype(str).tolist()
    )
    if not unit_ids:
        raise RuntimeError(
            f"no tied units were identified for {value_column}"
        )
    return unit_ids


def loto_results(
    matrices: dict[str, np.ndarray],
    *,
    n_years: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    metadata = scenario_metadata()
    raw_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    for scenario, matrix in matrices.items():
        n = matrix.shape[0]
        info = metadata[scenario]
        loadings = _scenario_loadings(scenario, matrix)
        full_weights = np.full(n, 1.0 / n)
        full_rbar = _mean_off_diagonal(matrix)
        full_eps = eps_value(full_rbar, n)
        full_dimension = effective_dimension(matrix)
        scenario_rows: list[dict[str, Any]] = []

        for omitted in range(n):
            retained = [index for index in range(n) if index != omitted]
            reduced_weights = _weight_vector(retained, n)
            metrics = _comparison_metrics(
                matrix, reduced_weights, full_weights, n_years=n_years
            )
            reduced_matrix = matrix[np.ix_(retained, retained)]
            reduced_rbar = _mean_off_diagonal(reduced_matrix)
            reduced_eps = eps_value(reduced_rbar, n - 1)
            reduced_dimension = effective_dimension(reduced_matrix)
            row = {
                "scenario": scenario,
                "omitted_unit_id": f"T{omitted + 1:02d}",
                "unit_class": info.unit_class[omitted],
                "loading": (
                    None if np.isnan(loadings[omitted]) else float(loadings[omitted])
                ),
                **metrics,
                "full_rbar": full_rbar,
                "reduced_rbar": reduced_rbar,
                "delta_rbar": reduced_rbar - full_rbar,
                "full_eps": full_eps,
                "reduced_eps": reduced_eps,
                "delta_eps": reduced_eps - full_eps,
                "full_effective_dimension": full_dimension,
                "reduced_effective_dimension": reduced_dimension,
                "delta_effective_dimension": (
                    reduced_dimension - full_dimension
                ),
            }
            scenario_rows.append(row)
            raw_rows.append(row)

        frame = pd.DataFrame(scenario_rows)
        max_rms_unit_ids = _tied_extreme_unit_ids(
            frame,
            "standardized_rms_difference",
        )
        max_abs_eps_unit_ids = _tied_extreme_unit_ids(
            frame,
            "delta_eps",
            absolute=True,
        )
        max_abs_dimension_unit_ids = _tied_extreme_unit_ids(
            frame,
            "delta_effective_dimension",
            absolute=True,
        )
        summary_rows.append(
            {
                "scenario": scenario,
                "n_omissions": len(frame),
                "correlation_min": float(frame["correlation"].min()),
                "correlation_max": float(frame["correlation"].max()),
                "standardized_rms_min": float(
                    frame["standardized_rms_difference"].min()
                ),
                "standardized_rms_max": float(
                    frame["standardized_rms_difference"].max()
                ),
                "standardized_rms_sd": float(
                    frame["standardized_rms_difference"].std(ddof=1)
                ),
                "delta_rbar_min": float(frame["delta_rbar"].min()),
                "delta_rbar_max": float(frame["delta_rbar"].max()),
                "delta_eps_min": float(frame["delta_eps"].min()),
                "delta_eps_max": float(frame["delta_eps"].max()),
                "delta_effective_dimension_min": float(
                    frame["delta_effective_dimension"].min()
                ),
                "delta_effective_dimension_max": float(
                    frame["delta_effective_dimension"].max()
                ),
                # Singular fields retain a deterministic representative for
                # compatibility; plural fields preserve the scientifically
                # relevant tied set.
                "max_rms_unit_id": max_rms_unit_ids[0],
                "max_rms_unit_ids": "|".join(max_rms_unit_ids),
                "max_rms_tie_count": len(max_rms_unit_ids),
                "max_abs_delta_eps_unit_id": max_abs_eps_unit_ids[0],
                "max_abs_delta_eps_unit_ids": "|".join(
                    max_abs_eps_unit_ids
                ),
                "max_abs_delta_eps_tie_count": len(
                    max_abs_eps_unit_ids
                ),
                "max_abs_delta_effective_dimension_unit_id": (
                    max_abs_dimension_unit_ids[0]
                ),
                "max_abs_delta_effective_dimension_unit_ids": "|".join(
                    max_abs_dimension_unit_ids
                ),
                "max_abs_delta_effective_dimension_tie_count": len(
                    max_abs_dimension_unit_ids
                ),
            }
        )
    return pd.DataFrame(raw_rows), pd.DataFrame(summary_rows)


def exact_long_data(
    matrices: dict[str, np.ndarray],
    *,
    n_years: int,
    start_year: int,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    years = np.arange(start_year, start_year + n_years)
    for scenario, matrix in matrices.items():
        values = exact_sample_from_correlation(matrix, n_years=n_years)
        for unit_index in range(matrix.shape[0]):
            unit = f"T{unit_index + 1:02d}"
            rows.extend(
                {
                    "scenario": scenario,
                    "year": int(year),
                    "series_id": unit,
                    "tree_id": unit,
                    "value": float(value),
                }
                for year, value in zip(
                    years, values[:, unit_index], strict=True
                )
            )
    return pd.DataFrame(rows)


def key_results_table(
    scenario: pd.DataFrame,
    subset: pd.DataFrame,
    split: pd.DataFrame,
    bootstrap: pd.DataFrame,
    loto: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for scenario_name in SCENARIOS:
        structure = scenario.loc[scenario["scenario"] == scenario_name].iloc[0]
        subset_nine = subset.loc[
            (subset["scenario"] == scenario_name)
            & (subset["subset_size"] == 9)
        ].iloc[0]
        split_row = split.loc[split["scenario"] == scenario_name].iloc[0]
        bootstrap_row = bootstrap.loc[
            bootstrap["scenario"] == scenario_name
        ].iloc[0]
        loto_row = loto.loc[loto["scenario"] == scenario_name].iloc[0]
        rows.append(
            {
                "scenario": scenario_name,
                "n_trees": int(structure["n_trees"]),
                "rbar": float(structure["mean_interseries_correlation"]),
                "eps": float(structure["eps"]),
                "effective_dimension": float(
                    structure["effective_dimension"]
                ),
                "latent_mean_reliability": float(
                    structure["latent_mean_reliability"]
                ),
                "eps_minus_latent_reliability": float(
                    structure["eps_minus_latent_reliability"]
                ),
                "off_diagonal_sd": float(structure["off_diagonal_sd"]),
                "subset_n9_theoretical_correlation": float(
                    subset_nine["theoretical_exchangeable_correlation"]
                ),
                "subset_n9_q05": float(subset_nine["correlation_q05"]),
                "subset_n9_median": float(
                    subset_nine["correlation_median"]
                ),
                "subset_n9_q95": float(subset_nine["correlation_q95"]),
                "subset_n9_min": float(subset_nine["correlation_min"]),
                "split_half_expected_correlation": float(
                    split_row["theoretical_exchangeable_half_eps"]
                ),
                "split_half_q05": float(split_row["correlation_q05"]),
                "split_half_median": float(
                    split_row["correlation_median"]
                ),
                "split_half_q95": float(split_row["correlation_q95"]),
                "split_half_min": float(split_row["correlation_min"]),
                "bootstrap_srms_q95": float(
                    bootstrap_row[
                        "standardized_rms_difference_q95"
                    ]
                ),
                "bootstrap_variance_ratio_q05": float(
                    bootstrap_row["variance_ratio_q05"]
                ),
                "bootstrap_variance_ratio_q95": float(
                    bootstrap_row["variance_ratio_q95"]
                ),
                "loto_correlation_min": float(
                    loto_row["correlation_min"]
                ),
                "loto_standardized_rms_max": float(
                    loto_row["standardized_rms_max"]
                ),
                "loto_delta_eps_min": float(loto_row["delta_eps_min"]),
                "loto_delta_eps_max": float(loto_row["delta_eps_max"]),
                "loto_delta_effective_dimension_min": float(
                    loto_row["delta_effective_dimension_min"]
                ),
                "loto_delta_effective_dimension_max": float(
                    loto_row["delta_effective_dimension_max"]
                ),
            }
        )
    return pd.DataFrame(rows)


def write_table(frame: pd.DataFrame, path: Path) -> None:
    """Write a deterministic UTF-8 CSV, optionally in reproducible gzip form."""
    if path.suffix == ".gz":
        with path.open("wb") as binary:
            with gzip.GzipFile(
                filename="",
                mode="wb",
                fileobj=binary,
                compresslevel=9,
                mtime=0,
            ) as compressed:
                with io.TextIOWrapper(
                    compressed, encoding="utf-8", newline=""
                ) as text_stream:
                    frame.to_csv(
                        text_stream,
                        index=False,
                        float_format="%.15g",
                        lineterminator="\n",
                    )
    else:
        frame.to_csv(
            path,
            index=False,
            float_format="%.15g",
            lineterminator="\n",
        )


def run_stage_a(
    *,
    fixture_dir: Path,
    output_dir: Path,
    n_years: int = DEFAULT_N_YEARS,
    start_year: int = DEFAULT_START_YEAR,
    bootstrap_replicates: int = DEFAULT_BOOTSTRAP_REPLICATES,
    bootstrap_seed: int = DEFAULT_BOOTSTRAP_SEED,
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    matrices = load_scenario_matrices(fixture_dir)
    scenario = scenario_summary(matrices, n_years=n_years)
    subset = exhaustive_subset_summary(matrices)
    split = exhaustive_balanced_split_summary(matrices)
    counts, multiplicities = generate_bootstrap_multiplicities(
        next(iter(matrices.values())).shape[0],
        n_replicates=bootstrap_replicates,
        seed=bootstrap_seed,
    )
    bootstrap_raw, bootstrap_summary, bootstrap_driver = bootstrap_results(
        matrices, counts, n_years=n_years
    )
    loto_raw, loto_summary = loto_results(matrices, n_years=n_years)
    exact = exact_long_data(
        matrices, n_years=n_years, start_year=start_year
    )
    key = key_results_table(
        scenario,
        subset,
        split,
        bootstrap_summary,
        loto_summary,
    )

    outputs = {
        "scenario_summary": output_dir / "scenario_summary_v01.csv",
        "subset_summary": output_dir / "subset_summary_v01.csv",
        "split_half_summary": output_dir / "split_half_summary_v01.csv",
        "bootstrap_multiplicities": (
            output_dir / "bootstrap_multiplicities_v01.csv.gz"
        ),
        "bootstrap_raw": output_dir / "bootstrap_results_v01.csv.gz",
        "bootstrap_summary": output_dir / "bootstrap_summary_v01.csv",
        "bootstrap_driver": output_dir / "bootstrap_driver_summary_v01.csv",
        "loto_raw": output_dir / "loto_results_v01.csv",
        "loto_summary": output_dir / "loto_summary_v01.csv",
        "exact_long_data": output_dir / "exact_same_eps_long_v01.csv.gz",
        "key_results": output_dir / "key_results_v01.csv",
    }
    tables = {
        "scenario_summary": scenario,
        "subset_summary": subset,
        "split_half_summary": split,
        "bootstrap_multiplicities": multiplicities,
        "bootstrap_raw": bootstrap_raw,
        "bootstrap_summary": bootstrap_summary,
        "bootstrap_driver": bootstrap_driver,
        "loto_raw": loto_raw,
        "loto_summary": loto_summary,
        "exact_long_data": exact,
        "key_results": key,
    }
    for name, path in outputs.items():
        write_table(tables[name], path)

    manifest = {
        "analysis": "same_eps_stage_a",
        "version": VERSION,
        "experiment_type": "deterministic controlled covariance experiment",
        "n_trees": int(next(iter(matrices.values())).shape[0]),
        "n_years": n_years,
        "start_year": start_year,
        "bootstrap_replicates": bootstrap_replicates,
        "bootstrap_seed": bootstrap_seed,
        "bootstrap_rng": "numpy.random.default_rng (PCG64)",
        "subset_design": "exhaustive for every subset size 1..N-1",
        "subset_composition_count_total": int(2 ** next(iter(matrices.values())).shape[0] - 2),
        "split_half_design": (
            "exhaustive unique balanced partitions; T01 fixed in group A"
        ),
        "split_half_partition_count": int(
            math.comb(next(iter(matrices.values())).shape[0] - 1,
                      next(iter(matrices.values())).shape[0] // 2 - 1)
        ),
        "chronology_method": "unweighted arithmetic mean",
        "extreme_unit_tie_policy": (
            "all units within SUMMARY_TIE_TOLERANCE of the extremum; "
            "singular representative is the lexicographically first tied ID"
        ),
        "extreme_unit_tie_tolerance": SUMMARY_TIE_TOLERANCE,
        "input_matrix_files": {
            scenario_name: {
                "path": _manifest_path(
                    fixture_dir / f"{scenario_name}_correlation.csv"
                ),
                "sha256_bytes": _file_sha256(
                    fixture_dir / f"{scenario_name}_correlation.csv"
                ),
                "sha256_canonical_numeric": _matrix_sha256(matrix),
            }
            for scenario_name, matrix in matrices.items()
        },
        "output_files": {
            name: {
                "path": _manifest_path(path),
                "sha256": _file_sha256(path),
            }
            for name, path in outputs.items()
        },
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
    }
    manifest_path = output_dir / "analysis_manifest_v01.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    outputs["manifest"] = manifest_path
    return outputs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fixture-dir",
        type=Path,
        default=REPO_ROOT / "shared_fixtures",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=(
            REPO_ROOT
            / "analysis"
            / "same_eps_stage_a_v01"
        ),
    )
    parser.add_argument("--n-years", type=int, default=DEFAULT_N_YEARS)
    parser.add_argument(
        "--start-year", type=int, default=DEFAULT_START_YEAR
    )
    parser.add_argument(
        "--bootstrap-replicates",
        type=int,
        default=DEFAULT_BOOTSTRAP_REPLICATES,
    )
    parser.add_argument(
        "--bootstrap-seed",
        type=int,
        default=DEFAULT_BOOTSTRAP_SEED,
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = run_stage_a(
        fixture_dir=args.fixture_dir,
        output_dir=args.output_dir,
        n_years=args.n_years,
        start_year=args.start_year,
        bootstrap_replicates=args.bootstrap_replicates,
        bootstrap_seed=args.bootstrap_seed,
    )
    print("same-EPS Stage A completed.")
    for name, path in outputs.items():
        print(f"{name}: {path}")


if __name__ == "__main__":
    main()
