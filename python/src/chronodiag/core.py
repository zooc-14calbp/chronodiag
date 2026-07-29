from __future__ import annotations

import math
from typing import Any

import numpy as np


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


def _finite_correlation(value: Any, name: str = "rbar") -> float:
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise TypeError(f"{name} must be a finite correlation") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    if result < -1.0 or result > 1.0:
        raise ValueError(f"{name} must lie in [-1, 1]")
    return result


def eps_value(rbar: float, n: int) -> float:
    """Return model-relative EPS for mean interseries correlation ``rbar``.

    Negative values are not clipped.  A non-positive denominator indicates that
    the equicorrelation model is inadmissible for the requested ``n``.
    """
    rho = _finite_correlation(rbar)
    n_int = _positive_integer(n, "n")
    denominator = 1.0 + (n_int - 1) * rho
    if denominator <= 0.0 or math.isclose(
        denominator, 0.0, abs_tol=np.finfo(float).eps, rel_tol=0.0
    ):
        raise ValueError(
            "rbar is outside the admissible equicorrelation range for n"
        )
    return n_int * rho / denominator


def sss_value(rbar: float, n_sub: int, n_full: int) -> float:
    """Return SSS for a subset of ``n_sub`` within a full set of ``n_full``.

    The direct formula is used rather than EPS(n_sub) / EPS(n_full), so the
    well-defined zero-correlation case returns ``n_sub / n_full``.
    """
    rho = _finite_correlation(rbar)
    n_sub_int = _positive_integer(n_sub, "n_sub")
    n_full_int = _positive_integer(n_full, "n_full")
    if n_sub_int > n_full_int:
        raise ValueError("n_sub must not exceed n_full")
    full_term = 1.0 + (n_full_int - 1) * rho
    sub_term = 1.0 + (n_sub_int - 1) * rho
    if full_term <= 0.0 or sub_term <= 0.0:
        raise ValueError(
            "rbar is outside the admissible equicorrelation range"
        )
    return n_sub_int * full_term / (n_full_int * sub_term)


def effective_dimension(
    correlation: np.ndarray, tolerance: float = 1e-10
) -> float:
    """Return the eigenvalue/Frobenius effective dimension of a correlation matrix."""
    if not math.isfinite(float(tolerance)) or tolerance < 0:
        raise ValueError("tolerance must be finite and non-negative")
    matrix = np.asarray(correlation, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1] or matrix.shape[0] < 1:
        raise ValueError("correlation must be a non-empty square matrix")
    if not np.isfinite(matrix).all():
        raise ValueError("correlation contains non-finite values")
    if not np.allclose(matrix, matrix.T, atol=tolerance, rtol=0.0):
        raise ValueError("correlation must be symmetric")
    if not np.allclose(np.diag(matrix), 1.0, atol=tolerance, rtol=0.0):
        raise ValueError("correlation diagonal must equal one")
    eigenvalues = np.linalg.eigvalsh(matrix)
    if float(eigenvalues.min()) < -tolerance:
        raise ValueError("correlation must be positive semidefinite")
    trace = float(np.trace(matrix))
    denominator = float(np.square(matrix).sum())
    if denominator <= 0:
        raise ValueError("correlation has zero Frobenius norm")
    result = trace * trace / denominator
    n = matrix.shape[0]
    if result < 1.0 - tolerance or result > n + tolerance:
        raise RuntimeError("effective dimension fell outside [1, n]")
    return min(float(n), max(1.0, result))
