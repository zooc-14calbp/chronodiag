from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from chronodiag import effective_dimension, eps_value

FIXTURES = Path(__file__).resolve().parents[2] / "shared_fixtures"


def _mean_off_diagonal(matrix: np.ndarray) -> float:
    n = matrix.shape[0]
    return float((matrix.sum() - np.trace(matrix)) / (n * (n - 1)))


@pytest.mark.parametrize(
    "scenario",
    ["exchangeable", "two_clusters", "unequal_loadings", "influential_trees"],
)
def test_same_eps_different_structure_fixtures(scenario: str) -> None:
    expected = pd.read_csv(FIXTURES / "expected_summaries.csv").set_index("scenario")
    matrix = np.loadtxt(FIXTURES / f"{scenario}_correlation.csv", delimiter=",")
    row = expected.loc[scenario]

    assert matrix.shape == (int(row["n_trees"]), int(row["n_trees"]))
    assert np.allclose(matrix, matrix.T, atol=1e-12, rtol=0.0)
    assert np.allclose(np.diag(matrix), 1.0, atol=1e-12, rtol=0.0)

    rbar = _mean_off_diagonal(matrix)
    assert rbar == pytest.approx(row["mean_interseries_correlation"], abs=1e-12)
    assert eps_value(rbar, matrix.shape[0]) == pytest.approx(row["eps"], abs=1e-12)
    assert effective_dimension(matrix) == pytest.approx(
        row["effective_dimension"], abs=1e-10
    )

    eigenvalues = np.linalg.eigvalsh(matrix)
    leading_fraction = float(eigenvalues[-1] / eigenvalues.sum())
    assert leading_fraction == pytest.approx(
        row["leading_eigenvalue_fraction"], abs=1e-10
    )
    assert float(eigenvalues[0]) == pytest.approx(row["minimum_eigenvalue"], abs=1e-10)
