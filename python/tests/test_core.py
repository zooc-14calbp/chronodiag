from __future__ import annotations

import numpy as np
import pytest

from chronodiag import effective_dimension, eps_value, sss_value


def test_eps_reference_value() -> None:
    assert eps_value(0.23943661971830986, 18) == pytest.approx(0.85, abs=1e-12)


def test_sss_zero_correlation_is_defined() -> None:
    assert sss_value(0.0, 5, 20) == pytest.approx(0.25, abs=1e-15)


def test_sss_full_subset_is_one() -> None:
    assert sss_value(0.3, 20, 20) == pytest.approx(1.0, abs=1e-15)


def test_eps_rejects_inadmissible_equicorrelation() -> None:
    with pytest.raises(ValueError, match="admissible"):
        eps_value(-0.2, 10)


def test_effective_dimension_identity() -> None:
    assert effective_dimension(np.eye(5)) == pytest.approx(5.0)


def test_effective_dimension_rejects_non_psd() -> None:
    matrix = np.array([[1.0, 1.2], [1.2, 1.0]])
    with pytest.raises(ValueError, match="positive semidefinite"):
        effective_dimension(matrix)
