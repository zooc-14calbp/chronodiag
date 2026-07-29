from __future__ import annotations

from pathlib import Path
import json

import numpy as np
import pandas as pd
import pytest

from chronodiag import leave_one_tree_out, tree_bootstrap

FIXTURES = Path(__file__).resolve().parents[2] / "shared_fixtures"


def _data() -> pd.DataFrame:
    return pd.read_csv(FIXTURES / "influence_example_long_v01.csv")


def _assert_expected(actual: pd.DataFrame, expected_name: str) -> None:
    expected = pd.read_csv(FIXTURES / expected_name)
    assert list(actual.columns) == list(expected.columns)
    assert actual.shape == expected.shape
    for column in actual.columns:
        if pd.api.types.is_numeric_dtype(expected[column]):
            assert np.allclose(
                pd.to_numeric(actual[column], errors="coerce").to_numpy(float),
                pd.to_numeric(expected[column], errors="coerce").to_numpy(float),
                atol=1e-12,
                rtol=0.0,
                equal_nan=True,
            ), column
        else:
            left = actual[column].fillna("<NA>").astype(str).tolist()
            right = expected[column].fillna("<NA>").astype(str).tolist()
            assert left == right, column


def test_tree_bootstrap_matches_frozen_fixture() -> None:
    multiplicities = pd.read_csv(FIXTURES / "bootstrap_multiplicities_v01.csv")
    actual = tree_bootstrap(_data(), multiplicities)
    _assert_expected(actual, "expected_tree_bootstrap_v01.csv")


def test_identity_bootstrap_has_zero_difference() -> None:
    multiplicities = pd.read_csv(FIXTURES / "bootstrap_multiplicities_v01.csv")
    result = tree_bootstrap(_data(), multiplicities).set_index("replicate_id")
    identity = result.loc["B02"]
    assert identity["all_units_represented"]
    assert identity["unique_unit_count"] == 6
    assert identity["kish_effective_unit_count"] == pytest.approx(6.0)
    assert identity["correlation"] == pytest.approx(1.0)
    assert identity["rms_difference"] == pytest.approx(0.0, abs=1e-15)


def test_bootstrap_mapping_represents_literal_draws() -> None:
    mapping = {"M01": ["T01", "T01", "T02", "T03", "T04", "T05"]}
    result = tree_bootstrap(_data(), mapping)
    row = result.iloc[0]
    assert row["bootstrap_ids"] == "T01*2|T02*1|T03*1|T04*1|T05*1"
    assert row["nominal_draw_count"] == 6
    assert row["unique_unit_count"] == 5


def test_bootstrap_requires_positive_combined_multiplicity_rows() -> None:
    duplicate = pd.DataFrame(
        {
            "replicate_id": ["B", "B"],
            "unit_id": ["T01", "T01"],
            "multiplicity": [1, 1],
        }
    )
    with pytest.raises(ValueError, match="combine"):
        tree_bootstrap(_data(), duplicate, require_full_size=False)

    zero = pd.DataFrame(
        {"replicate_id": ["B"], "unit_id": ["T01"], "multiplicity": [0]}
    )
    with pytest.raises(ValueError, match="positive"):
        tree_bootstrap(_data(), zero, require_full_size=False)


def test_bootstrap_size_contract_is_explicit() -> None:
    short = pd.DataFrame(
        {
            "replicate_id": ["B", "B"],
            "unit_id": ["T01", "T02"],
            "multiplicity": [1, 1],
        }
    )
    with pytest.raises(ValueError, match="expected 6"):
        tree_bootstrap(_data(), short)
    result = tree_bootstrap(_data(), short, require_full_size=False)
    assert result.loc[0, "nominal_draw_count"] == 2


def test_kish_count_is_bounded_by_unique_and_nominal_counts() -> None:
    multiplicities = pd.read_csv(FIXTURES / "bootstrap_multiplicities_v01.csv")
    result = tree_bootstrap(_data(), multiplicities)
    assert (result["kish_effective_unit_count"] >= 1).all()
    assert (
        result["kish_effective_unit_count"] <= result["unique_unit_count"]
    ).all()
    assert (result["unique_unit_count"] <= result["nominal_draw_count"]).all()


def test_weighted_median_expands_literal_bootstrap_draws() -> None:
    data = pd.DataFrame(
        {
            "year": [2000, 2001, 2002] * 3,
            "tree_id": ["T01"] * 3 + ["T02"] * 3 + ["T03"] * 3,
            "value": [0.0, 1.0, 2.0, 10.0, 11.0, 12.0, 20.0, 21.0, 22.0],
        }
    )
    draws = {"M": ["T01", "T01", "T02"]}
    result = tree_bootstrap(data, draws, chronology_method="median").iloc[0]
    assert result["mean_difference"] == pytest.approx(-10.0)
    assert result["rms_difference"] == pytest.approx(10.0)
    assert result["correlation"] == pytest.approx(1.0)


def test_leave_one_tree_out_matches_frozen_fixture() -> None:
    actual = leave_one_tree_out(_data())
    _assert_expected(actual, "expected_leave_one_tree_out_v01.csv")


def test_loto_identifies_influential_fixture_unit_by_direct_difference() -> None:
    result = leave_one_tree_out(_data()).set_index("omitted_unit_id")
    assert result["rms_difference"].idxmax() == "T06"
    assert result.loc["T06", "standardized_rms_difference"] > result.loc[
        "T01", "standardized_rms_difference"
    ]
    assert result.loc["T06", "shared_unit_count"] == 5
    assert result.loc["T06", "correlation"] > 0.99


def test_loto_direct_difference_matches_complete_mean_identity() -> None:
    data = _data()
    matrix = data.pivot(index="year", columns="tree_id", values="value").dropna()
    full = matrix.mean(axis=1)
    reduced = matrix.drop(columns="T06").mean(axis=1)
    identity_difference = (reduced - matrix["T06"]) / matrix.shape[1]
    observed_difference = reduced - full
    assert np.allclose(
        observed_difference.to_numpy(),
        identity_difference.to_numpy(),
        atol=1e-14,
        rtol=0.0,
    )

    result = leave_one_tree_out(data).set_index("omitted_unit_id").loc["T06"]
    assert result["mean_difference"] == pytest.approx(observed_difference.mean())
    assert result["rms_difference"] == pytest.approx(
        np.sqrt(np.mean(np.square(observed_difference)))
    )


def test_loto_structure_deltas_use_one_common_complete_period() -> None:
    result = leave_one_tree_out(_data())
    assert set(result["structure_status"]) == {"ok"}
    assert set(result["structure_n_years"]) == {21}
    assert result["full_rbar"].nunique() == 1
    assert result["full_eps"].nunique() == 1
    assert result["full_effective_dimension"].nunique() == 1


def test_loto_requires_at_least_two_units() -> None:
    one = _data().loc[lambda x: x["tree_id"] == "T01"]
    with pytest.raises(ValueError, match="at least two"):
        leave_one_tree_out(one)


def test_result_schema_is_valid_json() -> None:
    schema = json.loads(
        (FIXTURES.parents[0] / "spec" / "influence_result_schema_v01.json").read_text(
            encoding="utf-8"
        )
    )
    assert schema["title"] == "chronodiag influence diagnostic result row"
