from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from chronodiag import split_half_stability, subset_stability

FIXTURES = Path(__file__).resolve().parents[2] / "shared_fixtures"
_RESULT_NUMERIC = [
    "min_depth_fraction",
    "min_overlap",
    "full_pool_size",
    "group_a_size",
    "group_b_size",
    "shared_unit_count",
    "correlation",
    "squared_correlation",
    "n_years",
    "start_year",
    "end_year",
    "mean_depth_a",
    "min_depth_a",
    "mean_depth_b",
    "min_depth_b",
]


def _assert_frame_matches(result: pd.DataFrame, expected_name: str) -> None:
    expected = pd.read_csv(FIXTURES / expected_name)
    assert list(result.columns) == list(expected.columns)
    assert len(result) == len(expected)
    for column in _RESULT_NUMERIC:
        assert np.allclose(
            pd.to_numeric(result[column]).to_numpy(dtype=float),
            pd.to_numeric(expected[column]).to_numpy(dtype=float),
            atol=1e-12,
            rtol=0.0,
            equal_nan=True,
        )
    text_columns = [column for column in result.columns if column not in _RESULT_NUMERIC]
    pd.testing.assert_frame_equal(
        result[text_columns].reset_index(drop=True),
        expected[text_columns].reset_index(drop=True),
        check_dtype=False,
    )


def test_subset_stability_matches_frozen_fixture() -> None:
    data = pd.read_csv(FIXTURES / "composition_example_long_v01.csv")
    subsets = pd.read_csv(FIXTURES / "subset_partitions_v01.csv")
    result = subset_stability(
        data,
        subsets,
        min_overlap=10,
        min_depth_fraction=1.0,
    )
    _assert_frame_matches(result, "expected_subset_stability_v01.csv")
    assert (result["shared_unit_count"] == result["group_a_size"]).all()
    assert set(result["correlation_method"]) == {"pearson"}


def test_split_half_stability_matches_frozen_fixture() -> None:
    data = pd.read_csv(FIXTURES / "composition_example_long_v01.csv")
    partitions = pd.read_csv(FIXTURES / "split_half_partitions_v01.csv")
    result = split_half_stability(
        data,
        partitions,
        min_overlap=10,
        min_depth_fraction=1.0,
    )
    _assert_frame_matches(result, "expected_split_half_stability_v01.csv")
    assert (result["shared_unit_count"] == 0).all()
    assert set(result["correlation_method"]) == {"pearson"}


def test_depth_fraction_is_explicit_and_changes_available_years() -> None:
    data = pd.read_csv(FIXTURES / "composition_example_long_v01.csv")
    subsets = pd.DataFrame(
        {"replicate_id": ["S", "S"], "unit_id": ["T03", "T04"]}
    )
    complete = subset_stability(data, subsets, min_depth_fraction=1.0)
    permissive = subset_stability(data, subsets, min_depth_fraction=0.5)
    assert int(complete.loc[0, "n_years"]) == 21
    assert int(permissive.loc[0, "n_years"]) == 24
    assert float(permissive.loc[0, "mean_depth_b"]) < 8.0


def test_split_half_rejects_overlapping_or_non_exhaustive_partitions() -> None:
    data = pd.read_csv(FIXTURES / "composition_example_long_v01.csv")
    overlapping = pd.DataFrame(
        {
            "replicate_id": ["H", "H", "H"],
            "group": ["A", "B", "B"],
            "unit_id": ["T01", "T01", "T02"],
        }
    )
    with pytest.raises(ValueError, match="not disjoint"):
        split_half_stability(data, overlapping, require_exhaustive=False)

    partial = pd.DataFrame(
        {
            "replicate_id": ["H", "H"],
            "group": ["A", "B"],
            "unit_id": ["T01", "T02"],
        }
    )
    with pytest.raises(ValueError, match="not exhaustive"):
        split_half_stability(data, partial)


def test_duplicate_year_unit_and_all_missing_unit_are_rejected() -> None:
    data = pd.read_csv(FIXTURES / "composition_example_long_v01.csv")
    duplicate = pd.concat([data, data.iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError, match="one value per year"):
        subset_stability(duplicate, {"S": "T01"})

    empty = pd.concat(
        [
            data,
            pd.DataFrame(
                {
                    "year": [2000, 2001],
                    "series_id": ["EMPTY", "EMPTY"],
                    "tree_id": ["EMPTY", "EMPTY"],
                    "value": [np.nan, np.nan],
                }
            ),
        ],
        ignore_index=True,
    )
    with pytest.raises(ValueError, match="no finite observations"):
        subset_stability(empty, {"S": "T01"})


def test_status_rows_are_retained_for_insufficient_and_constant_comparisons() -> None:
    insufficient_data = pd.DataFrame(
        {
            "year": [2000, 2001, 2002, 2000, 2001, 2002],
            "tree_id": ["A", "A", "A", "B", "B", "B"],
            "value": [0.0, 1.0, 0.0, 0.2, 1.2, 0.2],
        }
    )
    insufficient = subset_stability(
        insufficient_data,
        {"short": "A"},
        min_overlap=4,
    )
    assert insufficient.loc[0, "status"] == "insufficient_overlap"
    assert np.isnan(float(insufficient.loc[0, "correlation"]))

    constant_data = pd.DataFrame(
        {
            "year": [2000, 2001, 2002, 2000, 2001, 2002],
            "tree_id": ["A", "A", "A", "B", "B", "B"],
            "value": [1.0, 1.0, 1.0, 0.0, 1.0, 2.0],
        }
    )
    constant = subset_stability(constant_data, {"constant": "A"})
    assert constant.loc[0, "status"] == "constant_chronology"
    assert np.isnan(float(constant.loc[0, "squared_correlation"]))


def test_all_unit_subset_is_identical_to_full_chronology() -> None:
    data = pd.read_csv(FIXTURES / "composition_example_long_v01.csv")
    units = sorted(data["tree_id"].unique())
    result = subset_stability(data, {"all": units}, min_overlap=10)
    assert result.loc[0, "correlation"] == pytest.approx(1.0, abs=1e-14)
    assert result.loc[0, "squared_correlation"] == pytest.approx(1.0, abs=1e-14)
    assert bool(result.loc[0, "comparison_exhaustive"])
    assert int(result.loc[0, "full_pool_size"]) == len(units)


def test_partial_disjoint_comparison_is_reported_as_non_exhaustive() -> None:
    data = pd.read_csv(FIXTURES / "composition_example_long_v01.csv")
    partitions = pd.DataFrame(
        {
            "replicate_id": ["P", "P"],
            "group": ["A", "B"],
            "unit_id": ["T01", "T02"],
        }
    )
    result = split_half_stability(
        data, partitions, require_exhaustive=False, min_overlap=10
    )
    assert not bool(result.loc[0, "comparison_exhaustive"])
    assert int(result.loc[0, "full_pool_size"]) == 8
    assert int(result.loc[0, "shared_unit_count"]) == 0


def test_split_half_requires_balanced_groups_by_default() -> None:
    data = pd.read_csv(FIXTURES / "composition_example_long_v01.csv")
    partitions = pd.DataFrame(
        {
            "replicate_id": ["U"] * 8,
            "group": ["A", "A", "B", "B", "B", "B", "B", "B"],
            "unit_id": [f"T{i:02d}" for i in range(1, 9)],
        }
    )
    with pytest.raises(ValueError, match="not balanced"):
        split_half_stability(data, partitions)

    result = split_half_stability(
        data, partitions, require_balanced=False, min_overlap=10
    )
    assert int(result.loc[0, "group_a_size"]) == 2
    assert int(result.loc[0, "group_b_size"]) == 6
    assert bool(result.loc[0, "comparison_exhaustive"])


def test_exchangeable_fixture_links_subset_and_split_half_to_eps_sss() -> None:
    from chronodiag import eps_value, sss_value

    data = pd.read_csv(FIXTURES / "composition_exchangeable_exact_v01.csv")
    subset = subset_stability(
        data,
        {"S": ["E01", "E02", "E03"]},
        min_overlap=20,
    )
    assert subset.loc[0, "squared_correlation"] == pytest.approx(
        sss_value(0.3, 3, 6), abs=1e-12
    )

    partitions = pd.DataFrame(
        {
            "replicate_id": ["H"] * 6,
            "group": ["A", "A", "A", "B", "B", "B"],
            "unit_id": ["E01", "E02", "E03", "E04", "E05", "E06"],
        }
    )
    split = split_half_stability(data, partitions, min_overlap=20)
    assert split.loc[0, "correlation"] == pytest.approx(
        eps_value(0.3, 3), abs=1e-12
    )
