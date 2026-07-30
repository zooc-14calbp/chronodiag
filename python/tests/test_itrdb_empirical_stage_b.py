from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]


def load_script(name: str):
    path = REPO_ROOT / "scripts" / name
    spec = importlib.util.spec_from_file_location(path.stem, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_empirical_contract_and_keymap_counts() -> None:
    contract = pd.read_csv(
        REPO_ROOT / "data" / "itrdb" / "empirical_dataset_contract_v01.csv"
    )
    expected = {
        "ca533": (17, 27, 1727, 1968),
        "germ012": (21, 42, 1873, 1944),
        "fran003": (25, 25, 1819, 1976),
    }
    assert set(contract["dataset"]) == set(expected)
    assert set(contract["rwl_format"].astype(str).str.casefold()) == {"tucson"}
    assert set(contract["tucson_long"].astype(str).str.casefold()) == {"false"}
    assert set(contract["preserve_zero"].astype(str).str.casefold()) == {"true"}
    expected_header = {
        "ca533": "false",
        "germ012": "true",
        "fran003": "true",
    }
    for row in contract.itertuples(index=False):
        assert str(row.tucson_header).casefold() == expected_header[row.dataset]
        n_trees, n_series, start, end = expected[row.dataset]
        assert int(row.n_trees_expected) == n_trees
        assert int(row.n_series_expected) == n_series
        assert int(row.analysis_start) == start
        assert int(row.analysis_end) == end
        keymap = pd.read_csv(
            REPO_ROOT
            / "data"
            / "itrdb"
            / "keymaps"
            / f"{row.dataset}_keymap_v01.csv",
            dtype={"tree_id": str, "series_id": str},
        )
        included = keymap[
            keymap["include_primary_cohort"].astype(str).str.casefold().eq("true")
        ]
        assert included["tree_id"].nunique() == n_trees
        assert included["series_id"].nunique() == n_series


def test_cohort_selection_rule_is_value_blind_and_deterministic() -> None:
    module = load_script("audit_itrdb_empirical_cohorts_v01.py")
    keymap = pd.DataFrame(
        {
            "tree_id": ["A", "B", "C", "D"],
            "series_id": ["A1", "B1", "C1", "D1"],
            "series_start": [1900, 1900, 1910, 1920],
            "series_end": [2000, 2000, 2000, 2000],
        }
    )
    best, candidates = module.select_cohort(
        keymap, minimum_years=50, minimum_years_per_tree_ratio=3.0
    )
    assert int(best["n_trees"]) == 4
    assert int(best["analysis_start"]) == 1920
    assert int(best["analysis_end"]) == 2000
    assert not candidates.empty


def test_membership_seed_and_combination_generation_are_reproducible() -> None:
    module = load_script("generate_itrdb_empirical_memberships_v01.py")
    units = tuple(f"T{i:02d}" for i in range(1, 11))
    first = module.sample_unique_combinations(
        units, 5, 20, design_key="test|subset"
    )
    second = module.sample_unique_combinations(
        units, 5, 20, design_key="test|subset"
    )
    assert first == second
    assert len(first) == 20
    assert len(set(first)) == 20
    assert all(len(value) == 5 for value in first)


def test_empirical_distribution_summary_handles_finite_values() -> None:
    module = load_script("run_itrdb_empirical_stage_b_v01.py")
    result = module.distribution_summary([1.0, 2.0, np.nan, 3.0], prefix="x_")
    assert result["x_n"] == 3
    assert result["x_median"] == pytest.approx(2.0)
    assert result["x_min"] == pytest.approx(1.0)
    assert result["x_max"] == pytest.approx(3.0)


def test_tie_summary_retains_all_equivalent_units() -> None:
    module = load_script("run_itrdb_empirical_stage_b_v01.py")
    frame = pd.DataFrame(
        {
            "omitted_unit_id": ["T03", "T01", "T02"],
            "metric": [1.0, 1.0 + 5e-13, 0.5],
        }
    )
    singular, plural, count, value = module.tied_ids(frame, "metric")
    assert singular == "T01"
    assert plural == "T01|T03"
    assert count == 2
    assert value == pytest.approx(1.0 + 5e-13)



def test_empirical_r_reader_uses_explicit_tucson_contract() -> None:
    script = (
        REPO_ROOT / "scripts" / "prepare_itrdb_empirical_stage_b_v01.R"
    ).read_text(encoding="utf-8")
    assert 'format = "auto"' not in script
    assert "format = rwl_format" in script
    assert "header = tucson_header" in script
    assert "long = tucson_long" in script
    assert "edge.zeros = preserve_zero" in script
