from __future__ import annotations

from hashlib import sha256
import importlib.util
import json
import math
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

from chronodiag import effective_dimension, eps_value

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_PATH = REPO_ROOT / "scripts" / "run_same_eps_stage_a_v01.py"
SPEC = importlib.util.spec_from_file_location("same_eps_stage_a_v01", SCRIPT_PATH)
assert SPEC is not None and SPEC.loader is not None
STAGE_A = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = STAGE_A
SPEC.loader.exec_module(STAGE_A)

FIXTURES = REPO_ROOT / "shared_fixtures"


def test_four_scenarios_share_n_rbar_and_eps_but_not_dimension() -> None:
    matrices = STAGE_A.load_scenario_matrices(FIXTURES)
    summaries = STAGE_A.scenario_summary(matrices, n_years=100)
    assert set(summaries["scenario"]) == set(STAGE_A.SCENARIOS)
    assert summaries["n_trees"].nunique() == 1
    assert summaries["n_trees"].iloc[0] == 18
    assert np.ptp(summaries["mean_interseries_correlation"]) < 1e-12
    assert np.ptp(summaries["eps"]) < 1e-12
    assert summaries["eps"].iloc[0] == pytest.approx(0.85, abs=1e-12)
    assert summaries["effective_dimension"].nunique() == 4
    assert summaries.set_index("scenario").loc[
        "two_clusters", "effective_dimension"
    ] == pytest.approx(5.99889812828725, abs=1e-10)


    indexed = summaries.set_index("scenario")
    assert indexed.loc["exchangeable", "latent_mean_reliability"] == pytest.approx(
        0.85, abs=1e-12
    )
    assert indexed.loc["unequal_loadings", "latent_mean_reliability"] == pytest.approx(
        0.8524702924390316, abs=1e-12
    )
    assert indexed.loc["influential_trees", "latent_mean_reliability"] == pytest.approx(
        0.858766870400748, abs=1e-12
    )
    assert math.isnan(indexed.loc["two_clusters", "latent_mean_reliability"])


@pytest.mark.parametrize("scenario", STAGE_A.SCENARIOS)
def test_exact_sample_reproduces_each_target_correlation(scenario: str) -> None:
    matrices = STAGE_A.load_scenario_matrices(FIXTURES)
    target = matrices[scenario]
    values = STAGE_A.exact_sample_from_correlation(target, n_years=100)
    observed = np.corrcoef(values, rowvar=False)
    assert observed == pytest.approx(target, abs=5e-12)


def test_exchangeable_subset_correlation_equals_sqrt_sss() -> None:
    matrix = STAGE_A.load_scenario_matrices(FIXTURES)["exchangeable"]
    n = matrix.shape[0]
    rbar = STAGE_A._mean_off_diagonal(matrix)
    full = np.full(n, 1.0 / n)
    for subset_size, units in [
        (1, [0]),
        (3, [0, 4, 17]),
        (9, list(range(0, 18, 2))),
        (17, list(range(17))),
    ]:
        observed = STAGE_A._weighted_correlation(
            matrix,
            STAGE_A._weight_vector(units, n),
            full,
        )
        expected = math.sqrt(
            eps_value(rbar, subset_size) / eps_value(rbar, n)
        )
        assert observed == pytest.approx(expected, abs=1e-12)


def test_cluster_aligned_split_is_nearly_unrelated_despite_eps_point_85() -> None:
    matrix = STAGE_A.load_scenario_matrices(FIXTURES)["two_clusters"]
    n = matrix.shape[0]
    observed = STAGE_A._weighted_correlation(
        matrix,
        STAGE_A._weight_vector(range(9), n),
        STAGE_A._weight_vector(range(9, 18), n),
    )
    assert eps_value(STAGE_A._mean_off_diagonal(matrix), n) == pytest.approx(
        0.85, abs=1e-12
    )
    assert observed == pytest.approx(0.01408450704225335, abs=1e-12)


def test_balanced_split_summary_has_expected_exact_counts_and_tails() -> None:
    matrices = STAGE_A.load_scenario_matrices(FIXTURES)
    summary = STAGE_A.exhaustive_balanced_split_summary(matrices).set_index(
        "scenario"
    )
    assert (summary["correlation_n"] == math.comb(17, 8)).all()
    assert summary.loc["exchangeable", "correlation_sd"] < 1e-12
    assert summary.loc["two_clusters", "correlation_min"] == pytest.approx(
        0.01408450704225335, abs=1e-12
    )
    assert summary.loc["two_clusters", "correlation_q05"] < 0.5
    assert summary.loc["unequal_loadings", "correlation_q05"] > 0.73


def test_bootstrap_membership_generation_is_frozen_and_auditable() -> None:
    counts, memberships = STAGE_A.generate_bootstrap_multiplicities(
        18, n_replicates=1000, seed=20260730
    )
    digest = sha256(counts.astype(np.int16).tobytes()).hexdigest()
    assert digest == "9dc52cb467b3d703f7c6350ebc9403fe2e8f7333a176c466046e68b8cab2ab35"
    assert counts.shape == (1000, 18)
    assert np.all(counts.sum(axis=1) == 18)
    assert memberships.groupby("replicate_id")["multiplicity"].sum().eq(18).all()
    first = memberships.loc[memberships["replicate_id"] == "B0001"]
    assert dict(zip(first["unit_id"], first["multiplicity"], strict=True)) == {
        "T02": 2,
        "T04": 2,
        "T07": 3,
        "T09": 1,
        "T10": 1,
        "T13": 3,
        "T14": 1,
        "T17": 3,
        "T18": 2,
    }


def test_loto_influence_is_metric_specific() -> None:
    matrices = STAGE_A.load_scenario_matrices(FIXTURES)
    raw, summary = STAGE_A.loto_results(matrices, n_years=100)
    influential = raw.loc[raw["scenario"] == "influential_trees"].set_index(
        "omitted_unit_id"
    )
    # High-loading units have the largest EPS and effective-dimension changes.
    assert influential.loc["T01", "delta_eps"] < influential.loc[
        "T18", "delta_eps"
    ]
    assert influential.loc["T01", "delta_effective_dimension"] > 0.9
    # Lower-loading, residual-rich units move the chronology trajectory more.
    assert influential.loc["T18", "standardized_rms_difference"] > influential.loc[
        "T01", "standardized_rms_difference"
    ]
    indexed_summary = summary.set_index("scenario")
    row = indexed_summary.loc["influential_trees"]
    assert row["correlation_min"] > 0.995
    assert row["delta_eps_min"] < -0.027

    # Exact symmetries must be represented as ties, not as arbitrary
    # platform-dependent idxmax selections.
    all_units = "|".join(f"T{index:02d}" for index in range(1, 19))
    exchangeable = indexed_summary.loc["exchangeable"]
    assert exchangeable["max_rms_unit_id"] == "T01"
    assert exchangeable["max_rms_unit_ids"] == all_units
    assert exchangeable["max_rms_tie_count"] == 18

    two_clusters = indexed_summary.loc["two_clusters"]
    assert two_clusters["max_rms_unit_id"] == "T01"
    assert two_clusters["max_rms_unit_ids"] == all_units
    assert two_clusters["max_rms_tie_count"] == 18

    assert row["max_rms_unit_id"] == "T04"
    assert row["max_rms_unit_ids"] == "|".join(
        f"T{index:02d}" for index in range(4, 19)
    )
    assert row["max_rms_tie_count"] == 15
    assert row["max_abs_delta_eps_unit_ids"] == "T01|T02|T03"
    assert row["max_abs_delta_eps_tie_count"] == 3
    assert (
        row["max_abs_delta_effective_dimension_unit_ids"]
        == "T01|T02|T03"
    )
    assert row["max_abs_delta_effective_dimension_tie_count"] == 3


def test_frozen_stage_a_outputs_are_hash_verified_and_readable() -> None:
    output_dir = REPO_ROOT / "analysis" / "same_eps_stage_a_v01"
    manifest_path = output_dir / "analysis_manifest_v01.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    for record in manifest["input_matrix_files"].values():
        path = REPO_ROOT / record["path"]
        assert sha256(path.read_bytes()).hexdigest() == record["sha256_bytes"]

    for record in manifest["output_files"].values():
        path = REPO_ROOT / record["path"]
        assert path.exists()
        assert sha256(path.read_bytes()).hexdigest() == record["sha256"]

    key = pd.read_csv(output_dir / "key_results_v01.csv").set_index("scenario")
    assert key.loc["exchangeable", "eps"] == pytest.approx(0.85, abs=1e-12)
    assert key.loc["two_clusters", "split_half_min"] == pytest.approx(
        0.0140845070422533, abs=1e-12
    )
    assert key.loc["influential_trees", "latent_mean_reliability"] == pytest.approx(
        0.858766870400748, abs=1e-12
    )

    memberships = pd.read_csv(
        output_dir / "bootstrap_multiplicities_v01.csv.gz"
    )
    assert memberships.groupby("replicate_id")["multiplicity"].sum().eq(18).all()
    exact = pd.read_csv(output_dir / "exact_same_eps_long_v01.csv.gz")
    assert exact.shape == (4 * 100 * 18, 5)
