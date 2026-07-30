from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]


def load_script(name: str):
    path = REPO_ROOT / "scripts" / name
    spec = importlib.util.spec_from_file_location(path.stem, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_yoshino_contract_and_keymap() -> None:
    contract = pd.read_csv(
        REPO_ROOT / "data" / "yoshino" / "yoshino_dataset_contract_v01.csv"
    ).iloc[0]
    keymap = pd.read_csv(
        REPO_ROOT / "data" / "yoshino" / "yoshino_keymap_v01.csv",
        dtype={
            "series_id": str,
            "tree_id": str,
            "specimen_id": str,
            "source_group": str,
            "site_id": str,
        },
    )
    assert contract["dataset"] == "yoshino"
    assert contract["input_level"] == "specimen_mean"
    assert int(contract["n_raw_series_expected"]) == 60
    assert int(contract["n_trees_expected"]) == 50
    assert int(contract["analysis_start"]) == 1850
    assert int(contract["analysis_end"]) == 2004
    assert len(keymap) == 60
    assert keymap["series_id"].nunique() == 60
    assert keymap["tree_id"].nunique() == 60
    assert (keymap["series_id"] == keymap["tree_id"]).all()
    assert set(keymap["source_group"]) == {"001", "002", "003"}
    assert set(keymap["input_level"]) == {"specimen_mean"}
    selected = keymap[
        keymap["include_primary_cohort"].astype(str).str.casefold().eq("true")
    ]
    assert selected["tree_id"].nunique() == 50


def test_yoshino_value_blind_cohort_is_reproduced() -> None:
    module = load_script("audit_yoshino_empirical_stage_c_v01.py")
    keymap = pd.read_csv(
        REPO_ROOT / "data" / "yoshino" / "yoshino_keymap_v01.csv",
        dtype={"series_id": str, "tree_id": str, "source_group": str},
    )
    best, candidates = module.select_cohort(
        keymap, minimum_years=50, minimum_years_per_tree_ratio=3.0
    )
    assert int(best["n_trees"]) == 50
    assert int(best["analysis_start"]) == 1850
    assert int(best["analysis_end"]) == 2004
    assert int(best["analysis_years"]) == 155
    assert not candidates.empty


def test_yoshino_r_preprocessing_uses_specimen_mean_tucson_input() -> None:
    script = (
        REPO_ROOT / "scripts" / "prepare_yoshino_empirical_stage_c_v01.R"
    ).read_text(encoding="utf-8")
    assert 'format = "tucson"' in script
    assert "header = FALSE" in script
    assert "edge.zeros = TRUE" in script
    assert "specimen_mean" in script
    assert "rowMeans(" not in script


def test_yoshino_manifest_paths_accept_external_temporary_directories(
    tmp_path: Path,
) -> None:
    external = tmp_path / "regenerated" / "output.csv.gz"
    external.parent.mkdir(parents=True)
    external.write_bytes(b"fixture")

    generator = load_script("generate_yoshino_empirical_memberships_v01.py")
    runner = load_script("run_yoshino_empirical_stage_c_v01.py")

    expected = external.resolve().as_posix()
    assert generator.stable_manifest_path(external, REPO_ROOT) == expected
    assert runner.stable_manifest_path(external, REPO_ROOT) == expected

    inside = REPO_ROOT / "analysis" / "fixture.csv"
    assert generator.stable_manifest_path(inside, REPO_ROOT) == "analysis/fixture.csv"
    assert runner.stable_manifest_path(inside, REPO_ROOT) == "analysis/fixture.csv"


def test_yoshino_raw_manifest_is_line_ending_invariant(tmp_path: Path) -> None:
    validator = load_script("validate_yoshino_empirical_stage_c_v01.py")
    manifest = pd.read_csv(
        REPO_ROOT / "data" / "yoshino" / "yoshino_raw_manifest_v01.csv",
        encoding="utf-8-sig",
    ).iloc[0]
    raw = REPO_ROOT / "data" / "yoshino" / "raw" / "nr_all.rwl"

    assert validator.normalized_lf_sha256(raw) == manifest["normalized_lf_sha256"]
    assert validator.normalized_lf_size(raw) == int(
        manifest["normalized_lf_size_bytes"]
    )

    crlf = tmp_path / "crlf.rwl"
    lf = tmp_path / "lf.rwl"
    crlf.write_bytes(b"AAA 1900 100\r\nAAA 1910 200\r\n")
    lf.write_bytes(b"AAA 1900 100\nAAA 1910 200\n")
    assert validator.normalized_lf_sha256(crlf) == validator.normalized_lf_sha256(lf)
    assert validator.normalized_lf_size(crlf) == validator.normalized_lf_size(lf)

