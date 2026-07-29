from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from chronodiag import read_chronodiag, read_tucson_rwl, wide_to_long

FIXTURES = Path(__file__).resolve().parents[2] / "shared_fixtures"


def test_wide_csv_adapter() -> None:
    source = pd.DataFrame({"year": [1900, 1901], "KEY_A": [1.0, 2.0], "KEY_B": [3.0, None]})
    result = wide_to_long(source)
    assert list(result.columns[:4]) == ["year", "series_id", "tree_id", "value"]
    assert len(result) == 3
    assert (result["series_id"] == result["tree_id"]).all()


def test_tucson_reader_handles_header_zero_missing_spaces_and_precisions() -> None:
    result = read_tucson_rwl(FIXTURES / "example_tucson_v02.rwl", drop_missing=False)
    assert len(result.attrs["tucson_header"]) == 3
    assert set(result["series_id"]) == {"ZERO01", "MISS01", "Q  525", "MICRO1"}
    zero = result[result["series_id"] == "ZERO01"].set_index("year")["value"]
    assert zero.loc[1900] == pytest.approx(1.0)
    assert zero.loc[1901] == pytest.approx(0.0)
    assert zero.loc[1902] == pytest.approx(2.0)
    missing = result[result["series_id"] == "MISS01"].set_index("year")["value"]
    assert pd.isna(missing.loc[1901])
    micro = result[result["series_id"] == "MICRO1"].set_index("year")["value"]
    assert micro.loc[1900] == pytest.approx(1.234)
    assert set(result["unit"]) == {"mm"}


def test_zero_policy_is_explicit() -> None:
    result = read_tucson_rwl(
        FIXTURES / "example_tucson_v02.rwl",
        preserve_zero=False,
        drop_missing=False,
    )
    zero = result[(result["series_id"] == "ZERO01") & (result["year"] == 1901)]
    assert zero["value"].isna().all()


def test_key_parser_is_optional_and_auditable(tmp_path: Path) -> None:
    path = tmp_path / "wide.csv"
    pd.DataFrame({"year": [1900], "KEY_A": [1.0], "KEY_B": [2.0]}).to_csv(path, index=False)

    def parser(key: str) -> dict[str, str]:
        return {"tree_id": key.split("_")[0], "site_id": "TEST"}

    result = read_chronodiag(path, key_parser=parser)
    assert set(result["tree_id"]) == {"KEY"}
    assert set(result["site_id"]) == {"TEST"}


def test_duplicate_long_observation_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "long.csv"
    pd.DataFrame(
        {"year": [1900, 1900], "series_id": ["A", "A"], "value": [1.0, 2.0]}
    ).to_csv(path, index=False)
    with pytest.raises(ValueError, match="duplicate"):
        read_chronodiag(path)


def test_existing_tree_ids_remain_aligned_after_parser_merge(tmp_path: Path) -> None:
    path = tmp_path / "long_interleaved.csv"
    pd.DataFrame(
        {
            "year": [1900, 1900, 1901, 1901],
            "series_id": ["B1", "A1", "B1", "A1"],
            "tree_id": ["TREE_B", "TREE_A", "TREE_B", "TREE_A"],
            "value": [1.0, 2.0, 3.0, 4.0],
        }
    ).to_csv(path, index=False)

    def parser(key: str) -> dict[str, str]:
        return {"tree_id": f"PARSED_{key}", "site_id": "TEST"}

    result = read_chronodiag(path, key_parser=parser)
    observed = result.groupby("series_id")["tree_id"].unique().to_dict()
    assert observed["A1"].tolist() == ["TREE_A"]
    assert observed["B1"].tolist() == ["TREE_B"]
    assert set(result["site_id"]) == {"TEST"}
