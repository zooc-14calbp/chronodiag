from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass
from pathlib import Path
import re
import warnings

import numpy as np
import pandas as pd

KeyParser = Callable[[str], Mapping[str, object] | None]
_STOP_MARKERS = {999: 0.01, -9999: 0.001}


@dataclass(frozen=True)
class _TucsonRow:
    series_id: str
    start_year: int
    values: tuple[int | None, ...]
    line_number: int


def _apply_key_parser(data: pd.DataFrame, key_parser: KeyParser | None) -> pd.DataFrame:
    result = data.copy()
    result["series_id"] = result["series_id"].astype(str)
    if key_parser is None:
        if "tree_id" not in result.columns:
            result["tree_id"] = result["series_id"]
        else:
            result["tree_id"] = result["tree_id"].where(
                result["tree_id"].notna(), result["series_id"]
            )
        return result

    rows: list[dict[str, object]] = []
    for series_id in result["series_id"].drop_duplicates():
        parsed = key_parser(series_id)
        row: dict[str, object] = {"series_id": series_id}
        if parsed is not None:
            row.update(dict(parsed))
        if str(row.get("series_id", series_id)) != series_id:
            raise ValueError("key_parser must not change series_id")
        rows.append(row)
    metadata = pd.DataFrame(rows)
    if metadata["series_id"].duplicated().any():
        raise ValueError("key_parser returned duplicate series_id values")
    if "tree_id" not in metadata.columns:
        metadata["tree_id"] = metadata["series_id"]
    else:
        metadata["tree_id"] = metadata["tree_id"].where(
            metadata["tree_id"].notna(), metadata["series_id"]
        )

    row_id = "__chronodiag_row_id"
    result[row_id] = np.arange(len(result), dtype=np.int64)
    original_tree_by_row: pd.Series | None = None
    if "tree_id" in result.columns:
        original_tree_by_row = pd.Series(
            result["tree_id"].to_numpy(copy=True),
            index=result[row_id].to_numpy(copy=True),
        )
        result = result.drop(columns=["tree_id"])

    result = result.merge(
        metadata,
        on="series_id",
        how="left",
        validate="many_to_one",
        sort=False,
    )
    result = result.sort_values(row_id, kind="stable").reset_index(drop=True)
    if original_tree_by_row is not None:
        original_tree = original_tree_by_row.reindex(result[row_id]).reset_index(drop=True)
        result["tree_id"] = original_tree.where(
            original_tree.notna(), result["tree_id"]
        )
    return result.drop(columns=[row_id])


def _validate_years(years: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(years, errors="raise")
    values = numeric.to_numpy(dtype=float)
    if not np.isfinite(values).all() or not np.equal(values, np.floor(values)).all():
        raise ValueError("year must contain finite integers")
    return numeric.astype(int)


def _validate_values(values: pd.Series, name: str = "value") -> pd.Series:
    converted = pd.to_numeric(values, errors="coerce")
    bad = values.notna() & converted.isna()
    if bad.any():
        examples = values.loc[bad].astype(str).head(3).tolist()
        raise ValueError(f"{name} contains non-numeric values: {examples}")
    numeric = converted.to_numpy(dtype=float)
    finite_or_missing = np.isfinite(numeric) | np.isnan(numeric)
    if not finite_or_missing.all():
        raise ValueError(f"{name} must contain finite numbers or missing values")
    return converted


def _validate_unique_observations(data: pd.DataFrame) -> None:
    duplicated = data.duplicated(["series_id", "year"], keep=False)
    if duplicated.any():
        example = data.loc[duplicated, ["series_id", "year"]].iloc[0]
        raise ValueError(
            f"duplicate observation for {example['series_id']}, {int(example['year'])}"
        )


def wide_to_long(
    data: pd.DataFrame,
    *,
    key_parser: KeyParser | None = None,
    drop_missing: bool = True,
) -> pd.DataFrame:
    """Convert a ``year, key1, key2, ...`` table to canonical long form."""
    if not isinstance(data, pd.DataFrame) or data.shape[1] < 2:
        raise ValueError("wide input must contain year and at least one series column")
    first = str(data.columns[0])
    if first.casefold() != "year":
        raise ValueError("the first wide-CSV column must be named 'year'")
    if pd.Index(data.columns).duplicated().any():
        raise ValueError("wide input contains duplicate column names")

    wide = data.copy()
    wide[first] = _validate_years(wide[first])
    if wide[first].duplicated().any():
        raise ValueError("wide input contains duplicate years")
    long = wide.melt(id_vars=[first], var_name="series_id", value_name="value")
    long = long.rename(columns={first: "year"})
    long["value"] = _validate_values(long["value"])
    if drop_missing:
        long = long.dropna(subset=["value"])
    long = _apply_key_parser(long, key_parser)
    _validate_unique_observations(long)
    columns = ["year", "series_id", "tree_id", "value"]
    columns.extend(c for c in long.columns if c not in columns)
    return long.loc[:, columns].sort_values(["series_id", "year"], kind="stable").reset_index(drop=True)


def read_wide_csv(
    path: str | Path,
    *,
    key_parser: KeyParser | None = None,
    drop_missing: bool = True,
    **read_csv_kwargs: object,
) -> pd.DataFrame:
    data = pd.read_csv(path, **read_csv_kwargs)
    return wide_to_long(data, key_parser=key_parser, drop_missing=drop_missing)


def read_long_csv(
    path: str | Path,
    *,
    key_parser: KeyParser | None = None,
    drop_missing: bool = True,
    **read_csv_kwargs: object,
) -> pd.DataFrame:
    data = pd.read_csv(path, **read_csv_kwargs)
    required = {"year", "series_id", "value"}
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"long CSV is missing required columns: {sorted(missing)}")
    result = data.copy()
    result["series_id"] = result["series_id"].astype(str)
    result["year"] = _validate_years(result["year"])
    result["value"] = _validate_values(result["value"])
    if drop_missing:
        result = result.dropna(subset=["value"])
    result = _apply_key_parser(result, key_parser)
    _validate_unique_observations(result)
    columns = ["year", "series_id", "tree_id", "value"]
    columns.extend(c for c in result.columns if c not in columns)
    return result.loc[:, columns].sort_values(["series_id", "year"], kind="stable").reset_index(drop=True)


def _parse_fixed_width_tucson_line(line: str, line_number: int) -> _TucsonRow | None:
    for id_width, year_width in ((8, 4), (7, 5)):
        prefix_width = id_width + year_width
        if len(line) < prefix_width + 1:
            continue
        series_id = line[:id_width].strip()
        year_text = line[id_width:prefix_width].strip()
        if not series_id or not re.fullmatch(r"[+-]?\d+", year_text):
            continue
        start_year = int(year_text)
        if start_year < -10000 or start_year > 10000:
            continue
        remainder = line[prefix_width:]
        chunks = [remainder[i:i + 6] for i in range(0, len(remainder), 6)]
        nonblank = [i for i, chunk in enumerate(chunks) if chunk.strip()]
        if not nonblank:
            continue
        last = nonblank[-1]
        parsed: list[int | None] = []
        failed = False
        for chunk in chunks[:last + 1]:
            text = chunk.strip()
            if text in {"", "."}:
                parsed.append(None)
                continue
            if not re.fullmatch(r"[+-]?\d+", text):
                failed = True
                break
            parsed.append(int(text))
        if not failed and parsed:
            return _TucsonRow(series_id, start_year, tuple(parsed), line_number)
    return None


def _parse_whitespace_tucson_line(line: str, line_number: int) -> _TucsonRow | None:
    parts = line.split()
    if len(parts) < 3:
        return None
    year_index = None
    # The ID can contain one embedded whitespace component (e.g. "Q 3608").
    for candidate in (1, 2):
        if candidate < len(parts) and re.fullmatch(r"[+-]?\d+", parts[candidate]):
            year = int(parts[candidate])
            if -10000 <= year <= 10000:
                year_index = candidate
                break
    if year_index is None:
        return None
    series_id = " ".join(parts[:year_index])
    values: list[int | None] = []
    for text in parts[year_index + 1:]:
        if text == ".":
            values.append(None)
        elif re.fullmatch(r"[+-]?\d+", text):
            values.append(int(text))
        else:
            return None
    if not values:
        return None
    return _TucsonRow(series_id, int(parts[year_index]), tuple(values), line_number)


def _parse_tucson_data_line(line: str, line_number: int) -> _TucsonRow | None:
    return _parse_fixed_width_tucson_line(line, line_number) or _parse_whitespace_tucson_line(line, line_number)


def _read_tucson_rows(path: str | Path, *, strict: bool = True) -> tuple[list[_TucsonRow], list[str]]:
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()
    rows: list[_TucsonRow] = []
    header: list[str] = []
    data_started = False
    for line_number, line in enumerate(lines, start=1):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        parsed = _parse_tucson_data_line(line, line_number)
        if parsed is None:
            if not data_started:
                header.append(line)
                continue
            message = f"unparseable Tucson line {line_number}: {line[:80]!r}"
            if strict:
                raise ValueError(message)
            warnings.warn(message, RuntimeWarning, stacklevel=2)
            continue
        data_started = True
        rows.append(parsed)
    if not rows:
        raise ValueError("no Tucson RWL data rows were detected")
    if len(header) not in {0, 3}:
        warnings.warn(
            f"detected {len(header)} non-data header lines; standard Tucson headers use three",
            RuntimeWarning,
            stacklevel=2,
        )
    return rows, header


def _infer_precision_by_series(
    rows: list[_TucsonRow], precision: str | float
) -> tuple[dict[str, float], set[tuple[int, int]]]:
    if precision != "auto":
        scale = float(precision)
        if scale not in {0.01, 0.001}:
            raise ValueError("precision must be 'auto', 0.01, or 0.001")
        return {row.series_id: scale for row in rows}, set()

    rows_by_series: dict[str, list[tuple[int, _TucsonRow]]] = {}
    for index, row in enumerate(rows):
        rows_by_series.setdefault(row.series_id, []).append((index, row))

    scale_by_series: dict[str, float] = {}
    marker_positions: set[tuple[int, int]] = set()
    unknown: list[str] = []
    for series_id, indexed_rows in rows_by_series.items():
        # Latest row by calendar start; file order breaks ties.
        row_index, final_row = max(indexed_rows, key=lambda item: (item[1].start_year, item[0]))
        nonmissing = [i for i, value in enumerate(final_row.values) if value is not None]
        marker = None
        if nonmissing:
            final_position = nonmissing[-1]
            final_value = final_row.values[final_position]
            if final_value in _STOP_MARKERS:
                marker = int(final_value)
                marker_positions.add((row_index, final_position))
        if marker is None:
            unknown.append(series_id)
        else:
            scale_by_series[series_id] = _STOP_MARKERS[marker]

    if unknown:
        preview = ", ".join(unknown[:5]) + (" ..." if len(unknown) > 5 else "")
        raise ValueError(
            "Tucson precision could not be inferred for series without a final "
            f"999/-9999 stop marker: {preview}; pass precision explicitly"
        )
    return scale_by_series, marker_positions


def read_tucson_rwl(
    path: str | Path,
    *,
    precision: str | float = "auto",
    key_parser: KeyParser | None = None,
    drop_missing: bool = True,
    preserve_zero: bool = True,
    strict: bool = True,
) -> pd.DataFrame:
    """Read a Tucson/ITRDB decadal RWL file into canonical long form.

    Ring widths are returned in millimetres.  Zero-width rings are preserved by
    default, matching ``dplR::read.rwl(..., edge.zeros = TRUE)``.  Negative
    non-marker values are treated as missing.
    """
    rows, header = _read_tucson_rows(path, strict=strict)
    scale_by_series, marker_positions = _infer_precision_by_series(rows, precision)

    observations: list[dict[str, object]] = []
    seen: set[tuple[str, int]] = set()
    for row_index, row in enumerate(rows):
        scale = scale_by_series[row.series_id]
        for offset, raw in enumerate(row.values):
            if (row_index, offset) in marker_positions:
                continue
            year = row.start_year + offset
            key = (row.series_id, year)
            if key in seen:
                raise ValueError(f"duplicate Tucson observation for {row.series_id}, {year}")
            seen.add(key)
            if raw is None or raw < 0 or (raw == 0 and not preserve_zero):
                value = np.nan
            else:
                value = raw * scale
            observations.append(
                {
                    "year": year,
                    "series_id": row.series_id,
                    "value": value,
                }
            )

    result = pd.DataFrame(observations)
    if drop_missing:
        result = result.dropna(subset=["value"])
    result = _apply_key_parser(result, key_parser)
    _validate_unique_observations(result)
    result["unit"] = "mm"
    result["source_format"] = "tucson_rwl"
    result.attrs["tucson_header"] = tuple(header)
    result.attrs["precision_by_series"] = dict(scale_by_series)
    columns = ["year", "series_id", "tree_id", "value", "unit", "source_format"]
    columns.extend(c for c in result.columns if c not in columns)
    return result.loc[:, columns].sort_values(["series_id", "year"], kind="stable").reset_index(drop=True)


def read_chronodiag(
    path: str | Path,
    *,
    format: str = "auto",
    key_parser: KeyParser | None = None,
    drop_missing: bool = True,
    **kwargs: object,
) -> pd.DataFrame:
    path_obj = Path(path)
    selected = format.casefold()
    if selected == "auto":
        suffix = path_obj.suffix.casefold()
        if suffix in {".rwl", ".tuc"}:
            selected = "rwl"
        elif suffix == ".csv":
            header = pd.read_csv(path_obj, nrows=0)
            folded = {str(column).casefold() for column in header.columns}
            selected = "csv_long" if {"series_id", "value"}.issubset(folded) else "csv_wide"
        else:
            raise ValueError("cannot auto-detect format; use csv_wide, csv_long, or rwl")

    if selected in {"csv", "csv_wide", "wide"}:
        return read_wide_csv(path_obj, key_parser=key_parser, drop_missing=drop_missing, **kwargs)
    if selected in {"csv_long", "long"}:
        return read_long_csv(path_obj, key_parser=key_parser, drop_missing=drop_missing, **kwargs)
    if selected in {"rwl", "tucson", "tucson_rwl"}:
        return read_tucson_rwl(path_obj, key_parser=key_parser, drop_missing=drop_missing, **kwargs)
    raise ValueError(f"unsupported format: {format}")
