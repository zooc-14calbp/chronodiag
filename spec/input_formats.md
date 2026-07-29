# Input formats and adapter contract v02

## 1. Canonical representation

All readers return long-form observations with at least:

- `year` — integer calendar year;
- `series_id` — exact input key;
- `tree_id` — biological resampling unit, defaulting to `series_id`;
- `value` — finite numeric value or missing value.

Optional parsed identifiers include `site_id`, `specimen_id`, and `core_id`.
Readers must not infer independent trees solely from different series keys.
Dataset-specific parsing and aggregation rules must be explicit and auditable.

Duplicate `series_id`–`year` observations are rejected.

## 2. Wide CSV

Required layout:

```text
year,key1,key2,...
```

Rules:

- the first column is named `year` (case-insensitive);
- remaining column names are preserved exactly as series keys;
- blank fields and declared NA tokens are missing observations;
- values are not rescaled;
- duplicate years and duplicate column names are rejected;
- non-numeric nonmissing values are rejected;
- `drop_missing = TRUE` omits missing cells from the returned table.

## 3. Long CSV

Minimum columns:

```text
year,series_id,value
```

`tree_id` is optional and defaults to `series_id`.  Existing nonmissing `tree_id`
values take precedence over parser-derived fallback values.  Additional columns
are preserved.

## 4. Tucson/ITRDB RWL

The RWL adapter returns ring width in millimetres.

The R implementation delegates low-level parsing to `dplR::read.rwl()` and maps
the user-facing `preserve_zero` argument to dplR's `edge.zeros` behavior.  The
Python implementation is native but follows the same public contract.

Supported v0.1 features:

- headerless files and standard three-line headers;
- fixed-width 8-character ID + 4-character year;
- long-year 7-character ID + 5-character year;
- series IDs containing spaces;
- CR, LF, CRLF, and files producing blank lines under CRCRLF decoding;
- `999` precision marker for 0.01 mm;
- `-9999` precision marker for 0.001 mm;
- internal blank or `.` missing fields without shifting calendar years.

Zero-width rings are valid observations and are preserved by default.  Setting
`preserve_zero = FALSE` converts them to missing values explicitly.  Negative
non-marker values are treated as missing.

In `precision = "auto"` mode, Python recognizes a precision marker only as the
final nonblank value of the final row of a series.  If precision cannot be
inferred, the reader stops and requires `0.01` or `0.001`; it does not silently
assume a scale.

Other formats supported by dplR (compact, Heidelberg, TRiDaS) are outside the
v0.1 cross-language parity contract.

## 5. Format auto-detection

`read_chronodiag(path, format = "auto")` applies:

1. `.rwl` or `.tuc` -> Tucson RWL;
2. `.csv` containing `series_id` and `value` -> long CSV;
3. other `.csv` -> wide CSV;
4. otherwise, explicit format required.

## 6. Key parsing

A key parser is applied once to each unique `series_id`.  It may return
`tree_id`, `specimen_id`, `site_id`, and `core_id`, but must not alter the raw
`series_id` or return duplicate IDs.

For Yoshino, a dataset adapter will parse the frozen raw measurement keycode.
The generic package will not embed Yoshino-specific substring positions.
