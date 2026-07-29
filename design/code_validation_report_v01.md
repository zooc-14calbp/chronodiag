# chronodiag code validation report v01

Date: 2026-07-29

## 1. Scope

The current Dropbox scaffold was inspected before GitHub publication.  The
initial validation focused on algebraic functions and R/Python input adapters,
using six real ITRDB Tucson RWL files:

- ca533
- germ011
- germ012
- germ6
- fran003
- brit042

## 2. Important issues found in the current Python baseline

### 2.1 Three-line Tucson headers were not safely skipped

The parser's whitespace fallback attempted to convert textual header fields to
integers.  As a result, `germ011`, `germ012`, `germ6`, and `brit042` failed before
data were returned.  `fran003` happened to pass because its fixed-width header did
not trigger the same failure path.

### 2.2 Zero-width rings were converted to missing values

The baseline used `raw <= 0` as the missing-value condition.  In ca533 this
converted 234 locally absent/zero-width rings to `NaN`.  The corrected behavior
preserves zero by default, matching the default `edge.zeros = TRUE` behavior of
`dplR::read.rwl()`.

### 2.3 SSS failed at zero mean correlation

The baseline computed SSS as `EPS(n_sub) / EPS(n_full)`, which becomes `0/0` at
`rbar = 0`.  The direct SSS formula is well defined and returns `n_sub/n_full`.

### 2.4 Blank fixed-width fields could shift calendar years

The baseline discarded every blank six-character field.  An internal blank must
retain its calendar position as a missing observation; only trailing blanks may be
ignored.

### 2.5 Precision inference was too permissive

Any occurrence of `999` or `-9999` was interpreted as a precision marker.  A
legitimate value can equal 999 in the integer representation.  The corrected
prototype recognizes a marker only as the final nonblank value in the final row of
a series, or requires explicit precision.

### 2.6 Long-CSV key parsing could misalign existing tree IDs

The Python baseline reassigned existing `tree_id` values after a merge using
integer-index alignment.  Merge reordering could attach a tree ID to the wrong
observation.  The corrected prototype preserves stable row identity.

## 3. Corrected local Python prototype

The local prototype now:

- detects and retains zero or three header lines;
- supports standard 8+4 and long-year 7+5 Tucson layouts;
- preserves series IDs containing spaces;
- handles CR, LF, CRLF, and CRCRLF line endings;
- preserves zero-width rings by default;
- distinguishes missing values from stop markers;
- requires auditable precision inference;
- rejects duplicate `series_id`–`year` observations;
- computes SSS directly at `rbar = 0`;
- validates the admissible equicorrelation range.

## 4. Test result

Repository-ready Python unit-test result:

```text
11 passed
```

These tests cover algebraic edge cases, zero-ring preservation, missing fixed-width
fields, embedded-space series IDs, three-line headers, both Tucson precisions, CSV
validation, and effective-dimension validation.  A separate integration audit then
read all six full ITRDB RWL files successfully and reproduced the counts reported
below.

R is not installed in the current execution environment.  R code can be revised
against the same specification, but execution must be performed on the user's
machine or in GitHub Actions before release.

## 5. Raw RWL audit summary

See `itrdb_raw_audit_v01.csv` for the machine-readable table.  Key findings are:

- ca533: 34 series, 626–1983, 23,276 observations, 234 zero rings;
- germ012: 44 core/radius series, with a key pattern consistent with 22 paired
  biological units; all 44 overlap for 41 years;
- fran003: all 25 series share 1819–1976 (158 years), making it the cleanest
  complete-case structural example;
- germ011 and germ6 are strongly changing historical assemblages rather than
  simple dense modern samples;
- brit042 provides a useful long record and a practical line-ending/header stress
  test.

## 6. Recommended next sequence

1. Replace the Python algebra and I/O files with the corrected versions.
2. Apply the same algebraic contract to R (`sss_value()` in particular).
3. Add the six-file local audit script and smaller distributable fixtures.
4. Run `R CMD check`, `pytest`, and cross-language parity tests in GitHub Actions.
5. Only after CI passes, create the private GitHub repository and upload the
   validated scaffold.
6. Update the manuscript and Supplement from the frozen mathematical
   specification rather than from implementation-specific prose.
