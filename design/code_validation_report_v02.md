# chronodiag pre-GitHub validation report v02

Date: 2026-07-29

## Scope

This validation round tested the Python package as source code and as a built
wheel, repeated the six-file ITRDB integration audit, expanded shared numerical
fixtures, and revised the R package metadata and GitHub Actions workflow in
preparation for a private GitHub repository.

## Python unit and parity tests

Command:

```text
cd python
python -m pytest -q
```

Result:

```text
16 passed
```

Coverage of the deterministic test suite includes:

- EPS reference algebra;
- direct SSS evaluation at `rbar = 0`;
- admissible equicorrelation boundaries;
- effective-dimension validation;
- wide and long CSV adapters;
- Tucson header, zero-ring, missing-field, embedded-space-ID, and dual-precision
  handling;
- duplicate observation rejection;
- all four same-EPS/different-structure correlation fixtures.

## Python packaging test

The package metadata was corrected to declare both runtime dependencies used by
the implementation:

```text
numpy >= 1.24
pandas >= 2.0
```

Explicit src-layout package discovery was also added. A wheel was built with
setuptools:

```text
chronodiag-0.0.1.dev0-py3-none-any.whl
```

Wheel SHA-256 in the validation environment:

```text
0abc1edf6f8d5d1ef0736e48f5564373715aaba49de99c1df0ac3cd217ec1589
```

The wheel was installed into a clean target directory and imported outside the
source tree. It returned:

```text
EPS(18, 0.23943661971830986) = 0.8500000000000001
SSS(5 | 20, 0) = 0.25
```

## Full ITRDB integration audit

The audit was rerun through the installed wheel on the six hash-verified raw RWL
files. The recomputed table exactly matched the frozen v01 audit for the shared
fields.

| Dataset | Series | Period | Observations | Zero rings | Full common period |
|---|---:|---:|---:|---:|---:|
| ca533 | 34 | 626–1983 | 23,276 | 234 | none across all 34 |
| germ011 | 66 | 1338–1967 | 8,520 | 0 | none across all 66 |
| germ012 | 44 | 1800–1969 | 5,290 | 0 | 1873–1913 (41 yr) |
| germ6 | 138 | 1376–1972 | 14,104 | 0 | none across all 138 |
| fran003 | 25 | 1531–1979 | 7,994 | 0 | 1819–1976 (158 yr) |
| brit042 | 28 | 1649–1992 | 5,720 | 0 | none across all 28 |

## R package corrections prepared

The R package metadata now:

- lists `dplR`, `roxygen2`, and `testthat` in `Suggests`;
- exports all algebra and CSV/RWL input functions;
- retains native R implementations rather than wrapping Python;
- has repository-level parity scripts for shared metric and RWL fixtures.

R was unavailable in the present execution environment, so the R code has not yet
been executed here. The revised GitHub Actions workflow will:

1. install R dependencies;
2. generate documentation with roxygen2;
3. run shared metric parity fixtures;
4. run the shared Tucson RWL fixture; and
5. run `R CMD check`, failing on warnings.

## GitHub Actions matrix

The revised workflow tests Python 3.10, 3.12, and 3.13 on Ubuntu and performs a
separate R package/parity/check job. The first private-repository CI run is the
remaining gate before composition diagnostics are implemented.

## Decision

The Python baseline is ready for the first private GitHub commit. The R baseline
is ready for CI validation but should not be described as fully validated until the
first successful GitHub Actions run or an equivalent local R check.
