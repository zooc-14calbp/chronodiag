# R4 public snapshot checks

Preparation date: 2026-11-02. The following checks were run locally on the
candidate public tree based on commit
`da165d0624937e71c4db016cb14cbd1a6c036de2`.

Environment: Windows, Python 3.12.14, NumPy 2.5.4, pandas 3.1.0, pytest 9.2.0,
Matplotlib 3.11.2, setuptools 84.2.1, and R 4.6.1. These are the observed check
versions, not additional minimum-version requirements.

| Check | Result |
| --- | --- |
| `python -m pytest python/tests -q` | 61 passed |
| `testthat::test_local("r/chronodiag", stop_on_failure=TRUE)` | 63 expectations passed; no test failures, warnings or skips |
| All five scripts in `r/parity_tests/` | Passed shared metrics, I/O, composition, influence and same-EPS fixtures |
| `R CMD build` and `R CMD check --no-manual` on the resulting source package | Status: OK |
| Python wheel build | Passed; wheel contains full GPL text and `License-Expression: GPL-3.0-or-later` |
| `scripts/validate_same_eps_stage_a_v01.py` | End-to-end controlled-analysis reproduction passed |
| `scripts/validate_itrdb_empirical_stage_b_v01.py` | Frozen empirical validation passed |
| `scripts/validate_yoshino_empirical_stage_c_v01.py` | Frozen empirical validation passed |
| `revision_support/reproduce_revision.py` | All three examples and 2,000 frozen bootstrap draws per example matched the numerical references; point-only figure generated |
| Existing accessible JSON input path/hash references | All 62 matched the distributed files |
| Two revised PowerShell wrappers | Syntax parsing passed |

The Stage A validator regenerates the controlled analysis. Stage B/C validators
check archived inputs/cohorts and regenerate memberships and diagnostic tables
from archived preprocessed series; R/dplR raw detrending was not rerun for these
checks. The R4 script uses the existing memberships and generates no new
resampling or surrogate series. Its numerical comparison uses absolute tolerance
`1e-11` and relative tolerance `1e-10`. Figure PDF bytes can vary with rendering
metadata, fonts and dependency versions.

Core calculation files, raw values, preprocessed values, membership streams and
numerical outputs retain their original Git content. Three CSV source-path fields
were made portable, and five existing input-hash references were corrected.
`public_release_changes.json` records those provenance-only changes. Wrapper
defaults, package metadata, documentation, license files and the added revision
support are release-preparation changes.

GitHub CI is configured to repeat the existing tests and frozen validators, with
the added R4 numerical check on Python 3.12. This record reports the local checks;
the GitHub run attached to the published commit is the separate CI evidence.
