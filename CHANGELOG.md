# Changelog

## Unreleased

- Added explicit `tree_bootstrap()` and `leave_one_tree_out()` implementations in
  native R and Python.
- Defined standard `N`-out-of-`N` tree bootstrap draws through auditable positive
  integer multiplicities, with an explicit opt-in for `m`-out-of-`N` designs.
- Added bootstrap membership metadata including unique, omitted, and duplicated
  units, maximum multiplicity, and Kish effective unit count.
- Made direct chronology changes—RMS, standardized RMS, mean absolute, maximum
  absolute, signed maximum-year difference, and variance ratio—the primary
  bootstrap/LOTO outputs; shared-data correlations remain descriptive.
- Added LOTO changes in mean correlation, EPS, and effective dimension evaluated
  on one fixed full complete period so unit influence is not confounded with an
  interval change.
- Added a language-neutral influence-diagnostics specification and result schema.
- Added deterministic bootstrap and LOTO fixtures, frozen R/Python expected
  outputs, a fixture generator, and CI parity checks.
- Added explicit `subset_stability()` and `split_half_stability()` implementations
  in native R and Python.
- Added a language-neutral composition-diagnostics specification and result schema.
- Added frozen subset and split-half partition fixtures with expected cross-language
  outputs.
- Added an exact exchangeable fixture linking squared subset-to-full correlation
  to SSS and equal-half disjoint correlation to EPS.
- Made shared-unit overlap, chronology constructor, correlation method, full-pool
  size, sample-depth rule, overlap requirement, comparison exhaustiveness, and
  realized chronology depths explicit in every composition result row.
- Required literal split halves to be balanced by default while allowing an
  explicitly requested unbalanced disjoint comparison.
- Kept random partition generation outside the v01 numerical engine so that R and
  Python can evaluate identical auditable compositions.
- Fixed the project, R package, and Python package working name as `chronodiag`.
- Added `pandas` as an explicit Python runtime dependency and explicit src-layout package discovery.
- Added four shared same-EPS/different-structure parity tests; the baseline Python suite contains 16 passing tests.
- Verified that the Python wheel builds and reproduces the six-file ITRDB audit.
- Expanded GitHub Actions to Python 3.10/3.12/3.13 and an R documentation, parity, and `R CMD check` job.
- Corrected the R package metadata and exports for CSV/RWL input functions.
- Restored physically blank six-character Tucson data fields to missing values in the R adapter while retaining genuine zero-width rings.
- Added native R and Python implementations governed by a language-neutral scientific contract.
- Added wide CSV, long CSV, and Tucson RWL input adapters.
- Corrected SSS to use its direct formula, including the `rbar = 0` case.
- Added admissible-range validation for EPS and SSS.
- Corrected Python Tucson parsing for standard three-line headers, embedded-space IDs, internal missing fields, and both precision markers.
- Preserved zero-width rings by default; ca533 retains 234 zero rings.
- Added a deterministic Tucson fixture covering headers, zeros, missing fields, embedded-space IDs, and 0.01/0.001-mm precision.
- Added a mathematical specification linking equations, implementations, and tests.
- Added an ITRDB audit script, candidate register, and six-file integration audit.
- Kept Yoshino metadata deliberately minimal: raw keycode plus the original paper as the scientific authority, with specimen means reconstructed by the workflow.
