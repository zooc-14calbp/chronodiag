# Changelog

## Unreleased

- Fixed the project, R package, and Python package working name as `chronodiag`.
- Added `pandas` as an explicit Python runtime dependency and explicit src-layout package discovery.
- Added four shared same-EPS/different-structure parity tests; the Python suite now contains 16 passing tests.
- Verified that the Python wheel builds and reproduces the six-file ITRDB audit.
- Expanded GitHub Actions to Python 3.10/3.12/3.13 and an R documentation, parity, and `R CMD check` job.
- Corrected the R package metadata and exports for CSV/RWL input functions.
- Restored physically blank six-character Tucson data fields to missing values in the R adapter while retaining genuine zero-width rings.
- Added native R and Python implementations governed by a language-neutral
  scientific contract.
- Added wide CSV, long CSV, and Tucson RWL input adapters.
- Corrected SSS to use its direct formula, including the `rbar = 0` case.
- Added admissible-range validation for EPS and SSS.
- Corrected Python Tucson parsing for standard three-line headers, embedded-space
  IDs, internal missing fields, and both precision markers.
- Preserved zero-width rings by default; ca533 retains 234 zero rings.
- Added a deterministic Tucson fixture covering headers, zeros, missing fields,
  embedded-space IDs, and 0.01/0.001-mm precision.
- Added a mathematical specification linking equations, implementations, and
  tests.
- Added an ITRDB audit script, candidate register, and six-file integration audit.
- Added revised R tests; R execution is pending local/GitHub Actions validation.
- Kept Yoshino metadata deliberately minimal: raw keycode plus the original paper
  as the scientific authority, with specimen means reconstructed by the workflow.
