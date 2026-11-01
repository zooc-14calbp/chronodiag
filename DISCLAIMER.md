# Research software scope

chronodiag and its scripts support the reported calculations and their
inspection. The documented checks include analytical reference cases,
deterministic R/Python fixtures, unit and package checks, raw-input hash and
format checks, and frozen-output comparisons. Stage A reproduction regenerates
the controlled analysis. Stage B/C reproduction checks archived inputs and
cohorts, then regenerates memberships and diagnostic outputs from archived
preprocessed series; it does not rerun R/dplR detrending from raw measurements.

The authors reviewed the mathematical derivations and reported outputs and
take responsibility for the analyses and interpretations. Checks concern the
specified inputs, outputs and implementations. They do not establish error-free
operation or constitute a comprehensive audit of third-party libraries.

The software is distributed under GPL-3.0-or-later without warranty, as stated
in LICENSE. Users should inspect inputs, check assumptions and evaluate
suitability for their application. Report defects with a reproducible example.
The Git tag and commit identify this snapshot. An equivalent open IZUMI deposit
is prepared; no archival DOI is asserted until assigned by the repository.
