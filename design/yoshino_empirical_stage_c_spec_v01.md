# Yoshino empirical Stage C specification v01

## 1. Purpose

Stage C applies the same composition-aware diagnostics used in the exact Stage A
experiment and public ITRDB Stage B cases to the provenance-screened Yoshino
Japanese cedar chronology.

The scientific question is narrow:

> Does a chronology assembled from a centuries-long, intensively managed and
> market-sourced plantation system show stable collective coherence under
> transparent perturbations of specimen composition?

Stage C does not reanalyse climate--growth relationships, management releases,
or reconstruction skill.

## 2. Input level

The source file `nr_all.rwl` contains 60 specimen-mean series. Project records
indicate 72 measured radii from 60 specimens; multiple radii were averaged in the
original Yoshino workflow. The present analysis therefore treats each RWL column
as one biological specimen and never labels the file as a 72-radius raw archive.

This distinction is essential for both scientific interpretation and public data
release.

## 3. Provenance and key convention

The six-character key is preserved. The first three digits identify project
source groups 001--003, mapped to NR01--NR03. The full key is `specimen_id` and
`tree_id`. The analysis input is already a specimen mean.

The project site register identifies NR01--03 collectively as the Yoshino
Forestry Industry Area, Japanese cedar (*Cryptomeria japonica*), with permission
for use.

## 4. Value-blind complete cohort

Cohort selection uses only series start/end years. Candidate intervals satisfy

\[
T\ge 50,\qquad T\ge 3N.
\]

Candidates are ranked by maximum specimen count, then maximum interval length,
then earliest start. The primary cohort is:

\[
N=50,\qquad 1850\text{--}2004,\qquad T=155.
\]

The ratio rule is an analysis-design constraint for estimating one complete
correlation matrix; it is not a universal chronology-quality threshold.

## 5. Preprocessing

Primary:

```text
dplR::detrend(method = "AgeDepSpline")
```

Sensitivity:

```text
ModNegExp
Spline, nyrs = 30, f = 0.5
```

Each specimen-mean series is detrended over its available span, cut to the fixed
complete interval, and standardized to zero mean and unit sample standard
deviation within that interval. The chronology is an unweighted arithmetic mean
of specimens.

## 6. Frozen perturbations

- subset sizes 2 through 48, at most 250 unique subsets per size;
- 1,000 unique balanced 25-versus-25 disjoint partitions;
- 2,000 explicit 50-out-of-50 specimen bootstrap draws;
- all 50 LOTO omissions.

Memberships and multiplicities are frozen before diagnostic evaluation.

## 7. Primary outputs

- mean interseries correlation and conventional plug-in EPS;
- effective dimension and leading eigenvalue fractions;
- subset-to-full distribution by subset size;
- disjoint split-half distribution;
- bootstrap direct differences and variance ratios;
- LOTO direct chronology change and changes in RBAR, EPS, and effective dimension;
- preprocessing sensitivity.

No composite score or pass--fail label is produced.

## 8. Scope and interpretation

The bootstrap is conditional on the observed 50-specimen cohort. The fixed
complete cohort is not the full 60-specimen archive. High internal coherence does
not establish climatic meaning. A stable result would nevertheless show that a
managed, market-assembled regional archive can provide a collective annual chronology whose
coherence is not dominated by one arbitrary specimen composition.

## 9. Portable validation manifests

Frozen production outputs are stored with repository-relative paths. During
deterministic validation, regenerated memberships and results are written to an
external temporary directory. Manifest serialization therefore falls back to an
absolute POSIX path when a file is outside the repository. Numerical and
membership equality is checked from table contents; temporary manifest paths are
not treated as scientific results.


## 10. Source-byte and repository-text identity

The Yoshino RWL is an ASCII fixed-width text file. Its original source-byte
SHA-256 is retained for provenance. Because the repository enforces LF line
endings for `*.rwl`, reproducibility across Windows and Linux is verified with a
second SHA-256 after CRLF and lone-CR normalization to LF. These hashes answer
different questions and are not substituted for one another.
