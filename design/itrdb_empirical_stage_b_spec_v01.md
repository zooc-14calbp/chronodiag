# ITRDB empirical composition-sensitivity analysis — Stage B v01

## 1. Purpose

Stage B moves from the exact covariance experiment to three public empirical
tree-ring datasets with distinct methodological roles:

1. **ca533** — a canonical, strongly coherent, multiple-core benchmark;
2. **GERM012** — a paired-core hierarchy with a shorter complete period;
3. **FRAN003** — a long, single-series-per-tree complete-case example.

The aim is not to find datasets that make EPS appear inadequate.  The portfolio
includes a positive control and progressively less exchangeable or less coherent
structures.  Stage B asks whether conventional plug-in EPS, covariance
architecture, and explicit biological-composition sensitivities agree or
separate in real data.

Climate-growth relationships and reconstruction skill remain outside this stage.

## 2. Data provenance

Raw Tucson RWL files remain outside the public repository and are acquired from
DOI-addressable ITRDB records.  The project-local copies are verified against
the frozen SHA-256 values in `empirical_dataset_contract_v01.csv`.

| Dataset | DOI | Role |
|---|---|---|
| ca533 | 10.25921/f65s-sx21 | canonical multiple-core benchmark |
| GERM012 | 10.25921/az18-vr59 | paired-core hierarchy |
| FRAN003 | 10.25921/g06k-g948 | long complete-case structural example |

Zero-width rings are valid observations and are retained.  Missing fields,
zero-width rings, and Tucson stop markers remain distinct.


### 2.1 Explicit Tucson input contract

All three raw files are declared explicitly as Tucson RWL input. The empirical
workflow does not use `dplR::read.rwl(format = "auto")`, because format
auto-detection is a rough convenience rather than a provenance rule. The
dataset contract freezes:

```text
rwl_format
tucson_header
tucson_long
preserve_zero
```

ca533 is headerless; GERM012 and FRAN003 have three-line Tucson headers.
GERM012's investigator header contains a comma, so a generic comma-based
auto-detector can misclassify the file as CSV even though its ring-width body is
valid Tucson fixed-width data. Explicit format and header settings keep input
interpretation independent of header punctuation.

## 3. Biological-unit parsing

### 3.1 ca533

The six-character key is parsed as:

```text
site = characters 1–3
tree = characters 1–5
core = character 6
```

This follows the documented `dplR::read.ids(..., stc = c(3, 2, 1))` convention
used for ca533.  The 34 series represent 21 trees.

### 3.2 GERM012

The six-character key is provisionally parsed as:

```text
tree = characters 1–4
core/radius = character 5
measurement code = character 6
```

The file contains 22 exact pairs and every current TRW series ends in `3`.
The mapping is structurally compelling but is labelled as inferred rather than
source-documented.  It must not be presented as archival fact until original
documentation confirms the key convention.

### 3.3 FRAN003

Each exact series key is treated as one biological unit.  Embedded spaces are
preserved.  No substring-based tree grouping is imposed.

## 4. Value-blind complete-cohort selection

A tree's strict available interval is the intersection of all core/radius spans
assigned to that tree.  Candidate intervals use observed strict starts and ends.
The frozen cohort is selected by:

1. requiring at least 50 complete years;
2. requiring `T >= 3N`, where `T` is interval length and `N` is tree count;
3. maximizing `N`;
4. among equal `N`, maximizing `T`;
5. using the earliest start as the deterministic final tie-break.

The rule uses only series coverage and biological identifiers, not correlation,
EPS, or any composition-diagnostic result.  It is an analysis-design rule for a
complete covariance matrix, not a universal dendrochronological threshold.

The selected cohorts are:

| Dataset | Trees | Core/radius series | Fixed interval | Years |
|---|---:|---:|---:|---:|
| ca533 | 17 | 27 | 1727–1968 | 242 |
| GERM012 | 21 | 42 | 1873–1944 | 72 |
| FRAN003 | 25 | 25 | 1819–1976 | 158 |

The excluded tree IDs and the top alternative coverage cohorts remain in the
machine-readable audit.

## 5. Detrending and tree-level construction

Detrending is deliberately upstream from `chronodiag`.

### 5.1 Primary preprocessing

The primary method is the age-dependent spline in `dplR`:

```r
detrend(rwl, method = "AgeDepSpline")
```

Each selected core/radius is detrended over its full available span before the
fixed cohort interval is extracted.  No prewhitening is applied.  Indices are
formed by division, not subtraction.

### 5.2 Sensitivity methods

Two preprocessing sensitivities are retained:

```text
ModNegExp
Spline with nyrs = 30 and f = 0.5
```

The full composition analysis is initially evaluated for AgeDepSpline.  The
alternative methods are compared through the complete correlation matrix,
conventional plug-in EPS, effective dimension, and full-chronology difference.
If these structural results differ materially, the full composition diagnostics
will be repeated for the affected method rather than assumed invariant.

### 5.3 Core-to-tree aggregation

Within each method and fixed interval, detrended cores/radii from one tree are
averaged with equal weight.  Every tree then contributes one biological series.
This prevents multiple cores from being counted as independent trees.

The unscaled tree RWI is retained.  For correlation-based diagnostics, every
tree mean is additionally z-scored over the common interval using its sample
mean and sample standard deviation.  This gives each biological unit unit
variance and aligns empirical calculations with the correlation-based EPS
contract.

The chronology constructor is the unweighted arithmetic mean.  Tukey's biweight
mean is not used in this stage because the exact EPS/SSS relationships are
derived for arithmetic means.

## 6. Frozen empirical perturbation designs

Membership generation is separate from diagnostic evaluation.  The master seed
is `20260730`; dataset- and design-specific seeds are derived by SHA-256.

### 6.1 Subset-to-full

Every subset size from 2 to `N-2` is included.  At most 250 unique compositions
are retained per size; all combinations are used when fewer than 250 exist.

### 6.2 Disjoint split-half

For odd `N`, group A has `floor(N/2)` trees and group B the remainder.  One
thousand unique exhaustive partitions are retained per dataset.

### 6.3 Tree bootstrap

Two thousand explicit `N`-out-of-`N` tree-bootstrap draws are retained per
dataset as positive integer multiplicities.

### 6.4 LOTO

Every tree is omitted once.  Ties in extremal influence measures use tolerance
`1e-12` and retain all tied identifiers.

The primary diagnostics use:

```text
chronology_method = mean
min_depth_fraction = 1
min_overlap = full fixed-period length
```

Because the selected grids are complete, every valid realization uses the same
calendar interval.

## 7. Outputs

For the primary method, Stage B writes:

- baseline `rbar`, conventional plug-in EPS, effective dimension, eigenvalues,
  and off-diagonal correlation dispersion;
- realization-level subset, split-half, bootstrap, and LOTO tables;
- compact quantile summaries;
- standardization-sensitivity summaries;
- manuscript-oriented key results;
- deterministic memberships and an analysis manifest;
- color figure prototypes.

The lower tail is reported together with the median and valid replicate count.
The sampled minimum for a non-exhaustive empirical design is labelled as such
and is not interpreted as the mathematical worst case.

## 8. Pre-specified empirical expectations

Stage B is descriptive rather than confirmatory, but the portfolio was chosen
with the following expectations:

1. ca533 should operate as a positive-control chronology with high conventional
   EPS and comparatively narrow composition distributions.
2. GERM012 should reveal whether a short but core-complete paired hierarchy
   remains compositionally stable after tree aggregation.
3. FRAN003 should provide the strongest opportunity for composition tails to
   differ from the central tendency because it combines lower replication
   coherence with a long complete period.
4. Effective dimension and composition lower tails need not rank datasets in
   exactly the same order as conventional plug-in EPS.
5. High subset-to-full or LOTO-to-full correlations will not be interpreted as
   independent replication because the compared chronologies share trees.

No result is converted into a new pass/fail rule or composite score.

## 9. Scope limitations

- Inference is conditional on the observed tree pools.
- GERM012's key hierarchy remains an explicit inference pending source
  confirmation.
- Detrending uncertainty is represented by three methods, not exhausted.
- The fixed complete-period design sacrifices some available series and years in
  exchange for one auditable covariance matrix.
- Frequency dependence, temporal moving windows, climate response, and
  reconstruction skill are separate future analyses.
- Yoshino is not included until its raw keycoded measurements are installed and
  the 72-radius to 60-specimen adapter is validated.
