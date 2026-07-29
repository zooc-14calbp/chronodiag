# Dual-language implementation specification v01

## 1. Canonical source of truth

Neither R nor Python is the canonical scientific definition. The canonical source
of truth is the language-neutral specification in `spec/`, together with frozen
numerical fixtures in `shared_fixtures/`.

The manuscript reproduction workflow may use R for figure production, while the
Python implementation must reproduce the same reference values and diagnostic
summaries within declared numerical tolerances.

## 2. Scope of v0.1

### Included

- explicit tree/core hierarchy;
- mean interseries correlation with documented overlap policy;
- EPS and SSS;
- effective covariance dimension;
- subset-to-full stability distributions;
- non-overlapping split-half stability distributions;
- tree bootstrap;
- leave-one-tree-out influence;
- non-overlapping period and moving-window coherence diagnostics;
- machine-readable diagnostic-profile output.

### Excluded from the core release

- automatic detrending or standardization;
- climate-response analysis;
- reconstruction calibration or verification;
- a composite chronology-strength score;
- universal thresholds or automatic PASS/FAIL labels;
- black-box machine learning.

The first public release should accept already standardized series and treat
standardization sensitivity as an upstream analysis layer.

## 3. Input contract

The core numerical representation is a `T x N` matrix of standardized tree-level
series, with years in rows and trees in columns. Convenience adapters may accept
long-format data or RWL-like wide tables.

Long-format minimum fields:

- `year` — integer;
- `series_id` — unique measurement-series identifier;
- `tree_id` — biological resampling unit;
- `value` — standardized value.

Optional fields:

- `core_id`;
- `site_id`;
- `measurement`;
- `standardization_id`.

Multiple cores must not silently become independent trees. The adapter must either:

1. aggregate cores to one tree-level series using a declared rule; or
2. preserve the hierarchy for a hierarchical resampling routine.

## 4. Correlation and overlap policy

Two distinct correlation products should be supported and labelled:

1. **common-period correlation matrix** — required by default for eigenvalue and
   covariance-architecture diagnostics;
2. **pairwise-overlap mean correlation** — optional descriptive input for EPS/SSS.

A pairwise-complete matrix may fail to be positive semidefinite because different
pairs use different years. Effective dimension must therefore not be computed from
such a matrix unless it is explicitly projected to the nearest valid correlation
matrix and that intervention is reported.

Required metadata:

- overlap rule;
- minimum pairwise overlap;
- analysis period;
- treatment of missing values;
- tree/core aggregation rule.

## 5. Core diagnostics

### 5.1 Model-relative reliability

- mean interseries correlation;
- EPS;
- SSS;
- uncertainty from tree-level resampling.

### 5.2 Covariance architecture

- complete-period correlation matrix;
- eigenvalue spectrum;
- effective dimension;
- leading-mode variance fraction;
- optional loading dispersion and cluster diagnostics.

### 5.3 Composition stability

- random subset-to-full correlation distribution, retained as the empirical
  counterpart of the original SSS target;
- subset-to-complement correlation distribution, which avoids shared-series
  inflation;
- repeated split-half correlation, both raw and optionally Spearman–Brown adjusted;
- tree-bootstrap distributions based on explicit multiplicities;
- leave-one-tree-out influence using direct differences, not only correlations.

### 5.4 Sample-depth and temporal stability

- empirical stability curve by subset size;
- non-overlapping period comparisons;
- moving-window summaries used descriptively;
- sample depth reported separately from annual SSS.

## 6. No composite score

`build_diagnostic_profile()` must return a structured collection of evidence, not
one weighted score. The software may summarize observed features in words, but it
must not hide the separate diagnostics or assign a universal quality class.

## 7. Cross-language parity

Exact parity tests should use frozen:

- input matrices;
- tree identifiers;
- subset indices;
- split-half partitions;
- bootstrap multiplicities; and
- expected output tables.

Native random-number generators may be used for large simulations. Cross-language
validation of simulations should then compare distributions and Monte Carlo
summaries within tolerance, rather than requiring identical random streams.

Suggested tolerances for deterministic fixtures:

- scalar algebraic results: `1e-12`;
- matrix/eigenvalue summaries: `1e-10`;
- frozen composition, bootstrap, and LOTO outputs: `1e-12`.

## 8. Provisional public APIs

| Task | R | Python |
|---|---|---|
| EPS | `eps_value()` | `eps_value()` |
| SSS | `sss_value()` | `sss_value()` |
| effective dimension | `effective_dimension()` | `effective_dimension()` |
| core aggregation | `aggregate_cores_to_trees()` | `aggregate_cores_to_trees()` |
| covariance profile | `covariance_profile()` | `covariance_profile()` |
| subset stability | `subset_stability()` | `subset_stability()` |
| repeated split-half | `split_half_stability()` | `split_half_stability()` |
| tree bootstrap | `tree_bootstrap()` | `tree_bootstrap()` |
| leave-one-tree-out | `leave_one_tree_out()` | `leave_one_tree_out()` |
| temporal stability | `temporal_stability()` | `temporal_stability()` |
| profile assembly | `build_diagnostic_profile()` | `build_diagnostic_profile()` |

## 9. Publication and repository boundary

The public GitHub repository should contain:

- source code;
- tests and frozen fixtures;
- simulation scripts;
- public-data acquisition scripts;
- machine-readable outputs needed to reproduce paper figures;
- documentation and software citation metadata.

The Yoshino raw ring-width data may be included as a public worked example after
coauthor approval. Its measurement keycode remains the authoritative identifier,
and specimen means must be reconstructed by the workflow rather than distributed
as a second source of truth. The original Yoshino publication remains the
scientific description of the material and chronology development.

## 10. Development stages

### Milestone 0.1 — algebra and specification — completed

- EPS/SSS;
- effective dimension;
- shared deterministic fixtures;
- R/Python parity tests.

### Milestone 0.2 — composition diagnostics — implemented, empirical application pending

- subset-to-full and disjoint/subset-to-complement evaluation;
- split-half;
- explicit tree bootstrap;
- leave-one-tree-out;
- frozen cross-language membership and influence fixtures.

### Milestone 0.3 — temporal diagnostics and reporting

- period/window analysis;
- diagnostic-profile object;
- HTML/PDF report generation;
- manuscript figure pipeline.

### Milestone 1.0 — paper-linked public release

- frozen manuscript analyses;
- tagged GitHub release;
- complete citation metadata;
- reproducible environment files for both languages.
