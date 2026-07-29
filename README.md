# chronodiag

`chronodiag` is a dual-language research-software project for **assistive
diagnostics of tree-ring chronology coherence**. The R and Python packages are
native implementations of the same language-neutral scientific specification.

The project is currently a validated development baseline, not a finished public
release.

## Scientific position

The software does not seek a universal replacement for EPS or SSS. It retains
EPS/SSS as model-relative reliability summaries and evaluates other evidence
domains separately:

- covariance architecture;
- tree-composition sensitivity;
- sample-depth dependence; and
- temporal stability of chronology coherence.

It does **not** assign a single chronology-strength score or automatic PASS/FAIL
classification. Climate–growth relationships and climate-reconstruction skill
are downstream questions and are outside the core scope of v0.1.

## User-facing input formats

`read_chronodiag()` accepts:

1. wide CSV: `year,key1,key2,...`;
2. long CSV: at least `year,series_id,value`; and
3. Tucson/ITRDB ring-width files (`.rwl` or `.tuc`).

All readers return canonical long-form observations with at least
`year`, `series_id`, `tree_id`, and `value`.

Series keys are preserved exactly. By default, `tree_id` equals `series_id`.
A dataset-specific parser may add `site_id`, `specimen_id`, `tree_id`, and
`core_id`; those rules remain outside the generic numerical core.

### Tucson RWL policy

- standard three-line headers and headerless files are supported;
- 8-character ID/4-character year and 7-character ID/5-character year layouts
  are supported;
- both `999` (0.01 mm) and `-9999` (0.001 mm) precision markers are supported;
- zero-width rings are preserved by default;
- missing fields, zero rings, and stop markers remain distinct;
- ambiguous precision requires an explicit user choice rather than a silent
  assumption.

## Composition-aware diagnostics

The first composition-aware API evaluates **explicit, auditable partitions**.
Random partition generation is kept separate so that the same tree compositions
can be evaluated exactly in R and Python.

### Subset versus full chronology

```r
subsets <- data.frame(
  replicate_id = c("S01", "S01", "S02", "S02"),
  unit_id = c("T01", "T02", "T03", "T04")
)
subset_results <- subset_stability(tree_level, subsets)
```

```python
from chronodiag import subset_stability

subset_results = subset_stability(tree_level, subsets)
```

The result reports both `correlation` and `squared_correlation`. The square is an
empirical counterpart of SSS only under the model conditions that define SSS.
`shared_unit_count` makes the subset/full overlap explicit.

### Disjoint split-half chronology comparison

```r
halves <- data.frame(
  replicate_id = rep("H01", 4),
  group = c("A", "A", "B", "B"),
  unit_id = c("T01", "T03", "T02", "T04")
)
split_results <- split_half_stability(tree_level, halves)
```

```python
from chronodiag import split_half_stability

split_results = split_half_stability(tree_level, halves)
```

Groups A and B must be disjoint. By default they must exhaust the biological-unit
pool and differ in size by at most one unit. `min_depth_fraction = 1` is the
complete-period default; lower values require an explicit choice and are recorded
in every output row.

The full contract is in
`design/composition_diagnostics_spec_v01.md`.

## Repository design

- `spec/`: language-neutral scientific, input, and output contracts;
- `shared_fixtures/`: deterministic numerical and file-format fixtures;
- `r/chronodiag/`: native R package;
- `python/`: native Python package;
- `scripts/`: data-audit and manuscript-reproduction scripts;
- `data/yoshino/`: public Yoshino worked-example location;
- `data/itrdb/`: DOI and hash manifest for public ITRDB validation datasets;
- `design/`: manuscript outline, mathematical specification, and simulation plan;
- `.github/workflows/`: cross-language CI.

## Yoshino example

The Yoshino raw ring-width dataset may be included as a public worked example.
Its measurement keycode is the authoritative identifier. The workflow will
reconstruct specimen-level biological units from that key rather than duplicate
extensive metadata or distribute precomputed specimen means as a second source of
truth. The original Yoshino publication remains the scientific description of
the material and chronology construction.

## Quick start

### R

```r
library(chronodiag)
wide <- read_chronodiag("measurements.csv")
rwl  <- read_chronodiag("measurements.rwl")
```

### Python

```python
from chronodiag import read_chronodiag

wide = read_chronodiag("measurements.csv")
rwl = read_chronodiag("measurements.rwl")
```

## Validation status

- The baseline Python suite passed on Python 3.10, 3.12, and 3.13 in GitHub Actions.
- The integrated Python suite contains 26 passing baseline and composition tests.
- A separate integration audit successfully read six full ITRDB files: ca533,
  germ011, germ012, germ6, fran003, and brit042.
- ca533 retained all 23,276 observations, including 234 zero-width rings.
- A Python wheel builds successfully and reproduces the six-file ITRDB audit when installed from the wheel.
- The baseline R package passed shared metric/I-O parity and `R CMD check` with
  0 errors, 0 warnings, and 0 notes.
- The composition branch adds ten Python tests plus frozen R/Python parity
  fixtures for subset-to-full and disjoint split-half diagnostics.

## Immediate next milestone

1. merge the explicit subset and split-half diagnostic contract after CI review;
2. add auditable partition generators without requiring identical R/Python RNG streams;
3. implement tree bootstrap and leave-one-tree-out influence using the same
   chronology-construction contract;
4. apply the diagnostics to the same-EPS simulations and selected empirical datasets;
5. update the manuscript and Supplement from the frozen mathematical contract.
