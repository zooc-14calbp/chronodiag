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

Composition functions evaluate **explicit, auditable memberships**. Random
partition or bootstrap generation is kept separate so that identical tree
compositions can be evaluated in R and Python.

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

### Explicit tree bootstrap

The bootstrap evaluator accepts positive integer multiplicities:

```r
boot <- data.frame(
  replicate_id = c("B01", "B01", "B01"),
  unit_id = c("T01", "T02", "T03"),
  multiplicity = c(2, 1, 1)
)
bootstrap_results <- tree_bootstrap(tree_level, boot, require_full_size = FALSE)
```

```python
from chronodiag import tree_bootstrap

bootstrap_results = tree_bootstrap(
    tree_level,
    {"B01": ["T01", "T01", "T02", "T03"]},
    require_full_size=False,
)
```

Ordinary `N`-out-of-`N` draws are required by default. Output retains the nominal
draw count, unique and omitted units, multiplicities, Kish effective unit count,
and direct chronology differences. Bootstrap/full correlation is descriptive
because the two chronologies share represented trees.

### Leave-one-tree-out influence

```r
loto_results <- leave_one_tree_out(tree_level)
```

```python
from chronodiag import leave_one_tree_out

loto_results = leave_one_tree_out(tree_level)
```

A reduced chronology shares `N - 1` trees with the full chronology, so its
correlation is mechanically high for moderate or large `N`. The primary outputs
are RMS, mean absolute, and maximum absolute chronology differences. On one fixed
complete period, the result also reports changes in mean interseries correlation,
EPS, and effective dimension; it never silently switches to pairwise-complete
covariance estimation.

The full contracts are in:

- `design/composition_diagnostics_spec_v01.md`;
- `design/influence_diagnostics_spec_v01.md`.

## Controlled same-EPS Stage A experiment

The first manuscript-linked experiment holds `N = 18`, mean interseries
correlation `0.2394366`, and the **conventional plug-in EPS = 0.85** constant
across four frozen covariance structures:

- exchangeable;
- unequal one-factor loadings;
- a high-loading minority;
- two equal clusters.

The script exhaustively evaluates every non-empty proper subset and every unique
balanced split half, then applies one frozen set of 1,000 explicit tree-bootstrap
draws and all 18 LOTO omissions:

```text
python scripts/run_same_eps_stage_a_v01.py
```

The end-to-end validator regenerates the analysis in a temporary directory and
compares the scientific outputs with the frozen tables:

```text
python scripts/validate_same_eps_stage_a_v01.py
```

The main controlled result is that identical `N`, mean correlation, and plug-in
EPS do not fix covariance architecture or compositional replication. Exchangeable
nine-versus-nine split correlations are identically 0.739, whereas the
two-cluster structure ranges from 0.014 to 0.803. Effective dimension declines
from 9.12 to 6.00 across those two structures.

Frozen tables, hashes, and figure prototypes are in:

```text
analysis/same_eps_stage_a_v01/
```

The numerical and interpretive contracts are:

- `design/same_eps_stage_a_spec_v01.md`;
- `design/same_eps_stage_a_results_v01.md`;
- `design/same_eps_stage_a_figure_captions_v01.md`.

## Repository design

- `spec/`: language-neutral scientific, input, and output contracts;
- `shared_fixtures/`: deterministic numerical and file-format fixtures;
- `r/chronodiag/`: native R package;
- `python/`: native Python package;
- `scripts/`: data-audit, analysis, and manuscript-figure scripts;
- `analysis/`: frozen machine-readable manuscript analyses;
- `data/yoshino/`: public Yoshino worked-example location;
- `data/itrdb/`: DOI and hash manifest for public ITRDB validation datasets;
- `design/`: manuscript outline, mathematical specifications, and analysis plans;
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

- The integrated Python suite contains 50 baseline, composition, bootstrap, LOTO,
  and same-EPS Stage A tests and is enforced on Python 3.10, 3.12, and 3.13.
- A separate integration audit successfully reads six full ITRDB files: ca533,
  germ011, germ012, germ6, fran003, and brit042.
- ca533 retains all 23,276 observations, including 234 zero-width rings.
- A Python wheel builds successfully and reproduces the six-file ITRDB audit when
  installed from the wheel.
- The R package passes shared metric, RWL I/O, subset/split-half,
  bootstrap/LOTO, and same-EPS Stage A parity checks.
- `R CMD check` is enforced with errors on warnings.
- Frozen R/Python fixtures cover subset-to-full, disjoint split-half, explicit
  tree bootstrap, LOTO influence, exact EPS/SSS linkage, and the controlled
  same-EPS experiment.

## Immediate next milestone

1. run and freeze the empirical Stage B profile for ca533, GERM012, and FRAN003;
2. distinguish positive-control agreement from structural/compositional
   discordance in those empirical cases;
3. install and validate the Yoshino raw keycoded measurements before analysis;
4. integrate the empirical results into the main text and Supplement; and
5. assess whether temporal or frequency-dependent extensions are required by
   the empirical findings.


## Empirical ITRDB Stage B

The first empirical composition-sensitivity workflow uses three DOI-addressable
ITRDB datasets with frozen biological-unit maps and value-blind complete-cohort
selection:

- `ca533`: 17 trees, 27 selected cores, 1727–1968;
- `GERM012`: 21 trees, 42 selected cores, 1873–1944;
- `FRAN003`: 25 one-series trees, 1819–1976.

Core/radius series are detrended in R with `dplR`; the primary method is
`AgeDepSpline`, with `ModNegExp` and a 30-year spline retained as preprocessing
sensitivities. Cores are averaged within tree, and tree means are z-scored over
one fixed complete period before composition diagnostics.

The workflow retains explicit subset, disjoint split-half, bootstrap, and LOTO
memberships and writes deterministic manifests:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass `
  -File ".\scripts\run_itrdb_empirical_stage_b_v01.ps1" `
  -ProjectRoot "D:\Dropbox\papers_2026\03f_eps"
```

The raw ITRDB files are not silently repackaged in the public repository. Their
DOIs and expected SHA-256 values are recorded under `data/itrdb/`.
