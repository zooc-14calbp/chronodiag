# chronodiag

`chronodiag` is a dual-language research-software project for **assistive
diagnostics of tree-ring chronology coherence**.  The R and Python packages are
native implementations of the same language-neutral scientific specification.

The project is currently a validated development baseline, not a finished public
release.

## Scientific position

The software does not seek a universal replacement for EPS or SSS.  It retains
EPS/SSS as model-relative reliability summaries and evaluates other evidence
domains separately:

- covariance architecture;
- tree-composition sensitivity;
- sample-depth dependence; and
- temporal stability of chronology coherence.

It does **not** assign a single chronology-strength score or automatic PASS/FAIL
classification.  Climate–growth relationships and climate-reconstruction skill
are downstream questions and are outside the core scope of v0.1.

## User-facing input formats

`read_chronodiag()` accepts:

1. wide CSV: `year,key1,key2,...`;
2. long CSV: at least `year,series_id,value`; and
3. Tucson/ITRDB ring-width files (`.rwl` or `.tuc`).

All readers return canonical long-form observations with at least
`year`, `series_id`, `tree_id`, and `value`.

Series keys are preserved exactly.  By default, `tree_id` equals `series_id`.
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
Its measurement keycode is the authoritative identifier.  The workflow will
reconstruct specimen-level biological units from that key rather than duplicate
extensive metadata or distribute precomputed specimen means as a second source of
truth.  The original Yoshino publication remains the scientific description of
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

- 16 repository-ready Python tests pass, including all four same-EPS covariance fixtures.
- A separate integration audit successfully read six full ITRDB files: ca533,
  germ011, germ012, germ6, fran003, and brit042.
- ca533 retained all 23,276 observations, including 234 zero-width rings.
- A Python wheel builds successfully and reproduces the six-file ITRDB audit when installed from the wheel.
- The revised R implementation follows the same contract; execution remains to
  be confirmed with `R CMD check` locally and in GitHub Actions.

## Immediate next milestone

1. run R and Python checks on the user's machine;
2. freeze cross-language input and algebra fixtures;
3. upload the validated scaffold to a private GitHub repository;
4. implement subset, split-half, bootstrap, and LOTO diagnostics; and
5. update the manuscript and Supplement from the frozen mathematical contract.
