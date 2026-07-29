# Tree-bootstrap and leave-one-tree-out diagnostics specification v01

## 1. Purpose

This specification defines two composition-aware diagnostics in `chronodiag`:

- `tree_bootstrap()` — evaluation of explicit tree-level bootstrap draws;
- `leave_one_tree_out()` — direct influence analysis obtained by omitting each
  biological unit in turn.

Both procedures provide **assistive evidence about chronology composition**. They
are not replacements for EPS or SSS, do not identify the environmental meaning of
a common signal, and do not produce a universal chronology-strength score or an
automatic pass/fail decision.

## 2. Biological unit and input contract

The default resampling unit is `tree_id`. A different unit column may be supplied,
but it must represent the biological replication unit whose composition is being
varied. Multiple cores or radii from the same tree/specimen must be aggregated
upstream by a documented rule, or handled by a later explicitly hierarchical
method. Treating cores as independent trees changes the estimand.

Input is canonical long-form data with one row per:

```text
year × biological unit
```

and at least:

```text
year, tree_id, value
```

Duplicate year–unit observations, non-finite years, empty unit identifiers, and
non-numeric nonmissing values are rejected. The functions do not detrend,
standardize, prewhiten, or aggregate cores internally. Those operations and the
analysis interval are part of the upstream scientific specification.

## 3. Explicit resampling objects as the numerical primitive

The numerical engines evaluate fixed memberships. They do not generate random
bootstrap draws internally in v01. This separates:

1. **resampling design**, including the random-number generator, seed, number of
   replicates, stratification, and draw size; from
2. **diagnostic evaluation**, which converts an explicit draw into a chronology
   and computes auditable outputs.

This separation allows the identical draws to be evaluated in native R and Python
without requiring their random-number generators to produce identical streams.
Large analyses may use language-native draw generators, but the generated
membership table and seed must be retained when exact reproducibility is required.

## 4. Tree bootstrap

### 4.1 Membership table

The primary input has one row per unique biological unit represented in a
bootstrap replicate:

```text
replicate_id, unit_id, multiplicity
```

`multiplicity` is a positive integer. Repeating a replicate–unit row is an error;
counts must be combined explicitly into one row. As a convenience, a named mapping
or list may contain the literal sequence of draws, which is collapsed to the same
multiplicity representation.

For replicate `b`, let the multiplicity of unit `i` be `m_bi`. The nominal number
of draws is

\[
M_b = \sum_i m_{bi}.
\]

By default, `require_full_size = TRUE`, so

\[
M_b = N,
\]

which defines the ordinary nonparametric `N`-out-of-`N` tree bootstrap. An
`m`-out-of-`N` design is accepted only when explicitly requested by setting
`require_full_size = FALSE`; the realized draw count remains in every result row.

### 4.2 Bootstrap chronology

For the arithmetic-mean constructor, the bootstrap chronology in year `t` is

\[
C_b(t)
=
\frac{\sum_i m_{bi}x_i(t)I_i(t)}
     {\sum_i m_{bi}I_i(t)},
\]

where `I_i(t)` is one when unit `i` has a finite observation and zero otherwise.
Thus repeated draws act as integer weights; they are not falsely labelled as
additional independent trees.

For the median constructor, observations are conceptually repeated according to
multiplicity before the median is evaluated. Model-defined EPS/SSS relationships
apply to arithmetic means, not automatically to medians.

### 4.3 Sample-depth rule

For a bootstrap sample with nominal draw count `M_b`, the algorithmic weighted
depth in year `t` is

\[
d_b(t)=\sum_i m_{bi}I_i(t).
\]

A value is retained only when

\[
d_b(t) \geq \lceil fM_b\rceil,
\]

where `f = min_depth_fraction`. The default `f = 1` requires every nominal draw to
be observed. Because duplicated draws do not add biological independence,
`mean_depth_bootstrap` and `min_depth_bootstrap` are algorithmic draw depths, not
counts of unique trees. The result therefore also reports unique, omitted, and
duplicated units and the Kish effective unit count.

### 4.4 Kish effective unit count

Multiplicity concentration is summarized by

\[
N_{\mathrm{Kish},b}
=
\frac{\left(\sum_i m_{bi}\right)^2}
     {\sum_i m_{bi}^2}.
\]

For positive multiplicities,

\[
1 \leq N_{\mathrm{Kish},b}
\leq N_{\mathrm{unique},b}
\leq M_b.
\]

It equals the number of unique represented units when all represented units have
equal multiplicity, and equals `N` for the identity draw in an `N`-out-of-`N`
bootstrap. It is a concentration descriptor, not a replacement for the biological
tree count or a reliability coefficient.

### 4.5 Comparison with the full chronology

Each bootstrap chronology is compared with the chronology constructed from every
available biological unit using the same constructor and sample-depth rule. The
bootstrap and full chronologies share all unique units represented in the draw.
Consequently, their correlation is descriptive and can remain very high even when
the direct chronology difference is non-negligible.

Primary outputs are therefore direct differences over common valid years:

\[
\Delta_b(t)=C_b(t)-C_N(t),
\]

including:

- mean difference;
- root-mean-square difference;
- RMS difference standardized by the sample standard deviation of the full
  chronology over the same years;
- mean absolute difference;
- maximum absolute difference, its calendar year, and its sign;
- bootstrap/full chronology variance ratio.

Correlation and squared correlation are retained with `shared_unit_count`, but are
not interpreted as independent replication or as bootstrap reliability.

### 4.6 Replicate-level rather than pre-summarized output

`tree_bootstrap()` returns one row per explicit replicate. It does not impose a
confidence level, percentile convention, or universal threshold. Distributional
summaries such as median, interquartile range, 5–95% range, or percentile
intervals are downstream analysis choices and must report the number and design of
replicates.

## 5. Leave-one-tree-out influence

### 5.1 Chronology comparison

For every biological unit `i`, `leave_one_tree_out()` constructs:

- the full chronology `C_N(t)`;
- the reduced chronology `C_{-i}(t)` from the remaining `N-1` units.

The two chronologies share `N-1` units. Their correlation is therefore expected to
approach one rapidly as `N` increases and is not, by itself, evidence that unit `i`
has negligible influence.

For complete arithmetic-mean chronologies,

\[
C_{-i}(t)-C_N(t)
=
\frac{C_{-i}(t)-x_i(t)}{N}.
\]

This identity shows both why the difference shrinks with `N` and why direct
chronology changes are more informative than the shared-data correlation.

The same direct metrics used for bootstrap comparisons are reported for every
omitted unit. If several years share the same maximum absolute difference, the
first calendar year in ascending order is recorded, giving deterministic output.

### 5.2 Covariance-structure deltas

LOTO also reports changes in:

- mean off-diagonal correlation (`rbar`);
- plug-in EPS;
- correlation-matrix effective dimension.

All full and reduced covariance summaries are evaluated on **one fixed complete
period containing every unit in the original pool**. The period is not allowed to
change after omitting a tree, because changing the observations would confound
unit influence with interval selection.

If the full complete period is shorter than `min_overlap`, a series is constant, or
the resulting correlation matrix is invalid, the structural metrics remain missing
and `structure_status` records the reason. The function never silently switches to
a pairwise-complete matrix or repairs a correlation matrix.

Because reduced EPS changes both `rbar` and the number of units, `delta_eps` is a
joint change in coherence and replication. It is not interpreted as a pure effect
on either component alone. Effective dimension is a covariance-architecture
descriptor, not reliability.

## 6. Chronology constructors and overlap

The supported constructors are:

- `mean` — arithmetic mean;
- `median` — median across available values or multiplicity-expanded draws.

The constructor is recorded in every row. A comparison is assigned:

```text
ok
insufficient_overlap
constant_chronology
```

`correlation`, `squared_correlation`, and variance-ratio outputs remain missing
when their prerequisites are not met. Direct metrics are computed whenever the
minimum overlap is met, including when one chronology is constant; this preserves
useful information instead of deleting the realization.

## 7. Output interpretation

### 7.1 Tree-bootstrap metadata

Every replicate records at least:

```text
nominal_draw_count
unique_unit_count
omitted_unit_count
duplicated_unit_count
maximum_multiplicity
kish_effective_unit_count
shared_unit_count
all_units_represented
bootstrap_ids
```

`duplicated_unit_count` is the number of distinct biological units whose
multiplicity exceeds one, not the number of extra draw positions.

### 7.2 LOTO metadata

Every omission records at least:

```text
omitted_unit_id
retained_unit_count
shared_unit_count
structure_n_years
structure_status
```

The full structural baseline is repeated in each row to keep each realization
self-describing and readily auditable in long-form output.

### 7.3 Scope boundary

Neither a narrow bootstrap distribution nor small LOTO changes establish:

- climate-growth validity;
- calibration or reconstruction skill;
- temporal stationarity;
- frequency-dependent fidelity;
- representativeness beyond the observed tree pool.

Tree bootstrap conditions on the observed pool. LOTO diagnoses influence within
that pool. External validity and temporal/frequency stability remain separate
evidence domains.

## 8. Deterministic cross-language fixtures

R and Python parity uses frozen:

```text
influence_example_long_v01.csv
bootstrap_multiplicities_v01.csv
expected_tree_bootstrap_v01.csv
expected_leave_one_tree_out_v01.csv
```

The fixture contains missing values and one deliberately influential unit. It
verifies that:

1. the identity bootstrap draw exactly reproduces the full chronology;
2. literal draw lists and multiplicity tables are equivalent;
3. Kish effective unit count and membership metadata agree across languages;
4. direct bootstrap differences agree to the frozen tolerance;
5. LOTO identifies the deliberately influential unit by RMS difference even while
   reduced/full correlation remains above 0.99;
6. every LOTO structural comparison uses the same full complete period.

Frozen deterministic numerical outputs use an absolute tolerance of `1e-12` in
the parity tests. Production simulations should assess Monte Carlo convergence and
need not share random-number streams across languages.

## 9. Required invariants

- all multiplicities are positive integers;
- replicate–unit rows are unique;
- an ordinary bootstrap replicate has nominal draw count `N`;
- `1 <= kish_effective_unit_count <= unique_unit_count <= nominal_draw_count`;
- an identity bootstrap has zero direct difference and correlation one;
- `shared_unit_count = unique_unit_count` for bootstrap-to-full comparisons;
- `shared_unit_count = N - 1` for LOTO;
- LOTO full structural quantities are identical across omitted units because one
  fixed complete period is used;
- no result is converted to a pass/fail label or combined score.
