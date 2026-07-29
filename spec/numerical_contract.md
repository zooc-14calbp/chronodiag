# Numerical contract v04

## Algebra

For finite `n >= 1`, a finite correlation `rbar` in `[-1, 1]`, and a positive
model denominator,

\[
EPS(n,r)=\frac{nr}{1+(n-1)r}.
\]

For `1 <= n_sub <= n_full`, SSS is computed directly:

\[
SSS(n_{sub}\mid n_{full},r)
=\frac{n_{sub}\{1+(n_{full}-1)r\}}
{n_{full}\{1+(n_{sub}-1)r\}}.
\]

The direct formula is mandatory because at `r = 0`, SSS is well defined as
`n_sub/n_full`, while the numerical ratio `EPS(n_sub)/EPS(n_full)` is `0/0`.

For a valid `N x N` correlation matrix `R`,

\[
D_{eff}(R)=\frac{(\mathrm{tr}R)^2}{\mathrm{tr}(R^2)}
          =\frac{N^2}{\sum_{ij}r_{ij}^2}.
\]

## Composition diagnostics

Composition functions operate on one value per calendar year and biological unit.
Duplicate year–unit observations are errors.

For group `G` containing `g` units and `min_depth_fraction = f`, the chronology is
available in year `t` only when

\[
n_t \ge \lceil fg\rceil.
\]

The default `f = 1` gives a complete-period comparison. Arithmetic mean and median
chronologies are supported and must be labelled in output.

`subset_stability()` reports correlation and squared correlation between an
explicit subset and the full chronology. Its `shared_unit_count` must equal the
subset size. The squared correlation is labelled an empirical counterpart of SSS
only under the exchangeable arithmetic-mean model on a fixed interval.

`split_half_stability()` reports correlation and squared correlation between
explicit disjoint A/B chronologies. Its `shared_unit_count` must equal zero. By
default A and B must exhaust the full biological-unit pool and their sizes may
differ by at most one unit.

Every syntactically valid comparison returns a row with one of:

- `ok`;
- `insufficient_overlap`;
- `constant_chronology`.

The last two states have missing correlation outputs; malformed memberships are
hard errors.

## Explicit tree bootstrap

A bootstrap replicate is represented by positive integer multiplicities `m_i`.
The nominal draw count is

\[
M=\sum_i m_i.
\]

The default contract requires `M = N`. An `m`-out-of-`N` design must be enabled
explicitly. For an arithmetic-mean bootstrap chronology, repeated draws are
integer weights:

\[
C_b(t)=
\frac{\sum_i m_i x_i(t)I_i(t)}{\sum_i m_iI_i(t)}.
\]

The weighted sample-depth rule uses the nominal-draw depth

\[
d_b(t)=\sum_i m_iI_i(t)
\]

and requires `d_b(t) >= ceil(f M)`. This is algorithmic bootstrap depth, not a
claim of additional biological independence.

Multiplicity concentration is reported as the Kish effective unit count,

\[
N_{Kish}=\frac{(\sum_i m_i)^2}{\sum_i m_i^2}.
\]

For positive integer multiplicities,

\[
1 \le N_{Kish}\le N_{unique}\le M.
\]

Bootstrap/full correlation is descriptive because the two chronologies share the
unique represented units. Direct chronology differences are the primary outputs:
mean difference, RMS difference, RMS standardized by the sample standard deviation
of the full chronology over the same years, mean and maximum absolute difference,
the year and sign of the maximum, and the variance ratio.

## Leave-one-tree-out influence

For complete arithmetic-mean chronologies,

\[
C_{-i}(t)-C_N(t)=\frac{C_{-i}(t)-x_i(t)}{N}.
\]

A LOTO chronology shares `N - 1` units with the full chronology. Its correlation
is therefore descriptive and must not be treated as the primary influence metric.
The direct difference outputs are the same as for bootstrap/full comparison.

Changes in `rbar`, EPS, and effective dimension are evaluated on one fixed period
where all units in the original pool are complete. The period must not change after
omission. If the complete period is insufficient, a series is constant, or the
correlation matrix is invalid, structural outputs remain missing and
`structure_status` records the reason. Pairwise-complete substitution and silent
matrix repair are prohibited.

## Required invariants

- `SSS(n_full, n_full, rbar) == 1`;
- `SSS(n_sub, n_full, 0) == n_sub/n_full`;
- `1 <= effective_dimension(R) <= N` for a valid correlation matrix;
- identity matrix gives `effective_dimension == N`;
- all-ones matrix gives `effective_dimension == 1`;
- column/row permutation does not change effective dimension;
- an exchangeable matrix gives
  `N / (1 + (N - 1) * rbar^2)` for effective dimension;
- for fixed `N` and mean off-diagonal correlation `rbar`,
  `effective_dimension <= N / (1 + (N - 1) * rbar^2)`;
- subset-to-full output has `shared_unit_count == group_a_size`;
- disjoint split-half output has `shared_unit_count == 0`;
- `squared_correlation == correlation^2` whenever `status == "ok"`;
- complete-period frozen composition fixtures match across R and Python;
- the exact exchangeable composition fixture satisfies
  `subset_to_full_r2 == SSS(3 | 6, 0.3)` and
  `equal_half_disjoint_r == EPS(3, 0.3)`;
- every bootstrap multiplicity is a positive integer and replicate–unit rows are
  unique;
- a default tree-bootstrap replicate has `nominal_draw_count == full_pool_size`;
- `1 <= kish_effective_unit_count <= unique_unit_count <= nominal_draw_count`;
- an identity bootstrap has correlation one and all direct differences equal zero;
- bootstrap/full output has `shared_unit_count == unique_unit_count`;
- LOTO output has `shared_unit_count == full_pool_size - 1`;
- every LOTO row uses the same `structure_n_years` and full structural baseline;
- the frozen influence fixture identifies T06 as the largest RMS LOTO influence
  while its reduced/full correlation remains above 0.99.

## Error handling

Functions must reject:

- non-finite scalar arguments;
- correlations outside `[-1, 1]`;
- `n < 1`;
- `n_sub > n_full`;
- non-positive equicorrelation denominators;
- non-square matrices;
- materially asymmetric matrices;
- matrices whose diagonal is not approximately one;
- matrices with materially negative eigenvalues when a correlation matrix is required;
- duplicate year–unit observations;
- unknown or empty biological-unit IDs;
- overlapping split-half groups;
- non-exhaustive split halves unless explicitly allowed;
- unbalanced split halves unless explicitly allowed;
- sample-depth fractions outside `(0, 1]`;
- missing, zero, negative, fractional, or non-finite bootstrap multiplicities;
- duplicate bootstrap replicate–unit rows;
- bootstrap draw counts different from `N` unless explicitly allowed;
- LOTO input containing fewer than two biological units.

Negative plug-in EPS values must not be clipped or re-labelled as reliability.
No function may silently repair a matrix. A future optional repair function may
be added, but its use must be explicit and recorded in output metadata.

## Deterministic tolerances

- scalar algebra: `1e-12`;
- matrix/eigenvalue summaries: `1e-10`;
- frozen-index composition outputs: `1e-12`;
- frozen bootstrap and LOTO outputs: `1e-12`;
- future stochastic distribution summaries: tolerance defined by the simulation
  contract rather than bitwise RNG equality.
