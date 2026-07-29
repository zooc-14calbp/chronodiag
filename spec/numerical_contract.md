# Numerical contract v03

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

Every syntactically valid replicate returns a row with one of:

- `ok`;
- `insufficient_overlap`;
- `constant_chronology`.

The last two states have missing correlation outputs; malformed memberships are
hard errors.

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
  `equal_half_disjoint_r == EPS(3, 0.3)`.

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
- sample-depth fractions outside `(0, 1]`.

Negative plug-in EPS values must not be clipped or re-labelled as reliability.
No function may silently repair a matrix. A future optional repair function may
be added, but its use must be explicit and recorded in output metadata.

## Deterministic tolerances

- scalar algebra: `1e-12`;
- matrix/eigenvalue summaries: `1e-10`;
- frozen-index composition outputs: `1e-12`;
- future stochastic distribution summaries: tolerance defined by the simulation
  contract rather than bitwise RNG equality.
