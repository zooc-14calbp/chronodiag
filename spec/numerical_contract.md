# Numerical contract v02

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
  `effective_dimension <= N / (1 + (N - 1) * rbar^2)`.

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
- matrices with materially negative eigenvalues when a correlation matrix is
  required.

Negative plug-in EPS values must not be clipped or re-labelled as reliability.
No function may silently repair a matrix.  A future optional repair function may
be added, but its use must be explicit and recorded in output metadata.

## Deterministic tolerances

- scalar algebra: `1e-12`;
- matrix/eigenvalue summaries: `1e-10`;
- frozen-index resampling outputs: `1e-8`.
