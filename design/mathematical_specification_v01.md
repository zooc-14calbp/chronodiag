# chronodiag mathematical specification v01

## 1. Purpose and scope

This document is the language-neutral mathematical contract for the R and Python
implementations of `chronodiag`.  It defines the quantities to be computed and
distinguishes model-relative reliability from structural, compositional, and
temporal diagnostics.  It does **not** define a universal chronology-strength
score or a pass/fail threshold.

The elementary derivations in Sections 2–8 require only means, variances,
covariances, correlations, and algebra.  Matrix notation is introduced only in
Sections 9–10.

## 2. Notation and replication unit

Let `t = 1, ..., T` index years and `i = 1, ..., N` index biological replication
units (normally trees or specimens, not individual radii/cores unless explicitly
stated).  A standardized series is denoted by

\[
x_{i,t}.
\]

For a set `A` of `n` trees,

\[
\bar x_{A,t}=\frac{1}{n}\sum_{i\in A}x_{i,t}.
\]

The full chronology is

\[
\bar x_{N,t}=\frac{1}{N}\sum_{i=1}^{N}x_{i,t}.
\]

All functions must record the unit used for resampling and averaging.  Multiple
cores must either be aggregated to a tree/specimen mean by a declared rule or
retained in an explicitly hierarchical analysis.

## 3. Exchangeable common-signal model

The elementary model is

\[
x_{i,t}=s_t+\varepsilon_{i,t},
\]

with

\[
\operatorname{Var}(s_t)=\sigma_s^2,\qquad
\operatorname{Var}(\varepsilon_{i,t})=\sigma_\varepsilon^2,
\]

\[
\operatorname{Cov}(\varepsilon_{i,t},\varepsilon_{j,t})=0\quad(i\ne j),
\qquad
\operatorname{Cov}(s_t,\varepsilon_{i,t})=0.
\]

All trees have equal loadings and equal residual variances.  The population
correlation between two distinct tree series is

\[
\rho=\frac{\sigma_s^2}{\sigma_s^2+\sigma_\varepsilon^2},
\qquad 0\le \rho\le 1.
\]

Temporal independence is not required for the population identities below, but
temporal autocorrelation affects finite-sample uncertainty.

## 4. EPS as reliability of the mean chronology

The full mean is

\[
\bar x_{N,t}=s_t+\bar\varepsilon_{N,t}.
\]

Because independent residual components are averaged,

\[
\operatorname{Var}(\bar\varepsilon_{N,t})=\frac{\sigma_\varepsilon^2}{N},
\]

and therefore

\[
\operatorname{Var}(\bar x_{N,t})
=\sigma_s^2+\frac{\sigma_\varepsilon^2}{N}.
\]

Reliability relative to the specified common signal is

\[
\operatorname{Rel}(\bar x_N)
=\frac{\sigma_s^2}
{\sigma_s^2+\sigma_\varepsilon^2/N}.
\]

Substituting

\[
\frac{\sigma_\varepsilon^2}{\sigma_s^2}=\frac{1-\rho}{\rho}
\]

gives

\[
\boxed{
\operatorname{EPS}_N(\rho)
=\frac{N\rho}{1+(N-1)\rho}
}.
\]

The empirical plug-in form replaces `rho` by an explicitly defined estimate
`rbar`.  Negative plug-in values must not be clipped to zero or interpreted as
negative reliability; they must be reported as evidence that the assumed
one-common-signal model is unsupported or that sampling/preprocessing effects are
material.

### 4.1 Inverting the EPS equation

For a chosen target `q` with `0 < q < 1`, the continuous sample size required to
reach that EPS is

\[
N^*(q,\rho)
=\frac{q(1-\rho)}{\rho(1-q)},
\]

with the operational minimum obtained by rounding upward.  This is a joint
function of correlation and replication, not a threshold on either alone.

## 5. SSS as subset-to-full squared correlation

Let `A` be a subset of `n` trees contained in the full set of `N` trees.  Under
exchangeability,

\[
\operatorname{Var}(\bar x_A)
=\frac{1+(n-1)\rho}{n},
\]

and

\[
\operatorname{Var}(\bar x_N)
=\frac{1+(N-1)\rho}{N}.
\]

Because the subset is contained in the full chronology,

\[
\operatorname{Cov}(\bar x_A,\bar x_N)
=\operatorname{Var}(\bar x_N).
\]

Hence

\[
\operatorname{Corr}^2(\bar x_A,\bar x_N)
=\frac{\operatorname{Var}(\bar x_N)}
{\operatorname{Var}(\bar x_A)},
\]

and

\[
\boxed{
\operatorname{SSS}_{n\mid N}(\rho)
=\frac{n\{1+(N-1)\rho\}}
{N\{1+(n-1)\rho\}}
}.
\]

For `rho > 0`,

\[
\operatorname{SSS}_{n\mid N}
=\frac{\operatorname{EPS}_n}{\operatorname{EPS}_N}.
\]

At `rho = 0`, the EPS ratio is formally `0/0`, but SSS itself remains well
defined by the direct correlation derivation:

\[
\operatorname{SSS}_{n\mid N}(0)=\frac{n}{N}.
\]

**Implementation requirement:** `sss_value()` must use the direct SSS formula,
not numerical division of two EPS values.

### 5.1 Inverting the SSS equation

For a chosen relative-reliability target `q`,

\[
\boxed{
n^*_{\min}(q;N,\rho)
=\frac{qN(1-\rho)}
{1+(N-1)\rho-qN\rho}
}.
\]

The threshold `q` is an application-dependent tolerance for loss relative to the
finite full chronology.

## 6. Overlapping and disjoint subset diagnostics

SSS compares a subset with a full chronology that contains the same subset.
This shared-series structure mechanically increases the correlation.  A disjoint
comparison answers a different question.

Let `A` and `B` be non-overlapping subsets of sizes `n` and `m`.  Then

\[
\operatorname{Cov}(\bar x_A,\bar x_B)=\rho,
\]

and

\[
\boxed{
\operatorname{Corr}(\bar x_A,\bar x_B)
=\frac{\rho\sqrt{nm}}
{\sqrt{\{1+(n-1)\rho\}\{1+(m-1)\rho\}}}
}.
\]

For `rho >= 0`, this can be written as

\[
\operatorname{Corr}(\bar x_A,\bar x_B)
=\sqrt{\operatorname{EPS}_n\operatorname{EPS}_m}.
\]

For equal halves (`n=m`),

\[
\operatorname{Corr}(\bar x_A,\bar x_B)=\operatorname{EPS}_n.
\]

The Spearman–Brown adjustment of an equal-half correlation recovers the
reliability expected for the combined `2n`-tree mean under the exchangeable
model:

\[
\frac{2\operatorname{EPS}_n}{1+\operatorname{EPS}_n}
=\operatorname{EPS}_{2n}.
\]

`subset_to_full` and `subset_to_complement` outputs must therefore be labelled
separately and must never be pooled as estimates of the same quantity.

## 7. Leave-one-tree-out correlation and influence

A leave-one-tree-out chronology shares `N-1` of `N` trees with the full
chronology.  Its squared correlation with the full chronology is

\[
\operatorname{Corr}^2(\bar x_{-i},\bar x_N)
=\operatorname{SSS}_{N-1\mid N}.
\]

This value approaches one rapidly as `N` increases and is therefore a weak
influence diagnostic by itself.

The exact difference is

\[
\bar x_{-i}-\bar x_N
=\frac{\bar x_{-i}-x_i}{N}.
\]

Consequently, `leave_one_tree_out()` must emphasize direct changes, such as
RMS difference, maximum absolute difference, changes in covariance summaries,
and changes in downstream diagnostic quantities.  The shared-data correlation
may be reported only as a descriptive quantity with an explicit overlap warning.

## 8. Correspondence with standardized coefficient alpha

For standardized series with mean off-diagonal correlation `rbar`, standardized
coefficient alpha is

\[
\alpha_{\mathrm{std}}
=\frac{N\bar r}{1+(N-1)\bar r}.
\]

This is algebraically identical to the conventional plug-in EPS formula for any
correlation matrix summarized by the same `rbar`.  The equality is algebraic;
its interpretation as reliability requires an appropriate measurement model.
Unequal loadings, correlated residuals, multidimensional structure, unequal
overlap, and estimated detrending functions can invalidate the parallel-test
interpretation without changing the numerical equality.

## 9. Effective covariance dimension

Let `R` be an `N x N` valid correlation matrix with eigenvalues
`lambda_1, ..., lambda_N`.  Define

\[
\boxed{
D_{\mathrm{eff}}(R)
=\frac{\left(\sum_k\lambda_k\right)^2}
{\sum_k\lambda_k^2}
=\frac{N^2}{\sum_{i,j}r_{ij}^2}
}.
\]

It satisfies

\[
1\le D_{\mathrm{eff}}\le N.
\]

For an exchangeable matrix,

\[
D_{\mathrm{eff}}
=\frac{N}{1+(N-1)\rho^2}.
\]

Writing the mean squared off-diagonal correlation as

\[
\overline{r^2}
=\frac{1}{N(N-1)}\sum_{i\ne j}r_{ij}^2,
\]

we have

\[
D_{\mathrm{eff}}
=\frac{N}{1+(N-1)\overline{r^2}}.
\]

For a fixed mean off-diagonal correlation `rbar`, Jensen's inequality gives

\[
\overline{r^2}\ge \bar r^2,
\]

and therefore

\[
\boxed{
D_{\mathrm{eff}}(R)
\le \frac{N}{1+(N-1)\bar r^2}
}.
\]

Equality holds when all off-diagonal correlations are equal.  Thus the
exchangeable matrix has the maximum effective dimension among matrices with the
same `N` and mean correlation; heterogeneous or clustered correlation structure
can retain the same EPS while lowering effective dimension.  This is a structural
contrast, not a quality ranking.

`effective_dimension()` may be evaluated only for a symmetric positive
semidefinite correlation matrix.  A pairwise-complete matrix assembled from
different overlaps must not be silently repaired or used without disclosure.

## 10. General non-exchangeable reliability model

For

\[
\boldsymbol{x}_t=\boldsymbol{\lambda}s_t+\boldsymbol{\varepsilon}_t,
\qquad
\operatorname{Cov}(\boldsymbol{\varepsilon}_t)=\boldsymbol{\Psi},
\]

and chronology

\[
c_t=\boldsymbol{w}^{\mathsf T}\boldsymbol{x}_t,
\]

reliability relative to `s_t` is

\[
\operatorname{Rel}(c)
=\frac{\sigma_s^2(\boldsymbol{w}^{\mathsf T}\boldsymbol{\lambda})^2}
{\sigma_s^2(\boldsymbol{w}^{\mathsf T}\boldsymbol{\lambda})^2
+\boldsymbol{w}^{\mathsf T}\boldsymbol{\Psi}\boldsymbol{w}}.
\]

`N` and mean correlation alone do not identify this quantity.  The same EPS can
therefore coexist with unequal loadings, clusters, influential trees, correlated
residuals, or temporal changes.

## 11. Numerical and reporting contract

1. Never clip negative `rbar` or negative plug-in EPS to zero.
2. Report model inadmissibility when `1+(N-1)rbar <= 0`.
3. Compute SSS directly so that `rbar = 0` is supported.
4. Preserve zero-width rings when reading Tucson RWL data unless the user
   explicitly requests edge-zero conversion.
5. Distinguish missing values, zero-width rings, and precision stop markers.
6. Use complete-period correlation matrices for eigenvalue diagnostics by
   default.
7. Report tree/core aggregation, analysis interval, overlap rule, and minimum
   overlap.
8. Do not return a universal composite score or automatic pass/fail label.

## 12. Equation–implementation–test traceability

| Concept | R function | Python function | Minimum deterministic test |
|---|---|---|---|
| EPS | `eps_value()` | `eps_value()` | `N=18`, EPS `=0.85` implies `rbar=0.2394366197...` |
| SSS | `sss_value()` | `sss_value()` | `rbar=0`, `n=5`, `N=20` gives `0.25` |
| Effective dimension | `effective_dimension()` | `effective_dimension()` | identity matrix gives `N`; same-EPS fixtures differ |
| Subset-to-full | `subset_stability()` | `subset_stability()` | fixed subset indices, theoretical exchangeable expectation |
| Disjoint split-half | `split_half_stability()` | `split_half_stability()` | fixed partitions; equal halves give EPS of one half |
| LOTO | `leave_one_tree_out()` | `leave_one_tree_out()` | fixed influential-tree fixture and direct-difference outputs |
| RWL input | `read_tucson_rwl()` | `read_tucson_rwl()` | ca533: 34 series, 23,276 observations, 234 zero rings |
