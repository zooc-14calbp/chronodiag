# Same-EPS/different-structure simulation specification — Stage A v01

## 1. Purpose

Stage A provides the first controlled demonstration for the manuscript:

> Identical tree count, mean interseries correlation, and conventional plug-in
> EPS do not imply identical covariance architecture or identical sensitivity to
> tree composition.

The experiment is deliberately narrower than a climate-reconstruction simulation.
It evaluates chronology coherence before any climate target is introduced.  The
four scenarios use standardized tree-level series, an unweighted arithmetic-mean
chronology, a fixed complete period, and one biological series per tree.  They
therefore isolate covariance architecture and composition while holding the
traditional EPS inputs fixed.

The analysis does not propose a replacement for EPS, a composite
chronology-strength score, or a universal pass/fail threshold.

## 2. Fixed design

| Quantity | Value |
|---|---:|
| Biological units | 18 trees |
| Analysis years | 100 |
| Mean off-diagonal correlation | 0.23943661971830982 |
| Conventional plug-in EPS | 0.85 |
| Chronology constructor | unweighted arithmetic mean |
| Missing values | none |
| Subset analysis | exhaustive |
| Balanced split-half analysis | exhaustive |
| Tree-bootstrap replicates | 1,000 |
| Bootstrap design | explicit 18-out-of-18 tree multiplicities |
| Bootstrap seed | 20260730 |
| LOTO analysis | all 18 omissions |

The target mean correlation is the inverse EPS-formula solution for \(N=18\)
and \(\widehat{\mathrm{EPS}}=0.85\).  Only the exchangeable scenario gives this
coefficient its exact parallel-model reliability interpretation.

## 3. Correlation scenarios

The frozen correlation matrices already stored in `shared_fixtures/` are used as
the scientific inputs.

### 3.1 Exchangeable

Every off-diagonal element equals

\[
\rho=0.23943661971830982.
\]

This is the positive-control case in which every subset of a given size is
equivalent.  The empirical subset-to-full correlation has zero composition
spread and agrees exactly with the exchangeable SSS expression.

### 3.2 Unequal loadings

The matrix follows a one-factor model,

\[
R_{ij}=\lambda_i\lambda_j \quad (i\ne j),
\]

with smoothly increasing positive loadings.  Residual variances retain unit
diagonals.  The mean off-diagonal correlation is held at the target value.

### 3.3 Two clusters

The 18 units form two equal groups of nine.  Within-cluster correlation is 0.50;
between-cluster correlation is 0.007824726134585196, chosen so that the overall
mean off-diagonal correlation remains the target value.

This case is particularly important because a balanced split may either mix both
clusters or align almost exactly with the latent grouping.  A single average
correlation cannot distinguish those partitions.

### 3.4 Influential-tree/loading-concentration case

Three units have loading 0.95 and fifteen have loading
0.400207481631523 in a one-factor correlation model.  The scenario represents
strong concentration of covariance loading in a small group.

The term *influence* is diagnostic-specific.  A high-loading tree can strongly
affect RBAR, EPS, or effective dimension when omitted, while a lower-loading,
residual-rich tree can move the mean chronology trajectory more.  Stage A
therefore reports influence components separately rather than ranking trees with
one score.

### 3.5 Model-relative latent reliability

For the exchangeable and positive one-factor scenarios, the equal-weight mean
reliability relative to the specified factor is also retained:

\[
\operatorname{Rel}(ar x)
=
rac{(ar\lambda)^2}
{(ar\lambda)^2+N^{-2}\sum_i(1-\lambda_i^2)}.
\]

It equals 0.85 in the parallel exchangeable case.  In the unequal-loading and
loading-concentrated cases it need not equal the conventional plug-in EPS,
despite the identical \(N\) and mean interseries correlation.  The two-cluster
case is deliberately not assigned a one-factor reliability because doing so
would require an additional scientific decision about whether global and
cluster-specific shared components belong to the target signal.

## 4. Exact finite sample

For each target correlation matrix \(R\), a centered \(T\times N\) matrix is
constructed as

\[
X=\sqrt{T-1}\,QL^\mathsf{T},
\]

where \(Q^\mathsf{T}Q=I\), every column of \(Q\) is centered, and
\(LL^\mathsf{T}=R\).  The same deterministic temporal basis \(Q\) is used in all
four scenarios.  Thus,

\[
\operatorname{cor}(X)=R
\]

to numerical tolerance, and scenario differences enter only through covariance
architecture.  This removes finite-sample drift in \(\bar r\) and EPS from the
Stage A comparison.

The exact long-form data are retained as
`analysis/same_eps_stage_a_v01/exact_same_eps_long_v01.csv.gz` for R/Python audit
checks.  The manuscript's distributional results are evaluated directly from
the known correlation matrices.

## 5. Weight-vector formulation

For chronology weights \(w_A\) and \(w_B\),

\[
\operatorname{cor}(C_A,C_B)=
\frac{w_A^\mathsf{T}Rw_B}
{\sqrt{(w_A^\mathsf{T}Rw_A)(w_B^\mathsf{T}Rw_B)}}.
\]

This expression is used for exhaustive subset and split-half calculations and
for bootstrap/full and LOTO/full correlations.

For the finite exact sample, the standardized RMS chronology difference is

\[
\mathrm{sRMS}(A,B)=
\frac{
\sqrt{\frac{T-1}{T}(w_A-w_B)^\mathsf{T}R(w_A-w_B)}
}{
\sqrt{w_B^\mathsf{T}Rw_B}
}.
\]

The factor \((T-1)/T\) matches the implemented `chronodiag` definition: RMS uses
the mean squared annual difference, while the reference chronology standard
deviation uses the usual sample denominator.

## 6. Exhaustive subset-to-full analysis

Every non-empty proper subset is evaluated.  For subset size \(n\), the number of
evaluated compositions is

\[
{18\choose n}.
\]

The analysis retains the mean, standard deviation, 1st, 5th, 25th, 50th, 75th,
95th, and 99th percentiles, and the full range of subset-to-full correlation.

The exchangeable reference is

\[
\sqrt{\mathrm{SSS}_{n\mid18}}
=
\sqrt{\frac{\mathrm{EPS}_n}{\mathrm{EPS}_{18}}}.
\]

For nonexchangeable scenarios this curve is a reference derived from the common
\(N\) and \(\bar r\), not an assertion that every empirical subset has that
correlation.

## 7. Exhaustive balanced split-half analysis

All unique 9-versus-9 partitions are evaluated.  Group-label reversal is removed
by fixing `T01` in group A, leaving

\[
{17\choose8}=24\,310
\]

unique partitions.

Under exchangeability the expected split-half correlation is

\[
\mathrm{EPS}_9=0.7391304347826088.
\]

The distribution in the other scenarios quantifies the dependence of independent
compositional replication on which trees form each half.

## 8. Tree bootstrap

The bootstrap uses one frozen set of 1,000 explicit 18-out-of-18 draws.  The
numerical primitive is

```text
replicate_id, unit_id, multiplicity
```

and is saved to
`analysis/same_eps_stage_a_v01/bootstrap_multiplicities_v01.csv.gz`.

For every replicate the analysis reports correlation with the full chronology,
standardized RMS difference, variance ratio, unique and omitted units,
multiplicity concentration, and Kish effective unit count.

Scenario-specific composition axes are retained only as explanatory simulation
metadata:

- absolute cluster draw imbalance;
- mean loading of the drawn units;
- number of high-loading draws.

They are not part of the generic `chronodiag` API.

## 9. Leave-one-tree-out

Every unit is omitted once.  The analysis reports direct chronology differences
and fixed-period changes in:

\[
\bar r,\qquad \mathrm{EPS},\qquad D_{\mathrm{eff}}.
\]

Because each reduced chronology shares 17 of 18 units with the full chronology,
LOTO/full correlation is explicitly treated as descriptive.

Exact symmetric structures can give mathematically identical LOTO extrema for
several or all units.  Stage A therefore treats values as tied when they differ
from the extremum by no more than \(10^{-12}\) on the reported metric scale.
The summary retains a lexicographically first singular unit ID only for
backward-compatible machine use, and additionally reports the full pipe-delimited
set of tied unit IDs and its count.  This prevents floating-point roundoff from
creating a spurious biologically ``most influential'' tree.

## 10. Primary hypotheses

Stage A is considered scientifically successful if the following pre-specified
patterns are observed.

1. All scenarios reproduce \(N=18\), the target \(\bar r\), and EPS = 0.85.
2. Effective dimension is lower in heterogeneous or clustered scenarios than in
   the exchangeable scenario.
3. Exchangeable subset and split-half distributions collapse to their theoretical
   values.
4. Two-cluster partitions show a pronounced lower tail, including a
   cluster-aligned split with very weak cross-half correlation.
5. Unequal/loading-concentrated scenarios show nonzero subset and LOTO influence
   dispersion despite identical EPS.
6. Bootstrap mean behavior may remain similar across scenarios, but tail width
   and variance-ratio sensitivity depend on composition.
7. LOTO influence differs by estimand; direct trajectory change and structural
   changes need not identify the same unit.

## 11. Output contract

`run_same_eps_stage_a_v01.py` writes:

```text
analysis/same_eps_stage_a_v01/
├─ analysis_manifest_v01.json
├─ scenario_summary_v01.csv
├─ subset_summary_v01.csv
├─ split_half_summary_v01.csv
├─ bootstrap_multiplicities_v01.csv.gz
├─ bootstrap_results_v01.csv.gz
├─ bootstrap_summary_v01.csv
├─ bootstrap_driver_summary_v01.csv
├─ loto_results_v01.csv
├─ loto_summary_v01.csv
├─ exact_same_eps_long_v01.csv.gz
└─ key_results_v01.csv
```

The manifest records input and output hashes, analysis settings, RNG metadata,
and software versions.  Large raw membership and chronology tables use
reproducible gzip encoding (`mtime = 0`); compact summaries remain plain CSV.
Figure code must consume these frozen tables rather than recalculate undocumented
variants.

## 12. Scope limits

Stage A does not address:

- climate-growth relationships or reconstruction skill;
- temporal changes in loadings;
- frequency-dependent coherence;
- missingness or declining sample depth;
- multiple cores per tree;
- uncertainty from estimating a correlation matrix from a short record.

Those extensions belong to later simulation stages or empirical analyses.  The
controlled Stage A result establishes only that EPS-equivalent covariance systems
can differ materially in composition stability.
