# Same-EPS Stage A results and manuscript interpretation v01

## 1. Status

The deterministic Stage A experiment satisfies its pre-specified same-
**conventional plug-in EPS** constraint.  All four scenarios have:

\[
N=18,\qquad
\bar r=0.23943661971831,\qquad
\widehat{\mathrm{EPS}}=0.85,
\]

and the variance of the full arithmetic-mean chronology is identical
(\(0.28169014\)) because it is determined by \(N\) and the mean off-diagonal
correlation for standardized series.

The differences below therefore arise from covariance architecture and tree
composition, not from tree count, mean correlation, the conventional plug-in EPS
coefficient, chronology weighting, sample depth, or analysis-period length.

Only the exchangeable case gives the coefficient its exact parallel-model
reliability interpretation.  In the two one-factor but nonparallel cases, the
equal-weight mean reliability relative to the specified latent factor is 0.8525
(unequal loadings) and 0.8588 (loading concentration), rather than exactly 0.85.
The two-cluster stress case has no unique one-common-factor reliability without an
additional choice about which shared components count as the target signal.  This
is a secondary but useful reminder that numerical equality of the plug-in EPS
formula does not identify the same latent estimand across covariance models.

## 2. Covariance architecture

Effective dimension differs substantially:

| Scenario | Effective dimension | Off-diagonal SD |
|---|---:|---:|
| Exchangeable | 9.116 | ~0 |
| Unequal loadings | 8.705 | 0.074 |
| Influential/high-loading group | 7.848 | 0.137 |
| Two clusters | 5.999 | 0.246 |

The two-cluster scenario has the same leading-eigenvalue fraction as the
exchangeable scenario (0.282), but almost the same amount of variance is
concentrated in a second leading mode (0.274).  This is why the leading mode alone
does not reveal the dimensional contrast, whereas effective dimension does.

**Manuscript interpretation:** identical conventional plug-in EPS can coexist
with covariance systems that range from exchangeable to strongly bimodal.
Effective dimension is complementary structural information, not a replacement
reliability coefficient.

## 3. Exhaustive subset-to-full results

For a nine-tree subset, the exchangeable reference correlation is

\[
\sqrt{\mathrm{SSS}_{9\mid18}}=0.932505.
\]

Because every exchangeable subset is equivalent, the complete distribution
collapses to this value.  The nonexchangeable scenarios differ:

| Scenario | 5th percentile | Median | 95th percentile | Minimum |
|---|---:|---:|---:|---:|
| Exchangeable | 0.9325 | 0.9325 | 0.9325 | 0.9325 |
| Unequal loadings | 0.9238 | 0.9335 | 0.9401 | 0.9067 |
| Influential/high-loading group | 0.9112 | 0.9352 | 0.9456 | 0.9112 |
| Two clusters | 0.8549 | 0.9494 | 0.9494 | 0.7121 |

The two-cluster median is higher than the exchangeable reference because most
nine-tree subsets mix the two groups and approximate the full mean well.  Its
lower tail is nevertheless much worse because cluster-unbalanced compositions
represent the full covariance system poorly.

**Manuscript interpretation:** a central tendency alone can be reassuring while
the lower-tail composition risk remains large.  The distribution, especially its
lower tail, is the relevant assistive diagnostic.

## 4. Exhaustive balanced split-half results

All 24,310 unique nine-versus-nine partitions were evaluated.  Under
exchangeability, every split has

\[
r_{\mathrm{split}}=\mathrm{EPS}_9=0.739130.
\]

The results are:

| Scenario | 5th percentile | Median | 95th percentile | Minimum |
|---|---:|---:|---:|---:|
| Exchangeable | 0.7391 | 0.7391 | 0.7391 | 0.7391 |
| Unequal loadings | 0.7356 | 0.7420 | 0.7429 | 0.7185 |
| Influential/high-loading group | 0.7277 | 0.7498 | 0.7498 | 0.7277 |
| Two clusters | 0.4617 | 0.8028 | 0.8028 | 0.0141 |

The minimum two-cluster split assigns one latent cluster to each half.  The
resulting cross-half correlation is only 0.014 despite full-chronology plug-in EPS = 0.85.
Most partitions instead mix both clusters and yield a relatively high
correlation, producing a strongly non-Gaussian, multimodal distribution.

This is the clearest Stage A demonstration of the paper's central claim.

**Manuscript interpretation:** under the exchangeable model, EPS describes the
reliability of the full mean.  Outside that model, the familiar formula remains a
summary of \(N\) and mean correlation but does not establish that two
non-overlapping biological samples reproduce the same collective signal.
Split-half replication answers that different question directly.

## 5. Tree-bootstrap results

The same 1,000 explicit 18-out-of-18 draws were evaluated in every scenario.
Median bootstrap/full correlations remain high (about 0.94) in all cases because
the compared chronologies share represented trees.

The 95th percentile of standardized RMS difference is:

| Scenario | sRMS 95th percentile |
|---|---:|
| Exchangeable | 0.481 |
| Unequal loadings | 0.484 |
| Influential/high-loading group | 0.489 |
| Two clusters | 0.541 |

Differences are clearer in the variance-ratio tails.  The 5th–95th percentile
range is:

| Scenario | Variance-ratio 5th–95th |
|---|---:|
| Exchangeable | 1.083–1.233 |
| Two clusters | 1.066–1.296 |
| Unequal loadings | 1.019–1.304 |
| Influential/high-loading group | 0.937–1.420 |

Scenario-specific audit variables confirm the mechanism rather than merely
showing wider distributions:

- two-cluster sRMS correlates strongly with absolute cluster draw imbalance
  (\(r=0.851\));
- variance ratio correlates with the number of high-loading draws in the
  loading-concentrated scenario (\(r=0.935\));
- variance ratio correlates with the mean loading represented by a bootstrap draw
  in the unequal-loading scenario (\(r=0.830\)).

**Manuscript interpretation:** bootstrap/full correlation is not an independent
replication test.  Direct chronology differences, variance ratios, multiplicity
metadata, and composition-linked tails are more informative.

## 6. Leave-one-tree-out results

LOTO/full correlations exceed 0.9949 in every scenario.  Their uniformly high
values follow from sharing 17 of 18 trees and should not be interpreted as evidence
that omission has no effect.

The exchangeable and symmetric two-cluster scenarios give identical LOTO results
for every tree.  Their extremal-unit summaries are therefore reported as complete
tie sets rather than allowing floating-point roundoff to designate an arbitrary
tree.  In the loading-concentrated scenario, the 15 lower-loading trees are tied
for the maximum direct standardized RMS effect, whereas the three high-loading
trees are tied for the largest absolute changes in EPS and effective dimension.
Unequal loadings produce graded changes:

\[
\Delta\bar r=-0.0101\ \text{to}\ +0.0104,
\]

\[
\Delta\mathrm{EPS}=-0.0151\ \text{to}\ -0.0001,
\]

and

\[
\Delta D_{\mathrm{eff}}=-0.531\ \text{to}\ +0.169.
\]

In the loading-concentrated scenario, omitting a high-loading tree changes EPS
most strongly (\(\Delta\mathrm{EPS}=-0.0275\)) and increases effective dimension
by about 0.993.  Omitting a lower-loading, residual-rich tree changes the
chronology trajectory more strongly in standardized RMS terms (0.0983 versus
0.0640).

**Manuscript interpretation:** there is no single self-evident meaning of
"influential tree."  Influence on the chronology trajectory, average covariance,
EPS, and covariance dimension are different estimands and can identify different
units.  This supports reporting an influence profile rather than a scalar rank.

## 7. Main-text result proposed

The strongest concise result for the main text is:

> Holding \(N=18\), \(\bar r=0.2394\), and the conventional plug-in EPS =
> 0.85 constant did not fix
> chronology composition stability.  Exchangeable split-half correlations were
> identically 0.739, whereas the two-cluster system ranged from 0.014 to 0.803
> across all unique balanced partitions.  Effective dimension simultaneously
> declined from 9.12 to 6.00.  Unequal- and concentrated-loading scenarios showed
> smaller but systematic subset, bootstrap, and LOTO heterogeneity.  Thus a
> single EPS value can be compatible with materially different covariance
> architecture and biological-sample reproducibility.

## 8. Figure strategy emerging from Stage A

### Main Figure 1

EPS surface and the common Stage A point:

\[
(N,\bar r,\mathrm{EPS})=(18,0.2394,0.85).
\]

### Main Figure 2

Four fixed correlation matrices, all plotted on the same scale, accompanied by
effective dimension.

### Main Figure 3

Composition distributions:

- subset-to-full 5–95% envelopes across subset size;
- balanced split-half distributions, emphasizing the two-cluster lower tail;
- bootstrap standardized RMS or variance-ratio tails;
- LOTO \(\Delta\mathrm{EPS}\) versus direct standardized RMS.

The matrices and diagnostics should not be collapsed into a composite score.

## 9. Interpretation limits and reviewer risks

1. **This is an exact covariance experiment, not sampling uncertainty.**  
   That is intentional for Stage A, but the manuscript must state that finite-
   sample estimation uncertainty is deferred to later simulation stages.

2. **The two-cluster case is deliberately strong.**  
   It should be described as a diagnostic stress case, not as a claim that all
   real chronologies have such a structure.

3. **Subset-to-full correlation shares trees with the full chronology.**  
   The manuscript must not present it as independent replication.

4. **Bootstrap results depend on the observed tree pool.**  
   They quantify composition uncertainty conditional on that pool, not population
   generalization beyond it.

5. **LOTO influence is multidimensional.**  
   The apparently counterintuitive high-/low-loading contrast is a result to
   explain, not a defect to hide.

6. **No climate-validity claim follows from Stage A.**  
   The results establish structural and compositional non-equivalence at fixed
   EPS; climate-growth and reconstruction validity remain downstream questions.

## 10. Decision

Stage A meets the go criterion.  The same-plug-in-EPS/different-structure result is
unambiguous and large enough to support the Original Article framing.  The next
scientific step is to convert these frozen outputs into manuscript figures and
then apply the same diagnostics to the positive-control and heterogeneous
empirical datasets.
