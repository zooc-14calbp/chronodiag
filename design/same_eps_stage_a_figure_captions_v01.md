# Same-EPS Stage A figure captions v01

## Figure 1. EPS surface and the controlled Stage A design

**Expressed population signal as a joint function of biological replication and
mean interseries correlation.** Colors show the conventional plug-in expressed
population signal,

\[
\widehat{\mathrm{EPS}}
=
\frac{N\bar r}{1+(N-1)\bar r},
\]

for \(N=1\)--60 and \(\bar r=0\)--0.60. Contours identify selected coefficient
values; the heavy contour marks 0.85. The star shows the common Stage A design
(\(N=18\), \(\bar r=0.2394\), \(\widehat{\mathrm{EPS}}=0.85\)). The same
coefficient can arise from many combinations of replication and average
correlation and does not, by itself, identify covariance architecture.

Source prototype: `fig01_eps_surface_v01.svg`.

## Figure 2. Same conventional EPS, different covariance architecture

**Four standardized correlation structures with identical \(N\), mean
interseries correlation, and conventional plug-in EPS.** Diagonal elements are
masked in gray; all panels use the same correlation scale. (a) Exchangeable
structure, in which every off-diagonal correlation is 0.2394. (b) Unequal
one-factor loadings. (c) A high-loading minority comprising three of 18 units.
(d) Two equal clusters with within-cluster correlation 0.50 and between-cluster
correlation 0.0078. The corresponding effective dimensions are 9.12, 8.71, 7.85,
and 6.00, respectively. Identical \(N\), \(\bar r\), and plug-in EPS therefore
do not imply the same covariance dimension or modal structure.

Source prototypes:

- `fig02_matrix_exchangeable_v01.svg`;
- `fig02_matrix_unequal_loadings_v01.svg`;
- `fig02_matrix_influential_trees_v01.svg`;
- `fig02_matrix_two_clusters_v01.svg`.

## Figure 3. Composition-aware diagnostics at fixed conventional EPS

**Composition stability for the four Stage A covariance structures.**  
(a) Median subset-to-full chronology correlation across every non-empty proper
subset; shaded bands show the 5th--95th percentiles. The dashed line is the
exchangeable \(\sqrt{\mathrm{SSS}_{n\mid18}}\) reference. Subset/full
comparisons share units and are not independent replication tests.  
(b) Correlations for all 24,310 unique balanced nine-versus-nine splits. Thin
segments show the full range, thick segments the 5th--95th percentiles, and
points the median. The dashed line is the exchangeable half-chronology value
\(\mathrm{EPS}_9=0.739\). The cluster-aligned split in the two-cluster scenario
has correlation 0.014 despite full-chronology plug-in EPS = 0.85.  
(c) Bootstrap/full chronology variance ratios from the same 1,000 explicit
18-out-of-18 tree draws in every scenario. Thin and thick segments show the
1st--99th and 5th--95th percentiles, respectively; points show the median.
Bootstrap/full correlations are not interpreted as independent replication
because represented units are shared.  
(d) Leave-one-tree-out change in plug-in EPS plotted against the standardized
RMS chronology difference. LOTO/full correlations exceed 0.9949 in every
scenario because 17 of 18 units are shared. In the high-loading-minority
scenario, the unit with the largest effect on EPS is not the unit with the
largest trajectory difference, showing that influence is estimand-specific.

Source prototypes:

- `fig03_subset_profiles_v01.svg`;
- `fig04_split_half_intervals_v01.svg`;
- `fig05_bootstrap_variance_ratio_v01.svg`;
- `fig06_loto_influence_v01.svg`.

## Production note

The source prototypes are separate, reproducibly generated figures. Final panel
assembly, panel lettering, font-size harmonization, and journal dimensions will
be applied in the manuscript figure pipeline after the empirical cases determine
the final allocation of main-text and supplementary panels.
