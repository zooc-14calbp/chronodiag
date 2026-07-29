# Manuscript v03 detailed outline v02

## Fixed publication metadata

- **Working title:** *What does the expressed population signal express? A reliability-theoretic and composition-aware framework for tree-ring chronology coherence*
- **Target journal:** *Dendrochronologia*
- **Article type:** Original Article
- **Authors:** Hitoshi Yonenobu; Yasuharu Hoshino; Motonari Ohyama
- **Corresponding author:** Hitoshi Yonenobu

The v02 manuscript remains the frozen theory baseline. Version v03 will be the first manuscript version that combines the three-author front matter, Original Article scope, mathematical Supplement, simulation, empirical datasets, and the R/Python implementation.

## 1. Introduction

1. Replication is fundamental to chronology construction, but internal coherence and environmental interpretation are distinct inferential stages.
2. EPS and SSS answer useful model-defined questions.
3. A single point estimate and the 0.85 convention cannot represent multidimensional robustness.
4. Central principle: before a chronology is interpreted climatically or ecologically, it should be shown to function as a stable collective measurement.
5. Objective: develop a reproducible diagnostic companion rather than a successor score.

## 2. EPS and SSS as model-relative reliability coefficients

- Wigley et al.'s original targets;
- exchangeable common-signal model;
- reliability of the mean chronology;
- correspondence with standardized coefficient alpha;
- SSS as subset-to-full squared correlation and relative reliability;
- direct SSS formula at zero mean correlation;
- 0.85 as an application-dependent uncertainty budget;
- tree versus core as the biological replication unit.

## 3. Why one average correlation is structurally incomplete

- unequal loadings;
- clustered covariance;
- influential trees;
- correlated residual structure;
- temporal change;
- same N, same mean correlation, same EPS, different covariance architecture;
- effective dimension as a complementary structural descriptor, not a replacement score.

## 4. Assistive chronology diagnostic profile

### 4.1 Model-relative reliability
EPS, SSS, RBAR, and tree-level resampling uncertainty.

### 4.2 Covariance architecture
Complete-period correlation matrix, eigenvalue spectrum, leading-mode fraction, and effective dimension.

### 4.3 Composition stability
Subset-to-full, subset-to-complement, repeated split-half, tree bootstrap, and leave-one-tree-out influence based on direct changes rather than shared-data correlation alone.

### 4.4 Sample-depth and temporal stability
Empirical subset-size curves, non-overlapping period comparisons, and descriptive moving-window diagnostics.

### 4.5 Scope boundary
The diagnostics establish robustness of the chronology as an object of inference. They do not identify the environmental origin of common variation or validate a climate reconstruction.

## 5. Mathematical Supplement

- notation and data hierarchy;
- elementary derivation of EPS;
- elementary derivation of SSS;
- overlapping versus disjoint subsets;
- leave-one-tree-out overlap;
- multiple cores and biological units;
- general non-exchangeable reliability model;
- effective dimension and its exchangeable upper bound;
- numerical contract linking equations, functions, and tests.

## 6. Simulation design

- Stage A: same expected EPS under exchangeable, unequal-loading, two-cluster, and influential-tree structures;
- Stage B: time-varying structure;
- Stage C: sample size, record length, autocorrelation, and missingness;
- outcomes: estimator bias, subset-distribution width, influence concentration, covariance dimension, interval coverage, and temporal drift.

## 7. Empirical datasets

- ca533 as the canonical dplR/dplPy benchmark;
- GERM012 as the tree/core hierarchy case;
- FRAN003 as the long complete-case structural example;
- one historical-composition stress case selected after common preprocessing;
- Yoshino as a provenance-constrained, intensively managed plantation case;
- no climate target required for inclusion.

The Yoshino chronology is a central empirical case. The chronology developers are therefore included as coauthors, while the present article remains focused on structural, compositional, sample-depth, and temporal coherence rather than climate-reconstruction performance.

## 8. R and Python implementations

- language-neutral scientific specification;
- native R and native Python implementations;
- shared deterministic fixtures;
- CI-based numerical parity;
- reproducible reports and machine-readable output;
- no integrated score or threshold-based classification.

## 9. Results

- exchangeable benchmark;
- same-EPS/different-structure simulation;
- cross-language agreement;
- empirical composition and temporal diagnostics;
- Yoshino chronology as a non-classical archive evaluated independently of climate attribution.

## 10. Implications for chronology reporting

Minimum reporting set:

1. trees and cores separately;
2. overlap and correlation definition;
3. EPS/SSS formula and implementation;
4. uncertainty interval;
5. sample depth separate from annual SSS;
6. at least one composition-stability diagnostic;
7. at least one temporal-stability diagnostic when record length permits;
8. covariance architecture where sufficient overlap exists;
9. no universal 0.85 pass/fail use;
10. climatic or ecological attribution assessed downstream and directly.

## 11. Discussion and conclusions

- EPS remains informative for its estimand;
- chronology robustness is an evidence profile rather than one number;
- stable coherence is a prerequisite for, not proof of, later dendroclimatic or dendroecological interpretation;
- diverse archives, including managed plantations, can be evaluated from data behavior rather than site category alone.
