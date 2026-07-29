# Manuscript v03 detailed outline v01

## Working title

**What does the expressed population signal express? A reliability-theoretic and
composition-aware diagnostic framework for tree-ring chronologies**

## 1. Introduction

1. Replication is fundamental to chronology construction, but internal coherence
   and environmental interpretation are different inferential stages.
2. EPS and SSS answer useful model-defined questions.
3. The practical problem is not that EPS is invalid, but that a single point estimate
   and the 0.85 convention are often asked to stand for multidimensional robustness.
4. Central principle:

   > Before a chronology is interpreted climatically or ecologically, it should be
   > shown to function as a stable collective measurement.

5. Objective: develop a reproducible diagnostic companion, not a successor score.

## 2. EPS and SSS as model-relative reliability coefficients

- original Wigley targets;
- exchangeable common-signal model;
- reliability of the mean chronology;
- correspondence with standardized coefficient alpha;
- SSS as relative reliability;
- 0.85 as an application-dependent uncertainty budget;
- tree versus core as the biological replication unit.

## 3. Why one average correlation is structurally incomplete

- unequal loadings;
- clustered covariance;
- influential trees;
- correlated residual structure;
- temporal change;
- same `N`, same mean correlation, same EPS, different covariance architecture.

## 4. An assistive chronology diagnostic profile

### 4.1 Model-relative reliability

EPS, SSS and resampling uncertainty.

### 4.2 Covariance architecture

Correlation matrix, eigenvalue spectrum and effective dimension.

### 4.3 Composition stability

Subset-to-full, subset-to-complement, repeated split-half, tree bootstrap and
leave-one-tree-out influence.

### 4.4 Sample-depth and temporal stability

Empirical subset-size curves, non-overlapping period comparisons and descriptive
moving-window diagnostics.

### 4.5 Scope boundary

These diagnostics establish robustness of the chronology as an object of inference.
They do not identify the environmental origin of the common signal or validate a
climate reconstruction.

## 5. Simulation design

- Stage A: same expected EPS under four covariance structures;
- Stage B: time-varying structure;
- Stage C: sample size, record length, autocorrelation and missingness factorial;
- outcomes: estimator bias, subset-distribution width, influence concentration,
  covariance dimension, interval coverage and temporal drift;
- latent-signal recovery is used only where a latent signal is defined by simulation.

## 6. Empirical cases

- public ITRDB datasets selected for contrasting replication and covariance
  structures;
- Yoshino as a provenance-constrained, intensively managed plantation test case;
- no climate target is required for dataset inclusion;
- public and private data remain physically separated.

## 7. R and Python reference implementations

- language-neutral scientific specification;
- native implementations with shared fixtures;
- reproducible reports and machine-readable output;
- no integrated score or threshold-based classification.

## 8. Results

- ideal exchangeable benchmark;
- same-EPS/different-structure simulation;
- composition and temporal diagnostics in empirical datasets;
- cross-language numerical agreement.

## 9. Implications for chronology reporting

Proposed minimum reporting set:

1. trees and cores separately;
2. overlap and correlation definition;
3. EPS/SSS formula and implementation;
4. uncertainty interval;
5. sample depth separate from annual SSS;
6. one composition-stability diagnostic;
7. one temporal-stability diagnostic when record length permits;
8. covariance architecture where sufficient overlap exists;
9. no universal 0.85 pass/fail use;
10. climate or ecological attribution assessed downstream and directly.

## 10. Discussion and conclusions

- EPS remains informative for its estimand;
- chronology robustness is an evidence profile, not one number;
- stable coherence is a prerequisite for, not proof of, later dendroclimatic or
  dendroecological interpretation;
- diverse archives, including managed plantations, can be evaluated by data
  behaviour rather than site category alone.
