# Composition-aware chronology diagnostics specification v01

## 1. Purpose

This specification defines the first composition-aware diagnostics in
`chronodiag`:

- `subset_stability()` — explicit subset versus full-chronology comparisons;
- `split_half_stability()` — explicit, disjoint A/B chronology comparisons.

These diagnostics are **assistive evidence**, not replacements for EPS/SSS and not
components of a universal chronology-strength score. Their purpose is to reveal
how chronology behavior depends on the identities of the biological replication
units used to construct it.

## 2. Biological unit and analysis input

The default resampling unit is `tree_id`. A different unit column may be supplied,
but it must identify the biological unit whose composition is being varied.
Individual cores or radii must not be treated as independent trees unless that is
the explicit scientific design.

Input is canonical long-form data with one row per:

```text
year × biological unit
```

and at least:

```text
year, tree_id, value
```

Duplicate year–unit observations are rejected. Multiple cores/radii must therefore
be aggregated upstream by a documented rule, or analyzed later with a specifically
hierarchical diagnostic.

`chronodiag` does not detrend or standardize within these functions. The input
series and the standardization method are part of the analysis specification.

## 3. Explicit partitions as the computational primitive

The numerical engines consume explicit membership tables. They do not generate
random subsets internally in v01.

This design separates two tasks:

1. **partition generation**, where a random or enumerative design is chosen; and
2. **diagnostic evaluation**, where fixed groups are converted to chronologies and
   compared.

The separation makes every realization auditable and permits exact R/Python
parity. Language-native random-number generators need not produce identical
streams because cross-language validation uses the same frozen partition files.

### 3.1 Subset membership table

```text
replicate_id, unit_id
```

Each `replicate_id` defines one subset. Repeated assignment of a unit within a
replicate is an error.

### 3.2 Split-half membership table

```text
replicate_id, group, unit_id
```

`group` is `A` or `B`. The two groups must be disjoint. By default their union must
exhaust the complete unit pool, although deliberately partial disjoint comparisons
may be requested explicitly.

## 4. Chronology construction

Two chronology constructors are supported:

- `mean` — arithmetic mean;
- `median` — median across available units.

The constructor is reported in every result row. The model-defined relationships
to EPS and SSS apply to arithmetic means, not automatically to medians.

### 4.1 Sample-depth rule

For a group containing `g` units, a chronology value is available in year `t` only
when

\[
n_t \ge \lceil f g \rceil,
\]

where `n_t` is the observed group depth and `f = min_depth_fraction`.

The default is

```text
min_depth_fraction = 1
```

so that each comparison is complete with respect to the units assigned to each
group. Lower fractions are allowed only by explicit request. Result rows report the
mean and minimum realized depths of both chronologies over their common years.

## 5. Subset-to-full stability

For subset `A_b` in replicate `b`, let

\[
C_{A_b}(t)
\]

be its chronology and let

\[
C_N(t)
\]

be the chronology constructed from the full biological-unit pool. The diagnostic
reports

\[
r^{\mathrm{full}}_b = \operatorname{cor}\{C_{A_b}, C_N\}
\]

and

\[
\left(r^{\mathrm{full}}_b\right)^2.
\]

Because `A_b` is contained in the full set, the two chronologies share units. The
result therefore describes **subset-to-full correspondence**, not independent
replication.

Under an exchangeable common-signal model, arithmetic-mean chronologies, and a
fixed analysis interval,

\[
\operatorname{cor}^2(C_{A_b}, C_N)
=
\operatorname{SSS}_{|A_b|\mid N}.
\]

Outside those conditions, the squared empirical correlation is not labelled SSS.
The output records `shared_unit_count = |A_b|` so that overlap is explicit.

## 6. Disjoint split-half stability

For non-overlapping groups `A_b` and `B_b`, the diagnostic reports

\[
r^{\mathrm{disjoint}}_b
=
\operatorname{cor}\{C_{A_b}, C_{B_b}\}
\]

and its square. Here,

```text
shared_unit_count = 0
```

by construction.

Under the exchangeable arithmetic-mean model, for group sizes `n` and `m`,

\[
\operatorname{cor}(C_A,C_B)
=
\frac{\rho\sqrt{nm}}
{\sqrt{\{1+(n-1)\rho\}\{1+(m-1)\rho\}}}.
\]

For equal halves (`n=m`), this correlation equals `EPS_n`. It is therefore a more
stringent composition-replication diagnostic than subset-to-full correspondence,
which benefits mechanically from shared units.

## 7. Correlation status

A result row is always retained for every syntactically valid replicate.

- `status = ok` — at least `min_overlap` years and nonconstant chronologies;
- `status = insufficient_overlap` — fewer than `min_overlap` common years;
- `status = constant_chronology` — sufficient overlap but zero variance in at
  least one chronology.

For the latter two states, `correlation` and `squared_correlation` are missing.
Malformed partitions, unknown unit IDs, duplicate observations, overlapping
split halves, and non-exhaustive halves under the default policy are hard errors.

## 8. Result contract

Each row contains:

```text
diagnostic
replicate_id
chronology_method
correlation_method
min_depth_fraction
min_overlap
full_pool_size
group_a_size
group_b_size
shared_unit_count
comparison_exhaustive
correlation
squared_correlation
n_years
start_year
end_year
mean_depth_a
min_depth_a
mean_depth_b
min_depth_b
group_a_ids
group_b_ids
status
```

The result is a realization-level table. Distribution summaries such as medians,
interquartile ranges, 5–95% intervals, and lower-tail quantiles are downstream
operations and must retain the number of valid and invalid realizations.

## 9. Frozen cross-language fixtures

The shared fixture set contains:

```text
composition_example_long_v01.csv
composition_exchangeable_exact_v01.csv
subset_partitions_v01.csv
split_half_partitions_v01.csv
expected_subset_stability_v01.csv
expected_split_half_stability_v01.csv
```

The synthetic example includes explicit missing observations. The frozen expected
outputs use:

```text
chronology_method = mean
min_depth_fraction = 1
min_overlap = 10
```

R and Python must reproduce all numerical columns to absolute tolerance `1e-12`
and all identifier/status columns exactly. The exact exchangeable fixture verifies
that squared subset-to-full correlation reduces to SSS and that equal, disjoint
half-chronology correlation reduces to the EPS of one half.

## 10. Scope boundary and next stage

v01 does not include:

- language-independent random partition generation;
- tree bootstrap;
- leave-one-tree-out influence;
- temporal moving-window diagnostics;
- frequency-domain diagnostics;
- automatic interpretation or pass/fail thresholds.

The next composition stage will reuse this chronology-construction and result
contract for tree bootstrap and leave-one-tree-out influence, while retaining their
separate estimands.
