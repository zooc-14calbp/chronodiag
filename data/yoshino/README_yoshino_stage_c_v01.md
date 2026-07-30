# Yoshino empirical Stage C data contract

## Scientific input

The frozen input used here is `nr_all.rwl`, copied byte-for-byte from:

```text
D:\Dropbox\papers_2026\01c_Dendro_2026_yoshino\floatsR0\data_raw\rwl\nr_all.rwl
```

The file contains **60 specimen-mean ring-width series**, one series per
biological specimen. Project records indicate that the original measurement
campaign comprised 72 measured radii from 60 specimens; multiple radii were
averaged upstream in the original Yoshino workflow. Stage C therefore does not
pretend that `nr_all.rwl` contains 72 independent raw radii.

The specimen-mean input is appropriate for composition diagnostics because the
biological resampling unit is the specimen. A future public raw-radii release may
be added after the proprietary measurement archive and radius-to-specimen map
have been curated separately.

## Key convention

The six-character series key is preserved exactly. In the project convention:

- the first three digits define source group 001, 002, or 003;
- `site_id` is recorded as NR01, NR02, or NR03;
- the full six-character key is both `specimen_id` and `tree_id`;
- `core_id` is recorded as `specimen_mean`, not as an original radius.

## Primary complete cohort

The fixed complete cohort is selected from series spans without using ring-width
values, correlations, EPS, or any composition diagnostic. Candidate intervals
must satisfy:

```text
T >= 50 years
T >= 3N
```

Candidates are ranked by maximum specimen count, then maximum interval length,
then earliest start. The resulting primary cohort is:

```text
N = 50 specimens
1850--2004
T = 155 years
```

This is an analysis-design constraint for covariance estimation, not a universal
chronology-quality threshold.
