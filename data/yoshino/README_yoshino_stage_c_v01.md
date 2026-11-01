# Yoshino empirical Stage C data contract

## Scientific input

The scientific source for the frozen input is `nr_all.rwl`, included at the repository-relative path:

```text
data/yoshino/raw/nr_all.rwl
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


## Dual-hash provenance and Git line endings

The original source-byte size and SHA-256 are retained in
`yoshino_raw_manifest_v01.csv`. The repository `.gitattributes` declares
`*.rwl text eol=lf`, so Git stores and checks out a canonical LF representation.
The same manifest therefore also records the size and SHA-256 after CRLF or
lone-CR normalization to LF.

Validation requires the LF-normalized SHA-256 to match on Windows and Linux.
The original source-byte SHA-256 remains provenance metadata and is not
incorrectly compared with a Git-normalized checkout.

Public-release provenance metadata replaces the original workstation path with
the distributed relative path. Original source-byte and LF-normalized hashes
remain unchanged. The analysis inputs and results are unchanged.
