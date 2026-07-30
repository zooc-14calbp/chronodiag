# ITRDB empirical Stage B contracts

This directory contains public metadata and biological-unit parsing contracts,
not the raw ITRDB RWL files.

## Files

- `empirical_dataset_contract_v01.csv` — DOI, raw hash, explicit Tucson
  reader settings, fixed cohort, and preprocessing design for ca533, GERM012,
  and FRAN003.
- `keymaps/*_keymap_v01.csv` — exact mapping from raw series IDs to biological
  tree and core/radius IDs, including inclusion/exclusion status.

## Source status

- ca533 parsing follows the documented dplR `read.ids(..., stc=c(3,2,1))`
  convention.
- GERM012 parsing is an explicit inference from the complete paired key pattern
  and parallel TRW/LWW identifiers. It is not labelled as source-documented.
- FRAN003 preserves each exact series ID as one biological unit.

## Raw data

Raw files are expected at the project level:

```text
<project-root>/data/itrdb_raw/
```

They are acquired from the DOI-addressable ITRDB records and verified against
the SHA-256 values in the contract. The public repository should distribute an
acquisition/verification workflow rather than silently replacing the source
archive.


## Explicit input format

The empirical workflow never auto-detects these archived files. Each dataset
declares `rwl_format = tucson`, whether a three-line header is present, whether
long-year parsing is needed, and whether zero-width rings are preserved. This
avoids allowing punctuation in a valid Tucson header to determine the parser.
