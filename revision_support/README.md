# R4 analytical comparisons and point-only figure

Run from the repository root with Python, NumPy and Matplotlib installed:

    python revision_support/reproduce_revision.py --output-dir /path/to/separate-output

Use a directory appropriate to your system. The --output-dir option is optional;
its default is revision_support/reproduced. Add --skip-figure for numerical checks
only. --repo-root allows an explicit repository-root location.

The script verifies the three published-data examples using the unchanged
archived standardized series and 2,000 frozen bootstrap memberships per dataset.
It checks the EPS/SSS population benchmarks, the exact conditional mean bootstrap
variance ratio, and the effective-dimension/pairwise-dispersion identity, and
recreates the point-only comparison figure from unrounded structure summaries.
It generates no new resampling or surrogate series and does not rerun R/dplR
detrending. Correlation benchmarks are not exact finite-record expected medians;
the conditional bootstrap expectation is a mean identity, not a tail test.

expected_results contains the numerical reference records and preparation-record
hashes. The script writes its actual repository-relative input hashes, dependency
versions, numerical comparisons and generated figure into the chosen output
folder. Figures may differ at the PDF byte level because of rendering metadata,
fonts or dependency versions; their data and plotting specifications are fixed.

This directory contains scientific reproduction support only. Submission letters,
reviewer comments and manuscript files are not part of the public source release.
