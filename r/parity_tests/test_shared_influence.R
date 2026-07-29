source(file.path("r", "chronodiag", "R", "core_metrics.R"))
source(file.path("r", "chronodiag", "R", "composition_diagnostics.R"))
source(file.path("r", "chronodiag", "R", "influence_diagnostics.R"))

fixture_dir <- "shared_fixtures"
data <- utils::read.csv(
  file.path(fixture_dir, "influence_example_long_v01.csv"),
  stringsAsFactors = FALSE
)
multiplicities <- utils::read.csv(
  file.path(fixture_dir, "bootstrap_multiplicities_v01.csv"),
  stringsAsFactors = FALSE
)

actual_bootstrap <- tree_bootstrap(data, multiplicities)
actual_loto <- leave_one_tree_out(data)
expected_bootstrap <- utils::read.csv(
  file.path(fixture_dir, "expected_tree_bootstrap_v01.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
expected_loto <- utils::read.csv(
  file.path(fixture_dir, "expected_leave_one_tree_out_v01.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

compare_frames <- function(actual, expected, label) {
  stopifnot(identical(names(actual), names(expected)))
  stopifnot(nrow(actual) == nrow(expected))
  logical_columns <- names(expected)[vapply(expected, is.logical, logical(1))]
  numeric_columns <- names(expected)[vapply(expected, is.numeric, logical(1))]
  text_columns <- setdiff(names(expected), c(logical_columns, numeric_columns))

  for (column in numeric_columns) {
    stopifnot(isTRUE(all.equal(
      as.numeric(actual[[column]]),
      as.numeric(expected[[column]]),
      tolerance = 1e-12,
      check.attributes = FALSE
    )))
  }
  for (column in logical_columns) {
    stopifnot(identical(as.logical(actual[[column]]), as.logical(expected[[column]])))
  }
  for (column in text_columns) {
    stopifnot(identical(as.character(actual[[column]]), as.character(expected[[column]])))
  }
  cat(label, "passed.\n")
}

compare_frames(
  actual_bootstrap,
  expected_bootstrap,
  "R tree-bootstrap parity fixture"
)
compare_frames(
  actual_loto,
  expected_loto,
  "R leave-one-tree-out parity fixture"
)

identity <- actual_bootstrap[actual_bootstrap$replicate_id == "B02", ]
stopifnot(isTRUE(all.equal(identity$rms_difference, 0, tolerance = 1e-14)))
stopifnot(identity$unique_unit_count == identity$full_pool_size)
influential <- actual_loto[which.max(actual_loto$rms_difference), ]
stopifnot(influential$omitted_unit_id == "T06")
stopifnot(influential$correlation > 0.99)
cat("R bootstrap/LOTO invariants passed.\n")
