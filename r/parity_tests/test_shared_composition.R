source(file.path("r", "chronodiag", "R", "core_metrics.R"))
source(file.path("r", "chronodiag", "R", "composition_diagnostics.R"))

fixture_dir <- "shared_fixtures"
data <- utils::read.csv(
  file.path(fixture_dir, "composition_example_long_v01.csv"),
  stringsAsFactors = FALSE
)
subsets <- utils::read.csv(
  file.path(fixture_dir, "subset_partitions_v01.csv"),
  stringsAsFactors = FALSE
)
splits <- utils::read.csv(
  file.path(fixture_dir, "split_half_partitions_v01.csv"),
  stringsAsFactors = FALSE
)

actual_subset <- subset_stability(
  data,
  subsets,
  min_overlap = 10,
  min_depth_fraction = 1
)
actual_split <- split_half_stability(
  data,
  splits,
  min_overlap = 10,
  min_depth_fraction = 1
)
expected_subset <- utils::read.csv(
  file.path(fixture_dir, "expected_subset_stability_v01.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
expected_split <- utils::read.csv(
  file.path(fixture_dir, "expected_split_half_stability_v01.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

compare_frames <- function(actual, expected, label) {
  stopifnot(identical(names(actual), names(expected)))
  stopifnot(nrow(actual) == nrow(expected))
  numeric_columns <- c(
    "min_depth_fraction", "min_overlap", "full_pool_size",
    "group_a_size", "group_b_size",
    "shared_unit_count", "correlation", "squared_correlation", "n_years",
    "start_year", "end_year", "mean_depth_a", "min_depth_a",
    "mean_depth_b", "min_depth_b"
  )
  for (column in numeric_columns) {
    stopifnot(isTRUE(all.equal(
      as.numeric(actual[[column]]),
      as.numeric(expected[[column]]),
      tolerance = 1e-12,
      check.attributes = FALSE
    )))
  }
  text_columns <- setdiff(names(actual), numeric_columns)
  for (column in text_columns) {
    stopifnot(identical(as.character(actual[[column]]), as.character(expected[[column]])))
  }
  cat(label, "passed.\n")
}

compare_frames(actual_subset, expected_subset, "R subset-stability parity fixture")
compare_frames(actual_split, expected_split, "R split-half parity fixture")

exchangeable <- utils::read.csv(
  file.path(fixture_dir, "composition_exchangeable_exact_v01.csv"),
  stringsAsFactors = FALSE
)
exchangeable_subset <- subset_stability(
  exchangeable,
  list(S = c("E01", "E02", "E03")),
  min_overlap = 20
)
stopifnot(isTRUE(all.equal(
  exchangeable_subset$squared_correlation,
  sss_value(0.3, 3, 6),
  tolerance = 1e-12
)))
exchangeable_partitions <- data.frame(
  replicate_id = rep("H", 6),
  group = c("A", "A", "A", "B", "B", "B"),
  unit_id = sprintf("E%02d", 1:6),
  stringsAsFactors = FALSE
)
exchangeable_split <- split_half_stability(
  exchangeable, exchangeable_partitions, min_overlap = 20
)
stopifnot(isTRUE(all.equal(
  exchangeable_split$correlation,
  eps_value(0.3, 3),
  tolerance = 1e-12
)))
cat("R exchangeable EPS/SSS linkage fixture passed.\n")
