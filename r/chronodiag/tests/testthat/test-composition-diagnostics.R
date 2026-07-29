.make_composition_test_data <- function() {
  years <- 2000:2007
  units <- sprintf("T%02d", 1:4)
  signal <- c(-1.0, -0.2, 0.8, 1.1, 0.4, -0.7, -1.2, 0.3)
  rows <- lapply(seq_along(units), function(i) {
    data.frame(
      year = years,
      tree_id = units[i],
      value = signal * c(1.0, 0.9, 1.1, 0.8)[i] + 0.03 * i * seq_along(years),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

test_that("subset stability reports shared units and squared correlation", {
  data <- .make_composition_test_data()
  subsets <- data.frame(
    replicate_id = c("S", "S"),
    unit_id = c("T01", "T02"),
    stringsAsFactors = FALSE
  )
  out <- subset_stability(data, subsets, min_overlap = 5)
  expect_equal(out$diagnostic, "subset_to_full")
  expect_equal(out$group_a_size, 2)
  expect_equal(out$group_b_size, 4)
  expect_equal(out$shared_unit_count, 2)
  expect_equal(out$squared_correlation, out$correlation^2, tolerance = 1e-14)
  expect_equal(out$status, "ok")
})

test_that("split-half stability enforces disjoint exhaustive groups", {
  data <- .make_composition_test_data()
  partitions <- data.frame(
    replicate_id = rep("H", 4),
    group = c("A", "A", "B", "B"),
    unit_id = c("T01", "T03", "T02", "T04"),
    stringsAsFactors = FALSE
  )
  out <- split_half_stability(data, partitions, min_overlap = 5)
  expect_equal(out$diagnostic, "disjoint_split_half")
  expect_equal(out$shared_unit_count, 0)
  expect_equal(out$group_a_size, 2)
  expect_equal(out$group_b_size, 2)

  overlapping <- rbind(
    partitions,
    data.frame(replicate_id = "H", group = "B", unit_id = "T01")
  )
  expect_error(
    split_half_stability(data, overlapping, require_exhaustive = FALSE),
    "not disjoint"
  )

  partial <- data.frame(
    replicate_id = c("P", "P"),
    group = c("A", "B"),
    unit_id = c("T01", "T02")
  )
  expect_error(split_half_stability(data, partial), "not exhaustive")
})

test_that("minimum depth fraction is explicit", {
  data <- .make_composition_test_data()
  data$value[data$tree_id == "T01" & data$year == 2003] <- NA_real_
  subsets <- list(S = c("T01", "T02"))
  complete <- subset_stability(data, subsets, min_depth_fraction = 1)
  permissive <- subset_stability(data, subsets, min_depth_fraction = 0.5)
  expect_equal(complete$n_years, 7)
  expect_equal(permissive$n_years, 8)
})

test_that("duplicate year-unit observations are rejected", {
  data <- .make_composition_test_data()
  data <- rbind(data, data[1, ])
  expect_error(
    subset_stability(data, list(S = "T01")),
    "one value per year"
  )
})

test_that("literal split halves are balanced by default", {
  data <- .make_composition_test_data()
  partitions <- data.frame(
    replicate_id = rep("U", 4),
    group = c("A", "B", "B", "B"),
    unit_id = c("T01", "T02", "T03", "T04"),
    stringsAsFactors = FALSE
  )
  expect_error(split_half_stability(data, partitions), "not balanced")
  out <- split_half_stability(data, partitions, require_balanced = FALSE)
  expect_equal(out$group_a_size, 1)
  expect_equal(out$group_b_size, 3)
})
