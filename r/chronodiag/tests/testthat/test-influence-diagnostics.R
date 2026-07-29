.make_influence_test_data <- function() {
  years <- 2000:2023
  t <- seq_along(years) - 1
  base <- sin(2 * pi * t / 7) + 0.03 * t
  values <- list(
    T01 = base + 0.10 * cos(2 * pi * t / 5),
    T02 = 0.90 * base - 0.08 * sin(2 * pi * t / 4),
    T03 = 1.05 * base + 0.06 * cos(2 * pi * t / 3),
    T04 = base + 0.12 * sin(2 * pi * t / 6),
    T05 = 0.85 * base + 0.15 * cos(2 * pi * t / 8),
    T06 = 1.25 * base + 0.45 * sin(2 * pi * t / 4) + 0.10
  )
  values$T03[6] <- NA_real_
  values$T05[13] <- NA_real_
  values$T06[19] <- NA_real_
  do.call(rbind, lapply(names(values), function(unit) {
    data.frame(
      year = years,
      tree_id = unit,
      value = values[[unit]],
      stringsAsFactors = FALSE
    )
  }))
}

.make_bootstrap_multiplicities <- function() {
  data.frame(
    replicate_id = c(
      rep("B01", 5), rep("B02", 6), rep("B03", 4), rep("B04", 3)
    ),
    unit_id = c(
      "T01", "T02", "T03", "T04", "T05",
      "T01", "T02", "T03", "T04", "T05", "T06",
      "T01", "T03", "T05", "T06",
      "T02", "T04", "T05"
    ),
    multiplicity = c(
      2, 1, 1, 1, 1,
      1, 1, 1, 1, 1, 1,
      1, 1, 1, 3,
      2, 2, 2
    ),
    stringsAsFactors = FALSE
  )
}

test_that("identity bootstrap has zero chronology difference", {
  out <- tree_bootstrap(
    .make_influence_test_data(), .make_bootstrap_multiplicities()
  )
  identity <- out[out$replicate_id == "B02", ]
  expect_true(identity$all_units_represented)
  expect_equal(identity$unique_unit_count, 6)
  expect_equal(identity$kish_effective_unit_count, 6, tolerance = 1e-14)
  expect_equal(identity$correlation, 1, tolerance = 1e-14)
  expect_equal(identity$rms_difference, 0, tolerance = 1e-14)
})

test_that("named bootstrap lists represent literal draws", {
  out <- tree_bootstrap(
    .make_influence_test_data(),
    list(M01 = c("T01", "T01", "T02", "T03", "T04", "T05"))
  )
  expect_equal(out$bootstrap_ids, "T01*2|T02*1|T03*1|T04*1|T05*1")
  expect_equal(out$nominal_draw_count, 6)
  expect_equal(out$unique_unit_count, 5)
})

test_that("bootstrap multiplicity and size contracts are explicit", {
  duplicate <- data.frame(
    replicate_id = c("B", "B"),
    unit_id = c("T01", "T01"),
    multiplicity = c(1, 1)
  )
  expect_error(
    tree_bootstrap(.make_influence_test_data(), duplicate, require_full_size = FALSE),
    "combine"
  )
  zero <- data.frame(replicate_id = "B", unit_id = "T01", multiplicity = 0)
  expect_error(
    tree_bootstrap(.make_influence_test_data(), zero, require_full_size = FALSE),
    "positive"
  )
  short <- data.frame(
    replicate_id = c("B", "B"),
    unit_id = c("T01", "T02"),
    multiplicity = c(1, 1)
  )
  expect_error(tree_bootstrap(.make_influence_test_data(), short), "expected 6")
  out <- tree_bootstrap(
    .make_influence_test_data(), short, require_full_size = FALSE
  )
  expect_equal(out$nominal_draw_count, 2)
})

test_that("Kish count is bounded by unique and nominal counts", {
  out <- tree_bootstrap(
    .make_influence_test_data(), .make_bootstrap_multiplicities()
  )
  expect_true(all(out$kish_effective_unit_count >= 1))
  expect_true(all(out$kish_effective_unit_count <= out$unique_unit_count))
  expect_true(all(out$unique_unit_count <= out$nominal_draw_count))
})

test_that("weighted median expands literal bootstrap draws", {
  data <- data.frame(
    year = rep(2000:2002, 3),
    tree_id = rep(c("T01", "T02", "T03"), each = 3),
    value = c(0, 1, 2, 10, 11, 12, 20, 21, 22),
    stringsAsFactors = FALSE
  )
  out <- tree_bootstrap(
    data,
    list(M = c("T01", "T01", "T02")),
    chronology_method = "median"
  )
  expect_equal(out$mean_difference, -10, tolerance = 1e-14)
  expect_equal(out$rms_difference, 10, tolerance = 1e-14)
  expect_equal(out$correlation, 1, tolerance = 1e-14)
})

test_that("LOTO emphasizes direct differences despite shared-data correlation", {
  out <- leave_one_tree_out(.make_influence_test_data())
  influential <- out[which.max(out$rms_difference), ]
  expect_equal(influential$omitted_unit_id, "T06")
  expect_gt(influential$correlation, 0.99)
  expect_equal(influential$shared_unit_count, 5)
  expect_gt(influential$standardized_rms_difference, 0)
})

test_that("LOTO direct difference follows the complete mean identity", {
  data <- .make_influence_test_data()
  wide <- reshape(
    data,
    idvar = "year",
    timevar = "tree_id",
    direction = "wide"
  )
  value_columns <- grep("^value\\.", names(wide), value = TRUE)
  complete <- stats::complete.cases(wide[value_columns])
  matrix_data <- as.matrix(wide[complete, value_columns, drop = FALSE])
  full <- rowMeans(matrix_data)
  t06_column <- grep("T06$", colnames(matrix_data), value = TRUE)
  reduced <- rowMeans(matrix_data[, setdiff(colnames(matrix_data), t06_column), drop = FALSE])
  identity_difference <- (reduced - matrix_data[, t06_column]) / ncol(matrix_data)
  observed_difference <- reduced - full
  expect_equal(observed_difference, as.numeric(identity_difference), tolerance = 1e-14)

  out <- leave_one_tree_out(data)
  row <- out[out$omitted_unit_id == "T06", ]
  expect_equal(row$mean_difference, mean(observed_difference), tolerance = 1e-12)
  expect_equal(row$rms_difference, sqrt(mean(observed_difference^2)), tolerance = 1e-12)
})

test_that("LOTO structure deltas use one complete period", {
  out <- leave_one_tree_out(.make_influence_test_data())
  expect_true(all(out$structure_status == "ok"))
  expect_true(all(out$structure_n_years == 21))
  expect_equal(length(unique(out$full_rbar)), 1)
  expect_equal(length(unique(out$full_eps)), 1)
  expect_equal(length(unique(out$full_effective_dimension)), 1)
})

test_that("LOTO requires at least two biological units", {
  one <- subset(.make_influence_test_data(), tree_id == "T01")
  expect_error(leave_one_tree_out(one), "at least two")
})
