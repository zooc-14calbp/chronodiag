test_that("wide CSV is converted to canonical long format", {
  x <- data.frame(year = 1900:1901, A = c(1, 2), B = c(3, NA), check.names = FALSE)
  out <- wide_to_long(x)
  expect_equal(names(out)[1:4], c("year", "series_id", "tree_id", "value"))
  expect_equal(nrow(out), 3)
  expect_true(all(out$series_id == out$tree_id))
})

test_that("long CSV duplicate observations are rejected", {
  f <- tempfile(fileext = ".csv")
  write.csv(data.frame(year = c(1900, 1900), series_id = c("A", "A"), value = c(1, 2)), f, row.names = FALSE)
  expect_error(read_long_csv(f), "duplicate")
})

test_that("existing tree ID has priority over parser fallback", {
  f <- tempfile(fileext = ".csv")
  write.csv(data.frame(year = 1900:1901, series_id = c("A1", "A1"), tree_id = c("TREE", "TREE"), value = c(1, 2)), f, row.names = FALSE)
  parser <- function(ids) data.frame(series_id = ids, tree_id = "PARSED", site_id = "TEST")
  out <- read_long_csv(f, key_parser = parser)
  expect_true(all(out$tree_id == "TREE"))
  expect_true(all(out$site_id == "TEST"))
})

test_that("Tucson blank missing cells remain distinct from zero-width rings", {
  skip_if_not_installed("dplR")
  f <- tempfile(fileext = ".rwl")
  writeLines(c(
    "TEST    1 CHRONODIAG TUCSON FIXTURE",
    "TEST    2 TEST REGION  TEST SPECIES                       1900 1905",
    "TEST    3 GENERATED FIXTURE",
    "ZERO01  1900   100     0   200   999",
    "MISS01  1900   100         200   999"
  ), f, useBytes = TRUE)

  out <- read_tucson_rwl(f, drop_missing = FALSE)
  expect_equal(out$value[out$series_id == "ZERO01" & out$year == 1901], 0)
  expect_true(is.na(out$value[out$series_id == "MISS01" & out$year == 1901]))
  expect_equal(attr(out, "explicit_missing_cells"), 1L)
})
