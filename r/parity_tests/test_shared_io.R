source(file.path("r", "chronodiag", "R", "io.R"))

if (!requireNamespace("dplR", quietly = TRUE)) {
  stop("dplR is required for the shared RWL parity test")
}

fixture <- file.path("shared_fixtures", "example_tucson_v02.rwl")
out <- read_tucson_rwl(fixture, drop_missing = FALSE)

stopifnot(setequal(unique(out$series_id), c("ZERO01", "MISS01", "Q  525", "MICRO1")))
stopifnot(out$value[out$series_id == "ZERO01" & out$year == 1900] == 1.0)
stopifnot(out$value[out$series_id == "ZERO01" & out$year == 1901] == 0.0)
stopifnot(out$value[out$series_id == "ZERO01" & out$year == 1902] == 2.0)
stopifnot(is.na(out$value[out$series_id == "MISS01" & out$year == 1901]))
stopifnot(abs(out$value[out$series_id == "MICRO1" & out$year == 1900] - 1.234) < 1e-12)

cat("R shared-I/O parity fixture passed.\n")
