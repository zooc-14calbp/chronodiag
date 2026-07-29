source(file.path("r", "chronodiag", "R", "core_metrics.R"))

fixture_dir <- "shared_fixtures"
expected <- read.csv(file.path(fixture_dir, "expected_summaries.csv"),
                     stringsAsFactors = FALSE)

mean_off_diagonal <- function(x) {
  n <- nrow(x)
  (sum(x) - sum(diag(x))) / (n * (n - 1))
}

for (i in seq_len(nrow(expected))) {
  scenario <- expected$scenario[i]
  matrix <- as.matrix(read.csv(
    file.path(fixture_dir, paste0(scenario, "_correlation.csv")),
    header = FALSE,
    check.names = FALSE
  ))
  storage.mode(matrix) <- "double"
  rho <- mean_off_diagonal(matrix)

  stopifnot(nrow(matrix) == expected$n_trees[i])
  stopifnot(max(abs(matrix - t(matrix))) < 1e-12)
  stopifnot(max(abs(diag(matrix) - 1)) < 1e-12)
  stopifnot(abs(rho - expected$mean_interseries_correlation[i]) < 1e-12)
  stopifnot(abs(eps_value(rho, nrow(matrix)) - expected$eps[i]) < 1e-12)
  stopifnot(abs(effective_dimension(matrix) - expected$effective_dimension[i]) < 1e-10)

  ev <- eigen(matrix, symmetric = TRUE, only.values = TRUE)$values
  leading_fraction <- max(ev) / sum(ev)
  stopifnot(abs(leading_fraction - expected$leading_eigenvalue_fraction[i]) < 1e-10)
  stopifnot(abs(min(ev) - expected$minimum_eigenvalue[i]) < 1e-10)
}

cat("R shared-metric parity fixtures passed.\n")
