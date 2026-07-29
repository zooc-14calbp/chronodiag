.validate_count <- function(x, name) {
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x < 1 || x != floor(x)) {
    stop(name, " must be one positive integer")
  }
  as.integer(x)
}

.validate_rbar <- function(rbar) {
  if (length(rbar) != 1L || is.na(rbar) || !is.finite(rbar)) {
    stop("rbar must be one finite number")
  }
  if (rbar < -1 || rbar > 1) stop("rbar must lie in [-1, 1]")
  as.numeric(rbar)
}

#' Expressed population signal
#'
#' @param rbar Mean interseries correlation.
#' @param n Number of independent biological replication units (trees).
#' @return Numeric scalar. Negative values are not clipped.
#' @export
eps_value <- function(rbar, n) {
  rho <- .validate_rbar(rbar)
  n <- .validate_count(n, "n")
  denom <- 1 + (n - 1) * rho
  if (denom <= 0 || isTRUE(all.equal(denom, 0, tolerance = .Machine$double.eps))) {
    stop("rbar is outside the admissible equicorrelation range for n")
  }
  n * rho / denom
}

#' Subsample signal strength
#'
#' @param rbar Mean interseries correlation.
#' @param n_sub Subsample tree count.
#' @param n_full Full tree count.
#' @return Numeric scalar.
#' @export
sss_value <- function(rbar, n_sub, n_full) {
  rho <- .validate_rbar(rbar)
  n_sub <- .validate_count(n_sub, "n_sub")
  n_full <- .validate_count(n_full, "n_full")
  if (n_sub > n_full) stop("n_sub must not exceed n_full")
  full_term <- 1 + (n_full - 1) * rho
  sub_term <- 1 + (n_sub - 1) * rho
  if (full_term <= 0 || sub_term <= 0) {
    stop("rbar is outside the admissible equicorrelation range")
  }
  n_sub * full_term / (n_full * sub_term)
}

#' Effective covariance dimension
#'
#' @param correlation A finite, symmetric correlation matrix.
#' @param tolerance Numerical tolerance for validation.
#' @return Numeric scalar.
#' @export
effective_dimension <- function(correlation, tolerance = 1e-10) {
  if (length(tolerance) != 1L || is.na(tolerance) || !is.finite(tolerance) || tolerance < 0) {
    stop("tolerance must be finite and non-negative")
  }
  x <- as.matrix(correlation)
  if (length(dim(x)) != 2L || nrow(x) != ncol(x) || nrow(x) < 1L) {
    stop("correlation must be a non-empty square matrix")
  }
  if (any(!is.finite(x))) stop("correlation contains non-finite values")
  if (max(abs(x - t(x))) > tolerance) stop("correlation must be symmetric")
  if (max(abs(diag(x) - 1)) > tolerance) stop("correlation diagonal must equal one")
  ev <- eigen(x, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) < -tolerance) stop("correlation must be positive semidefinite")
  tr <- sum(diag(x))
  denom <- sum(x^2)
  if (denom <= 0) stop("correlation has zero Frobenius norm")
  result <- tr^2 / denom
  n <- nrow(x)
  if (result < 1 - tolerance || result > n + tolerance) {
    stop("effective dimension fell outside [1, n]")
  }
  min(n, max(1, result))
}
