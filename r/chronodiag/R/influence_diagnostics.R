.influence_normalize_bootstrap <- function(multiplicities) {
  if (is.list(multiplicities) && !is.data.frame(multiplicities)) {
    if (!length(multiplicities)) {
      stop("multiplicities must contain at least one replicate")
    }
    if (is.null(names(multiplicities)) || any(nchar(trimws(names(multiplicities))) == 0L)) {
      stop("a bootstrap list must be named by replicate_id")
    }
    if (anyDuplicated(names(multiplicities))) {
      stop("a bootstrap list contains duplicate replicate_id names")
    }
    rows <- lapply(seq_along(multiplicities), function(i) {
      draws <- multiplicities[[i]]
      if (!length(draws)) {
        stop("bootstrap replicate ", sQuote(names(multiplicities)[i]), " contains no draws")
      }
      draws <- .composition_validate_identifiers(draws, "multiplicities")
      counts <- table(draws)
      data.frame(
        replicate_id = rep(names(multiplicities)[i], length(counts)),
        unit_id = names(counts),
        multiplicity = as.integer(counts),
        stringsAsFactors = FALSE
      )
    })
    out <- do.call(rbind, rows)
  } else if (is.data.frame(multiplicities)) {
    required <- c("replicate_id", "unit_id", "multiplicity")
    missing <- setdiff(required, names(multiplicities))
    if (length(missing)) {
      stop("multiplicities is missing required columns: ", paste(missing, collapse = ", "))
    }
    out <- multiplicities[required]
    out$replicate_id <- .composition_validate_identifiers(
      out$replicate_id, "multiplicities"
    )
    out$unit_id <- .composition_validate_identifiers(out$unit_id, "multiplicities")
    raw_multiplicity <- out$multiplicity
    if (is.factor(raw_multiplicity)) raw_multiplicity <- as.character(raw_multiplicity)
    numeric_multiplicity <- suppressWarnings(as.numeric(raw_multiplicity))
    bad <- is.na(numeric_multiplicity) | !is.finite(numeric_multiplicity) |
      numeric_multiplicity < 1 | numeric_multiplicity != floor(numeric_multiplicity)
    if (any(bad)) stop("multiplicity must contain positive integers")
    out$multiplicity <- as.integer(numeric_multiplicity)
  } else {
    stop("multiplicities must be a data.frame or named list")
  }

  if (!nrow(out)) stop("multiplicities must contain at least one replicate")
  if (anyDuplicated(paste(out$replicate_id, out$unit_id, sep = "\r"))) {
    stop(
      "multiplicities contains duplicate replicate/unit rows; combine them ",
      "into one positive multiplicity"
    )
  }
  out <- out[order(out$replicate_id, out$unit_id), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.influence_validate_bootstrap_units <- function(block, available, replicate_id) {
  units <- as.character(block$unit_id)
  missing <- setdiff(units, available)
  if (length(missing)) {
    stop(
      "bootstrap replicate ", sQuote(replicate_id),
      " refers to units absent from data: ", paste(sort(missing), collapse = ", ")
    )
  }
  order_index <- order(units)
  list(
    units = units[order_index],
    weights = as.integer(block$multiplicity[order_index])
  )
}

.influence_build_weighted_chronology <- function(
    matrix_data, units, weights, method, min_depth_fraction) {
  selected <- matrix_data[, units, drop = FALSE]
  observed <- !is.na(selected)
  weighted_depth <- as.integer(observed %*% as.integer(weights))
  nominal_draw_count <- sum(weights)
  required_depth <- ceiling(min_depth_fraction * nominal_draw_count)

  chronology <- rep(NA_real_, nrow(selected))
  if (method == "mean") {
    weighted_values <- sweep(selected, 2L, weights, `*`)
    numerator <- rowSums(weighted_values, na.rm = TRUE)
    valid_depth <- weighted_depth > 0L
    chronology[valid_depth] <- numerator[valid_depth] / weighted_depth[valid_depth]
  } else {
    for (i in seq_len(nrow(selected))) {
      keep <- !is.na(selected[i, ])
      if (any(keep)) {
        expanded <- rep(as.numeric(selected[i, keep]), times = weights[keep])
        chronology[i] <- stats::median(expanded)
      }
    }
  }
  chronology[weighted_depth < required_depth] <- NA_real_
  list(chronology = chronology, depth = weighted_depth)
}

.influence_compare <- function(
    chronology_a,
    depth_a,
    chronology_b,
    depth_b,
    min_overlap,
    depth_a_label,
    depth_b_label) {
  valid <- !is.na(chronology_a) & !is.na(chronology_b)
  n_years <- sum(valid)
  years <- as.integer(names(chronology_a))[valid]

  out <- list(
    correlation = NA_real_,
    squared_correlation = NA_real_,
    n_years = as.integer(n_years),
    start_year = if (n_years) min(years) else NA_integer_,
    end_year = if (n_years) max(years) else NA_integer_,
    mean_difference = NA_real_,
    rms_difference = NA_real_,
    standardized_rms_difference = NA_real_,
    mean_absolute_difference = NA_real_,
    max_absolute_difference = NA_real_,
    max_abs_year = NA_integer_,
    signed_difference_at_max = NA_real_,
    variance_ratio = NA_real_,
    status = "insufficient_overlap"
  )
  out[[paste0("mean_depth_", depth_a_label)]] <- if (n_years) mean(depth_a[valid]) else NA_real_
  out[[paste0("min_depth_", depth_a_label)]] <- if (n_years) min(depth_a[valid]) else 0L
  out[[paste0("mean_depth_", depth_b_label)]] <- if (n_years) mean(depth_b[valid]) else NA_real_
  out[[paste0("min_depth_", depth_b_label)]] <- if (n_years) min(depth_b[valid]) else 0L

  if (n_years < min_overlap) return(out)

  a <- chronology_a[valid]
  b <- chronology_b[valid]
  difference <- a - b
  absolute <- abs(difference)
  maximum_index <- which.max(absolute)
  variance_a <- stats::var(a)
  variance_b <- stats::var(b)
  rms <- sqrt(mean(difference^2))

  out$mean_difference <- mean(difference)
  out$rms_difference <- rms
  out$standardized_rms_difference <- if (is.finite(variance_b) && variance_b > 0) {
    rms / sqrt(variance_b)
  } else {
    NA_real_
  }
  out$mean_absolute_difference <- mean(absolute)
  out$max_absolute_difference <- absolute[maximum_index]
  out$max_abs_year <- years[maximum_index]
  out$signed_difference_at_max <- difference[maximum_index]
  out$variance_ratio <- if (is.finite(variance_b) && variance_b > 0) {
    variance_a / variance_b
  } else {
    NA_real_
  }

  if (diff(range(a)) == 0 || diff(range(b)) == 0) {
    out$status <- "constant_chronology"
    return(out)
  }
  correlation <- .composition_bound_correlation(stats::cor(a, b, method = "pearson"))
  out$correlation <- correlation
  out$squared_correlation <- correlation^2
  out$status <- "ok"
  out
}

.influence_mean_off_diagonal <- function(correlation) {
  n <- nrow(correlation)
  if (n < 2L) stop("at least two units are required")
  (sum(correlation) - sum(diag(correlation))) / (n * (n - 1))
}

.influence_loto_structure <- function(matrix_data, omitted_unit, min_overlap) {
  complete <- rowSums(is.na(matrix_data)) == 0L
  n_years <- sum(complete)
  empty <- list(
    structure_n_years = as.integer(n_years),
    full_rbar = NA_real_,
    reduced_rbar = NA_real_,
    delta_rbar = NA_real_,
    full_eps = NA_real_,
    reduced_eps = NA_real_,
    delta_eps = NA_real_,
    full_effective_dimension = NA_real_,
    reduced_effective_dimension = NA_real_,
    delta_effective_dimension = NA_real_,
    structure_status = "insufficient_complete_period"
  )
  if (n_years < min_overlap || ncol(matrix_data) < 3L) return(empty)

  full_values <- matrix_data[complete, , drop = FALSE]
  if (any(apply(full_values, 2L, function(x) diff(range(x)) == 0))) {
    empty$structure_status <- "constant_series"
    return(empty)
  }

  reduced_columns <- setdiff(colnames(matrix_data), omitted_unit)
  calculated <- tryCatch({
    full_correlation <- stats::cor(full_values, use = "everything")
    reduced_correlation <- stats::cor(
      full_values[, reduced_columns, drop = FALSE], use = "everything"
    )
    full_rbar <- .influence_mean_off_diagonal(full_correlation)
    reduced_rbar <- .influence_mean_off_diagonal(reduced_correlation)
    full_eps <- eps_value(full_rbar, ncol(matrix_data))
    reduced_eps <- eps_value(reduced_rbar, length(reduced_columns))
    full_dimension <- effective_dimension(full_correlation)
    reduced_dimension <- effective_dimension(reduced_correlation)
    list(
      structure_n_years = as.integer(n_years),
      full_rbar = full_rbar,
      reduced_rbar = reduced_rbar,
      delta_rbar = reduced_rbar - full_rbar,
      full_eps = full_eps,
      reduced_eps = reduced_eps,
      delta_eps = reduced_eps - full_eps,
      full_effective_dimension = full_dimension,
      reduced_effective_dimension = reduced_dimension,
      delta_effective_dimension = reduced_dimension - full_dimension,
      structure_status = "ok"
    )
  }, error = function(e) NULL)

  if (is.null(calculated)) {
    empty$structure_status <- "invalid_structure"
    return(empty)
  }
  calculated
}

#' Explicit tree-level bootstrap chronology diagnostics
#'
#' Evaluate bootstrap draws supplied as positive integer multiplicities. Random
#' draws are intentionally generated outside this numerical engine so that the
#' same bootstrap compositions can be audited in R and Python.
#'
#' @param data Long-form data containing one value per year and biological unit.
#' @param multiplicities A data frame with `replicate_id`, `unit_id`, and positive
#'   integer `multiplicity`, or a named list containing literal unit draws.
#' @param year_col Name of the calendar-year column.
#' @param unit_col Name of the biological-unit column, normally `tree_id`.
#' @param value_col Name of the numeric value column.
#' @param chronology_method Arithmetic `"mean"` or multiplicity-expanded
#'   `"median"` chronology.
#' @param min_depth_fraction Minimum fraction of nominal draws required in a year.
#' @param min_overlap Minimum number of common chronology years.
#' @param require_full_size Require each replicate to contain `N` nominal draws,
#'   where `N` is the full number of biological units.
#' @return One auditable result row per bootstrap replicate. Correlation with the
#'   full chronology is descriptive because the chronologies share units; direct
#'   difference metrics are the primary influence outputs.
#' @export
tree_bootstrap <- function(
    data,
    multiplicities,
    year_col = "year",
    unit_col = "tree_id",
    value_col = "value",
    chronology_method = "mean",
    min_depth_fraction = 1,
    min_overlap = 3L,
    require_full_size = TRUE) {
  if (length(require_full_size) != 1L || is.na(require_full_size) ||
      !is.logical(require_full_size)) {
    stop("require_full_size must be TRUE or FALSE")
  }
  method <- .composition_validate_method(chronology_method)
  fraction <- .composition_validate_fraction(min_depth_fraction)
  overlap <- .composition_validate_count(min_overlap, "min_overlap")
  matrix_data <- .composition_prepare_matrix(data, year_col, unit_col, value_col)
  table <- .influence_normalize_bootstrap(multiplicities)
  available <- colnames(matrix_data)
  full_units <- sort(available)
  full <- .composition_build_chronology(
    matrix_data, full_units, method, fraction
  )
  names(full$chronology) <- attr(matrix_data, "years")

  replicate_ids <- sort(unique(table$replicate_id))
  rows <- lapply(replicate_ids, function(replicate_id) {
    block <- table[table$replicate_id == replicate_id, , drop = FALSE]
    membership <- .influence_validate_bootstrap_units(
      block, available, replicate_id
    )
    nominal_draw_count <- sum(membership$weights)
    if (require_full_size && nominal_draw_count != length(full_units)) {
      stop(
        "bootstrap replicate ", sQuote(replicate_id), " has ",
        nominal_draw_count, " draws; expected ", length(full_units)
      )
    }
    bootstrap <- .influence_build_weighted_chronology(
      matrix_data, membership$units, membership$weights, method, fraction
    )
    names(bootstrap$chronology) <- attr(matrix_data, "years")
    metrics <- .influence_compare(
      bootstrap$chronology,
      bootstrap$depth,
      full$chronology,
      full$depth,
      overlap,
      "bootstrap",
      "full"
    )
    unique_count <- length(membership$units)
    weight_sum <- sum(membership$weights)
    kish <- weight_sum^2 / sum(membership$weights^2)

    data.frame(
      diagnostic = "tree_bootstrap_to_full",
      replicate_id = replicate_id,
      chronology_method = method,
      correlation_method = "pearson",
      min_depth_fraction = fraction,
      min_overlap = overlap,
      full_pool_size = length(full_units),
      nominal_draw_count = nominal_draw_count,
      unique_unit_count = unique_count,
      omitted_unit_count = length(full_units) - unique_count,
      duplicated_unit_count = sum(membership$weights > 1L),
      maximum_multiplicity = max(membership$weights),
      kish_effective_unit_count = kish,
      shared_unit_count = unique_count,
      all_units_represented = setequal(membership$units, available),
      correlation = metrics$correlation,
      squared_correlation = metrics$squared_correlation,
      n_years = metrics$n_years,
      start_year = metrics$start_year,
      end_year = metrics$end_year,
      mean_depth_bootstrap = metrics$mean_depth_bootstrap,
      min_depth_bootstrap = metrics$min_depth_bootstrap,
      mean_depth_full = metrics$mean_depth_full,
      min_depth_full = metrics$min_depth_full,
      mean_difference = metrics$mean_difference,
      rms_difference = metrics$rms_difference,
      standardized_rms_difference = metrics$standardized_rms_difference,
      mean_absolute_difference = metrics$mean_absolute_difference,
      max_absolute_difference = metrics$max_absolute_difference,
      max_abs_year = metrics$max_abs_year,
      signed_difference_at_max = metrics$signed_difference_at_max,
      variance_ratio = metrics$variance_ratio,
      bootstrap_ids = paste0(
        membership$units, "*", membership$weights, collapse = "|"
      ),
      status = metrics$status,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Leave-one-tree-out chronology influence diagnostics
#'
#' Omit each biological unit in turn and compare the reduced chronology with the
#' full chronology. Because the two chronologies share `N - 1` units, correlation
#' is descriptive; direct chronology differences and common-period changes in
#' covariance summaries are the primary influence outputs.
#'
#' @param data Long-form data containing one value per year and biological unit.
#' @param year_col Name of the calendar-year column.
#' @param unit_col Name of the biological-unit column, normally `tree_id`.
#' @param value_col Name of the numeric value column.
#' @param chronology_method Arithmetic `"mean"` or `"median"` chronology.
#' @param min_depth_fraction Minimum fraction of assigned units required in a year.
#' @param min_overlap Minimum number of common chronology years.
#' @return One auditable result row per omitted unit, including direct chronology
#'   differences and changes in mean correlation, EPS, and effective dimension.
#' @export
leave_one_tree_out <- function(
    data,
    year_col = "year",
    unit_col = "tree_id",
    value_col = "value",
    chronology_method = "mean",
    min_depth_fraction = 1,
    min_overlap = 3L) {
  method <- .composition_validate_method(chronology_method)
  fraction <- .composition_validate_fraction(min_depth_fraction)
  overlap <- .composition_validate_count(min_overlap, "min_overlap")
  matrix_data <- .composition_prepare_matrix(data, year_col, unit_col, value_col)
  full_units <- sort(colnames(matrix_data))
  if (length(full_units) < 2L) {
    stop("leave-one-tree-out requires at least two biological units")
  }
  full <- .composition_build_chronology(
    matrix_data, full_units, method, fraction
  )
  names(full$chronology) <- attr(matrix_data, "years")

  rows <- lapply(full_units, function(omitted_unit) {
    retained <- setdiff(full_units, omitted_unit)
    reduced <- .composition_build_chronology(
      matrix_data, retained, method, fraction
    )
    names(reduced$chronology) <- attr(matrix_data, "years")
    metrics <- .influence_compare(
      reduced$chronology,
      reduced$depth,
      full$chronology,
      full$depth,
      overlap,
      "reduced",
      "full"
    )
    structure <- .influence_loto_structure(
      matrix_data, omitted_unit, overlap
    )

    data.frame(
      diagnostic = "leave_one_tree_out",
      omitted_unit_id = omitted_unit,
      chronology_method = method,
      correlation_method = "pearson",
      min_depth_fraction = fraction,
      min_overlap = overlap,
      full_pool_size = length(full_units),
      retained_unit_count = length(retained),
      shared_unit_count = length(retained),
      correlation = metrics$correlation,
      squared_correlation = metrics$squared_correlation,
      n_years = metrics$n_years,
      start_year = metrics$start_year,
      end_year = metrics$end_year,
      mean_depth_reduced = metrics$mean_depth_reduced,
      min_depth_reduced = metrics$min_depth_reduced,
      mean_depth_full = metrics$mean_depth_full,
      min_depth_full = metrics$min_depth_full,
      mean_difference = metrics$mean_difference,
      rms_difference = metrics$rms_difference,
      standardized_rms_difference = metrics$standardized_rms_difference,
      mean_absolute_difference = metrics$mean_absolute_difference,
      max_absolute_difference = metrics$max_absolute_difference,
      max_abs_year = metrics$max_abs_year,
      signed_difference_at_max = metrics$signed_difference_at_max,
      variance_ratio = metrics$variance_ratio,
      structure_n_years = structure$structure_n_years,
      full_rbar = structure$full_rbar,
      reduced_rbar = structure$reduced_rbar,
      delta_rbar = structure$delta_rbar,
      full_eps = structure$full_eps,
      reduced_eps = structure$reduced_eps,
      delta_eps = structure$delta_eps,
      full_effective_dimension = structure$full_effective_dimension,
      reduced_effective_dimension = structure$reduced_effective_dimension,
      delta_effective_dimension = structure$delta_effective_dimension,
      structure_status = structure$structure_status,
      status = metrics$status,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
