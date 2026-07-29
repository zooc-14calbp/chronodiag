.composition_validate_count <- function(x, name) {
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x < 1 || x != floor(x)) {
    stop(name, " must be one positive integer")
  }
  as.integer(x)
}

.composition_validate_fraction <- function(x) {
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x <= 0 || x > 1) {
    stop("min_depth_fraction must lie in (0, 1]")
  }
  as.numeric(x)
}

.composition_validate_method <- function(x) {
  if (length(x) != 1L || is.na(x)) {
    stop("chronology_method must be 'mean' or 'median'")
  }
  method <- tolower(as.character(x))
  if (!method %in% c("mean", "median")) {
    stop("chronology_method must be 'mean' or 'median'")
  }
  method
}

.composition_as_numeric <- function(x, name) {
  original_missing <- is.na(x)
  if (is.factor(x)) x <- as.character(x)
  out <- suppressWarnings(as.numeric(x))
  bad <- !original_missing & is.na(out)
  if (any(bad)) stop(name, " contains non-numeric values")
  if (any(!is.finite(out) & !is.na(out))) {
    stop(name, " must contain finite numbers or missing values")
  }
  out
}

.composition_validate_identifiers <- function(x, name) {
  if (anyNA(x)) stop(name, " contains missing identifiers")
  out <- as.character(x)
  if (any(nchar(trimws(out)) == 0L)) stop(name, " contains empty identifiers")
  out
}

.composition_prepare_matrix <- function(data, year_col, unit_col, value_col) {
  if (!is.data.frame(data)) stop("data must be a data.frame")
  required <- c(year_col, unit_col, value_col)
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop("data is missing required columns: ", paste(missing, collapse = ", "))
  }

  units <- .composition_validate_identifiers(data[[unit_col]], unit_col)
  years_num <- .composition_as_numeric(data[[year_col]], year_col)
  if (anyNA(years_num) || any(years_num != floor(years_num))) {
    stop(year_col, " must contain finite integer years")
  }
  years <- as.integer(years_num)
  values <- .composition_as_numeric(data[[value_col]], value_col)

  key <- paste(years, units, sep = "\r")
  if (anyDuplicated(key)) {
    index <- which(duplicated(key) | duplicated(key, fromLast = TRUE))[1L]
    stop(
      "data must contain one value per year and biological unit; duplicate found for ",
      sQuote(units[index]), " in ", years[index]
    )
  }

  year_levels <- sort(unique(years))
  unit_levels <- sort(unique(units))
  matrix_out <- matrix(
    NA_real_,
    nrow = length(year_levels),
    ncol = length(unit_levels),
    dimnames = list(as.character(year_levels), unit_levels)
  )
  index <- cbind(match(years, year_levels), match(units, unit_levels))
  matrix_out[index] <- values
  empty_units <- colnames(matrix_out)[colSums(!is.na(matrix_out)) == 0L]
  if (length(empty_units)) {
    stop(
      "biological units contain no finite observations: ",
      paste(empty_units, collapse = ", ")
    )
  }
  attr(matrix_out, "years") <- year_levels
  matrix_out
}

.composition_normalize_subsets <- function(subsets) {
  if (is.list(subsets) && !is.data.frame(subsets)) {
    if (!length(subsets)) stop("subsets must contain at least one replicate")
    if (is.null(names(subsets)) || any(nchar(trimws(names(subsets))) == 0L)) {
      stop("a subset list must be named by replicate_id")
    }
    if (anyDuplicated(names(subsets))) {
      stop("a subset list contains duplicate replicate_id names")
    }
    if (any(lengths(subsets) == 0L)) {
      empty <- names(subsets)[lengths(subsets) == 0L][1L]
      stop("subset replicate ", sQuote(empty), " contains no biological units")
    }
    rows <- lapply(seq_along(subsets), function(i) {
      data.frame(
        replicate_id = rep(names(subsets)[i], length(subsets[[i]])),
        unit_id = .composition_validate_identifiers(subsets[[i]], "subsets"),
        stringsAsFactors = FALSE
      )
    })
    out <- do.call(rbind, rows)
  } else if (is.data.frame(subsets)) {
    required <- c("replicate_id", "unit_id")
    missing <- setdiff(required, names(subsets))
    if (length(missing)) {
      stop("subsets is missing required columns: ", paste(missing, collapse = ", "))
    }
    out <- subsets[required]
    out$replicate_id <- .composition_validate_identifiers(
      out$replicate_id, "subsets"
    )
    out$unit_id <- .composition_validate_identifiers(out$unit_id, "subsets")
  } else {
    stop("subsets must be a data.frame or named list")
  }

  if (!nrow(out)) stop("subsets must contain at least one replicate")
  if (anyDuplicated(paste(out$replicate_id, out$unit_id, sep = "\r"))) {
    stop("subsets contains duplicate unit assignments")
  }
  rownames(out) <- NULL
  out
}

.composition_normalize_splits <- function(partitions) {
  if (!is.data.frame(partitions)) stop("partitions must be a data.frame")
  required <- c("replicate_id", "group", "unit_id")
  missing <- setdiff(required, names(partitions))
  if (length(missing)) {
    stop("partitions is missing required columns: ", paste(missing, collapse = ", "))
  }
  out <- partitions[required]
  out$replicate_id <- .composition_validate_identifiers(
    out$replicate_id, "partitions"
  )
  out$group <- toupper(.composition_validate_identifiers(out$group, "partitions"))
  out$unit_id <- .composition_validate_identifiers(out$unit_id, "partitions")
  if (!nrow(out)) stop("partitions must contain at least one replicate")
  if (any(!out$group %in% c("A", "B"))) {
    stop("partition group labels must be A or B")
  }
  if (anyDuplicated(paste(out$replicate_id, out$group, out$unit_id, sep = "\r"))) {
    stop("partitions contains duplicate unit assignments")
  }
  rownames(out) <- NULL
  out
}

.composition_validate_units <- function(requested, available, name) {
  units <- sort(unique(as.character(requested)))
  if (!length(units)) stop(name, " contains no biological units")
  missing <- setdiff(units, available)
  if (length(missing)) {
    stop(name, " refers to units absent from data: ", paste(missing, collapse = ", "))
  }
  units
}

.composition_build_chronology <- function(matrix_data, units, method, min_depth_fraction) {
  selected <- matrix_data[, units, drop = FALSE]
  depth <- rowSums(!is.na(selected))
  required_depth <- ceiling(min_depth_fraction * length(units))
  if (method == "mean") {
    chronology <- rowMeans(selected, na.rm = TRUE)
  } else {
    chronology <- apply(selected, 1L, stats::median, na.rm = TRUE)
  }
  chronology[depth == 0L] <- NA_real_
  chronology[depth < required_depth] <- NA_real_
  list(chronology = as.numeric(chronology), depth = as.integer(depth))
}

.composition_bound_correlation <- function(x) {
  tolerance <- 64 * .Machine$double.eps
  if (!is.finite(x)) stop("Pearson correlation is not finite")
  if (x > 1 && x <= 1 + tolerance) return(1)
  if (x < -1 && x >= -1 - tolerance) return(-1)
  if (x < -1 || x > 1) stop("Pearson correlation fell outside [-1, 1]")
  x
}

.composition_compare <- function(
    diagnostic,
    replicate_id,
    matrix_data,
    group_a,
    group_b,
    chronology_method,
    min_depth_fraction,
    min_overlap,
    full_pool_size) {
  a <- .composition_build_chronology(
    matrix_data, group_a, chronology_method, min_depth_fraction
  )
  b <- .composition_build_chronology(
    matrix_data, group_b, chronology_method, min_depth_fraction
  )
  valid <- !is.na(a$chronology) & !is.na(b$chronology)
  n_years <- sum(valid)
  years <- attr(matrix_data, "years")[valid]

  correlation <- NA_real_
  status <- "insufficient_overlap"
  if (n_years >= min_overlap) {
    range_a <- range(a$chronology[valid])
    range_b <- range(b$chronology[valid])
    if (diff(range_a) == 0 || diff(range_b) == 0) {
      status <- "constant_chronology"
    } else {
      correlation <- .composition_bound_correlation(
        stats::cor(a$chronology[valid], b$chronology[valid], method = "pearson")
      )
      status <- "ok"
    }
  }

  if (n_years > 0L) {
    mean_depth_a <- mean(a$depth[valid])
    min_depth_a <- min(a$depth[valid])
    mean_depth_b <- mean(b$depth[valid])
    min_depth_b <- min(b$depth[valid])
    start_year <- min(years)
    end_year <- max(years)
  } else {
    mean_depth_a <- NA_real_
    min_depth_a <- 0L
    mean_depth_b <- NA_real_
    min_depth_b <- 0L
    start_year <- NA_integer_
    end_year <- NA_integer_
  }

  shared <- length(intersect(group_a, group_b))
  data.frame(
    diagnostic = diagnostic,
    replicate_id = as.character(replicate_id),
    chronology_method = chronology_method,
    correlation_method = "pearson",
    min_depth_fraction = min_depth_fraction,
    min_overlap = min_overlap,
    full_pool_size = full_pool_size,
    group_a_size = length(group_a),
    group_b_size = length(group_b),
    shared_unit_count = shared,
    comparison_exhaustive = length(union(group_a, group_b)) == full_pool_size,
    correlation = correlation,
    squared_correlation = correlation^2,
    n_years = n_years,
    start_year = start_year,
    end_year = end_year,
    mean_depth_a = mean_depth_a,
    min_depth_a = min_depth_a,
    mean_depth_b = mean_depth_b,
    min_depth_b = min_depth_b,
    group_a_ids = paste(group_a, collapse = "|"),
    group_b_ids = paste(group_b, collapse = "|"),
    status = status,
    stringsAsFactors = FALSE
  )
}

#' Explicit subset-to-full chronology stability
#'
#' Compare user-supplied biological-unit subsets with the chronology formed from
#' the complete unit pool. Random subsets are intentionally generated outside
#' this function so that identical compositions can be audited in R and Python.
#'
#' @param data Long-form data containing one value per year and biological unit.
#' @param subsets A data frame with `replicate_id` and `unit_id`, or a named list
#'   whose elements contain unit IDs.
#' @param year_col Name of the calendar-year column.
#' @param unit_col Name of the biological-unit column, normally `tree_id`.
#' @param value_col Name of the numeric value column.
#' @param chronology_method Arithmetic `"mean"` or `"median"` chronology.
#' @param min_depth_fraction Minimum fraction of each group's units required in a
#'   year. The default 1 uses a complete-period comparison.
#' @param min_overlap Minimum number of common chronology years.
#' @return A data frame containing one auditable result row per replicate. Pearson
#'   correlation and squared correlation are reported together with overlap and
#'   depth information.
#' @export
subset_stability <- function(
    data,
    subsets,
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
  partition_table <- .composition_normalize_subsets(subsets)
  available <- colnames(matrix_data)
  full_units <- sort(available)

  replicate_ids <- sort(unique(partition_table$replicate_id))
  rows <- lapply(replicate_ids, function(replicate_id) {
    block <- partition_table[
      partition_table$replicate_id == replicate_id, , drop = FALSE
    ]
    subset_units <- .composition_validate_units(
      block$unit_id,
      available,
      paste0("subset replicate ", sQuote(replicate_id))
    )
    .composition_compare(
      "subset_to_full",
      replicate_id,
      matrix_data,
      subset_units,
      full_units,
      method,
      fraction,
      overlap,
      length(full_units)
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Explicit disjoint split-half chronology stability
#'
#' Compare user-supplied, non-overlapping A/B biological-unit partitions. By
#' default every replicate must assign the complete unit pool exactly once.
#'
#' @param data Long-form data containing one value per year and biological unit.
#' @param partitions A data frame with `replicate_id`, `group` (`A` or `B`), and
#'   `unit_id`.
#' @param year_col Name of the calendar-year column.
#' @param unit_col Name of the biological-unit column, normally `tree_id`.
#' @param value_col Name of the numeric value column.
#' @param chronology_method Arithmetic `"mean"` or `"median"` chronology.
#' @param min_depth_fraction Minimum fraction of each group's units required in a
#'   year. The default 1 uses a complete-period comparison.
#' @param min_overlap Minimum number of common chronology years.
#' @param require_exhaustive Require groups A and B together to contain every unit
#'   in `data` exactly once.
#' @param require_balanced Require group sizes to differ by at most one unit, as
#'   expected for a literal split half.
#' @return A data frame containing one auditable result row per partition. Pearson
#'   correlation and squared correlation are reported with overlap and depth data.
#' @export
split_half_stability <- function(
    data,
    partitions,
    year_col = "year",
    unit_col = "tree_id",
    value_col = "value",
    chronology_method = "mean",
    min_depth_fraction = 1,
    min_overlap = 3L,
    require_exhaustive = TRUE,
    require_balanced = TRUE) {
  if (length(require_exhaustive) != 1L || is.na(require_exhaustive) ||
      !is.logical(require_exhaustive)) {
    stop("require_exhaustive must be TRUE or FALSE")
  }
  if (length(require_balanced) != 1L || is.na(require_balanced) ||
      !is.logical(require_balanced)) {
    stop("require_balanced must be TRUE or FALSE")
  }
  method <- .composition_validate_method(chronology_method)
  fraction <- .composition_validate_fraction(min_depth_fraction)
  overlap <- .composition_validate_count(min_overlap, "min_overlap")
  matrix_data <- .composition_prepare_matrix(data, year_col, unit_col, value_col)
  partition_table <- .composition_normalize_splits(partitions)
  available <- colnames(matrix_data)

  replicate_ids <- sort(unique(partition_table$replicate_id))
  rows <- lapply(replicate_ids, function(replicate_id) {
    block <- partition_table[
      partition_table$replicate_id == replicate_id, , drop = FALSE
    ]
    if (!setequal(unique(block$group), c("A", "B"))) {
      stop("split-half replicate ", sQuote(replicate_id), " must contain groups A and B")
    }
    group_a <- .composition_validate_units(
      block$unit_id[block$group == "A"],
      available,
      paste0("split-half replicate ", sQuote(replicate_id), " group A")
    )
    group_b <- .composition_validate_units(
      block$unit_id[block$group == "B"],
      available,
      paste0("split-half replicate ", sQuote(replicate_id), " group B")
    )
    shared <- intersect(group_a, group_b)
    if (length(shared)) {
      stop(
        "split-half replicate ", sQuote(replicate_id),
        " is not disjoint: ", paste(shared, collapse = ", ")
      )
    }
    if (require_balanced && abs(length(group_a) - length(group_b)) > 1L) {
      stop(
        "split-half replicate ", sQuote(replicate_id),
        " is not balanced; group sizes are ", length(group_a),
        " and ", length(group_b)
      )
    }
    assigned <- union(group_a, group_b)
    if (require_exhaustive && !setequal(assigned, available)) {
      omitted <- setdiff(available, assigned)
      extra <- setdiff(assigned, available)
      stop(
        "split-half replicate ", sQuote(replicate_id),
        " is not exhaustive; omitted=", paste(omitted, collapse = ","),
        ", extra=", paste(extra, collapse = ",")
      )
    }
    .composition_compare(
      "disjoint_split_half",
      replicate_id,
      matrix_data,
      group_a,
      group_b,
      method,
      fraction,
      overlap,
      length(available)
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
