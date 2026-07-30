#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
options(warn = 1)

parse_args <- function(args) {
  result <- list()
  index <- 1L
  while (index <= length(args)) {
    token <- args[[index]]
    if (!startsWith(token, "--")) stop("Unexpected argument: ", token)
    key <- sub("^--", "", token)
    if (index == length(args)) stop("Missing value for --", key)
    result[[key]] <- args[[index + 1L]]
    index <- index + 2L
  }
  result
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
if (length(file_arg) != 1L) stop("Could not resolve script path")
script_path <- normalizePath(
  sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE
)
repo_root <- if (!is.null(args[["repo-root"]])) {
  normalizePath(args[["repo-root"]], winslash = "/", mustWork = TRUE)
} else {
  normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
}
output_dir <- if (!is.null(args[["output-dir"]])) {
  normalizePath(args[["output-dir"]], winslash = "/", mustWork = FALSE)
} else {
  file.path(repo_root, "analysis", "yoshino_empirical_stage_c_v01", "preprocessed")
}
raw_path <- if (!is.null(args[["raw-path"]])) {
  normalizePath(args[["raw-path"]], winslash = "/", mustWork = TRUE)
} else {
  normalizePath(
    file.path(repo_root, "data", "yoshino", "raw", "nr_all.rwl"),
    winslash = "/", mustWork = TRUE
  )
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("dplR", quietly = TRUE)) {
  stop("Yoshino empirical Stage C requires the R package 'dplR'.")
}

contract_path <- file.path(
  repo_root, "data", "yoshino", "yoshino_dataset_contract_v01.csv"
)
keymap_path <- file.path(
  repo_root, "data", "yoshino", "yoshino_keymap_v01.csv"
)
if (!file.exists(contract_path)) stop("Missing contract: ", contract_path)
if (!file.exists(keymap_path)) stop("Missing key map: ", keymap_path)

contract <- utils::read.csv(contract_path, check.names = FALSE)[1L, , drop = FALSE]
keymap <- utils::read.csv(
  keymap_path, check.names = FALSE,
  colClasses = c(
    series_id = "character", tree_id = "character",
    specimen_id = "character", source_group = "character",
    site_id = "character"
  )
)
keymap$include_primary_cohort <-
  tolower(as.character(keymap$include_primary_cohort)) == "true"

if (!identical(unique(keymap$input_level), "specimen_mean")) {
  stop("Yoshino key map must contain one specimen-mean series per specimen")
}
if (anyDuplicated(keymap$series_id) || anyDuplicated(keymap$tree_id)) {
  stop("Yoshino specimen-mean IDs must be one-to-one")
}

message(
  "Reading Yoshino specimen means from ", raw_path,
  " [format=tucson, header=FALSE, long=FALSE, preserve_zero=TRUE]"
)
rwl <- dplR::read.rwl(
  fname = raw_path,
  format = "tucson",
  header = FALSE,
  long = FALSE,
  edge.zeros = TRUE
)

if (!setequal(colnames(rwl), keymap$series_id)) {
  stop("Raw Yoshino RWL and key map contain different series IDs")
}
if (ncol(rwl) != as.integer(contract$n_raw_series_expected)) {
  stop("Unexpected number of Yoshino specimen-mean series")
}

selected_map <- keymap[keymap$include_primary_cohort, , drop = FALSE]
selected_series <- sort(selected_map$series_id)
selected_trees <- sort(selected_map$tree_id)
if (length(selected_trees) != as.integer(contract$n_trees_expected)) {
  stop("Selected Yoshino specimen count differs from contract")
}
if (!identical(selected_series, selected_trees)) {
  stop("Yoshino Stage C expects one specimen-mean series per tree_id")
}

selected_rwl <- rwl[, selected_series, drop = FALSE]
raw_years <- as.integer(rownames(selected_rwl))
analysis_start <- as.integer(contract$analysis_start)
analysis_end <- as.integer(contract$analysis_end)
expected_years <- seq.int(analysis_start, analysis_end)
period_mask <- raw_years >= analysis_start & raw_years <= analysis_end
if (!identical(raw_years[period_mask], expected_years)) {
  stop("Yoshino analysis years are not contiguous")
}
if (anyNA(selected_rwl[period_mask, , drop = FALSE])) {
  stop("Yoshino primary interval contains missing specimen-mean widths")
}

method_specs <- list(
  AgeDepSpline = list(method = "AgeDepSpline", nyrs = NULL, role = "primary"),
  ModNegExp = list(method = "ModNegExp", nyrs = NULL, role = "sensitivity"),
  Spline30 = list(method = "Spline", nyrs = 30, role = "sensitivity")
)

detrend_with_warnings <- function(rwl_input, method_spec) {
  warning_messages <- character()
  value <- withCallingHandlers(
    {
      if (is.null(method_spec$nyrs)) {
        dplR::detrend(
          rwl = rwl_input,
          method = method_spec$method,
          make.plot = FALSE,
          verbose = FALSE
        )
      } else {
        dplR::detrend(
          rwl = rwl_input,
          method = method_spec$method,
          nyrs = method_spec$nyrs,
          f = 0.5,
          make.plot = FALSE,
          verbose = FALSE
        )
      }
    },
    warning = function(w) {
      warning_messages <<- c(warning_messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = unique(warning_messages))
}

tree_results <- list()
specimen_results <- list()
summary_results <- list()
warning_results <- list()

for (method_label in names(method_specs)) {
  method_spec <- method_specs[[method_label]]
  message("Detrending Yoshino specimen means with ", method_label)
  detrended <- tryCatch(
    detrend_with_warnings(selected_rwl, method_spec),
    error = function(e) e
  )

  if (inherits(detrended, "error")) {
    if (identical(method_spec$role, "primary")) {
      stop("Yoshino primary detrending failed: ", conditionMessage(detrended))
    }
    warning(method_label, " sensitivity failed: ", conditionMessage(detrended))
    summary_results[[length(summary_results) + 1L]] <- data.frame(
      dataset = "yoshino",
      detrending_method = method_label,
      method_role = method_spec$role,
      status = "detrending_failed",
      failure_message = conditionMessage(detrended),
      analysis_start = analysis_start,
      analysis_end = analysis_end,
      analysis_years = length(expected_years),
      n_trees = length(selected_trees),
      n_series = length(selected_series),
      input_level = "specimen_mean",
      warning_count = NA_integer_,
      stringsAsFactors = FALSE
    )
    next
  }

  rwi <- as.data.frame(detrended$value, check.names = FALSE)
  if (!identical(colnames(rwi), colnames(selected_rwl))) {
    stop("Detrending changed Yoshino specimen IDs")
  }
  rwi_years <- as.integer(rownames(rwi))
  rwi_period <- rwi[
    rwi_years >= analysis_start & rwi_years <= analysis_end,
    selected_series,
    drop = FALSE
  ]
  if (!identical(as.integer(rownames(rwi_period)), expected_years)) {
    stop("Detrended Yoshino years do not match the fixed interval")
  }
  if (anyNA(rwi_period) || any(!is.finite(as.matrix(rwi_period)))) {
    stop("Detrended Yoshino interval contains invalid values")
  }

  specimen_matrix <- as.matrix(rwi_period)
  specimen_means <- colMeans(specimen_matrix)
  specimen_sds <- apply(specimen_matrix, 2L, stats::sd)
  if (any(!is.finite(specimen_sds)) || any(specimen_sds <= 0)) {
    stop("One or more Yoshino specimen-mean series are constant")
  }
  specimen_z <- sweep(specimen_matrix, 2L, specimen_means, FUN = "-")
  specimen_z <- sweep(specimen_z, 2L, specimen_sds, FUN = "/")

  tree_rows <- vector("list", length(selected_trees))
  specimen_rows <- vector("list", length(selected_trees))
  for (position in seq_along(selected_trees)) {
    tree_id <- selected_trees[[position]]
    meta <- selected_map[selected_map$tree_id == tree_id, , drop = FALSE]
    tree_rows[[position]] <- data.frame(
      dataset = "yoshino",
      detrending_method = method_label,
      method_role = method_spec$role,
      year = expected_years,
      tree_id = tree_id,
      specimen_id = tree_id,
      site_id = as.character(meta$site_id),
      source_group = as.character(meta$source_group),
      value_rwi = as.numeric(specimen_matrix[, position]),
      value_z = as.numeric(specimen_z[, position]),
      n_input_series = 1L,
      input_level = "specimen_mean",
      series_ids = tree_id,
      stringsAsFactors = FALSE
    )
    specimen_rows[[position]] <- data.frame(
      dataset = "yoshino",
      detrending_method = method_label,
      method_role = method_spec$role,
      year = expected_years,
      series_id = tree_id,
      tree_id = tree_id,
      specimen_id = tree_id,
      site_id = as.character(meta$site_id),
      source_group = as.character(meta$source_group),
      value_rwi = as.numeric(specimen_matrix[, position]),
      input_level = "specimen_mean",
      stringsAsFactors = FALSE
    )
  }

  tree_results[[length(tree_results) + 1L]] <- do.call(rbind, tree_rows)
  specimen_results[[length(specimen_results) + 1L]] <- do.call(rbind, specimen_rows)
  summary_results[[length(summary_results) + 1L]] <- data.frame(
    dataset = "yoshino",
    detrending_method = method_label,
    method_role = method_spec$role,
    status = "ok",
    failure_message = "",
    analysis_start = analysis_start,
    analysis_end = analysis_end,
    analysis_years = length(expected_years),
    n_trees = length(selected_trees),
    n_series = length(selected_series),
    input_level = "specimen_mean",
    n_zero_rwi = sum(specimen_matrix == 0),
    n_missing_tree_values = sum(is.na(specimen_matrix)),
    minimum_tree_sd = min(specimen_sds),
    maximum_tree_sd = max(specimen_sds),
    warning_count = length(detrended$warnings),
    stringsAsFactors = FALSE
  )
  if (length(detrended$warnings)) {
    warning_results[[length(warning_results) + 1L]] <- data.frame(
      dataset = "yoshino",
      detrending_method = method_label,
      warning_message = detrended$warnings,
      stringsAsFactors = FALSE
    )
  }
}

tree_table <- do.call(rbind, tree_results)
specimen_table <- do.call(rbind, specimen_results)
summary_table <- do.call(rbind, summary_results)
warning_table <- if (length(warning_results)) {
  do.call(rbind, warning_results)
} else {
  data.frame(
    dataset = character(),
    detrending_method = character(),
    warning_message = character(),
    stringsAsFactors = FALSE
  )
}

utils::write.csv(
  tree_table,
  file.path(output_dir, "tree_level_all_methods_v01.csv"),
  row.names = FALSE, na = "", eol = "\n"
)
utils::write.csv(
  specimen_table,
  file.path(output_dir, "specimen_level_all_methods_v01.csv"),
  row.names = FALSE, na = "", eol = "\n"
)
utils::write.csv(
  summary_table,
  file.path(output_dir, "preprocessing_summary_v01.csv"),
  row.names = FALSE, na = "", eol = "\n"
)
utils::write.csv(
  warning_table,
  file.path(output_dir, "detrending_warnings_v01.csv"),
  row.names = FALSE, na = "", eol = "\n"
)
utils::write.csv(
  data.frame(
    item = c("R", "dplR"),
    version = c(
      paste(R.version$major, R.version$minor, sep = "."),
      as.character(utils::packageVersion("dplR"))
    ),
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "preprocessing_environment_v01.csv"),
  row.names = FALSE, eol = "\n"
)

message("")
message("Yoshino empirical Stage C preprocessing completed.")
message("Tree-level table: ", file.path(output_dir, "tree_level_all_methods_v01.csv"))
message("Summary: ", file.path(output_dir, "preprocessing_summary_v01.csv"))
