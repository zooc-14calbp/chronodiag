#!/usr/bin/env Rscript

# Prepare frozen tree-level ITRDB data for empirical Stage B.

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
default_repo <- normalizePath(
  file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE
)

repo_root <- if (!is.null(args[["repo-root"]])) {
  normalizePath(args[["repo-root"]], winslash = "/", mustWork = TRUE)
} else {
  default_repo
}
project_root <- if (!is.null(args[["project-root"]])) {
  normalizePath(args[["project-root"]], winslash = "/", mustWork = TRUE)
} else {
  normalizePath(
    file.path(repo_root, "..", ".."), winslash = "/", mustWork = TRUE
  )
}
output_dir <- if (!is.null(args[["output-dir"]])) {
  normalizePath(args[["output-dir"]], winslash = "/", mustWork = FALSE)
} else {
  file.path(
    repo_root, "analysis", "itrdb_empirical_stage_b_v01", "preprocessed"
  )
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("dplR", quietly = TRUE)) {
  stop("The empirical preprocessing requires the R package 'dplR'.")
}

contract_path <- file.path(
  repo_root, "data", "itrdb", "empirical_dataset_contract_v01.csv"
)
if (!file.exists(contract_path)) {
  stop("Dataset contract not found: ", contract_path)
}
contracts <- utils::read.csv(contract_path, check.names = FALSE)

method_specs <- list(
  AgeDepSpline = list(method = "AgeDepSpline", nyrs = NULL, role = "primary"),
  ModNegExp = list(method = "ModNegExp", nyrs = NULL, role = "sensitivity"),
  Spline30 = list(method = "Spline", nyrs = 30, role = "sensitivity")
)


parse_contract_bool <- function(value, field, dataset) {
  text <- tolower(trimws(as.character(value)))
  if (length(text) != 1L || is.na(text) || !text %in% c("true", "false")) {
    stop(
      dataset, ": contract field '", field,
      "' must be exactly true or false"
    )
  }
  identical(text, "true")
}

read_empirical_rwl <- function(raw_path, contract) {
  dataset <- as.character(contract$dataset)
  rwl_format <- tolower(trimws(as.character(contract$rwl_format)))
  if (length(rwl_format) != 1L || is.na(rwl_format) ||
      !identical(rwl_format, "tucson")) {
    stop(
      dataset,
      ": empirical Stage B currently requires rwl_format='tucson'"
    )
  }

  tucson_header <- parse_contract_bool(
    contract$tucson_header, "tucson_header", dataset
  )
  tucson_long <- parse_contract_bool(
    contract$tucson_long, "tucson_long", dataset
  )
  preserve_zero <- parse_contract_bool(
    contract$preserve_zero, "preserve_zero", dataset
  )

  message(
    "Reading ", dataset, " from ", raw_path,
    " [format=", rwl_format,
    ", header=", tucson_header,
    ", long=", tucson_long,
    ", preserve_zero=", preserve_zero, "]"
  )

  value <- dplR::read.rwl(
    fname = raw_path,
    format = rwl_format,
    header = tucson_header,
    long = tucson_long,
    edge.zeros = preserve_zero
  )

  attr(value, "chronodiag_input_format") <- rwl_format
  attr(value, "chronodiag_tucson_header") <- tucson_header
  attr(value, "chronodiag_tucson_long") <- tucson_long
  attr(value, "chronodiag_preserve_zero") <- preserve_zero
  value
}

detrend_with_warnings <- function(rwl, method_spec) {
  warning_messages <- character()
  value <- withCallingHandlers(
    {
      if (is.null(method_spec$nyrs)) {
        dplR::detrend(
          rwl = rwl,
          method = method_spec$method,
          make.plot = FALSE,
          verbose = FALSE
        )
      } else {
        dplR::detrend(
          rwl = rwl,
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

process_method <- function(
    dataset,
    selected_rwl,
    selected_map,
    selected_trees,
    selected_series,
    analysis_start,
    analysis_end,
    method_label,
    method_spec,
    input_format,
    tucson_header,
    tucson_long,
    preserve_zero) {

  detrended <- detrend_with_warnings(selected_rwl, method_spec)
  rwi <- as.data.frame(detrended$value, check.names = FALSE)
  if (!identical(colnames(rwi), colnames(selected_rwl))) {
    stop("detrending changed series IDs")
  }

  expected_years <- seq.int(analysis_start, analysis_end)
  rwi_years <- as.integer(rownames(rwi))
  rwi_period <- rwi[
    rwi_years >= analysis_start & rwi_years <= analysis_end,
    selected_series,
    drop = FALSE
  ]
  if (!identical(as.integer(rownames(rwi_period)), expected_years)) {
    stop("detrended years do not match the fixed interval")
  }
  if (anyNA(rwi_period)) {
    stop("selected interval contains missing RWI values")
  }
  if (any(!is.finite(as.matrix(rwi_period)))) {
    stop("selected interval contains non-finite RWI values")
  }

  tree_matrix <- matrix(
    NA_real_,
    nrow = nrow(rwi_period),
    ncol = length(selected_trees),
    dimnames = list(rownames(rwi_period), selected_trees)
  )
  tree_metadata <- vector("list", length(selected_trees))

  for (tree_position in seq_along(selected_trees)) {
    tree_id <- selected_trees[[tree_position]]
    series_ids <- sort(
      selected_map$series_id[selected_map$tree_id == tree_id]
    )
    tree_matrix[, tree_position] <- rowMeans(
      rwi_period[, series_ids, drop = FALSE],
      na.rm = FALSE
    )
    tree_metadata[[tree_position]] <- data.frame(
      tree_id = tree_id,
      n_cores = length(series_ids),
      series_ids = paste(series_ids, collapse = "|"),
      stringsAsFactors = FALSE
    )
  }
  tree_metadata <- do.call(rbind, tree_metadata)

  if (anyNA(tree_matrix) || any(!is.finite(tree_matrix))) {
    stop("tree aggregation produced invalid values")
  }
  tree_means <- colMeans(tree_matrix)
  tree_sds <- apply(tree_matrix, 2L, stats::sd)
  if (any(!is.finite(tree_sds)) || any(tree_sds <= 0)) {
    stop("one or more tree means are constant")
  }
  tree_z <- sweep(tree_matrix, 2L, tree_means, FUN = "-")
  tree_z <- sweep(tree_z, 2L, tree_sds, FUN = "/")

  tree_rows <- vector("list", length(selected_trees))
  for (tree_position in seq_along(selected_trees)) {
    tree_id <- selected_trees[[tree_position]]
    meta <- tree_metadata[
      tree_metadata$tree_id == tree_id, , drop = FALSE
    ]
    tree_rows[[tree_position]] <- data.frame(
      dataset = dataset,
      detrending_method = method_label,
      method_role = method_spec$role,
      year = expected_years,
      tree_id = tree_id,
      value_rwi = as.numeric(tree_matrix[, tree_position]),
      value_z = as.numeric(tree_z[, tree_position]),
      n_cores = as.integer(meta$n_cores),
      series_ids = as.character(meta$series_ids),
      stringsAsFactors = FALSE
    )
  }

  core_rows <- vector("list", length(selected_series))
  for (series_position in seq_along(selected_series)) {
    series_id <- selected_series[[series_position]]
    map_row <- selected_map[
      selected_map$series_id == series_id, , drop = FALSE
    ]
    core_rows[[series_position]] <- data.frame(
      dataset = dataset,
      detrending_method = method_label,
      method_role = method_spec$role,
      year = expected_years,
      series_id = series_id,
      tree_id = as.character(map_row$tree_id),
      core_id = as.character(map_row$core_id),
      value_rwi = as.numeric(rwi_period[, series_id]),
      stringsAsFactors = FALSE
    )
  }

  warning_rows <- if (length(detrended$warnings)) {
    data.frame(
      dataset = dataset,
      detrending_method = method_label,
      warning_message = detrended$warnings,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      dataset = character(),
      detrending_method = character(),
      warning_message = character(),
      stringsAsFactors = FALSE
    )
  }

  summary <- data.frame(
    dataset = dataset,
    detrending_method = method_label,
    method_role = method_spec$role,
    status = "ok",
    failure_message = "",
    analysis_start = analysis_start,
    analysis_end = analysis_end,
    analysis_years = length(expected_years),
    n_trees = length(selected_trees),
    n_series = length(selected_series),
    input_format = input_format,
    tucson_header = tucson_header,
    tucson_long = tucson_long,
    preserve_zero = preserve_zero,
    n_zero_rwi = sum(as.matrix(rwi_period) == 0),
    n_missing_tree_values = sum(is.na(tree_matrix)),
    minimum_tree_sd = min(tree_sds),
    maximum_tree_sd = max(tree_sds),
    warning_count = length(detrended$warnings),
    stringsAsFactors = FALSE
  )

  list(
    tree = do.call(rbind, tree_rows),
    core = do.call(rbind, core_rows),
    summary = summary,
    warnings = warning_rows
  )
}

tree_results <- list()
core_results <- list()
summary_results <- list()
warning_results <- list()

for (contract_index in seq_len(nrow(contracts))) {
  contract <- contracts[contract_index, , drop = FALSE]
  dataset <- as.character(contract$dataset)
  raw_path <- file.path(
    project_root, "data", "itrdb_raw", as.character(contract$raw_filename)
  )
  keymap_path <- file.path(
    repo_root, "data", "itrdb", "keymaps",
    paste0(dataset, "_keymap_v01.csv")
  )

  if (!file.exists(raw_path)) stop(dataset, ": raw RWL not found: ", raw_path)
  if (!file.exists(keymap_path)) stop(dataset, ": key map not found: ", keymap_path)

  keymap <- utils::read.csv(
    keymap_path,
    check.names = FALSE,
    colClasses = "character"
  )
  keymap$include_primary_cohort <-
    tolower(keymap$include_primary_cohort) == "true"
  selected_map <- keymap[keymap$include_primary_cohort, , drop = FALSE]
  selected_series <- sort(unique(selected_map$series_id))
  selected_trees <- sort(unique(selected_map$tree_id))

  if (length(selected_trees) != as.integer(contract$n_trees_expected)) {
    stop(dataset, ": selected tree count differs from contract")
  }
  if (length(selected_series) != as.integer(contract$n_series_expected)) {
    stop(dataset, ": selected series count differs from contract")
  }

  rwl <- read_empirical_rwl(raw_path, contract)
  input_format <- as.character(attr(rwl, "chronodiag_input_format"))
  tucson_header <- isTRUE(attr(rwl, "chronodiag_tucson_header"))
  tucson_long <- isTRUE(attr(rwl, "chronodiag_tucson_long"))
  preserve_zero <- isTRUE(attr(rwl, "chronodiag_preserve_zero"))
  missing_series <- setdiff(selected_series, colnames(rwl))
  unmapped_series <- setdiff(colnames(rwl), keymap$series_id)
  if (length(missing_series) > 0L || length(unmapped_series) > 0L) {
    stop(
      dataset, ": raw/key-map mismatch; missing selected series=",
      paste(missing_series, collapse = "|"),
      "; unmapped raw series=", paste(unmapped_series, collapse = "|")
    )
  }

  selected_rwl <- rwl[, selected_series, drop = FALSE]
  raw_years <- as.integer(rownames(selected_rwl))
  analysis_start <- as.integer(contract$analysis_start)
  analysis_end <- as.integer(contract$analysis_end)
  period_mask <- raw_years >= analysis_start & raw_years <= analysis_end
  expected_years <- seq.int(analysis_start, analysis_end)
  if (!identical(raw_years[period_mask], expected_years)) {
    stop(dataset, ": selected analysis years are not contiguous")
  }
  if (anyNA(selected_rwl[period_mask, , drop = FALSE])) {
    stop(dataset, ": selected complete-core interval contains missing raw values")
  }

  for (method_label in names(method_specs)) {
    method_spec <- method_specs[[method_label]]
    message("Detrending ", dataset, " with ", method_label)
    result <- tryCatch(
      process_method(
        dataset = dataset,
        selected_rwl = selected_rwl,
        selected_map = selected_map,
        selected_trees = selected_trees,
        selected_series = selected_series,
        analysis_start = analysis_start,
        analysis_end = analysis_end,
        method_label = method_label,
        method_spec = method_spec,
        input_format = input_format,
        tucson_header = tucson_header,
        tucson_long = tucson_long,
        preserve_zero = preserve_zero
      ),
      error = function(e) e
    )

    if (inherits(result, "error")) {
      if (identical(method_spec$role, "primary")) {
        stop(
          dataset, " primary detrending failed: ",
          conditionMessage(result)
        )
      }
      warning(
        dataset, " ", method_label, " sensitivity failed: ",
        conditionMessage(result)
      )
      summary_results[[length(summary_results) + 1L]] <- data.frame(
        dataset = dataset,
        detrending_method = method_label,
        method_role = method_spec$role,
        status = "detrending_failed",
        failure_message = conditionMessage(result),
        analysis_start = analysis_start,
        analysis_end = analysis_end,
        analysis_years = length(expected_years),
        n_trees = length(selected_trees),
        n_series = length(selected_series),
        input_format = input_format,
        tucson_header = tucson_header,
        tucson_long = tucson_long,
        preserve_zero = preserve_zero,
        n_zero_rwi = NA_integer_,
        n_missing_tree_values = NA_integer_,
        minimum_tree_sd = NA_real_,
        maximum_tree_sd = NA_real_,
        warning_count = NA_integer_,
        stringsAsFactors = FALSE
      )
      next
    }

    tree_results[[length(tree_results) + 1L]] <- result$tree
    core_results[[length(core_results) + 1L]] <- result$core
    summary_results[[length(summary_results) + 1L]] <- result$summary
    if (nrow(result$warnings) > 0L) {
      warning_results[[length(warning_results) + 1L]] <- result$warnings
    }
  }
}

tree_table <- do.call(rbind, tree_results)
core_table <- do.call(rbind, core_results)
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

tree_path <- file.path(output_dir, "tree_level_all_methods_v01.csv")
core_path <- file.path(output_dir, "core_level_all_methods_v01.csv")
summary_path <- file.path(output_dir, "preprocessing_summary_v01.csv")
warning_path <- file.path(output_dir, "detrending_warnings_v01.csv")
environment_path <- file.path(output_dir, "preprocessing_environment_v01.csv")

utils::write.csv(tree_table, tree_path, row.names = FALSE, na = "", eol = "\n")
utils::write.csv(core_table, core_path, row.names = FALSE, na = "", eol = "\n")
utils::write.csv(summary_table, summary_path, row.names = FALSE, na = "", eol = "\n")
utils::write.csv(warning_table, warning_path, row.names = FALSE, na = "", eol = "\n")
utils::write.csv(
  data.frame(
    item = c("R", "dplR"),
    version = c(
      paste(R.version$major, R.version$minor, sep = "."),
      as.character(utils::packageVersion("dplR"))
    ),
    stringsAsFactors = FALSE
  ),
  environment_path,
  row.names = FALSE,
  eol = "\n"
)

message("")
message("ITRDB empirical preprocessing completed.")
message("Tree-level table: ", tree_path)
message("Core-level table: ", core_path)
message("Summary: ", summary_path)
message("Warnings: ", warning_path)
message("Environment: ", environment_path)
