.apply_key_parser <- function(x, key_parser = NULL) {
  x$series_id <- as.character(x$series_id)
  if (is.null(key_parser)) {
    if (!"tree_id" %in% names(x)) x$tree_id <- x$series_id
    x$tree_id[is.na(x$tree_id)] <- x$series_id[is.na(x$tree_id)]
    return(x)
  }

  ids <- unique(x$series_id)
  parsed <- key_parser(ids)
  if (!is.data.frame(parsed) || !"series_id" %in% names(parsed)) {
    stop("key_parser must return a data.frame containing series_id")
  }
  parsed$series_id <- as.character(parsed$series_id)
  if (anyDuplicated(parsed$series_id)) stop("key_parser returned duplicate series_id values")
  if (!setequal(parsed$series_id, ids)) {
    stop("key_parser must return exactly one row for every input series_id")
  }
  if (!"tree_id" %in% names(parsed)) parsed$tree_id <- parsed$series_id
  parsed$tree_id[is.na(parsed$tree_id)] <- parsed$series_id[is.na(parsed$tree_id)]

  x$.chronodiag_row_id <- seq_len(nrow(x))
  original_tree <- if ("tree_id" %in% names(x)) x$tree_id else NULL
  if (!is.null(original_tree)) x$tree_id <- NULL
  out <- merge(x, parsed, by = "series_id", all.x = TRUE, sort = FALSE)
  out <- out[order(out$.chronodiag_row_id), , drop = FALSE]
  if (!is.null(original_tree)) {
    out$tree_id <- ifelse(is.na(original_tree), out$tree_id, original_tree)
  }
  out$.chronodiag_row_id <- NULL
  rownames(out) <- NULL
  out
}

.validate_year <- function(x) {
  out <- suppressWarnings(as.numeric(x))
  if (any(!is.finite(out)) || any(out != floor(out))) stop("year must contain finite integers")
  as.integer(out)
}

.validate_value <- function(x) {
  out <- suppressWarnings(as.numeric(x))
  bad <- !is.na(x) & is.na(out)
  if (any(bad)) stop("value contains non-numeric entries")
  if (any(!is.finite(out) & !is.na(out))) stop("value must contain finite numbers or missing values")
  out
}

.validate_unique_observations <- function(x) {
  key <- paste(x$series_id, x$year, sep = "\r")
  if (anyDuplicated(key)) stop("duplicate series_id-year observations detected")
  invisible(TRUE)
}

#' Convert wide year-by-series input to canonical long form
#'
#' @param data A data frame whose first column is named `year` and remaining
#'   columns are series keys.
#' @param key_parser Optional dataset-specific key parser.
#' @param drop_missing Drop missing observations from the long table.
#' @return Data frame with at least year, series_id, tree_id, and value.
#' @export
wide_to_long <- function(data, key_parser = NULL, drop_missing = TRUE) {
  if (!is.data.frame(data) || ncol(data) < 2L) {
    stop("wide input must contain year and at least one series column")
  }
  if (tolower(names(data)[1L]) != "year") stop("the first wide-CSV column must be named 'year'")
  if (anyDuplicated(names(data))) stop("wide input contains duplicate column names")
  names(data)[1L] <- "year"
  data$year <- .validate_year(data$year)
  if (anyDuplicated(data$year)) stop("wide input contains duplicate years")

  ids <- names(data)[-1L]
  values <- unlist(data[-1L], use.names = FALSE)
  out <- data.frame(
    year = rep(data$year, times = length(ids)),
    series_id = rep(ids, each = nrow(data)),
    value = .validate_value(values),
    stringsAsFactors = FALSE
  )
  if (drop_missing) out <- out[!is.na(out$value), , drop = FALSE]
  out <- .apply_key_parser(out, key_parser)
  .validate_unique_observations(out)
  front <- c("year", "series_id", "tree_id", "value")
  out <- out[c(front, setdiff(names(out), front))]
  out <- out[order(out$series_id, out$year), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Read a wide CSV file
#'
#' @param path Path to a CSV file whose first column is `year`.
#' @param key_parser Optional dataset-specific key parser.
#' @param drop_missing Drop missing observations from the returned long table.
#' @param ... Additional arguments passed to [utils::read.csv()].
#' @return Canonical long-form data frame.
#' @export
read_wide_csv <- function(path, key_parser = NULL, drop_missing = TRUE, ...) {
  x <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, ...)
  wide_to_long(x, key_parser = key_parser, drop_missing = drop_missing)
}

#' Read a long CSV file
#'
#' @param path Path to a CSV file containing at least `year`, `series_id`, and `value`.
#' @param key_parser Optional dataset-specific key parser.
#' @param drop_missing Drop missing observations from the returned long table.
#' @param ... Additional arguments passed to [utils::read.csv()].
#' @return Canonical long-form data frame.
#' @export
read_long_csv <- function(path, key_parser = NULL, drop_missing = TRUE, ...) {
  x <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, ...)
  needed <- c("year", "series_id", "value")
  missing <- setdiff(needed, names(x))
  if (length(missing)) stop("long CSV is missing required columns: ", paste(missing, collapse = ", "))
  x$series_id <- as.character(x$series_id)
  x$year <- .validate_year(x$year)
  x$value <- .validate_value(x$value)
  if (drop_missing) x <- x[!is.na(x$value), , drop = FALSE]
  x <- .apply_key_parser(x, key_parser)
  .validate_unique_observations(x)
  front <- c("year", "series_id", "tree_id", "value")
  x <- x[c(front, setdiff(names(x), front))]
  x <- x[order(x$series_id, x$year), , drop = FALSE]
  rownames(x) <- NULL
  x
}

.tucson_explicit_missing_cells <- function(path) {
  lines <- readLines(path, warn = FALSE)
  if (!length(lines)) {
    return(data.frame(series_id = character(), year = integer(),
                      line_index = integer(), field_start = integer(),
                      field_end = integer(), normalize = logical(),
                      stringsAsFactors = FALSE))
  }

  missing_series <- character()
  missing_year <- integer()
  missing_line <- integer()
  missing_start <- integer()
  missing_end <- integer()
  normalize <- logical()
  before_999 <- logical()
  before_9999 <- logical()
  final_year <- list()
  final_marker <- list()

  for (line_index in seq_along(lines)) {
    line <- lines[[line_index]]
    if (!nzchar(trimws(line)) || startsWith(trimws(line), "#")) next

    parsed <- FALSE
    series_id <- NA_character_
    start_year <- NA_integer_

    # A minus sign in column 8 indicates the 7-character-ID/5-character-year
    # Tucson layout.  Otherwise prefer the standard 8+4 layout.
    long_first <- nchar(line, type = "chars") >= 12L && identical(substr(line, 8L, 8L), "-")
    layouts <- if (long_first) list(c(7L, 5L), c(8L, 4L)) else list(c(8L, 4L), c(7L, 5L))

    for (layout in layouts) {
      if (parsed) break
      id_width <- layout[1L]
      year_width <- layout[2L]
      prefix_width <- id_width + year_width
      if (nchar(line, type = "chars") < prefix_width) next
      year_text <- trimws(substr(line, id_width + 1L, prefix_width))
      if (!grepl("^[+-]?[0-9]+$", year_text)) next
      candidate <- suppressWarnings(as.integer(year_text))
      if (is.na(candidate) || candidate < -10000L || candidate > 10000L) next
      series_id <- trimws(substr(line, 1L, id_width))
      start_year <- candidate
      parsed <- nzchar(series_id)
    }

    if (!parsed || nchar(line, type = "chars") <= 12L) next

    n_fields <- ceiling((nchar(line, type = "chars") - 12L) / 6L)
    content_end <- nchar(sub("[[:blank:]]+$", "", line), type = "chars")
    # A 999 value is a genuine 0.999 mm ring when the series ends in -9999.
    # Match the declared precision of the final dated record in each series.
    last_start <- 13L + 6L * floor((content_end - 13L) / 6L)
    last_token <- trimws(substr(line, last_start, content_end))
    if (is.null(final_year[[series_id]]) || start_year >= final_year[[series_id]]) {
      final_year[[series_id]] <- start_year
      final_marker[[series_id]] <- if (identical(last_token, "-9999")) "-9999" else "999"
    }
    terminated_999 <- FALSE
    terminated_9999 <- FALSE
    for (j in seq_len(n_fields)) {
      field_start <- 13L + 6L * (j - 1L)
      field_end <- min(field_start + 5L, nchar(line, type = "chars"))
      field <- substr(line, field_start, field_end)
      token <- trimws(field)
      if (identical(token, "999")) terminated_999 <- TRUE
      if (identical(token, "-9999")) terminated_9999 <- TRUE
      if (!nzchar(token) || identical(token, ".")) {
        missing_series <- c(missing_series, series_id)
        missing_year <- c(missing_year, start_year + j - 1L)
        missing_line <- c(missing_line, line_index)
        missing_start <- c(missing_start, field_start)
        missing_end <- c(missing_end, field_end)
        normalize <- c(normalize, field_start <= content_end)
        before_999 <- c(before_999, !terminated_999)
        before_9999 <- c(before_9999, !terminated_9999)
      }
    }
  }

  markers <- vapply(missing_series, function(id) final_marker[[id]], character(1))
  normalize <- normalize & ifelse(markers == "-9999", before_9999, before_999)
  data.frame(
    series_id = missing_series,
    year = missing_year,
    line_index = missing_line,
    field_start = missing_start,
    field_end = missing_end,
    normalize = normalize,
    stringsAsFactors = FALSE
  )
}

#' Read a Tucson/ITRDB RWL file
#'
#' @param path File path.
#' @param key_parser Optional dataset-specific key parser.
#' @param drop_missing Drop missing observations.
#' @param rwl_format Format passed to `dplR::read.rwl()`; default `"auto"`.
#' @param preserve_zero Preserve zero-width rings. This maps to dplR's
#'   `edge.zeros` argument and defaults to `TRUE`.
#' @param ... Additional arguments passed to `dplR::read.rwl()`.
#' @return Canonical long-form data frame with values in millimetres.
#' @export
read_tucson_rwl <- function(path, key_parser = NULL, drop_missing = TRUE,
                             rwl_format = "auto", preserve_zero = TRUE, ...) {
  if (!requireNamespace("dplR", quietly = TRUE)) {
    stop("Reading RWL files requires the suggested package 'dplR'")
  }
  if (length(preserve_zero) != 1L || is.na(preserve_zero) || !is.logical(preserve_zero)) {
    stop("preserve_zero must be TRUE or FALSE")
  }

  # Preserve explicit missing fields across dplR readers: older versions may
  # turn blanks into zero, while 1.8.0 rejects an otherwise aligned blank field
  # during its fixed-width/whitespace comparison. Use temporary numeric
  # placeholders only for these known missing cells, then restore their NA
  # values. Genuine zeros, stop markers and post-marker padding are untouched.
  explicit_missing <- .tucson_explicit_missing_cells(path)
  read_path <- path
  cells <- which(explicit_missing$normalize)
  if (length(cells)) {
    prepared_lines <- readLines(path, warn = FALSE)
    for (i in cells) {
      row <- explicit_missing$line_index[i]
      first <- explicit_missing$field_start[i]
      last <- explicit_missing$field_end[i]
      substr(prepared_lines[row], first, last) <- sprintf("%*d", last - first + 1L, 0L)
    }
    read_path <- tempfile(fileext = ".rwl")
    on.exit(unlink(read_path), add = TRUE)
    writeLines(prepared_lines, read_path, useBytes = TRUE)
  }

  wide <- dplR::read.rwl(
    read_path, format = rwl_format, edge.zeros = preserve_zero, ...
  )
  if (!is.data.frame(wide)) stop("dplR::read.rwl() did not return a wide rwl data frame")
  years <- suppressWarnings(as.integer(rownames(wide)))
  if (anyNA(years)) stop("RWL row names could not be converted to calendar years")
  tab <- data.frame(year = years, wide, check.names = FALSE)

  # Convert without dropping values so that explicit blank cells can be masked
  # before the caller's drop_missing policy is applied.
  out <- wide_to_long(tab, key_parser = key_parser, drop_missing = FALSE)
  if (nrow(explicit_missing)) {
    out_key <- paste(out$series_id, out$year, sep = "\r")
    missing_key <- paste(explicit_missing$series_id, explicit_missing$year, sep = "\r")
    out$value[out_key %in% missing_key] <- NA_real_
  }
  if (drop_missing) out <- out[!is.na(out$value), , drop = FALSE]

  out$unit <- "mm"
  out$source_format <- "tucson_rwl"
  front <- c("year", "series_id", "tree_id", "value", "unit", "source_format")
  out <- out[c(front, setdiff(names(out), front))]
  out <- out[order(out$series_id, out$year), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "explicit_missing_cells") <- nrow(explicit_missing)
  out
}

#' Read chronodiag input with format auto-detection
#'
#' @param path Input file path.
#' @param format One of `"auto"`, `"csv_wide"`, `"csv_long"`, or `"rwl"`.
#' @param key_parser Optional dataset-specific key parser.
#' @param drop_missing Drop missing observations from the returned long table.
#' @param ... Additional arguments passed to the selected reader.
#' @return Canonical long-form data frame.
#' @export
read_chronodiag <- function(path, format = "auto", key_parser = NULL,
                            drop_missing = TRUE, ...) {
  selected <- tolower(format)
  if (selected == "auto") {
    ext <- tolower(tools::file_ext(path))
    if (ext %in% c("rwl", "tuc")) {
      selected <- "rwl"
    } else if (ext == "csv") {
      header <- utils::read.csv(path, nrows = 1L, check.names = FALSE,
                                stringsAsFactors = FALSE)
      selected <- if (all(c("series_id", "value") %in% names(header))) "csv_long" else "csv_wide"
    } else {
      stop("cannot auto-detect format; use csv_wide, csv_long, or rwl")
    }
  }
  if (selected %in% c("csv", "csv_wide", "wide")) {
    return(read_wide_csv(path, key_parser = key_parser, drop_missing = drop_missing, ...))
  }
  if (selected %in% c("csv_long", "long")) {
    return(read_long_csv(path, key_parser = key_parser, drop_missing = drop_missing, ...))
  }
  if (selected %in% c("rwl", "tucson", "tucson_rwl")) {
    return(read_tucson_rwl(path, key_parser = key_parser, drop_missing = drop_missing, ...))
  }
  stop("unsupported format: ", format)
}
