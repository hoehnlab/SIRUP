# RAD and pRESTO sequence-header contracts ----------------------------------

.sirap_header_capture <- function(x, pattern) {
  matched <- regexec(pattern, x, perl = TRUE, ignore.case = TRUE)
  pieces <- regmatches(x, matched)
  vapply(pieces, function(piece) {
    if (length(piece) < 2L) NA_character_ else piece[[2L]]
  }, character(1))
}

.sirap_duplicate_presto_fields <- function(header) {
  token <- sub("[[:space:]].*$", "", sub("^[@>]", "", header))
  parts <- strsplit(token, "|", fixed = TRUE)[[1L]]
  if (length(parts) < 2L) {
    return(character())
  }
  fields <- sub("=.*$", "", parts[-1L])
  fields <- toupper(fields[grepl("=", parts[-1L], fixed = TRUE)])
  unique(fields[duplicated(fields)])
}

.sirap_validate_presto_component <- function(x, name, allow_missing = FALSE) {
  x <- as.character(x)
  invalid_missing <- is.na(x) | !nzchar(x)
  if (!allow_missing && any(invalid_missing)) {
    rlang::abort(paste0("`", name, "` cannot contain missing or empty values."))
  }
  check <- !invalid_missing
  invalid_reserved <- check & grepl("[|=,[:space:]]", x)
  if (any(invalid_reserved)) {
    first <- which(invalid_reserved)[[1L]]
    rlang::abort(paste0(
      "`", name, "` contains a pRESTO delimiter or whitespace at row ",
      first, ": ", encodeString(x[[first]], quote = "\""), "."
    ))
  }
  invisible(x)
}

#' Parse RAD metadata from sequence identifiers
#'
#' Parses the canonical RAD/pRESTO identifier written for IgBLAST,
#' `read_id|BARCODE=value|UMI=value`. For compatibility at input boundaries,
#' legacy `|barcode:value|umi:value` fields and native RAD SAM-style comments
#' such as `CB:Z:value UB:Z:value` are also recognized. Newly written SIRAP
#' files always use the canonical form.
#'
#' @param airr_df AIRR-style data frame.
#' @param id Column containing sequence identifiers or full headers.
#' @param strip Replace `id` with the bare read identifier before its first
#'   annotation field.
#' @param require_barcode Require one non-empty spatial barcode per row.
#' @return A data frame with `rad_read_id`, `barcode`, `umi`, `raw_barcode`,
#'   `rad_cb`, `rad_bc`, `rad_cr`, `rad_ub`, `rad_ts`, and `rad_rg` columns.
#' @export
parseRadSequenceIds <- function(airr_df,
                                id = "sequence_id",
                                strip = TRUE,
                                require_barcode = FALSE) {
  .check_required_cols(airr_df, id, "airr_df")
  out <- .as_data_frame(airr_df)
  headers <- as.character(out[[id]])
  if (anyNA(headers) || any(!nzchar(headers))) {
    rlang::abort(paste0("`", id, "` cannot contain missing or empty values."))
  }

  duplicates <- lapply(headers, .sirap_duplicate_presto_fields)
  bad <- which(lengths(duplicates) > 0L)
  if (length(bad)) {
    rlang::abort(paste0(
      "Duplicate pRESTO field in `", id, "` at row ", bad[[1L]], ": ",
      paste(duplicates[[bad[[1L]]]], collapse = ", "), "."
    ))
  }

  clean <- sub("^[@>]", "", headers)
  token <- sub("[[:space:]].*$", "", clean)
  bare_id <- sub("\\|.*$", "", token)

  canonical_barcode <- .sirap_header_capture(
    token, "(?:^|\\|)BARCODE=([^|]*)"
  )
  legacy_barcode <- .sirap_header_capture(
    token, "(?:^|\\|)barcode:([^|]*)"
  )
  canonical_umi <- .sirap_header_capture(token, "(?:^|\\|)UMI=([^|]*)")
  legacy_umi <- .sirap_header_capture(token, "(?:^|\\|)umi:([^|]*)")

  sam_capture <- function(tag) {
    .sirap_header_capture(
      clean,
      paste0("(?:^|[\\t ])", tag, ":[A-Za-z]:([^\\t ]*)")
    )
  }
  rad_cb <- sam_capture("CB")
  rad_bc <- sam_capture("BC")
  rad_cr <- sam_capture("CR")
  rad_ub <- sam_capture("UB")
  rad_ts <- sam_capture("TS")
  rad_rg <- sam_capture("RG")

  first_present <- function(...) {
    values <- list(...)
    answer <- rep(NA_character_, length(headers))
    for (value in values) {
      take <- is.na(answer) & !is.na(value) & nzchar(value)
      answer[take] <- value[take]
    }
    answer
  }

  existing_barcode <- if ("barcode" %in% colnames(out)) {
    as.character(out$barcode)
  } else {
    rep(NA_character_, length(headers))
  }
  existing_umi <- if ("umi" %in% colnames(out)) {
    as.character(out$umi)
  } else {
    rep(NA_character_, length(headers))
  }

  # BARCODE is the explicit interchange field. BC is the spatial bin emitted
  # by coordinate reformatting; otherwise the corrected RAD barcode is CB.
  # Already-tabulated AIRR columns are retained as the final fallback.
  barcode <- first_present(
    canonical_barcode, legacy_barcode, rad_bc, rad_cb, existing_barcode
  )
  umi <- first_present(canonical_umi, legacy_umi, rad_ub, existing_umi)

  if (isTRUE(require_barcode) && any(is.na(barcode) | !nzchar(barcode))) {
    missing <- which(is.na(barcode) | !nzchar(barcode))
    rlang::abort(paste0(
      "Missing spatial barcode in `", id, "` at row ", missing[[1L]],
      ". Expected BARCODE, BC, or CB metadata."
    ))
  }

  out$rad_read_id <- bare_id
  out$barcode <- barcode
  out$umi <- umi
  out$raw_barcode <- rad_cr
  out$rad_cb <- rad_cb
  out$rad_bc <- rad_bc
  out$rad_cr <- rad_cr
  out$rad_ub <- rad_ub
  out$rad_ts <- rad_ts
  out$rad_rg <- rad_rg
  if (isTRUE(strip)) {
    out[[id]] <- bare_id
  }
  out
}

.sirap_format_spatial_headers <- function(id, barcode, umi = NULL,
                                          barcode_key = "BARCODE") {
  id <- as.character(id)
  barcode <- as.character(barcode)
  if (length(id) != length(barcode)) {
    rlang::abort("`id` and `barcode` must have the same length.")
  }
  if (is.null(umi)) {
    umi <- rep(NA_character_, length(id))
  } else {
    umi <- as.character(umi)
    if (length(umi) != length(id)) {
      rlang::abort("`umi` must be NULL or have the same length as `id`.")
    }
  }

  .sirap_validate_presto_component(id, "id")
  .sirap_validate_presto_component(barcode, "barcode")
  .sirap_validate_presto_component(umi, "umi", allow_missing = TRUE)
  if (!is.character(barcode_key) || length(barcode_key) != 1L ||
      is.na(barcode_key) || !grepl("^[A-Za-z][A-Za-z0-9_]*$", barcode_key)) {
    rlang::abort(
      "`barcode_key` must be one pRESTO field name using letters, numbers, or underscores."
    )
  }
  barcode_key <- toupper(barcode_key)
  if (anyDuplicated(id)) {
    duplicate <- id[duplicated(id)][[1L]]
    rlang::abort(paste0(
      "Sequence identifiers must be unique before pRESTO; duplicate: ",
      encodeString(duplicate, quote = "\""), "."
    ))
  }

  header <- paste0(id, "|", barcode_key, "=", barcode)
  has_umi <- !is.na(umi) & nzchar(umi)
  header[has_umi] <- paste0(header[has_umi], "|UMI=", umi[has_umi])
  header
}

.sirap_write_spatial_fasta <- function(airr_df,
                                       path,
                                       id = "sequence_id",
                                       barcode = "barcode",
                                       umi = "umi",
                                       sequence = "sequence_alignment",
                                       barcode_key = "BARCODE") {
  required <- c(id, barcode, sequence)
  .check_required_cols(airr_df, required, "airr_df")
  umi_values <- if (!is.null(umi) && umi %in% colnames(airr_df)) {
    airr_df[[umi]]
  } else {
    NULL
  }
  headers <- .sirap_format_spatial_headers(
    airr_df[[id]], airr_df[[barcode]], umi_values,
    barcode_key = barcode_key
  )
  sequences <- as.character(airr_df[[sequence]])
  if (anyNA(sequences) || any(!nzchar(sequences))) {
    bad <- which(is.na(sequences) | !nzchar(sequences))[[1L]]
    rlang::abort(paste0(
      "`", sequence, "` contains a missing or empty sequence at row ",
      bad, "."
    ))
  }
  writeFastqas(
    data.frame(id = headers, seq = sequences, stringsAsFactors = FALSE),
    fn = path,
    type = "fa",
    header_col = "id"
  )
  normalizePath(path, mustWork = TRUE)
}
