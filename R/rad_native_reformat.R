# Native RAD reformat facade -------------------------------------------------

.sirap_rad_reformat_identity_fields <- c(
  "core_version", "core_api_version", "source_commit",
  "embedded_source_digest", "resource_schema", "resource_digest"
)

.sirap_validate_rad_reformat_core <- function(core, to_fasta = FALSE,
                                              header_format = NULL) {
  .sirap_validate_rad_core(core)
  required <- "reformat"
  if (isTRUE(to_fasta)) {
    required <- c(required, "reformat-fasta")
  }
  presto_requested <- is.character(header_format) &&
    length(header_format) == 1L && !is.na(header_format) &&
    identical(tolower(header_format), "presto")
  if (presto_requested) {
    required <- c(required, "reformat-presto-header")
  }
  if (!is.character(core$features) || anyNA(core$features)) {
    .sirap_rad_incompatible(
      "The embedded RAD core lacks the required reformat capability."
    )
  }
  missing <- setdiff(required, core$features)
  if (length(missing)) {
    .sirap_rad_incompatible(paste0(
      "The embedded RAD core lacks required reformat feature",
      if (length(missing) == 1L) " " else "s ",
      paste(missing, collapse = ", "), "."
    ))
  }
  invisible(core)
}

.sirap_rad_reformat_identity_mismatch <- function(reference, candidate) {
  fields <- .sirap_rad_reformat_identity_fields
  fields[vapply(fields, function(field) {
    !identical(
      as.character(reference[[field]]),
      as.character(candidate[[field]])
    )
  }, logical(1))]
}

.sirap_assert_rad_reformat_identity <- function(reference, candidate,
                                                context) {
  mismatch <- .sirap_rad_reformat_identity_mismatch(reference, candidate)
  if (length(mismatch)) {
    .sirap_rad_incompatible(paste0(
      "RAD changed ", context, " (", paste(mismatch, collapse = ", "), ")."
    ))
  }
  invisible(candidate)
}

.sirap_validate_native_rad_reformat_result <- function(result,
                                                       expected_core,
                                                       to_fasta = FALSE,
                                                       header_format = NULL) {
  required <- c(
    "status", "backend", "core", "input", "output", "config",
    "artifacts", "stats", "log"
  )
  if (!is.list(result) || !inherits(result, "rad_reformat_result") ||
      is.null(names(result))) {
    .sirap_rad_incompatible(
      "`rrad::rad_reformat()` did not return a `rad_reformat_result` object."
    )
  }
  missing_fields <- setdiff(required, names(result))
  if (length(missing_fields)) {
    .sirap_rad_incompatible(paste0(
      "The native RAD reformat result is missing: ",
      paste(missing_fields, collapse = ", "), "."
    ))
  }
  if (!identical(result$status, "success")) {
    rlang::abort(
      paste0("The embedded RAD reformat run returned status '",
             .sirap_rad_value(result$status), "'."),
      class = "sirap_rad_run_error"
    )
  }
  if (!identical(result$backend, "embedded")) {
    .sirap_rad_incompatible(paste0(
      "Expected rrad's embedded reformat backend, but the result reported '",
      .sirap_rad_value(result$backend), "'."
    ))
  }
  if (!is.character(result$input) || length(result$input) != 1L ||
      is.na(result$input) || !nzchar(result$input) ||
      !is.character(result$output) || length(result$output) != 1L ||
      is.na(result$output) || !nzchar(result$output) ||
      !is.list(result$config) || !is.list(result$artifacts) ||
      !is.list(result$stats) || !is.character(result$log) ||
      length(result$log) != 1L || is.na(result$log)) {
    .sirap_rad_incompatible(
      "The native RAD reformat result has invalid metadata."
    )
  }

  .sirap_validate_rad_reformat_core(
    result$core, to_fasta = to_fasta, header_format = header_format
  )
  .sirap_assert_rad_reformat_identity(
    expected_core,
    result$core,
    "between the preflight handshake and the reformat result"
  )
  invisible(result)
}

.rrad_reformat <- function(...) {
  rrad::rad_reformat(...)
}

#' Reformat RAD-tagged reads with SIRAP's native RAD engine
#'
#' `rad_reformat()` is SIRAP's native facade for [rrad::rad_reformat()]. It
#' collapses RAD barcode tags into read names, splits records into one file per
#' barcode, converts FASTQ to FASTA, emits canonical pRESTO/IgBLAST headers, or
#' translates supported spatial barcode pairs into coordinate headers without
#' launching an external RAD executable. SIRAP verifies the native engine
#' immediately before and after the run and records the exact engine identity
#' in the returned provenance.
#'
#' @param fastq Input FASTQ, FASTA, gzip-compressed FASTQ/FASTA path, or a
#'   `rad_demux_result` returned by [rad_demux()].
#' @param out_dir Output directory for barcode-specific files. Required when
#'   `split_bc = TRUE` and unused otherwise.
#' @param split_bc Write one gzip file per observed `CB:Z` barcode.
#' @param reformat_header Collapse the read name, cell barcode, and UMI to
#'   `QNAME<delimiter>CB<delimiter>UB`.
#' @param collapsed_input Parse tagless names as already-collapsed
#'   `QNAME<delimiter>CB<delimiter>UB` values.
#' @param delimiter Single-byte delimiter used for collapsed headers.
#' @param coordinate Optional coordinate mode. `NULL` disables coordinate
#'   rewriting; `"vizHD-v1"` maps supported Visium HD barcode pairs.
#' @param bin_size Positive spatial bin size in microns.
#' @param output_fastq Optional aggregate output path for a non-splitting run.
#'   When omitted, a header or coordinate rewrite replaces `fastq` atomically
#'   in place. This is also the FASTA destination when `to_fasta = TRUE`.
#' @param threads Positive requested worker-thread count.
#' @param chunk_size Positive number of records per streaming chunk.
#' @param verbose Print native progress captured during the completed run.
#' @param to_fasta Convert FASTQ to FASTA in the native streaming pass.
#' @param header_format Optional header contract: `"preserve"`, legacy
#'   `"collapsed"`, or canonical `"presto"`. The pRESTO form is
#'   `QNAME|BARCODE=value|UMI=UB` and is safe for the first IgBLAST pass.
#'   `BARCODE` uses a coordinate-reformatted `BC:Z` value when present and
#'   otherwise the corrected `CB:Z` value.
#'
#' @return A result inheriting from `sirap_rad_reformat_result`,
#'   `sirap_rad_result`, and `rad_reformat_result`. The original rrad result is
#'   retained and augmented with SIRAP version, timing, backend, and core
#'   handshake provenance.
#' @seealso [rad_demux()], [radReformat()], [rrad::rad_reformat()]
#' @export
rad_reformat <- function(fastq,
                         out_dir = NULL,
                         split_bc = FALSE,
                         reformat_header = FALSE,
                         collapsed_input = FALSE,
                         delimiter = "_",
                         coordinate = NULL,
                         bin_size = 2L,
                         output_fastq = NULL,
                         threads = 1L,
                         chunk_size = 5000L,
                         verbose = FALSE,
                         to_fasta = FALSE,
                         header_format = NULL) {
  started <- Sys.time()
  preflight_core <- rad_core_info("native")
  .sirap_validate_rad_reformat_core(
    preflight_core, to_fasta = to_fasta, header_format = header_format
  )

  result <- .rrad_reformat(
    fastq = fastq,
    out_dir = out_dir,
    split_bc = split_bc,
    reformat_header = reformat_header,
    collapsed_input = collapsed_input,
    delimiter = delimiter,
    coordinate = coordinate,
    bin_size = bin_size,
    output_fastq = output_fastq,
    threads = threads,
    chunk_size = chunk_size,
    verbose = verbose,
    to_fasta = to_fasta,
    header_format = header_format
  )

  postflight_core <- rad_core_info("native")
  .sirap_validate_rad_reformat_core(
    postflight_core, to_fasta = to_fasta, header_format = header_format
  )
  .sirap_validate_native_rad_reformat_result(
    result, preflight_core,
    to_fasta = to_fasta,
    header_format = header_format
  )
  .sirap_assert_rad_reformat_identity(
    preflight_core,
    postflight_core,
    "between the preflight and postflight handshakes"
  )

  result <- .sirap_finish_rad_result(
    result,
    backend = "native",
    core = preflight_core,
    started = started
  )
  class(result) <- unique(c("sirap_rad_reformat_result", class(result)))
  result
}

#' @export
print.sirap_rad_reformat_result <- function(x, ...) {
  core_version <- x$core$core_version %||% "unknown"
  api_version <- x$core$core_api_version %||% NA_integer_
  api_label <- if (length(api_version) == 1L && !is.na(api_version)) {
    paste0("core API ", api_version)
  } else {
    "unknown core API"
  }

  cat("<SIRAP RAD reformat result>\n", sep = "")
  cat("  Status:  ", x$status %||% "unknown", "\n", sep = "")
  cat("  Backend: ", x$facade_backend %||% "native", "\n", sep = "")
  cat("  RAD:     ", core_version, " (", api_label, ")\n", sep = "")
  if (is.character(x$input) && length(x$input) == 1L) {
    cat("  Input:   ", x$input, "\n", sep = "")
  }
  if (is.character(x$output) && length(x$output) == 1L) {
    cat("  Output:  ", x$output, "\n", sep = "")
  }
  invisible(x)
}
