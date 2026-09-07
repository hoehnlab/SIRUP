# Native RAD whitelist-scan facade ------------------------------------------

.sirap_rad_scan_wl_identity_fields <- c(
  "core_version", "core_api_version", "source_commit",
  "embedded_source_digest", "resource_schema", "resource_digest"
)

.sirap_validate_rad_scan_wl_core <- function(core) {
  .sirap_validate_rad_core(core)
  if (!is.character(core$features) || anyNA(core$features) ||
      !("scan-wl" %in% core$features)) {
    .sirap_rad_incompatible(
      "The embedded RAD core lacks the required scan-wl capability."
    )
  }
  invisible(core)
}

.sirap_rad_scan_wl_identity_mismatch <- function(reference, candidate) {
  fields <- .sirap_rad_scan_wl_identity_fields
  fields[vapply(fields, function(field) {
    !identical(
      as.character(reference[[field]]),
      as.character(candidate[[field]])
    )
  }, logical(1))]
}

.sirap_assert_rad_scan_wl_identity <- function(reference, candidate,
                                               context) {
  mismatch <- .sirap_rad_scan_wl_identity_mismatch(reference, candidate)
  if (length(mismatch)) {
    .sirap_rad_incompatible(paste0(
      "RAD changed ", context, " (", paste(mismatch, collapse = ", "), ")."
    ))
  }
  invisible(candidate)
}

.sirap_validate_native_rad_scan_wl_run <- function(result, expected_core) {
  required <- c(
    "status", "backend", "core", "input", "output_prefix", "config",
    "artifacts", "stats", "log"
  )
  if (!is.list(result) || !inherits(result, "rad_scan_wl_result") ||
      is.null(names(result))) {
    .sirap_rad_incompatible(
      "`rrad::rad_scan_wl()` did not return a `rad_scan_wl_result` object."
    )
  }
  missing_fields <- setdiff(required, names(result))
  if (length(missing_fields)) {
    .sirap_rad_incompatible(paste0(
      "The native RAD whitelist-scan result is missing: ",
      paste(missing_fields, collapse = ", "), "."
    ))
  }
  if (!identical(result$status, "success")) {
    rlang::abort(
      paste0(
        "The embedded RAD whitelist scan returned status '",
        .sirap_rad_value(result$status), "'."
      ),
      class = "sirap_rad_run_error"
    )
  }
  if (!identical(result$backend, "embedded")) {
    .sirap_rad_incompatible(paste0(
      "Expected rrad's embedded whitelist-scan backend, but the result ",
      "reported '", .sirap_rad_value(result$backend), "'."
    ))
  }
  if (!is.character(result$input) || length(result$input) != 1L ||
      is.na(result$input) || !nzchar(result$input) ||
      !is.character(result$output_prefix) ||
      length(result$output_prefix) != 1L || is.na(result$output_prefix) ||
      !nzchar(result$output_prefix) || !is.list(result$config) ||
      !is.list(result$artifacts) || !is.list(result$stats) ||
      !is.character(result$log) || length(result$log) != 1L ||
      is.na(result$log)) {
    .sirap_rad_incompatible(
      "The native RAD whitelist-scan result has invalid metadata."
    )
  }

  .sirap_validate_rad_scan_wl_core(result$core)
  .sirap_assert_rad_scan_wl_identity(
    expected_core,
    result$core,
    "between the preflight handshake and the whitelist-scan result"
  )
  invisible(result)
}

.sirap_validate_native_rad_scan_wl_batch <- function(result, expected_core) {
  required <- c(
    "status", "backend", "core", "batch_csv", "config", "runs"
  )
  if (!is.list(result) || !inherits(result, "rad_scan_wl_batch_result") ||
      is.null(names(result))) {
    .sirap_rad_incompatible(paste0(
      "`rrad::rad_scan_wl()` did not return a ",
      "`rad_scan_wl_batch_result` object."
    ))
  }
  missing_fields <- setdiff(required, names(result))
  if (length(missing_fields)) {
    .sirap_rad_incompatible(paste0(
      "The native RAD whitelist-scan batch result is missing: ",
      paste(missing_fields, collapse = ", "), "."
    ))
  }
  if (!identical(result$status, "success")) {
    rlang::abort(
      paste0(
        "The embedded RAD whitelist-scan batch returned status '",
        .sirap_rad_value(result$status), "'."
      ),
      class = "sirap_rad_run_error"
    )
  }
  if (!identical(result$backend, "embedded")) {
    .sirap_rad_incompatible(paste0(
      "Expected rrad's embedded whitelist-scan backend, but the batch ",
      "result reported '", .sirap_rad_value(result$backend), "'."
    ))
  }
  if (!is.character(result$batch_csv) || length(result$batch_csv) != 1L ||
      is.na(result$batch_csv) || !nzchar(result$batch_csv) ||
      !is.list(result$config) || !is.list(result$runs) ||
      !length(result$runs)) {
    .sirap_rad_incompatible(
      "The native RAD whitelist-scan batch result has invalid metadata."
    )
  }

  .sirap_validate_rad_scan_wl_core(result$core)
  .sirap_assert_rad_scan_wl_identity(
    expected_core,
    result$core,
    "between the preflight handshake and the whitelist-scan batch result"
  )
  for (run in result$runs) {
    .sirap_validate_native_rad_scan_wl_run(run, expected_core)
  }
  invisible(result)
}

.sirap_validate_native_rad_scan_wl_result <- function(result,
                                                      expected_core) {
  if (inherits(result, "rad_scan_wl_batch_result")) {
    return(.sirap_validate_native_rad_scan_wl_batch(result, expected_core))
  }
  if (inherits(result, "rad_scan_wl_result")) {
    return(.sirap_validate_native_rad_scan_wl_run(result, expected_core))
  }
  .sirap_rad_incompatible(paste0(
    "`rrad::rad_scan_wl()` did not return a `rad_scan_wl_result` or ",
    "`rad_scan_wl_batch_result` object."
  ))
}

.rrad_scan_wl <- function(...) {
  rrad::rad_scan_wl(...)
}

#' Discover a barcode whitelist with SIRAP's native RAD engine
#'
#' `rad_scan_wl()` is SIRAP's zero-setup facade for
#' [rrad::rad_scan_wl()]. It scans FASTQ or FASTA reads for barcodes after an
#' adapter, optionally restricts observations to a reference whitelist, and
#' writes RAD's selected whitelist and scan statistics. Statistical rescans,
#' two-part BC1/BC2 scans, and batch input use the same packaged native engine.
#'
#' SIRAP verifies the native core and its packaged resource identity before
#' and after every call, independently checks the core recorded in the result,
#' and adds SIRAP timing and version provenance. Large built-in reference
#' whitelists retain rrad's download-on-first-use, checksum verification, and
#' user-cache behavior; no standalone RAD executable is needed.
#'
#' @param input Input FASTQ, FASTA, or gzip-compressed FASTQ/FASTA path. With
#'   `rescan = TRUE`, a CSV, TSV, or whitespace-delimited barcode/count table
#'   instead. Supply exactly one of `input` and `batch_csv`.
#' @param output_prefix Prefix for generated files. The default is
#'   `"barcodes"` for sequence scans and `<input>_rescan` for statistical
#'   rescans.
#' @param batch_csv Optional two-column CSV containing an input sequence-file
#'   path and output prefix per row. Every row is validated by rrad before the
#'   first native scan starts.
#' @param adapter_seq Adapter or primer immediately upstream of the barcode.
#'   Required for sequence scans and unused for statistical rescans.
#' @param barcode_length Number of barcode bases to extract, from 1 through 32.
#' @param left_margin,right_margin Non-negative RAD extraction offsets.
#' @param max_reads Optional positive maximum number of reads to inspect;
#'   `NULL` scans the entire input.
#' @param max_error Maximum adapter edit-distance ratio, from 0 through 1.
#' @param whitelist Optional built-in whitelist key or whitelist file used to
#'   restrict and validate single-barcode observations.
#' @param threads Positive requested worker-thread count. A build without
#'   OpenMP warns and uses one effective thread.
#' @param chunk_size Positive number of reads per processing chunk.
#' @param selection Which barcode set is written: RAD's
#'   `"high_specificity"` final calls or every `"above_floor"` barcode.
#' @param rescan Re-run only RAD's statistical selector on barcode/count data;
#'   no sequence records are scanned.
#' @param bc1_whitelist,bc2_whitelist Optional built-in keys or files for a
#'   two-part barcode scan. Supply both or neither.
#' @param bc1_length,bc2_length Optional packed-whitelist barcode lengths,
#'   each from 1 through 32.
#' @param umi_length Non-negative UMI length between the adapter and BC1 in
#'   two-part mode.
#' @param offset_min,offset_max Inclusive BC1 search-offset range.
#' @param verbose Print native progress captured during each completed run.
#'
#' @return A result inheriting from `sirap_rad_scan_wl_result`,
#'   `sirap_rad_result`, and rrad's `rad_scan_wl_result`. Batch mode also adds
#'   `sirap_rad_scan_wl_batch_result` and retains
#'   `rad_scan_wl_batch_result`. The original rrad result is augmented with
#'   SIRAP version, timing, backend, and verified core-handshake provenance.
#' @seealso [rad_demux()], [radScanWl()], [rrad::rad_scan_wl()]
#' @export
rad_scan_wl <- function(input = NULL,
                        output_prefix = NULL,
                        batch_csv = NULL,
                        adapter_seq = NULL,
                        barcode_length = 16L,
                        left_margin = 0L,
                        right_margin = 0L,
                        max_reads = NULL,
                        max_error = 0.3,
                        whitelist = NULL,
                        threads = 1L,
                        chunk_size = 10000L,
                        selection = c("high_specificity", "above_floor"),
                        rescan = FALSE,
                        bc1_whitelist = NULL,
                        bc2_whitelist = NULL,
                        bc1_length = NULL,
                        bc2_length = NULL,
                        umi_length = 9L,
                        offset_min = 0L,
                        offset_max = 3L,
                        verbose = FALSE) {
  started <- Sys.time()
  preflight_core <- rad_core_info("native")
  .sirap_validate_rad_scan_wl_core(preflight_core)

  config <- list(
    input = input,
    output_prefix = output_prefix,
    batch_csv = batch_csv,
    adapter_seq = adapter_seq,
    barcode_length = barcode_length,
    left_margin = left_margin,
    right_margin = right_margin,
    max_reads = max_reads,
    max_error = max_error,
    whitelist = whitelist,
    threads = threads,
    chunk_size = chunk_size,
    selection = selection,
    rescan = rescan,
    bc1_whitelist = bc1_whitelist,
    bc2_whitelist = bc2_whitelist,
    bc1_length = bc1_length,
    bc2_length = bc2_length,
    umi_length = umi_length,
    offset_min = offset_min,
    offset_max = offset_max,
    verbose = verbose
  )
  result <- do.call(.rrad_scan_wl, config)

  postflight_core <- rad_core_info("native")
  .sirap_validate_rad_scan_wl_core(postflight_core)
  .sirap_validate_native_rad_scan_wl_result(result, preflight_core)
  .sirap_assert_rad_scan_wl_identity(
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
  if (inherits(result, "rad_scan_wl_batch_result")) {
    class(result) <- unique(c(
      "sirap_rad_scan_wl_batch_result",
      "sirap_rad_scan_wl_result",
      class(result)
    ))
  } else {
    class(result) <- unique(c("sirap_rad_scan_wl_result", class(result)))
  }
  result
}

#' @param x A SIRAP native RAD whitelist-scan result.
#' @param ... Additional arguments reserved for print-method compatibility.
#' @rdname rad_scan_wl
#' @export
print.sirap_rad_scan_wl_result <- function(x, ...) {
  core_version <- x$core$core_version %||% "unknown"
  api_version <- x$core$core_api_version %||% NA_integer_
  api_label <- if (length(api_version) == 1L && !is.na(api_version)) {
    paste0("core API ", api_version)
  } else {
    "unknown core API"
  }

  is_batch <- inherits(x, "rad_scan_wl_batch_result")
  label <- if (is_batch) {
    "<SIRAP RAD whitelist-scan batch result>"
  } else {
    "<SIRAP RAD whitelist-scan result>"
  }
  cat(label, "\n", sep = "")
  cat("  Status:  ", x$status %||% "unknown", "\n", sep = "")
  cat("  Backend: ", x$facade_backend %||% "native", "\n", sep = "")
  cat("  RAD:     ", core_version, " (", api_label, ")\n", sep = "")
  if (is_batch) {
    cat("  Batch:   ", x$batch_csv, "\n", sep = "")
    cat("  Runs:    ", length(x$runs), "\n", sep = "")
  } else {
    cat("  Input:   ", x$input, "\n", sep = "")
    cat("  Output:  ", x$output_prefix, "\n", sep = "")
  }
  invisible(x)
}
