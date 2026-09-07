# rad/rrad compatibility contract used by the SIRAP facade. These values are
# intentionally independent of package versions: a patch or feature release is
# compatible when it preserves the core API and resource schema.
.sirap_rad_core_api_version <- 1L
.sirap_rad_resource_schema <- 1L
.sirap_rad_required_features <- "demux"
.sirap_rad_supported_cli_versions <- "1.0.2"

.sirap_rad_configuration_error <- function(message) {
  rlang::abort(message, class = "sirap_rad_configuration_error")
}

.sirap_rad_scalar_character <- function(x, name, allow_null = FALSE) {
  if (allow_null && is.null(x)) {
    return(NULL)
  }
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    .sirap_rad_configuration_error(
      paste0("`", name, "` must be one non-empty character value.")
    )
  }
  enc2utf8(x)
}

.sirap_rad_scalar_logical <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    .sirap_rad_configuration_error(
      paste0("`", name, "` must be TRUE or FALSE.")
    )
  }
  x
}

.sirap_rad_scalar_integer <- function(x, name, minimum = 0L,
                                      maximum = .Machine$integer.max,
                                      allow_null = FALSE) {
  if (allow_null && is.null(x)) {
    return(NULL)
  }
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x != trunc(x) || x < minimum || x > maximum) {
    .sirap_rad_configuration_error(paste0(
      "`", name, "` must be one whole number between ", minimum,
      " and ", maximum, "."
    ))
  }
  as.integer(x)
}

.sirap_rad_scalar_number <- function(x, name, minimum = 0,
                                     maximum = Inf, allow_null = FALSE) {
  if (allow_null && is.null(x)) {
    return(NULL)
  }
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < minimum || x > maximum) {
    requirement <- if (is.finite(maximum)) {
      paste0("between ", minimum, " and ", maximum)
    } else {
      paste0("greater than or equal to ", minimum)
    }
    .sirap_rad_configuration_error(paste0(
      "`", name, "` must be one finite number ", requirement, "."
    ))
  }
  as.numeric(x)
}

.sirap_rad_existing_file <- function(x, name) {
  x <- .sirap_rad_scalar_character(x, name)
  if (!file.exists(x) || dir.exists(x)) {
    .sirap_rad_configuration_error(
      paste0("`", name, "` does not identify an existing file: ", x, ".")
    )
  }
  normalizePath(x, mustWork = TRUE)
}

.sirap_rad_optional_existing_file <- function(x, name) {
  if (is.null(x)) {
    return(NULL)
  }
  .sirap_rad_existing_file(x, name)
}

.sirap_rad_choice <- function(x, name, choices) {
  x <- tolower(.sirap_rad_scalar_character(x, name))
  if (!(x %in% choices)) {
    .sirap_rad_configuration_error(paste0(
      "`", name, "` must be one of ",
      paste0("\"", choices, "\"", collapse = ", "), "."
    ))
  }
  x
}

.sirap_normalize_rad_config <- function(config) {
  .sirap_validate_whitelist_config(config)

  layout <- .sirap_rad_scalar_character(config$layout, "layout")
  if (file.exists(layout)) {
    if (dir.exists(layout)) {
      .sirap_rad_configuration_error(
        "`layout` must identify a layout file or a built-in layout key."
      )
    }
    layout <- normalizePath(layout, mustWork = TRUE)
  } else if (grepl("[/\\\\]", layout) ||
             grepl("\\.(csv|ya?ml)$", layout, ignore.case = TRUE)) {
    .sirap_rad_configuration_error(
      paste0("`layout` file does not exist: ", layout, ".")
    )
  }

  fastq <- .sirap_rad_existing_file(config$fastq, "fastq")
  if (!grepl("\\.(fastq|fq)(\\.gz)?$", fastq, ignore.case = TRUE)) {
    .sirap_rad_configuration_error(
      "`fastq` must end in .fastq, .fq, .fastq.gz, or .fq.gz."
    )
  }

  out_dir <- .sirap_rad_scalar_character(config$out_dir, "out_dir")
  if (file.exists(out_dir) && !dir.exists(out_dir)) {
    .sirap_rad_configuration_error(
      paste0("`out_dir` identifies a file, not a directory: ", out_dir, ".")
    )
  }
  if (!dir.exists(out_dir) &&
      !dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)) {
    .sirap_rad_configuration_error(
      paste0("Could not create `out_dir`: ", out_dir, ".")
    )
  }
  out_dir <- normalizePath(out_dir, mustWork = TRUE)

  output <- .sirap_rad_scalar_character(config$output, "output")
  if (output %in% c(".", "..") || grepl("[/\\\\]", output)) {
    .sirap_rad_configuration_error(
      "`output` must be a filename prefix within `out_dir`."
    )
  }

  config$layout <- layout
  config$fastq <- fastq
  config$out_dir <- out_dir
  config$output <- output
  config["kit"] <- list(
    .sirap_rad_scalar_character(config$kit, "kit", allow_null = TRUE)
  )
  if (!is.null(config$kit) && grepl(":", config$kit, fixed = TRUE)) {
    .sirap_rad_configuration_error(
      "`kit` must be one built-in kit key without ':'."
    )
  }
  config["global_whitelist"] <- list(
    .sirap_rad_optional_existing_file(
      config$global_whitelist, "global_whitelist"
    )
  )
  config["custom_whitelist"] <- list(
    .sirap_rad_optional_existing_file(
      config$custom_whitelist, "custom_whitelist"
    )
  )
  config$threads <- .sirap_rad_scalar_integer(
    config$threads, "threads", minimum = 1L
  )
  config$chunk_size <- .sirap_rad_scalar_integer(
    config$chunk_size, "chunk_size", minimum = 1L
  )
  config["max_reads"] <- list(
    .sirap_rad_scalar_integer(
      config$max_reads, "max_reads", minimum = 1L, allow_null = TRUE
    )
  )
  for (name in c(
    "auto_whitelist", "rc_umi", "write_debug", "verbose", "reuse_cache"
  )) {
    config[[name]] <- .sirap_rad_scalar_logical(config[[name]], name)
  }
  config["scan_adapter"] <- list(
    .sirap_rad_scalar_character(
      config$scan_adapter, "scan_adapter", allow_null = TRUE
    )
  )
  config["scan_barcode_length"] <- list(
    .sirap_rad_scalar_integer(
      config$scan_barcode_length, "scan_barcode_length",
      minimum = 1L, maximum = 32L, allow_null = TRUE
    )
  )
  config$scan_max_error <- .sirap_rad_scalar_number(
    config$scan_max_error, "scan_max_error", minimum = 0, maximum = 1
  )
  config["scan_max_reads"] <- list(
    .sirap_rad_scalar_integer(
      config$scan_max_reads, "scan_max_reads",
      minimum = 1L, allow_null = TRUE
    )
  )
  config["scan_chunk_size"] <- list(
    .sirap_rad_scalar_integer(
      config$scan_chunk_size, "scan_chunk_size",
      minimum = 1L, allow_null = TRUE
    )
  )
  config["scan_threads"] <- list(
    .sirap_rad_scalar_integer(
      config$scan_threads, "scan_threads",
      minimum = 1L, allow_null = TRUE
    )
  )
  config$scan_selection <- .sirap_rad_choice(
    config$scan_selection,
    "scan_selection",
    c("high_specificity", "above_floor")
  )
  config["whitelist_mutation"] <- list(
    .sirap_rad_scalar_integer(
      config$whitelist_mutation, "whitelist_mutation",
      minimum = 0L, allow_null = TRUE
    )
  )
  config["generated_mutation"] <- list(
    .sirap_rad_scalar_integer(
      config$generated_mutation, "generated_mutation",
      minimum = 0L, allow_null = TRUE
    )
  )
  config$bc_correction_mode <- .sirap_rad_choice(
    config$bc_correction_mode, "bc_correction_mode", c("offensive", "defensive")
  )
  config$joint_bc_mode <- .sirap_rad_choice(
    config$joint_bc_mode, "joint_bc_mode", c("default", "strict")
  )
  config["min_read_length"] <- list(
    .sirap_rad_scalar_integer(
      config$min_read_length, "min_read_length", minimum = 0L,
      allow_null = TRUE
    )
  )
  config
}

#' Inspect the RAD engine available to SIRAP
#'
#' `rad_core_info()` performs the same compatibility handshake used by
#' [rad_demux()]. The default `"native"` backend reports the RAD core compiled
#' into the required `rrad` dependency. The explicit `"cli"` backend reports
#' the external executable selected by [radPath()]. That compatibility adapter
#' requires RAD 1.0.2; legacy CLI builds provide only a version/binary-digest
#' handshake and are marked as such.
#'
#' @param backend Execution backend. `"native"` uses the RAD core bundled with
#'   `rrad`; `"cli"` uses the legacy external executable.
#' @param rad Optional path to an external RAD executable, used only when
#'   `backend = "cli"`.
#'
#' @return A named list inheriting from `sirap_rad_core_info`, including the
#'   core version/API, upstream source commit, patched embedded-source digest,
#'   resource schema/digest, capabilities, and build metadata.
#' @export
rad_core_info <- function(backend = c("native", "cli"), rad = NULL) {
  backend <- match.arg(backend)

  if (identical(backend, "native")) {
    info <- .rrad_core_info()
    .sirap_validate_rad_core(info)
    return(structure(info, class = unique(c("sirap_rad_core_info", class(info)))))
  }

  info <- .sirap_cli_core_info(rad)
  .sirap_validate_cli_core(info)
  structure(info, class = c("sirap_rad_core_info", "list"))
}

#' Demultiplex reads with the RAD engine bundled for R
#'
#' `rad_demux()` is SIRAP's zero-setup demultiplexing entry point. By default it
#' delegates to the native RAD core packaged by `rrad`, after checking the core
#' API, resource schema, and required capabilities. It returns RAD's structured
#' run result with a SIRAP provenance record attached.
#' When a named kit needs a large whitelist that is not bundled in the package,
#' `rrad` downloads it automatically on first use, verifies its pinned size and
#' SHA-256 digest, and reuses the per-user cached copy on later runs. Use
#' [rad_download_whitelists()] only to prepare an offline or shared cache.
#'
#' Set `backend = "cli"` only to use an existing external RAD executable. This
#' compatibility path delegates to [radDemux()] and normalizes its result into
#' the same SIRAP result class, but an older CLI can provide only a version-level
#' handshake. `kit` and `global_whitelist` are alternative global/reference
#' barcode sources and cannot be supplied together. Likewise,
#' `auto_whitelist = TRUE` cannot be combined with `custom_whitelist`, because
#' both select the true/sample barcode set. Current CLI semantics additionally
#' require `global_whitelist` and `custom_whitelist` together unless the global
#' whitelist is the reference for automatic discovery; the native backend also
#' supports either file on its own. For CLI runs, SIRAP reports only recognized
#' outputs created or updated by that invocation and treats a status-zero run
#' without a valid new primary FASTQ as a failure. It validates the complete
#' input and output FASTQ streams because legacy RAD 1.0.2 does not reliably
#' surface parser or gzip-writer failures. Legacy CLI debug mode is rejected
#' because that build can suppress primary output when debugging is enabled.
#'
#' @param layout Built-in RAD layout key or path to a custom layout file.
#' @param fastq Input FASTQ or gzip-compressed FASTQ path. The CLI fallback
#'   performs a complete streaming structural/integrity pass before execution.
#' @param out_dir Directory in which RAD should write its artifacts.
#' @param output Output filename prefix.
#' @param kit Optional built-in whitelist kit.
#' @param global_whitelist Optional global whitelist CSV path.
#' @param custom_whitelist Optional true-cell/custom whitelist CSV path.
#' @param threads Requested number of worker threads. Current macOS builds run
#'   serially; `rrad` warns and records requested/effective counts when it must
#'   reduce this value to one.
#' @param chunk_size Number of reads processed per chunk.
#' @param max_reads Optional maximum number of reads to process.
#' @param auto_whitelist Whether to discover the sample whitelist before
#'   demultiplexing. This currently supports single-barcode layouts.
#' @param scan_adapter Optional adapter sequence used for automatic whitelist
#'   discovery. By default RAD derives it from the selected layout.
#' @param scan_barcode_length Optional barcode length from 1 to 32 for
#'   automatic whitelist discovery. By default RAD derives it from the layout.
#' @param scan_max_error Maximum adapter edit-distance ratio from 0 to 1 used
#'   during automatic whitelist discovery.
#' @param scan_max_reads Optional positive scan-specific read limit. `NULL`
#'   inherits `max_reads`.
#' @param scan_chunk_size Optional positive scan-specific chunk size. `NULL`
#'   inherits `chunk_size`.
#' @param scan_threads Optional positive scan-specific thread count. `NULL`
#'   inherits `threads`.
#' @param scan_selection Barcode-selection policy for automatic discovery:
#'   `"high_specificity"` or `"above_floor"`.
#' @param whitelist_mutation Edit-distance search space for whitelist entries.
#' @param generated_mutation Edit-distance search space for generated entries.
#' @param bc_correction_mode Barcode correction mode, either `"offensive"` or
#'   `"defensive"`.
#' @param joint_bc_mode Joint barcode policy, either `"default"` or `"strict"`.
#' @param rc_umi Whether reverse-strand UMIs should be reverse-complemented.
#' @param min_read_length Optional minimum retained cDNA length.
#' @param write_debug Whether to write RAD diagnostic artifacts in addition to
#'   primary output. Supported by the native backend; rejected for the legacy
#'   CLI backend because its debug writer can suppress primary data.
#' @param verbose Whether to enable verbose RAD logging.
#' @param reuse_cache Whether the native backend may reuse existing layout and
#'   position-map files at the same output prefix. The safe default, `FALSE`,
#'   recalibrates them for the requested input. The legacy CLI backend retains
#'   the external command's own cache behavior, which is reported as
#'   `"legacy-cli-uncontrolled"` rather than as the requested value.
#' @param backend Execution backend. The default `"native"` requires no
#'   separately installed RAD executable; `"cli"` explicitly selects the
#'   legacy [radDemux()] adapter.
#' @param rad Optional external RAD executable path for `backend = "cli"`.
#' @param capture Whether the CLI fallback should capture its console output.
#'
#' @return A structured result inheriting from `sirap_rad_result`. Native
#'   results also retain the `rad_demux_result` class supplied by `rrad`.
#' @seealso [rad_core_info()], [rad_download_whitelists()], [radDemux()]
#' @export
rad_demux <- function(layout,
                      fastq,
                      out_dir = ".",
                      output = "demux",
                      kit = NULL,
                      global_whitelist = NULL,
                      custom_whitelist = NULL,
                      threads = 1L,
                      chunk_size = 5000L,
                      max_reads = NULL,
                      auto_whitelist = FALSE,
                      scan_adapter = NULL,
                      scan_barcode_length = NULL,
                      scan_max_error = 0.3,
                      scan_max_reads = NULL,
                      scan_chunk_size = NULL,
                      scan_threads = NULL,
                      scan_selection = c("high_specificity", "above_floor"),
                      whitelist_mutation = NULL,
                      generated_mutation = NULL,
                      bc_correction_mode = "offensive",
                      joint_bc_mode = "default",
                      rc_umi = TRUE,
                      min_read_length = NULL,
                      write_debug = FALSE,
                      verbose = FALSE,
                      reuse_cache = FALSE,
                      backend = c("native", "cli"),
                      rad = NULL,
                      capture = TRUE) {
  backend <- match.arg(backend)
  if (missing(scan_selection)) {
    scan_selection <- "high_specificity"
  }
  started <- Sys.time()

  config <- list(
    layout = layout,
    fastq = fastq,
    out_dir = out_dir,
    output = output,
    kit = kit,
    global_whitelist = global_whitelist,
    custom_whitelist = custom_whitelist,
    threads = threads,
    chunk_size = chunk_size,
    max_reads = max_reads,
    auto_whitelist = auto_whitelist,
    scan_adapter = scan_adapter,
    scan_barcode_length = scan_barcode_length,
    scan_max_error = scan_max_error,
    scan_max_reads = scan_max_reads,
    scan_chunk_size = scan_chunk_size,
    scan_threads = scan_threads,
    scan_selection = scan_selection,
    whitelist_mutation = whitelist_mutation,
    generated_mutation = generated_mutation,
    bc_correction_mode = bc_correction_mode,
    joint_bc_mode = joint_bc_mode,
    rc_umi = rc_umi,
    min_read_length = min_read_length,
    write_debug = write_debug,
    verbose = verbose,
    reuse_cache = reuse_cache
  )
  config <- .sirap_normalize_rad_config(config)
  capture <- .sirap_rad_scalar_logical(capture, "capture")

  if (identical(backend, "native")) {
    core <- rad_core_info("native")
    .sirap_validate_rad_demux_features(core, config)
    result <- do.call(.rrad_demux, config)
    .sirap_validate_native_rad_result(result, expected_core = core)
    completed_core <- rad_core_info("native")
    .sirap_assert_native_rad_identity(
      core,
      completed_core,
      "between the preflight and postflight handshakes"
    )
    return(.sirap_finish_rad_result(result, "native", core, started))
  }

  .sirap_validate_cli_whitelist_config(config)
  core <- rad_core_info("cli", rad = rad)
  .sirap_validate_fastq_file(
    config$fastq,
    label = "legacy CLI input FASTQ",
    error_class = "sirap_rad_configuration_error"
  )
  before <- .sirap_snapshot_cli_outputs(config)
  cli_result <- .sirap_cli_demux(config, rad = rad, capture = capture)
  completed_core <- rad_core_info("cli", rad = rad)
  .sirap_validate_cli_identity(core, completed_core)
  result <- .sirap_normalize_cli_result(cli_result, config, core, before)
  .sirap_finish_rad_result(result, "cli", core, started)
}

#' @export
print.sirap_rad_result <- function(x, ...) {
  core_version <- x$core$core_version %||% "unknown"
  api_version <- x$core$core_api_version %||% NA_integer_
  facade_backend <- x$facade_backend %||% x$provenance$backend %||%
    x$backend %||% "unknown"
  api_label <- if (length(api_version) == 1L && !is.na(api_version)) {
    paste0("core API ", api_version)
  } else {
    "version-only handshake"
  }

  cat("<SIRAP RAD demultiplexing result>\n", sep = "")
  cat("  Status:  ", x$status %||% "unknown", "\n", sep = "")
  cat("  Backend: ", facade_backend, "\n", sep = "")
  if (!is.null(x$backend) && !identical(x$backend, facade_backend)) {
    cat("  Engine:  ", x$backend, "\n", sep = "")
  }
  cat("  RAD:     ", core_version, " (", api_label, ")\n", sep = "")
  if (!is.null(x$input$fastq)) {
    cat("  Input:   ", x$input$fastq, "\n", sep = "")
  }
  invisible(x)
}

.rrad_core_info <- function() {
  rrad::rad_core_info()
}

.rrad_demux <- function(...) {
  rrad::rad_demux(...)
}

.sirap_validate_rad_demux_features <- function(core, config) {
  required <- "demux"
  if (isTRUE(config$auto_whitelist)) {
    required <- c(required, "auto-whitelist", "scan-wl")
  }
  if (isTRUE(config$write_debug)) {
    required <- c(required, "debug-output")
  }
  missing_features <- setdiff(required, core$features)
  if (length(missing_features)) {
    .sirap_rad_incompatible(paste0(
      "The embedded RAD core lacks capabilities required by this run: ",
      paste(missing_features, collapse = ", "), "."
    ))
  }
  invisible(core)
}

.sirap_validate_rad_core <- function(info) {
  required <- c(
    "core_version", "core_api_version", "source_commit",
    "embedded_source_digest", "resource_schema", "resource_digest",
    "features", "build", "vendored"
  )
  if (!is.list(info) || is.null(names(info))) {
    .sirap_rad_incompatible("`rrad::rad_core_info()` did not return a named list.")
  }
  missing_fields <- setdiff(required, names(info))
  if (length(missing_fields)) {
    .sirap_rad_incompatible(paste0(
      "The embedded RAD handshake is missing: ",
      paste(missing_fields, collapse = ", "), "."
    ))
  }

  identity_fields <- c(
    "core_version", "source_commit", "embedded_source_digest",
    "resource_digest"
  )
  invalid_identity <- vapply(identity_fields, function(field) {
    value <- info[[field]]
    !is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)
  }, logical(1))
  if (any(invalid_identity)) {
    .sirap_rad_incompatible(paste0(
      "The embedded RAD handshake has an invalid ",
      paste(identity_fields[invalid_identity], collapse = ", "), "."
    ))
  }
  if (!grepl("^[[:xdigit:]]{64}$", info$resource_digest)) {
    .sirap_rad_incompatible(
      "The embedded RAD handshake has an invalid resource SHA-256 digest."
    )
  }
  if (!grepl("^[[:xdigit:]]{64}$", info$embedded_source_digest)) {
    .sirap_rad_incompatible(
      "The embedded RAD handshake has an invalid embedded-source SHA-256 digest."
    )
  }
  if (!grepl("^[[:xdigit:]]{40}$", info$source_commit)) {
    .sirap_rad_incompatible(
      "The embedded RAD handshake has an invalid 40-character source commit."
    )
  }
  if (!is.list(info$build)) {
    .sirap_rad_incompatible("The embedded RAD handshake has invalid build metadata.")
  }
  if (!is.logical(info$build$resources_available) ||
      length(info$build$resources_available) != 1L ||
      !isTRUE(info$build$resources_available) ||
      !is.character(info$build$resource_dir) ||
      length(info$build$resource_dir) != 1L ||
      is.na(info$build$resource_dir) || !nzchar(info$build$resource_dir)) {
    .sirap_rad_incompatible(
      "The embedded RAD handshake reports unavailable packaged resources."
    )
  }
  if (!isTRUE(info$vendored)) {
    .sirap_rad_incompatible(
      "The selected rrad build does not report a vendored RAD core."
    )
  }

  api <- info$core_api_version
  if (!is.numeric(api) || length(api) != 1L || is.na(api) ||
      !is.finite(api) || api != trunc(api) ||
      as.integer(api) != .sirap_rad_core_api_version) {
    .sirap_rad_incompatible(paste0(
      "SIRAP requires RAD core API ", .sirap_rad_core_api_version,
      ", but rrad reported ", .sirap_rad_value(info$core_api_version), "."
    ))
  }

  schema <- info$resource_schema
  if (!is.numeric(schema) || length(schema) != 1L || is.na(schema) ||
      !is.finite(schema) || schema != trunc(schema) ||
      as.integer(schema) != .sirap_rad_resource_schema) {
    .sirap_rad_incompatible(paste0(
      "SIRAP requires RAD resource schema ", .sirap_rad_resource_schema,
      ", but rrad reported ", .sirap_rad_value(info$resource_schema), "."
    ))
  }

  if (!is.character(info$features) || anyNA(info$features)) {
    .sirap_rad_incompatible("The embedded RAD handshake has an invalid capability list.")
  }
  features <- info$features
  missing_features <- setdiff(.sirap_rad_required_features, features)
  if (length(missing_features)) {
    .sirap_rad_incompatible(paste0(
      "The embedded RAD core lacks required capability: ",
      paste(missing_features, collapse = ", "), "."
    ))
  }

  invisible(info)
}

.sirap_assert_native_rad_identity <- function(reference, candidate, context) {
  .sirap_validate_rad_core(candidate)
  keys <- c(
    "core_version", "core_api_version", "source_commit",
    "embedded_source_digest", "resource_schema", "resource_digest"
  )
  mismatch <- vapply(
    keys,
    function(key) {
      !identical(as.character(candidate[[key]]), as.character(reference[[key]]))
    },
    logical(1)
  )
  if (any(mismatch)) {
    .sirap_rad_incompatible(paste0(
      "RAD changed ", context, " (",
      paste(keys[mismatch], collapse = ", "), ")."
    ))
  }
  invisible(candidate)
}

.sirap_validate_native_rad_result <- function(result, expected_core) {
  required <- c(
    "status", "backend", "core", "input", "config", "artifacts", "stats",
    "scan", "log"
  )
  if (!is.list(result) || !inherits(result, "rad_demux_result") ||
      is.null(names(result))) {
    .sirap_rad_incompatible(
      "`rrad::rad_demux()` did not return a `rad_demux_result` object."
    )
  }
  missing_fields <- setdiff(required, names(result))
  if (length(missing_fields)) {
    .sirap_rad_incompatible(paste0(
      "The native RAD result is missing: ",
      paste(missing_fields, collapse = ", "), "."
    ))
  }
  if (!identical(result$status, "success")) {
    rlang::abort(
      paste0("The embedded RAD run returned status '", result$status, "'."),
      class = "sirap_rad_run_error"
    )
  }
  if (!identical(result$backend, "embedded")) {
    .sirap_rad_incompatible(paste0(
      "Expected rrad's embedded backend, but the result reported '",
      result$backend, "'."
    ))
  }

  valid_string <- function(value) {
    is.character(value) && length(value) == 1L && !is.na(value) &&
      nzchar(value)
  }
  if (!is.list(result$input) ||
      !valid_string(result$input$fastq) ||
      !valid_string(result$input$layout) ||
      !is.list(result$config) || is.null(names(result$config)) ||
      !is.logical(result$config$auto_whitelist) ||
      length(result$config$auto_whitelist) != 1L ||
      is.na(result$config$auto_whitelist) ||
      !is.list(result$artifacts) || is.null(names(result$artifacts)) ||
      !is.list(result$stats) || is.null(names(result$stats)) ||
      !is.character(result$log) || length(result$log) != 1L ||
      is.na(result$log)) {
    .sirap_rad_incompatible(
      "The native RAD demultiplexing result has invalid metadata."
    )
  }

  .sirap_assert_native_rad_identity(
    expected_core,
    result$core,
    "between the preflight handshake and the completed run"
  )
  if (isTRUE(result$config$auto_whitelist)) {
    if (is.null(result$scan)) {
      .sirap_rad_incompatible(
        "The native RAD result is missing its automatic whitelist scan."
      )
    }
    .sirap_validate_native_rad_scan_wl_run(result$scan, expected_core)
  } else if (!is.null(result$scan)) {
    .sirap_rad_incompatible(
      "The native RAD result contains an unexpected automatic whitelist scan."
    )
  }
  invisible(result)
}

.sirap_cli_demux <- function(config, rad, capture) {
  rad <- if (is.null(rad)) radPath() else radPath(rad)
  extra <- c(
    if (isTRUE(config$auto_whitelist)) "--auto-wl" else character(),
    if (!isTRUE(config$rc_umi)) "--no-umi-rc" else character(),
    .rad_flag("--min-read-length", config$min_read_length)
  )

  radDemux(
    layout = config$layout,
    fastq = config$fastq,
    kit = config$kit,
    global_whitelist = config$global_whitelist,
    custom_whitelist = config$custom_whitelist,
    bc_correction_mode = config$bc_correction_mode,
    joint_bc_mode = config$joint_bc_mode,
    whitelist_mutation = config$whitelist_mutation,
    generated_mutation = config$generated_mutation,
    max_reads = config$max_reads,
    chunk_size = config$chunk_size,
    output = config$output,
    out_dir = config$out_dir,
    write_dbg = config$write_debug,
    threads = config$threads,
    verbose = config$verbose,
    extra,
    capture = capture,
    rad = rad
  )
}

.sirap_validate_cli_whitelist_config <- function(config) {
  .sirap_validate_whitelist_config(config)

  has_kit <- !is.null(config$kit)
  has_global <- !is.null(config$global_whitelist)
  has_custom <- !is.null(config$custom_whitelist)
  auto_whitelist <- isTRUE(config$auto_whitelist)

  if (has_custom && !has_kit && !has_global) {
    rlang::abort(
      paste0(
        "The legacy RAD CLI ignores an unpaired custom whitelist. ",
        "Supply `kit` or `global_whitelist`, or use `backend = \"native\"`."
      ),
      class = "sirap_rad_cli_compatibility_error"
    )
  }
  if (has_global && !has_custom && !has_kit && !auto_whitelist) {
    rlang::abort(
      paste0(
        "The legacy RAD CLI ignores an unpaired global whitelist. ",
        "Supply `custom_whitelist`, enable `auto_whitelist`, or use ",
        "`backend = \"native\"`."
      ),
      class = "sirap_rad_cli_compatibility_error"
    )
  }
  if (isTRUE(config$write_debug)) {
    rlang::abort(
      paste0(
        "The legacy RAD CLI's debug mode can suppress the primary FASTQ. ",
        "Use the default native backend for `write_debug = TRUE`."
      ),
      class = "sirap_rad_cli_compatibility_error"
    )
  }
  scan_defaults <- list(
    scan_adapter = NULL,
    scan_barcode_length = NULL,
    scan_max_error = 0.3,
    scan_max_reads = NULL,
    scan_chunk_size = NULL,
    scan_threads = NULL,
    scan_selection = "high_specificity"
  )
  changed_scan_options <- names(scan_defaults)[vapply(
    names(scan_defaults),
    function(name) {
      if (!(name %in% names(config))) {
        return(FALSE)
      }
      value <- config[[name]]
      default <- scan_defaults[[name]]
      if (is.null(value) && is.null(default)) {
        return(FALSE)
      }
      !identical(value, default)
    },
    logical(1)
  )]
  if (length(changed_scan_options)) {
    rlang::abort(
      paste0(
        "The legacy RAD CLI adapter cannot guarantee the requested automatic-",
        "whitelist scan controls: ",
        paste0("`", changed_scan_options, "`", collapse = ", "),
        ". Use the default native backend."
      ),
      class = "sirap_rad_cli_compatibility_error"
    )
  }
  invisible(config)
}

.sirap_validate_fastq_file <- function(path, label, error_class) {
  connection <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
  closed <- FALSE
  on.exit({
    if (!closed) {
      try(suppressWarnings(close(connection)), silent = TRUE)
    }
  }, add = TRUE)

  fail <- function(message) {
    rlang::abort(
      paste0("Invalid ", label, " at ", path, ": ", message),
      class = error_class
    )
  }
  record_offset <- 0L
  repeat {
    lines <- tryCatch(
      withCallingHandlers(
        readLines(connection, n = 32768L, warn = FALSE),
        warning = function(w) fail(conditionMessage(w))
      ),
      error = function(e) {
        if (inherits(e, error_class)) {
          stop(e)
        }
        fail(conditionMessage(e))
      }
    )
    if (!length(lines)) {
      break
    }
    if (length(lines) %% 4L != 0L) {
      fail(paste0(
        "the file ends with an incomplete four-line record near record ",
        record_offset + (length(lines) %/% 4L) + 1L, "."
      ))
    }

    header_idx <- seq.int(1L, length(lines), by = 4L)
    sequence_idx <- header_idx + 1L
    plus_idx <- header_idx + 2L
    quality_idx <- header_idx + 3L
    bad_header <- which(substr(lines[header_idx], 1L, 1L) != "@")[1L]
    if (!is.na(bad_header)) {
      fail(paste0("record ", record_offset + bad_header,
                  " does not begin with '@'."))
    }
    bad_plus <- which(substr(lines[plus_idx], 1L, 1L) != "+")[1L]
    if (!is.na(bad_plus)) {
      fail(paste0("record ", record_offset + bad_plus,
                  " has no '+' separator."))
    }
    sequence_lengths <- nchar(lines[sequence_idx], type = "bytes")
    quality_lengths <- nchar(lines[quality_idx], type = "bytes")
    bad_length <- which(sequence_lengths != quality_lengths)[1L]
    if (!is.na(bad_length)) {
      fail(paste0(
        "record ", record_offset + bad_length,
        " has different sequence and quality lengths."
      ))
    }
    record_offset <- record_offset + length(header_idx)
  }

  close_result <- tryCatch(
    withCallingHandlers(
      close(connection),
      warning = function(w) fail(conditionMessage(w))
    ),
    error = function(e) {
      if (inherits(e, error_class)) {
        stop(e)
      }
      fail(conditionMessage(e))
    }
  )
  closed <- TRUE
  invisible(close_result)
}

.sirap_validate_whitelist_config <- function(config) {
  has_kit <- !is.null(config$kit)
  has_global <- !is.null(config$global_whitelist)
  has_custom <- !is.null(config$custom_whitelist)
  auto_whitelist <- isTRUE(config$auto_whitelist)

  if (has_kit && has_global) {
    rlang::abort(
      paste0(
        "`kit` and `global_whitelist` are alternative global/reference ",
        "barcode sources; supply only one."
      ),
      class = "sirap_rad_configuration_error"
    )
  }
  if (auto_whitelist && has_custom) {
    rlang::abort(
      paste0(
        "`auto_whitelist = TRUE` and `custom_whitelist` both select the ",
        "true/sample barcode set; supply only one."
      ),
      class = "sirap_rad_configuration_error"
    )
  }
  invisible(config)
}

.sirap_cli_artifact_paths <- function(config) {
  output_dir <- normalizePath(config$out_dir, mustWork = FALSE)
  output_prefix <- file.path(output_dir, config$output)
  paths <- c(
    fastq = paste0(output_prefix, ".fq.gz"),
    layout = paste0(output_prefix, "_layout.csv"),
    position_map = paste0(output_prefix, "_position_map.csv"),
    whitelist_true = paste0(output_prefix, "_whitelist_true.csv"),
    whitelist_global = paste0(output_prefix, "_whitelist_global.csv"),
    demux_log = paste0(output_prefix, "_demux.log")
  )
  if (isTRUE(config$auto_whitelist)) {
    paths <- c(
      paths,
      scan_whitelist_csv = paste0(output_prefix, "_scanwl.csv"),
      scan_whitelist_txt = paste0(output_prefix, "_scanwl.txt"),
      scan_whitelist_log = paste0(output_prefix, "_scan_wl.log")
    )
  }
  if (isTRUE(config$write_debug)) {
    paths <- c(
      paths,
      debug_signature = paste0(output_prefix, "_dbg.sig.gz"),
      debug_csv = paste0(output_prefix, "_dbg.csv.gz"),
      debug_fastq = paste0(output_prefix, "_dbg.fq.gz"),
      metrics = paste0(output_prefix, ".metrics.tsv")
    )
  }
  paths
}

.sirap_snapshot_cli_outputs <- function(config) {
  paths <- .sirap_cli_artifact_paths(config)
  signatures <- lapply(paths, function(path) {
    info <- file.info(path)
    exists <- !is.na(info$size[[1L]]) && !isTRUE(info$isdir[[1L]])
    list(
      path = path,
      exists = exists,
      size = if (exists) unname(info$size[[1L]]) else NA_real_,
      mtime = if (exists) as.numeric(info$mtime[[1L]]) else NA_real_,
      ctime = if (exists) as.numeric(info$ctime[[1L]]) else NA_real_
    )
  })
  names(signatures) <- names(paths)
  signatures
}

.sirap_cli_output_changed <- function(before, after) {
  if (!isTRUE(after$exists)) {
    return(FALSE)
  }
  if (is.null(before) || !isTRUE(before$exists)) {
    return(TRUE)
  }
  !identical(before$size, after$size) ||
    !identical(before$mtime, after$mtime) ||
    !identical(before$ctime, after$ctime)
}

.sirap_normalize_cli_result <- function(cli_result, config, core, before) {
  raw_status <- if (is.list(cli_result)) cli_result$status else NULL
  valid_status <- is.numeric(raw_status) && length(raw_status) == 1L &&
    !is.na(raw_status) && is.finite(raw_status) &&
    raw_status == trunc(raw_status) && raw_status >= 0 &&
    raw_status <= .Machine$integer.max
  status_code <- if (valid_status) as.integer(raw_status) else NA_integer_
  if (!valid_status || status_code != 0L) {
    status <- if (is.list(cli_result)) .sirap_rad_value(cli_result$status) else "<invalid>"
    rlang::abort(
      paste0("The external RAD CLI returned status ", status, "."),
      class = "sirap_rad_run_error"
    )
  }
  output_prefix <- file.path(
    normalizePath(config$out_dir, mustWork = FALSE),
    config$output
  )
  after <- .sirap_snapshot_cli_outputs(config)
  changed <- vapply(names(after), function(name) {
    .sirap_cli_output_changed(before[[name]], after[[name]])
  }, logical(1))

  if (!isTRUE(changed[["fastq"]]) || after$fastq$size <= 0) {
    rlang::abort(
      paste0(
        "The external RAD CLI reported success but did not create or update ",
        "the required primary FASTQ output: ", after$fastq$path, "."
      ),
      class = "sirap_rad_run_error"
    )
  }
  .sirap_validate_fastq_file(
    after$fastq$path,
    label = "RAD CLI primary output",
    error_class = "sirap_rad_run_error"
  )
  artifact_paths <- vapply(after[changed], `[[`, character(1), "path")
  artifact_paths <- stats::setNames(
    normalizePath(unname(artifact_paths), mustWork = FALSE),
    names(artifact_paths)
  )
  artifacts <- c(
    list(output_prefix = output_prefix),
    as.list(artifact_paths),
    list(files = unname(artifact_paths))
  )

  effective_config <- config
  effective_config$reuse_cache_requested <- config$reuse_cache
  effective_config$reuse_cache <- NA
  effective_config$cache_policy <- "legacy-cli-uncontrolled"

  structure(
    list(
      status = "success",
      backend = "cli",
      core = unclass(core),
      input = list(
        fastq = normalizePath(config$fastq, mustWork = FALSE),
        layout = config$layout
      ),
      config = effective_config,
      artifacts = artifacts,
      stats = list(),
      log = list(
        status = status_code,
        stdout = cli_result$stdout,
        stderr = cli_result$stderr,
        command = cli_result$command
      )
    ),
    class = c("sirap_rad_cli_result", "list")
  )
}

.sirap_finish_rad_result <- function(result, backend, core, started) {
  finished <- Sys.time()
  engine_backend <- result$backend
  engine_provenance <- result$provenance
  result$facade_backend <- backend
  result$provenance <- list(
    schema_version = 1L,
    sirap_version = .sirap_package_version(),
    rrad_version = .sirap_rrad_version(),
    backend = backend,
    started_at = format(started, tz = "UTC", usetz = TRUE),
    finished_at = format(finished, tz = "UTC", usetz = TRUE),
    elapsed_seconds = unname(as.numeric(difftime(finished, started, units = "secs"))),
    core_handshake = unclass(core),
    engine_backend = engine_backend,
    engine = engine_provenance
  )
  class(result) <- unique(c("sirap_rad_result", class(result)))
  result
}

.sirap_cli_core_info <- function(rad = NULL) {
  path <- if (is.null(rad)) radPath() else radPath(rad)
  version_result <- .rad_call(
    "--version", rad = path, capture = TRUE, echo = FALSE, check = TRUE
  )
  version_line <- trimws(strsplit(version_result$stdout, "\n", fixed = TRUE)[[1L]][1L])
  version <- sub("^rad[[:space:]]+", "", version_line)
  if (!nzchar(version) || identical(version, version_line)) {
    version <- "unknown"
  }

  list(
    core_version = version,
    core_api_version = NA_integer_,
    source_commit = NA_character_,
    embedded_source_digest = NA_character_,
    resource_schema = NA_integer_,
    resource_digest = NA_character_,
    features = "demux",
    build = list(
      executable = path,
      md5 = unname(as.character(tools::md5sum(path)))
    ),
    vendored = FALSE,
    handshake = "legacy-cli-version-only"
  )
}

.sirap_validate_cli_core <- function(info) {
  version <- if (is.list(info)) info$core_version else NULL
  if (!is.character(version) || length(version) != 1L || is.na(version) ||
      !(version %in% .sirap_rad_supported_cli_versions)) {
    rlang::abort(
      paste0(
        "SIRAP's legacy CLI adapter supports RAD ",
        paste(.sirap_rad_supported_cli_versions, collapse = ", "),
        "; the selected executable reported ", .sirap_rad_value(version),
        ". Use the default native backend for a fully verified engine."
      ),
      class = "sirap_rad_cli_compatibility_error"
    )
  }
  invisible(info)
}

.sirap_validate_cli_identity <- function(before, after) {
  fields <- list(
    core_version = c(before$core_version, after$core_version),
    executable = c(before$build$executable, after$build$executable),
    md5 = c(before$build$md5, after$build$md5)
  )
  changed <- names(fields)[vapply(fields, function(values) {
    length(values) != 2L || anyNA(values) || !identical(values[[1L]], values[[2L]])
  }, logical(1))]
  if (length(changed)) {
    rlang::abort(
      paste0(
        "The RAD CLI identity changed between preflight and completion (",
        paste(changed, collapse = ", "), ")."
      ),
      class = "sirap_rad_cli_compatibility_error"
    )
  }
  invisible(after)
}

.sirap_rad_incompatible <- function(message) {
  rlang::abort(
    paste0(message, " Install a SIRAP-compatible release of `rrad`."),
    class = "sirap_rad_incompatible_error"
  )
}

.sirap_rad_value <- function(value) {
  if (is.null(value) || !length(value)) return("<missing>")
  paste(as.character(value), collapse = ", ")
}

.sirap_package_version <- function() {
  tryCatch(
    as.character(utils::packageVersion("SIRAP")),
    error = function(e) NA_character_
  )
}

.sirap_rrad_version <- function() {
  tryCatch(
    as.character(utils::packageVersion("rrad")),
    error = function(e) NA_character_
  )
}
