.sirap_test_scan_core <- function(
    commit = paste(rep("a", 40L), collapse = ""),
    features = c("demux", "scan-wl", "file-fastq", "gzip")) {
  list(
    core_version = "1.1.0",
    core_api_version = 1L,
    source_commit = commit,
    embedded_source_digest = paste(rep("c", 64L), collapse = ""),
    resource_schema = 1L,
    resource_digest = paste(rep("b", 64L), collapse = ""),
    features = features,
    build = list(
      platform = "test",
      openmp = TRUE,
      resources_available = TRUE,
      resource_dir = "/mock/rrad/resources"
    ),
    vendored = TRUE
  )
}

.sirap_test_scan_result <- function(args, core = .sirap_test_scan_core()) {
  structure(
    list(
      status = "success",
      backend = "embedded",
      core = core,
      input = args$input %||% "reads.fastq",
      output_prefix = args$output_prefix %||% "cells",
      config = args,
      artifacts = list(
        whitelist = paste0(args$output_prefix %||% "cells", ".txt"),
        statistics = paste0(args$output_prefix %||% "cells", ".csv"),
        scan_log = paste0(args$output_prefix %||% "cells", "_scan_wl.log")
      ),
      stats = list(reads_processed = 2L, selected_barcodes = 2L),
      log = "mock scan",
      provenance = list(source = "mock native scan")
    ),
    class = c("rad_scan_wl_result", "list")
  )
}

.sirap_test_scan_batch_result <- function(args,
                                          core = .sirap_test_scan_core()) {
  run_one <- .sirap_test_scan_result(
    list(input = "one.fastq", output_prefix = "out/one"), core
  )
  run_two <- .sirap_test_scan_result(
    list(input = "two.fastq", output_prefix = "out/two"), core
  )
  structure(
    list(
      status = "success",
      backend = "embedded",
      core = core,
      batch_csv = args$batch_csv %||% "batch.csv",
      config = args,
      runs = list(one = run_one, two = run_two),
      provenance = list(source = "mock native batch")
    ),
    class = c("rad_scan_wl_batch_result", "list")
  )
}

.sirap_test_scan_call <- function(...) {
  rad_scan_wl(
    input = "reads.fastq",
    output_prefix = "cells",
    adapter_seq = "ACGT",
    ...
  )
}

test_that("rad_scan_wl mirrors the complete rrad interface", {
  sirap_formals <- formals(rad_scan_wl)
  rrad_formals <- formals(rrad::rad_scan_wl)

  expect_identical(names(sirap_formals), names(rrad_formals))
  expect_identical(
    vapply(sirap_formals, function(value) {
      paste(deparse(value), collapse = " ")
    }, character(1)),
    vapply(rrad_formals, function(value) {
      paste(deparse(value), collapse = " ")
    }, character(1))
  )
})

test_that("rad_scan_wl delegates every native option and adds provenance", {
  core <- .sirap_test_scan_core()
  observed <- new.env(parent = emptyenv())
  core_calls <- 0L
  input <- tempfile(fileext = ".fastq")
  bc1 <- tempfile(fileext = ".txt")
  bc2 <- tempfile(fileext = ".txt")
  writeLines(c("@read", "ACGTAACC", "+", "IIIIIIII"), input)
  writeLines("AACC", bc1)
  writeLines("TTGA", bc2)
  expected <- list(
    input = input,
    output_prefix = tempfile(pattern = "sirap-scan-cells-"),
    batch_csv = NULL,
    adapter_seq = "ACGT",
    barcode_length = 12L,
    left_margin = 1L,
    right_margin = 2L,
    max_reads = 250L,
    max_error = 0.2,
    whitelist = NULL,
    threads = 3L,
    chunk_size = 64L,
    selection = "above_floor",
    rescan = FALSE,
    bc1_whitelist = bc1,
    bc2_whitelist = bc2,
    bc1_length = 8L,
    bc2_length = 8L,
    umi_length = 10L,
    offset_min = 1L,
    offset_max = 4L,
    verbose = TRUE
  )
  testthat::local_mocked_bindings(
    .rrad_core_info = function() {
      core_calls <<- core_calls + 1L
      core
    },
    .rrad_scan_wl = function(...) {
      observed$args <- list(...)
      .sirap_test_scan_result(observed$args, core)
    },
    .package = "SIRAP"
  )

  result <- do.call(rad_scan_wl, expected)

  expect_identical(observed$args, expected)
  expect_identical(core_calls, 2L)
  expect_s3_class(result, "sirap_rad_scan_wl_result")
  expect_s3_class(result, "sirap_rad_result")
  expect_s3_class(result, "rad_scan_wl_result")
  expect_identical(result$status, "success")
  expect_identical(result$backend, "embedded")
  expect_identical(result$facade_backend, "native")
  expect_identical(result$provenance$schema_version, 1L)
  expect_identical(result$provenance$backend, "native")
  expect_identical(result$provenance$engine_backend, "embedded")
  expect_identical(
    result$provenance$engine,
    list(source = "mock native scan")
  )
  expect_identical(
    result$provenance$core_handshake$source_commit,
    core$source_commit
  )
  expect_true(is.numeric(result$provenance$elapsed_seconds))
  expect_output(
    print(result),
    "<SIRAP RAD whitelist-scan result>",
    fixed = TRUE
  )
})

test_that("rad_scan_wl preserves built-in whitelist and rescan arguments", {
  core <- .sirap_test_scan_core()
  calls <- list()
  testthat::local_mocked_bindings(
    .rrad_core_info = function() core,
    .rrad_scan_wl = function(...) {
      calls[[length(calls) + 1L]] <<- list(...)
      .sirap_test_scan_result(calls[[length(calls)]], core)
    },
    .package = "SIRAP"
  )

  whitelist_result <- rad_scan_wl(
    input = "reads.fastq",
    output_prefix = "cells",
    adapter_seq = "ACGT",
    whitelist = "10x_3v3"
  )
  rescan_result <- rad_scan_wl(
    input = "counts.csv",
    output_prefix = "recalled",
    rescan = TRUE,
    selection = "high_specificity"
  )

  expect_s3_class(whitelist_result, "sirap_rad_scan_wl_result")
  expect_identical(calls[[1L]]$whitelist, "10x_3v3")
  expect_null(calls[[1L]]$bc1_whitelist)
  expect_s3_class(rescan_result, "sirap_rad_scan_wl_result")
  expect_true(calls[[2L]]$rescan)
  expect_null(calls[[2L]]$adapter_seq)
  expect_identical(calls[[2L]]$selection, "high_specificity")
})

test_that("rad_scan_wl preserves batch results and validates every run", {
  core <- .sirap_test_scan_core()
  observed <- NULL
  testthat::local_mocked_bindings(
    .rrad_core_info = function() core,
    .rrad_scan_wl = function(...) {
      observed <<- list(...)
      .sirap_test_scan_batch_result(observed, core)
    },
    .package = "SIRAP"
  )

  result <- rad_scan_wl(batch_csv = "batch.csv", adapter_seq = "ACGT")

  expect_identical(observed$batch_csv, "batch.csv")
  expect_null(observed$input)
  expect_null(observed$output_prefix)
  expect_s3_class(result, "sirap_rad_scan_wl_batch_result")
  expect_s3_class(result, "sirap_rad_scan_wl_result")
  expect_s3_class(result, "sirap_rad_result")
  expect_s3_class(result, "rad_scan_wl_batch_result")
  expect_length(result$runs, 2L)
  expect_true(all(vapply(
    result$runs, inherits, logical(1), what = "rad_scan_wl_result"
  )))
  expect_false(any(vapply(
    result$runs, inherits, logical(1), what = "sirap_rad_result"
  )))
  expect_identical(
    result$provenance$engine,
    list(source = "mock native batch")
  )
  expect_output(
    print(result),
    "<SIRAP RAD whitelist-scan batch result>",
    fixed = TRUE
  )
})

test_that("rad_scan_wl requires the native scan-wl capability", {
  core <- .sirap_test_scan_core(features = "demux")
  scan_calls <- 0L
  testthat::local_mocked_bindings(
    .rrad_core_info = function() core,
    .rrad_scan_wl = function(...) {
      scan_calls <<- scan_calls + 1L
      .sirap_test_scan_result(list(...), core)
    },
    .package = "SIRAP"
  )

  expect_error(
    .sirap_test_scan_call(),
    "scan-wl capability",
    class = "sirap_rad_incompatible_error"
  )
  expect_identical(scan_calls, 0L)
})

test_that("rad_scan_wl detects result and postflight identity changes", {
  preflight <- .sirap_test_scan_core(
    commit = paste(rep("a", 40L), collapse = "")
  )
  changed <- .sirap_test_scan_core(
    commit = paste(rep("d", 40L), collapse = "")
  )
  testthat::local_mocked_bindings(
    .rrad_core_info = function() preflight,
    .rrad_scan_wl = function(...) {
      .sirap_test_scan_result(list(...), changed)
    },
    .package = "SIRAP"
  )
  expect_error(
    .sirap_test_scan_call(),
    "whitelist-scan result",
    class = "sirap_rad_incompatible_error"
  )

  core_calls <- 0L
  testthat::local_mocked_bindings(
    .rrad_core_info = function() {
      core_calls <<- core_calls + 1L
      if (core_calls == 1L) preflight else changed
    },
    .rrad_scan_wl = function(...) {
      .sirap_test_scan_result(list(...), preflight)
    },
    .package = "SIRAP"
  )
  expect_error(
    .sirap_test_scan_call(),
    "preflight and postflight handshakes",
    class = "sirap_rad_incompatible_error"
  )
  expect_identical(core_calls, 2L)
})

test_that("rad_scan_wl rejects malformed or unsuccessful native results", {
  core <- .sirap_test_scan_core()
  response <- new.env(parent = emptyenv())
  response$value <- .sirap_test_scan_result(
    list(input = "reads.fastq", output_prefix = "cells"), core
  )
  testthat::local_mocked_bindings(
    .rrad_core_info = function() core,
    .rrad_scan_wl = function(...) response$value,
    .package = "SIRAP"
  )

  wrong_class <- unclass(response$value)
  response$value <- wrong_class
  expect_error(
    .sirap_test_scan_call(),
    "did not return",
    class = "sirap_rad_incompatible_error"
  )

  response$value <- .sirap_test_scan_result(
    list(input = "reads.fastq", output_prefix = "cells"), core
  )
  response$value$artifacts <- NULL
  expect_error(
    .sirap_test_scan_call(),
    "missing: artifacts",
    class = "sirap_rad_incompatible_error"
  )

  response$value <- .sirap_test_scan_result(
    list(input = "reads.fastq", output_prefix = "cells"), core
  )
  response$value$status <- "failed"
  expect_error(
    .sirap_test_scan_call(),
    "returned status 'failed'",
    class = "sirap_rad_run_error"
  )

  response$value <- .sirap_test_scan_result(
    list(input = "reads.fastq", output_prefix = "cells"), core
  )
  response$value$backend <- "cli"
  expect_error(
    .sirap_test_scan_call(),
    "embedded whitelist-scan backend",
    class = "sirap_rad_incompatible_error"
  )

  response$value <- .sirap_test_scan_result(
    list(input = "reads.fastq", output_prefix = "cells"), core
  )
  response$value$input <- character()
  expect_error(
    .sirap_test_scan_call(),
    "invalid metadata",
    class = "sirap_rad_incompatible_error"
  )
})

test_that("rad_scan_wl rejects a malformed nested batch run", {
  core <- .sirap_test_scan_core()
  testthat::local_mocked_bindings(
    .rrad_core_info = function() core,
    .rrad_scan_wl = function(...) {
      result <- .sirap_test_scan_batch_result(list(...), core)
      result$runs[[2L]]$core$resource_digest <- paste(
        rep("d", 64L), collapse = ""
      )
      result
    },
    .package = "SIRAP"
  )

  expect_error(
    rad_scan_wl(batch_csv = "batch.csv", adapter_seq = "ACGT"),
    "whitelist-scan result",
    class = "sirap_rad_incompatible_error"
  )
})

test_that("SIRAP runs a real whitelist scan through the packaged RAD engine", {
  work <- withr::local_tempdir(pattern = "sirap-rad-native-scan-")
  input <- file.path(work, "reads.fastq")
  writeLines(c(
    "@read-one", "ACGTAACC", "+", "IIIIIIII",
    "@read-two", "ACGTTTGA", "+", "IIIIIIII"
  ), input)

  result <- rad_scan_wl(
    input = input,
    output_prefix = file.path(work, "cells"),
    adapter_seq = "ACGT",
    barcode_length = 4L,
    max_error = 0,
    threads = 1L,
    chunk_size = 1L,
    selection = "above_floor"
  )

  expect_s3_class(result, "sirap_rad_scan_wl_result")
  expect_s3_class(result, "sirap_rad_result")
  expect_s3_class(result, "rad_scan_wl_result")
  expect_identical(result$status, "success")
  expect_identical(result$backend, "embedded")
  expect_identical(result$facade_backend, "native")
  expect_equal(result$stats$reads_processed, 2)
  expect_equal(result$stats$selected_barcodes, 2)
  expect_true(all(file.exists(unlist(result$artifacts[c(
    "whitelist", "statistics", "scan_log"
  )]))))
  expect_setequal(
    readLines(result$artifacts$whitelist, warn = FALSE),
    c("AACC", "TTGA")
  )
  expect_identical(
    result$provenance$core_handshake$resource_digest,
    result$core$resource_digest
  )
  expect_output(print(result), "Input:", fixed = TRUE)
})
