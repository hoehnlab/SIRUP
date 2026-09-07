sirap_reformat_mock_core <- function(
    commit = paste(rep("a", 40L), collapse = ""),
    source_digest = paste(rep("c", 64L), collapse = ""),
    resource_digest = paste(rep("b", 64L), collapse = ""),
    features = c(
      "demux", "reformat", "reformat-fasta", "reformat-presto-header",
      "file-fastq", "gzip"
    )) {
  list(
    core_version = "1.1.0",
    core_api_version = 1L,
    source_commit = commit,
    embedded_source_digest = source_digest,
    resource_schema = 1L,
    resource_digest = resource_digest,
    features = features,
    build = list(
      platform = "test",
      resources_available = TRUE,
      resource_dir = "/mock/rrad/resources"
    ),
    vendored = TRUE
  )
}

sirap_reformat_mock_result <- function(core, input = "/input.fastq.gz",
                                       output = "/output/by-cell",
                                       status = "success",
                                       backend = "embedded") {
  structure(
    list(
      status = status,
      backend = backend,
      core = core,
      input = input,
      output = output,
      config = list(split_bc = TRUE),
      artifacts = list(split_fastq = character()),
      stats = list(records_read = 2, records_written = 2),
      log = "",
      provenance = list(owner = "rrad")
    ),
    class = c("rad_reformat_result", "list")
  )
}

test_that("SIRAP exposes the complete native rrad reformat interface", {
  expect_identical(formals(rad_reformat), formals(rrad::rad_reformat))
})

test_that("rad_reformat delegates every option and adds SIRAP provenance", {
  core <- sirap_reformat_mock_core()
  core_calls <- 0L
  observed <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    .rrad_core_info = function() {
      core_calls <<- core_calls + 1L
      core
    },
    .rrad_reformat = function(...) {
      observed$args <- list(...)
      sirap_reformat_mock_result(
        core,
        input = observed$args$fastq,
        output = observed$args$out_dir
      )
    },
    .package = "SIRAP"
  )

  result <- rad_reformat(
    fastq = "/input.fastq.gz",
    out_dir = "/output/by-cell",
    split_bc = TRUE,
    reformat_header = TRUE,
    collapsed_input = TRUE,
    delimiter = "-",
    coordinate = NULL,
    bin_size = 8L,
    output_fastq = NULL,
    threads = 3L,
    chunk_size = 128L,
    verbose = TRUE
  )

  expect_identical(core_calls, 2L)
  expect_named(observed$args, names(formals(rrad::rad_reformat)))
  expect_identical(observed$args$fastq, "/input.fastq.gz")
  expect_identical(observed$args$out_dir, "/output/by-cell")
  expect_true(observed$args$split_bc)
  expect_true(observed$args$reformat_header)
  expect_true(observed$args$collapsed_input)
  expect_identical(observed$args$delimiter, "-")
  expect_null(observed$args$coordinate)
  expect_identical(observed$args$bin_size, 8L)
  expect_null(observed$args$output_fastq)
  expect_identical(observed$args$threads, 3L)
  expect_identical(observed$args$chunk_size, 128L)
  expect_true(observed$args$verbose)
  expect_false(observed$args$to_fasta)
  expect_null(observed$args$header_format)

  expect_s3_class(result, "sirap_rad_reformat_result")
  expect_s3_class(result, "sirap_rad_result")
  expect_s3_class(result, "rad_reformat_result")
  expect_identical(result$status, "success")
  expect_identical(result$backend, "embedded")
  expect_identical(result$facade_backend, "native")
  expect_identical(result$provenance$schema_version, 1L)
  expect_identical(result$provenance$backend, "native")
  expect_identical(result$provenance$engine_backend, "embedded")
  expect_identical(result$provenance$engine$owner, "rrad")
  expect_identical(
    result$provenance$core_handshake$resource_digest,
    core$resource_digest
  )
  expect_true(is.numeric(result$provenance$elapsed_seconds))
  expect_output(
    print.sirap_rad_reformat_result(result),
    "SIRAP RAD reformat result"
  )
})

test_that("rad_reformat requires the native reformat capability", {
  called <- FALSE
  core <- sirap_reformat_mock_core(features = "demux")
  testthat::local_mocked_bindings(
    .rrad_core_info = function() core,
    .rrad_reformat = function(...) {
      called <<- TRUE
      stop("must not run")
    },
    .package = "SIRAP"
  )

  expect_error(
    rad_reformat("input.fastq", reformat_header = TRUE),
    "reformat feature",
    class = "sirap_rad_incompatible_error"
  )
  expect_false(called)
})

test_that("rad_reformat negotiates FASTA and pRESTO features", {
  called <- FALSE
  core <- sirap_reformat_mock_core(features = c("demux", "reformat"))
  testthat::local_mocked_bindings(
    .rrad_core_info = function() core,
    .rrad_reformat = function(...) {
      called <<- TRUE
      stop("must not run")
    },
    .package = "SIRAP"
  )

  expect_error(
    rad_reformat(
      "input.fastq", output_fastq = "input.fa", to_fasta = TRUE
    ),
    "reformat-fasta",
    class = "sirap_rad_incompatible_error"
  )
  expect_false(called)

  core <- sirap_reformat_mock_core(
    features = c("demux", "reformat", "reformat-fasta")
  )
  expect_error(
    rad_reformat(
      "input.fastq", output_fastq = "input.fa", to_fasta = TRUE,
      header_format = "presto"
    ),
    "reformat-presto-header",
    class = "sirap_rad_incompatible_error"
  )
  expect_false(called)
})

test_that("rad_reformat rejects a result from a different native core", {
  preflight <- sirap_reformat_mock_core()
  result_core <- sirap_reformat_mock_core(
    source_digest = paste(rep("d", 64L), collapse = "")
  )
  testthat::local_mocked_bindings(
    .rrad_core_info = function() preflight,
    .rrad_reformat = function(...) {
      sirap_reformat_mock_result(result_core)
    },
    .package = "SIRAP"
  )

  expect_error(
    rad_reformat("input.fastq", reformat_header = TRUE),
    "preflight handshake and the reformat result",
    class = "sirap_rad_incompatible_error"
  )
})

test_that("rad_reformat rejects a postflight engine identity change", {
  preflight <- sirap_reformat_mock_core()
  postflight <- sirap_reformat_mock_core(
    resource_digest = paste(rep("e", 64L), collapse = "")
  )
  core_calls <- 0L
  testthat::local_mocked_bindings(
    .rrad_core_info = function() {
      core_calls <<- core_calls + 1L
      if (core_calls == 1L) preflight else postflight
    },
    .rrad_reformat = function(...) {
      sirap_reformat_mock_result(preflight)
    },
    .package = "SIRAP"
  )

  expect_error(
    rad_reformat("input.fastq", reformat_header = TRUE),
    "preflight and postflight handshakes",
    class = "sirap_rad_incompatible_error"
  )
  expect_identical(core_calls, 2L)
})

test_that("rad_reformat rejects malformed and unsuccessful native results", {
  core <- sirap_reformat_mock_core()
  testthat::local_mocked_bindings(
    .rrad_core_info = function() core,
    .rrad_reformat = function(...) list(status = "success"),
    .package = "SIRAP"
  )
  expect_error(
    rad_reformat("input.fastq", reformat_header = TRUE),
    "did not return a `rad_reformat_result`",
    class = "sirap_rad_incompatible_error"
  )

  testthat::local_mocked_bindings(
    .rrad_core_info = function() core,
    .rrad_reformat = function(...) {
      sirap_reformat_mock_result(core, status = "failed")
    },
    .package = "SIRAP"
  )
  expect_error(
    rad_reformat("input.fastq", reformat_header = TRUE),
    "returned status 'failed'",
    class = "sirap_rad_run_error"
  )
})

test_that("SIRAP reformats a tagged FASTQ through the real native core", {
  work <- withr::local_tempdir(pattern = "sirap-rad-reformat-")
  input <- file.path(work, "tagged.fastq")
  output <- file.path(work, "collapsed.fastq.gz")
  writeLines(
    c(
      "@read-one CB:Z:AACC UB:Z:TGCA extra:Z:retained-only-without-collapse",
      "ACGT",
      "+read-one",
      "!!!!"
    ),
    input
  )

  result <- rad_reformat(
    fastq = input,
    reformat_header = TRUE,
    output_fastq = output,
    threads = 1L,
    chunk_size = 1L
  )

  connection <- gzfile(output, open = "rt")
  on.exit(close(connection), add = TRUE)
  lines <- readLines(connection, warn = FALSE)

  expect_s3_class(result, "sirap_rad_reformat_result")
  expect_s3_class(result, "sirap_rad_result")
  expect_s3_class(result, "rad_reformat_result")
  expect_identical(result$status, "success")
  expect_identical(result$facade_backend, "native")
  expect_identical(lines[[1L]], "@read-one_AACC_TGCA")
  expect_identical(lines[[3L]], "+read-one")
  expect_identical(result$stats$records_read, 1)
  expect_identical(result$stats$records_written, 1)
  expect_identical(
    result$provenance$core_handshake$resource_digest,
    result$core$resource_digest
  )
})

test_that("SIRAP streams RAD FASTQ to canonical IgBLAST FASTA", {
  work <- withr::local_tempdir(pattern = "sirap-rad-igblast-fasta-")
  input <- file.path(work, "tagged.fastq.gz")
  output <- file.path(work, "igblast.fasta")
  connection <- gzfile(input, open = "wt")
  writeLines(
    c(
      "@read-F TS:A:+ UB:Z:TGCA CB:Z:AACC BC:Z:s_002um_001_002-1 CR:Z:raw",
      "ACGT",
      "+read-F",
      "IIII"
    ),
    connection
  )
  close(connection)

  result <- rad_reformat(
    input,
    output_fastq = output,
    threads = 1L,
    to_fasta = TRUE,
    header_format = "presto"
  )

  expect_identical(
    readLines(output, warn = FALSE),
    c(">read-F|BARCODE=s_002um_001_002-1|UMI=TGCA", "ACGT")
  )
  expect_identical(result$artifacts$fasta, normalizePath(output))
  expect_identical(result$config$header_format, "presto")
  expect_identical(result$stats$records_converted_to_fasta, 1)
})
