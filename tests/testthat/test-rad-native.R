mock_rad_core <- function(commit = paste(rep("a", 40L), collapse = ""),
                          api = 1L, schema = 1L,
                          features = c(
                            "demux", "file-fastq", "gzip", "scan-wl",
                            "auto-whitelist", "debug-output"
                          )) {
  list(
    core_version = "1.1.0",
    core_api_version = api,
    source_commit = commit,
    embedded_source_digest = paste(rep("c", 64L), collapse = ""),
    resource_schema = schema,
    resource_digest = paste(rep("b", 64L), collapse = ""),
    features = features,
    build = list(
      platform = "test",
      resources_available = TRUE,
      resource_dir = "/mock/rrad/resources"
    ),
    vendored = TRUE
  )
}

mock_rad_scan_result <- function(core = mock_rad_core(),
                                 fastq = "reads.fq.gz") {
  structure(
    list(
      status = "success",
      backend = "embedded",
      core = core,
      input = fastq,
      output_prefix = "run/demo_scanwl",
      config = list(selection = "high_specificity"),
      artifacts = list(whitelist = "run/demo_scanwl.txt"),
      stats = list(reads_processed = 12L),
      log = ""
    ),
    class = c("rad_scan_wl_result", "list")
  )
}

mock_rad_result <- function(core = mock_rad_core(), fastq = "reads.fq.gz",
                            auto_whitelist = FALSE) {
  structure(
    list(
      status = "success",
      backend = "embedded",
      core = core,
      input = list(fastq = fastq, layout = "five_prime"),
      config = list(
        layout = "five_prime",
        auto_whitelist = auto_whitelist
      ),
      artifacts = list(fastq = "run/demo.fq.gz"),
      stats = list(reads_seen = 12L),
      scan = if (auto_whitelist) mock_rad_scan_result(core, fastq) else NULL,
      log = ""
    ),
    class = c("rad_demux_result", "list")
  )
}

write_test_fastq_gz <- function(path, records = c("@read", "ACGT", "+", "!!!!")) {
  connection <- gzfile(path, open = "wt")
  on.exit(close(connection), add = TRUE)
  writeLines(records, connection)
  invisible(path)
}

test_that("SIRAP runs one read through the packaged RAD engine", {
  fastq <- test_path("fixtures", "rad-native", "nanopore-smoke.fastq")
  out_dir <- withr::local_tempdir(pattern = "sirap-rad-native-")

  result <- rad_demux(
    layout = "nanopore_rapid_bc",
    fastq = fastq,
    out_dir = out_dir,
    output = "smoke",
    threads = 1L,
    chunk_size = 1L,
    max_reads = 1L,
    min_read_length = 0L
  )

  expect_s3_class(result, "sirap_rad_result")
  expect_s3_class(result, "rad_demux_result")
  expect_identical(result$status, "success")
  expect_identical(result$backend, "embedded")
  expect_identical(result$facade_backend, "native")
  expect_equal(result$stats$total_reads, 1)
  expect_equal(result$stats$reads_demultiplexed, 1)
  expect_true(file.exists(result$artifacts$fastq))
  expect_identical(result$provenance$engine_backend, "embedded")
  expect_identical(
    result$provenance$core_handshake$resource_digest,
    result$core$resource_digest
  )
})

test_that("rad_core_info enforces the native compatibility handshake", {
  core <- mock_rad_core()
  testthat::local_mocked_bindings(
    .rrad_core_info = function() core,
    .package = "SIRAP"
  )

  info <- rad_core_info()

  expect_s3_class(info, "sirap_rad_core_info")
  expect_identical(info$core_api_version, 1L)
  expect_identical(info$resource_schema, 1L)
  expect_match(info$resource_digest, "^[[:xdigit:]]{64}$")
  expect_match(info$embedded_source_digest, "^[[:xdigit:]]{64}$")
  expect_true("demux" %in% info$features)
})

test_that("rad_core_info rejects incompatible native engines", {
  testthat::local_mocked_bindings(
    .rrad_core_info = function() mock_rad_core(api = 2L),
    .package = "SIRAP"
  )
  expect_error(
    rad_core_info(),
    class = "sirap_rad_incompatible_error"
  )

  fractional_api <- mock_rad_core()
  fractional_api$core_api_version <- 1.5
  testthat::local_mocked_bindings(
    .rrad_core_info = function() fractional_api,
    .package = "SIRAP"
  )
  expect_error(
    rad_core_info(),
    class = "sirap_rad_incompatible_error"
  )

  unavailable_resources <- mock_rad_core()
  unavailable_resources$build$resources_available <- FALSE
  testthat::local_mocked_bindings(
    .rrad_core_info = function() unavailable_resources,
    .package = "SIRAP"
  )
  expect_error(
    rad_core_info(),
    "packaged resources",
    class = "sirap_rad_incompatible_error"
  )

  testthat::local_mocked_bindings(
    .rrad_core_info = function() mock_rad_core(features = "gzip"),
    .package = "SIRAP"
  )
  expect_error(
    rad_core_info(),
    class = "sirap_rad_incompatible_error"
  )

  testthat::local_mocked_bindings(
    .rrad_core_info = function() mock_rad_core(schema = 2L),
    .package = "SIRAP"
  )
  expect_error(
    rad_core_info(),
    class = "sirap_rad_incompatible_error"
  )

  bad_digest <- mock_rad_core()
  bad_digest$resource_digest <- "not-a-sha256"
  testthat::local_mocked_bindings(
    .rrad_core_info = function() bad_digest,
    .package = "SIRAP"
  )
  expect_error(
    rad_core_info(),
    "resource SHA-256",
    class = "sirap_rad_incompatible_error"
  )

  bad_source_digest <- mock_rad_core()
  bad_source_digest$embedded_source_digest <- "not-a-sha256"
  testthat::local_mocked_bindings(
    .rrad_core_info = function() bad_source_digest,
    .package = "SIRAP"
  )
  expect_error(
    rad_core_info(),
    "embedded-source SHA-256",
    class = "sirap_rad_incompatible_error"
  )

  bad_commit <- mock_rad_core(commit = "abc1234")
  testthat::local_mocked_bindings(
    .rrad_core_info = function() bad_commit,
    .package = "SIRAP"
  )
  expect_error(
    rad_core_info(),
    "source commit",
    class = "sirap_rad_incompatible_error"
  )
})

test_that("rad_demux delegates every native option and adds provenance", {
  core <- mock_rad_core()
  observed <- new.env(parent = emptyenv())
  fastq <- test_path("fixtures", "rad-native", "nanopore-smoke.fastq")
  out_dir <- withr::local_tempdir(pattern = "sirap-rad-delegate-")
  testthat::local_mocked_bindings(
    .rrad_core_info = function() core,
    .rrad_demux = function(...) {
      observed$args <- list(...)
      mock_rad_result(
        core = core,
        fastq = observed$args$fastq,
        auto_whitelist = observed$args$auto_whitelist
      )
    },
    .package = "SIRAP"
  )

  result <- rad_demux(
    layout = "five_prime",
    fastq = fastq,
    out_dir = out_dir,
    output = "sample",
    kit = "10x_5v1",
    threads = 4L,
    chunk_size = 1000L,
    max_reads = 2500L,
    auto_whitelist = TRUE,
    scan_adapter = "ACGTACGT",
    scan_barcode_length = 17L,
    scan_max_error = 0.2,
    scan_max_reads = 2000L,
    scan_chunk_size = 800L,
    scan_threads = 3L,
    scan_selection = "above_floor",
    whitelist_mutation = 1L,
    generated_mutation = 2L,
    bc_correction_mode = "defensive",
    joint_bc_mode = "strict",
    rc_umi = FALSE,
    min_read_length = 100L,
    write_debug = TRUE,
    verbose = TRUE,
    reuse_cache = TRUE
  )

  expect_s3_class(result, "sirap_rad_result")
  expect_s3_class(result, "rad_demux_result")
  expect_identical(result$status, "success")
  expect_identical(result$backend, "embedded")
  expect_identical(result$facade_backend, "native")
  expect_identical(observed$args$layout, "five_prime")
  expect_named(
    observed$args,
    c(
      "layout", "fastq", "out_dir", "output", "kit", "global_whitelist",
      "custom_whitelist", "threads", "chunk_size", "max_reads",
      "auto_whitelist", "scan_adapter", "scan_barcode_length",
      "scan_max_error", "scan_max_reads", "scan_chunk_size",
      "scan_threads", "scan_selection", "whitelist_mutation",
      "generated_mutation",
      "bc_correction_mode", "joint_bc_mode", "rc_umi", "min_read_length",
      "write_debug", "verbose", "reuse_cache"
    )
  )
  expect_identical(observed$args$fastq, normalizePath(fastq))
  expect_identical(observed$args$out_dir, normalizePath(out_dir))
  expect_identical(observed$args$output, "sample")
  expect_identical(observed$args$threads, 4L)
  expect_true(observed$args$auto_whitelist)
  expect_identical(observed$args$scan_adapter, "ACGTACGT")
  expect_identical(observed$args$scan_barcode_length, 17L)
  expect_identical(observed$args$scan_max_error, 0.2)
  expect_identical(observed$args$scan_max_reads, 2000L)
  expect_identical(observed$args$scan_chunk_size, 800L)
  expect_identical(observed$args$scan_threads, 3L)
  expect_identical(observed$args$scan_selection, "above_floor")
  expect_false(observed$args$rc_umi)
  expect_true(observed$args$write_debug)
  expect_true(observed$args$reuse_cache)
  expect_identical(result$provenance$schema_version, 1L)
  expect_identical(result$provenance$backend, "native")
  expect_true(is.character(result$provenance$rrad_version))
  expect_identical(result$provenance$engine_backend, "embedded")
  expect_identical(
    result$provenance$core_handshake$source_commit,
    paste(rep("a", 40L), collapse = "")
  )
  expect_identical(
    result$provenance$core_handshake$embedded_source_digest,
    paste(rep("c", 64L), collapse = "")
  )
  expect_identical(
    result$provenance$core_handshake$resource_digest,
    paste(rep("b", 64L), collapse = "")
  )
  expect_true(is.numeric(result$provenance$elapsed_seconds))
})

test_that("rad_demux preserves the complete native rrad interface", {
  rrad_formals <- formals(rrad::rad_demux)
  sirap_formals <- formals(rad_demux)

  expect_identical(
    as.list(sirap_formals)[names(rrad_formals)],
    as.list(rrad_formals)
  )
  expect_identical(
    setdiff(names(sirap_formals), names(rrad_formals)),
    c("backend", "rad", "capture")
  )
})

test_that("rad_demux checks configuration-sensitive native capabilities", {
  fastq <- test_path("fixtures", "rad-native", "nanopore-smoke.fastq")
  cases <- list(
    list(
      args = list(auto_whitelist = TRUE),
      missing = c("auto-whitelist", "scan-wl")
    ),
    list(
      args = list(write_debug = TRUE),
      missing = "debug-output"
    )
  )

  for (case in cases) {
    core <- mock_rad_core(
      features = setdiff(mock_rad_core()$features, case$missing)
    )
    called <- FALSE
    testthat::local_mocked_bindings(
      .rrad_core_info = function() core,
      .rrad_demux = function(...) {
        called <<- TRUE
        mock_rad_result(core)
      },
      .package = "SIRAP"
    )
    args <- c(
      list(
        layout = "five_prime",
        fastq = fastq,
        out_dir = withr::local_tempdir(pattern = "sirap-rad-feature-")
      ),
      case$args
    )

    expect_error(
      do.call(rad_demux, args),
      paste(case$missing, collapse = ".*"),
      class = "sirap_rad_incompatible_error"
    )
    expect_false(called)
  }
})

test_that("rad_demux validates the complete native result contract", {
  core <- mock_rad_core()
  fastq <- test_path("fixtures", "rad-native", "nanopore-smoke.fastq")
  out_dir <- withr::local_tempdir(pattern = "sirap-rad-result-contract-")

  malformed <- list(
    missing_scan_field = local({
      value <- mock_rad_result(core)
      value$scan <- NULL
      value$config$auto_whitelist <- TRUE
      value
    }),
    invalid_input = local({
      value <- mock_rad_result(core)
      value$input <- "reads.fq.gz"
      value
    }),
    unexpected_scan = local({
      value <- mock_rad_result(core)
      value$scan <- mock_rad_scan_result(core)
      value
    }),
    mismatched_scan_core = local({
      value <- mock_rad_result(core, auto_whitelist = TRUE)
      value$scan$core$source_commit <- paste(rep("d", 40L), collapse = "")
      value
    })
  )

  for (name in names(malformed)) {
    result <- malformed[[name]]
    testthat::local_mocked_bindings(
      .rrad_core_info = function() core,
      .rrad_demux = function(...) result,
      .package = "SIRAP"
    )
    expect_error(
      rad_demux("five_prime", fastq, out_dir = out_dir),
      class = "sirap_rad_incompatible_error",
      info = name
    )
  }
})

test_that("rad_demux rejects a result produced by a different core", {
  preflight <- mock_rad_core(commit = paste(rep("b", 40L), collapse = ""))
  completed <- mock_rad_core(commit = paste(rep("c", 40L), collapse = ""))
  testthat::local_mocked_bindings(
    .rrad_core_info = function() preflight,
    .rrad_demux = function(...) mock_rad_result(completed),
    .package = "SIRAP"
  )

  fastq <- test_path("fixtures", "rad-native", "nanopore-smoke.fastq")

  expect_error(
    rad_demux(
      "five_prime", fastq,
      out_dir = withr::local_tempdir(pattern = "sirap-rad-mismatch-")
    ),
    "changed between the preflight handshake",
    class = "sirap_rad_incompatible_error"
  )
})

test_that("rad_demux rejects a postflight native identity change", {
  preflight <- mock_rad_core(commit = paste(rep("b", 40L), collapse = ""))
  postflight <- mock_rad_core(commit = paste(rep("c", 40L), collapse = ""))
  calls <- 0L
  testthat::local_mocked_bindings(
    .rrad_core_info = function() {
      calls <<- calls + 1L
      if (calls == 1L) preflight else postflight
    },
    .rrad_demux = function(...) mock_rad_result(preflight),
    .package = "SIRAP"
  )

  fastq <- test_path("fixtures", "rad-native", "nanopore-smoke.fastq")
  expect_error(
    rad_demux(
      "five_prime", fastq,
      out_dir = withr::local_tempdir(pattern = "sirap-rad-postflight-")
    ),
    "preflight and postflight handshakes",
    class = "sirap_rad_incompatible_error"
  )
})

test_that("the CLI backend is explicit and normalized", {
  core <- list(
    core_version = "1.0.2",
    core_api_version = NA_integer_,
    source_commit = NA_character_,
    embedded_source_digest = NA_character_,
    resource_schema = NA_integer_,
    resource_digest = NA_character_,
    features = "demux",
    build = list(executable = "/mock/rad", md5 = "123"),
    vendored = FALSE,
    handshake = "legacy-cli-version-only"
  )
  cli_call <- list(
    status = 0L,
    stdout = "done",
    stderr = NA_character_,
    command = c("/mock/rad", "demux")
  )
  observed <- new.env(parent = emptyenv())
  out_dir <- withr::local_tempdir(pattern = "sirap-rad-cli-")
  fastq <- test_path("fixtures", "rad-native", "nanopore-smoke.fastq")
  testthat::local_mocked_bindings(
    .sirap_cli_core_info = function(rad = NULL) core,
    .sirap_cli_demux = function(config, rad, capture) {
      observed$config <- config
      observed$rad <- rad
      observed$capture <- capture
      write_test_fastq_gz(
        file.path(config$out_dir, paste0(config$output, ".fq.gz"))
      )
      cli_call
    },
    .package = "SIRAP"
  )

  result <- rad_demux(
    "five_prime", fastq,
    out_dir = out_dir,
    backend = "cli", rad = "/mock/rad", capture = FALSE,
    auto_whitelist = TRUE, rc_umi = FALSE, min_read_length = 42L
  )

  expect_s3_class(result, "sirap_rad_result")
  expect_s3_class(result, "sirap_rad_cli_result")
  expect_identical(result$status, "success")
  expect_identical(result$backend, "cli")
  expect_identical(result$facade_backend, "cli")
  expect_identical(result$log$status, 0L)
  expect_identical(result$core$handshake, "legacy-cli-version-only")
  expect_identical(observed$rad, "/mock/rad")
  expect_false(observed$capture)
  expect_true(observed$config$auto_whitelist)
  expect_false(observed$config$rc_umi)
  expect_identical(observed$config$min_read_length, 42L)
  expect_true(is.na(result$config$reuse_cache))
  expect_false(result$config$reuse_cache_requested)
  expect_identical(result$config$cache_policy, "legacy-cli-uncontrolled")
  expect_identical(
    result$artifacts$fastq,
    normalizePath(file.path(out_dir, "demux.fq.gz"))
  )
})

test_that("the CLI adapter reuses radDemux", {
  config <- list(
    layout = "five_prime", fastq = "reads.fq.gz", out_dir = ".",
    output = "demo", kit = NULL, global_whitelist = NULL,
    custom_whitelist = NULL, threads = 2L, chunk_size = 5000L,
    max_reads = NULL, auto_whitelist = TRUE, whitelist_mutation = NULL,
    scan_adapter = NULL, scan_barcode_length = NULL, scan_max_error = 0.3,
    scan_max_reads = NULL, scan_chunk_size = NULL, scan_threads = NULL,
    scan_selection = "high_specificity",
    generated_mutation = NULL, bc_correction_mode = "offensive",
    joint_bc_mode = "default", rc_umi = FALSE, min_read_length = 10L,
    write_debug = TRUE, verbose = FALSE, reuse_cache = FALSE
  )
  observed <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    radPath = function(path = NULL) "/mock/rad",
    radDemux = function(...) {
      observed$args <- list(...)
      list(status = 0L, stdout = "done", stderr = NULL, command = "rad")
    },
    .package = "SIRAP"
  )

  result <- SIRAP:::.sirap_cli_demux(config, rad = "/mock/rad", capture = TRUE)

  expect_identical(result$status, 0L)
  expect_true(observed$args$write_dbg)
  expect_identical(observed$args$threads, 2L)
  raw_flags <- unlist(observed$args[names(observed$args) == ""], use.names = FALSE)
  expect_true("--auto-wl" %in% raw_flags)
  expect_true("--no-umi-rc" %in% raw_flags)
  expect_true("--min-read-length" %in% raw_flags)
  expect_true("10" %in% raw_flags)
})

test_that("the CLI execution boundary preserves literal arguments", {
  script_dir <- withr::local_tempdir(pattern = "sirap fake rad ")
  script <- file.path(script_dir, "fake rad")
  writeLines(
    c("#!/bin/sh", "for arg in \"$@\"; do printf '<%s>\\n' \"$arg\"; done"),
    script
  )
  Sys.chmod(script, mode = "0755")

  result <- SIRAP:::.rad_call(
    c("alpha beta", "literal; /usr/bin/printf '[INJECTED]'"),
    rad = script,
    capture = TRUE,
    echo = FALSE
  )

  expect_identical(result$status, 0L)
  expect_identical(
    result$stdout,
    paste0(
      "<alpha beta>\n",
      "<literal; /usr/bin/printf '[INJECTED]'>"
    )
  )
  expect_identical(
    result$command,
    c(script, "alpha beta", "literal; /usr/bin/printf '[INJECTED]'")
  )
})

test_that("CLI whitelist combinations cannot silently change native semantics", {
  combinations <- expand.grid(
    kit = c(FALSE, TRUE),
    global = c(FALSE, TRUE),
    custom = c(FALSE, TRUE),
    KEEP.OUT.ATTRS = FALSE
  )

  for (i in seq_len(nrow(combinations))) {
    combo <- combinations[i, ]
    config <- list(
      kit = if (combo$kit) "kit" else NULL,
      global_whitelist = if (combo$global) "global.csv" else NULL,
      custom_whitelist = if (combo$custom) "custom.csv" else NULL,
      auto_whitelist = FALSE
    )
    rejected <- (combo$kit && combo$global) ||
      (!combo$kit && combo$global && !combo$custom) ||
      (!combo$kit && !combo$global && combo$custom)

    if (rejected) {
      expect_error(
        SIRAP:::.sirap_validate_cli_whitelist_config(config),
        class = if (combo$kit && combo$global) {
          "sirap_rad_configuration_error"
        } else {
          "sirap_rad_cli_compatibility_error"
        },
        info = paste(unlist(combo), collapse = "/")
      )
    } else {
      expect_invisible(SIRAP:::.sirap_validate_cli_whitelist_config(config))
    }
  }

  expect_invisible(SIRAP:::.sirap_validate_cli_whitelist_config(list(
    kit = NULL,
    global_whitelist = "global.csv",
    custom_whitelist = NULL,
    auto_whitelist = TRUE
  )))
})

test_that("shared demux arguments are validated before backend selection", {
  fastq <- test_path("fixtures", "rad-native", "nanopore-smoke.fastq")
  invalid <- list(
    list(auto_whitelist = 1),
    list(rc_umi = NA),
    list(write_debug = 1),
    list(verbose = logical()),
    list(reuse_cache = "false"),
    list(threads = 0L),
    list(chunk_size = 1.5),
    list(max_reads = Inf),
    list(scan_adapter = ""),
    list(scan_barcode_length = 0L),
    list(scan_barcode_length = 33L),
    list(scan_max_error = -0.1),
    list(scan_max_error = 1.01),
    list(scan_max_reads = 0L),
    list(scan_chunk_size = 0L),
    list(scan_threads = 0L),
    list(scan_selection = "unknown"),
    list(min_read_length = -1L),
    list(bc_correction_mode = "unknown"),
    list(joint_bc_mode = "unknown")
  )

  for (backend in c("native", "cli")) {
    for (case in invalid) {
      args <- c(
        list(
          layout = "five_prime",
          fastq = fastq,
          out_dir = withr::local_tempdir(pattern = "sirap-rad-invalid-"),
          backend = backend
        ),
        case
      )
      expect_error(
        do.call(rad_demux, args),
        class = "sirap_rad_configuration_error",
        info = paste(backend, names(case))
      )
    }

    expect_error(
      rad_demux(
        "five_prime", fastq,
        out_dir = withr::local_tempdir(pattern = "sirap-rad-capture-"),
        backend = backend, capture = 1
      ),
      class = "sirap_rad_configuration_error",
      info = paste(backend, "capture")
    )
  }
})

test_that("legacy CLI rejects native-only automatic scan controls", {
  fastq <- test_path("fixtures", "rad-native", "nanopore-smoke.fastq")
  cases <- list(
    scan_adapter = "ACGT",
    scan_barcode_length = 16L,
    scan_max_error = 0.1,
    scan_max_reads = 10L,
    scan_chunk_size = 10L,
    scan_threads = 2L,
    scan_selection = "above_floor"
  )

  for (name in names(cases)) {
    args <- list(
      layout = "five_prime",
      fastq = fastq,
      out_dir = withr::local_tempdir(pattern = "sirap-rad-cli-scan-"),
      backend = "cli",
      auto_whitelist = TRUE
    )
    args[[name]] <- cases[[name]]
    expect_error(
      do.call(rad_demux, args),
      name,
      class = "sirap_rad_cli_compatibility_error",
      info = name
    )
  }
})

test_that("legacy CLI debug mode is rejected before it can drop primary data", {
  expect_error(
    rad_demux(
      "five_prime",
      test_path("fixtures", "rad-native", "nanopore-smoke.fastq"),
      out_dir = withr::local_tempdir(pattern = "sirap-rad-cli-debug-"),
      write_debug = TRUE,
      backend = "cli"
    ),
    "suppress the primary FASTQ",
    class = "sirap_rad_cli_compatibility_error"
  )

  expect_error(
    radDemux(
      "five_prime", "reads.fq.gz", write_dbg = TRUE,
      rad = "/not/evaluated"
    ),
    "can suppress the primary FASTQ"
  )
})

test_that("legacy CLI input is structurally and gzip validated", {
  out_dir <- withr::local_tempdir(pattern = "sirap-rad-cli-input-")
  testthat::local_mocked_bindings(
    .sirap_cli_core_info = function(rad = NULL) list(core_version = "1.0.2"),
    .package = "SIRAP"
  )
  malformed <- file.path(out_dir, "malformed.fastq")
  writeLines(c("@read", "ACGT", "+"), malformed)
  expect_error(
    rad_demux("five_prime", malformed, out_dir = out_dir, backend = "cli"),
    "incomplete four-line record",
    class = "sirap_rad_configuration_error"
  )

  compressed <- file.path(out_dir, "truncated.fastq.gz")
  write_test_fastq_gz(compressed)
  bytes <- readBin(compressed, "raw", n = file.info(compressed)$size)
  writeBin(bytes[seq_len(length(bytes) - 5L)], compressed)
  expect_error(
    rad_demux("five_prime", compressed, out_dir = out_dir, backend = "cli"),
    class = "sirap_rad_configuration_error"
  )

  unsupported <- file.path(out_dir, "reads.txt")
  writeLines(c("@read", "ACGT", "+", "!!!!"), unsupported)
  expect_error(
    rad_demux("five_prime", unsupported, out_dir = out_dir, backend = "cli"),
    "must end in",
    class = "sirap_rad_configuration_error"
  )
})

test_that("CLI version handshake rejects unknown adapter contracts", {
  supported <- list(core_version = "1.0.2")
  testthat::local_mocked_bindings(
    .sirap_cli_core_info = function(rad = NULL) supported,
    .package = "SIRAP"
  )
  expect_s3_class(rad_core_info("cli"), "sirap_rad_core_info")

  for (version in c("0.9.0", "unknown", "2.0.0")) {
    unsupported <- list(core_version = version)
    testthat::local_mocked_bindings(
      .sirap_cli_core_info = function(rad = NULL) unsupported,
      .package = "SIRAP"
    )
    expect_error(
      rad_core_info("cli"),
      "supports RAD 1.0.2",
      class = "sirap_rad_cli_compatibility_error",
      info = version
    )
  }
})

test_that("CLI identity is stable across an execution", {
  core_calls <- 0L
  out_dir <- withr::local_tempdir(pattern = "sirap-rad-cli-identity-")
  fastq <- test_path("fixtures", "rad-native", "nanopore-smoke.fastq")
  testthat::local_mocked_bindings(
    .sirap_cli_core_info = function(rad = NULL) {
      core_calls <<- core_calls + 1L
      list(
        core_version = "1.0.2",
        build = list(
          executable = "/mock/rad",
          md5 = if (core_calls == 1L) "before" else "after"
        )
      )
    },
    .sirap_cli_demux = function(config, rad, capture) {
      write_test_fastq_gz(
        file.path(config$out_dir, paste0(config$output, ".fq.gz"))
      )
      list(status = 0L, stdout = "", stderr = "", command = "rad")
    },
    .package = "SIRAP"
  )

  expect_error(
    rad_demux(
      "five_prime", fastq, out_dir = out_dir,
      backend = "cli", rad = "/mock/rad"
    ),
    "identity changed",
    class = "sirap_rad_cli_compatibility_error"
  )
})

test_that("rad_demux rejects mutually ignored whitelist settings", {
  for (backend in c("native", "cli")) {
    expect_error(
      rad_demux(
        "five_prime", "reads.fq.gz",
        kit = "10x_5v1", global_whitelist = "global.csv",
        backend = backend
      ),
      "alternative global/reference",
      class = "sirap_rad_configuration_error",
      info = backend
    )

    expect_error(
      rad_demux(
        "five_prime", "reads.fq.gz",
        custom_whitelist = "custom.csv", auto_whitelist = TRUE,
        backend = backend
      ),
      "true/sample barcode set",
      class = "sirap_rad_configuration_error",
      info = backend
    )
  }
})

test_that("CLI failures cannot be normalized as successful runs", {
  config <- list(
    out_dir = withr::local_tempdir(pattern = "sirap-rad-cli-failure-"),
    output = "demo", fastq = "reads.fq.gz", layout = "five_prime",
    auto_whitelist = FALSE, write_debug = FALSE
  )
  before <- SIRAP:::.sirap_snapshot_cli_outputs(config)
  expect_error(
    SIRAP:::.sirap_normalize_cli_result(
      list(status = 3L),
      config,
      list(core_version = "legacy"),
      before
    ),
    "status 3",
    class = "sirap_rad_run_error"
  )

  for (status in list("0", 0.5, NA_real_, logical())) {
    expect_error(
      SIRAP:::.sirap_normalize_cli_result(
        list(status = status), config, list(core_version = "legacy"), before
      ),
      class = "sirap_rad_run_error"
    )
  }
})

test_that("CLI success requires a newly produced primary artifact", {
  out_dir <- withr::local_tempdir(pattern = "sirap-rad-cli-artifacts-")
  config <- list(
    out_dir = out_dir, output = "demo", fastq = "reads.fq.gz",
    layout = "five_prime", auto_whitelist = FALSE, write_debug = FALSE
  )
  stale_primary <- file.path(out_dir, "demo.fq.gz")
  unrelated <- file.path(out_dir, "demo_old_notes.txt")
  writeBin(charToRaw("stale"), stale_primary)
  writeLines("unrelated", unrelated)
  before <- SIRAP:::.sirap_snapshot_cli_outputs(config)
  successful_status <- list(
    status = 0L, stdout = "", stderr = "", command = "rad"
  )

  expect_error(
    SIRAP:::.sirap_normalize_cli_result(
      successful_status, config, list(core_version = "legacy"), before
    ),
    "did not create or update",
    class = "sirap_rad_run_error"
  )

  Sys.sleep(0.01)
  write_test_fastq_gz(stale_primary)
  result <- SIRAP:::.sirap_normalize_cli_result(
    successful_status, config, list(core_version = "legacy"), before
  )
  expect_identical(result$artifacts$fastq, normalizePath(stale_primary))
  expect_false(unrelated %in% result$artifacts$files)
})

test_that("CLI success cannot normalize a corrupt gzip artifact", {
  out_dir <- withr::local_tempdir(pattern = "sirap-rad-cli-corrupt-output-")
  config <- list(
    out_dir = out_dir, output = "demo", fastq = "reads.fq.gz",
    layout = "five_prime", auto_whitelist = FALSE, write_debug = FALSE,
    reuse_cache = FALSE
  )
  before <- SIRAP:::.sirap_snapshot_cli_outputs(config)
  primary <- file.path(out_dir, "demo.fq.gz")
  write_test_fastq_gz(primary)
  bytes <- readBin(primary, "raw", n = file.info(primary)$size)
  writeBin(bytes[seq_len(length(bytes) - 5L)], primary)

  expect_error(
    SIRAP:::.sirap_normalize_cli_result(
      list(status = 0L, stdout = "", stderr = "", command = "rad"),
      config, list(core_version = "1.0.2"), before
    ),
    "Invalid RAD CLI primary output",
    class = "sirap_rad_run_error"
  )
})

test_that("sirap RAD results print their handshake", {
  core <- mock_rad_core()
  result <- mock_rad_result(core)
  result <- SIRAP:::.sirap_finish_rad_result(result, "native", core, Sys.time())

  expect_output(print(result), "core API 1")
  expect_output(print(result), "Backend: native")
  expect_output(print(result), "Engine:  embedded")
})
