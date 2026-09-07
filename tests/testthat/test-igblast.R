sirap_igblast_test_database <- function(directory,
                                        receptor = "IG",
                                        include_c = FALSE) {
  segments <- c("V", "D", "J")
  if (include_c) {
    segments <- c(segments, "C")
  }
  prefixes <- stats::setNames(
    file.path(directory, paste0(receptor, "_", segments)),
    tolower(segments)
  )
  for (prefix in prefixes) {
    for (path in paste0(prefix, c(".nhr", ".nin", ".nsq"))) {
      writeBin(as.raw(1L), path)
      expect_gt(file.info(path)$size, 0)
    }
  }
  prefixes
}

sirap_igblast_test_executable <- function(
    directory,
    airr_lines = c(
      "sequence_id\tlocus\tv_call\td_call\tj_call",
      "stub-read\tIGH\tIGHV1\tIGHD1\tIGHJ1"
    ),
    name = "igblastn-stub",
    exit_status = 0L,
    write_output = TRUE) {
  executable <- file.path(directory, name)
  script <- c(
    "#!/bin/sh",
    "out=''",
    "while [ \"$#\" -gt 0 ]; do",
    "  if [ \"$1\" = '-out' ]; then out=$2; shift 2; else shift; fi",
    "done",
    "printf '%s\\n' 'stub IgBLAST stdout'",
    "printf '%s\\n' 'stub IgBLAST stderr' >&2",
    "printf 'IGDATA:%s\\n' \"$IGDATA\""
  )
  if (write_output) {
    quoted_lines <- paste(vapply(airr_lines, shQuote, character(1)), collapse = " ")
    script <- c(script, paste0("printf '%s\\n' ", quoted_lines, " > \"$out\""))
  }
  script <- c(script, paste("exit", as.integer(exit_status)))
  writeLines(script, executable)
  Sys.chmod(executable, mode = "0755")
  normalizePath(executable, mustWork = TRUE)
}

sirap_igblast_query_aware_executable <- function(
    directory,
    name = "igblastn-query-aware",
    fail_on = NULL) {
  executable <- file.path(directory, name)
  script <- c(
    "#!/bin/sh",
    "query=''",
    "out=''",
    "threads=''",
    "while [ \"$#\" -gt 0 ]; do",
    "  case \"$1\" in",
    "    -query) query=$2; shift 2 ;;",
    "    -out) out=$2; shift 2 ;;",
    "    -num_threads) threads=$2; shift 2 ;;",
    "    *) shift ;;",
    "  esac",
    "done"
  )
  if (!is.null(fail_on)) {
    script <- c(script, paste0(
      "if grep -q '^>", fail_on, "' \"$query\"; then ",
      "printf '%s\\n' 'deliberate chunk failure' >&2; exit 7; fi"
    ))
  }
  script <- c(
    script,
    "printf 'query:%s\\nthreads:%s\\n' \"$query\" \"$threads\"",
    "printf '%s\\n' 'sequence_id\tlocus\tv_call\td_call\tj_call' > \"$out\"",
    paste0(
      "awk 'BEGIN { OFS=\"\\t\" } /^>/ { id=substr($0, 2); ",
      "sub(/[ \\t].*$/, \"\", id); print id, \"IGH\", \"IGHV1\", ",
      "\"IGHD1\", \"IGHJ1\" }' \"$query\" >> \"$out\""
    ),
    "exit 0"
  )
  writeLines(script, executable)
  Sys.chmod(executable, mode = "0755")
  normalizePath(executable, mustWork = TRUE)
}

sirap_igblast_test_installation <- function(
    directory,
    organism = "human",
    receptor = "ig",
    include_c = TRUE,
    airr_lines = c(
      "sequence_id\tlocus\tv_call\td_call\tj_call",
      "stub-read\tIGH\tIGHV1\tIGHD1\tIGHJ1"
    )) {
  base <- file.path(directory, "igblast installation")
  database <- file.path(base, "database")
  optional <- file.path(base, "optional_file")
  bin <- file.path(base, "bin")
  expect_true(dir.create(database, recursive = TRUE))
  expect_true(dir.create(optional))
  expect_true(dir.create(bin))

  segments <- c("v", "d", "j")
  if (include_c) {
    segments <- c(segments, "c")
  }
  prefixes <- stats::setNames(
    file.path(
      database,
      paste0("imgt_", organism, "_", receptor, "_", segments)
    ),
    segments
  )
  for (prefix in prefixes) {
    for (path in paste0(prefix, c(".nhr", ".nin", ".nsq"))) {
      writeBin(as.raw(1L), path)
    }
  }
  auxiliary <- file.path(optional, paste0(organism, "_gl.aux"))
  writeLines("auxiliary fixture", auxiliary)
  executable <- sirap_igblast_test_executable(
    bin,
    airr_lines = airr_lines,
    name = "igblastn"
  )
  list(
    base = normalizePath(base),
    databases = prefixes,
    auxiliary = normalizePath(auxiliary),
    executable = executable
  )
}

sirap_igblast_arg <- function(args, flag) {
  position <- match(flag, args)
  if (is.na(position) || position == length(args)) {
    return(NULL)
  }
  args[[position + 1L]]
}

test_that("IgBLAST database prefixes are resolved and validated", {
  directory <- withr::local_tempdir(pattern = "sirap-igblast-db-")
  prefixes <- sirap_igblast_test_database(directory, include_c = TRUE)

  dbs <- .sirap_igblast_db_paths(directory, receptor = "ig")
  resolved_prefix <- function(prefix) {
    file.path(normalizePath(dirname(prefix)), basename(prefix))
  }
  expect_identical(dbs$v, resolved_prefix(prefixes[["v"]]))
  expect_identical(dbs$d, resolved_prefix(prefixes[["d"]]))
  expect_identical(dbs$j, resolved_prefix(prefixes[["j"]]))
  expect_identical(dbs$c, resolved_prefix(prefixes[["c"]]))

  tr_prefixes <- sirap_igblast_test_database(directory, receptor = "TR")
  tr_dbs <- .sirap_igblast_db_paths(directory, receptor = "tr")
  expect_identical(tr_dbs$v, resolved_prefix(tr_prefixes[["v"]]))
  expect_identical(tr_dbs$d, resolved_prefix(tr_prefixes[["d"]]))
  expect_identical(tr_dbs$j, resolved_prefix(tr_prefixes[["j"]]))
  expect_null(tr_dbs$c)

  custom_v <- file.path(directory, "custom-human-V")
  for (path in paste0(custom_v, c(".nhr", ".nin", ".nsq"))) {
    writeBin(as.raw(1L), path)
  }
  overridden <- .sirap_igblast_db_paths(
    directory,
    receptor = "ig",
    germline_db_v = custom_v
  )
  expect_identical(overridden$v, resolved_prefix(custom_v))
  expect_identical(overridden$d, resolved_prefix(prefixes[["d"]]))
  expect_error(
    .sirap_validate_igblast_db_prefix(
      paste0(prefixes[["v"]], ".nhr"), "germline_db_v"
    ),
    "prefix, without a database-file extension",
    class = "sirap_igblast_error"
  )

  incomplete_volume <- file.path(directory, "incomplete-volume")
  for (path in paste0(incomplete_volume, ".00", c(".nhr", ".nin", ".nsq"))) {
    writeBin(as.raw(1L), path)
  }
  writeBin(as.raw(1L), paste0(incomplete_volume, ".01.nhr"))
  expect_false(.sirap_igblast_db_prefix_exists(incomplete_volume))
  empty_alias <- file.path(directory, "empty-alias")
  expect_true(file.create(paste0(empty_alias, ".nal")))
  expect_false(.sirap_igblast_db_prefix_exists(empty_alias))

  unlink(paste0(prefixes[["d"]], ".nsq"))
  expect_error(
    .sirap_igblast_db_paths(directory, receptor = "ig"),
    "germline_db_d.*prefix",
    class = "sirap_igblast_error"
  )
})

test_that("Change-O IgBLAST installation paths resolve by organism and receptor", {
  directory <- withr::local_tempdir(pattern = "sirap-igblast-install-")
  installation <- sirap_igblast_test_installation(
    directory,
    organism = "mouse",
    receptor = "tr",
    include_c = FALSE
  )
  dbs <- .sirap_igblast_db_paths(
    receptor = "tr",
    igblast_base = installation$base,
    organism = "mouse"
  )
  expected <- lapply(installation$databases, function(prefix) {
    file.path(normalizePath(dirname(prefix)), basename(prefix))
  })
  expect_identical(dbs$v, expected$v)
  expect_identical(dbs$d, expected$d)
  expect_identical(dbs$j, expected$j)
  expect_null(dbs$c)
  expect_identical(
    .sirap_igblast_aux_file(NULL, installation$base, "mouse"),
    installation$auxiliary
  )

  override_aux <- file.path(directory, "explicit auxiliary.aux")
  writeLines("explicit", override_aux)
  expect_identical(
    .sirap_igblast_aux_file(override_aux, installation$base, "mouse"),
    normalizePath(override_aux)
  )
})

test_that("IgBLAST arguments enforce AIRR output and preserve extra options", {
  dbs <- list(v = "/db/V", d = "/db/D", j = "/db/J", c = "/db/C")
  args <- .sirap_igblast_args(
    input_fasta = "/reads/query.fasta",
    output_tsv = "/results/query.tsv",
    receptor = "ig",
    organism = "mouse",
    dbs = dbs,
    aux_file = "/aux/mouse_gl.aux",
    num_threads = 4L,
    extra_args = c("-word_size", "9")
  )

  expect_identical(sirap_igblast_arg(args, "-query"), "/reads/query.fasta")
  expect_identical(sirap_igblast_arg(args, "-ig_seqtype"), "Ig")
  expect_identical(sirap_igblast_arg(args, "-organism"), "mouse")
  expect_identical(sirap_igblast_arg(args, "-germline_db_D"), "/db/D")
  expect_identical(sirap_igblast_arg(args, "-c_region_db"), "/db/C")
  expect_identical(sirap_igblast_arg(args, "-auxiliary_data"), "/aux/mouse_gl.aux")
  expect_identical(sirap_igblast_arg(args, "-domain_system"), "imgt")
  expect_identical(sirap_igblast_arg(args, "-outfmt"), "19")
  expect_identical(sirap_igblast_arg(args, "-num_threads"), "4")
  expect_identical(sirap_igblast_arg(args, "-out"), "/results/query.tsv")
  expect_identical(utils::tail(args, 2L), c("-word_size", "9"))

  tr_args <- .sirap_igblast_args(
    input_fasta = "query.fasta",
    output_tsv = "query.tsv",
    receptor = "tr",
    dbs = list(v = "TR_V", d = "TR_D", j = "TR_J", c = NULL)
  )
  expect_identical(sirap_igblast_arg(tr_args, "-ig_seqtype"), "TCR")
  expect_false("-c_region_db" %in% tr_args)

  expect_error(
    .sirap_igblast_args(
      "query.fasta", "query.tsv", dbs = dbs,
      extra_args = c("-outfmt", "7")
    ),
    "cannot override.*-outfmt",
    class = "sirap_igblast_error"
  )
})

test_that("IgBLAST FASTA chunks preserve records and balance sequence bases", {
  directory <- withr::local_tempdir(pattern = "sirap-igblast-chunks-")
  fasta_path <- file.path(directory, "wrapped.fasta")
  fasta_lines <- c(
    ">read-1 description", "AAAA", "AAAA",
    ">read-2", "CC",
    ">read-3", "GGGG", "GG",
    ">read-4", "TTTT",
    ">read-5", "NNNN", "NNNN",
    ">read-6", "ACGT"
  )
  writeLines(fasta_lines, fasta_path)
  fasta <- .sirap_validate_igblast_fasta(fasta_path)
  automatic_plan <- .sirap_igblast_process_plan(
    fasta, num_threads = 12L, max_processes = 4L
  )
  expect_identical(automatic_plan$desired_chunks, 1L)
  plan <- .sirap_igblast_process_plan(
    fasta, num_threads = 12L, max_processes = 4L, chunk_bases = 11L
  )
  chunks <- .sirap_split_igblast_fasta(
    fasta, file.path(directory, "chunks"), plan
  )

  expect_length(chunks, 3L)
  expect_identical(
    unlist(lapply(chunks, function(chunk) readLines(chunk$path)),
           use.names = FALSE),
    fasta_lines
  )
  expect_equal(sum(vapply(chunks, `[[`, numeric(1), "queries")), 6)
  expect_equal(sum(vapply(chunks, `[[`, numeric(1), "bases")), 32)
  expect_true(all(vapply(chunks, function(chunk) {
    first <- readLines(chunk$path, n = 1L)
    length(first) == 1L && startsWith(first, ">")
  }, logical(1))))

  explicit_plan <- .sirap_igblast_process_plan(
    fasta, num_threads = 12L, max_processes = 4L, chunk_bases = 7L
  )
  expect_equal(explicit_plan$desired_chunks, ceiling(fasta$bases / 7))
  monolithic_plan <- .sirap_igblast_process_plan(
    fasta, num_threads = 12L, max_processes = 1L, chunk_bases = 7L
  )
  expect_identical(monolithic_plan$desired_chunks, 1L)
  expect_error(
    .sirap_igblast_process_plan(fasta, 12L, 0L),
    "max_processes.*positive whole number",
    class = "sirap_igblast_error"
  )
})

test_that("IgBLAST chunk merge rejects schema drift", {
  directory <- withr::local_tempdir(pattern = "sirap-igblast-schema-")
  first <- file.path(directory, "first.tsv")
  second <- file.path(directory, "second.tsv")
  writeLines(c(
    "sequence_id\tlocus\tv_call\td_call\tj_call",
    "read-1\tIGH\tIGHV1\tIGHD1\tIGHJ1"
  ), first)
  writeLines(c(
    "sequence_id\tlocus\tv_call\td_call\tj_call\tproductive",
    "read-2\tIGH\tIGHV1\tIGHD1\tIGHJ1\tT"
  ), second)
  results <- list(
    list(index = 1L, airr = .sirap_inspect_airr_output(first)),
    list(index = 2L, airr = .sirap_inspect_airr_output(second))
  )

  expect_error(
    .sirap_merge_igblast_airr(results, file.path(directory, "merged.tsv")),
    "chunk 2.*schema",
    class = "sirap_igblast_error"
  )
})

test_that("runIgblast performs one FASTA-to-AIRR pass with a stubbed runner", {
  temporary <- withr::local_tempdir(pattern = "sirap-igblast-run-")
  directory <- file.path(temporary, "work with spaces")
  expect_true(dir.create(directory))
  installation <- sirap_igblast_test_installation(
    directory,
    airr_lines = c(
      "sequence_id\tlocus\tv_call\td_call\tj_call",
      "read-1|key=one\tIGH\tIGHV1\tIGHD1\tIGHJ1",
      "read-2|key=two\tIGH\tIGHV2\tIGHD2\tIGHJ2"
    )
  )
  fasta <- file.path(directory, "reads.fasta")
  writeLines(c(
    ">read-1|key=one",
    "ACGT",
    "NN",
    ">read-2|key=two",
    "TTAA"
  ), fasta)
  output <- file.path(directory, "nested", "reads-airr.tsv")
  original_igdata <- Sys.getenv("IGDATA", unset = NA_character_)

  result <- runIgblast(
    input_fasta = fasta,
    output_tsv = output,
    receptor = "ig",
    organism = "human",
    igblast_base = installation$base,
    num_threads = 3L,
    extra_args = c("-word_size", "9")
  )

  expect_s3_class(result, "sirap_igblast_result")
  expect_identical(result$status, "success")
  expect_identical(result$stage, "read")
  expect_identical(result$artifacts$query_fasta, normalizePath(fasta))
  expect_identical(result$artifacts$airr, normalizePath(output))
  expect_identical(
    result$artifacts$log,
    file.path(dirname(normalizePath(output)), "reads-airr.log")
  )
  expect_true(all(file.exists(unlist(result$artifacts, use.names = FALSE))))
  log_text <- paste(readLines(result$artifacts$log), collapse = "\n")
  expect_match(log_text, "stub IgBLAST stdout")
  expect_match(log_text, "stub IgBLAST stderr")
  expect_match(log_text, paste0("IGDATA:", installation$base), fixed = TRUE)
  expect_identical(Sys.getenv("IGDATA", unset = NA_character_), original_igdata)
  expect_identical(sirap_igblast_arg(result$command$args, "-query"), normalizePath(fasta))
  expect_identical(sirap_igblast_arg(result$command$args, "-outfmt"), "19")
  expect_identical(sirap_igblast_arg(result$command$args, "-num_threads"), "3")
  expect_identical(sirap_igblast_arg(result$command$args, "-out"), normalizePath(output))
  expect_false(identical(
    sirap_igblast_arg(result$command$execution_args, "-out"),
    normalizePath(output)
  ))
  expect_identical(result$config$stage, "read")
  expect_identical(result$config$receptor, "ig")
  expect_identical(result$config$organism, "human")
  expect_identical(result$config$igblast_base, installation$base)
  expect_null(result$config$igblast_db_dir)
  expect_identical(result$config$threads, 3L)
  expect_identical(result$config$thread_budget, 3L)
  expect_identical(result$config$processes, 1L)
  expect_identical(result$config$threads_per_process, 3L)
  expect_identical(result$config$executable, installation$executable)
  expect_identical(result$config$auxiliary_data, installation$auxiliary)
  expect_named(
    result$config$databases,
    c("v", "d", "j", "c"),
    ignore.order = FALSE
  )
  expect_identical(result$command$executable, installation$executable)
  expect_length(result$command$executions, 1L)
  expect_identical(result$command$executions[[1L]]$threads, 3L)
  expect_identical(result$command$env, c(IGDATA = installation$base))
  expect_identical(
    sirap_igblast_arg(result$command$execution_args, "-outfmt"),
    "19"
  )
  expect_equal(result$stats$query_count, 2)
  expect_equal(result$stats$queries, 2)
  expect_equal(result$stats$query_bases, 10)
  expect_equal(result$stats$airr_rows, 2)
  expect_equal(result$stats$annotated_queries, 2)
  expect_equal(result$stats$unannotated_queries, 0)
  expect_equal(result$stats$annotation_rate, 1)
  expect_equal(result$stats$chunks, 1)
  expect_equal(result$stats$chunk_queries, 2)
  expect_equal(result$stats$chunk_bases, 10)
  expect_gt(result$stats$airr_bytes, 0)
  expect_gt(result$stats$log_bytes, 0)
  expect_s3_class(result$timing$started, "POSIXct")
  expect_s3_class(result$timing$completed, "POSIXct")
  expect_gte(result$timing$elapsed_seconds, 0)

  airr <- .sirap_read_airr(result$artifacts$airr)
  expect_identical(airr$sequence_id, c("read-1|key=one", "read-2|key=two"))
  expect_identical(airr$v_call, c("IGHV1", "IGHV2"))
})

test_that("runIgblast parallelizes base-balanced chunks and merges in order", {
  skip_on_os("windows")
  skip_on_cran()
  directory <- withr::local_tempdir(pattern = "sirap-igblast-parallel-")
  sirap_igblast_test_database(directory)
  executable <- sirap_igblast_query_aware_executable(directory)
  fasta <- file.path(directory, "unequal reads.fasta")
  identifiers <- sprintf("read-%02d", seq_len(8L))
  lengths <- c(
    40000L, 2000L, 2000L, 2000L, 24000L, 22000L, 20000L, 18000L
  )
  fasta_lines <- unlist(Map(function(identifier, length) {
    c(paste0(">", identifier), paste(rep("A", length), collapse = ""))
  }, identifiers, lengths), use.names = FALSE)
  writeLines(fasta_lines, fasta)

  result <- runIgblast(
    fasta,
    file.path(directory, "parallel.tsv"),
    igblast_db_dir = directory,
    igblast_exec = executable,
    num_threads = 12L,
    max_processes = 4L
  )

  expect_identical(result$config$thread_budget, 12L)
  expect_identical(result$config$max_processes, 4L)
  expect_identical(result$config$processes, 4L)
  expect_identical(result$config$threads_per_process, 3L)
  expect_equal(result$stats$chunks, 4)
  expect_equal(sum(result$stats$chunk_queries), length(identifiers))
  expect_equal(sum(result$stats$chunk_bases), sum(lengths))
  expect_length(result$command$executions, 4L)
  expect_true(all(vapply(result$command$executions, function(command) {
    identical(sirap_igblast_arg(command$args, "-num_threads"), "3")
  }, logical(1))))
  expect_null(result$command$execution_args)

  airr <- data.table::fread(result$artifacts$airr)
  expect_identical(airr$sequence_id, identifiers)
  expect_equal(
    sum(startsWith(readLines(result$artifacts$airr), "sequence_id\t")),
    1
  )
  expect_equal(
    sum(startsWith(readLines(result$artifacts$log), "===== IgBLAST chunk")),
    4
  )
  temporary_queries <- vapply(result$command$executions, function(command) {
    sirap_igblast_arg(command$args, "-query")
  }, character(1))
  expect_false(any(file.exists(temporary_queries)))
})

test_that("explicit databases and auxiliary data override an installation", {
  directory <- withr::local_tempdir(pattern = "sirap-igblast-overrides-")
  installation <- sirap_igblast_test_installation(directory)
  custom_v <- file.path(directory, "custom V database")
  for (path in paste0(custom_v, c(".nhr", ".nin", ".nsq"))) {
    writeBin(as.raw(1L), path)
  }
  custom_aux <- file.path(directory, "custom auxiliary.aux")
  writeLines("custom auxiliary", custom_aux)
  fasta <- file.path(directory, "reads.fasta")
  writeLines(c(">stub-read", "ACGT"), fasta)

  result <- runIgblast(
    fasta,
    file.path(directory, "overridden.tsv"),
    igblast_base = installation$base,
    germline_db_v = custom_v,
    aux_file = custom_aux
  )
  expect_identical(
    result$config$databases$v,
    file.path(normalizePath(dirname(custom_v)), basename(custom_v))
  )
  expect_identical(
    result$config$databases$d,
    file.path(
      installation$base,
      "database",
      "imgt_human_ig_d"
    )
  )
  expect_identical(result$config$auxiliary_data, normalizePath(custom_aux))
  expect_identical(result$command$env, c(IGDATA = installation$base))
})

test_that("runIgblast records an explicit consensus pass", {
  directory <- withr::local_tempdir(pattern = "sirap-igblast-stage-")
  sirap_igblast_test_database(directory)
  executable <- sirap_igblast_test_executable(
    directory,
    airr_lines = c(
      "sequence_id\tlocus\tv_call\td_call\tj_call",
      "consensus-1\tIGH\tIGHV1\tIGHD1\tIGHJ1"
    )
  )
  fasta <- file.path(directory, "consensus.fa")
  writeLines(c(">consensus-1", "ACGT"), fasta)

  result <- runIgblast(
    fasta,
    file.path(directory, "consensus.tsv"),
    igblast_db_dir = directory,
    igblast_exec = executable,
    stage = "consensus"
  )
  expect_identical(result$stage, "consensus")
  expect_identical(result$config$stage, "consensus")
})

test_that("runIgblast never accepts stale or malformed staged output", {
  directory <- withr::local_tempdir(pattern = "sirap-igblast-staging-")
  sirap_igblast_test_database(directory)
  fasta <- file.path(directory, "reads.fasta")
  writeLines(c(">new-read", "ACGT"), fasta)

  stale <- file.path(directory, "stale.tsv")
  stale_lines <- c(
    "sequence_id\tlocus\tv_call\td_call\tj_call",
    "OLD\tIGH\tIGHV1\tIGHD1\tIGHJ1"
  )
  writeLines(stale_lines, stale)
  no_output <- sirap_igblast_test_executable(
    directory,
    name = "igblastn-no-output",
    write_output = FALSE
  )
  expect_error(
    runIgblast(
      fasta, stale,
      igblast_db_dir = directory,
      igblast_exec = no_output
    ),
    "does not identify an existing file",
    class = "sirap_igblast_error"
  )
  expect_identical(readLines(stale), stale_lines)

  replacement <- c(
    "sequence_id\tlocus\tv_call\td_call\tj_call",
    "new-read\tIGH\tIGHV2\tIGHD2\tIGHJ2"
  )
  valid <- sirap_igblast_test_executable(
    directory,
    airr_lines = replacement,
    name = "igblastn-valid-output"
  )
  replaced <- runIgblast(
    fasta, stale,
    igblast_db_dir = directory,
    igblast_exec = valid
  )
  expect_identical(readLines(replaced$artifacts$airr), replacement)

  malformed_output <- sirap_igblast_test_executable(
    directory,
    airr_lines = c("sequence_id\tlocus", "new-read\tIGH"),
    name = "igblastn-malformed-output"
  )
  malformed_path <- file.path(directory, "malformed-output.tsv")
  expect_error(
    runIgblast(
      fasta, malformed_path,
      igblast_db_dir = directory,
      igblast_exec = malformed_output
    ),
    "missing required AIRR fields",
    class = "sirap_igblast_error"
  )
  expect_false(file.exists(malformed_path))

  malformed_row <- sirap_igblast_test_executable(
    directory,
    airr_lines = c(
      "sequence_id\tlocus\tv_call\td_call\tj_call",
      "new-read\tIGH"
    ),
    name = "igblastn-malformed-row"
  )
  malformed_row_path <- file.path(directory, "malformed-row.tsv")
  expect_error(
    runIgblast(
      fasta, malformed_row_path,
      igblast_db_dir = directory,
      igblast_exec = malformed_row
    ),
    "AIRR data row 1 has 2 fields; the header has 5",
    class = "sirap_igblast_error"
  )
  expect_false(file.exists(malformed_row_path))

  no_hit <- sirap_igblast_test_executable(
    directory,
    airr_lines = "sequence_id\tlocus\tv_call\td_call\tj_call",
    name = "igblastn-no-hit"
  )
  no_hit_result <- runIgblast(
    fasta, file.path(directory, "no-hit.tsv"),
    igblast_db_dir = directory,
    igblast_exec = no_hit,
    stage = "read"
  )
  expect_equal(no_hit_result$stats$query_count, 1)
  expect_equal(no_hit_result$stats$airr_rows, 0)
  expect_equal(no_hit_result$stats$annotated_queries, 0)
  expect_equal(no_hit_result$stats$unannotated_queries, 1)
  expect_equal(no_hit_result$stats$annotation_rate, 0)

  blank_call <- sirap_igblast_test_executable(
    directory,
    airr_lines = c(
      "sequence_id\tlocus\tv_call\td_call\tj_call",
      "new-read\t\t\t\t"
    ),
    name = "igblastn-blank-call"
  )
  blank_call_result <- runIgblast(
    fasta, file.path(directory, "blank-call.tsv"),
    igblast_db_dir = directory,
    igblast_exec = blank_call,
    stage = "read"
  )
  expect_equal(blank_call_result$stats$airr_rows, 1)
  expect_equal(blank_call_result$stats$annotated_queries, 0)
  expect_equal(blank_call_result$stats$unannotated_queries, 1)

  too_many <- sirap_igblast_test_executable(
    directory,
    airr_lines = c(
      "sequence_id\tlocus\tv_call\td_call\tj_call",
      "new-read\tIGH\tIGHV1\tIGHD1\tIGHJ1",
      "impossible-extra\tIGH\tIGHV2\tIGHD2\tIGHJ2"
    ),
    name = "igblastn-too-many"
  )
  too_many_path <- file.path(directory, "too-many.tsv")
  expect_error(
    runIgblast(
      fasta, too_many_path,
      igblast_db_dir = directory,
      igblast_exec = too_many
    ),
    "more data rows than input queries",
    class = "sirap_igblast_error"
  )
  expect_false(file.exists(too_many_path))
})

test_that("a failed parallel IgBLAST chunk publishes no partial artifacts", {
  skip_on_os("windows")
  directory <- withr::local_tempdir(pattern = "sirap-igblast-failed-chunk-")
  sirap_igblast_test_database(directory)
  executable <- sirap_igblast_query_aware_executable(
    directory, name = "igblastn-failed-chunk", fail_on = "read-3"
  )
  fasta <- file.path(directory, "reads.fasta")
  writeLines(c(
    ">read-1", paste(rep("A", 100L), collapse = ""),
    ">read-2", paste(rep("C", 100L), collapse = ""),
    ">read-3", paste(rep("G", 100L), collapse = ""),
    ">read-4", paste(rep("T", 100L), collapse = "")
  ), fasta)
  output <- file.path(directory, "existing.tsv")
  log <- file.path(directory, "existing.log")
  original_output <- c("old header", "old result")
  original_log <- "old log"
  writeLines(original_output, output)
  writeLines(original_log, log)

  expect_error(
    runIgblast(
      fasta, output,
      igblast_db_dir = directory,
      igblast_exec = executable,
      num_threads = 4L,
      max_processes = 2L,
      chunk_bases = 200L
    ),
    "IgBLAST chunk 2/2 failed.*status 7",
    class = "sirap_igblast_command_error"
  )
  expect_identical(readLines(output), original_output)
  expect_identical(readLines(log), original_log)
})

test_that("runIgblast rejects output and log symlinks to its FASTA input", {
  skip_on_os("windows")
  directory <- withr::local_tempdir(pattern = "sirap-igblast-collision-")
  fasta <- file.path(directory, "reads.fasta")
  original <- c(">read-1", "ACGT")
  writeLines(original, fasta)

  output_alias <- file.path(directory, "output-alias.tsv")
  expect_true(file.symlink(fasta, output_alias))
  expect_error(
    runIgblast(fasta, output_alias),
    "different from `input_fasta`",
    class = "sirap_igblast_error"
  )
  expect_identical(readLines(fasta), original)

  output <- file.path(directory, "log-alias.tsv")
  log_alias <- file.path(directory, "log-alias.log")
  expect_true(file.symlink(fasta, log_alias))
  expect_error(
    runIgblast(fasta, output),
    "different from `input_fasta`",
    class = "sirap_igblast_error"
  )
  expect_identical(readLines(fasta), original)
})

test_that("the external command runner quotes arguments and reports failures", {
  skip_on_os("windows")
  temporary <- withr::local_tempdir(pattern = "sirap-command-")
  directory <- file.path(temporary, "runner with spaces")
  expect_true(dir.create(directory))
  executable <- file.path(directory, "echo argument")
  writeLines(c(
    "#!/bin/sh",
    "printf '%s\\n' \"$1\"",
    "printf '%s\\n' \"$IGDATA\""
  ), executable)
  Sys.chmod(executable, "0755")
  marker <- file.path(temporary, "must-not-exist")
  env_marker <- file.path(temporary, "env-must-not-exist")
  payload <- paste0("alpha ' beta $(touch ", shQuote(marker), ")")
  env_value <- paste0(
    file.path(directory, "IGDATA with spaces"),
    " $(touch ", shQuote(env_marker), ")"
  )
  log <- file.path(directory, "quoted output.log")
  original_igdata <- Sys.getenv("IGDATA", unset = NA_character_)

  expect_identical(
    .sirap_run_command(
      executable, payload, log = log, env = c(IGDATA = env_value)
    ),
    0L
  )
  expect_identical(readLines(log), c(payload, env_value))
  expect_false(file.exists(marker))
  expect_false(file.exists(env_marker))
  expect_identical(Sys.getenv("IGDATA", unset = NA_character_), original_igdata)

  failure <- sirap_igblast_test_executable(
    directory,
    name = "igblastn-fail",
    exit_status = 7L,
    write_output = FALSE
  )
  failure_log <- file.path(directory, "failure.log")
  expect_error(
    .sirap_run_command(failure, character(), log = failure_log),
    "External command failed with status 7.*igblastn-fail",
    class = "sirap_igblast_command_error"
  )
  expect_match(paste(readLines(failure_log), collapse = "\n"), "stub IgBLAST stderr")
})

test_that("runIgblast rejects FASTQ, gzip, malformed FASTA, and bad threads", {
  directory <- withr::local_tempdir(pattern = "sirap-igblast-input-")
  fastq <- file.path(directory, "reads.fastq")
  writeLines(c("@read-1", "ACGT", "+", "IIII"), fastq)
  expect_error(
    .sirap_validate_igblast_fasta(fastq),
    "does not begin with a FASTA header",
    class = "sirap_igblast_error"
  )
  expect_error(
    runIgblast(fastq, file.path(directory, "fastq.tsv")),
    "does not begin with a FASTA header",
    class = "sirap_igblast_error"
  )

  extensionless <- file.path(directory, "extensionless")
  writeLines(c(">read-1", "ACGT"), extensionless)
  expect_equal(.sirap_validate_igblast_fasta(extensionless)$queries, 1)
  ffn <- file.path(directory, "coding.ffn")
  writeLines(c(">read-1", "ACGT"), ffn)
  expect_equal(.sirap_validate_igblast_fasta(ffn)$queries, 1)

  compressed <- file.path(directory, "compressed.fasta")
  connection <- gzfile(compressed, open = "wt")
  writeLines(c(">read-1", "ACGT"), connection)
  close(connection)
  expect_error(
    .sirap_validate_igblast_fasta(compressed),
    "gzip-compressed",
    class = "sirap_igblast_error"
  )
  expect_error(
    runIgblast(compressed, file.path(directory, "gzip.tsv")),
    "gzip-compressed",
    class = "sirap_igblast_error"
  )

  malformed <- file.path(directory, "malformed.fasta")
  writeLines(c(">read-1", "ACGT", ">read-2"), malformed)
  expect_error(
    .sirap_validate_igblast_fasta(malformed),
    "has no sequence",
    class = "sirap_igblast_error"
  )
  expect_error(
    runIgblast(malformed, file.path(directory, "malformed.tsv")),
    "has no sequence",
    class = "sirap_igblast_error"
  )
  expect_error(
    .sirap_igblast_positive_integer(0, "num_threads"),
    "positive whole number",
    class = "sirap_igblast_error"
  )
  expect_error(
    .sirap_igblast_positive_integer(1.5, "num_threads"),
    "positive whole number",
    class = "sirap_igblast_error"
  )
  expect_error(
    runIgblast(
      malformed, file.path(directory, "processes.tsv"),
      max_processes = 0
    ),
    "max_processes.*positive whole number",
    class = "sirap_igblast_error"
  )
  expect_error(
    runIgblast(
      malformed, file.path(directory, "chunks.tsv"),
      chunk_bases = 1.5
    ),
    "chunk_bases.*positive whole number",
    class = "sirap_igblast_error"
  )
  expect_error(
    runIgblast(malformed, file.path(directory, "threads.tsv"), num_threads = 0),
    "positive whole number",
    class = "sirap_igblast_error"
  )
  expect_error(
    runIgblast(malformed, file.path(directory, "stage.tsv"), stage = "molecule"),
    "should be one of"
  )
})

test_that("IgBLAST executable validation returns a resolved path", {
  directory <- withr::local_tempdir(pattern = "sirap-igblast-tools-")
  executable <- sirap_igblast_test_executable(directory)
  tools <- .sirap_check_tools(c(igblastn = executable))
  expect_identical(unname(tools), executable)
  expect_identical(names(tools), "igblastn")

  path_executable <- sirap_igblast_test_executable(
    directory, name = "igblastn-on-path"
  )
  withr::local_envvar(c(
    PATH = paste(directory, Sys.getenv("PATH"), sep = .Platform$path.sep)
  ))
  path_tool <- .sirap_check_tools(c(igblastn = "igblastn-on-path"))
  expect_identical(unname(path_tool), path_executable)

  not_executable <- file.path(directory, "not-executable")
  writeLines("not executable", not_executable)
  Sys.chmod(not_executable, mode = "0644")
  expect_error(
    .sirap_check_tools(c(igblastn = not_executable)),
    "not found or not executable",
    class = "sirap_igblast_tool_error"
  )
  expect_error(
    .sirap_check_tools(c(igblastn = file.path(directory, "missing"))),
    "not found or not executable",
    class = "sirap_igblast_tool_error"
  )
})
