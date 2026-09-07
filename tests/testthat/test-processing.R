test_that("gap length helpers work", {
  expect_equal(.sirap_max_gap_length("AA...TT..G"), 3L)
  expect_equal(.sirap_opening_gap_length("...A.."), 3L)
  expect_equal(.sirap_opening_gap_length("A.."), 0L)
})

test_that("filterAirr filters by locus and lengths", {
  airr <- data.frame(
    locus = c("IGH", "IGK", "IGH"),
    v_sequence_start = c(1, 1, 1),
    v_sequence_end = c(120, 120, 50),
    j_sequence_start = c(1, 1, 1),
    j_sequence_end = c(20, 20, 20),
    junction_length = c(12, 12, 12),
    productive = c(TRUE, TRUE, TRUE)
  )
  out <- filterAirr(airr, locus_val = "IGH", productive_only = TRUE)
  expect_equal(nrow(out), 1L)
  expect_equal(out$locus, "IGH")
})

test_that("filterAirr can require complete per-read VDJ assignments", {
  airr <- data.frame(
    v_sequence_start = rep(1, 3),
    v_sequence_end = rep(120, 3),
    j_sequence_start = rep(1, 3),
    j_sequence_end = rep(20, 3),
    junction_length = rep(12, 3),
    productive = c("T", "T", "F"),
    complete_vdj = c("T", "F", "T")
  )

  out <- filterAirr(
    airr, productive_only = TRUE, complete_vdj = TRUE,
    min_consensus_count = NULL
  )

  expect_equal(nrow(out), 1L)
  expect_identical(out$complete_vdj, "T")

  expect_error(
    filterAirr(
      airr[, setdiff(names(airr), "complete_vdj")],
      productive_only = TRUE,
      complete_vdj = TRUE,
      min_consensus_count = NULL
    ),
    "complete_vdj"
  )
})

test_that("filterAirr uses canonical pRESTO spatial identifiers", {
  airr <- data.frame(
    sequence_id = c("read-1", "read-2"),
    barcode = c("AACC", "TTGA"),
    v_sequence_start = c(1, 1),
    v_sequence_end = c(120, 120),
    j_sequence_start = c(1, 1),
    j_sequence_end = c(20, 20),
    junction_length = c(12, 12)
  )

  out <- filterAirr(
    airr,
    productive_only = FALSE,
    min_consensus_count = NULL,
    spatial_id = "barcode",
    add_spatial_sequence_id = TRUE,
    spatial_id_key = "spot"
  )

  expect_identical(
    out$spatial_sequence_id,
    c("read-1|SPOT=AACC", "read-2|SPOT=TTGA")
  )
})

test_that("pRESTO commands keep consensus IgBLAST distinct and copy barcode", {
  paths <- .presto_paths(
    "/tmp/path with spaces/reads.fasta", 0.975, 0.51,
    "/tmp/path with spaces"
  )
  commands <- .sirap_presto_commands(
    paths = paths,
    receptor = "ig",
    organism = "human",
    identity = 0.975,
    consensus_freq = 0.51,
    germline_files = c("/db/IGHV.fa", "/db/IGHD.fa", "/db/IGHJ.fa"),
    igblast_base = "/opt/igblast",
    igblast_exec = "/opt/igblast/bin/igblastn",
    muscle_exec = "/opt/muscle",
    nproc = 4L,
    cluster_field = "BARCODE",
    delimiters = c("|", "=", ","),
    parse_fields = c("BARCODE", "CLUSTER"),
    parse_key = "BARCODE_CLUSTER",
    parse_act = "cat",
    consensus_min_count = 2L,
    consensus_copy_fields = "BARCODE",
    consensus_copy_action = "majority",
    align_log = NULL,
    consensus_log = NULL,
    igblast_dbs = list(
      v = "/opt/igblast/database/custom_v",
      d = "/opt/igblast/database/custom_d",
      j = "/opt/igblast/database/custom_j",
      c = NULL
    )
  )

  expect_identical(
    vapply(commands, `[[`, character(1), "stage"),
    c(
      "cluster", "merge_headers", "align", "consensus",
      "consensus_igblast", "make_db"
    )
  )
  expect_identical(
    vapply(commands, `[[`, character(1), "output"),
    unname(unlist(paths[c(
      "cluster", "reheader", "aligned", "consensus", "fmt7", "db_pass"
    )]))
  )
  consensus <- commands[[4L]]$args
  expect_identical(consensus[match("--cf", consensus) + 1L], "BARCODE")
  second_pass <- commands[[5L]]$args
  expect_identical(
    second_pass[match("-s", second_pass) + 1L],
    paths$consensus
  )
  expect_false(paths$input %in% second_pass)
  expect_identical(
    second_pass[match("--vdb", second_pass) + 1L],
    "/opt/igblast/database/custom_v"
  )
})

test_that("bare FASTA defaults to the mandatory read-level IgBLAST pass", {
  expect_identical(
    .sirap_process_input_type("reads.fasta", "auto"),
    "rad_fasta"
  )
  expect_identical(
    .sirap_process_input_type("reads.fasta", "spatial_fasta"),
    "spatial_fasta"
  )
})

test_that("compressed RAD FASTA is normalized before IgBLAST", {
  work <- withr::local_tempdir(pattern = "sirap-rad-fasta-")
  input <- file.path(work, "reads.fasta.gz")
  connection <- gzfile(input, open = "wt")
  writeLines(c(">read-F|BARCODE=AACC|UMI=TGCA", "ACGT"), connection)
  close(connection)

  prepared <- .sirap_prepare_rad_fasta(
    input, work, sample_name = "reads", threads = 1L
  )

  expect_s3_class(prepared$result, "sirap_rad_reformat_result")
  expect_false(grepl("\\.gz$", prepared$path))
  expect_identical(
    readLines(prepared$path, warn = FALSE),
    c(">read-F|BARCODE=AACC|UMI=TGCA", "ACGT")
  )
})

test_that("RAD FASTA SAM tags become pRESTO fields before IgBLAST", {
  work <- withr::local_tempdir(pattern = "sirap-rad-fasta-tags-")
  input <- file.path(work, "reads.fasta")
  writeLines(c(">read-F CB:Z:AACC UB:Z:TGCA", "ACGT"), input)

  prepared <- .sirap_prepare_rad_fasta(
    input, work, sample_name = "reads", threads = 1L
  )

  expect_s3_class(prepared$result, "sirap_rad_reformat_result")
  expect_identical(
    readLines(prepared$path, warn = FALSE),
    c(">read-F|BARCODE=AACC|UMI=TGCA", "ACGT")
  )
})

test_that("processSpatialSample runs read and consensus IgBLAST in order", {
  work <- withr::local_tempdir(pattern = "sirap-two-pass-")
  input <- file.path(work, "rad reads.fastq.gz")
  writeLines("stub", input)
  germline <- file.path(work, "imgt.fa")
  writeLines(c(">IGHV1", "ACGT"), germline)
  igblast_root <- file.path(work, "mock igblast")
  dir.create(igblast_root)
  calls <- character()
  observed <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    rad_reformat = function(fastq, output_fastq, threads, to_fasta,
                            header_format, ...) {
      calls <<- c(calls, "rad_reformat")
      expect_identical(fastq, normalizePath(input))
      expect_true(to_fasta)
      expect_identical(header_format, "presto")
      writeLines(c(
        ">read-F|BARCODE=AACC|UMI=TGCA", "ACGTACGT",
        ">read-R|BARCODE=TTGA", "TTTTAAAA"
      ), output_fastq)
      structure(
        list(artifacts = list(fasta = normalizePath(output_fastq))),
        class = c("sirap_rad_reformat_result", "rad_reformat_result", "list")
      )
    },
    runIgblast = function(input_fasta, output_tsv, stage, igblast_base,
                          num_threads, max_processes, ...) {
      calls <<- c(calls, paste0(stage, "_igblast"))
      expect_identical(stage, "read")
      expect_identical(igblast_base, normalizePath(igblast_root))
      expect_identical(num_threads, 1L)
      expect_identical(max_processes, 4L)
      expect_match(readLines(input_fasta, n = 1L), "BARCODE=AACC")
      data.table::fwrite(data.frame(
        sequence_id = c(
          "read-F|BARCODE=AACC|UMI=TGCA",
          "read-R|BARCODE=TTGA"
        ),
        locus = c("IGH", "IGH"),
        v_call = c("IGHV1", "IGHV2"),
        d_call = c("IGHD1", "IGHD2"),
        j_call = c("IGHJ1", "IGHJ2"),
        v_sequence_start = c(1L, 1L),
        v_sequence_end = c(120L, 120L),
        j_sequence_start = c(1L, 1L),
        j_sequence_end = c(20L, 20L),
        junction_length = c(12L, 12L),
        productive = c("T", "T"),
        complete_vdj = c("T", "T"),
        sequence_alignment = c("ACGTACGT", "TTTTAAAA")
      ), output_tsv, sep = "\t")
      structure(
        list(
          stage = stage,
          artifacts = list(airr = normalizePath(output_tsv)),
          config = list(
            igblast_base = normalizePath(igblast_root),
            databases = list(
              v = file.path(igblast_root, "database/imgt_human_ig_v"),
              d = file.path(igblast_root, "database/imgt_human_ig_d"),
              j = file.path(igblast_root, "database/imgt_human_ig_j"),
              c = NULL
            )
          )
        ),
        class = c("sirap_igblast_result", "list")
      )
    },
    .sirap_check_tools = function(execs) execs,
    .sirap_run_command = function(exec, args, log = NULL) {
      stage <- switch(
        exec,
        "ClusterSets.py" = "cluster",
        "ParseHeaders.py" = "merge_headers",
        "AlignSets.py" = "align",
        "BuildConsensus.py" = "consensus",
        "AssignGenes.py" = "consensus_igblast",
        "MakeDb.py" = "make_db"
      )
      calls <<- c(calls, stage)
      if (identical(exec, "BuildConsensus.py")) observed$consensus_args <- args
      if (identical(exec, "AssignGenes.py")) observed$second_args <- args
      output_position <- match("-o", args)
      expect_false(is.na(output_position))
      writeLines(paste("artifact", stage), args[[output_position + 1L]])
      invisible(0L)
    },
    .package = "SIRAP"
  )

  result <- processSpatialSample(
    input,
    output_dir = file.path(work, "output with spaces"),
    igblast_db_dir = work,
    germline_glob = germline,
    igblast_base = igblast_root,
    igblast_exec = "/mock/igblastn",
    muscle_exec = "/mock/muscle",
    nproc = 1L
  )

  expect_identical(calls, c(
    "rad_reformat", "read_igblast", "cluster", "merge_headers", "align",
    "consensus", "consensus_igblast", "make_db"
  ))
  expect_s3_class(result, "sirap_spatial_process_result")
  expect_identical(result$stages$read_igblast$stage, "read")
  expect_equal(result$n_raw, 2L)
  expect_equal(result$n_filtered, 2L)
  expect_true(all(file.exists(unlist(result$artifacts[c(
    "clustered_fasta", "aligned_fasta", "consensus_fasta",
    "consensus_igblast", "final_airr"
  )], use.names = FALSE))))
  annotated <- data.table::fread(result$artifacts$annotated_read_airr)
  expect_identical(annotated$barcode, c("AACC", "TTGA"))
  expect_identical(annotated$v_call, c("IGHV1", "IGHV2"))
  expect_identical(annotated$d_call, c("IGHD1", "IGHD2"))
  expect_identical(annotated$j_call, c("IGHJ1", "IGHJ2"))
  expect_identical(readLines(result$artifacts$spatial_fasta), c(
    ">read-F|BARCODE=AACC|UMI=TGCA", "ACGTACGT",
    ">read-R|BARCODE=TTGA", "TTTTAAAA"
  ))
  expect_identical(
    observed$second_args[match("-s", observed$second_args) + 1L],
    result$artifacts$consensus_fasta
  )
  expect_identical(
    observed$second_args[match("--vdb", observed$second_args) + 1L],
    file.path(igblast_root, "database/imgt_human_ig_v")
  )
  copy_position <- match("--cf", observed$consensus_args)
  expect_identical(observed$consensus_args[copy_position + 1L], "BARCODE")
})

test_that("external stages require a freshly produced artifact", {
  work <- withr::local_tempdir(pattern = "sirap-stage-artifact-")
  output <- file.path(work, "stage.fasta")
  writeLines("previous valid artifact", output)
  command <- list(
    stage = "example",
    exec = "example-tool",
    args = c("-s", "input.fasta", "-o", output),
    output = output
  )

  testthat::local_mocked_bindings(
    .sirap_run_command = function(exec, args, log = NULL) invisible(0L),
    .package = "SIRAP"
  )
  expect_error(
    .sirap_run_stage_command(command),
    "produced no non-empty output",
    class = "sirap_igblast_command_error"
  )
  expect_identical(readLines(output), "previous valid artifact")
})
