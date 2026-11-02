test_that(".sirup_detect_input_type classifies inputs", {
  expect_equal(.sirup_detect_input_type(data.frame(x = 1)), "airr_df")
  expect_equal(.sirup_detect_input_type("x.tsv"), "airr_tsv")
  expect_equal(.sirup_detect_input_type("x.tsv.gz"), "airr_tsv")
  expect_equal(.sirup_detect_input_type("x.txt"), "airr_tsv")
  expect_equal(.sirup_detect_input_type("x.fa"), "raw_reads")
  expect_equal(.sirup_detect_input_type("x.fasta"), "raw_reads")
  expect_equal(.sirup_detect_input_type("x.fastq.gz"), "raw_reads")
  expect_equal(.sirup_detect_input_type("x.fq"), "raw_reads")
  expect_error(.sirup_detect_input_type("x.bam"), "Cannot detect")
  expect_error(.sirup_detect_input_type(c("a.fa", "b.fa")), "single file path")
})

test_that(".sirup_prepare_query_fasta streams fastq.gz to a plain fasta", {
  skip_on_os("windows")
  dir <- withr::local_tempdir()
  fq <- file.path(dir, "reads.fastq")
  writeLines(c(
    "@read1|umi:AACG|barcode:BC01",
    "ACGTACGT",
    "+",
    "IIIIIIII",
    "@read2",
    "TTTTAAAA",
    "+",
    "IIIIIIII"
  ), fq)
  fq_gz <- file.path(dir, "reads.fastq.gz")
  con <- gzfile(fq_gz, "w")
  writeLines(readLines(fq), con)
  close(con)

  out <- .sirup_prepare_query_fasta(fq_gz, dir, "reads")
  expect_true(file.exists(out))
  lines <- readLines(out)
  expect_equal(lines, c(
    ">read1|umi:AACG|barcode:BC01",
    "ACGTACGT",
    ">read2",
    "TTTTAAAA"
  ))
  expect_false(file.exists(file.path(dir, "reads_input.fastq")))

  # A plain fasta passes through with no copy.
  fa <- file.path(dir, "reads.fasta")
  writeLines(c(">read1", "ACGT"), fa)
  expect_equal(.sirup_prepare_query_fasta(fa, dir, "x"), fa)
})

test_that("parseRadSequenceIds splits RRAD tags", {
  airr <- data.frame(
    sequence_id = c(
      "read1|umi:AACG|barcode:BC01",
      "read2|barcode:BC02",
      "read3"
    ),
    stringsAsFactors = FALSE
  )
  out <- parseRadSequenceIds(airr)
  expect_equal(out$sequence_id, c("read1", "read2", "read3"))
  expect_equal(out$umi, c("AACG", NA, NA))
  expect_equal(out$barcode, c("BC01", "BC02", NA))

  out_keep <- parseRadSequenceIds(airr, strip = FALSE)
  expect_equal(out_keep$sequence_id, airr$sequence_id)
})

test_that(".sirup_run_command aborts on failure with the log tail", {
  skip_on_os("windows")
  dir <- withr::local_tempdir()
  log <- file.path(dir, "cmd.log")
  expect_error(
    .sirup_run_command("false", character(0), log = log),
    "Command failed with status"
  )
  expect_equal(.sirup_run_command("true", character(0)), 0L)
})

test_that("runIgblast passes the plain FASTA and locus to airrblast::assignGenes", {
  dir <- withr::local_tempdir()
  fa <- file.path(dir, "reads.fasta")
  writeLines(c(">read1", "ACGT"), fa)
  out_tsv <- file.path(dir, "out", "reads_igblast.tsv")
  seen <- list()
  testthat::local_mocked_bindings(
    assignGenes = function(file, igblast, refs, organism, outfile, nproc, locus, ...) {
      seen <<- list(file = file, igblast = igblast, refs = refs, organism = organism,
                    outfile = outfile, nproc = nproc, locus = locus)
      writeLines(c("sequence_id\tlocus", "read1\tIGH"), outfile)
      invisible(NULL)
    },
    .package = "airrblast"
  )
  result <- runIgblast(
    input = fa, output_tsv = out_tsv, receptor = "tr",
    igblast_base = "/ig", num_threads = 3
  )
  expect_equal(result, out_tsv)
  expect_equal(seen$file, fa)
  expect_equal(seen$igblast, "/ig")
  expect_equal(seen$refs, "/ig")
  expect_equal(seen$locus, "TR")
  expect_equal(seen$nproc, 3)
  expect_equal(.sirup_read_airr(out_tsv)$sequence_id, "read1")
})

test_that(".sirup_presto_commands builds the four pRESTO commands", {
  paths <- .presto_paths("sample_spatial.fasta", 0.975, 0.51, "out")
  commands <- .sirup_presto_commands(
    paths = paths,
    receptor = "ig",
    organism = "human",
    identity = 0.975,
    consensus_freq = 0.51,
    muscle_exec = "/software/muscle",
    nproc = 2,
    cluster_field = "barcode",
    delimiters = c("|", "=", ","),
    parse_fields = c("barcode", "cluster"),
    parse_key = "barcode_cluster",
    parse_act = "cat",
    consensus_min_count = 1,
    consensus_copy_fields = NULL,
    consensus_copy_action = "majority",
    align_log = NULL,
    consensus_log = NULL
  )
  execs <- vapply(commands, function(cmd) cmd$exec, character(1))
  expect_equal(execs, c(
    "ClusterSets.py", "ParseHeaders.py", "AlignSets.py", "BuildConsensus.py"
  ))

  cluster_args <- commands[[1]]$args
  expect_equal(cluster_args[1], "set")
  expect_equal(cluster_args[which(cluster_args == "-f") + 1L], "barcode")
  expect_equal(cluster_args[which(cluster_args == "--ident") + 1L], "0.975")
  expect_true(all(c("|", "=", ",") %in% cluster_args))

  align_args <- commands[[3]]$args
  expect_false("--log" %in% align_args)
  expect_equal(align_args[which(align_args == "--exec") + 1L], "/software/muscle")

  consensus_args <- commands[[4]]$args
  expect_equal(consensus_args[which(consensus_args == "--freq") + 1L], "0.51")
  expect_false("--cf" %in% consensus_args)
})
