test_that(".sirap_detect_input_type classifies inputs", {
  expect_equal(.sirap_detect_input_type(data.frame(x = 1)), "airr_df")
  expect_equal(.sirap_detect_input_type("x.tsv"), "airr_tsv")
  expect_equal(.sirap_detect_input_type("x.tsv.gz"), "airr_tsv")
  expect_equal(.sirap_detect_input_type("x.txt"), "airr_tsv")
  expect_equal(.sirap_detect_input_type("x.fa"), "raw_reads")
  expect_equal(.sirap_detect_input_type("x.fasta"), "raw_reads")
  expect_equal(.sirap_detect_input_type("x.fastq.gz"), "raw_reads")
  expect_equal(.sirap_detect_input_type("x.fq"), "raw_reads")
  expect_error(.sirap_detect_input_type("x.bam"), "Cannot detect")
  expect_error(.sirap_detect_input_type(c("a.fa", "b.fa")), "single file path")
})

test_that(".sirap_igblast_db_paths applies the naming convention and overrides", {
  dbs <- .sirap_igblast_db_paths("/db", receptor = "ig")
  expect_equal(dbs$v, file.path("/db", "IG_V"))
  expect_equal(dbs$d, file.path("/db", "IG_D"))
  expect_equal(dbs$j, file.path("/db", "IG_J"))
  expect_equal(dbs$c, file.path("/db", "IG_C"))

  dbs <- .sirap_igblast_db_paths("/db", receptor = "tr")
  expect_equal(dbs$v, file.path("/db", "TR_V"))
  expect_null(dbs$c)

  dbs <- .sirap_igblast_db_paths("/db", receptor = "ig", germline_db_v = "/other/V")
  expect_equal(dbs$v, "/other/V")
  expect_equal(dbs$d, file.path("/db", "IG_D"))

  expect_error(
    .sirap_igblast_db_paths(NULL, receptor = "ig"),
    "Cannot resolve germline database"
  )
})

test_that(".sirap_igblast_args builds the igblastn argument vector", {
  dbs <- .sirap_igblast_db_paths("/db", receptor = "ig")
  args <- .sirap_igblast_args(
    query_fasta = "in.fasta",
    output_tsv = "out.tsv",
    receptor = "ig",
    organism = "human",
    dbs = dbs,
    aux_file = "/aux/human_gl.aux",
    num_threads = 4
  )
  expect_equal(args[which(args == "-ig_seqtype") + 1L], "Ig")
  expect_equal(args[which(args == "-c_region_db") + 1L], file.path("/db", "IG_C"))
  expect_equal(args[which(args == "-auxiliary_data") + 1L], "/aux/human_gl.aux")
  expect_equal(args[which(args == "-outfmt") + 1L], "19")
  expect_equal(args[which(args == "-num_threads") + 1L], "4")
  expect_equal(args[which(args == "-out") + 1L], "out.tsv")

  dbs_tr <- .sirap_igblast_db_paths("/db", receptor = "tr")
  args_tr <- .sirap_igblast_args(
    query_fasta = "in.fasta",
    output_tsv = "out.tsv",
    receptor = "tr",
    dbs = dbs_tr,
    extra_args = c("-word_size", "9")
  )
  expect_equal(args_tr[which(args_tr == "-ig_seqtype") + 1L], "TCR")
  expect_false("-c_region_db" %in% args_tr)
  expect_false("-auxiliary_data" %in% args_tr)
  expect_equal(utils::tail(args_tr, 2L), c("-word_size", "9"))
})

test_that(".sirap_igblast_args passes and clamps the length thresholds", {
  dbs <- .sirap_igblast_db_paths("/db", receptor = "ig")
  args <- .sirap_igblast_args(
    query_fasta = "in.fasta",
    output_tsv = "out.tsv",
    receptor = "ig",
    dbs = dbs,
    min_v_length = 100,
    min_j_length = 10
  )
  expect_equal(args[which(args == "-min_V_length") + 1L], "100")
  expect_equal(args[which(args == "-min_J_length") + 1L], "10")

  # igblastn rejects -min_V_length below 9; small values are raised.
  clamped <- .sirap_igblast_args(
    query_fasta = "in.fasta",
    output_tsv = "out.tsv",
    receptor = "ig",
    dbs = dbs,
    min_v_length = 1
  )
  expect_equal(clamped[which(clamped == "-min_V_length") + 1L], "9")
  expect_false("-min_J_length" %in% clamped)

  # NULL omits the flags entirely (igblastn defaults apply).
  omitted <- .sirap_igblast_args(
    query_fasta = "in.fasta",
    output_tsv = "out.tsv",
    receptor = "ig",
    dbs = dbs
  )
  expect_false("-min_V_length" %in% omitted)
  expect_false("-min_J_length" %in% omitted)
})

test_that(".sirap_prepare_query_fasta streams fastq.gz to a plain fasta", {
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

  out <- .sirap_prepare_query_fasta(fq_gz, dir, "reads")
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
  expect_equal(.sirap_prepare_query_fasta(fa, dir, "x"), fa)
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

test_that(".sirap_run_command aborts on failure with the log tail", {
  skip_on_os("windows")
  dir <- withr::local_tempdir()
  log <- file.path(dir, "cmd.log")
  expect_error(
    .sirap_run_command("false", character(0), log = log),
    "Command failed with status"
  )
  expect_equal(.sirap_run_command("true", character(0)), 0L)
})

test_that("runIgblast runs a stub executable and returns the output path", {
  skip_on_os("windows")
  dir <- withr::local_tempdir()
  stub <- file.path(dir, "igblastn_stub.sh")
  writeLines(c(
    "#!/bin/sh",
    "out=\"\"",
    "while [ $# -gt 0 ]; do",
    "  if [ \"$1\" = \"-out\" ]; then out=\"$2\"; shift; fi",
    "  shift",
    "done",
    "printf 'sequence_id\\tlocus\\nread1\\tIGH\\n' > \"$out\""
  ), stub)
  Sys.chmod(stub, "0755")

  fa <- file.path(dir, "reads.fasta")
  writeLines(c(">read1", "ACGT"), fa)
  out_tsv <- file.path(dir, "out", "reads_igblast.tsv")

  result <- runIgblast(
    input = fa,
    output_tsv = out_tsv,
    receptor = "ig",
    igblast_db_dir = "/db",
    igblast_exec = stub
  )
  expect_equal(result, out_tsv)
  airr <- .sirap_read_airr(out_tsv)
  expect_equal(airr$sequence_id, "read1")
  expect_equal(airr$locus, "IGH")
})

test_that(".sirap_presto_commands builds the six pipeline commands", {
  paths <- .presto_paths("sample_spatial.fasta", 0.975, 0.51, "out")
  commands <- .sirap_presto_commands(
    paths = paths,
    receptor = "ig",
    organism = "human",
    identity = 0.975,
    consensus_freq = 0.51,
    germline_files = c("/g/imgt_human_IGHV.fasta", "/g/imgt_human_IGHJ.fasta"),
    igblast_base = "/software/igblast",
    igblast_exec = "/software/igblast/bin/igblastn",
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
    "ClusterSets.py", "ParseHeaders.py", "AlignSets.py",
    "BuildConsensus.py", "AssignGenes.py", "MakeDb.py"
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

  makedb_args <- commands[[6]]$args
  r_pos <- which(makedb_args == "-r")
  expect_equal(makedb_args[r_pos + 1L], "/g/imgt_human_IGHV.fasta")
  expect_equal(makedb_args[r_pos + 2L], "/g/imgt_human_IGHJ.fasta")
  expect_true("--extended" %in% makedb_args)
})
