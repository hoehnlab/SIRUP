.presto_paths <- function(input_fasta, identity, consensus_freq, output_dir = NULL) {
  if (is.null(output_dir)) {
    output_dir <- dirname(input_fasta)
  }
  stem <- tools::file_path_sans_ext(basename(input_fasta))
  cluster <- file.path(
    output_dir,
    sprintf("%s_ident_%s_barcode_cluster_cluster-pass.fasta", stem, identity)
  )
  reheader <- sub("\\.fasta$", "_reheader.fasta", cluster)
  aligned <- sub("\\.fasta$", "_align-pass.fasta", reheader)
  consensus <- file.path(
    output_dir,
    sprintf("%s_ident_%s_barcode_cluster_freq_%s_cluster-pass_reheader_align-pass_consensus-pass.fasta", stem, identity, consensus_freq)
  )
  fmt7 <- sub("\\.fasta$", "_igblast.fmt7", consensus)
  db_pass <- sub("_igblast\\.fmt7$", "_igblast_db-pass.tsv", fmt7)
  list(
    input = input_fasta,
    cluster = cluster,
    reheader = reheader,
    aligned = aligned,
    consensus = consensus,
    fmt7 = fmt7,
    db_pass = db_pass
  )
}

#' Run pRESTO/IgBLAST consensus calling for one spatial FASTA
#'
#' This keeps users inside R while running the external Immcantation/pRESTO
#' command-line tools. Output from each command is printed directly to the R
#' console via `system()`.
#'
#' @param input_fasta Spatial FASTA written from preprocessed receptor reads.
#' @param output_dir Directory for pRESTO/IgBLAST outputs. Defaults to the
#'   input FASTA directory.
#' @param receptor Either `ig` for BCR or `tr` for TCR IgBLAST loci.
#' @param identity Cluster identity threshold.
#' @param consensus_freq Consensus frequency threshold.
#' @param germline_glob Germline FASTA glob for `MakeDb.py`.
#' @param igblast_base Path passed to `AssignGenes.py igblast -b`.
#' @param igblast_exec Optional path passed to `AssignGenes.py --exec`. When
#'   `NULL`, defaults to `file.path(igblast_base, "bin", "igblastn")`.
#' @param muscle_exec Path passed to `AlignSets.py --exec`.
#' @param nproc Number of processors for pRESTO and AssignGenes steps.
#' @param cluster_field FASTA header field used for within-spot clustering.
#' @param delimiters String of pRESTO delimiter values passed after `--delim`.
#'   Defaults to `"'|' '=' ','"`.
#' @param parse_fields Header fields passed to `ParseHeaders.py merge -f`.
#' @param parse_key Header field name created by `ParseHeaders.py merge -k`.
#' @param parse_act Action passed to `ParseHeaders.py merge --act`.
#' @param consensus_min_count Minimum sequence count passed to
#'   `BuildConsensus.py -n`.
#' @param consensus_copy_fields Optional fields passed to `BuildConsensus.py
#'   --cf`.
#' @param consensus_copy_action Action passed to `BuildConsensus.py --act` when
#'   `consensus_copy_fields` is provided.
#' @param align_log Optional log path passed to `AlignSets.py --log`.
#' @param consensus_log Optional log path passed to `BuildConsensus.py --log`.
#' @return A named list of output paths, invisibly.
#' @export
processSpatialSample <- function(input_fasta,
                                      output_dir = NULL,
                                      receptor = c("ig", "tr"),
                                      identity = 0.975,
                                      consensus_freq = 0.51,
                                      germline_glob,
                                      igblast_base = "~/software/igblast",
                                      igblast_exec = NULL,
                                      muscle_exec = "~/software/muscle3.8.31_i86linux64",
                                      nproc = 8,
                                      cluster_field = "barcode2um",
                                      delimiters = "'|' '=' ','",
                                      parse_fields = c(cluster_field, "cluster"),
                                      parse_key = paste0(cluster_field, "_cluster"),
                                      parse_act = "set",
                                      consensus_min_count = 1,
                                      consensus_copy_fields = NULL,
                                      consensus_copy_action = "majority",
                                      align_log = NULL,
                                      consensus_log = NULL) {
  receptor <- match.arg(receptor)
  paths <- .presto_paths(input_fasta, identity, consensus_freq, output_dir)
  if (!dir.exists(dirname(paths$cluster))) {
    dir.create(dirname(paths$cluster), recursive = TRUE, showWarnings = FALSE)
  }

  quote_path <- function(x) shQuote(path.expand(x))
  if (is.null(igblast_exec)) {
    igblast_exec <- file.path(igblast_base, "bin", "igblastn")
  }
  delim_args <- paste("--delim", delimiters)
  align_log_args <- if (is.null(align_log)) "" else paste("--log", quote_path(align_log))
  consensus_log_args <- if (is.null(consensus_log)) "" else paste("--log", quote_path(consensus_log))
  consensus_copy_args <- ""
  if (!is.null(consensus_copy_fields) && length(consensus_copy_fields) > 0L) {
    consensus_copy_args <- sprintf(
      "--cf %s --act %s",
      paste(consensus_copy_fields, collapse = " "),
      consensus_copy_action
    )
  }

  commands <- c(
    sprintf(
      "ClusterSets.py set -s %s -f %s -k cluster --ident %s --cluster vsearch %s --nproc %s -o %s",
      quote_path(paths$input),
      cluster_field,
      identity,
      delim_args,
      nproc,
      quote_path(paths$cluster)
    ),
    sprintf(
      "ParseHeaders.py merge -s %s -f %s -k %s %s --act %s",
      quote_path(paths$cluster),
      paste(parse_fields, collapse = " "),
      parse_key,
      delim_args,
      parse_act
    ),
    sprintf(
      "AlignSets.py muscle -s %s %s %s --nproc %s --bf %s --exec %s",
      quote_path(paths$reheader),
      align_log_args,
      delim_args,
      nproc,
      parse_key,
      quote_path(muscle_exec)
    )
  )

  commands <- c(
    commands,
    sprintf(
      "BuildConsensus.py -s %s -o %s %s --nproc %s -n %s -q 20 --freq %s --bf %s --maxerror 0.1 --maxgap 0.3 --fasta %s %s",
      quote_path(paths$aligned),
      quote_path(paths$consensus),
      consensus_log_args,
      nproc,
      consensus_min_count,
      consensus_freq,
      parse_key,
      delim_args,
      consensus_copy_args
    ),
    sprintf(
      "AssignGenes.py igblast -s %s -b %s --organism human --loci %s --format blast --exec %s --nproc %s",
      quote_path(paths$consensus),
      quote_path(igblast_base),
      receptor,
      quote_path(igblast_exec),
      nproc
    ),
    sprintf(
      "MakeDb.py igblast -i %s -s %s -r %s --extended",
      quote_path(paths$fmt7),
      quote_path(paths$consensus),
      path.expand(germline_glob)
    )
  )

  for (command in commands) {
    message("$ ", command)
    status <- system(command)
    if (!identical(status, 0L)) {
      rlang::abort(paste0("Command failed with status ", status, ": ", command))
    }
  }

  invisible(paths)
}

#' Run pRESTO/IgBLAST consensus calling for multiple spatial FASTA files
#'
#' @param input_fastas Character vector of spatial FASTA files.
#' @param output_dirs Optional vector of output directories, one per FASTA.
#' @param ... Passed to `processSpatialSample()`.
#' @return A list of output path lists, invisibly.
#' @export
processSpatialSamples <- function(input_fastas,
                                  output_dirs = NULL,
                                  ...) {
  if (is.null(output_dirs)) {
    output_dirs <- rep(NA_character_, length(input_fastas))
  }
  if (length(output_dirs) != length(input_fastas)) {
    rlang::abort("output_dirs must be NULL or the same length as input_fastas.")
  }
  results <- vector("list", length(input_fastas))
  names(results) <- basename(input_fastas)
  for (i in seq_along(input_fastas)) {
    output_dir <- output_dirs[[i]]
    if (is.na(output_dir)) {
      output_dir <- NULL
    }
    message("Processing spatial FASTA ", i, " of ", length(input_fastas), ": ", input_fastas[[i]])
    results[[i]] <- processSpatialSample(
      input_fasta = input_fastas[[i]],
      output_dir = output_dir,
      ...
    )
  }
  invisible(results)
}
