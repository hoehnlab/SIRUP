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
#' @param igblast_exec Path passed to `AssignGenes.py --exec`.
#' @param muscle_exec Path passed to `AlignSets.py --exec`.
#' @param nproc Number of processors for pRESTO steps.
#' @param cluster_field FASTA header field used for within-spot clustering.
#' @return A named list of output paths, invisibly.
#' @export
processSpatialSample <- function(input_fasta,
                                      output_dir = NULL,
                                      receptor = c("ig", "tr"),
                                      identity = 0.975,
                                      consensus_freq = 0.51,
                                      germline_glob,
                                      igblast_base = "~/software/igblast",
                                      igblast_exec = "~/software/igblast/bin/igblastn",
                                      muscle_exec = "~/software/muscle3.8.31_i86linux64",
                                      nproc = 8,
                                      cluster_field = "barcode2um") {
  receptor <- match.arg(receptor)
  paths <- .presto_paths(input_fasta, identity, consensus_freq, output_dir)
  if (!dir.exists(dirname(paths$cluster))) {
    dir.create(dirname(paths$cluster), recursive = TRUE, showWarnings = FALSE)
  }

  quote_path <- function(x) shQuote(path.expand(x))
  delim_args <- paste("--delim", shQuote("|"), shQuote("="), shQuote(","))
  half_nproc <- max(1, floor(nproc / 2))

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
      "ParseHeaders.py merge -s %s -f %s cluster -k %s %s --act set",
      quote_path(paths$cluster),
      cluster_field,
      paste0(cluster_field, "_cluster"),
      delim_args
    ),
    sprintf(
      "AlignSets.py muscle -s %s %s --nproc %s --bf %s --exec %s",
      quote_path(paths$reheader),
      delim_args,
      half_nproc,
      paste0(cluster_field, "_cluster"),
      quote_path(muscle_exec)
    ),
    sprintf(
      "BuildConsensus.py -s %s -o %s --nproc %s -n 1 -q 20 --freq %s --bf %s --maxerror 0.1 --maxgap 0.3 --fasta",
      quote_path(paths$aligned),
      quote_path(paths$consensus),
      nproc,
      consensus_freq,
      paste0(cluster_field, "_cluster")
    ),
    sprintf(
      "AssignGenes.py igblast -s %s -b %s --organism human --loci %s --format blast --exec %s",
      quote_path(paths$consensus),
      quote_path(igblast_base),
      receptor,
      quote_path(igblast_exec)
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
