# Helpers and exported functions for the raw-read IgBLAST step of the
# spatial processing pipeline. Raw reads come from the RRAD package as
# FASTA/FASTQ files, optionally gzipped, with headers like
# `id|umi:XXX|barcode:YYY`.

# Check that each external tool exists. `execs` is a named character vector;
# names are used in the error message. A tool passes if its expanded path is
# an existing file or if it resolves on the PATH.
.sirap_check_tools <- function(execs) {
  found <- vapply(execs, function(x) {
    file.exists(path.expand(x)) || nzchar(Sys.which(x))
  }, logical(1))
  if (any(!found)) {
    missing <- execs[!found]
    rlang::abort(paste0(
      "External tool", if (length(missing) == 1L) "" else "s", " not found: ",
      paste0(names(missing), " (", missing, ")", collapse = ", "),
      ". Install the tool or pass its full path."
    ))
  }
  invisible(TRUE)
}

# Run one external command with system2. Every argument is shell-quoted here,
# so callers pass plain strings. `log` redirects stdout and stderr to one
# file. `stdin`/`stdout` take file paths for redirection. Aborts on a
# non-zero exit status and includes the log tail when a log exists.
.sirap_run_command <- function(exec, args, log = NULL, stdin = "", stdout = "") {
  exec_path <- path.expand(exec)
  if (!file.exists(exec_path)) {
    exec_path <- exec
  }
  args <- as.character(args)
  shown <- paste(c(exec, args), collapse = " ")
  message("$ ", shown)
  if (!is.null(log)) {
    stdout <- log
    stderr <- log
  } else {
    stderr <- ""
  }
  status <- system2(exec_path, shQuote(args), stdout = stdout, stderr = stderr, stdin = stdin)
  if (!identical(status, 0L)) {
    detail <- ""
    if (!is.null(log) && file.exists(log)) {
      tail_lines <- utils::tail(readLines(log, warn = FALSE), 10L)
      detail <- paste0("\nLog tail (", log, "):\n", paste(tail_lines, collapse = "\n"))
    }
    rlang::abort(paste0("Command failed with status ", status, ": ", shown, detail))
  }
  invisible(status)
}

# Classify a pipeline input. Returns "airr_df" for a data frame, "airr_tsv"
# for a tabular file path, or "raw_reads" for a FASTA/FASTQ path (optionally
# gzipped). Aborts for anything else.
.sirap_detect_input_type <- function(input) {
  if (is.data.frame(input)) {
    return("airr_df")
  }
  if (!is.character(input) || length(input) != 1L || is.na(input)) {
    rlang::abort("`input` must be a data frame or a single file path.")
  }
  base <- sub("\\.gz$", "", input, ignore.case = TRUE)
  ext <- tolower(tools::file_ext(base))
  if (ext %in% c("tsv", "tab", "txt")) {
    return("airr_tsv")
  }
  if (ext %in% c("fa", "fasta", "fna", "fq", "fastq")) {
    return("raw_reads")
  }
  rlang::abort(paste0(
    "Cannot detect the input type from extension '.", ext, "'. ",
    "Accepted: .tsv/.tab/.txt for AIRR tables, ",
    ".fa/.fasta/.fna/.fq/.fastq for reads (each optionally .gz)."
  ))
}

# Read an AIRR TSV into a data frame. fread handles gzipped files natively.
.sirap_read_airr <- function(path) {
  if (!file.exists(path.expand(path))) {
    rlang::abort(paste0("AIRR file not found: ", path))
  }
  out <- data.table::fread(path.expand(path), sep = "\t", header = TRUE, data.table = FALSE)
  out
}

# Resolve IgBLAST germline database paths from a database directory and the
# receptor type, using the naming convention IG_V/IG_D/IG_J/IG_C and
# TR_V/TR_D/TR_J. Any non-NULL override wins. The C-region database is only
# used for `ig`.
.sirap_igblast_db_paths <- function(igblast_db_dir = NULL,
                                    receptor = c("ig", "tr"),
                                    germline_db_v = NULL,
                                    germline_db_d = NULL,
                                    germline_db_j = NULL,
                                    c_region_db = NULL) {
  receptor <- match.arg(receptor)
  prefix <- toupper(receptor)
  from_dir <- function(segment) {
    if (is.null(igblast_db_dir)) {
      return(NULL)
    }
    file.path(path.expand(igblast_db_dir), paste0(prefix, "_", segment))
  }
  dbs <- list(
    v = if (!is.null(germline_db_v)) path.expand(germline_db_v) else from_dir("V"),
    d = if (!is.null(germline_db_d)) path.expand(germline_db_d) else from_dir("D"),
    j = if (!is.null(germline_db_j)) path.expand(germline_db_j) else from_dir("J"),
    c = if (receptor == "ig") {
      if (!is.null(c_region_db)) path.expand(c_region_db) else from_dir("C")
    } else {
      NULL
    }
  )
  needed <- c("v", "d", "j")
  missing <- needed[vapply(dbs[needed], is.null, logical(1))]
  if (length(missing) > 0L) {
    rlang::abort(paste0(
      "Cannot resolve germline database", if (length(missing) == 1L) "" else "s",
      " for segment", if (length(missing) == 1L) "" else "s", " ",
      paste(toupper(missing), collapse = ", "),
      ". Provide `igblast_db_dir` or the per-segment paths."
    ))
  }
  dbs
}

# Build the igblastn argument vector for one raw-read annotation run.
# Pure function: it does not touch the filesystem.
.sirap_igblast_args <- function(query_fasta,
                                output_tsv,
                                receptor = c("ig", "tr"),
                                organism = "human",
                                dbs,
                                aux_file = NULL,
                                num_threads = 8,
                                extra_args = NULL) {
  receptor <- match.arg(receptor)
  seqtype <- if (receptor == "ig") "Ig" else "TCR"
  args <- c(
    "-query", path.expand(query_fasta),
    "-ig_seqtype", seqtype,
    "-organism", organism,
    "-germline_db_V", dbs$v,
    "-germline_db_D", dbs$d,
    "-germline_db_J", dbs$j
  )
  if (!is.null(dbs$c)) {
    args <- c(args, "-c_region_db", dbs$c)
  }
  if (!is.null(aux_file)) {
    args <- c(args, "-auxiliary_data", path.expand(aux_file))
  }
  args <- c(
    args,
    "-outfmt", "19",
    "-num_threads", as.character(num_threads),
    "-out", path.expand(output_tsv)
  )
  c(args, extra_args)
}

# Normalize raw reads to a plain FASTA that igblastn can read. igblastn
# cannot read gzip and only accepts FASTA, so this streams on disk with
# gunzip and awk; it never loads sequences into R. A plain FASTA passes
# through unchanged. FASTQ quality strings are dropped.
.sirap_prepare_query_fasta <- function(input, output_dir, sample_name) {
  input <- path.expand(input)
  if (!file.exists(input)) {
    rlang::abort(paste0("Input reads file not found: ", input))
  }
  is_gz <- grepl("\\.gz$", input, ignore.case = TRUE)
  base <- sub("\\.gz$", "", input, ignore.case = TRUE)
  ext <- tolower(tools::file_ext(base))
  is_fastq <- ext %in% c("fq", "fastq")
  if (!is_gz && !is_fastq) {
    return(input)
  }
  out_fasta <- file.path(output_dir, paste0(sample_name, "_input.fasta"))
  if (!is_fastq) {
    # Gzipped FASTA: stream-decompress only.
    .sirap_run_command("gunzip", c("-c", input), stdout = out_fasta)
    return(out_fasta)
  }
  fastq <- input
  tmp_fastq <- NULL
  if (is_gz) {
    tmp_fastq <- file.path(output_dir, paste0(sample_name, "_input.fastq"))
    .sirap_run_command("gunzip", c("-c", input), stdout = tmp_fastq)
    fastq <- tmp_fastq
  }
  # FASTQ to FASTA in one streaming pass; qualities are dropped.
  .sirap_run_command(
    "awk",
    "NR % 4 == 1 { sub(/^@/, \">\"); print } NR % 4 == 2 { print }",
    stdin = fastq,
    stdout = out_fasta
  )
  if (!is.null(tmp_fastq)) {
    unlink(tmp_fastq)
  }
  out_fasta
}

#' Annotate raw reads with IgBLAST
#'
#' Runs `igblastn` on raw reads and writes an AIRR rearrangement table
#' (`-outfmt 19`). The input can be FASTA or FASTQ, optionally gzipped.
#' Gzipped and FASTQ inputs are first normalized to a plain FASTA next to
#' `output_tsv` (FASTQ quality strings are dropped, because igblastn only
#' accepts FASTA). One receptor locus group runs per call: `ig` uses the
#' `IG_V`/`IG_D`/`IG_J`/`IG_C` databases, `tr` uses `TR_V`/`TR_D`/`TR_J`.
#'
#' Database paths come from `igblast_db_dir` plus the naming convention
#' above. Per-segment arguments override the convention. igblastn output is
#' logged to `<output_tsv stem>.log`.
#'
#' @param input Path to a FASTA/FASTQ file, optionally gzipped.
#' @param output_tsv Path for the AIRR TSV that igblastn writes.
#' @param receptor Either `ig` for BCR or `tr` for TCR.
#' @param organism Organism passed to `igblastn -organism`.
#' @param igblast_db_dir Directory that holds the germline databases named
#'   `IG_V`, `IG_D`, `IG_J`, `IG_C`, `TR_V`, `TR_D`, `TR_J`.
#' @param igblast_exec Path to the `igblastn` executable.
#' @param aux_file Optional path passed to `-auxiliary_data`.
#' @param germline_db_v,germline_db_d,germline_db_j Optional per-segment
#'   database paths. Each overrides the `igblast_db_dir` convention.
#' @param c_region_db Optional C-region database path (`ig` only).
#' @param num_threads Threads passed to `-num_threads`.
#' @param extra_args Optional character vector of extra igblastn arguments.
#' @return The path to `output_tsv`, invisibly.
#' @export
runIgblast <- function(input,
                       output_tsv,
                       receptor = c("ig", "tr"),
                       organism = "human",
                       igblast_db_dir = NULL,
                       igblast_exec = "~/software/igblast/bin/igblastn",
                       aux_file = NULL,
                       germline_db_v = NULL,
                       germline_db_d = NULL,
                       germline_db_j = NULL,
                       c_region_db = NULL,
                       num_threads = 8,
                       extra_args = NULL) {
  receptor <- match.arg(receptor)
  .sirap_check_tools(c(igblastn = igblast_exec))
  dbs <- .sirap_igblast_db_paths(
    igblast_db_dir = igblast_db_dir,
    receptor = receptor,
    germline_db_v = germline_db_v,
    germline_db_d = germline_db_d,
    germline_db_j = germline_db_j,
    c_region_db = c_region_db
  )
  output_tsv <- path.expand(output_tsv)
  output_dir <- dirname(output_tsv)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  sample_name <- tools::file_path_sans_ext(basename(output_tsv))
  query_fasta <- .sirap_prepare_query_fasta(input, output_dir, sample_name)
  args <- .sirap_igblast_args(
    query_fasta = query_fasta,
    output_tsv = output_tsv,
    receptor = receptor,
    organism = organism,
    dbs = dbs,
    aux_file = aux_file,
    num_threads = num_threads,
    extra_args = extra_args
  )
  log <- paste0(tools::file_path_sans_ext(output_tsv), ".log")
  .sirap_run_command(igblast_exec, args, log = log)
  invisible(output_tsv)
}

#' Parse RRAD tags out of AIRR sequence IDs
#'
#' RRAD writes read headers like `id|umi:XXX|barcode:YYY`. This function
#' splits that format: it adds `umi` and `barcode` columns and, when
#' `strip = TRUE`, replaces the ID column with the bare read ID. Rows
#' without a tag get `NA` in the matching column. For other header formats,
#' see [extractSpatialIds()].
#'
#' @param airr_df AIRR-style data frame.
#' @param id Column that contains the read identifiers.
#' @param strip Replace the ID column with the part before the first `|`.
#' @return The data frame with `umi` and `barcode` columns added.
#' @export
parseRadSequenceIds <- function(airr_df, id = "sequence_id", strip = TRUE) {
  .check_required_cols(airr_df, id, "airr_df")
  out <- .as_data_frame(airr_df)
  ids <- as.character(out[[id]])
  out$umi <- stringr::str_match(ids, "\\|umi:([^|]*)")[, 2]
  out$barcode <- stringr::str_match(ids, "\\|barcode:([^|]*)")[, 2]
  if (strip) {
    out[[id]] <- sub("\\|.*$", "", ids)
  }
  out
}
