# Resolve shell-script paths in both source tree and installed package layouts.
.rad_resolve_script_path <- function(script, script_dir = NULL) {
  candidates <- character()

  if (!is.null(script_dir)) {
    candidates <- c(candidates, file.path(script_dir, script))
  }

  pkg_script <- system.file("scripts", script, package = "rad")
  if (nzchar(pkg_script)) {
    candidates <- c(candidates, pkg_script)
  }

  wd <- getwd()
  candidates <- c(
    candidates,
    file.path(wd, "scripts", script),
    file.path(wd, "inst", "scripts", script)
  )

  existing <- unique(candidates[file.exists(candidates)])
  if (length(existing) == 0) {
    stop(
      "Could not resolve script '",
      script,
      "'. Set `script_dir` explicitly or run from the package root."
    )
  }

  normalizePath(existing[[1]], mustWork = TRUE)
}

#' Run a shell script from R
#'
#' @param script script filename (for example, `"extract_rad_headers.sh"`)
#' @param args character vector of CLI arguments passed to the script
#' @param script_dir optional directory containing scripts
#' @param env optional environment variables (`c("KEY=VALUE")`)
#' @param intern whether to capture combined stdout/stderr and return it
#' @return If `intern = FALSE`, an invisible list with script path and status. If `intern = TRUE`, captured output lines.
#' @export
run_shell_script <- function(script,
                             args = character(),
                             script_dir = NULL,
                             env = character(),
                             intern = FALSE) {
  script_path <- .rad_resolve_script_path(script = script, script_dir = script_dir)

  output <- system2(
    command = "bash",
    args = c(script_path, args),
    env = env,
    stdout = if (intern) TRUE else "",
    stderr = if (intern) TRUE else ""
  )

  if (intern) {
    status <- attr(output, "status")
    if (is.null(status)) {
      status <- 0L
    }
    if (!identical(as.integer(status), 0L)) {
      stop("Script failed with status ", status, ": ", script_path)
    }
    return(output)
  }

  if (!identical(as.integer(output), 0L)) {
    stop("Script failed with status ", output, ": ", script_path)
  }

  invisible(list(script = script_path, status = 0L))
}

#' Convert PAF to CSV natively in R
#'
#' This is a native-R replacement for `scripts/paf2csv.sh` and writes the same two
#' outputs: `output.csv` and `output_processed.csv`.
#'
#' @param input_paf path to `.paf` or `.paf.gz`
#' @param output_dir output directory
#' @param include_rest whether to include a `rest` column (full original PAF line) in `output.csv`
#' @param nthreads `data.table::fread` thread count
#' @return Invisible list with output paths and row count
#' @export
paf2csv_r <- function(input_paf,
                      output_dir,
                      include_rest = TRUE,
                      nthreads = 1L) {
  if (!file.exists(input_paf)) {
    stop("Input PAF does not exist: ", input_paf)
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  paf <- data.table::fread(
    input = input_paf,
    sep = "\t",
    header = FALSE,
    fill = TRUE,
    data.table = TRUE,
    quote = "",
    nThread = nthreads
  )

  if (ncol(paf) < 12) {
    stop("PAF must have at least 12 tab-separated fields.")
  }

  out <- paf[, 1:12]
  data.table::setnames(
    out,
    c(
      "qname", "qlen", "qstart", "qend", "strand", "tname",
      "tlen", "tstart", "tend", "nmatch", "alen", "mapq"
    )
  )

  if (isTRUE(include_rest)) {
    con <- if (grepl("\\.gz$", input_paf, ignore.case = TRUE)) {
      gzfile(input_paf, open = "rt")
    } else {
      file(input_paf, open = "rt")
    }
    on.exit(close(con), add = TRUE)

    rest_lines <- readLines(con = con, warn = FALSE)
    if (length(rest_lines) == nrow(out)) {
      out[, rest := rest_lines]
    } else {
      warning("Could not align raw PAF lines to parsed rows; `rest` column omitted.")
    }
  }

  qname <- as.character(out$qname)
  id <- character(length(qname))
  umi <- character(length(qname))
  barcode <- character(length(qname))

  for (i in seq_along(qname)) {
    token <- qname[[i]]
    if (is.na(token) || !nzchar(token)) {
      next
    }

    parts <- strsplit(token, split = "|", fixed = TRUE)[[1]]
    if (length(parts) == 0) {
      next
    }

    id[[i]] <- parts[[1]]

    umi_hit <- parts[startsWith(parts, "umi:")]
    bc_hit <- parts[startsWith(parts, "barcode:")]

    if (length(umi_hit) > 0) {
      umi[[i]] <- sub("^umi:", "", umi_hit[[1]])
    }
    if (length(bc_hit) > 0) {
      barcode[[i]] <- sub("^barcode:", "", bc_hit[[1]])
    }
  }

  processed <- data.table::data.table(id = id, umi = umi, barcode = barcode)
  processed <- cbind(processed, out)

  csv_output <- file.path(output_dir, "output.csv")
  processed_output <- file.path(output_dir, "output_processed.csv")

  data.table::fwrite(out, file = csv_output)
  data.table::fwrite(processed, file = processed_output)

  invisible(list(
    output_csv = normalizePath(csv_output, mustWork = TRUE),
    output_processed_csv = normalizePath(processed_output, mustWork = TRUE),
    n_rows = nrow(out)
  ))
}

#' Convert PAF to CSV via native R or shell script
#'
#' @param input_paf path to `.paf` or `.paf.gz`
#' @param output_dir output directory
#' @param engine either `"R"` (native parser) or `"shell"` (runs `paf2csv.sh`)
#' @param script_dir optional script directory when `engine = "shell"`
#' @param include_rest passed to `paf2csv_r()` when `engine = "R"`
#' @param nthreads passed to `paf2csv_r()` when `engine = "R"`
#' @return Invisible list of output paths
#' @export
paf2csv <- function(input_paf,
                    output_dir,
                    engine = c("R", "shell"),
                    script_dir = NULL,
                    include_rest = TRUE,
                    nthreads = 1L) {
  engine <- match.arg(engine)

  if (engine == "R") {
    return(paf2csv_r(
      input_paf = input_paf,
      output_dir = output_dir,
      include_rest = include_rest,
      nthreads = nthreads
    ))
  }

  run_shell_script(
    script = "paf2csv.sh",
    args = c(input_paf, output_dir),
    script_dir = script_dir
  )

  invisible(list(
    output_csv = normalizePath(file.path(output_dir, "output.csv"), mustWork = FALSE),
    output_processed_csv = normalizePath(file.path(output_dir, "output_processed.csv"), mustWork = FALSE)
  ))
}

#' Wrapper for `extract_rad_headers.sh`
#'
#' @param input_fastq FASTQ/FASTQ.gz path
#' @param output_csv output CSV path
#' @param include_length include read-length column (`--length`)
#' @param script_dir optional script directory
#' @return Invisible list with script path and status
#' @export
extract_rad_headers <- function(input_fastq,
                                output_csv,
                                include_length = FALSE,
                                script_dir = NULL) {
  args <- c(input_fastq, output_csv)
  if (isTRUE(include_length)) {
    args <- c(args, "--length")
  }
  run_shell_script("extract_rad_headers.sh", args = args, script_dir = script_dir)
}

#' Wrapper for `extract_blaze_headers.sh`
#'
#' @param input_fastq FASTQ/FASTQ.gz path
#' @param output_csv output CSV path
#' @param script_dir optional script directory
#' @return Invisible list with script path and status
#' @export
extract_blaze_headers <- function(input_fastq,
                                  output_csv,
                                  script_dir = NULL) {
  run_shell_script(
    script = "extract_blaze_headers.sh",
    args = c(input_fastq, output_csv),
    script_dir = script_dir
  )
}

#' Wrapper for `extract_flexiplex_headers.sh`
#'
#' @param input_fastq FASTQ/FASTQ.gz path
#' @param output_csv output CSV path
#' @param script_dir optional script directory
#' @return Invisible list with script path and status
#' @export
extract_flexiplex_headers <- function(input_fastq,
                                      output_csv,
                                      script_dir = NULL) {
  run_shell_script(
    script = "extract_flexiplex_headers.sh",
    args = c(input_fastq, output_csv),
    script_dir = script_dir
  )
}

#' Wrapper for `extract_bam_headers.sh`
#'
#' @param input_bam BAM path
#' @param output_csv output CSV path
#' @param script_dir optional script directory
#' @return Invisible list with script path and status
#' @export
extract_bam_headers <- function(input_bam,
                                output_csv,
                                script_dir = NULL) {
  run_shell_script(
    script = "extract_bam_headers.sh",
    args = c(input_bam, output_csv),
    script_dir = script_dir
  )
}
