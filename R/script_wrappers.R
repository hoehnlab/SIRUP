#' Convert PAF to CSV
#'
#' Writes two outputs: `output.csv` and `output_processed.csv`.
#'
#' @param input_paf path to `.paf` or `.paf.gz`
#' @param output_dir output directory
#' @param include_rest whether to include a `rest` column (full original PAF line) in `output.csv`
#' @param nthreads `data.table::fread` thread count
#' @return Invisible list with output paths and row count
#' @export
paf2csv <- function(input_paf,
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

  metadata <- parseRadSequenceIds(
    data.frame(
      sequence_id = as.character(out$qname),
      stringsAsFactors = FALSE
    ),
    require_barcode = FALSE
  )
  processed <- data.table::data.table(
    id = metadata$rad_read_id,
    umi = metadata$umi,
    barcode = metadata$barcode
  )
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
