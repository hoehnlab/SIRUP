# Adds edlib HW edit distance per row for barcode correction analysis.
# Expects a data.table with columns for final/original barcodes and optional strand.

#' Add edlib HW edit distance per row
#' @param DT data.table containing barcode columns
#' @param query_col column name for corrected/final barcode
#' @param target_col column name for original barcode
#' @param strand_col column name for strand (optional)
#' @param corrected_col if present, RC only when corrected==TRUE
#' @param incorrect_col unused (kept for compatibility)
#' @param out_col output column for edit distance
#' @param pad string to pad target before alignment
#' @param add_eq additional equalities for edlibR
#' @return DT with new column
#' @export

dt_add_hw_edist <- function(
    DT,
    query_col     = "final_barcode",
    target_col    = "original_barcode",
    strand_col    = "strand",      # if missing in DT, RC is skipped
    corrected_col = "corrected",   # optional: if present, RC only where corrected==TRUE
    incorrect_col = "incorrect",
    out_col       = "ed_hw",
    pad           = "NN",
    add_eq        = list(c('N','A'), c('N','C'), c('N','G'), c('N','T'))
) {
  stopifnot(is.data.table(DT))
  if (!requireNamespace("edlibR", quietly = TRUE)) {
    stop("edlibR is required for dt_add_hw_edist()")
  }
  req <- c(query_col, target_col)
  miss <- setdiff(req, names(DT))
  if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "))

  # prepare target (optional RC), then pad once
  tgt <- DT[[target_col]]

  has_strand    <- !is.null(strand_col)    && strand_col    %in% names(DT)
  has_corrected <- !is.null(corrected_col) && corrected_col %in% names(DT)

  if (has_strand) {
    s <- DT[[strand_col]]
    if (has_corrected) {
      cor <- DT[[corrected_col]]
      neg <- which(!is.na(s) & s == "-" & !is.na(cor) & cor)
    } else {
      neg <- which(!is.na(s) & s == "-")
    }
    if (length(neg)) tgt[neg] <- revcomp(tgt[neg])
  }

  if (!identical(pad, "") && !is.na(pad)) tgt <- paste0(pad, tgt, pad)

  # edlib distance per row
  qv <- DT[[query_col]]
  dists <- vapply(seq_len(nrow(DT)), function(i) {
    qi <- qv[i]; ti <- tgt[i]
    if (is.na(qi) || is.na(ti)) return(NA_integer_)
    edlibR::align(
      query = qi,
      target = ti,
      mode = "HW",
      task = "distance",
      additionalEqualities = add_eq
    )$editDistance
  }, integer(1L))

  DT[, (out_col) := dists][]
}
