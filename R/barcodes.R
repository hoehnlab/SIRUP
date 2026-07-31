#' Read a Visium HD barcode mapping file
#'
#' Reads a barcode mapping from CSV/TSV or parquet. Parquet support requires
#' the optional `arrow` package.
#'
#' @param path Path to a barcode mapping file.
#' @return A data frame.
#' @export
readBarcodeMapping <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "parquet") {
    .require_namespace("arrow")
    return(.as_data_frame(arrow::read_parquet(path)))
  }
  if (ext %in% c("tsv", "tab")) {
    return(.as_data_frame(data.table::fread(path, sep = "\t")))
  }
  .as_data_frame(data.table::fread(path))
}

#' Extract spatial barcodes from sequence IDs
#'
#' @param airr_df AIRR-style data frame.
#' @param sequence_id_col Column containing read or consensus identifiers.
#' @param barcode_col Name of the output barcode column.
#' @param regex Regular expression used to extract the barcode. If the regex
#'   contains a capture group, the first capture group is used. Otherwise the
#'   full match is used.
#' @param prefix,suffix Optional text added before or after the extracted match.
#' @return Input data with `barcode_col` added.
#' @export
extractSpatialBarcodes <- function(airr_df,
                                    sequence_id_col = "sequence_id",
                                    barcode_col = "barcode",
                                    regex,
                                    prefix = "",
                                    suffix = "") {
  .check_required_cols(airr_df, sequence_id_col, "airr_df")
  out <- .as_data_frame(airr_df)
  ids <- as.character(out[[sequence_id_col]])

  # str_match makes a matrix with the full match in the first column and capture groups in subsequent columns
  match_matrix <- stringr::str_match(ids, regex)
  if (ncol(match_matrix) > 1) {
    barcode <- match_matrix[, 2]
  } else {
    barcode <- match_matrix[, 1]
  }
  out[[barcode_col]] <- ifelse(is.na(barcode), NA_character_, paste0(prefix, barcode, suffix))
  out
}

#' Add tissue-position coordinates to airr.tsv style data frames
#'
#' @param airr_df AIRR-style data frame containing a barcode column.
#' @param positions Tissue positions data frame or CSV/parquet path.
#' @param barcode_col Column in `airr_df` to match.
#' @param position_barcode_col Barcode column in `positions`.
#' @param x_col,y_col Coordinate columns to copy or rename to `x` and `y`.
#' @param keep_position_cols Whether to retain all position columns.
#' @return Data frame with `x` and `y` columns.
#' @export
addTissuePositions <- function(airr_df,
                                 positions,
                                 barcode_col = "barcode_8um",
                                 position_barcode_col = "barcode",
                                 x_col = "pxl_col_in_fullres",
                                 y_col = "pxl_row_in_fullres",
                                 keep_position_cols = TRUE) {
  .check_required_cols(airr_df, barcode_col, "airr_df")
  if (is.character(positions) && length(positions) == 1) {
    positions <- readBarcodeMapping(positions)
  }
  .check_required_cols(positions, c(position_barcode_col, x_col, y_col), "positions")

  out <- dplyr::left_join(
    .as_data_frame(airr_df),
    .as_data_frame(positions),
    by = stats::setNames(position_barcode_col, barcode_col)
  )
  out$x <- out[[x_col]]
  out$y <- out[[y_col]]

  if (!keep_position_cols) {
    keep <- union(colnames(airr_df), c("x", "y"))
    out <- out[, keep, drop = FALSE]
  }
  out
}
