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
#' @param id Column containing read or consensus identifiers.
#' @param barcode Name of the output spatial barcode column.
#' @param regex Regular expression used to extract the barcode. If the regex
#'   contains a capture group, the first capture group is used. Otherwise the
#'   full match is used.
#' @param prefix,suffix Optional text added before or after the extracted match.
#' @return Input data with `barcode_col` added.
#' @export
extractSpatialBarcodes <- function(airr_df,
                                    id = "sequence_id",
                                    barcode = "spatial_barcode",
                                    regex,
                                    prefix = "",
                                    suffix = "") {
  .check_required_cols(airr_df, id, "airr_df")
  out <- .as_data_frame(airr_df)
  ids <- as.character(out[[id]])

  # str_match makes a matrix with the full match in the first column and capture groups in subsequent columns
  match_matrix <- stringr::str_match(ids, regex)
  if (ncol(match_matrix) > 1) {
    extracted_barcode <- match_matrix[, 2]
  } else {
    extracted_barcode <- match_matrix[, 1]
  }
  out[[barcode]] <- ifelse(is.na(extracted_barcode), NA_character_, paste0(prefix, extracted_barcode, suffix))
  out
}

#' Add tissue-position coordinates to airr.tsv style data frames
#'
#' @param airr_df AIRR-style data frame containing a barcode column.
#' @param positions Tissue positions data frame or CSV/parquet path.
#' @param barcode Column in `airr_df` to match.
#' @param position_barcode Barcode column in `positions`.
#' @param x,y Coordinate columns to copy or rename to `x` and `y`.
#' @param keep_position_cols Whether to retain all position columns.
#' @return Data frame with `x` and `y` columns.
#' @export
addTissuePositions <- function(airr_df,
                                 positions,
                                 barcode = "spatial_barcode_8um",
                                 position_barcode = "barcode",
                                 x = "pxl_col_in_fullres",
                                 y = "pxl_row_in_fullres",
                                 keep_position_cols = TRUE) {
  .check_required_cols(airr_df, barcode, "airr_df")
  if (is.character(positions) && length(positions) == 1) {
    positions <- readBarcodeMapping(positions)
  }
  .check_required_cols(positions, c(position_barcode, x, y), "positions")

  out <- dplyr::left_join(
    .as_data_frame(airr_df),
    .as_data_frame(positions),
    by = stats::setNames(position_barcode, barcode)
  )
  out$x <- out[[x]]
  out$y <- out[[y]]

  if (!keep_position_cols) {
    keep <- union(colnames(airr_df), c("x", "y"))
    out <- out[, keep, drop = FALSE]
  }
  out
}
