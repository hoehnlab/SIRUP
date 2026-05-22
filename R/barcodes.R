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

#' Swap row/column order in a Visium HD square barcode
#'
#' Converts barcodes such as `s_002um_12_34-1` to `s_002um_34_12-1`.
#' This mirrors the post-consensus correction used in the original scripts.
#'
#' @param x Character vector of Visium HD square barcodes.
#' @return Character vector with row/column components swapped.
#' @export
swapVisiumHdBarcode <- function(x) {
  sub("^(.*_)(\\d+)_(\\d+)(-\\d+)$", "\\1\\3_\\2\\4", x)
}

#' Add Visium HD spatial barcodes to AIRR records
#'
#' @param airr_df AIRR-style data frame.
#' @param barcode_map Data frame with square barcode mapping columns.
#' @param sequence_id_col Column containing read or consensus identifiers.
#' @param mode How to extract the 2 um barcode. `raw_suffix` extracts the last
#'   `row_col` suffix from `sequence_id` and wraps it as `s_002um_<row_col>-1`.
#'   `presto_header` extracts the value after the first comma, matching the
#'   consensus IDs produced by the pRESTO header merge step.
#' @param map_from,map_to_8um,map_to_16um Column names in `barcode_map`.
#' @param swap Whether to swap row/column order before matching to `barcode_map`.
#' @return Input data with `barcode_2um`, optional swapped barcode, `barcode_8um`,
#'   and `barcode_16um` columns.
#' @export
addSpatialBarcodes <- function(airr_df,
                                 barcode_map,
                                 sequence_id_col = "sequence_id",
                                 mode = c("raw_suffix", "presto_header"),
                                 map_from = "square_002um",
                                 map_to_8um = "square_008um",
                                 map_to_16um = "square_016um",
                                 swap = FALSE) {
  mode <- match.arg(mode)
  .check_required_cols(airr_df, sequence_id_col, "airr_df")
  .check_required_cols(barcode_map, map_from, "barcode_map")

  out <- .as_data_frame(airr_df)
  barcode_map <- .as_data_frame(barcode_map)
  ids <- as.character(out[[sequence_id_col]])

  if (mode == "raw_suffix") {
    suffix <- stringr::str_extract(ids, "[0-9]+_[0-9]+$")
    out$barcode_2um <- paste0("s_002um_", suffix, "-1")
  } else {
    split_ids <- stringr::str_split_fixed(ids, ",", 2)
    out$barcode_2um <- split_ids[, 2]
    out$barcode_2um[out$barcode_2um == ""] <- NA_character_
  }

  match_key <- out$barcode_2um
  if (swap) {
    out$barcode_2um_swapped <- swapVisiumHdBarcode(out$barcode_2um)
    match_key <- out$barcode_2um_swapped
  }

  keep <- !is.na(match_key) & match_key %in% barcode_map[[map_from]]
  out <- out[keep, , drop = FALSE]
  match_key <- match_key[keep]
  idx <- match(match_key, barcode_map[[map_from]])

  if (map_to_8um %in% colnames(barcode_map)) {
    out$barcode_8um <- barcode_map[[map_to_8um]][idx]
  }
  if (map_to_16um %in% colnames(barcode_map)) {
    out$barcode_16um <- barcode_map[[map_to_16um]][idx]
  }

  out
}

#' Add tissue-position coordinates to spatial AIRR records
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
