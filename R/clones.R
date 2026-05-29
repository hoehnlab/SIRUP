#' Build a clone-by-spatial-barcode count matrix
#'
#' @param clones Clone-assigned AIRR data frame.
#' @param clone_col,barcode_col Column names.
#' @param sparse Return a sparse Matrix when the optional `Matrix` package exists.
#' @return Clone-by-barcode count matrix.
#' @export
buildCloneCountMatrix <- function(clones,
                                     clone_col = "clone_id",
                                     barcode_col = "barcode",
                                     sparse = TRUE) {
  .check_required_cols(clones, c(clone_col, barcode_col), "clones")
  x <- clones[!is.na(clones[[clone_col]]) & !is.na(clones[[barcode_col]]), , drop = FALSE]
  mat <- table(x[[clone_col]], x[[barcode_col]]) %>% as.matrix()
  if (sparse && requireNamespace("Matrix", quietly = TRUE)) {
    mat <- Matrix::Matrix(mat, sparse = TRUE)
  }
  mat
}

#' Assign stable colors within selected clones
#'
#' Adds a per-tip color identifier and color value that can be used consistently
#' in tree plots and matched spatial plots.
#'
#' @param clones Clone-assigned AIRR data frame.
#' @param clone_ids Optional clone IDs to retain before assigning colors.
#' @param clone_col,id_col Column names for clone and sequence/tip identifiers.
#' @param palette Optional color vector. Defaults to `pals::glasbey()` when
#'   available, otherwise base R qualitative colors.
#' @return Data frame with `tip_color_id` and `tip_color` columns.
#' @export
assignTipColors <- function(clones,
                              clone_ids = NULL,
                              clone_col = "clone_id",
                              id_col = "sequence_id",
                              palette = NULL) {
  .check_required_cols(clones, c(clone_col, id_col), "clones")
  out <- .as_data_frame(clones)
  if (!is.null(clone_ids)) {
    out <- out[out[[clone_col]] %in% clone_ids, , drop = FALSE]
  }
  out <- out[!is.na(out[[clone_col]]) & !is.na(out[[id_col]]), , drop = FALSE]
  out <- out %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(clone_col))) %>%
    dplyr::mutate(tip_color_id = as.character(dplyr::row_number())) %>%
    dplyr::ungroup()

  n_colors <- max(as.integer(out$tip_color_id), na.rm = TRUE)
  if (is.null(palette)) {
    if (requireNamespace("pals", quietly = TRUE)) {
      palette <- pals::glasbey(max(32, n_colors))
    } else {
      palette <- grDevices::hcl.colors(max(12, n_colors), palette = "Dark 3")
    }
  }
  out$tip_color <- palette[as.integer(out$tip_color_id)]
  out
}
