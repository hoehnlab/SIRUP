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
#' @param shared_color_col Optional column whose shared values should receive
#'   the same tip color within each clone, for example a spatial barcode column.
#' @param palette Optional color vector. Defaults to `pals::glasbey()` when
#'   available, otherwise base R qualitative colors.
#' @return Data frame with `tip_color_id` and `tip_color` columns.
#' @export
assignTipColors <- function(clones,
                              clone_ids = NULL,
                              clone_col = "clone_id",
                              id_col = "sequence_id",
                              shared_color_col = NULL,
                              palette = NULL) {
  .check_required_cols(clones, c(clone_col, id_col), "clones")
  if (!is.null(shared_color_col)) {
    .check_required_cols(clones, shared_color_col, "clones")
  }
  out <- .as_data_frame(clones)
  if (!is.null(clone_ids)) {
    out <- out[out[[clone_col]] %in% clone_ids, , drop = FALSE]
  }
  out <- out[!is.na(out[[clone_col]]) & !is.na(out[[id_col]]), , drop = FALSE]
  if (is.null(shared_color_col)) {
    out <- out %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(clone_col))) %>%
      dplyr::arrange(.data[[id_col]], .by_group = TRUE) %>%
      dplyr::mutate(tip_color_id = as.character(dplyr::row_number())) %>%
      dplyr::ungroup()
  } else {
    out <- out %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(clone_col))) %>%
      dplyr::arrange(.data[[shared_color_col]], .data[[id_col]], .by_group = TRUE) %>%
      dplyr::mutate(tip_color_id = as.character(dplyr::dense_rank(.data[[shared_color_col]]))) %>%
      dplyr::ungroup()
  }

  n_colors <- max(as.integer(out$tip_color_id), na.rm = TRUE)
  if (is.null(palette)) {
    if (requireNamespace("pals", quietly = TRUE)) {
      palette <- .sirap_extend_palette(pals::glasbey(32), n_colors)
    } else {
      palette <- grDevices::hcl.colors(max(12, n_colors), palette = "Dark 3")
    }
  } else if (length(palette) < n_colors) {
    palette <- .sirap_extend_palette(palette, n_colors)
  }
  out$tip_color <- palette[as.integer(out$tip_color_id)]
  out
}

.sirap_extend_palette <- function(palette, n) {
  palette <- unique(as.character(stats::na.omit(palette)))
  if (length(palette) == 0L) {
    palette <- grDevices::hcl.colors(max(2L, n), palette = "Dark 3")
  } else if (length(palette) == 1L && n > 1L) {
    palette <- unique(c(palette, grDevices::hcl.colors(max(2L, n), palette = "Dark 3")))
  }
  if (length(palette) >= n) {
    return(palette[seq_len(n)])
  }

  rgb <- grDevices::col2rgb(palette)
  averaged <- character()
  for (i in seq_len(ncol(rgb))) {
    for (j in seq_len(ncol(rgb))) {
      if (i >= j) {
        next
      }
      avg <- round(rowMeans(cbind(rgb[, i], rgb[, j])))
      averaged <- c(averaged, grDevices::rgb(avg[1], avg[2], avg[3], maxColorValue = 255))
      if (length(unique(c(palette, averaged))) >= n) {
        return(unique(c(palette, averaged))[seq_len(n)])
      }
    }
  }

  grDevices::colorRampPalette(palette)(n)
}
