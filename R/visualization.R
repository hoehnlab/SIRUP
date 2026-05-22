#' Plot V/J/junction length QC histograms
#'
#' @param airr_df AIRR-style data frame.
#' @param title Optional plot title prefix.
#' @return A ggplot object.
#' @export
plotLengthQc <- function(airr_df, title = NULL) {
  .check_required_cols(airr_df, c("v_sequence_end", "v_sequence_start", "j_sequence_end", "j_sequence_start", "junction_length"), "airr_df")
  out <- .as_data_frame(airr_df)
  out$vlength <- out$v_sequence_end - out$v_sequence_start
  out$jlength <- out$j_sequence_end - out$j_sequence_start
  long <- tidyr::pivot_longer(
    out,
    cols = dplyr::all_of(c("vlength", "jlength", "junction_length")),
    names_to = "metric",
    values_to = "length"
  )
  ggplot2::ggplot(long, ggplot2::aes(x = .data$length)) +
    ggplot2::geom_histogram(binwidth = 5, na.rm = TRUE) +
    ggplot2::facet_wrap(~metric, scales = "free_x") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title, x = "Length", y = "Count")
}

#' Plot clone locations from AIRR coordinates
#'
#' This is a lightweight ggplot alternative to Seurat `SpatialFeaturePlot()`.
#'
#' @param clones Clone-assigned AIRR data frame with coordinates.
#' @param clone_id Clone to highlight; `NULL` plots all clone-assigned rows.
#' @param x_col,y_col Coordinate columns.
#' @param clone_col Clone ID column.
#' @param color_col Optional color/shape annotation column.
#' @return A ggplot object.
#' @export
plotCloneSpatialPoints <- function(clones,
                                      clone_id = NULL,
                                      x_col = "x",
                                      y_col = "y",
                                      clone_col = "clone_id",
                                      color_col = NULL) {
  .check_required_cols(clones, c(x_col, y_col, clone_col), "clones")
  out <- .as_data_frame(clones)
  if (!is.null(clone_id)) {
    out <- out[out[[clone_col]] %in% clone_id, , drop = FALSE]
  } else {
    out <- out[!is.na(out[[clone_col]]), , drop = FALSE]
  }
  mapping <- ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]])
  if (!is.null(color_col)) {
    .check_required_cols(out, color_col, "clones")
    mapping$colour <- rlang::expr(.data[[!!color_col]])
  } else {
    mapping$colour <- rlang::expr(.data[[!!clone_col]])
  }
  ggplot2::ggplot(out, mapping) +
    ggplot2::geom_point(alpha = 0.85, size = 1.8) +
    ggplot2::coord_fixed() +
    ggplot2::scale_y_reverse() +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = x_col, y = y_col, colour = ifelse(is.null(color_col), clone_col, color_col))
}

#' Add a clone assay to a Seurat spatial object
#'
#' @param seurat_obj Seurat object.
#' @param clone_matrix Clone-by-barcode count matrix.
#' @param assay_name Name for the new assay.
#' @param default_assay Whether to make the clone assay the default.
#' @return Modified Seurat object.
#' @export
addCloneAssay <- function(seurat_obj,
                            clone_matrix,
                            assay_name = "VDJ_Clonal",
                            default_assay = TRUE) {
  .require_namespace("SeuratObject")
  seurat_obj[[assay_name]] <- SeuratObject::CreateAssay5Object(counts = clone_matrix)
  if (default_assay) {
    SeuratObject::DefaultAssay(seurat_obj) <- assay_name
  }
  seurat_obj
}

#' Plot tissue coordinates with receptor locations overlaid
#'
#' @param tissue_positions Data frame of tissue/barcode coordinates. If
#'   `gene_col` is supplied, points are colored by expression.
#' @param receptors Data frame of BCR/TCR locations.
#' @param gene_col Optional column in `tissue_positions` to color as expression.
#' @param x_col,y_col Coordinate columns in `tissue_positions`.
#' @param receptor_x_col,receptor_y_col Coordinate columns in `receptors`.
#' @param receptor_color Outline color for receptor circles.
#' @param receptor_size Size for receptor circles.
#' @param title Optional plot title.
#' @return A ggplot object.
#' @export
plotReceptorOverlay <- function(tissue_positions,
                                  receptors,
                                  gene_col = NULL,
                                  x_col = "x",
                                  y_col = "y",
                                  receptor_x_col = x_col,
                                  receptor_y_col = y_col,
                                  receptor_color = "black",
                                  receptor_size = 1.6,
                                  title = NULL) {
  .check_required_cols(tissue_positions, c(x_col, y_col), "tissue_positions")
  .check_required_cols(receptors, c(receptor_x_col, receptor_y_col), "receptors")
  tissue_positions <- .as_data_frame(tissue_positions)
  receptors <- .as_data_frame(receptors)

  if (!is.null(gene_col)) {
    .check_required_cols(tissue_positions, gene_col, "tissue_positions")
    base <- ggplot2::ggplot(tissue_positions, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]], color = .data[[gene_col]])) +
      ggplot2::geom_point(size = 0.7, alpha = 0.9) +
      ggplot2::scale_color_viridis_c(option = "magma")
  } else {
    base <- ggplot2::ggplot(tissue_positions, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]])) +
      ggplot2::geom_point(color = "grey85", size = 0.7, alpha = 0.8)
  }

  base +
    ggplot2::geom_point(
      data = receptors,
      ggplot2::aes(x = .data[[receptor_x_col]], y = .data[[receptor_y_col]]),
      inherit.aes = FALSE,
      shape = 21,
      fill = NA,
      color = receptor_color,
      stroke = 0.45,
      size = receptor_size
    ) +
    ggplot2::coord_fixed() +
    ggplot2::scale_y_reverse() +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title, x = x_col, y = y_col, color = gene_col)
}

#' Plot a clone with tip colors on tissue coordinates
#'
#' @param colored_clones Output from `assignTipColors()`.
#' @param clone_id Clone ID to plot.
#' @param tissue_positions Optional background tissue coordinates.
#' @param x_col,y_col Coordinate columns.
#' @param clone_col Clone ID column.
#' @return A ggplot object.
#' @export
plotColoredCloneSpatial <- function(colored_clones,
                                       clone_id,
                                       tissue_positions = NULL,
                                       x_col = "x",
                                       y_col = "y",
                                       clone_col = "clone_id") {
  .check_required_cols(colored_clones, c(clone_col, x_col, y_col, "tip_color"), "colored_clones")
  clone_df <- .as_data_frame(colored_clones)
  clone_df <- clone_df[clone_df[[clone_col]] == clone_id, , drop = FALSE]

  if (!is.null(tissue_positions)) {
    .check_required_cols(tissue_positions, c(x_col, y_col), "tissue_positions")
    p <- ggplot2::ggplot(.as_data_frame(tissue_positions), ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]])) +
      ggplot2::geom_point(color = "grey88", size = 0.55, alpha = 0.7)
  } else {
    p <- ggplot2::ggplot()
  }

  p +
    ggplot2::geom_point(
      data = clone_df,
      ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]], fill = .data$tip_color),
      shape = 21,
      color = "black",
      stroke = 0.25,
      size = 2.4
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::coord_fixed() +
    ggplot2::scale_y_reverse() +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = paste("Clone", clone_id), x = x_col, y = y_col)
}
