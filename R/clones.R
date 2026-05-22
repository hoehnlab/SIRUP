#' Assign receptor clones with scoper
#'
#' @param airr_df AIRR-style data frame.
#' @param threshold Distance threshold for `scoper::hierarchicalClones()`.
#' @param productive_only Whether to keep productive rows before clone calling.
#' @param ... Passed to `scoper::hierarchicalClones()`.
#' @return Data frame with clone assignments.
#' @export
assignSpatialClones <- function(airr_df,
                                  threshold = 0.1,
                                  productive_only = TRUE,
                                  ...) {
  .require_namespace("scoper")
  out <- .as_data_frame(airr_df)
  if (productive_only && "productive" %in% colnames(out)) {
    out <- out[out$productive %in% c(TRUE, "TRUE", "T"), , drop = FALSE]
  }
  scoper::hierarchicalClones(out, threshold = threshold, summarize_clones = FALSE, ...)
}

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
  mat <- table(x[[clone_col]], x[[barcode_col]]) |> as.matrix()
  if (sparse && requireNamespace("Matrix", quietly = TRUE)) {
    mat <- Matrix::Matrix(mat, sparse = TRUE)
  }
  mat
}

#' Simplify common human IGH constant calls
#'
#' @param c_call Character vector of constant-region allele calls.
#' @return Character vector such as `IgM`, `IgG1`, or `unknown`.
#' @export
simplifyIghConstantCall <- function(c_call) {
  out <- dplyr::case_when(
    stringr::str_detect(c_call, "IGHM") ~ "IgM",
    stringr::str_detect(c_call, "IGHD") ~ "IgD",
    stringr::str_detect(c_call, "IGHG1") ~ "IgG1",
    stringr::str_detect(c_call, "IGHG2") ~ "IgG2",
    stringr::str_detect(c_call, "IGHG3") ~ "IgG3",
    stringr::str_detect(c_call, "IGHG4") ~ "IgG4",
    stringr::str_detect(c_call, "IGHA1") ~ "IgA1",
    stringr::str_detect(c_call, "IGHA2") ~ "IgA2",
    stringr::str_detect(c_call, "IGHE") ~ "IgE",
    stringr::str_detect(c_call, "IGKC") ~ "kappa",
    stringr::str_detect(c_call, "IGLC") ~ "lambda",
    TRUE ~ "unknown"
  )
  out[is.na(out)] <- "unknown"
  out
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
  out <- out |>
    dplyr::group_by(dplyr::across(dplyr::all_of(clone_col))) |>
    dplyr::mutate(tip_color_id = as.character(dplyr::row_number())) |>
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

#' Build B-cell lineage trees with dowser
#'
#' Convenience wrapper around `dowser::formatClones()` and `dowser::getTrees()`.
#' The user-requested `method = "pml"` is passed to dowser as `build = "pml"`.
#'
#' @param clones Clone-assigned AIRR data frame.
#' @param clone_ids Optional clone IDs to keep.
#' @param id Sequence/tip ID column.
#' @param traits Trait columns retained for plotting.
#' @param columns Extra columns retained in formatted clones.
#' @param method Tree-building method passed to `dowser::getTrees(build = ...)`.
#' @param nproc Number of processors.
#' @param ... Passed to `dowser::getTrees()`.
#' @return A list with `formatted_clones` and `trees`.
#' @export
buildDowserTrees <- function(clones,
                               clone_ids = NULL,
                               id = "sequence_id",
                               traits = c("tip_color_id", "tip_color"),
                               columns = NULL,
                               method = "pml",
                               nproc = 1,
                               ...) {
  .require_namespace("dowser")
  out <- .as_data_frame(clones)
  if (!is.null(clone_ids)) {
    out <- out[out$clone_id %in% clone_ids, , drop = FALSE]
  }
  keep_traits <- intersect(traits, colnames(out))
  keep_columns <- unique(c(intersect(columns, colnames(out)), keep_traits))
  formatted <- dowser::formatClones(
    out,
    id = id,
    traits = keep_traits,
    columns = keep_columns,
    filterstop = FALSE,
    collapse = TRUE,
    add_count = "consensus_count" %in% colnames(out),
    num_fields = if ("consensus_count" %in% colnames(out)) "consensus_count" else NULL
  )
  trees <- dowser::getTrees(formatted, build = method, nproc = nproc, ...)
  list(formatted_clones = formatted, trees = trees)
}
