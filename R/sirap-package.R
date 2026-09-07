#' SIRAP: Spatial Immune Receptor Analysis and Processing
#'
#' A native R interface to the RAD demultiplexer through `rrad`, with optional
#' legacy command-line wrappers, plus a consolidated toolkit for long-read
#' spatial immune receptor sequencing: AIRR-style
#' rearrangement tables with spatial IDs and coordinates, Visium HD
#' spatial ID mapping, pRESTO/IgBLAST shell pipeline generation, clonal
#' summaries, lightweight SIRAP containers pairing spatial Seurat data with
#' tibble-backed receptor tables, tree plotting,
#' int64-packed DNA primitives, and PAF parsing utilities.
#'
#' @keywords internal
#' @useDynLib SIRAP, .registration = TRUE
#' @import bit64
#' @import ggalign
#' @import pbapply
#' @importFrom grDevices hcl.colors
#' @importFrom Rcpp evalCpp
#' @importFrom data.table :=
#' @importFrom stats end setNames start
#' @importFrom utils head
"_PACKAGE"

utils::globalVariables(c(
  ".",
  "..curv..",
  ".N",
  ".SD",
  ".x",
  "barcode",
  "bin",
  "count",
  "count_16um",
  "count_8um",
  "d",
  "gene_i",
  "gene_j",
  "grp",
  "i",
  "i.end",
  "i.start",
  "in_sr_filtered_16um",
  "in_sr_filtered_2um",
  "in_sr_filtered_8um",
  "in_tissue",
  "j",
  "label",
  "log1p_ncpm",
  "log1p_ncpm_16um",
  "log1p_ncpm_8um",
  "n_spots_16um",
  "ncpm",
  "ncpm_16um",
  "ncpm_8um",
  "ovlen",
  "rad_valid_16um",
  "rest",
  "spaceranger_16um",
  "spaceranger_2um",
  "spaceranger_8um",
  "value",
  "w",
  "x",
  "x_16um",
  "x_2um",
  "x_8um",
  "x_bc",
  "xend",
  "xmax",
  "xmin",
  "y",
  "y_16um",
  "y_2um",
  "y_8um",
  "y_bc"
))
