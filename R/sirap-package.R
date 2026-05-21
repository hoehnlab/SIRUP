#' SIRAP: Spatial Immune Receptor Analysis and Processing
#'
#' An R interface to the rad command-line demultiplexer plus a consolidated
#' toolkit for long-read spatial immune receptor sequencing: AIRR-style
#' rearrangement tables with spatial barcodes and coordinates, Visium HD
#' barcode mapping, pRESTO/IgBLAST shell pipeline generation, clonal
#' summaries, Seurat clone assays, tree plotting, int64-packed DNA
#' primitives, and PAF parsing utilities.
#'
#' @keywords internal
#' @useDynLib SIRAP, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom data.table :=
"_PACKAGE"
