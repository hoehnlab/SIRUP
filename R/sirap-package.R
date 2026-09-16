#' SIRAP: Spatial Immune Receptor Analysis and Processing
#'
#' A consolidated toolkit for long-read spatial immune receptor sequencing.
#' SIRAP takes demultiplexed reads from the companion RRAD package and runs
#' an end-to-end IgBLAST and pRESTO pipeline to produce AIRR-style
#' rearrangement tables with spatial IDs and coordinates. It also provides
#' Visium HD spatial ID mapping, clonal summaries, lightweight SIRAP
#' containers pairing spatial Seurat data with tibble-backed receptor tables,
#' tree plotting, polygon annotation, and spatial statistics.
#'
#' @keywords internal
#' @import pbapply
#' @importFrom grDevices hcl.colors
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
  "d",
  "grp",
  "i",
  "i.end",
  "i.start",
  "j",
  "label",
  "value",
  "w",
  "x",
  "xend",
  "xmax",
  "xmin",
  "y"
))
