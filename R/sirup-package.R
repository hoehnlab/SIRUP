#' SIRUP: Spatial Immune Receptor Utilities Package
#'
#' SIRUP analyzes B cell and T cell receptor sequences from long-read spatial
#' sequencing. It turns demultiplexed reads into standard-formatted AIRR data with
#' spatial IDs, joins the AIRR data to Seurat spatial data, and includes tools for 
#' examining and analyzing clones, lineage trees, tissue regions, and more in *space*.
#'
#' @details
#' SIRUP takes demultiplexed reads from the companion RRAD package. The main
#' steps are:
#' \enumerate{
#'   \item Process the reads. [processSpatialSample()] intakes demultiplexed reads
#'     and writes an AIRR table. 
#'   \item Build a SIRUP object. [createSIRUP()] joins a Seurat object and the
#'     AIRR table. Use [saveSIRUP()] and [readSIRUP()] to keep the object in
#'     a file.
#'   \item Draw and compare tissue regions. [selectSpatialPolygons()] and
#'     [addPolygons()] store polygons. [spatialCloneDominance()] and
#'     [polygonDiversity()] summarize the receptors in each polygon.
#'   \item Test spatial patterns. [spatialMoran()], [spatialCloneJoinCount()],
#'     and [spatialRipleyCrossK()] test if genes or clones cluster in space.
#'   \item Work with AIRR data as you typically would using [updateReceptors()]. 
#'     Do the same with the spatial gene expression data with [updateSpatial()]. 
#'     Then work with both together.
#' }
#' The vignettes show each step: `vignette("SIRUP", package = "SIRUP")`,
#' `vignette("processing", package = "SIRUP")`,
#' `vignette("clones-and-trees", package = "SIRUP")`,
#' `vignette("regions", package = "SIRUP")`,
#' `vignette("spatial-statistics", package = "SIRUP")`, and
#' `vignette("multisample", package = "SIRUP")`.
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
