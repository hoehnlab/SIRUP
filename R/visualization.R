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
