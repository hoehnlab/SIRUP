sirup_test_seurat <- function(counts, ...) {
  testthat::skip_if_not_installed("SeuratObject")
  testthat::skip_if_not_installed("Matrix")
  if (!inherits(counts, "Matrix")) {
    counts <- Matrix::Matrix(counts, sparse = TRUE)
  }
  counts <- methods::as(methods::as(counts, "generalMatrix"), "CsparseMatrix")
  SeuratObject::CreateSeuratObject(counts = counts, ...)
}
