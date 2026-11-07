sirup_test_seurat <- function(counts, ...) {
  if (!inherits(counts, "Matrix")) {
    counts <- Matrix::Matrix(counts, sparse = TRUE)
  }
  counts <- methods::as(methods::as(counts, "generalMatrix"), "CsparseMatrix")
  SeuratObject::CreateSeuratObject(counts = counts, ...)
}
