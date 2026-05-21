# Build a simple consensus sequence from a vector of aligned sequences (MSA strings).
# ambiguity: character to use when no single base is dominant at a position.

#' Generate a consensus sequence from aligned strings
#' @param seqs character vector of equal-length aligned sequences
#' @param ambiguity character used when no consensus base (default "N")
#' @return consensus string
#' @export
msa_consensus <- function(seqs, ambiguity = "N") {
  if (length(seqs) == 0) return("")
  n <- nchar(seqs[1])
  if (any(nchar(seqs) != n)) {
    stop("All sequences must be the same length (aligned).")
  }
  mat <- do.call(rbind, strsplit(seqs, ""))
  apply(mat, 2, function(col) {
    tab <- sort(table(col), decreasing = TRUE)
    if (length(tab) == 0) return(ambiguity)
    # if top is unique, use it; else ambiguity
    if (length(tab) > 1 && tab[1] == tab[2]) ambiguity else names(tab)[1]
  }) |> paste0(collapse = "")
}
