# Assert that a data-frame-like object contains named columns before a helper
# tries to subset it, so errors point to the public argument that is wrong.
.check_required_cols <- function(x, cols, data_name = deparse(substitute(x))) {
  missing <- setdiff(cols, colnames(x))
  if (length(missing) > 0) {
    rlang::abort(
      paste0(
        data_name,
        " is missing required column(s): ",
        paste(missing, collapse = ", ")
      )
    )
  }
  invisible(x)
}

#' @importFrom rlang .data
#' @importFrom magrittr %>%
NULL

# Normalize tabular inputs used by internal helpers. Keep ordinary data frames
# as base data frames, but turn data.table inputs into tibbles so data.table's
# by-reference behavior and printing do not leak into SIRAP internals.
.as_data_frame <- function(x) {
  if (data.table::is.data.table(x)) {
    tibble::as_tibble(x)
  } else {
    as.data.frame(x, stringsAsFactors = FALSE)
  }
}

# Require optional dependencies only at the functions that need them.
.require_namespace <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    rlang::abort(paste0("Package '", pkg, "' is required for this operation."))
  }
  invisible(TRUE)
}

# Null-coalescing helper used for optional arguments and metadata fallbacks.
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
