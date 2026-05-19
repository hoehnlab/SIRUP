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
NULL

.as_data_frame <- function(x) {
  if (data.table::is.data.table(x)) {
    tibble::as_tibble(x)
  } else {
    as.data.frame(x, stringsAsFactors = FALSE)
  }
}

.require_namespace <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    rlang::abort(paste0("Package '", pkg, "' is required for this operation."))
  }
  invisible(TRUE)
}

.ensure_dir <- function(path) {
  dir <- dirname(path)
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}
