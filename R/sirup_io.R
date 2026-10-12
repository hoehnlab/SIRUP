#' Save and load SIRUP objects
#'
#' `saveSIRUP()` writes a validated SIRUP object to a single RDS file, and
#' `readSIRUP()` reads it back. Unlike a bare [saveRDS()]/[readRDS()] round
#' trip, `saveSIRUP()` validates the object and records the SIRUP and
#' 'SeuratObject' package versions used to write it, and `readSIRUP()` upgrades
#' older objects to the current schema, warns when the writing versions differ
#' from the versions currently installed, and re-validates before returning.
#'
#' @param x A SIRUP object created by [createSIRUP()].
#' @param file Path to the RDS file. A `.sirup.rds` extension is recommended
#'   but any path is accepted.
#' @param ... Additional arguments passed to [saveRDS()] (for example
#'   `compress`).
#'
#' @return `saveSIRUP()` returns `file` invisibly. `readSIRUP()` returns the
#'   validated SIRUP object.
#'
#' @examples
#' \dontrun{
#' sirup <- createSIRUP(spatial = seurat_object, receptors = cloned_receptors)
#' path <- saveSIRUP(sirup, tempfile(fileext = ".sirup.rds"))
#' sirup2 <- readSIRUP(path)
#' }
#' @name sirup-io
NULL

#' @rdname sirup-io
#' @export
saveSIRUP <- function(x, file, ...) {
  validateSIRUP(x)
  if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file)) {
    rlang::abort("`file` must be a single non-empty file path.")
  }
  attr(x, "sirup_pkg_version") <- as.character(utils::packageVersion("SIRUP"))
  attr(x, "sirup_seurat_version") <- .sirup_seurat_version()
  saveRDS(x, file = file, ...)
  invisible(file)
}

#' @rdname sirup-io
#' @export
readSIRUP <- function(file) {
  if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file)) {
    rlang::abort("`file` must be a single non-empty file path.")
  }
  if (!file.exists(file)) {
    rlang::abort(paste0("`file` does not exist: ", file))
  }
  x <- readRDS(file)
  x <- .sirup_from_legacy_sirap(x)
  if (!isSIRUP(x)) {
    rlang::abort("`file` does not contain a SIRUP object.")
  }
  saved_pkg <- attr(x, "sirup_pkg_version")
  saved_seurat <- attr(x, "sirup_seurat_version")

  x <- .sirup_upgrade(x)

  current_pkg <- as.character(utils::packageVersion("SIRUP"))
  if (!is.null(saved_pkg) && !identical(as.character(saved_pkg), current_pkg)) {
    cli::cli_warn(c(
      "This object was saved with SIRUP {saved_pkg} but SIRUP {current_pkg} is installed.",
      "i" = "It has been upgraded to the current schema; review results if behavior differs."
    ))
  }
  current_seurat <- .sirup_seurat_version()
  if (!is.null(saved_seurat) && !is.na(saved_seurat) &&
      !is.na(current_seurat) && !identical(as.character(saved_seurat), current_seurat)) {
    cli::cli_warn(
      "This object's Seurat data was saved with SeuratObject {saved_seurat} but SeuratObject {current_seurat} is installed."
    )
  }

  validateSIRUP(x)
  x
}

# The package was called SIRAP before it was renamed to SIRUP. Objects saved
# under the old name have class "sirap" and "sirap_*" version attributes.
# Stored analyses can have "sirap_*" classes too. Convert all of them to the
# current names. Return other objects unchanged.
.sirup_from_legacy_sirap <- function(x) {
  if (!inherits(x, "sirap") || inherits(x, "sirup")) {
    return(x)
  }
  class(x)[class(x) == "sirap"] <- "sirup"
  if (is.list(x$analyses)) {
    for (i in seq_along(x$analyses)) {
      analysis_class <- class(x$analyses[[i]])
      if (any(startsWith(analysis_class, "sirap_"))) {
        class(x$analyses[[i]]) <- sub("^sirap_", "sirup_", analysis_class)
      }
    }
  }
  for (field in c("pkg_version", "seurat_version")) {
    old_attr <- paste0("sirap_", field)
    if (!is.null(attr(x, old_attr))) {
      attr(x, paste0("sirup_", field)) <- attr(x, old_attr)
      attr(x, old_attr) <- NULL
    }
  }
  x
}

# Installed SeuratObject version, or NA when SeuratObject is unavailable.
.sirup_seurat_version <- function() {
  tryCatch(
    as.character(utils::packageVersion("SeuratObject")),
    error = function(e) NA_character_
  )
}
