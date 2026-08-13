#' Save and load SIRAP objects
#'
#' `saveSIRAP()` writes a validated SIRAP object to a single RDS file, and
#' `readSIRAP()` reads it back. Unlike a bare [saveRDS()]/[readRDS()] round
#' trip, `saveSIRAP()` validates the object and records the SIRAP and
#' 'SeuratObject' package versions used to write it, and `readSIRAP()` upgrades
#' older objects to the current schema, warns when the writing versions differ
#' from the versions currently installed, and re-validates before returning.
#'
#' @param x A SIRAP object created by [createSIRAP()].
#' @param file Path to the RDS file. A `.sirap.rds` extension is recommended
#'   but any path is accepted.
#' @param ... Additional arguments passed to [saveRDS()] (for example
#'   `compress`).
#'
#' @return `saveSIRAP()` returns `file` invisibly. `readSIRAP()` returns the
#'   validated SIRAP object.
#'
#' @examples
#' \dontrun{
#' sirap <- createSIRAP(spatial = seurat_object, receptors = cloned_receptors)
#' path <- saveSIRAP(sirap, tempfile(fileext = ".sirap.rds"))
#' sirap2 <- readSIRAP(path)
#' }
#' @name sirap-io
NULL

#' @rdname sirap-io
#' @export
saveSIRAP <- function(x, file, ...) {
  validateSIRAP(x)
  if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file)) {
    rlang::abort("`file` must be a single non-empty file path.")
  }
  attr(x, "sirap_pkg_version") <- as.character(utils::packageVersion("SIRAP"))
  attr(x, "sirap_seurat_version") <- .sirap_seurat_version()
  saveRDS(x, file = file, ...)
  invisible(file)
}

#' @rdname sirap-io
#' @export
readSIRAP <- function(file) {
  if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file)) {
    rlang::abort("`file` must be a single non-empty file path.")
  }
  if (!file.exists(file)) {
    rlang::abort(paste0("`file` does not exist: ", file))
  }
  x <- readRDS(file)
  if (!isSIRAP(x)) {
    rlang::abort("`file` does not contain a SIRAP object.")
  }
  saved_pkg <- attr(x, "sirap_pkg_version")
  saved_seurat <- attr(x, "sirap_seurat_version")

  x <- .sirap_upgrade(x)

  current_pkg <- as.character(utils::packageVersion("SIRAP"))
  if (!is.null(saved_pkg) && !identical(as.character(saved_pkg), current_pkg)) {
    cli::cli_warn(c(
      "This object was saved with SIRAP {saved_pkg} but SIRAP {current_pkg} is installed.",
      "i" = "It has been upgraded to the current schema; review results if behavior differs."
    ))
  }
  current_seurat <- .sirap_seurat_version()
  if (!is.null(saved_seurat) && !is.na(saved_seurat) &&
      !is.na(current_seurat) && !identical(as.character(saved_seurat), current_seurat)) {
    cli::cli_warn(
      "This object's Seurat data was saved with SeuratObject {saved_seurat} but SeuratObject {current_seurat} is installed."
    )
  }

  validateSIRAP(x)
  x
}

# Installed SeuratObject version, or NA when SeuratObject is unavailable.
.sirap_seurat_version <- function() {
  tryCatch(
    as.character(utils::packageVersion("SeuratObject")),
    error = function(e) NA_character_
  )
}
