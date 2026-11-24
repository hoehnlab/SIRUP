#' Save and read SIRUP objects
#'
#' `saveSIRUP` writes a SIRUP object to one RDS file. `readSIRUP` reads the file
#' back.
#'
#' @details
#' `saveSIRUP` validates the object. It stores the versions of SIRUP and
#' SeuratObject that write the file. 
#'
#' `readSIRUP` does these steps:
#' \enumerate{
#'   \item It converts an object that was saved under the old package name,
#'     SIRAP (class `sirap`), to the current names.
#'   \item It upgrades an older object to the current schema.
#'   \item It gives a warning if the versions of SIRUP or SeuratObject that
#'     wrote the file differ from the installed versions. Version 0.0.0.9001
#'     is the development version that became release 1.0.0. It gives no
#'     warning for that pair.
#'   \item It validates the object, then returns it.
#' }
#'
#' @param x A SIRUP object from [createSIRUP()].
#' @param file Path of the RDS file. A `.sirup.rds` extension is a good choice,
#'   but any path works.
#' @param ... Other arguments for [saveRDS()], for example `compress`.
#'
#' @return `saveSIRUP` returns `file`, invisibly. `readSIRUP` returns the
#'   validated SIRUP object.
#' @seealso [createSIRUP()], [validateSIRUP()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#'
#' # Save the object to a temporary file, then read it back.
#' path <- saveSIRUP(x, tempfile(fileext = ".sirup.rds"))
#' x2 <- readSIRUP(path)
#' x2
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
  if (!is.null(saved_pkg) && !.sirup_same_release(saved_pkg, current_pkg)) {
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

# TRUE when two SIRUP versions are the same code. Version 0.0.0.9001 is the
# development version that was released as 1.0.0, so the pair counts as equal.
.sirup_same_release <- function(saved, current) {
  saved <- as.character(saved)
  current <- as.character(current)
  identical(saved, current) ||
    (identical(saved, "0.0.0.9001") && identical(current, "1.0.0"))
}

# Installed SeuratObject version, or NA when SeuratObject is unavailable.
.sirup_seurat_version <- function() {
  tryCatch(
    as.character(utils::packageVersion("SeuratObject")),
    error = function(e) NA_character_
  )
}
