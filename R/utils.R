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
# by-reference behavior and printing do not leak into SIRUP internals.
.as_data_frame <- function(x) {
  if (tibble::is_tibble(x)) {
    x
  } else if (data.table::is.data.table(x)) {
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

# Example files for the guides on the SIRUP website. The files are not part of
# the package. Each entry is a file name with its MD5 sum. The function checks
# the sum after each download, so a broken download never looks like good data.
.sirup_example_files <- c(
  mg_thymus_bulk_blood_month0_1.tsv.gz = "32a2b5082c36db2c8f01a91dbb895a55",
  mg_thymus_bulk_blood_month0_2.tsv.gz = "4c6c8fcb9050363857073bbe34406544",
  mg_thymus_bulk_blood_month12_1.tsv.gz = "c002fe40968a6772dd90e623797dacb4",
  mg_thymus_bulk_blood_month12_2.tsv.gz = "a2dd5a1d9663c4064d541e596a6d9c1d",
  mg_thymus_bulk_thymus_1.tsv.gz = "7b98b5f6db3fd313d80e69c83ee82971",
  mg_thymus_bulk_thymus_2.tsv.gz = "e9cb1ad9e56659d6d622dac94d219b1b",
  mg_thymus_cell_types.csv.gz = "aff1d673873863ffe599fe09d285cbd0",
  mg_thymus_slice1_bcr.tsv.gz = "8160318f04f9e4d429e2df7715050c4a",
  mg_thymus_slice1_gc_polygons.csv = "dd81ab1de77af3e4acde62783fba00ea",
  mg_thymus_slice1_spatial.rds = "7ba68b7aa2b15d79f307d6c10fbe2444",
  mg_thymus_slice2_bcr.tsv.gz = "5e321c97072f5385404e035993a2b3b1",
  mg_thymus_slice2_gc_polygons.csv = "645f0c556f6aa69a3adb735a6333cb43",
  mg_thymus_slice2_qupath_gcs.geojson = "015da4f75651f2beb49a0a6fa08fa6be",
  mg_thymus_slice2_reads.tsv.gz = "3fb3a4dfd97170ce7c779fd169ec8d52",
  mg_thymus_slice2_spatial.rds = "70af91a05c16c267e21bb98b021818dc",
  mg_thymus_slice2_tcr.tsv.gz = "3ac507853a38c88363116c062df1a4ba"
)

.sirup_example_url <- "https://raw.githubusercontent.com/hoehnlab/SIRUP/main/example-data/"

#' Download the SIRUP example data
#'
#' Download the myasthenia gravis (MG) thymus example files that the guides on
#' the SIRUP website use. The files are not part of the package, to keep the
#' package small. The function skips files that are already in `dest` and
#' checks every download against a stored MD5 sum.
#'
#' The files are:
#' * `mg_thymus_slice1_spatial.rds`, `mg_thymus_slice2_spatial.rds`: spatial
#'   Seurat objects for the two tissue sections.
#' * `mg_thymus_slice1_bcr.tsv.gz`, `mg_thymus_slice2_bcr.tsv.gz`: AIRR
#'   rearrangement tables of IGH sequences with spatial IDs.
#' * `mg_thymus_slice2_tcr.tsv.gz`, `mg_thymus_slice2_reads.tsv.gz`: TCR
#'   sequences and the demultiplexed reads of section 2.
#' * `mg_thymus_slice1_gc_polygons.csv`, `mg_thymus_slice2_gc_polygons.csv`:
#'   annotated polygons of the ectopic germinal centers.
#' * `mg_thymus_slice2_qupath_gcs.geojson`: six germinal center outlines that
#'   a user drew on section 2 in QuPath. Use it with [addQuPathPolygons()].
#' * `mg_thymus_cell_types.csv.gz`: cell-type proportions for each spot.
#' * `mg_thymus_bulk_*.tsv.gz`: bulk repertoire tables from thymus and blood.
#'
#' @param dest Folder to save the files in. The function creates it if needed.
#'   The default is a folder in the R session's temporary directory, so the
#'   files are lost when you close R. Give a permanent folder to keep them.
#' @param files Names of the files to download. If `NULL` (default), the
#'   function downloads all files. Use this to download only what a guide needs.
#' @param overwrite Download a file again even if it is already in `dest`
#'   (default `FALSE`)?
#' @param base_url Web address of the folder that holds the files. Change it only
#'   to use a mirror. A `file://` address of a local folder also works. The
#'   default is the option `SIRUP.example_url` if you set it.
#' @param quiet Hide the download progress messages (default `FALSE`)?
#'
#' @return The path of `dest`, invisibly. Use `file.path()` with the path and a
#'   file name to read a file.
#' @seealso [readSIRUP()]
#'
#' @examplesIf interactive()
#' extdata <- downloadSIRUPExample(
#'   files = c("mg_thymus_slice2_spatial.rds", "mg_thymus_slice2_bcr.tsv.gz")
#' )
#' thymus <- readRDS(file.path(extdata, "mg_thymus_slice2_spatial.rds"))
#' @export
downloadSIRUPExample <- function(dest = file.path(tempdir(), "SIRUP_example"),
                                 files = NULL,
                                 overwrite = FALSE,
                                 base_url = getOption("SIRUP.example_url", .sirup_example_url),
                                 quiet = FALSE) {
  if (is.null(files)) {
    files <- names(.sirup_example_files)
  }
  unknown <- setdiff(files, names(.sirup_example_files))
  if (length(unknown) > 0) {
    rlang::abort(c(
      paste0("Unknown example file(s): ", paste(unknown, collapse = ", ")),
      i = paste0("Valid names: ", paste(names(.sirup_example_files), collapse = ", "))
    ))
  }
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(dest)) {
    rlang::abort(paste0("Could not create the folder: ", dest))
  }
  if (!endsWith(base_url, "/")) {
    base_url <- paste0(base_url, "/")
  }

  old_timeout <- options(timeout = max(600, getOption("timeout")))
  on.exit(options(old_timeout), add = TRUE)

  for (file in unique(files)) {
    path <- file.path(dest, file)
    expected <- .sirup_example_files[[file]]
    if (!overwrite && file.exists(path) && unname(tools::md5sum(path)) == expected) {
      if (!quiet) cli::cli_alert_info("{.file {file}} is already in {.path {dest}}.")
      next
    }
    part <- paste0(path, ".part")
    on.exit(unlink(part), add = TRUE)
    url <- paste0(base_url, file)
    if (!quiet) cli::cli_alert_info("Downloading {.file {file}}")
    status <- tryCatch(
      suppressWarnings(utils::download.file(url, part, mode = "wb", quiet = quiet)),
      error = function(e) conditionMessage(e)
    )
    if (!identical(status, 0L) || !file.exists(part)) {
      rlang::abort(c(
        paste0("Could not download ", file, "."),
        i = paste0("URL: ", url),
        i = "Check your internet connection, then try again."
      ))
    }
    if (unname(tools::md5sum(part)) != expected) {
      rlang::abort(c(
        paste0("The downloaded file ", file, " is damaged (MD5 sum does not match)."),
        i = "Try again with `overwrite = TRUE`."
      ))
    }
    file.rename(part, path)
  }
  invisible(dest)
}
