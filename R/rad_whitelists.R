# Native RAD whitelist resources -------------------------------------------

.rrad_download_whitelists <- function(...) {
  rrad::rad_download_whitelists(...)
}

.sirap_validate_rad_resource_core <- function(core) {
  if (!is.character(core$features) || anyNA(core$features) ||
      !("remote-resource-cache" %in% core$features)) {
    .sirap_rad_incompatible(
      "The embedded RAD core lacks the required remote-resource-cache capability."
    )
  }
  invisible(core)
}

#' Prepare RAD whitelists for native SIRAP runs
#'
#' Native [rad_demux()] and [rad_scan_wl()] runs automatically download a
#' missing large built-in whitelist the first time it is needed. The file is
#' stored in R's standard per-user cache and reused on later runs after size
#' and SHA-256 verification. No whitelist is downloaded during package
#' installation or attachment.
#'
#' Most users therefore do not need to call `rad_download_whitelists()`.
#' It is provided for explicitly preparing an offline workstation, compute
#' node, or shared cache before a run. SIRAP verifies its `rrad` compatibility
#' handshake before delegating the download.
#'
#' @param kits Optional character vector of built-in whitelist names or kit
#'   aliases. Downloadable names are `"10x_3v3"`, `"10x_3v4"`, and
#'   `"10x_5v3"`; `"10x_3v3.1"` and `"10x_3HTv3.1"` resolve to
#'   `"10x_3v3"`. `NULL` prepares every downloadable large whitelist.
#' @param force Re-download and verify the requested whitelists even when valid
#'   cached copies already exist. A failed refresh leaves a valid cached copy
#'   untouched.
#' @param quiet Suppress download progress messages.
#' @param cache_dir Optional cache base directory. The default is
#'   `getOption("rrad.cache_dir")`, then `RRAD_CACHE_DIR`, when set, otherwise
#'   `tools::R_user_dir("rrad", "cache")`.
#'
#' @return A named character vector containing the verified whitelist paths.
#' @seealso [rad_demux()], [rad_scan_wl()],
#'   [rrad::rad_download_whitelists()]
#' @export
rad_download_whitelists <- function(kits = NULL, force = FALSE, quiet = FALSE,
                                    cache_dir = NULL) {
  core <- rad_core_info("native")
  .sirap_validate_rad_resource_core(core)
  .rrad_download_whitelists(
    kits = kits,
    force = force,
    quiet = quiet,
    cache_dir = cache_dir
  )
}
