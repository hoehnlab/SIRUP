.mock_whitelist_core <- function(api = 1L) {
  list(
    core_version = "1.1.0",
    core_api_version = api,
    source_commit = paste(rep("a", 40L), collapse = ""),
    embedded_source_digest = paste(rep("c", 64L), collapse = ""),
    resource_schema = 1L,
    resource_digest = paste(rep("b", 64L), collapse = ""),
    features = c(
      "demux", "scan-wl", "reformat", "remote-resource-cache"
    ),
    build = list(
      platform = "test",
      resources_available = TRUE,
      resource_dir = "/mock/rrad/resources"
    ),
    vendored = TRUE
  )
}

test_that("rad_download_whitelists delegates after the native handshake", {
  core <- .mock_whitelist_core()
  observed <- new.env(parent = emptyenv())
  expected <- stats::setNames(
    c("/cache/3v3.txt.gz", "/cache/5v3.txt.gz"),
    c("10x_3v3", "10x_5v3")
  )
  testthat::local_mocked_bindings(
    .rrad_core_info = function() core,
    .rrad_download_whitelists = function(...) {
      observed$args <- list(...)
      expected
    },
    .package = "SIRAP"
  )

  result <- rad_download_whitelists(
    kits = c("10x_3v3", "10x_5v3"),
    force = TRUE,
    quiet = TRUE,
    cache_dir = "/cache"
  )

  expect_identical(result, expected)
  expect_identical(observed$args, list(
    kits = c("10x_3v3", "10x_5v3"),
    force = TRUE,
    quiet = TRUE,
    cache_dir = "/cache"
  ))
})

test_that("rad_download_whitelists refuses an incompatible native core", {
  called <- FALSE
  testthat::local_mocked_bindings(
    .rrad_core_info = function() .mock_whitelist_core(api = 2L),
    .rrad_download_whitelists = function(...) {
      called <<- TRUE
      character()
    },
    .package = "SIRAP"
  )

  expect_error(
    rad_download_whitelists("10x_3v3"),
    class = "sirap_rad_incompatible_error"
  )
  expect_false(called)
})

test_that("rad_download_whitelists requires the native resource cache", {
  core <- .mock_whitelist_core()
  core$features <- setdiff(core$features, "remote-resource-cache")
  called <- FALSE
  testthat::local_mocked_bindings(
    .rrad_core_info = function() core,
    .rrad_download_whitelists = function(...) {
      called <<- TRUE
      character()
    },
    .package = "SIRAP"
  )

  expect_error(
    rad_download_whitelists("10x_3v3"),
    "remote-resource-cache",
    class = "sirap_rad_incompatible_error"
  )
  expect_false(called)
})

test_that("rad_download_whitelists has rrad argument parity", {
  expect_identical(
    formals(rad_download_whitelists),
    formals(rrad::rad_download_whitelists)
  )
})
