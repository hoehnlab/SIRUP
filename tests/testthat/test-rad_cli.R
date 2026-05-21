skip_if_not_unix <- function() {
  if (.Platform$OS.type != "unix") testthat::skip("POSIX shell mock not portable.")
}

make_mock_rad <- function(log_path, exit_status = 0L) {
  bin <- tempfile(fileext = ".sh")
  lines <- c(
    "#!/bin/sh",
    sprintf("printf '%%s\\n' \"$*\" >> \"%s\"", log_path),
    sprintf("exit %d", exit_status)
  )
  writeLines(lines, bin)
  Sys.chmod(bin, "0755")
  bin
}

# ---- rad_path() resolution -----------------------------------------------

test_that("rad_path() honours explicit path argument", {
  bin <- tempfile()
  file.create(bin)
  on.exit(unlink(bin), add = TRUE)
  expect_equal(rad_path(bin), normalizePath(bin))
})

test_that("rad_path() errors when explicit path is missing", {
  expect_error(rad_path("/definitely/not/here/rad"), "not found")
})

test_that("rad_path() falls back to option, then env, then Sys.which", {
  skip_if_not_unix()
  bin_opt <- tempfile(); file.create(bin_opt)
  bin_env <- tempfile(); file.create(bin_env)
  on.exit({ unlink(c(bin_opt, bin_env)) }, add = TRUE)

  old_opt <- getOption("SIRAP.rad_path")
  old_env <- Sys.getenv("SIRAP_RAD_PATH", unset = NA)
  on.exit({
    options(SIRAP.rad_path = old_opt)
    if (is.na(old_env)) Sys.unsetenv("SIRAP_RAD_PATH") else Sys.setenv(SIRAP_RAD_PATH = old_env)
  }, add = TRUE)

  options(SIRAP.rad_path = bin_opt)
  expect_equal(rad_path(), normalizePath(bin_opt))

  options(SIRAP.rad_path = NULL)
  Sys.setenv(SIRAP_RAD_PATH = bin_env)
  expect_equal(rad_path(), normalizePath(bin_env))
})

# ---- rad_prep argument assembly ------------------------------------------

test_that("rad_prep requires at least one mode", {
  expect_error(
    rad_prep(layout = "five_prime", read_layout = FALSE, position_map = FALSE,
             rad = tempfile()),
    "read_layout"
  )
})

test_that("rad_prep validates position_map prerequisites", {
  expect_error(
    rad_prep(layout = "five_prime", read_layout = FALSE, position_map = TRUE,
             rad = tempfile()),
    "requires"
  )
})

test_that("rad_prep assembles the expected argv", {
  skip_if_not_unix()
  log <- tempfile(); file.create(log); on.exit(unlink(log), add = TRUE)
  mock <- make_mock_rad(log); on.exit(unlink(mock), add = TRUE)

  rad_prep(layout = "five_prime", threads = 4L, verbose = TRUE,
           rad = mock, capture = TRUE)
  got <- paste(readLines(log), collapse = " ")
  expect_match(got, "prep")
  expect_match(got, "--layout five_prime")
  expect_match(got, "--read-layout")
  expect_match(got, "--threads 4")
  expect_match(got, "--verbose")
})

test_that("rad_prep position_map mode includes fastq + output", {
  skip_if_not_unix()
  log <- tempfile(); file.create(log); on.exit(unlink(log), add = TRUE)
  mock <- make_mock_rad(log); on.exit(unlink(mock), add = TRUE)

  rad_prep(layout = "five_prime", read_layout = FALSE, position_map = TRUE,
           fastq = "reads.fq.gz", output = "out/run1",
           rad = mock, capture = TRUE)
  got <- paste(readLines(log), collapse = " ")
  expect_match(got, "--position-map")
  expect_match(got, "--fastq reads.fq.gz")
  expect_match(got, "--output out/run1")
  expect_false(grepl("--read-layout", got))
})

# ---- rad_demux ------------------------------------------------------------

test_that("rad_demux maps named args to the right CLI flags", {
  skip_if_not_unix()
  log <- tempfile(); file.create(log); on.exit(unlink(log), add = TRUE)
  mock <- make_mock_rad(log); on.exit(unlink(mock), add = TRUE)

  rad_demux(layout = "five_prime", fastq = "reads.fq.gz",
            kit = "10xv3", out_dir = "results", output = "demo",
            threads = 8L, bc_split = TRUE,
            rad = mock, capture = TRUE)
  got <- paste(readLines(log), collapse = " ")
  expect_match(got, "^demux")
  expect_match(got, "--kit 10xv3")
  expect_match(got, "--dir results")
  expect_match(got, "--output demo")
  expect_match(got, "--bc-split")
  expect_match(got, "--threads 8")
})

# ---- rad_reformat ---------------------------------------------------------

test_that("rad_reformat split_bc requires out_dir", {
  expect_error(
    rad_reformat(fastq = "x.fq.gz", split_bc = TRUE, rad = tempfile()),
    "split_bc"
  )
})

test_that("rad_reformat coordinate accepts logical and string", {
  skip_if_not_unix()
  log <- tempfile(); file.create(log); on.exit(unlink(log), add = TRUE)
  mock <- make_mock_rad(log); on.exit(unlink(mock), add = TRUE)

  rad_reformat(fastq = "x.fq.gz", coordinate = TRUE,
               rad = mock, capture = TRUE)
  rad_reformat(fastq = "x.fq.gz", coordinate = "vizHD-v1", bin_size = 8L,
               rad = mock, capture = TRUE)

  got <- readLines(log)
  expect_match(got[[1]], " --coordinate( |$)")
  expect_match(got[[2]], "--coordinate=vizHD-v1")
  expect_match(got[[2]], "--bin-size 8")
})

# ---- rad_modify validation ------------------------------------------------

test_that("rad_modify enforces exactly one target and one action", {
  expect_error(
    rad_modify(rad = tempfile()),
    "exactly one of `layout_key`"
  )
  expect_error(
    rad_modify(layout_key = "x", whitelist_kit = "y", rad = tempfile()),
    "exactly one of `layout_key`"
  )
  expect_error(
    rad_modify(layout_key = "x", rad = tempfile()),
    "exactly one of `add ="
  )
  expect_error(
    rad_modify(layout_key = "x", add = "p", remove = TRUE, rad = tempfile()),
    "exactly one of `add ="
  )
})

# ---- rad_scan_wl validation -----------------------------------------------

test_that("rad_scan_wl requires either input or batch_csv", {
  expect_error(
    rad_scan_wl(adapter_seq = "ACGT", rad = tempfile()),
    "input"
  )
})

# ---- non-zero exit propagates --------------------------------------------

test_that(".rad_call surfaces non-zero exit codes by default", {
  skip_if_not_unix()
  log <- tempfile(); file.create(log); on.exit(unlink(log), add = TRUE)
  mock <- make_mock_rad(log, exit_status = 2L); on.exit(unlink(mock), add = TRUE)

  expect_error(
    rad_prep(layout = "five_prime", rad = mock, capture = TRUE),
    "status 2"
  )
})
