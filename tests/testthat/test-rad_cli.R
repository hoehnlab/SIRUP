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
            kit = "10x_3v3", out_dir = "results", output = "demo",
            threads = 8L, bc_split = TRUE,
            rad = mock, capture = TRUE)
  got <- paste(readLines(log), collapse = " ")
  expect_match(got, "^demux")
  expect_match(got, "--kit 10x_3v3")
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

# ---- rad_src_path() resolution -------------------------------------------

make_fake_rad_src <- function() {
  d <- tempfile("rad_src_")
  dir.create(file.path(d, ".git"), recursive = TRUE)
  file.create(file.path(d, "CMakeLists.txt"))
  normalizePath(d)
}

test_that("rad_src_path() honours explicit path argument", {
  d <- make_fake_rad_src(); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  expect_equal(rad_src_path(d), d)
})

test_that("rad_src_path() rejects a non-clone path", {
  d <- tempfile("not_rad_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  expect_error(rad_src_path(d), "not a rad git clone")
})

test_that("rad_src_path() rejects a missing path", {
  expect_error(rad_src_path("/definitely/not/here"), "not found")
})

test_that("rad_src_path() finds the source by walking up from rad_path()", {
  skip_if_not_unix()
  src <- make_fake_rad_src(); on.exit(unlink(src, recursive = TRUE), add = TRUE)
  build_dir <- file.path(src, "build"); dir.create(build_dir)
  bin <- file.path(build_dir, "rad"); file.create(bin); Sys.chmod(bin, "0755")

  old_opt <- getOption("SIRAP.rad_path")
  on.exit(options(SIRAP.rad_path = old_opt), add = TRUE)
  options(SIRAP.rad_path = bin)

  expect_equal(rad_src_path(), src)
})

test_that("rad_src_path() resolves a symlinked binary", {
  skip_if_not_unix()
  src <- make_fake_rad_src(); on.exit(unlink(src, recursive = TRUE), add = TRUE)
  build_dir <- file.path(src, "build"); dir.create(build_dir)
  real_bin <- file.path(build_dir, "rad"); file.create(real_bin); Sys.chmod(real_bin, "0755")

  link_dir <- tempfile("bin_"); dir.create(link_dir); on.exit(unlink(link_dir, recursive = TRUE), add = TRUE)
  link_bin <- file.path(link_dir, "rad")
  file.symlink(real_bin, link_bin)

  old_opt <- getOption("SIRAP.rad_path")
  on.exit(options(SIRAP.rad_path = old_opt), add = TRUE)
  options(SIRAP.rad_path = link_bin)

  expect_equal(rad_src_path(), src)
})

# ---- rad_run end-to-end orchestrator -------------------------------------

test_that("rad_run requires out_dir", {
  expect_error(
    rad_run(layout = "five_prime", fastq = "x.fq.gz", rad = tempfile()),
    "out_dir"
  )
})

test_that("rad_run chains prep + demux by default (position_map = TRUE)", {
  skip_if_not_unix()
  log <- tempfile(); file.create(log); on.exit(unlink(log), add = TRUE)
  mock <- make_mock_rad(log); on.exit(unlink(mock), add = TRUE)
  out_dir <- tempfile(); dir.create(out_dir); on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  rad_run(
    layout = "five_prime", fastq = "reads.fq.gz",
    out_dir = out_dir, output = "demo",
    threads = 4L, rad = mock, capture = TRUE
  )

  got <- readLines(log)
  expect_length(got, 2)                       # prep + demux, no reformat
  expect_match(got[[1]], "^prep")
  expect_match(got[[1]], "--position-map")
  expect_match(got[[1]], "--fastq reads.fq.gz")
  expect_match(got[[1]], paste0("--output ", out_dir, "/demo_prep"))
  expect_match(got[[2]], "^demux")
  expect_match(got[[2]], paste0("--dir ", out_dir))
  expect_match(got[[2]], "--output demo")
})

test_that("rad_run with position_map = FALSE runs a lightweight prep", {
  skip_if_not_unix()
  log <- tempfile(); file.create(log); on.exit(unlink(log), add = TRUE)
  mock <- make_mock_rad(log); on.exit(unlink(mock), add = TRUE)
  out_dir <- tempfile(); dir.create(out_dir); on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  rad_run(
    layout = "five_prime", fastq = "reads.fq.gz",
    out_dir = out_dir, position_map = FALSE,
    rad = mock, capture = TRUE
  )

  got <- readLines(log)
  expect_length(got, 2)
  expect_match(got[[1]], "^prep")
  expect_match(got[[1]], "--read-layout")
  expect_false(grepl("--position-map", got[[1]]))
  expect_false(grepl("--fastq",        got[[1]]))
})

test_that("rad_run forwards `...` to demux", {
  skip_if_not_unix()
  log <- tempfile(); file.create(log); on.exit(unlink(log), add = TRUE)
  mock <- make_mock_rad(log); on.exit(unlink(mock), add = TRUE)
  out_dir <- tempfile(); dir.create(out_dir); on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  rad_run(
    layout = "five_prime", fastq = "reads.fq.gz", out_dir = out_dir,
    kit = "10x_3v3", bc_correction_mode = "offensive",
    rad = mock, capture = TRUE
  )

  got <- readLines(log)
  expect_match(got[[2]], "^demux")
  expect_match(got[[2]], "--kit 10x_3v3")
  expect_match(got[[2]], "--bc-correction-mode offensive")
})

test_that("rad_run runs all three steps when do_reformat = TRUE and demux output exists", {
  skip_if_not_unix()
  log <- tempfile(); file.create(log); on.exit(unlink(log), add = TRUE)
  mock <- make_mock_rad(log); on.exit(unlink(mock), add = TRUE)
  out_dir <- tempfile(); dir.create(out_dir); on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  # Touch the expected demux output so the reformat step proceeds.
  file.create(file.path(out_dir, "demo.fq.gz"))

  rad_run(
    layout = "five_prime", fastq = "reads.fq.gz",
    out_dir = out_dir, output = "demo",
    do_reformat      = TRUE,
    split_bc         = TRUE,
    reformat_out_dir = file.path(out_dir, "by_bc"),
    rad = mock, capture = TRUE
  )

  got <- readLines(log)
  expect_length(got, 3)
  expect_match(got[[1]], "^prep")
  expect_match(got[[2]], "^demux")
  expect_match(got[[3]], "^reformat")
  expect_match(got[[3]], "--split-bc")
  expect_match(got[[3]], paste0("--outdir ", out_dir, "/by_bc"))
  expect_match(got[[3]], paste0("--fastq ", out_dir, "/demo.fq.gz"))
})

test_that("rad_run reformat-specific flags propagate", {
  skip_if_not_unix()
  log <- tempfile(); file.create(log); on.exit(unlink(log), add = TRUE)
  mock <- make_mock_rad(log); on.exit(unlink(mock), add = TRUE)
  out_dir <- tempfile(); dir.create(out_dir); on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)
  file.create(file.path(out_dir, "demo.fq.gz"))

  rad_run(
    layout = "five_prime", fastq = "reads.fq.gz",
    out_dir = out_dir, output = "demo",
    do_reformat     = TRUE,
    reformat_header = TRUE,
    reformat_delim  = "_",
    coordinate      = "vizHD-v1",
    bin_size        = 8L,
    rad = mock, capture = TRUE
  )

  rfmt <- readLines(log)[[3]]
  expect_match(rfmt, "--reformat-header")
  expect_match(rfmt, "--reformat-delim _")
  expect_match(rfmt, "--coordinate=vizHD-v1")
  expect_match(rfmt, "--bin-size 8")
})

test_that("rad_run warns and skips reformat when demux output is missing", {
  skip_if_not_unix()
  log <- tempfile(); file.create(log); on.exit(unlink(log), add = TRUE)
  mock <- make_mock_rad(log); on.exit(unlink(mock), add = TRUE)
  out_dir <- tempfile(); dir.create(out_dir); on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  expect_warning(
    rad_run(
      layout = "five_prime", fastq = "reads.fq.gz",
      out_dir = out_dir,
      do_reformat      = TRUE,
      split_bc         = TRUE,
      reformat_out_dir = out_dir,
      rad = mock, capture = TRUE
    ),
    "skipping reformat"
  )

  got <- readLines(log)
  expect_length(got, 2)              # prep + demux only; reformat skipped
})

test_that("rad_run honours reformat_fastq override", {
  skip_if_not_unix()
  log <- tempfile(); file.create(log); on.exit(unlink(log), add = TRUE)
  mock <- make_mock_rad(log); on.exit(unlink(mock), add = TRUE)
  out_dir <- tempfile(); dir.create(out_dir); on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)
  custom_fq <- file.path(out_dir, "elsewhere.fq.gz")
  file.create(custom_fq)

  rad_run(
    layout = "five_prime", fastq = "reads.fq.gz",
    out_dir = out_dir,
    do_reformat    = TRUE,
    reformat_fastq = custom_fq,
    rad = mock, capture = TRUE
  )

  rfmt <- readLines(log)[[3]]
  expect_match(rfmt, paste0("--fastq ", custom_fq))
})
