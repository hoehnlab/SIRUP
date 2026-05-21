# End-to-end tests for build_annotated_pairs() using a tiny synthetic
# Visium HD fixture. The fixture has 2 distinct bc1 values and 2 bc2
# values, giving 4 unique RAD pairs. Coord ref maps them onto a small
# region; the 2um grid is downsampled 4x to 8um and 2x further to 16um
# (the real Visium HD ratios).

make_coords <- function() {
  # Construct directly as a data.table; column names with digits need
  # backticks at construction time. We rely on build_annotated_pairs to
  # rename them internally (matching the rad/resources/wl/ schema).
  data.table::data.table(
    x_bc    = c("BC_A_seq", "BC_B_seq"),
    y_bc    = c("BC_X_seq", "BC_Y_seq"),
    `002um_x` = c(100L, 200L),
    `002um_y` = c(50L, 150L),
    `008um_x` = c(25L,  50L),
    `008um_y` = c(12L,  37L),
    `016um_x` = c(12L,  25L),
    `016um_y` = c(6L,   18L)
  )
}

make_valid_pairs <- function() {
  data.table::data.table(
    bc1   = c("BC_A_seq", "BC_A_seq", "BC_B_seq", "BC_B_seq"),
    bc2   = c("BC_X_seq", "BC_Y_seq", "BC_X_seq", "BC_Y_seq"),
    pair  = c("A_X", "A_Y", "B_X", "B_Y"),
    count = c(10L, 5L, 20L, 15L)
  )
}

# ---- structural / annotation tests ---------------------------------------

test_that("build_annotated_pairs returns one row per pair with the documented columns", {
  out <- build_annotated_pairs(make_valid_pairs(), make_coords(), verbose = FALSE)

  expect_s3_class(out, "data.table")
  expect_equal(nrow(out), 4L)

  expect_identical(
    names(out),
    c("bc1","bc2","pair","count","ncpm","log1p_ncpm",
      "x_2um","y_2um","spaceranger_2um","in_tissue","in_sr_filtered_2um",
      "x_8um","y_8um","spaceranger_8um","in_sr_filtered_8um",
      "count_8um","n_spots_8um","ncpm_8um","log1p_ncpm_8um",
      "x_16um","y_16um","spaceranger_16um","in_sr_filtered_16um",
      "count_16um","n_spots_16um","ncpm_16um","log1p_ncpm_16um","rad_valid_16um")
  )
})

test_that("coordinate lookup maps bc1 -> x columns and bc2 -> y columns", {
  out <- build_annotated_pairs(make_valid_pairs(), make_coords(), verbose = FALSE)

  # bc1 == BC_A_seq -> x_2um = 100; bc1 == BC_B_seq -> x_2um = 200
  expect_equal(unique(out[bc1 == "BC_A_seq", x_2um]), 100L)
  expect_equal(unique(out[bc1 == "BC_B_seq", x_2um]), 200L)
  expect_equal(unique(out[bc1 == "BC_A_seq", x_8um]), 25L)
  expect_equal(unique(out[bc1 == "BC_B_seq", x_16um]), 25L)

  # bc2 == BC_X_seq -> y_2um = 50; bc2 == BC_Y_seq -> y_2um = 150
  expect_equal(unique(out[bc2 == "BC_X_seq", y_2um]), 50L)
  expect_equal(unique(out[bc2 == "BC_Y_seq", y_2um]), 150L)
})

test_that("spaceranger barcode strings follow s_<res>_<row>_<col>-1 with 5-digit padding", {
  out <- build_annotated_pairs(make_valid_pairs(), make_coords(), verbose = FALSE)

  # (A, X): 2um coords (100, 50) — y is row, x is col -> s_002um_00050_00100-1
  expect_equal(out[bc1 == "BC_A_seq" & bc2 == "BC_X_seq", spaceranger_2um],
               "s_002um_00050_00100-1")
  expect_equal(out[bc1 == "BC_A_seq" & bc2 == "BC_X_seq", spaceranger_8um],
               "s_008um_00012_00025-1")
  expect_equal(out[bc1 == "BC_A_seq" & bc2 == "BC_X_seq", spaceranger_16um],
               "s_016um_00006_00012-1")

  # (B, Y): 16um coords (25, 18) -> s_016um_00018_00025-1
  expect_equal(out[bc1 == "BC_B_seq" & bc2 == "BC_Y_seq", spaceranger_16um],
               "s_016um_00018_00025-1")
})

test_that("missing tissue + filtered inputs leave all flags at 0", {
  out <- build_annotated_pairs(make_valid_pairs(), make_coords(), verbose = FALSE)
  expect_true(all(out$in_tissue          == 0L))
  expect_true(all(out$in_sr_filtered_2um == 0L))
  expect_true(all(out$in_sr_filtered_8um == 0L))
  expect_true(all(out$in_sr_filtered_16um == 0L))
})

# ---- bin aggregation -----------------------------------------------------

test_that("8 um and 16 um aggregations sum counts and count spots correctly", {
  # Build a fixture where two pairs land in the same 8 um bin so we can
  # verify the sum. (BC_A, BC_X) and (BC_A, BC_X) again — but our base
  # fixture has unique (bc1, bc2) pairs, so duplicate at the bc level by
  # using two bc2 values that map to the same 8um row.
  coords <- data.table::data.table(
    x_bc    = c("X1", "X2"),
    y_bc    = c("Y1", "Y2"),
    `002um_x` = c(100L, 200L),
    `002um_y` = c(50L, 150L),
    # Both 2um rows map to the SAME 8um and 16um row (Y1, Y2 -> y_8um=12)
    `008um_x` = c(25L,  50L),
    `008um_y` = c(12L,  12L),
    `016um_x` = c(12L,  25L),
    `016um_y` = c( 6L,   6L)
  )
  pairs <- data.table::data.table(
    bc1 = c("X1","X1","X2","X2"),
    bc2 = c("Y1","Y2","Y1","Y2"),
    pair = c("a","b","c","d"),
    count = c(10L, 5L, 20L, 15L)
  )
  out <- build_annotated_pairs(pairs, coords, verbose = FALSE)

  # Each 8um bin holds two 2um pairs that share the same (x_8um, y_8um).
  # Pair (X1,Y1) -> 8um (25,12); (X1,Y2) -> 8um (25,12) too. So both share
  # spaceranger_8um == s_008um_00012_00025-1; their counts sum to 10+5=15.
  bin_sumX1 <- out[bc1 == "X1", unique(count_8um)]
  expect_equal(bin_sumX1, 10L + 5L)
  bin_sumX2 <- out[bc1 == "X2", unique(count_8um)]
  expect_equal(bin_sumX2, 20L + 15L)

  # n_spots_8um is the count of 2um pairs in that 8um bin, so 2 each
  expect_equal(unique(out[bc1 == "X1", n_spots_8um]), 2L)
  expect_equal(unique(out[bc1 == "X2", n_spots_8um]), 2L)
})

test_that("rad_valid_16um is 1 iff n_spots_16um >= 32", {
  out <- build_annotated_pairs(make_valid_pairs(), make_coords(), verbose = FALSE)
  # Only 4 rows total -> n_spots_16um is at most 4 in any bin -> all 0
  expect_true(all(out$rad_valid_16um == 0L))

  # Build a fixture with 32+ pairs landing in the same 16um bin
  coords <- data.table::data.table(
    x_bc    = paste0("X", 1:8),
    y_bc    = paste0("Y", 1:8),
    `002um_x` = seq.int(8) * 10L,
    `002um_y` = seq.int(8) * 10L,
    `008um_x` = rep(1L, 8),    # all in 8um bin (1, 1)
    `008um_y` = rep(1L, 8),
    `016um_x` = rep(1L, 8),    # all in 16um bin (1, 1) too
    `016um_y` = rep(1L, 8)
  )
  pairs <- data.table::CJ(bc1 = paste0("X", 1:8), bc2 = paste0("Y", 1:8))
  pairs[, pair  := paste0(bc1, "_", bc2)]
  pairs[, count := 1L]
  out2 <- build_annotated_pairs(pairs, coords, verbose = FALSE)
  # 8 * 8 = 64 pairs all in the same 16um bin -> all rows valid
  expect_equal(unique(out2$n_spots_16um), 64L)
  expect_true(all(out2$rad_valid_16um == 1L))
})

# ---- normalization -------------------------------------------------------

test_that("ncpm = count * 1e6 / sum(count); log1p_ncpm = log1p(ncpm)", {
  pairs <- make_valid_pairs()
  out <- build_annotated_pairs(pairs, make_coords(), verbose = FALSE)

  total <- sum(pairs$count)  # 50
  expect_equal(out$ncpm, pairs$count * 1e6 / total)
  expect_equal(out$log1p_ncpm, log1p(out$ncpm))

  # ncpm_8um and ncpm_16um use the SAME scale factor (per the Python port —
  # not a re-normalization at each level)
  expect_equal(out$ncpm_8um,  out$count_8um  * 1e6 / total)
  expect_equal(out$ncpm_16um, out$count_16um * 1e6 / total)
})

# ---- previous-run column drop -------------------------------------------

test_that("pre-computed columns from a previous run are dropped before re-derive", {
  pairs <- make_valid_pairs()
  pairs[, rad_valid_16um := 99L]      # bogus stale column
  pairs[, n_spots_16um   := 999L]

  out <- build_annotated_pairs(pairs, make_coords(), verbose = FALSE)
  # The re-derived rad_valid_16um is 0 (only 4 rows), so the bogus 99 is gone.
  expect_true(all(out$rad_valid_16um == 0L))
  expect_false(any(out$n_spots_16um == 999L))
})

# ---- Space Ranger flag annotation ---------------------------------------

test_that("in_tissue + in_sr_filtered_* flags fire when the corresponding files exist", {
  tmp_dir <- tempfile("sr_"); dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  # tissue_positions.csv: column names per Space Ranger
  tp_path <- file.path(tmp_dir, "tissue_positions.csv")
  data.table::fwrite(data.table::data.table(
    barcode = c("s_002um_00050_00100-1",     # (BC_A, BC_X)
                "s_002um_00150_00200-1"),    # (BC_B, BC_Y)
    in_tissue = c(1L, 1L)
  ), tp_path)

  # filtered barcodes at each resolution (gzipped .tsv, one per line)
  filt_2 <- file.path(tmp_dir, "filtered_2um.tsv.gz")
  filt_8 <- file.path(tmp_dir, "filtered_8um.tsv.gz")
  filt_16 <- file.path(tmp_dir, "filtered_16um.tsv.gz")
  writeLines("s_002um_00050_00100-1", gzfile(filt_2));  close.connection <- NULL
  writeLines("s_008um_00012_00025-1", gzfile(filt_8))
  writeLines("s_016um_00018_00025-1", gzfile(filt_16))

  out <- build_annotated_pairs(
    make_valid_pairs(), make_coords(),
    tissue_positions  = tp_path,
    sr_filtered_2um   = filt_2,
    sr_filtered_8um   = filt_8,
    sr_filtered_16um  = filt_16,
    verbose = FALSE
  )

  # (A, X) and (B, Y) are in-tissue; (A, Y) and (B, X) are not.
  expect_equal(out[bc1 == "BC_A_seq" & bc2 == "BC_X_seq", in_tissue], 1L)
  expect_equal(out[bc1 == "BC_B_seq" & bc2 == "BC_Y_seq", in_tissue], 1L)
  expect_equal(out[bc1 == "BC_A_seq" & bc2 == "BC_Y_seq", in_tissue], 0L)
  expect_equal(out[bc1 == "BC_B_seq" & bc2 == "BC_X_seq", in_tissue], 0L)

  # 2um filtered: only (A, X) is in the filtered set
  expect_equal(out[bc1 == "BC_A_seq" & bc2 == "BC_X_seq", in_sr_filtered_2um], 1L)
  expect_equal(sum(out$in_sr_filtered_2um), 1L)

  # 8um filtered: in make_coords() each bc2 maps to a distinct y_8um, so
  # all four 2um pairs land in DIFFERENT 8um bins. Only (A, X)'s bin
  # s_008um_00012_00025-1 was added to the filtered set.
  expect_equal(out[bc1 == "BC_A_seq" & bc2 == "BC_X_seq", in_sr_filtered_8um], 1L)
  expect_equal(sum(out$in_sr_filtered_8um), 1L)

  # 16um filtered: (B, Y) maps to s_016um_00018_00025-1
  expect_equal(out[bc1 == "BC_B_seq" & bc2 == "BC_Y_seq", in_sr_filtered_16um], 1L)
  expect_equal(sum(out$in_sr_filtered_16um), 1L)
})

# ---- file I/O round-trip -------------------------------------------------

test_that("build_annotated_pairs writes a gzipped CSV that round-trips identically", {
  out_path <- tempfile(fileext = ".csv.gz")
  on.exit(unlink(out_path), add = TRUE)

  result <- build_annotated_pairs(
    make_valid_pairs(), make_coords(),
    out_csv = out_path, verbose = FALSE
  )
  expect_true(file.exists(out_path))

  reread <- data.table::fread(out_path)
  expect_equal(nrow(reread), nrow(result))
  expect_equal(reread$count,         result$count)
  expect_equal(reread$spaceranger_2um, result$spaceranger_2um)
  expect_equal(reread$ncpm_16um,     result$ncpm_16um)
})

# ---- input file vs in-memory parity --------------------------------------

test_that("file-path and data.table inputs produce identical output", {
  pairs_dt  <- make_valid_pairs()
  coords_dt <- make_coords()

  pairs_path  <- tempfile(fileext = ".csv"); on.exit(unlink(pairs_path),  add = TRUE)
  coords_path <- tempfile(fileext = ".csv"); on.exit(unlink(coords_path), add = TRUE)
  data.table::fwrite(pairs_dt,  pairs_path)
  data.table::fwrite(coords_dt, coords_path)

  out_dt   <- build_annotated_pairs(pairs_dt,   coords_dt,   verbose = FALSE)
  out_file <- build_annotated_pairs(pairs_path, coords_path, verbose = FALSE)

  expect_equal(out_dt, out_file)
})

# ---- validation ----------------------------------------------------------

test_that("missing required columns surface a clear error", {
  bad_pairs <- data.table::data.table(bc1 = "a", count = 1L)   # no bc2
  expect_error(
    build_annotated_pairs(bad_pairs, make_coords(), verbose = FALSE),
    "bc2"
  )
  bad_coords <- make_coords()[, .(x_bc, y_bc)]                 # no x/y_2um etc.
  expect_error(
    build_annotated_pairs(make_valid_pairs(), bad_coords, verbose = FALSE),
    "rad_coords"
  )
})
