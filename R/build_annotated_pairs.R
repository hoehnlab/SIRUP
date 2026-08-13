#' Annotate RAD whitelist pairs with Visium HD coordinates and Space Ranger flags
#'
#' Takes RAD's whitelist output (valid barcode pairs with counts) and
#' annotates each pair with:
#'   - spatial grid coordinates at 2 um / 8 um / 16 um resolutions, looked
#'     up from rad's Visium HD coordinate reference
#'     (`rad/resources/wl/visium_hd_coordinates.csv.gz`)
#'   - Space Ranger style barcode strings
#'     (`s_<res>_<row>_<col>-1`, e.g. `s_002um_01468_00608-1`)
#'   - `in_tissue` flag from Space Ranger `tissue_positions.parquet`
#'   - `in_sr_filtered_{2,8,16}um` flags from each Space Ranger
#'     `filtered_feature_bc_matrix/barcodes.tsv.gz`
#'   - 8 um and 16 um bin-level aggregated counts and spot counts
#'   - `rad_valid_16um` (1 if the 16 um bin contains >= 32 of 64
#'     possible 2 um sub-bins with any RAD signal)
#'   - Normalized counts: `ncpm`, `log1p_ncpm` at the 2 / 8 / 16 um levels.
#'
#' R port of `rad_paper/visium_hd_mapping/buildAnnotatedPairs.py`.
#'
#' @param valid_pairs Path to RAD's whitelist output CSV (or a
#'   `data.table`). Must contain columns `bc1`, `bc2`, `pair`, `count`.
#'   Any pre-computed columns from a previous run of this function are
#'   dropped before re-deriving them.
#' @param rad_coords Path to the Visium HD coordinate reference
#'   (`visium_hd_coordinates.csv.gz` shipped in `rad/resources/wl/`) or a
#'   `data.table` with columns `x_bc`, `y_bc`, `002um_x`, `002um_y`,
#'   `008um_x`, `008um_y`, `016um_x`, `016um_y`.
#' @param spaceranger_dir Directory containing the Space Ranger output
#'   tree: it should have `square_002um/`, `square_008um/`,
#'   `square_016um/` subdirectories. If supplied, the four
#'   `tissue_positions` / `sr_filtered_*um` paths below are inferred from
#'   it unless explicitly overridden.
#' @param tissue_positions Path to Space Ranger
#'   `tissue_positions.parquet` (or `.csv` / `.csv.gz`). Reading
#'   `.parquet` requires the `arrow` package. If `NULL` and
#'   `spaceranger_dir` is set, inferred as
#'   `<spaceranger_dir>/square_002um/spatial/tissue_positions.parquet`
#'   (falling back to `.csv` / `.csv.gz` if the parquet is missing).
#' @param sr_filtered_2um Path to Space Ranger filtered 2um barcodes
#'   (`barcodes.tsv.gz`). If `NULL` and `spaceranger_dir` is set, inferred from
#'   that root. A path that doesn't exist is treated as an empty filtered set.
#' @param sr_filtered_8um Path to Space Ranger filtered 8um barcodes
#'   (`barcodes.tsv.gz`). If `NULL` and `spaceranger_dir` is set, inferred from
#'   that root. A path that doesn't exist is treated as an empty filtered set.
#' @param sr_filtered_16um Path to Space Ranger filtered 16um barcodes
#'   (`barcodes.tsv.gz`). If `NULL` and `spaceranger_dir` is set, inferred from
#'   that root. A path that doesn't exist is treated as an empty filtered set.
#' @param out_csv Output CSV path. Gzipped when the path ends in `.gz`.
#'   `NULL` (default) only returns the result.
#' @param verbose Print progress and a summary block.
#' @return Invisibly, a `data.table` with columns
#'   `bc1, bc2, pair, count, ncpm, log1p_ncpm,
#'   x_2um, y_2um, spaceranger_2um, in_tissue, in_sr_filtered_2um,
#'   x_8um, y_8um, spaceranger_8um, in_sr_filtered_8um,
#'   count_8um, n_spots_8um, ncpm_8um, log1p_ncpm_8um,
#'   x_16um, y_16um, spaceranger_16um, in_sr_filtered_16um,
#'   count_16um, n_spots_16um, ncpm_16um, log1p_ncpm_16um,
#'   rad_valid_16um`.
#' @export
buildAnnotatedPairs <- function(valid_pairs,
                                  rad_coords,
                                  spaceranger_dir = NULL,
                                  tissue_positions = NULL,
                                  sr_filtered_2um = NULL,
                                  sr_filtered_8um = NULL,
                                  sr_filtered_16um = NULL,
                                  out_csv = NULL,
                                  verbose = TRUE) {
  log_msg <- function(...) if (isTRUE(verbose)) message(...)

  # --- Resolve paths from spaceranger_dir when individual ones aren't set --
  if (!is.null(spaceranger_dir)) {
    if (is.null(tissue_positions)) {
      candidates <- c(
        file.path(spaceranger_dir, "square_002um", "spatial", "tissue_positions.parquet"),
        file.path(spaceranger_dir, "square_002um", "spatial", "tissue_positions.csv"),
        file.path(spaceranger_dir, "square_002um", "spatial", "tissue_positions.csv.gz")
      )
      hits <- candidates[file.exists(candidates)]
      if (length(hits) > 0L) tissue_positions <- hits[[1]]
    }
    if (is.null(sr_filtered_2um))  sr_filtered_2um  <- file.path(spaceranger_dir, "square_002um", "filtered_feature_bc_matrix", "barcodes.tsv.gz")
    if (is.null(sr_filtered_8um))  sr_filtered_8um  <- file.path(spaceranger_dir, "square_008um", "filtered_feature_bc_matrix", "barcodes.tsv.gz")
    if (is.null(sr_filtered_16um)) sr_filtered_16um <- file.path(spaceranger_dir, "square_016um", "filtered_feature_bc_matrix", "barcodes.tsv.gz")
  }

  # --- Step 1: load RAD valid pairs ---------------------------------------
  log_msg("Loading RAD whitelist...")
  rad <- if (data.table::is.data.table(valid_pairs)) {
    data.table::copy(valid_pairs)
  } else {
    data.table::fread(valid_pairs)
  }
  for (col in c("bc1", "bc2", "count")) {
    if (!col %in% names(rad)) {
      stop("`valid_pairs` is missing required column: ", col)
    }
  }
  log_msg(sprintf("  %s barcode pairs", format(nrow(rad), big.mark = ",")))

  # Drop columns from a previous run so we re-derive them downstream.
  drop_cols <- intersect(
    c("n_2um_spots_16","rad_valid_16um","n_2um_spots","n_spots_8um","n_spots_16um"),
    names(rad)
  )
  if (length(drop_cols) > 0L) rad[, (drop_cols) := NULL]

  # --- Step 2: load coordinate reference + build lookups ------------------
  log_msg("Loading coordinate reference...")
  coords <- if (data.table::is.data.table(rad_coords)) {
    data.table::copy(rad_coords)
  } else {
    data.table::fread(rad_coords)
  }
  required_coord_cols <- c("x_bc", "y_bc",
                           "002um_x", "002um_y",
                           "008um_x", "008um_y",
                           "016um_x", "016um_y")
  missing_cc <- setdiff(required_coord_cols, names(coords))
  if (length(missing_cc) > 0L) {
    stop("`rad_coords` is missing columns: ", paste(missing_cc, collapse = ", "))
  }

  # Rename digit-prefixed columns up front so we don't need backticks
  # everywhere. We mirror the Python script's final names.
  data.table::setnames(
    coords,
    old = c("002um_x","002um_y","008um_x","008um_y","016um_x","016um_y"),
    new = c("x_2um", "y_2um", "x_8um", "y_8um", "x_16um", "y_16um")
  )

  x_lookup <- unique(coords[, .(x_bc, x_2um, x_8um, x_16um)], by = "x_bc")
  y_lookup <- unique(coords[, .(y_bc, y_2um, y_8um, y_16um)], by = "y_bc")

  log_msg("Mapping bc1 -> X indices...")
  rad <- merge(rad, x_lookup, by.x = "bc1", by.y = "x_bc", all.x = TRUE, sort = FALSE)
  log_msg("Mapping bc2 -> Y indices...")
  rad <- merge(rad, y_lookup, by.x = "bc2", by.y = "y_bc", all.x = TRUE, sort = FALSE)

  mapped <- sum(!is.na(rad$x_2um))
  if (mapped < nrow(rad)) {
    log_msg(sprintf("  Mapped: %s / %s (%.2f%%)",
                    format(mapped,     big.mark = ","),
                    format(nrow(rad),  big.mark = ","),
                    100 * mapped / nrow(rad)))
  } else {
    log_msg(sprintf("  Mapped: %s / %s (100%%)",
                    format(mapped,    big.mark = ","),
                    format(nrow(rad), big.mark = ",")))
  }

  # --- Step 3: build Space Ranger barcode strings -------------------------
  # Visium HD convention: s_<res>_<array_row>_<array_col>-1
  # In the rad coord ref, <res>_x is the array column and <res>_y is the
  # array row, so y is first in the string.
  pad5 <- function(v) sprintf("%05d", as.integer(v))
  rad[, spaceranger_2um  := paste0("s_002um_", pad5(y_2um),  "_", pad5(x_2um),  "-1")]
  rad[, spaceranger_8um  := paste0("s_008um_", pad5(y_8um),  "_", pad5(x_8um),  "-1")]
  rad[, spaceranger_16um := paste0("s_016um_", pad5(y_16um), "_", pad5(x_16um), "-1")]

  # --- Step 4: annotate with tissue + filtered flags ----------------------
  load_barcode_set <- function(path, label) {
    if (is.null(path) || !file.exists(path)) {
      if (!is.null(path)) {
        log_msg(sprintf("  SKIP %s: %s missing", label, path))
      }
      return(character(0))
    }
    dt <- data.table::fread(path, header = FALSE)
    dt[[1]]
  }

  load_tissue_positions_dt <- function(path) {
    if (is.null(path) || !file.exists(path)) return(NULL)
    ext <- tools::file_ext(path)
    if (ext == "parquet") {
      if (!requireNamespace("arrow", quietly = TRUE)) {
        stop(
          "Reading tissue_positions.parquet requires the 'arrow' package. ",
          "Install it (install.packages('arrow')) or supply a .csv / .csv.gz path."
        )
      }
      return(data.table::as.data.table(arrow::read_parquet(path)))
    }
    data.table::fread(path)
  }

  log_msg("Loading tissue_positions...")
  tp <- load_tissue_positions_dt(tissue_positions)
  tissue_set <- if (!is.null(tp) && all(c("in_tissue", "barcode") %in% names(tp))) {
    tp[in_tissue == 1L, barcode]
  } else if (!is.null(tp)) {
    log_msg("  SKIP tissue: missing in_tissue/barcode columns")
    character(0)
  } else {
    character(0)
  }
  log_msg(sprintf("  %s spots in tissue", format(length(tissue_set), big.mark = ",")))

  log_msg("Loading Space Ranger filtered barcodes...")
  filt_2um  <- load_barcode_set(sr_filtered_2um,  "2um")
  filt_8um  <- load_barcode_set(sr_filtered_8um,  "8um")
  filt_16um <- load_barcode_set(sr_filtered_16um, "16um")
  log_msg(sprintf("  Filtered: %s (2um), %s (8um), %s (16um)",
                  format(length(filt_2um),  big.mark = ","),
                  format(length(filt_8um),  big.mark = ","),
                  format(length(filt_16um), big.mark = ",")))

  rad[, in_tissue          := as.integer(spaceranger_2um  %in% tissue_set)]
  rad[, in_sr_filtered_2um := as.integer(spaceranger_2um  %in% filt_2um)]
  rad[, in_sr_filtered_8um := as.integer(spaceranger_8um  %in% filt_8um)]
  rad[, in_sr_filtered_16um := as.integer(spaceranger_16um %in% filt_16um)]

  # --- Step 5: aggregate bin-level counts ---------------------------------
  log_msg("Aggregating 8um bins...")
  agg_8um  <- rad[, .(count_8um  = sum(count), n_spots_8um  = .N),
                  by = spaceranger_8um]

  log_msg("Aggregating 16um bins...")
  agg_16um <- rad[, .(count_16um = sum(count), n_spots_16um = .N),
                  by = spaceranger_16um]

  rad <- merge(rad, agg_8um,  by = "spaceranger_8um",  all.x = TRUE, sort = FALSE)
  rad <- merge(rad, agg_16um, by = "spaceranger_16um", all.x = TRUE, sort = FALSE)

  # rad_valid_16um: at least half of the 64 possible 2 um sub-bins lit up.
  rad[, rad_valid_16um := as.integer(n_spots_16um >= 32L)]

  # --- Step 6: normalized counts ------------------------------------------
  total <- sum(rad$count)
  scale_factor <- if (total > 0) 1e6 / total else 0
  rad[, ncpm            := count      * scale_factor]
  rad[, log1p_ncpm      := log1p(ncpm)]
  rad[, ncpm_8um        := count_8um  * scale_factor]
  rad[, log1p_ncpm_8um  := log1p(ncpm_8um)]
  rad[, ncpm_16um       := count_16um * scale_factor]
  rad[, log1p_ncpm_16um := log1p(ncpm_16um)]

  # --- Step 7: final column order -----------------------------------------
  out_cols <- c(
    "bc1","bc2","pair","count","ncpm","log1p_ncpm",
    "x_2um","y_2um","spaceranger_2um","in_tissue","in_sr_filtered_2um",
    "x_8um","y_8um","spaceranger_8um","in_sr_filtered_8um",
    "count_8um","n_spots_8um","ncpm_8um","log1p_ncpm_8um",
    "x_16um","y_16um","spaceranger_16um","in_sr_filtered_16um",
    "count_16um","n_spots_16um","ncpm_16um","log1p_ncpm_16um","rad_valid_16um"
  )
  # If `pair` was missing from the input, allow that.
  out_cols <- intersect(out_cols, names(rad))
  rad <- rad[, out_cols, with = FALSE]

  # --- Summary ------------------------------------------------------------
  if (isTRUE(verbose)) {
    pct <- function(x) sprintf("%s (%.2f%%)",
                               format(sum(x), big.mark = ","), 100 * mean(x))
    message("\n=== SUMMARY ===")
    message(sprintf("Total RAD barcode pairs:       %s", format(nrow(rad), big.mark = ",")))
    message(sprintf("In tissue (2um):               %s", pct(rad$in_tissue)))
    message(sprintf("In SR filtered (2um):          %s", pct(rad$in_sr_filtered_2um)))
    message(sprintf("In SR filtered (8um):          %s", pct(rad$in_sr_filtered_8um)))
    message(sprintf("In SR filtered (16um):         %s", pct(rad$in_sr_filtered_16um)))
    message(sprintf("Unique 8um bins:               %s",
                    format(data.table::uniqueN(rad$spaceranger_8um),  big.mark = ",")))
    message(sprintf("Unique 16um bins:              %s",
                    format(data.table::uniqueN(rad$spaceranger_16um), big.mark = ",")))
    n_valid_bins <- data.table::uniqueN(rad[rad_valid_16um == 1L, spaceranger_16um])
    message(sprintf("RAD valid 16um bins (>=32):    %s / %s",
                    format(n_valid_bins, big.mark = ","),
                    format(data.table::uniqueN(rad$spaceranger_16um), big.mark = ",")))
    message(sprintf("  2um pairs in valid bins:     %s", pct(rad$rad_valid_16um)))
  }

  # --- Write --------------------------------------------------------------
  if (!is.null(out_csv)) {
    data.table::fwrite(rad, out_csv)   # auto-gzips when path ends in .gz
    log_msg(sprintf("\nWritten to: %s", out_csv))
  }

  invisible(rad[])
}
