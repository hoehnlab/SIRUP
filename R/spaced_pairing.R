#' Build a SPaCeD-style clone spatial table
#'
#' SPaCeD pairs heavy and light chain clones over a weighted combination of
#' spatial evidence and a REPAIR model matrix. These helpers implement that
#' pairing surface while using sparse native spatial-bin operations that scale
#' better for Visium HD-sized data.
#'
#' The method follows the SPaCeD repository by H. J. Melton et al.
#' (https://github.com/hjm2217/spaced) and its use of a REPAIR model matrix
#' (https://github.com/almaan/star-repair; modified SPaCeD notes reference
#' https://github.com/Floraliu7/star-repair). The original SPaCeD spatial term
#' was computed with `spatstat::pppdist(type = "mat")`; SIRUP's default spatial
#' term uses sparse clone-by-barcode matrices. For Visium HD, use the Space
#' Ranger 8um barcodes as the native spatial bins rather than adding another
#' coordinate downsampling step.
#'
#' @param data Clone-assigned data frame.
#' @param positions Optional spatial positions data frame. If supplied, it is
#'   joined to `data` by barcode.
#' @param barcode Barcode column in `data`.
#' @param position_barcode Barcode column in `positions`.
#' @param clone Clone identifier column in `data`.
#' @param x Spatial x-coordinate column. If `positions` is supplied, this is
#'   read after the join.
#' @param y Spatial y-coordinate column. If `positions` is supplied, this is
#'   read after the join.
#' @param count Optional count column. When `NULL`, each row contributes 1.
#' @param clone_prefix Optional prefix added to clone IDs, for example
#'   `"heavy_"` or `"light_kappa_"`.
#' @param drop_missing_coords Drop rows missing spatial coordinates.
#' @return Data frame with `clone_id`, `barcode`, `x`, `y`, and `n`.
#' @export
buildSpacedCloneTable <- function(data,
                                  positions = NULL,
                                  barcode = "barcode",
                                  clone = "clone_id",
                                  x = "x",
                                  y = "y",
                                  count = NULL,
                                  position_barcode = "barcode",
                                  clone_prefix = NULL,
                                  drop_missing_coords = TRUE) {
  .check_required_cols(data, c(barcode, clone), "data")
  data <- .as_data_frame(data)

  if (!is.null(positions)) {
    .check_required_cols(positions, c(position_barcode, x, y), "positions")
    data <- data[, setdiff(colnames(data), c(x, y)), drop = FALSE]
    pos <- .as_data_frame(positions)
    join_cols <- stats::setNames(position_barcode, barcode)
    data <- dplyr::left_join(data, pos, by = join_cols)
  }

  .check_required_cols(data, c(barcode, clone, x, y), "data")
  if (!is.null(count)) {
    .check_required_cols(data, count, "data")
  }

  clone_id <- as.character(data[[clone]])
  if (!is.null(clone_prefix)) {
    clone_id <- paste0(clone_prefix, clone_id)
  }

  counts <- if (is.null(count)) {
    rep(1, nrow(data))
  } else {
    as.numeric(data[[count]])
  }
  counts[is.na(counts)] <- 0

  out <- data.frame(
    clone_id = clone_id,
    barcode = as.character(data[[barcode]]),
    x = as.numeric(data[[x]]),
    y = as.numeric(data[[y]]),
    n = counts,
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$clone_id) & out$clone_id != "" & !is.na(out$barcode), , drop = FALSE]
  if (drop_missing_coords) {
    out <- out[!is.na(out$x) & !is.na(out$y), , drop = FALSE]
  }
  out <- out[out$n > 0, , drop = FALSE]

  out <- stats::aggregate(
    n ~ clone_id + barcode + x + y,
    data = out,
    FUN = sum
  )
  out[order(out$clone_id, out$barcode), c("clone_id", "barcode", "x", "y", "n")]
}

#' Convert a wide SPaCeD expression table to clone spatial points
#'
#' @param df Wide table with barcode, coordinate, and clone count columns.
#' @param barcode Barcode column name.
#' @param x Spatial x-coordinate column name.
#' @param y Spatial y-coordinate column name.
#' @param clone Optional clone columns. Defaults to all non-barcode and
#'   non-coordinate columns.
#' @param clone_prefix Optional prefix added to clone IDs.
#' @return Data frame with `clone_id`, `barcode`, `x`, `y`, and `n`.
#' @export
spacedCloneTableFromWide <- function(df,
                                     barcode = colnames(df)[[1]],
                                     x = "x",
                                     y = "y",
                                     clone = NULL,
                                     clone_prefix = NULL) {
  df <- .as_data_frame(df)
  .check_required_cols(df, c(barcode, x, y), "df")
  if (is.null(clone)) {
    clone <- setdiff(colnames(df), c(barcode, x, y))
  }
  .check_required_cols(df, clone, "df")

  pieces <- lapply(clone, function(clone_id_col) {
    counts <- as.numeric(df[[clone_id_col]])
    keep <- !is.na(counts) & counts > 0
    if (!any(keep)) {
      return(NULL)
    }
    clone_id <- if (is.null(clone_prefix)) clone_id_col else paste0(clone_prefix, clone_id_col)
    data.frame(
      clone_id = clone_id,
      barcode = as.character(df[[barcode]][keep]),
      x = as.numeric(df[[x]][keep]),
      y = as.numeric(df[[y]][keep]),
      n = counts[keep],
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, pieces[!vapply(pieces, is.null, logical(1))])
  if (is.null(out)) {
    out <- data.frame(clone_id = character(), barcode = character(), x = numeric(), y = numeric(), n = numeric())
  }
  rownames(out) <- NULL
  out
}

#' Compute scalable SPaCeD spatial affinity
#'
#' Builds clone-by-spatial-feature matrices and compares heavy/light clones by a
#' spatial method. The default `method = "weighted_cosim"` is a scalable
#' distance-aware overlap score: nearby 8um barcodes receive partial credit
#' through a sparse Gaussian spatial kernel. `method = "ot"` uses the original
#' SPaCeD-style point-pattern optimal-transport distance through
#' `spatstat::pppdist()`.
#'
#' @param heavy_points Heavy-chain output from `buildSpacedCloneTable()` or any
#'   data frame with `clone_id`, `barcode`, and `n`.
#' @param light_points Light-chain output from `buildSpacedCloneTable()` or any
#'   data frame with `clone_id`, `barcode`, and `n`.
#' @param method Spatial scoring method. `"weighted_cosim"` is the scalable
#'   SIRUP default, and `"ot"` is the original SPaCeD optimal-transport
#'   point-pattern distance.
#' @param feature Column containing the native spatial feature/bin ID.
#' @param kernel_sigma Gaussian kernel scale for `method = "weighted_cosim"`,
#'   in the same units as `x` and `y`. When `NULL`, SIRUP estimates the
#'   median nearest-neighbor feature spacing and uses twice that value.
#' @param kernel_radius Maximum feature-feature distance retained in the sparse
#'   spatial kernel. Defaults to `3 * kernel_sigma`.
#' @param kernel_max_neighbors Maximum nearby features retained per feature when
#'   building the sparse kernel.
#' @param x Spatial x-coordinate column.
#' @param y Spatial y-coordinate column.
#' @param ot_cutoff Maximum transport distance passed to `pppdist()` when
#'   `method = "ot"`.
#' @param ot_scale Distance scale used by `transformSpacedDistances()` for
#'   `method = "ot"`.
#' @param nproc Number of workers used by the slow `method = "ot"` clone-pair
#'   distance loop. The sparse methods ignore this argument.
#' @param top_n Optional number of top light candidates retained per heavy
#'   clone. Non-retained affinities are set to zero.
#' @param min_affinity Optional lower threshold for retained affinities.
#' @param min_counts Optional minimum total clone count. Mirrors REPAIR's
#'   `--min_counts` filter.
#' @param min_x_obs Minimum number of spatial features where a clone has count
#'   greater than or equal to `filter_x`. Mirrors REPAIR's `--min_x_obs`
#'   filter.
#' @param filter_x Count threshold used with `min_x_obs`.
#' @return Heavy-by-light numeric affinity matrix with rows normalized to sum to
#'   one when a row has nonzero affinity.
#' @export
spacedSpatialAffinity <- function(heavy_points,
                                  light_points,
                                  method = c("weighted_cosim", "ot"),
                                  feature = "barcode",
                                  kernel_sigma = NULL,
                                  kernel_radius = NULL,
                                  kernel_max_neighbors = 500L,
                                  x = "x",
                                  y = "y",
                                  ot_cutoff = 1e6,
                                  ot_scale = 1,
                                  nproc = 1L,
                                  top_n = NULL,
                                  min_affinity = 0,
                                  min_counts = NULL,
                                  min_x_obs = 2,
                                  filter_x = 1) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    rlang::abort("Package 'Matrix' is required for `spacedSpatialAffinity()`.")
  }
  method <- match.arg(method)
  if (!is.character(feature) || length(feature) != 1L || is.na(feature) || feature == "") {
    rlang::abort("`feature` must be a single non-empty column name.")
  }
  heavy_points <- .sirup_spaced_filter_clone_points(
    heavy_points,
    name = "heavy_points",
    feature_col = feature,
    min_counts = min_counts,
    min_x_obs = min_x_obs,
    filter_x = filter_x
  )
  light_points <- .sirup_spaced_filter_clone_points(
    light_points,
    name = "light_points",
    feature_col = feature,
    min_counts = min_counts,
    min_x_obs = min_x_obs,
    filter_x = filter_x
  )

  sim <- switch(
    method,
    weighted_cosim = .sirup_spaced_weighted_cosim_affinity(
      heavy_points = heavy_points,
      light_points = light_points,
      feature_col = feature,
      kernel_sigma = kernel_sigma,
      kernel_radius = kernel_radius,
      kernel_max_neighbors = kernel_max_neighbors,
      x_col = x,
      y_col = y
    ),
    ot = .sirup_spaced_ot_affinity(
      heavy_points = heavy_points,
      light_points = light_points,
      x_col = x,
      y_col = y,
      cutoff = ot_cutoff,
      scale = ot_scale,
      nproc = nproc
    )
  )

  if (!is.null(top_n)) {
    top_n <- as.integer(top_n)
    if (is.na(top_n) || top_n < 1L) {
      rlang::abort("`top_n` must be a positive integer when supplied.")
    }
    sim <- .sirup_spaced_prune_top_n(sim, top_n)
  }
  if (!is.null(min_affinity) && min_affinity > 0) {
    sim[sim < min_affinity] <- 0
  }

  .sirup_spaced_row_normalize(sim)
}

#' Transform SPaCeD spatial distances into affinities
#'
#' @param dist_matrix Heavy-by-light distance matrix.
#' @param scale Optional divisor applied before exponentiation. Use a larger
#'   value when distances are in pixel units.
#' @return Row-normalized heavy-by-light affinity matrix.
#' @export
transformSpacedDistances <- function(dist_matrix, scale = 1) {
  dist_matrix <- as.matrix(dist_matrix)
  if (!is.numeric(dist_matrix)) {
    storage.mode(dist_matrix) <- "numeric"
  }
  if (!is.numeric(scale) || length(scale) != 1L || is.na(scale) || scale <= 0) {
    rlang::abort("`scale` must be a positive number.")
  }
  affinity <- exp(-dist_matrix / scale)
  affinity[!is.finite(affinity)] <- 0
  .sirup_spaced_row_normalize(affinity)
}

#' Pair heavy and light clones with SPaCeD's weighted scores
#'
#' Returns the top light-chain candidates for each heavy clone. Light clones are
#' allowed to appear as candidates for multiple heavy clones, which is useful
#' when public light chains are plausible.
#'
#' @param spatial_affinity Heavy-by-light spatial affinity matrix, data frame,
#'   or path to a delimited matrix file.
#' @param mapping_matrix Optional heavy-by-light REPAIR model matrix, data
#'   frame, or path to a delimited matrix file. If `NULL`, only spatial affinity
#'   is used. Unnamed matrices are aligned positionally to `spatial_affinity`,
#'   matching the current SPaCeD behavior.
#' @param weights Numeric vector of omega values. `1` is spatial-only; `0` is
#'   mapping-only.
#' @param top_n_light_chains Number of top light-chain candidates to report for
#'   each heavy clone.
#' @return If one weight is supplied, a data frame with up to
#'   `top_n_light_chains` candidate light clones per heavy clone. If multiple
#'   weights are supplied, a named list of those data frames. Each data frame
#'   includes spatial score (cosine similarity/optimal transport score based on position, cite SPaCeD paper),
#'   REPAIR score (score from REPAIR (cite Engblom '23)), combined score, row-normalized pairing score,
#'   and light-chain specificity score (normalized pairing score for this pair across all pairs for this chain).
#' @export
pairHeavyLightSpaced <- function(spatial_affinity,
                                 mapping_matrix = NULL,
                                 weights = seq(0, 1, by = 0.1),
                                 top_n_light_chains = 5L) {
  spatial_affinity <- .sirup_spaced_matrix_arg(spatial_affinity, "spatial_affinity")
  if (is.null(rownames(spatial_affinity)) || is.null(colnames(spatial_affinity))) {
    rlang::abort("`spatial_affinity` must have rownames for heavy clones and colnames for light clones.")
  }
  if (anyDuplicated(rownames(spatial_affinity)) || anyDuplicated(colnames(spatial_affinity))) {
    rlang::abort("`spatial_affinity` rownames and colnames must be unique clone IDs.")
  }
  mapping_matrix <- .sirup_spaced_prepare_mapping(mapping_matrix, spatial_affinity)
  weights <- as.numeric(weights)
  if (any(is.na(weights) | weights < 0 | weights > 1)) {
    rlang::abort("`weights` must contain values between 0 and 1.")
  }
  top_n_light_chains <- as.integer(top_n_light_chains)
  if (length(top_n_light_chains) != 1L || is.na(top_n_light_chains) || top_n_light_chains < 1L) {
    rlang::abort("`top_n_light_chains` must be a positive integer.")
  }
  top_n_light_chains <- min(top_n_light_chains, ncol(spatial_affinity))

  out <- lapply(weights, function(w) {
    score <- (w * spatial_affinity) + ((1 - w) * mapping_matrix)
    normalized_pairing <- .sirup_spaced_row_normalize(score)
    light_specificity <- .sirup_spaced_col_normalize(score)
    .sirup_spaced_top_light_candidates(
      spatial_affinity = spatial_affinity,
      mapping_matrix = mapping_matrix,
      score = score,
      normalized_pairing = normalized_pairing,
      light_specificity = light_specificity,
      weight = w,
      top_n_light_chains = top_n_light_chains
    )
  })
  names(out) <- paste0("w=", format(weights, trim = TRUE, scientific = FALSE))
  if (length(out) == 1L) {
    return(out[[1L]])
  }
  out
}

#' Build REPAIR input for SPaCeD
#'
#' REPAIR expects a spatial observation-by-clone expression table with heavy and
#' light clone columns distinguishable by regex patterns such as `"^heavy"` and
#' `"^light"`. This helper builds that matrix directly from SPaCeD/SIRUP clone
#' point tables.
#'
#' @param heavy_points Heavy-chain clone spatial table with `clone_id`,
#'   `barcode`, and `n`.
#' @param light_points Light-chain clone spatial table with `clone_id`,
#'   `barcode`, and `n`.
#' @param sparse Return a sparse `Matrix` instead of a dense data frame.
#' @return Barcode-by-clone count table. Dense output includes a `barcode`
#'   column; sparse output has barcodes as row names and clones as column names.
#' @export
buildRepairInputMatrix <- function(heavy_points,
                                   light_points,
                                   sparse = FALSE) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    rlang::abort("Package 'Matrix' is required for `buildRepairInputMatrix()`.")
  }
  heavy_points <- .sirup_spaced_check_repair_points(heavy_points, "heavy_points")
  light_points <- .sirup_spaced_check_repair_points(light_points, "light_points")
  points <- rbind(heavy_points, light_points)

  barcodes <- sort(unique(points$barcode))
  clones <- sort(unique(points$clone_id))
  mat <- Matrix::sparseMatrix(
    i = match(points$barcode, barcodes),
    j = match(points$clone_id, clones),
    x = points$n,
    dims = c(length(barcodes), length(clones)),
    dimnames = list(barcodes, clones)
  )
  mat <- Matrix::drop0(mat)
  if (sparse) {
    return(mat)
  }
  out <- as.data.frame(as.matrix(mat), check.names = FALSE)
  out <- cbind(barcode = rownames(out), out)
  rownames(out) <- NULL
  out
}

#' Write REPAIR input matrix
#'
#' @param repair_input Output from `buildRepairInputMatrix()`.
#' @param path Output TSV path.
#' @param include_barcode Include row/barcode IDs as the first column.
#' @return Invisibly, `path`.
#' @export
writeRepairInputMatrix <- function(repair_input,
                                   path,
                                   include_barcode = TRUE) {
  out_dir <- dirname(path)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (inherits(repair_input, "Matrix")) {
    df <- as.data.frame(as.matrix(repair_input), check.names = FALSE)
    if (include_barcode) {
      df <- cbind(barcode = rownames(df), df)
    }
  } else {
    df <- as.data.frame(repair_input, check.names = FALSE)
    if (!include_barcode && "barcode" %in% colnames(df)) {
      rownames(df) <- df$barcode
      df$barcode <- NULL
    }
  }
  utils::write.table(
    df,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
  invisible(path)
}

#' Construct a REPAIR analyze command for SPaCeD
#'
#' @param input_path REPAIR input TSV path.
#' @param output_dir REPAIR output directory.
#' @param a_pattern Heavy clone column regex.
#' @param b_pattern Light clone column regex.
#' @param mxo Maximum number of output candidates passed to `repair analyze`.
#' @param x REPAIR `x` argument passed to `repair analyze`.
#' @param b REPAIR `b` argument passed to `repair analyze`.
#' @param nb REPAIR `nb` argument passed to `repair analyze`.
#' @param min_counts Minimum count threshold passed to `repair analyze`.
#' @param repair_bin REPAIR executable.
#' @return Shell command string.
#' @export
repairAnalyzeCommand <- function(input_path,
                                 output_dir,
                                 a_pattern = "^heavy",
                                 b_pattern = "^light",
                                 mxo = 2,
                                 x = 1,
                                 b = 100,
                                 nb = 1,
                                 min_counts = NULL,
                                 repair_bin = "repair") {
  command <- c(
    shQuote(repair_bin),
    "analyze",
    "-i", shQuote(input_path),
    "-o", shQuote(output_dir),
    "-ap", shQuote(a_pattern),
    "-bp", shQuote(b_pattern),
    "-mxo", shQuote(as.character(mxo)),
    "-x", shQuote(as.character(x)),
    "-b", shQuote(as.character(b)),
    "-nb", shQuote(as.character(nb))
  )
  if (!is.null(min_counts)) {
    command <- c(command, "-mc", shQuote(as.character(min_counts)))
  }
  paste(command, collapse = " ")
}

#' Check for the SPaCeD-compatible REPAIR fork
#'
#' SPaCeD needs a modified REPAIR installation that writes a
#' `*_model_matrix.tsv` file. This helper checks the Python REPAIR source in a
#' conda environment and fails early if the installed package does not appear to
#' be the modified fork used by SPaCeD.
#'
#' @param conda_env Conda environment name. Set to `NULL` or `""` to run
#'   `python_bin` directly without `conda run`.
#' @param conda_exe Conda executable.
#' @param python_bin Python executable to run inside the conda environment.
#' @param verbose Print the detected REPAIR path.
#' @return Invisibly, a list with the detected REPAIR analyze path.
#' @export
checkSpacedRepair <- function(conda_env = Sys.getenv("REPAIR_CONDA_ENV", "repair"),
                              conda_exe = Sys.getenv("CONDA_EXE", "conda"),
                              python_bin = "python",
                              verbose = TRUE) {
  python_code <- paste(
    "import inspect",
    "import pathlib",
    "import repair.analyze as analyze",
    "path = pathlib.Path(inspect.getfile(analyze))",
    "text = path.read_text()",
    "print('SIRUP_REPAIR_ANALYZE_PATH=' + str(path))",
    "print('SIRUP_REPAIR_EXPORTS_MODEL_MATRIX=' + str('_model_matrix.tsv' in text))",
    sep = "; "
  )

  if (is.null(conda_env) || length(conda_env) != 1L || is.na(conda_env) || conda_env == "") {
    command <- python_bin
    args <- c("-c", python_code)
  } else {
    command <- conda_exe
    args <- c("run", "-n", conda_env, python_bin, "-c", python_code)
  }

  output <- tryCatch(
    system2(command, shQuote(args), stdout = TRUE, stderr = TRUE),
    error = function(e) e
  )
  if (inherits(output, "error")) {
    rlang::abort(paste0(
      "Could not run REPAIR check using `", command, "`. ",
      "Install the SPaCeD-compatible REPAIR fork or check `conda_exe`, ",
      "`conda_env`, and `python_bin`. Original error: ", conditionMessage(output)
    ))
  }

  status <- attr(output, "status")
  if (is.null(status)) {
    status <- 0L
  }
  if (!identical(as.integer(status), 0L)) {
    rlang::abort(paste0(
      "Could not import REPAIR from the requested environment. ",
      "Install the SPaCeD-compatible REPAIR fork in `", conda_env, "`. Output:\n",
      paste(output, collapse = "\n")
    ))
  }

  path_line <- grep("^SIRUP_REPAIR_ANALYZE_PATH=", output, value = TRUE)
  matrix_line <- grep("^SIRUP_REPAIR_EXPORTS_MODEL_MATRIX=", output, value = TRUE)
  if (length(path_line) != 1L || length(matrix_line) != 1L) {
    rlang::abort(paste0(
      "REPAIR check did not return the expected markers. Output:\n",
      paste(output, collapse = "\n")
    ))
  }

  analyze_path <- sub("^SIRUP_REPAIR_ANALYZE_PATH=", "", path_line)
  exports_model_matrix <- identical(sub("^SIRUP_REPAIR_EXPORTS_MODEL_MATRIX=", "", matrix_line), "True")
  if (!exports_model_matrix) {
    rlang::abort(paste0(
      "The installed REPAIR package does not appear to export `*_model_matrix.tsv`. ",
      "SIRUP's SPaCeD pairing workflow needs the modified REPAIR fork from ",
      "https://github.com/Floraliu7/star-repair. Detected REPAIR analyze file: ",
      analyze_path
    ))
  }

  if (isTRUE(verbose)) {
    message("Found SPaCeD-compatible REPAIR at ", analyze_path)
  }
  invisible(list(analyze_path = analyze_path, exports_model_matrix = TRUE))
}

.sirup_spaced_check_repair_points <- function(points, name) {
  points <- .as_data_frame(points)
  .check_required_cols(points, c("clone_id", "barcode", "n"), name)
  points$clone_id <- as.character(points$clone_id)
  points$barcode <- as.character(points$barcode)
  points$n <- as.numeric(points$n)
  points <- points[
    !is.na(points$clone_id) &
      !is.na(points$barcode) &
      !is.na(points$n) &
      points$n > 0,
    ,
    drop = FALSE
  ]
  if (nrow(points) == 0L) {
    rlang::abort(paste0("`", name, "` has no usable clone counts."))
  }
  stats::aggregate(n ~ clone_id + barcode, data = points, FUN = sum)
}

.sirup_spaced_matrix_arg <- function(x, name) {
  if (is.character(x) && length(x) == 1L) {
    if (!file.exists(x)) {
      rlang::abort(paste0("`", name, "` path does not exist: ", x))
    }
    x <- .sirup_spaced_read_matrix_file(x)
  } else if (is.data.frame(x)) {
    x <- .sirup_spaced_data_frame_to_matrix(x)
  } else {
    x <- as.matrix(x)
  }
  if (!is.numeric(x)) {
    storage.mode(x) <- "numeric"
  }
  if (any(!is.finite(x))) {
    x[!is.finite(x)] <- 0
  }
  x
}

.sirup_spaced_read_matrix_file <- function(path) {
  first_line <- readLines(path, n = 1L, warn = FALSE)
  if (length(first_line) == 0L) {
    rlang::abort(paste0("Matrix file is empty: ", path))
  }
  fields <- strsplit(first_line, "\t", fixed = TRUE)[[1]]
  fields_numeric <- suppressWarnings(!is.na(as.numeric(fields)))
  header <- !all(fields != "" & fields_numeric)
  df <- utils::read.delim(path, header = header, check.names = FALSE)
  .sirup_spaced_data_frame_to_matrix(df)
}

.sirup_spaced_data_frame_to_matrix <- function(df) {
  df <- as.data.frame(df, check.names = FALSE)
  if (ncol(df) > 1L) {
    first_col_numeric <- suppressWarnings(!any(is.na(as.numeric(df[[1L]]))))
    if (!first_col_numeric) {
      row_names <- as.character(df[[1L]])
      df <- df[-1L]
      mat <- as.matrix(df)
      rownames(mat) <- row_names
      return(mat)
    }
  }
  as.matrix(df)
}

.sirup_spaced_has_matching_dimnames <- function(mat, target) {
  !is.null(rownames(mat)) &&
    !is.null(colnames(mat)) &&
    setequal(rownames(target), rownames(mat)) &&
    setequal(colnames(target), colnames(mat))
}

.sirup_spaced_filter_clone_points <- function(points,
                                              name,
                                              feature_col,
                                              min_counts,
                                              min_x_obs,
                                              filter_x) {
  points <- .as_data_frame(points)
  .check_required_cols(points, c("clone_id", feature_col, "n"), name)
  points$clone_id <- as.character(points$clone_id)
  points$feature <- as.character(points[[feature_col]])
  points$n <- as.numeric(points$n)
  points <- points[
    !is.na(points$clone_id) &
      points$clone_id != "" &
      !is.na(points$feature) &
      points$feature != "" &
      !is.na(points$n) &
      points$n > 0,
    ,
    drop = FALSE
  ]
  if (nrow(points) == 0L) {
    rlang::abort(paste0("`", name, "` has no usable clone counts in `", feature_col, "`."))
  }
  if (!is.null(min_counts)) {
    min_counts <- as.numeric(min_counts)
    if (length(min_counts) != 1L || is.na(min_counts) || min_counts < 0) {
      rlang::abort("`min_counts` must be NULL or a non-negative number.")
    }
  }
  if (!is.null(min_x_obs)) {
    min_x_obs <- as.numeric(min_x_obs)
    if (length(min_x_obs) != 1L || is.na(min_x_obs) || min_x_obs < 0) {
      rlang::abort("`min_x_obs` must be NULL or a non-negative number.")
    }
    filter_x <- as.numeric(filter_x)
    if (length(filter_x) != 1L || is.na(filter_x) || filter_x < 0) {
      rlang::abort("`filter_x` must be a non-negative number.")
    }
  }

  counts <- stats::aggregate(n ~ clone_id + feature, data = points, FUN = sum)
  total_counts <- stats::aggregate(n ~ clone_id, data = counts, FUN = sum)
  keep <- total_counts$clone_id
  if (!is.null(min_counts)) {
    keep <- intersect(keep, total_counts$clone_id[total_counts$n > min_counts])
  }
  if (!is.null(min_x_obs)) {
    observed <- stats::aggregate(
      keep_feature ~ clone_id,
      data = data.frame(
        clone_id = counts$clone_id,
        keep_feature = counts$n >= filter_x,
        stringsAsFactors = FALSE
      ),
      FUN = sum
    )
    keep <- intersect(keep, observed$clone_id[observed$keep_feature >= min_x_obs])
  }
  points <- points[points$clone_id %in% keep, , drop = FALSE]
  if (nrow(points) == 0L) {
    rlang::abort(paste0(
      "`", name, "` has no clones after applying `min_counts`/`min_x_obs` filters."
    ))
  }
  points$feature <- NULL
  points
}

.sirup_spaced_native_features <- function(points, name, feature_col) {
  points <- .as_data_frame(points)
  .check_required_cols(points, c("clone_id", feature_col, "n"), name)
  points$clone_id <- as.character(points$clone_id)
  points$feature <- as.character(points[[feature_col]])
  points$n <- as.numeric(points$n)
  points <- points[
    !is.na(points$clone_id) &
      points$clone_id != "" &
      !is.na(points$feature) &
      points$feature != "" &
      !is.na(points$n) &
      points$n > 0,
    ,
    drop = FALSE
  ]
  if (nrow(points) == 0L) {
    rlang::abort(paste0("`", name, "` has no usable clone counts in `", feature_col, "`."))
  }
  stats::aggregate(n ~ clone_id + feature, data = points, FUN = sum) |>
    stats::setNames(c("clone_id", "feature", "value"))
}

.sirup_spaced_spatial_features <- function(points, name, feature_col, x_col, y_col) {
  points <- .as_data_frame(points)
  .check_required_cols(points, c("clone_id", feature_col, "n", x_col, y_col), name)
  points$clone_id <- as.character(points$clone_id)
  points$feature <- as.character(points[[feature_col]])
  points$n <- as.numeric(points$n)
  points$x <- as.numeric(points[[x_col]])
  points$y <- as.numeric(points[[y_col]])
  points <- points[
    !is.na(points$clone_id) &
      points$clone_id != "" &
      !is.na(points$feature) &
      points$feature != "" &
      !is.na(points$n) &
      points$n > 0 &
      is.finite(points$x) &
      is.finite(points$y),
    ,
    drop = FALSE
  ]
  if (nrow(points) == 0L) {
    rlang::abort(paste0("`", name, "` has no usable clone counts with finite spatial coordinates."))
  }
  stats::aggregate(
    n ~ clone_id + feature + x + y,
    data = points,
    FUN = sum
  ) |>
    stats::setNames(c("clone_id", "feature", "x", "y", "value"))
}

.sirup_spaced_sparse_clone_matrix <- function(features, feature_levels) {
  clone_levels <- sort(unique(features$clone_id))
  i <- match(features$clone_id, clone_levels)
  j <- match(features$feature, feature_levels)
  mat <- Matrix::sparseMatrix(
    i = i,
    j = j,
    x = features$value,
    dims = c(length(clone_levels), length(feature_levels)),
    dimnames = list(clone_levels, feature_levels)
  )
  mat
}

.sirup_spaced_weighted_cosim_affinity <- function(heavy_points,
                                                    light_points,
                                                    feature_col,
                                                    kernel_sigma,
                                                    kernel_radius,
                                                    kernel_max_neighbors,
                                                    x_col,
                                                    y_col) {
  heavy_features <- .sirup_spaced_spatial_features(heavy_points, "heavy_points", feature_col, x_col, y_col)
  light_features <- .sirup_spaced_spatial_features(light_points, "light_points", feature_col, x_col, y_col)
  feature_levels <- union(heavy_features$feature, light_features$feature)
  coords <- .sirup_spaced_feature_coordinates(rbind(heavy_features, light_features), feature_levels)

  if (is.null(kernel_sigma)) {
    kernel_sigma <- 2 * .sirup_spaced_median_nn_distance(coords)
  }
  kernel_sigma <- as.numeric(kernel_sigma)
  if (!is.finite(kernel_sigma) || kernel_sigma <= 0) {
    rlang::abort("`kernel_sigma` must be a positive finite number.")
  }
  if (is.null(kernel_radius)) {
    kernel_radius <- 3 * kernel_sigma
  }
  kernel_radius <- as.numeric(kernel_radius)
  if (!is.finite(kernel_radius) || kernel_radius <= 0) {
    rlang::abort("`kernel_radius` must be a positive finite number.")
  }

  heavy_mat <- .sirup_spaced_sparse_clone_matrix(heavy_features, feature_levels)
  light_mat <- .sirup_spaced_sparse_clone_matrix(light_features, feature_levels)
  kernel <- .sirup_spaced_spatial_kernel(
    coords = coords,
    sigma = kernel_sigma,
    radius = kernel_radius,
    max_neighbors = kernel_max_neighbors
  )
  .sirup_spaced_kernel_similarity(heavy_mat, light_mat, kernel)
}

.sirup_spaced_kernel_similarity <- function(heavy_mat, light_mat, kernel) {
  heavy_kernel <- heavy_mat %*% kernel
  light_kernel <- light_mat %*% kernel
  sim <- as.matrix(heavy_kernel %*% Matrix::t(light_mat))
  row_norm <- sqrt(Matrix::rowSums(heavy_kernel * heavy_mat))
  col_norm <- sqrt(Matrix::rowSums(light_kernel * light_mat))
  denom <- as.numeric(row_norm) %o% as.numeric(col_norm)
  sim[denom > 0] <- sim[denom > 0] / denom[denom > 0]
  sim[denom <= 0] <- 0
  sim[!is.finite(sim)] <- 0
  sim
}

.sirup_spaced_feature_coordinates <- function(features, feature_levels) {
  coords <- unique(features[, c("feature", "x", "y")])
  duplicated_features <- coords$feature[duplicated(coords$feature)]
  if (length(duplicated_features) > 0L) {
    bad <- unique(duplicated_features)
    rlang::abort(paste0(
      "Spatial feature coordinates are not unique for: ",
      paste(utils::head(bad, 5), collapse = ", "),
      if (length(bad) > 5L) ", ..." else "",
      ". Use one coordinate per `feature_col` value."
    ))
  }
  coords <- coords[match(feature_levels, coords$feature), , drop = FALSE]
  if (anyNA(coords$x) || anyNA(coords$y)) {
    rlang::abort("Could not find finite coordinates for every spatial feature.")
  }
  rownames(coords) <- coords$feature
  coords[, c("x", "y"), drop = FALSE]
}

.sirup_spaced_median_nn_distance <- function(coords) {
  if (nrow(coords) < 2L) {
    rlang::abort("At least two spatial features are required to estimate `kernel_sigma`; supply `kernel_sigma` directly.")
  }
  if (!requireNamespace("RANN", quietly = TRUE)) {
    rlang::abort("Package 'RANN' is required to estimate `kernel_sigma`.")
  }
  nn <- RANN::nn2(as.matrix(coords), k = 2L)
  d <- nn$nn.dists[, 2L]
  d <- d[is.finite(d) & d > 0]
  if (length(d) == 0L) {
    rlang::abort("Could not estimate spatial feature spacing; supply `kernel_sigma` directly.")
  }
  stats::median(d)
}

.sirup_spaced_spatial_kernel <- function(coords, sigma, radius, max_neighbors) {
  if (!requireNamespace("RANN", quietly = TRUE)) {
    rlang::abort("Package 'RANN' is required for `method = \"weighted_cosim\"`.")
  }
  max_neighbors <- as.integer(max_neighbors)
  if (is.na(max_neighbors) || max_neighbors < 1L) {
    rlang::abort("`kernel_max_neighbors` must be a positive integer.")
  }
  coords_mat <- as.matrix(coords)
  nn <- RANN::nn2(
    data = coords_mat,
    query = coords_mat,
    k = min(max_neighbors, nrow(coords_mat)),
    searchtype = "radius",
    radius = radius
  )
  i <- rep(seq_len(nrow(coords_mat)), ncol(nn$nn.idx))
  j <- as.vector(nn$nn.idx)
  d <- as.vector(nn$nn.dists)
  keep <- j > 0L & is.finite(d) & d <= radius
  i <- i[keep]
  j <- j[keep]
  d <- d[keep]
  weights <- exp(-(d ^ 2) / (2 * sigma ^ 2))
  Matrix::sparseMatrix(
    i = i,
    j = j,
    x = weights,
    dims = c(nrow(coords_mat), nrow(coords_mat)),
    dimnames = list(rownames(coords), rownames(coords))
  )
}

.sirup_spaced_ot_affinity <- function(heavy_points,
                                      light_points,
                                      x_col,
                                      y_col,
                                      cutoff,
                                      scale,
                                      nproc) {
  pppdist <- .sirup_spaced_pppdist()
  heavy <- .sirup_spaced_ot_points(heavy_points, "heavy_points", x_col, y_col)
  light <- .sirup_spaced_ot_points(light_points, "light_points", x_col, y_col)
  heavy_ids <- sort(unique(heavy$clone_id))
  light_ids <- sort(unique(light$clone_id))
  win <- .sirup_spaced_ppp_window(rbind(heavy, light))
  nproc <- .sirup_spaced_check_nproc(nproc)

  distance_row <- function(i) {
    heavy_id <- heavy_ids[[i]]
    heavy_ppp <- .sirup_spaced_as_ppp(heavy[heavy$clone_id == heavy_id, , drop = FALSE], win)
    row <- numeric(length(light_ids))
    for (j in seq_along(light_ids)) {
      light_id <- light_ids[[j]]
      light_ppp <- .sirup_spaced_as_ppp(light[light$clone_id == light_id, , drop = FALSE], win)
      row[[j]] <- tryCatch(
        pppdist(heavy_ppp, light_ppp, type = "mat", matching = FALSE, auction = FALSE, cutoff = cutoff),
        error = function(e) {
          rlang::abort(
            paste0("`pppdist()` failed for heavy clone ", heavy_id, " and light clone ", light_id, ": ", conditionMessage(e)),
            parent = e
          )
        }
      )
    }
    row
  }

  rows <- .sirup_spaced_lapply(seq_along(heavy_ids), distance_row, nproc)
  dist_matrix <- do.call(rbind, rows)
  dimnames(dist_matrix) <- list(heavy_ids, light_ids)
  transformSpacedDistances(dist_matrix, scale = scale)
}

.sirup_spaced_check_nproc <- function(nproc) {
  nproc <- as.integer(nproc)
  if (is.na(nproc) || nproc < 1L) {
    rlang::abort("`nproc` must be a positive integer.")
  }
  nproc
}

.sirup_spaced_lapply <- function(x, fun, nproc) {
  if (nproc == 1L || length(x) <= 1L) {
    return(lapply(x, fun))
  }
  if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(nproc)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    return(parallel::parLapply(cl, x, fun))
  }
  parallel::mclapply(x, fun, mc.cores = nproc)
}

.sirup_spaced_ot_points <- function(points, name, x_col, y_col) {
  points <- .as_data_frame(points)
  .check_required_cols(points, c("clone_id", x_col, y_col), name)
  points$clone_id <- as.character(points$clone_id)
  points$x <- as.numeric(points[[x_col]])
  points$y <- as.numeric(points[[y_col]])
  points <- unique(points[
    !is.na(points$clone_id) &
      points$clone_id != "" &
      is.finite(points$x) &
      is.finite(points$y),
    c("clone_id", "x", "y"),
    drop = FALSE
  ])
  if (nrow(points) == 0L) {
    rlang::abort(paste0("`", name, "` has no usable point-pattern coordinates."))
  }
  points
}

.sirup_spaced_pppdist <- function() {
  if (!requireNamespace("spatstat.geom", quietly = TRUE)) {
    rlang::abort("Package 'spatstat.geom' is required for `method = \"ot\"`.")
  }
  if (exists("pppdist", envir = asNamespace("spatstat.geom"), inherits = FALSE)) {
    return(get("pppdist", envir = asNamespace("spatstat.geom"), inherits = FALSE))
  }
  if (requireNamespace("spatstat.explore", quietly = TRUE) &&
      exists("pppdist", envir = asNamespace("spatstat.explore"), inherits = FALSE)) {
    return(get("pppdist", envir = asNamespace("spatstat.explore"), inherits = FALSE))
  }
  rlang::abort("Could not find `pppdist()` in installed spatstat packages.")
}

.sirup_spaced_ppp_window <- function(points) {
  xrange <- range(points$x, finite = TRUE)
  yrange <- range(points$y, finite = TRUE)
  if (diff(xrange) == 0) {
    xrange <- xrange + c(-0.5, 0.5)
  }
  if (diff(yrange) == 0) {
    yrange <- yrange + c(-0.5, 0.5)
  }
  spatstat.geom::owin(xrange = xrange, yrange = yrange)
}

.sirup_spaced_as_ppp <- function(points, window) {
  spatstat.geom::ppp(x = points$x, y = points$y, window = window)
}

.sirup_spaced_prune_top_n <- function(mat, top_n) {
  if (top_n >= ncol(mat)) {
    return(mat)
  }
  out <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  for (i in seq_len(nrow(mat))) {
    row <- mat[i, ]
    positive <- which(row > 0)
    if (length(positive) == 0L) {
      next
    }
    keep <- positive[order(row[positive], decreasing = TRUE)[seq_len(min(top_n, length(positive)))]]
    out[i, keep] <- row[keep]
  }
  out
}

.sirup_spaced_row_normalize <- function(mat) {
  rows <- rowSums(mat)
  keep <- rows > 0
  if (any(keep)) {
    mat[keep, ] <- mat[keep, , drop = FALSE] / rows[keep]
  }
  mat
}

.sirup_spaced_col_normalize <- function(mat) {
  cols <- colSums(mat)
  keep <- cols > 0
  if (any(keep)) {
    mat[, keep] <- t(t(mat[, keep, drop = FALSE]) / cols[keep])
  }
  mat
}

.sirup_spaced_top_light_candidates <- function(spatial_affinity,
                                               mapping_matrix,
                                               score,
                                               normalized_pairing,
                                               light_specificity,
                                               weight,
                                               top_n_light_chains) {
  out <- vector("list", nrow(score))
  for (i in seq_len(nrow(score))) {
    ord <- order(score[i, ], decreasing = TRUE)
    keep <- ord[seq_len(top_n_light_chains)]
    index <- cbind(rep(i, length(keep)), keep)
    out[[i]] <- data.frame(
      weight = weight,
      heavy_id = rownames(score)[i],
      light_id = colnames(score)[keep],
      spatial_affinity_score = spatial_affinity[index],
      repair_score = mapping_matrix[index],
      combined_score = score[index],
      normalized_pairing_score = normalized_pairing[index],
      light_specificity_score = light_specificity[index],
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

.sirup_spaced_prepare_mapping <- function(mapping_matrix, spatial_affinity) {
  if (is.null(mapping_matrix)) {
    mapping_matrix <- matrix(
      0,
      nrow = nrow(spatial_affinity),
      ncol = ncol(spatial_affinity),
      dimnames = dimnames(spatial_affinity)
    )
    return(mapping_matrix)
  }
  mapping_matrix <- .sirup_spaced_matrix_arg(mapping_matrix, "mapping_matrix")
  if (!identical(dim(mapping_matrix), dim(spatial_affinity))) {
    rlang::abort(paste0(
      "`mapping_matrix` dimensions must match `spatial_affinity`. Found ",
      nrow(mapping_matrix),
      "x",
      ncol(mapping_matrix),
      " mapping and ",
      nrow(spatial_affinity),
      "x",
      ncol(spatial_affinity),
      " spatial."
    ))
  }
  if (.sirup_spaced_has_matching_dimnames(mapping_matrix, spatial_affinity)) {
    if (anyDuplicated(rownames(mapping_matrix)) || anyDuplicated(colnames(mapping_matrix))) {
      rlang::abort("`mapping_matrix` rownames and colnames must be unique clone IDs.")
    }
    mapping_matrix <- mapping_matrix[rownames(spatial_affinity), colnames(spatial_affinity), drop = FALSE]
  } else {
    # TODO: replace SPaCeD's positional convention with mandatory named model
    # matrices once the REPAIR export is stable across environments.
    dimnames(mapping_matrix) <- dimnames(spatial_affinity)
  }
  .sirup_spaced_row_normalize(mapping_matrix)
}
