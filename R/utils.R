#' Determine a distance threshold via kernel density estimation to consider points as neighbors
#' 
#' @param    tissue_positions               Spatial tissue positions data frame containing x and y columns and loc_names as row names
#' @param    num_neighbors_non_boundary     Number of nearby points needed to be considered non-boundary 
#' @param    x                              Column name for x coordinates in tissue_positions
#' @param    y                              Column name for y coordinates in tissue_positions
#' 
#' @return  Distance threshold to consider points as neighbors
#' 
#' @export
calculateNeighborDistanceThreshold <- function(tissue_positions, num_neighbors_non_boundary, x = "x", y = "y"){
  # Calculate distance matrix
  dist_matrix <- tissue_positions %>%
    dplyr::select(dplyr::all_of(c(x, y))) %>%
    stats::dist(method = "euclidean") %>%
    as.matrix()

  # For each point, get the distance to its num_neighbors_non_boundary-th nearest neighbor
  neighbor_distances <- apply(dist_matrix, 1, function(row) {
    sorted_distances <- sort(row[row > 0])  # Exclude distance to self
    return(sorted_distances[num_neighbors_non_boundary])
  })
  
  # Remove any NAs, NaNs, or Infs
  neighbor_distances <- neighbor_distances[!is.na(neighbor_distances) & !is.nan(neighbor_distances) & !is.infinite(neighbor_distances)]

  # Need at least 3 unique distances for density estimation
  unique_distances <- unique(neighbor_distances)
  if (length(unique_distances) < 3 ) {
    stop("The automated distance threshold calculation requires at least 3 unique distance values.\n",
         "Found distance values: ", paste(unique_distances, collapse=", "))
  }
  
  # Estimate ideal bandwidth
  bandwidth <- kedd::h.ucv(unique_distances, 4)$h
  
  # Kernel density estimation
  dens <- KernSmooth::bkde(neighbor_distances, bandwidth = bandwidth, canonical = TRUE)
  xdens <- dens$x
  ydens <- dens$y
  
  # Find distance threshold
  tryCatch(threshold <- xdens[which(diff(sign(diff(ydens))) == 2)[1] + 1], 
           error = function(e) {
             warning('No minimum was found between two modes.')
             return(NULL) })
  
  return(threshold)
}

#' Collapse spatially contiguous spots with identical sequences
#'
#' For each unique BCR sequence, finds connected components of neighboring spots
#' sharing that sequence, then collapses each component to a single representative
#' spot (highest consensus count -> closest to centroid -> random tiebreak).
#'
#' @param data airr.tsv style data frame
#' @param tissue_positions Data frame with spatial coordinates for all spots
#' @param barcode Column name for spot barcodes (must be present in both data frames)
#' @param sequence Column name for sequence in data
#' @param consensus_count Column name for consensus count in data
#' @param x Column name for x coordinates in tissue_positions
#' @param y Column name for y coordinates in tissue_positions
#' @param num_neighbors_non_boundary Number of nearby points needed to be considered non-boundary
#'
#' @return bcr_data with bleedover rows removed, retaining only the representative
#'         spot for each contiguous group of identical sequences
#' @export
collapseSpatial <- function(
    data,
    tissue_positions,
    barcode                    = "barcode",
    sequence                   = "sequence",
    consensus_count            = "consensus_count",
    x                          = "x",
    y                          = "y",
    num_neighbors_non_boundary = 6L
) {

  # --- Step 0: Collapse to one row per barcode/sequence, summing consensus counts and averaging x and y --- #
  # This can happen if we used a higher-res barcode to build consensus sequences earlier
  cli::cli_progress_step("Ensuring one row per barcode/sequence...")
  data <- data %>%
  dplyr::group_by(dplyr::across(dplyr::all_of(c(barcode, sequence)))) %>%
  dplyr::summarise(
    dplyr::across(dplyr::all_of(consensus_count), sum),
    dplyr::across(dplyr::all_of(x), mean),
    dplyr::across(dplyr::all_of(y), mean),
    dplyr::across(-dplyr::all_of(consensus_count), dplyr::first),
    .groups = "drop"
  )
  cli::cli_progress_done()


  # --- Step 1: Estimate distance threshold from full tissue positions ---
  # ~1 min for 750k spots
  cli::cli_progress_step("Estimating distance threshold from full tissue positions...")

  all_coords  <- as.matrix(tissue_positions[, c(x, y)])
  n_all       <- nrow(all_coords)
  k_search    <- min(num_neighbors_non_boundary + 1L, n_all)
  sample_size <- min(2000L, n_all)
  sample_idx  <- sample(n_all, sample_size)

  knn_sample <- RANN::nn2(
    data       = all_coords,
    query      = all_coords[sample_idx, , drop = FALSE],
    k          = k_search,
    treetype   = "kd",
    searchtype = "standard"
  )

  neighbor_distances <- knn_sample$nn.dists[, k_search]
  neighbor_distances <- neighbor_distances[
    !is.na(neighbor_distances) &
    !is.nan(neighbor_distances) &
    !is.infinite(neighbor_distances) &
    neighbor_distances > 0
  ]

  unique_distances <- unique(neighbor_distances)
  if (length(unique_distances) < 3) {
    stop(
      "The automated distance threshold calculation requires at least 3 unique distance values.\n",
      "Found distance values: ", paste(unique_distances, collapse = ", ")
    )
  }

  bandwidth <- kedd::h.ucv(unique_distances, 4)$h
  dens      <- KernSmooth::bkde(neighbor_distances, bandwidth = bandwidth, canonical = TRUE)
  xdens     <- dens$x
  ydens     <- dens$y

  threshold <- tryCatch(
    xdens[which(diff(sign(diff(ydens))) == 2)[1] + 1],
    error = function(e) {
      warning("No minimum was found between two modes.")
      NULL
    }
  )

  if (is.null(threshold) || is.na(threshold)) {
    stop("Could not determine a valid distance threshold from the density estimate.")
  }

  cli::cli_progress_done()
  cli::cli_alert_info("Distance threshold: {round(threshold, 4)}")

  # --- Step 2: Restrict to seq-positive spots ---
  cli::cli_progress_step("Filtering to spots with observed sequences...")

  barcodes         <- unique(data[[barcode]])
  tissue_data      <- tissue_positions[tissue_positions[[barcode]] %in% barcodes, ]
  coords           <- as.matrix(tissue_data[, c(x, y)])
  barcodes_ordered <- tissue_data[[barcode]]
  n                <- nrow(tissue_data)
  n_pretty         <- scales::comma(n)
  n_all_pretty     <- scales::comma(n_all)
  cli::cli_progress_done()
  cli::cli_alert_info("Sequence-positive spots: {n_pretty} (of {n_all_pretty} total)")

  # --- Step 3: Build sparse adjacency matrix over BCR-positive spots only ---
  cli::cli_progress_step("Building adjacency matrix over spots with observed sequences...")

  knn_full <- RANN::nn2(
    data       = coords,
    query      = coords,
    k          = k_search,
    treetype   = "kd",
    searchtype = "standard"
  )

  row_idx <- rep(seq_len(n), times = k_search)
  col_idx <- as.vector(knn_full$nn.idx)
  dists   <- as.vector(knn_full$nn.dists)

  valid   <- col_idx > 0 & dists <= threshold & dists > 0
  row_idx <- row_idx[valid]
  col_idx <- col_idx[valid]

  adj_sparse <- Matrix::sparseMatrix(
    i        = c(row_idx, col_idx),
    j        = c(col_idx, row_idx),
    x        = TRUE,
    dims     = c(n, n),
    dimnames = list(barcodes_ordered, barcodes_ordered)
  )

  cli::cli_progress_done()

  # --- Step 4: Find connected components per unique sequence ---
  barcode_to_idx   <- stats::setNames(seq_len(n), barcodes_ordered)
  unique_sequences <- unique(data[[sequence]])
  n_seqs           <- length(unique_sequences)
  component_labels <- vector("list", n_seqs)

  cli::cli_progress_bar(
    name   = "Finding connected components",
    total  = n_seqs,
    format = "{cli::pb_name} {cli::pb_bar} {cli::pb_current}/{cli::pb_total} | ETA: {cli::pb_eta}"
  )

  for (si in seq_len(n_seqs)) {
    seq_i        <- unique_sequences[[si]]
    seq_barcodes <- data[[barcode]][data[[sequence]] == seq_i]
    seq_idx      <- barcode_to_idx[seq_barcodes]

    if (length(seq_idx) == 1L) {
      component_labels[[si]] <- data.frame(
        barcode   = seq_barcodes,
        sequence  = seq_i,
        component = 1L
      )
    } else {
      sub_adj    <- adj_sparse[seq_idx, seq_idx, drop = FALSE]
      g          <- igraph::graph_from_adjacency_matrix(sub_adj, mode = "undirected", diag = FALSE)
      components <- igraph::components(g)

      component_labels[[si]] <- data.frame(
        barcode   = seq_barcodes,
        sequence  = seq_i,
        component = components$membership
      )
    }

    cli::cli_progress_update()
  }

  cli::cli_progress_done()

  component_df <- do.call(rbind, component_labels)

  # --- Step 5: Collapse each component to its representative spot ---
  cli::cli_progress_step("Collapsing components to representative spots...")

  coords_df <- data.frame(barcodes_ordered, cx = coords[, 1], cy = coords[, 2])
  colnames(coords_df)[1] <- c(barcode)

  colnames(component_df)[colnames(component_df) == "barcode"] <- barcode
  colnames(component_df)[colnames(component_df) == "sequence"] <- sequence
  component_df <- merge(component_df, coords_df, by = barcode)
  component_df <- merge(
    component_df,
    data[, c(barcode, sequence, consensus_count)],
    by = c(barcode, sequence)
  )

  representative_rows <- lapply(
    split(component_df, list(component_df[[sequence]], component_df$component), drop = TRUE),
    function(grp) {
      if (nrow(grp) == 1L) return(grp)

      # Tiebreaker 1: highest consensus count
      max_count <- max(grp[[consensus_count]])
      grp       <- grp[grp[[consensus_count]] == max_count, , drop = FALSE]
      if (nrow(grp) == 1L) return(grp)

      # Tiebreaker 2: closest to centroid
      centroid_x    <- mean(grp$cx)
      centroid_y    <- mean(grp$cy)
      dist_centroid <- sqrt((grp$cx - centroid_x)^2 + (grp$cy - centroid_y)^2)
      grp           <- grp[dist_centroid == min(dist_centroid), , drop = FALSE]
      if (nrow(grp) == 1L) return(grp)

      # Tiebreaker 3: random
      grp[sample(nrow(grp), 1L), , drop = FALSE]
    }
  )

  representative_df <- do.call(rbind, representative_rows)
  cli::cli_progress_done()

  # --- Step 6: Filter data to only representative rows ---
  cli::cli_progress_step("Filtering data to representative spots...")

  keep_keys <- paste(representative_df[[barcode]], representative_df[[sequence]], sep = "___")
  all_keys  <- paste(data[[barcode]], data[[sequence]], sep = "___")
  result    <- data[all_keys %in% keep_keys, ]
  n_result_pretty <- scales::comma(nrow(result))
  n_data_pretty <- scales::comma(nrow(data))
  n_removed_pretty <- scales::comma(nrow(data) - nrow(result))

  cli::cli_progress_done()
  cli::cli_alert_success(
    "Done! Retained {n_result_pretty} rows from {n_data_pretty} ({n_removed_pretty} bleedover rows removed)"
  )

  return(result)
}


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

.as_data_frame <- function(x) {
  if (data.table::is.data.table(x)) {
    tibble::as_tibble(x)
  } else {
    as.data.frame(x, stringsAsFactors = FALSE)
  }
}

.require_namespace <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    rlang::abort(paste0("Package '", pkg, "' is required for this operation."))
  }
  invisible(TRUE)
}

.ensure_dir <- function(path) {
  dir <- dirname(path)
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}
