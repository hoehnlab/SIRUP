#' Process spatial AIRR seq data by filtering for productivity, v/j/junction length, a minimum 
#' consensus count, optional max opening gap length, and optionally adding a full_id combining sequence_id and
#' spatial barcode.
#'
#' @param airr_df AIRR-style data frame.
#' @param loci Character vector of loci to keep. Use `NULL` to skip locus
#'   filtering.
#' @param min_v_length Minimum accepted V gene length.
#' @param min_j_length Minimum accepted J gene length.
#' @param min_junction_length Minimum accepted junction length.
#' @param productive_only Keep only productive rearrangements.
#' @param min_consensus_count Optional minimum consensus_count.
#' @param max_opening_gap Optional maximum opening gap in `sequence_alignment`.
#' @param sequence_alignment_col Sequence alignment column for gap filtering.
#' @param sequence_id_col Column containing read or consensus identifiers.
#' @param barcode_col Optional name of a spatial barcode column that is already
#'   present in `airr_df`.
#' @param add_full_id Add a combined id with the original sequence id 
#' and the spatial barcode. barcode_col must be provided if TRUE.
#' @param full_id_barcode_key Barcode field name to use in `full_id`.
#' @return Preprocessed data frame with `vlength`, `jlength`, and optional
#'   `full_id` columns.
#' @export
processSpatialAirr <- function(airr_df,
                                 loci = NULL,
                                 min_v_length = 100,
                                 min_j_length = 10,
                                 min_junction_length = 1,
                                 productive_only = TRUE,
                                 min_consensus_count = 2,
                                 max_opening_gap = NULL,
                                 sequence_alignment_col = "sequence_alignment",
                                 sequence_id_col = "sequence_id",
                                 barcode_col = NULL,
                                 add_full_id = FALSE,
                                 full_id_barcode_key = "barcode2um") {
  required <- c(
    "v_sequence_end", "v_sequence_start",
    "j_sequence_end", "j_sequence_start",
    "junction_length"
  )
  .check_required_cols(airr_df, required, "airr_df")
  out <- .as_data_frame(airr_df)
  n_input <- nrow(out)
  out$vlength <- out$v_sequence_end - out$v_sequence_start
  out$jlength <- out$j_sequence_end - out$j_sequence_start

  if (!is.null(loci)) {
    .check_required_cols(out, "locus", "airr_df")
    out <- out[out$locus %in% loci, , drop = FALSE]
  }

  out <- out[
    !is.na(out$vlength) &
      !is.na(out$jlength) &
      !is.na(out$junction_length) &
      out$vlength >= min_v_length &
      out$jlength >= min_j_length &
      out$junction_length >= min_junction_length,
    ,
    drop = FALSE
  ]

  if (productive_only && "productive" %in% colnames(out)) {
    out <- out[isTRUE(out$productive) | out$productive %in% c(TRUE, "TRUE", "T"), , drop = FALSE]
  }
  if (!is.null(min_consensus_count) && "consensus_count" %in% colnames(out)) {
    out <- out[!is.na(out$consensus_count) & out$consensus_count >= min_consensus_count, , drop = FALSE]
  }
  if (!is.null(max_opening_gap)) {
    .check_required_cols(out, sequence_alignment_col, "airr_df")
    gaps <- vapply(out[[sequence_alignment_col]], openingGapLength, integer(1))
    out <- out[gaps <= max_opening_gap, , drop = FALSE]
  }

  if (!is.null(barcode_col)) {
    .check_required_cols(out, barcode_col, "airr_df")
  }

  if (add_full_id) {
    if (is.null(barcode_col)) {
      rlang::abort("barcode_col is required when add_full_id is TRUE.")
    }
    .check_required_cols(out, c(sequence_id_col, barcode_col), "airr_df")
    out$full_id <- paste0(
      "sequence_id=", out[[sequence_id_col]],
      "|", full_id_barcode_key, "=", out[[barcode_col]]
    )
  }

  message(
    "processSpatialAirr kept ", nrow(out), " of ", n_input,
    " rows (dropped ", n_input - nrow(out), ")."
  )
  out
}

#' Check the maximum internal gap run length
#'
#' @param seq Character sequence, typically `sequence_alignment`.
#' @return Integer length of the longest run of periods.
#' @export
maxGapLength <- function(seq) {
  seq <- as.character(seq)
  if (is.na(seq)) {
    return(NA_integer_)
  }
  gaps <- gregexpr("\\.+", seq)
  matches <- regmatches(seq, gaps)[[1]]
  if (length(matches) == 0 || identical(matches, character(0))) {
    return(0L)
  }
  max(nchar(matches))
}

#' Check the opening gap run length
#'
#' @param seq Character sequence, typically `sequence_alignment`.
#' @return Integer length of the leading run of periods.
#' @export
openingGapLength <- function(seq) {
  seq <- as.character(seq)
  if (is.na(seq)) {
    return(NA_integer_)
  }
  gap_match <- regexpr("^\\.+", seq)
  if (gap_match[[1]] == -1) {
    return(0L)
  }
  attr(gap_match, "match.length")
}

#' Collapse to the most abundant consensus sequence per barcode and locus
#'
#' Adds `prop_total_consensus`, `non_N_nucleotide_prop`, `max_gap_length`, and
#' `opening_gap_length`, then keeps the highest-consensus row per barcode/locus.
#'
#' @param airr_df AIRR-style data frame.
#' @param barcode_col,locus_col,consensus_col Column names for grouping/ranking.
#' @param sequence_col Sequence column for non-N proportion.
#' @param alignment_col Alignment column for gap metrics.
#' @param ties Keep all tied top rows or only the first tied row.
#' @return Collapsed data frame.
#' @export
collapseTopConsensusPerBarcode <- function(airr_df,
                                               barcode_col = "barcode_2um",
                                               locus_col = "locus",
                                               consensus_col = "consensus_count",
                                               sequence_col = "sequence",
                                               alignment_col = "sequence_alignment",
                                               ties = c("all", "first")) {
  ties <- match.arg(ties)
  .check_required_cols(airr_df, c(barcode_col, locus_col, consensus_col), "airr_df")
  out <- .as_data_frame(airr_df)
  tie_summary <- out %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(barcode_col, locus_col)))) %>%
    dplyr::summarise(
      max_consensus = max(.data[[consensus_col]], na.rm = TRUE),
      n_top_consensus = sum(.data[[consensus_col]] == max(.data[[consensus_col]], na.rm = TRUE), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(.data$n_top_consensus > 1)
  if (nrow(tie_summary) > 0L) {
    message(
      "collapseTopConsensusPerBarcode found ", nrow(tie_summary),
      " barcode/locus groups with tied top consensus rows."
    )
  }
  out <- out %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(barcode_col, locus_col)))) %>%
    dplyr::mutate(prop_total_consensus = .data[[consensus_col]] / sum(.data[[consensus_col]], na.rm = TRUE)) %>%
    dplyr::filter(.data[[consensus_col]] == max(.data[[consensus_col]], na.rm = TRUE))

  if (ties == "first") {
    out <- out %>% dplyr::slice_head(n = 1)
  }

  out <- out %>% dplyr::ungroup()
  if (sequence_col %in% colnames(out)) {
    out$non_N_nucleotide_prop <- stringr::str_count(out[[sequence_col]], "[ACGT]") / nchar(out[[sequence_col]])
  }
  if (alignment_col %in% colnames(out)) {
    out$max_gap_length <- vapply(out[[alignment_col]], maxGapLength, integer(1))
    out$opening_gap_length <- vapply(out[[alignment_col]], openingGapLength, integer(1))
  }
  out
}

#' Summarize a spatial receptor table
#'
#' @param airr_df AIRR-style data frame.
#' @param barcode_col Barcode column to count.
#' @return One-row tibble of basic counts.
#' @export
spatialReceptorSummary <- function(airr_df, barcode_col = "barcode_2um") {
  .check_required_cols(airr_df, barcode_col, "airr_df")
  tibble::tibble(
    n_sequences = nrow(airr_df),
    n_barcodes = dplyr::n_distinct(airr_df[[barcode_col]], na.rm = TRUE),
    n_loci = if ("locus" %in% colnames(airr_df)) dplyr::n_distinct(airr_df$locus, na.rm = TRUE) else NA_integer_,
    n_productive = if ("productive" %in% colnames(airr_df)) sum(airr_df$productive %in% c(TRUE, "TRUE", "T"), na.rm = TRUE) else NA_integer_
  )
}

#' Determine a distance threshold via kernel density estimation to consider points as neighbors
#'
#' @param tissue_positions Spatial tissue positions data frame containing x and y columns.
#' @param num_neighbors_non_boundary Number of nearby points needed to be
#'   considered non-boundary.
#' @param x Column name for x coordinates in tissue_positions.
#' @param y Column name for y coordinates in tissue_positions.
#' @return Distance threshold to consider points as neighbors.
#' @export
calculateNeighborDistanceThreshold <- function(tissue_positions, num_neighbors_non_boundary, x = "x", y = "y") {
  dist_matrix <- tissue_positions %>%
    dplyr::select(dplyr::all_of(c(x, y))) %>%
    stats::dist(method = "euclidean") %>%
    as.matrix()

  neighbor_distances <- apply(dist_matrix, 1, function(row) {
    sorted_distances <- sort(row[row > 0])
    sorted_distances[num_neighbors_non_boundary]
  })

  neighbor_distances <- neighbor_distances[
    !is.na(neighbor_distances) &
      !is.nan(neighbor_distances) &
      !is.infinite(neighbor_distances)
  ]

  unique_distances <- unique(neighbor_distances)
  if (length(unique_distances) < 3) {
    stop(
      "The automated distance threshold calculation requires at least 3 unique distance values.\n",
      "Found distance values: ", paste(unique_distances, collapse = ", ")
    )
  }

  bandwidth <- kedd::h.ucv(unique_distances, 4)$h
  dens <- KernSmooth::bkde(neighbor_distances, bandwidth = bandwidth, canonical = TRUE)
  xdens <- dens$x
  ydens <- dens$y

  threshold <- tryCatch(
    xdens[which(diff(sign(diff(ydens))) == 2)[1] + 1],
    error = function(e) {
      warning("No minimum was found between two modes.")
      NULL
    }
  )

  threshold
}

#' Collapse spatially contiguous spots with identical sequences
#'
#' For each unique receptor sequence, finds connected components of neighboring
#' spots sharing that sequence, then collapses each component to a single
#' representative spot.
#'
#' @param data AIRR-style data frame.
#' @param tissue_positions Data frame with spatial coordinates for all spots.
#' @param barcode Column name for spot barcodes in both data frames.
#' @param sequence Column name for sequence in data.
#' @param consensus_count Column name for consensus count in data.
#' @param x Column name for x coordinates.
#' @param y Column name for y coordinates.
#' @param num_neighbors_non_boundary Number of nearby points needed to be
#'   considered non-boundary.
#' @return `data` with bleedover rows removed, retaining only the representative
#'   spot for each contiguous group of identical sequences.
#' @export
collapseIdenticalNeighbors <- function(data,
                                       tissue_positions,
                                       barcode = "barcode",
                                       sequence = "sequence",
                                       consensus_count = "consensus_count",
                                       x = "x",
                                       y = "y",
                                       num_neighbors_non_boundary = 6L) {
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

  cli::cli_progress_step("Estimating distance threshold from full tissue positions...")

  all_coords <- as.matrix(tissue_positions[, c(x, y)])
  n_all <- nrow(all_coords)
  k_search <- min(num_neighbors_non_boundary + 1L, n_all)
  sample_size <- min(2000L, n_all)
  sample_idx <- sample(n_all, sample_size)

  knn_sample <- RANN::nn2(
    data = all_coords,
    query = all_coords[sample_idx, , drop = FALSE],
    k = k_search,
    treetype = "kd",
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
  dens <- KernSmooth::bkde(neighbor_distances, bandwidth = bandwidth, canonical = TRUE)
  xdens <- dens$x
  ydens <- dens$y

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

  cli::cli_progress_step("Filtering to spots with observed sequences...")

  barcodes <- unique(data[[barcode]])
  tissue_data <- tissue_positions[tissue_positions[[barcode]] %in% barcodes, ]
  coords <- as.matrix(tissue_data[, c(x, y)])
  barcodes_ordered <- tissue_data[[barcode]]
  n <- nrow(tissue_data)
  n_pretty <- scales::comma(n)
  n_all_pretty <- scales::comma(n_all)
  cli::cli_progress_done()
  cli::cli_alert_info("Sequence-positive spots: {n_pretty} (of {n_all_pretty} total)")

  cli::cli_progress_step("Building adjacency matrix over spots with observed sequences...")

  knn_full <- RANN::nn2(
    data = coords,
    query = coords,
    k = k_search,
    treetype = "kd",
    searchtype = "standard"
  )

  row_idx <- rep(seq_len(n), times = k_search)
  col_idx <- as.vector(knn_full$nn.idx)
  dists <- as.vector(knn_full$nn.dists)

  valid <- col_idx > 0 & dists <= threshold & dists > 0
  row_idx <- row_idx[valid]
  col_idx <- col_idx[valid]

  adj_sparse <- Matrix::sparseMatrix(
    i = c(row_idx, col_idx),
    j = c(col_idx, row_idx),
    x = TRUE,
    dims = c(n, n),
    dimnames = list(barcodes_ordered, barcodes_ordered)
  )

  cli::cli_progress_done()

  barcode_to_idx <- stats::setNames(seq_len(n), barcodes_ordered)
  unique_sequences <- unique(data[[sequence]])
  n_seqs <- length(unique_sequences)
  component_labels <- vector("list", n_seqs)

  cli::cli_progress_bar(
    name = "Finding connected components",
    total = n_seqs,
    format = "{cli::pb_name} {cli::pb_bar} {cli::pb_current}/{cli::pb_total} | ETA: {cli::pb_eta}"
  )

  for (si in seq_len(n_seqs)) {
    seq_i <- unique_sequences[[si]]
    seq_barcodes <- data[[barcode]][data[[sequence]] == seq_i]
    seq_idx <- barcode_to_idx[seq_barcodes]

    if (length(seq_idx) == 1L) {
      component_labels[[si]] <- data.frame(
        barcode = seq_barcodes,
        sequence = seq_i,
        component = 1L
      )
    } else {
      sub_adj <- adj_sparse[seq_idx, seq_idx, drop = FALSE]
      g <- igraph::graph_from_adjacency_matrix(sub_adj, mode = "undirected", diag = FALSE)
      components <- igraph::components(g)

      component_labels[[si]] <- data.frame(
        barcode = seq_barcodes,
        sequence = seq_i,
        component = components$membership
      )
    }

    cli::cli_progress_update()
  }

  cli::cli_progress_done()

  component_df <- do.call(rbind, component_labels)

  cli::cli_progress_step("Collapsing components to representative spots...")

  coords_df <- data.frame(barcodes_ordered, cx = coords[, 1], cy = coords[, 2])
  colnames(coords_df)[1] <- barcode

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
      if (nrow(grp) == 1L) {
        return(grp)
      }

      max_count <- max(grp[[consensus_count]])
      grp <- grp[grp[[consensus_count]] == max_count, , drop = FALSE]
      if (nrow(grp) == 1L) {
        return(grp)
      }

      centroid_x <- mean(grp$cx)
      centroid_y <- mean(grp$cy)
      dist_centroid <- sqrt((grp$cx - centroid_x)^2 + (grp$cy - centroid_y)^2)
      grp <- grp[dist_centroid == min(dist_centroid), , drop = FALSE]
      if (nrow(grp) == 1L) {
        return(grp)
      }

      grp[sample(nrow(grp), 1L), , drop = FALSE]
    }
  )

  representative_df <- do.call(rbind, representative_rows)
  cli::cli_progress_done()

  cli::cli_progress_step("Filtering data to representative spots...")

  keep_keys <- paste(representative_df[[barcode]], representative_df[[sequence]], sep = "___")
  all_keys <- paste(data[[barcode]], data[[sequence]], sep = "___")
  result <- data[all_keys %in% keep_keys, ]
  n_result_pretty <- scales::comma(nrow(result))
  n_data_pretty <- scales::comma(nrow(data))
  n_removed_pretty <- scales::comma(nrow(data) - nrow(result))

  cli::cli_progress_done()
  cli::cli_alert_success(
    "Done! Retained {n_result_pretty} rows from {n_data_pretty} ({n_removed_pretty} bleedover sequences removed)"
  )

  result
}
