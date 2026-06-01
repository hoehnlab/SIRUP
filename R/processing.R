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
