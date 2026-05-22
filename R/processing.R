#' Filter receptor reads by locus and alignment quality
#'
#' @param airr_df AIRR-style data frame.
#' @param loci Character vector of loci to keep. Use `NULL` to skip locus
#'   filtering, or set `exclude_loci = TRUE` to remove these loci.
#' @param exclude_loci Whether `loci` should be excluded rather than retained.
#' @param min_v_length,min_j_length,min_junction_length Minimum accepted lengths.
#' @param productive_only Keep only productive rearrangements when a `productive`
#'   column is present.
#' @param min_consensus_count Optional minimum `consensus_count`.
#' @param max_opening_gap Optional maximum opening gap in `sequence_alignment`.
#' @param sequence_alignment_col Sequence alignment column for gap filtering.
#' @return Filtered data frame with `vlength` and `jlength` columns.
#' @export
filterReceptorReads <- function(airr_df,
                                  loci = NULL,
                                  exclude_loci = FALSE,
                                  min_v_length = 100,
                                  min_j_length = 10,
                                  min_junction_length = 1,
                                  productive_only = FALSE,
                                  min_consensus_count = NULL,
                                  max_opening_gap = NULL,
                                  sequence_alignment_col = "sequence_alignment") {
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
    if (exclude_loci) {
      out <- out[!(out$locus %in% loci), , drop = FALSE]
    } else {
      out <- out[out$locus %in% loci, , drop = FALSE]
    }
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
  out
}

#' Maximum internal gap run length
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

#' Opening gap run length
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

#' Collapse to the most abundant consensus per barcode and locus
#'
#' Adds `prop_total_consensus`, `non_N_nucleotide_prop`, `maxGapLength`, and
#' `openingGapLength`, then keeps the highest-consensus row per barcode/locus.
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
  out <- out |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(barcode_col, locus_col)))) |>
    dplyr::mutate(prop_total_consensus = .data[[consensus_col]] / sum(.data[[consensus_col]], na.rm = TRUE)) |>
    dplyr::filter(.data[[consensus_col]] == max(.data[[consensus_col]], na.rm = TRUE))

  if (ties == "first") {
    out <- out |> dplyr::slice_head(n = 1)
  }

  out <- out |> dplyr::ungroup()
  if (sequence_col %in% colnames(out)) {
    out$non_N_nucleotide_prop <- stringr::str_count(out[[sequence_col]], "[ACGT]") / nchar(out[[sequence_col]])
  }
  if (alignment_col %in% colnames(out)) {
    out$maxGapLength <- vapply(out[[alignment_col]], maxGapLength, integer(1))
    out$openingGapLength <- vapply(out[[alignment_col]], openingGapLength, integer(1))
  }
  out
}

#' Preprocess an IgBLAST AIRR table into spatial AIRR and FASTA outputs
#'
#' This wraps the first step of the original processing script: locus/length QC,
#' Visium HD barcode mapping, and writing AIRR TSV plus FASTA for pRESTO.
#'
#' @param input AIRR data frame or path to TSV.
#' @param barcode_map Barcode mapping data frame or path.
#' @param loci Loci to keep, for example `"IGH"` or `"TRB"`.
#' @param exclude_loci Whether to exclude `loci`.
#' @param out_tsv,out_fasta Optional output paths.
#' @param min_v_length,min_j_length,min_junction_length Length cutoffs.
#' @param sequence_col Column used for FASTA sequences.
#' @return Processed spatial AIRR data frame.
#' @export
preprocessSpatialAirr <- function(input,
                                    barcode_map,
                                    loci,
                                    exclude_loci = FALSE,
                                    out_tsv = NULL,
                                    out_fasta = NULL,
                                    min_v_length = 100,
                                    min_j_length = 10,
                                    min_junction_length = 1,
                                    sequence_col = "sequence_alignment") {
  airr_df <- if (is.character(input) && length(input) == 1) {
    readSpatialAirr(input)
  } else {
    input
  }
  barcode_map <- if (is.character(barcode_map) && length(barcode_map) == 1) {
    readBarcodeMapping(barcode_map)
  } else {
    barcode_map
  }

  out <- filterReceptorReads(
    airr_df,
    loci = loci,
    exclude_loci = exclude_loci,
    min_v_length = min_v_length,
    min_j_length = min_j_length,
    min_junction_length = min_junction_length
  )
  out <- addSpatialBarcodes(out, barcode_map, mode = "raw_suffix")
  out$full_id <- paste0(
    "sequence_id=", out$sequence_id,
    "|barcode2um=", out$barcode_2um,
    "|barcode8um=", out$barcode_8um
  )

  writeSpatialOutputs(out, out_tsv = out_tsv, out_fasta = out_fasta, sequence_col = sequence_col)
  out
}

#' Read an AIRR TSV file
#'
#' Uses `airr::read_rearrangement()` when available, falling back to
#' `data.table::fread()`.
#'
#' @param path Path to AIRR TSV.
#' @return Data frame.
#' @export
readSpatialAirr <- function(path) {
  if (requireNamespace("airr", quietly = TRUE)) {
    return(.as_data_frame(airr::read_rearrangement(path)))
  }
  .as_data_frame(data.table::fread(path, sep = "\t"))
}

#' Write spatial AIRR and FASTA outputs
#'
#' @param airr_df AIRR-style data frame.
#' @param out_tsv,out_fasta Optional output paths.
#' @param id_col FASTA identifier column.
#' @param sequence_col FASTA sequence column.
#' @return Invisibly returns `airr_df`.
#' @export
writeSpatialOutputs <- function(airr_df,
                                  out_tsv = NULL,
                                  out_fasta = NULL,
                                  id_col = "full_id",
                                  sequence_col = "sequence_alignment") {
  if (!is.null(out_tsv)) {
    .ensure_dir(out_tsv)
    if (requireNamespace("airr", quietly = TRUE)) {
      airr::write_rearrangement(airr_df, out_tsv)
    } else {
      data.table::fwrite(airr_df, out_tsv, sep = "\t")
    }
  }
  if (!is.null(out_fasta)) {
    .check_required_cols(airr_df, c(id_col, sequence_col), "airr_df")
    .ensure_dir(out_fasta)
    .require_namespace("seqinr")
    seqinr::write.fasta(
      sequences = as.list(airr_df[[sequence_col]]),
      names = airr_df[[id_col]],
      file.out = out_fasta
    )
  }
  invisible(airr_df)
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
