#' Filter AIRR-seq data
#'
#' Filters AIRR-seq data for productivity, V/J/junction length, a minimum
#' consensus count, optional max opening gap length, and optionally adding a
#' spatial sequence ID combining `sequence_id` and a spatial ID.
#'
#' @param airr_df AIRR-style data frame.
#' @param locus_val Character vector of locus values to keep. Use `NULL` to skip
#'   filtering.
#' @param min_v_length Minimum accepted V gene length.
#' @param min_j_length Minimum accepted J gene length.
#' @param min_junction_length Minimum accepted junction length.
#' @param productive_only Keep only productive rearrangements.
#' @param complete_vdj Keep only rows with a complete V(D)J. IgBLAST writes
#'   this column as `"T"`/`"F"`. Skipped when the column is absent.
#' @param min_consensus_count Optional minimum consensus_count.
#' @param max_opening_gap Optional maximum opening gap in `sequence_alignment`.
#' @param max_gap_length Optional maximum gap run in `sequence_alignment`. A real
#'   somatic-hypermutation indel is a multiple of three and rarely longer than
#'   12 nt, so a longer run points to a chimeric or truncated consensus. The
#'   metric is measured on the IMGT-numbered alignment, where the padding is
#'   locus-specific: IGHV is padded at 6/9/12 nt and TRBV at 18/21 nt, so a
#'   ceiling tuned for IGH rejects every TRB sequence. Set it per locus.
#' @param sequence_alignment Sequence alignment column for gap filtering.
#' @param id Column containing read or consensus identifiers.
#' @param spatial_id Optional name of a spatial ID column that is already
#'   present in `airr_df`.
#' @param add_spatial_sequence_id Add a combined ID with the original sequence ID
#'   and spatial ID. `spatial_id` must be provided if `TRUE`.
#' @param spatial_sequence_id Name of the combined spatial sequence ID column.
#' @param spatial_id_key Spatial ID field name to use in the combined spatial
#'   sequence ID.
#' @return Preprocessed data frame with `v_sequence_length`,
#'   `j_sequence_length`, and optional spatial sequence ID columns.
#' @export
filterAirr <- function(airr_df,
                       locus_val = NULL,
                       min_v_length = 100,
                       min_j_length = 10,
                       min_junction_length = 1,
                       productive_only = TRUE,
                       complete_vdj = FALSE,
                       min_consensus_count = 2,
                       max_opening_gap = NULL,
                       max_gap_length = NULL,
                       sequence_alignment = "sequence_alignment",
                       id = "sequence_id",
                       spatial_id = NULL,
                       add_spatial_sequence_id = FALSE,
                       spatial_sequence_id = "spatial_sequence_id",
                       spatial_id_key = "spatial_id") {
  required <- c(
    "v_sequence_end", "v_sequence_start",
    "j_sequence_end", "j_sequence_start",
    "junction_length"
  )
  .check_required_cols(airr_df, required, "airr_df")
  # Whether the caller asked for these, as opposed to inheriting the default.
  asked_productive <- !missing(productive_only)
  asked_complete <- !missing(complete_vdj)
  asked_consensus <- !missing(min_consensus_count)
  out <- .as_data_frame(airr_df)
  n_input <- nrow(out)
  out$v_sequence_length <- out$v_sequence_end - out$v_sequence_start + 1L
  out$j_sequence_length <- out$j_sequence_end - out$j_sequence_start + 1L

  if (!is.null(locus_val)) {
    .check_required_cols(out, "locus", "airr_df")
    out <- out[out$locus %in% locus_val, , drop = FALSE]
  }

  out <- out[
    !is.na(out$v_sequence_length) &
      !is.na(out$j_sequence_length) &
      !is.na(out$junction_length) &
      out$v_sequence_length >= min_v_length &
      out$j_sequence_length >= min_j_length &
      out$junction_length >= min_junction_length,
    ,
    drop = FALSE
  ]

  if (productive_only) {
    if ("productive" %in% colnames(out)) {
      out <- out[out$productive %in% c(TRUE, "TRUE", "T"), , drop = FALSE]
    } else {
      if (asked_productive) .sirup_warn_missing_filter_col("productive_only", "productive")
    }
  }
  if (complete_vdj) {
    if ("complete_vdj" %in% colnames(out)) {
      out <- out[out$complete_vdj %in% c(TRUE, "TRUE", "T"), , drop = FALSE]
    } else {
      if (asked_complete) .sirup_warn_missing_filter_col("complete_vdj", "complete_vdj")
    }
  }
  if (!is.null(min_consensus_count)) {
    if ("consensus_count" %in% colnames(out)) {
      out <- out[!is.na(out$consensus_count) & out$consensus_count >= min_consensus_count, , drop = FALSE]
    } else {
      if (asked_consensus) .sirup_warn_missing_filter_col("min_consensus_count", "consensus_count")
    }
  }
  if (!is.null(max_opening_gap)) {
    .check_required_cols(out, sequence_alignment, "airr_df")
    gaps <- .sirup_opening_gap_length(out[[sequence_alignment]])
    out <- out[gaps <= max_opening_gap, , drop = FALSE]
  }
  if (!is.null(max_gap_length)) {
    .check_required_cols(out, sequence_alignment, "airr_df")
    inner_gaps <- .sirup_max_gap_length(out[[sequence_alignment]])
    out <- out[inner_gaps <= max_gap_length, , drop = FALSE]
  }

  if (!is.null(spatial_id)) {
    .check_required_cols(out, spatial_id, "airr_df")
  }

  if (add_spatial_sequence_id) {
    if (is.null(spatial_id)) {
      rlang::abort("`spatial_id` is required when `add_spatial_sequence_id` is TRUE.")
    }
    .check_required_cols(out, c(id, spatial_id), "airr_df")
    # paste0() treats zero-length inputs as "", so guard the 0-row case.
    out[[spatial_sequence_id]] <- if (nrow(out) == 0L) {
      character(0)
    } else {
      paste0(
        "sequence_id=", out[[id]],
        "|", spatial_id_key, "=", out[[spatial_id]]
      )
    }
  }

  message(
    "filterAirr kept ", nrow(out), " of ", n_input,
    " rows (dropped ", n_input - nrow(out), ")."
  )
  out
}

#' Process one spatial sample from raw reads to an AIRR table
#'
#' This function runs the full external processing pipeline for one sample
#' and one receptor locus group (`ig` or `tr`). The `input` can be one of:
#'
#' * Raw reads: a FASTA or FASTQ file, optionally gzipped, as written by the
#'   rrad package. Export rrad FASTA with `to_fasta = TRUE` and
#'   `header_format = "presto"` so barcode/UMI tags remain in the IgBLAST
#'   identifier. The function first annotates the reads with
#'   [runIgblast()] (`igblastn -outfmt 19`), which needs `igblast_db_dir`
#'   (or the per-segment database paths).
#' * An AIRR table: a `.tsv` path or a data frame. The raw IgBLAST step is
#'   skipped. A data frame input requires `output_dir`.
#' * A pRESTO-ready spatial FASTA (`input_type = "spatial_fasta"`): the
#'   IgBLAST and filter steps are both skipped and the FASTA goes straight
#'   to the consensus pipeline below.
#'
#' For raw reads and AIRR input, the function then:
#'
#' * splits rrad `|BARCODE=...|UMI=...` tags (or legacy colon tags) out of
#'   `sequence_id` with
#'   [parseRadSequenceIds()];
#' * removes junk sequences with [filterAirr()]. Pipeline defaults:
#'   `complete_vdj = TRUE`, `spatial_id = "barcode"`, and a combined header
#'   ID keyed by `cluster_field`. Override with `filter_args`;
#' * writes the raw and filtered AIRR tables to `output_dir` for
#'   inspection;
#' * writes the filtered sequences to a spatial FASTA whose headers carry
#'   the spatial ID.
#'
#' The spatial FASTA then goes through the external Immcantation/pRESTO
#' consensus workflow:
#'
#' * `ClusterSets.py set` groups reads within each spatial ID using `cluster_field`,
#'   the requested identity threshold `identity`, and vsearch.
#' * `ParseHeaders.py merge` is used to write a new field into the headers that combines
#'   the spatial ID and the cluster id (provided with parse_fields). These are
#'   combined using `parse_act` (default set) and written to a new header field `parse_key`.
#'   Uses `delimiters` to specify the FASTA header delimiters.
#' * `AlignSets.py muscle` aligns reads within each `parse_key` group created by `ParseHeaders`.
#'   Takes `muscle_exec` as the path to the muscle aligner (https://github.com/rcedgar/muscle).
#'   Optionally logs results to `align_log`. Uses `delimiters` to specify the FASTA header delimiters.
#' * `BuildConsensus.py` creates consensus sequences within each `parse_key` group. `consensus_min_count`
#'   specifies the minimum number of sequences needed to define a valid consensus sequence, and
#'   `consensus_freq` specifies the fraction of character occurrences under which an ambiguous
#'   character is assigned. Optionally logs results to `consensus_log`. If there are additional
#'   fields in the FASTA header that should be carried forward and copied to the consensus file,
#'   `consensus_copy_fields` specifies the field names and `consensus_copy_action` specifies how
#'   to copy them (majority, first, last, or all).
#' * `AssignGenes.py igblast` assigns V(D)J gene annotations using IgBLAST. Takes either `ig` or `tr`
#'   as the argument to `receptor` (for bcr/tcr respectively). Takes `organism` as the argument to
#'   specify the species (e.g., human, mouse, etc.). `igblast_base` specifies the path to the IgBLAST
#'   installation, and `igblast_exec` can be provided if the igblast executable is not located in
#'   the `bin/igblastn` path beneath `igblast_base`.
#' * `MakeDb.py igblast` converts IgBLAST output to an AIRR-style TSV using the specified IMGT
#'   germline reference `germline_glob`.
#'
#' Each command is printed before it is run, and any non-zero exit status stops
#' the function with the failing command in the error message. All commands
#' run through `system2()` without a shell. The final output is the
#' `db_pass` AIRR TSV in the returned path list.
#'
#' @param input Raw reads (FASTA/FASTQ, optionally gzipped), an AIRR TSV
#'   path, an AIRR data frame, or a pRESTO-ready spatial FASTA.
#' @param output_dir Directory for all outputs. Defaults to the input file's
#'   directory. Required when `input` is a data frame.
#' @param receptor Either `ig` for BCR or `tr` for TCR IgBLAST loci.
#' @param organism The organism for which to perform the analysis.
#' @param input_type How to treat `input`. The default `"auto"` detects the
#'   type from the object class and file extension: FASTA/FASTQ paths become
#'   `"raw_reads"`. Pass `"spatial_fasta"` to feed a pRESTO-ready FASTA
#'   straight to the consensus pipeline.
#'   Newly demultiplexed rrad FASTA should use `"auto"` or `"raw_reads"`:
#'   `"spatial_fasta"` skips the first read-level V(D)J annotation and filtering.
#' @param igblast_db_dir Directory with the raw-step germline databases
#'   (`IG_V`, `IG_D`, `IG_J`, `IG_C`, `TR_V`, `TR_D`, `TR_J`). See
#'   [runIgblast()].
#' @param aux_file Optional path passed to `igblastn -auxiliary_data`.
#' @param germline_db_v,germline_db_d,germline_db_j,c_region_db Optional
#'   per-segment database paths that override the `igblast_db_dir`
#'   convention.
#' @param filter_args Named list of [filterAirr()] overrides. Pipeline
#'   defaults: `complete_vdj = TRUE`, `min_consensus_count = NULL`,
#'   `spatial_id = "barcode"`, `add_spatial_sequence_id = TRUE`,
#'   `spatial_id_key = cluster_field`. An explicit `NULL` disables an
#'   optional filter (`locus_val`, `min_consensus_count`, `max_opening_gap`,
#'   or `spatial_id`); other parameters do not accept `NULL`.
#' @param invalid_spatial_ids What to do with filtered rows whose spatial
#'   barcode is missing or contains a delimiter: `"abort"` (default) stops
#'   with a count, `"drop"` removes those rows with a warning.
#' @param fasta_sequence_col AIRR column written to the spatial FASTA.
#'   Defaults to `sequence_alignment` (the V(D)J-aligned region), which
#'   keeps flanking read junk out of clustering and consensus calling. Use
#'   `sequence` for the full read.
#' @param sample_name Stem for intermediate file names. Defaults to the
#'   input file name without extension, or `"sample"` for data frames.
#' @param identity Cluster identity threshold passed to `ClusterSets.py`.
#' @param consensus_freq Consensus frequency threshold passed to
#'   `BuildConsensus.py`.
#' @param germline_glob IMGT germline glob passed to `MakeDb.py -r`.
#' @param igblast_base Path passed to `AssignGenes.py igblast -b`.
#' @param igblast_exec Optional path passed to `AssignGenes.py --exec`. When
#'   `NULL`, defaults to `file.path(igblast_base, "bin", "igblastn")`.
#' @param muscle_exec Path passed to `AlignSets.py --exec`.
#' @param nproc Number of processors for the raw IgBLAST step, pRESTO, and
#'   AssignGenes.
#' @param cluster_field FASTA header field used for within-spot clustering.
#'   The pipeline writes headers keyed by this name; the default `"barcode"`
#'   uses the parsed rrad `BARCODE` field. Override it for other spatial ID fields, for
#'   example `barcode_8um`.
#' @param delimiters Character vector of pRESTO delimiter values passed
#'   after `--delim`. Defaults to `c("|", "=", ",")`.
#' @param parse_fields Header fields passed to `ParseHeaders.py merge -f`.
#' @param parse_key Header field name created by `ParseHeaders.py merge -k`.
#' @param parse_act Action passed to `ParseHeaders.py merge --act`. The default
#'   `"set"` joins the spatial ID and the cluster ID with a comma, so the spatial
#'   ID stays recoverable from the consensus header. `"cat"` joins them with no
#'   separator, which makes the spatial ID unrecoverable for any variable-length
#'   barcode. Both produce the same grouping, so prefer `"set"`.
#' @param consensus_min_count Minimum sequence count passed to
#'   `BuildConsensus.py -n`. Recommended to be at least 2 for non-single-cell resolution data.
#'   This gates which groups are written, not the consensus of a group that
#'   passes, so a run at `1` can be filtered to any higher value afterwards with
#'   `filterAirr(min_consensus_count = )`. Note `BuildConsensus.py -q` is fixed at
#'   20 and has no effect here, because the pipeline input is FASTA and carries
#'   no base qualities.
#' @param umi_support Group reads by spatial barcode AND UMI before the barcode
#'   step. Reads of one molecule that differ by a sequencing error cluster apart
#'   at `identity = 1.0` and are then dropped as singletons; this stage joins
#'   them first, at `umi_identity`, so they support each other. Each
#'   barcode-UMI cluster yields one consensus carrying `cc_bc_umi` (its reads,
#'   counted before the productive subset) and `n_productive_reads` (the reads
#'   the consensus was built from; a cluster with none is discarded). The
#'   barcode step then runs with `BuildConsensus -n 1` and sums both fields, and
#'   `min_umi_support` replaces the read floor. Needs an AIRR table with `umi`
#'   and `productive` columns.
#' @param umi_field Column holding the UMI, used with `umi_support`.
#' @param umi_identity Cluster identity for the barcode-UMI step. Looser than
#'   the barcode step on purpose, so error-bearing reads of one molecule group.
#' @param min_umi_support Lowest summed `cc_bc_umi` a barcode consensus needs.
#'   Applied after `MakeDb.py`, writing a filtered copy beside `db_pass`, so a
#'   run can be re-thresholded without reprocessing.
#' @param consensus_copy_fields Optional fields passed to `BuildConsensus.py
#'   --cf`. With `umi_support`, `CC_BC_UMI` and `N_PRODUCTIVE_READS` are
#'   prepended and summed.
#' @param consensus_copy_action Action passed to `BuildConsensus.py --act` when
#'   `consensus_copy_fields` is provided.
#' @param consensus_max_error Maximum error rate passed to
#'   `BuildConsensus.py --maxerror`.
#' @param consensus_max_gap Maximum gap fraction passed to
#'   `BuildConsensus.py --maxgap`.
#' @param skip_existing Skip a pipeline command when its output already exists
#'   and is not empty. Use this to restart an interrupted run, or to reuse one
#'   clustering and alignment pass across several `consensus_freq` values,
#'   since those steps depend only on `identity`.
#' @param align_log Optional log path passed to `AlignSets.py --log`.
#' @param consensus_log Optional log path passed to `BuildConsensus.py --log`.
#' @return A named list, invisibly: `raw_airr`, `filtered_airr`,
#'   `spatial_fasta`, `n_raw`, `n_filtered` (each `NULL` for steps that did
#'   not run), plus the pRESTO/IgBLAST output paths. `db_pass` is the final
#'   AIRR TSV. With `umi_support`, `db_pass_umi` is that table filtered to
#'   `min_umi_support`.
#' @export
processSpatialSample <- function(input,
                                 output_dir = NULL,
                                 receptor = c("ig", "tr"),
                                 organism = "human",
                                 input_type = c("auto", "raw_reads", "airr", "spatial_fasta"),
                                 igblast_db_dir = NULL,
                                 aux_file = NULL,
                                 germline_db_v = NULL,
                                 germline_db_d = NULL,
                                 germline_db_j = NULL,
                                 c_region_db = NULL,
                                 filter_args = list(),
                                 invalid_spatial_ids = c("abort", "drop"),
                                 fasta_sequence_col = "sequence_alignment",
                                 sample_name = NULL,
                                 identity = 0.975,
                                 consensus_freq = 0.51,
                                 germline_glob,
                                 igblast_base = "~/software/igblast",
                                 igblast_exec = NULL,
                                 muscle_exec = "~/software/muscle3.8.31_i86linux64",
                                 nproc = 8,
                                 cluster_field = "barcode",
                                 delimiters = c("|", "=", ","),
                                 parse_fields = c(cluster_field, "cluster"),
                                 parse_key = paste0(cluster_field, "_cluster"),
                                 parse_act = "set",
                                 consensus_min_count = 1,
                                 umi_support = FALSE,
                                 umi_field = "umi",
                                 umi_identity = 0.95,
                                 min_umi_support = 2,
                                 consensus_copy_fields = NULL,
                                 consensus_copy_action = "majority",
                                 consensus_max_error = 0.1,
                                 consensus_max_gap = 0.3,
                                 skip_existing = FALSE,
                                 align_log = NULL,
                                 consensus_log = NULL) {
  receptor <- match.arg(receptor)
  input_type <- match.arg(input_type)
  invalid_spatial_ids <- match.arg(invalid_spatial_ids)
  if (length(cluster_field) != 1L || is.na(cluster_field) || !nzchar(cluster_field)) {
    rlang::abort("`cluster_field` must name one FASTA header field, for example `barcode` or `barcode_8um`.")
  }
  if (missing(germline_glob)) {
    rlang::abort("`germline_glob` is required. It names the IMGT germline FASTA files for `MakeDb.py -r`.")
  }
  if (is.null(igblast_exec)) {
    igblast_exec <- file.path(igblast_base, "bin", "igblastn")
  }

  if (isTRUE(umi_support) && identical(input_type, "spatial_fasta")) {
    rlang::abort("`umi_support = TRUE` needs an AIRR table, because the UMI and productive columns are read from it.")
  }

  if (input_type == "auto") {
    input_type <- .sirup_detect_input_type(input)
  } else if (input_type == "airr") {
    input_type <- if (is.data.frame(input)) "airr_df" else "airr_tsv"
  } else if (input_type == "spatial_fasta" && is.data.frame(input)) {
    rlang::abort("`input_type = \"spatial_fasta\"` needs a FASTA file path, not a data frame.")
  }

  if (is.data.frame(input)) {
    if (is.null(output_dir)) {
      rlang::abort("`output_dir` is required when `input` is a data frame.")
    }
  } else if (is.null(output_dir)) {
    output_dir <- dirname(input)
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (is.null(sample_name)) {
    sample_name <- if (is.character(input)) {
      tools::file_path_sans_ext(sub("\\.gz$", "", basename(input), ignore.case = TRUE))
    } else {
      "sample"
    }
  }

  raw_airr_tsv <- NULL
  filtered_tsv <- NULL
  n_raw <- NULL
  n_filtered <- NULL

  if (input_type == "spatial_fasta") {
    spatial_fasta <- input
  } else {
    if (input_type == "raw_reads") {
      raw_airr_tsv <- file.path(output_dir, sprintf("%s_%s_igblast_raw.tsv", sample_name, receptor))
      runIgblast(
        input = input,
        output_tsv = raw_airr_tsv,
        receptor = receptor,
        organism = organism,
        igblast_db_dir = igblast_db_dir,
        igblast_exec = igblast_exec,
        aux_file = aux_file,
        germline_db_v = germline_db_v,
        germline_db_d = germline_db_d,
        germline_db_j = germline_db_j,
        c_region_db = c_region_db,
        num_threads = nproc
      )
      airr <- .sirup_read_airr(raw_airr_tsv)
    } else if (input_type == "airr_tsv") {
      airr <- .sirup_read_airr(input)
    } else {
      airr <- .as_data_frame(input)
    }
    n_raw <- nrow(airr)

    airr <- parseRadSequenceIds(airr)

    filter_defaults <- list(
      complete_vdj = TRUE,
      min_consensus_count = NULL,
      spatial_id = "barcode",
      add_spatial_sequence_id = TRUE,
      spatial_id_key = cluster_field
    )
    # Explicit NULL disables an optional filter; other filterAirr parameters
    # cannot accept NULL, so reject those instead of passing them through.
    args <- utils::modifyList(filter_defaults, filter_args, keep.null = TRUE)
    nullable <- c("locus_val", "min_consensus_count", "max_opening_gap", "spatial_id")
    null_args <- names(args)[vapply(args, is.null, logical(1))]
    bad_nulls <- setdiff(null_args, nullable)
    if (length(bad_nulls)) {
      rlang::abort(paste0(
        "`filter_args` sets NULL for parameter",
        if (length(bad_nulls) == 1L) "" else "s", " that cannot be NULL: ",
        paste0("`", bad_nulls, "`", collapse = ", "),
        ". Only ", paste0("`", nullable, "`", collapse = ", "),
        " accept NULL. Omit the entry to use the default."
      ))
    }
    spatial_col <- args$spatial_id
    filtered <- do.call(filterAirr, c(list(airr_df = airr), args))
    n_filtered <- nrow(filtered)
    if (n_filtered == 0L) {
      rlang::abort("No sequences passed `filterAirr()`. Check the input and `filter_args`.")
    }
    if (!is.null(spatial_col)) {
      .check_required_cols(filtered, spatial_col, "filtered AIRR table")
      spatial_ids <- as.character(filtered[[spatial_col]])
      # A barcode that contains a pRESTO delimiter or whitespace corrupts the
      # FASTA header fields, so validate against the same `delimiters` the
      # pipeline passes to pRESTO.
      invalid <- which(is.na(spatial_ids) | !nzchar(spatial_ids) |
                         grepl(.sirup_delims_regex_class(delimiters), spatial_ids))
      if (length(invalid)) {
        detail <- paste0(
          length(invalid), " of ", n_filtered, " filtered rows have a missing ",
          "or invalid spatial barcode in `", spatial_col, "` (first at row ",
          invalid[[1L]], ")."
        )
        if (identical(invalid_spatial_ids, "drop")) {
          rlang::warn(paste0(detail, " Dropping those rows."))
          filtered <- filtered[-invalid, , drop = FALSE]
          n_filtered <- nrow(filtered)
          if (n_filtered == 0L) {
            rlang::abort("No filtered rows have a valid spatial barcode.")
          }
        } else {
          rlang::abort(paste0(
            detail,
            " For rrad reads, export FASTA with `to_fasta = TRUE, ",
            "header_format = \"presto\"` so `|BARCODE=...` survives IgBLAST. ",
            "For AIRR tables, add spatial IDs first (see `extractSpatialIds()`), ",
            "or set `invalid_spatial_ids = \"drop\"` to remove those rows."
          ))
        }
      }
    }
    filtered_tsv <- file.path(output_dir, sprintf("%s_%s_filtered.tsv", sample_name, receptor))
    data.table::fwrite(filtered, filtered_tsv, sep = "\t")

    header_col <- args$spatial_sequence_id
    if (is.null(header_col)) {
      header_col <- "spatial_sequence_id"
    }
    .check_required_cols(filtered, c(header_col, fasta_sequence_col), "filtered AIRR table")
    dup <- duplicated(filtered[[header_col]])
    if (any(dup)) {
      rlang::abort(paste0(
        sum(dup), " duplicate spatial FASTA header",
        if (sum(dup) == 1L) "" else "s", " in `", header_col, "` (first: '",
        filtered[[header_col]][dup][[1L]], "'). Sequence IDs must be unique ",
        "after tag stripping; IDs that differ only after a `|` or whitespace ",
        "collapse together."
      ))
    }
    if (isTRUE(umi_support)) {
      # One consensus per barcode-UMI cluster, carrying its read support, so the
      # barcode stage below weighs molecules rather than reads.
      umi_out <- .sirup_umi_consensus(
        df = filtered,
        output_dir = output_dir,
        sample_name = sample_name,
        receptor = receptor,
        cluster_field = cluster_field,
        umi_field = umi_field,
        umi_identity = umi_identity,
        consensus_freq = consensus_freq,
        consensus_min_count = consensus_min_count,
        consensus_max_error = consensus_max_error,
        consensus_max_gap = consensus_max_gap,
        delimiters = delimiters,
        muscle_exec = muscle_exec,
        nproc = nproc,
        fasta_sequence_col = fasta_sequence_col,
        # filterAirr's `id` default, unless the caller overrode it.
        id_col = if (is.null(args$id)) "sequence_id" else args$id,
        align_log = if (is.null(align_log)) NULL else sub("\\.log$", "_umi.log", align_log),
        consensus_log = if (is.null(consensus_log)) NULL else sub("\\.log$", "_umi.log", consensus_log)
      )
      spatial_fasta <- umi_out$spatial_fasta
    } else {
      spatial_fasta <- file.path(output_dir, sprintf("%s_%s_spatial.fasta", sample_name, receptor))
      writeFastqas(
        data.frame(
          id = filtered[[header_col]],
          seq = filtered[[fasta_sequence_col]],
          stringsAsFactors = FALSE
        ),
        fn = spatial_fasta,
        type = "fa",
        header_col = "id"
      )
      message("Wrote ", n_filtered, " filtered sequences to ", spatial_fasta)
    }
  }

  germline_files <- Sys.glob(path.expand(germline_glob))
  if (length(germline_files) == 0L) {
    rlang::abort(paste0("`germline_glob` matched no files: ", germline_glob))
  }

  # With UMI support the read floor moves off BuildConsensus: every molecule
  # consensus is kept (-n 1) and the summed `cc_bc_umi` decides afterwards.
  barcode_min_count <- consensus_min_count
  if (isTRUE(umi_support)) {
    barcode_min_count <- 1
    support_fields <- c("CC_BC_UMI", "N_PRODUCTIVE_READS")
    consensus_copy_action <- c(
      rep("sum", length(support_fields)),
      rep(consensus_copy_action, length.out = length(consensus_copy_fields))
    )
    consensus_copy_fields <- c(support_fields, consensus_copy_fields)
  }

  paths <- .presto_paths(
    spatial_fasta, identity, consensus_freq, output_dir,
    cluster_field = cluster_field,
    consensus_min_count = barcode_min_count
  )
  commands <- .sirup_presto_commands(
    paths = paths,
    receptor = receptor,
    organism = organism,
    identity = identity,
    consensus_freq = consensus_freq,
    germline_files = germline_files,
    igblast_base = igblast_base,
    igblast_exec = igblast_exec,
    muscle_exec = muscle_exec,
    nproc = nproc,
    cluster_field = cluster_field,
    delimiters = delimiters,
    parse_fields = parse_fields,
    parse_key = parse_key,
    parse_act = parse_act,
    consensus_min_count = barcode_min_count,
    consensus_copy_fields = consensus_copy_fields,
    consensus_copy_action = consensus_copy_action,
    consensus_max_error = consensus_max_error,
    consensus_max_gap = consensus_max_gap,
    align_log = align_log,
    consensus_log = consensus_log
  )

  tool_execs <- unique(vapply(commands, function(cmd) cmd$exec, character(1)))
  .sirup_check_tools(c(
    stats::setNames(tool_execs, basename(tool_execs)),
    muscle = muscle_exec,
    igblastn = igblast_exec
  ))
  for (cmd in commands) {
    if (skip_existing && !is.null(cmd$out) && file.exists(cmd$out) &&
        file.info(cmd$out)$size > 0) {
      message("Skipping ", cmd$exec, ": output exists (", cmd$out, ").")
      next
    }
    .sirup_run_command(cmd$exec, cmd$args)
  }

  db_pass_umi <- NULL
  if (isTRUE(umi_support) && file.exists(paths$db_pass)) {
    db <- .sirup_read_airr(paths$db_pass)
    sup_col <- intersect(c("cc_bc_umi", "CC_BC_UMI"), colnames(db))
    if (length(sup_col) == 0L) {
      rlang::warn("MakeDb wrote no `cc_bc_umi` column, so the `min_umi_support` floor was not applied.")
    } else {
      support <- suppressWarnings(as.numeric(db[[sup_col[[1L]]]]))
      keep <- !is.na(support) & support >= min_umi_support
      db_pass_umi <- sub("\\.tsv$", sprintf("_umi_support_%s.tsv", min_umi_support), paths$db_pass)
      data.table::fwrite(db[keep, , drop = FALSE], db_pass_umi, sep = "\t", quote = FALSE, na = "")
      message("UMI support floor (>= ", min_umi_support, "): kept ", sum(keep),
              " of ", nrow(db), " consensus sequences -> ", db_pass_umi)
    }
  }

  invisible(c(
    list(
      raw_airr = raw_airr_tsv,
      filtered_airr = filtered_tsv,
      spatial_fasta = spatial_fasta,
      db_pass_umi = db_pass_umi,
      n_raw = n_raw,
      n_filtered = n_filtered
    ),
    paths
  ))
}

#' Process multiple spatial samples
#'
#' @param inputs List or character vector of inputs, one per sample. Each is
#'   passed to `processSpatialSample()` and can be a reads file, an AIRR TSV
#'   path, or an AIRR data frame.
#' @param output_dirs Optional vector of output directories, one per input.
#' @param sample_names Optional vector of sample names, one per input.
#' @param ... Passed to `processSpatialSample()`.
#' @return A list of output path lists, invisibly.
#' @export
processSpatialSamples <- function(inputs,
                                  output_dirs = NULL,
                                  sample_names = NULL,
                                  ...) {
  if (!is.list(inputs)) {
    inputs <- as.list(inputs)
  }
  if (is.null(output_dirs)) {
    output_dirs <- rep(NA_character_, length(inputs))
  }
  if (length(output_dirs) != length(inputs)) {
    rlang::abort("output_dirs must be NULL or the same length as inputs.")
  }
  if (!is.null(sample_names) && length(sample_names) != length(inputs)) {
    rlang::abort("sample_names must be NULL or the same length as inputs.")
  }
  result_names <- vapply(seq_along(inputs), function(i) {
    if (!is.null(sample_names)) {
      sample_names[[i]]
    } else if (is.character(inputs[[i]])) {
      basename(inputs[[i]])
    } else {
      paste0("sample_", i)
    }
  }, character(1))
  results <- vector("list", length(inputs))
  names(results) <- result_names
  for (i in seq_along(inputs)) {
    output_dir <- output_dirs[[i]]
    if (is.na(output_dir)) {
      output_dir <- NULL
    }
    message("Processing sample ", i, " of ", length(inputs), ": ", result_names[[i]])
    results[[i]] <- processSpatialSample(
      input = inputs[[i]],
      output_dir = output_dir,
      sample_name = if (is.null(sample_names)) NULL else sample_names[[i]],
      ...
    )
  }
  invisible(results)
}

#' Collapse to the most abundant consensus sequence per spatial ID and locus
#'
#' Adds `prop_total_consensus`, `non_N_nucleotide_prop`, `max_gap_length`, and
#' `opening_gap_length`, then keeps the highest-consensus row per spatial ID/locus.
#'
#' @param airr_df AIRR-style data frame.
#' @param spatial_id Spatial ID column name used for grouping.
#' @param locus Locus column name used for grouping.
#' @param consensus Consensus-count column name used for ranking rows.
#' @param sequence Sequence column for non-N proportion.
#' @param alignment Alignment column for gap metrics.
#' @param ties Keep all tied top rows or only the first tied row.
#' @return Collapsed data frame.
#' @export
collapseTopConsensusPerSpatialId <- function(airr_df,
                                           spatial_id = "spatial_id_2um",
                                           locus = "locus",
                                           consensus = "consensus_count",
                                           sequence = "sequence",
                                           alignment = "sequence_alignment",
                                           ties = c("all", "first")) {
  ties <- match.arg(ties)
  .check_required_cols(airr_df, c(spatial_id, locus, consensus), "airr_df")
  out <- .as_data_frame(airr_df)
  tie_summary <- out %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(spatial_id, locus)))) %>%
    dplyr::summarise(
      max_consensus = max(.data[[consensus]], na.rm = TRUE),
      n_top_consensus = sum(.data[[consensus]] == max(.data[[consensus]], na.rm = TRUE), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(.data$n_top_consensus > 1)
  if (nrow(tie_summary) > 0L) {
    message(
      "collapseTopConsensusPerSpatialId found ", nrow(tie_summary),
      " spatial ID/locus groups with tied top consensus rows."
    )
  }
  out <- out %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(spatial_id, locus)))) %>%
    dplyr::mutate(prop_total_consensus = .data[[consensus]] / sum(.data[[consensus]], na.rm = TRUE)) %>%
    dplyr::filter(.data[[consensus]] == max(.data[[consensus]], na.rm = TRUE))

  if (ties == "first") {
    out <- out %>% dplyr::slice_head(n = 1)
  }

  out <- out %>% dplyr::ungroup()
  addSequenceQualityMetrics(out, sequence = sequence, alignment = alignment)
}

#' Add per-sequence quality metrics to an AIRR table
#'
#' Adds the columns used to judge consensus quality, without dropping any rows.
#' Use this to inspect the metric distributions before choosing thresholds for
#' [filterAirr()]. Previously these columns could only be obtained as a side
#' effect of [collapseTopConsensusPerSpatialId()], which also discards rows.
#'
#' Added columns, when the source column is present:
#'
#' * `non_N_nucleotide_prop` -- the share of called bases (`A`, `C`, `G`, `T`)
#'   in `sequence`. A low value means an ambiguous consensus. This must be taken
#'   from the ungapped sequence; using the alignment would put IMGT gap
#'   characters into the denominator.
#' * `max_gap_length` -- the longest gap run in `alignment`.
#' * `opening_gap_length` -- the leading gap run in `alignment`.
#'
#' @param airr_df AIRR-style data frame.
#' @param sequence Sequence column for the called-base proportion.
#' @param alignment Alignment column for the gap metrics.
#' @return `airr_df` with the metric columns added.
#' @export
addSequenceQualityMetrics <- function(airr_df,
                                      sequence = "sequence",
                                      alignment = "sequence_alignment") {
  out <- .as_data_frame(airr_df)
  if (sequence %in% colnames(out)) {
    out$non_N_nucleotide_prop <- stringr::str_count(out[[sequence]], "[ACGT]") / nchar(out[[sequence]])
  }
  if (alignment %in% colnames(out)) {
    out$max_gap_length <- .sirup_max_gap_length(out[[alignment]])
    out$opening_gap_length <- .sirup_opening_gap_length(out[[alignment]])
  }
  out
}

#' Summarize a spatial receptor table
#'
#' @param airr_df AIRR-style data frame or a SIRUP object.
#' @param spatial_id Spatial ID column to count. When omitted for a SIRUP object,
#'   the spatial ID field recorded in its metadata is used.
#' @return One-row tibble of basic counts.
#' @export
spatialReceptorSummary <- function(airr_df, spatial_id = "spatial_id_2um") {
  spatial_id_missing <- missing(spatial_id)
  inputs <- .sirup_receptor_input(
    airr_df,
    if (isSIRUP(airr_df) && spatial_id_missing) NULL else spatial_id
  )
  airr_df <- inputs$receptors
  spatial_id <- inputs$spatial_id_col %||% "spatial_id_2um"
  .check_required_cols(airr_df, spatial_id, "airr_df")
  sample_col <- if (is.list(inputs$receptor_fields)) inputs$receptor_fields$sample else NULL
  n_spatial_ids <- if (!is.null(sample_col) && sample_col %in% colnames(airr_df)) {
    sample_values <- as.character(airr_df[[sample_col]])
    spatial_values <- as.character(airr_df[[spatial_id]])
    keep <- !is.na(sample_values) & nzchar(sample_values) &
      !is.na(spatial_values) & nzchar(spatial_values)
    length(unique(paste(sample_values[keep], spatial_values[keep], sep = "\r")))
  } else {
    dplyr::n_distinct(airr_df[[spatial_id]], na.rm = TRUE)
  }
  tibble::tibble(
    n_sequences = nrow(airr_df),
    n_spatial_ids = n_spatial_ids,
    n_loci = if ("locus" %in% colnames(airr_df)) dplyr::n_distinct(airr_df$locus, na.rm = TRUE) else NA_integer_,
    n_productive = if ("productive" %in% colnames(airr_df)) sum(airr_df$productive %in% c(TRUE, "TRUE", "T"), na.rm = TRUE) else NA_integer_
  )
}

#' Collapse spatially contiguous spots with identical receptor sequences
#'
#' For each unique receptor sequence, finds connected components of neighboring
#' spots sharing that sequence, then collapses each component to a single
#' representative spot. A SIRUP object is the preferred input: receptor IDs
#' are read from the SIRUP field mapping and spatial coordinates are read from
#' the stored Seurat object. If you are not using a SIRUP object, pass an
#' AIRR-style receptor table as `object` and the matching Seurat object as
#' `seurat_object`.
#'
#' When `object` is a SIRUP object, the returned value is the updated SIRUP
#' object with its receptor table replaced. When `object` is a receptor table,
#' the returned value is the collapsed receptor table.
#'
#' @param object A SIRUP object, or an AIRR-style receptor data frame.
#' @param seurat_object Optional Seurat object used only when `object` is a
#'   receptor data frame.
#' @param spatial_id Column in the receptor table matching Seurat spot/cell
#'   spatial IDs. When omitted for a SIRUP object, the recorded SIRUP spatial ID field
#'   is used. When omitted for a receptor table, a spatial ID-like column is
#'   guessed.
#' @param sequence Column containing receptor sequence strings.
#' @param consensus_count Column containing consensus counts used to choose
#'   the representative spot within a connected component.
#' @param images Optional Seurat image names to use for spatial coordinates.
#' @param image_scale Image scale passed to
#'   `SeuratObject::GetTissueCoordinates()`.
#' @param swap_coords Swap x/y coordinates after reading Seurat tissue
#'   coordinates.
#' @param num_neighbors_non_boundary Number of nearby points needed to be
#'   considered non-boundary.
#' @param validate Whether to validate the updated SIRUP object before
#'   returning it. Only used when `object` is a SIRUP object.
#' @return A SIRUP object or receptor data frame with bleedover rows removed,
#'   retaining only the representative spot for each contiguous group of
#'   identical sequences.
#' @export
collapseIdenticalNeighbors <- function(object,
                                       seurat_object = NULL,
                                       spatial_id = NULL,
                                       sequence = "sequence",
                                       consensus_count = "consensus_count",
                                       images = NULL,
                                       image_scale = NULL,
                                       swap_coords = NULL,
                                       num_neighbors_non_boundary = 6L,
                                       validate = TRUE) {
  is_sirup <- isSIRUP(object)
  if (is_sirup) {
    if (!is.null(seurat_object)) {
      rlang::abort("`seurat_object` should not be supplied when `object` is a SIRUP object.")
    }
    validateSIRUP(object)
    receptor_fields <- .sirup_receptor_fields(object)
    receptors <- object$receptors
    spatial_id <- spatial_id %||% receptor_fields$spatial_id
    sample_col <- .sirup_receptor_sample_col(receptors, receptor_fields, object$spatial)
    seurat_object <- object$spatial
  } else {
    if (is.null(seurat_object)) {
      rlang::abort("`seurat_object` is required when `object` is a receptor table.")
    }
    receptors <- .as_data_frame(object)
    spatial_id <- .guess_spatial_id_col(receptors, spatial_id)
    sample_col <- .sirup_receptor_sample_col(receptors, spatial = seurat_object)
  }
  options <- .sirup_spatial_options(
    object = if (is_sirup) object else seurat_object,
    image_scale = image_scale,
    swap_coords = swap_coords
  )
  image_scale <- options$image_scale
  swap_coords <- options$swap_coords

  coords <- .sirup_spatial_coordinate_table(
    seurat_object,
    images = images,
    image_scale = image_scale,
    swap_coords = swap_coords
  )
  collapsed <- .collapse_identical_neighbors_table(
    receptors = receptors,
    tissue_positions = coords,
    spatial_id_col = spatial_id,
    sample_col = sample_col,
    sequence_col = sequence,
    consensus_count_col = consensus_count,
    num_neighbors_non_boundary = num_neighbors_non_boundary
  )

  if (!is_sirup) {
    return(collapsed)
  }

  object$receptors <- tibble::as_tibble(collapsed)
  if (isTRUE(validate)) {
    validateSIRUP(object)
  }
  object
}

# Collapse identical neighboring receptors after SIRUP/Seurat inputs have been
# resolved to a receptor table and a coordinate table.
.collapse_identical_neighbors_table <- function(receptors,
                                                tissue_positions,
                                                spatial_id_col,
                                                sample_col = NULL,
                                                sequence_col,
                                                consensus_count_col,
                                                num_neighbors_non_boundary) {
  data <- .as_data_frame(receptors)
  .check_required_cols(data, c(spatial_id_col, sequence_col, consensus_count_col), "receptors")
  .check_required_cols(tissue_positions, c("cell", "x", "y"), "spatial coordinates")

  if (!is.null(sample_col) && sample_col %in% colnames(data) && "sample" %in% colnames(tissue_positions)) {
    samples <- intersect(
      unique(stats::na.omit(as.character(data[[sample_col]]))),
      unique(stats::na.omit(as.character(tissue_positions$sample)))
    )
    if (length(samples) == 0L) {
      rlang::abort("No receptor sample values matched names in the spatial coordinate table.")
    }
    out <- lapply(samples, function(sample) {
      sample_data <- data[as.character(data[[sample_col]]) == sample, , drop = FALSE]
      sample_tissue <- tissue_positions[as.character(tissue_positions$sample) == sample, , drop = FALSE]
      sample_tissue$sample <- NULL
      .collapse_identical_neighbors_table(
        receptors = sample_data,
        tissue_positions = sample_tissue,
        spatial_id_col = spatial_id_col,
        sample_col = NULL,
        sequence_col = sequence_col,
        consensus_count_col = consensus_count_col,
        num_neighbors_non_boundary = num_neighbors_non_boundary
      )
    })
    return(dplyr::bind_rows(out))
  }

  temp <- .collapse_identical_neighbors_temp_cols(data)
  data[[temp$spatial_id]] <- as.character(data[[spatial_id_col]])
  if (!is.null(sample_col) && sample_col %in% colnames(data) && "sample" %in% colnames(tissue_positions)) {
    data[[temp$sample]] <- as.character(data[[sample_col]])
  }
  data[[temp$sequence]] <- as.character(data[[sequence_col]])
  data[[temp$row]] <- seq_len(nrow(data))

  tissue_positions <- .as_data_frame(tissue_positions)
  tissue_positions <- tissue_positions[is.finite(tissue_positions$x) & is.finite(tissue_positions$y), , drop = FALSE]
  sample_keep <- if (temp$sample %in% colnames(data)) "sample" else character()
  tissue_keep <- unique(c(intersect(sample_keep, colnames(tissue_positions)), "cell", "x", "y"))
  tissue_positions <- tissue_positions[!duplicated(tissue_positions[, intersect(c("sample", "cell"), tissue_keep), drop = FALSE]), tissue_keep, drop = FALSE]
  colnames(tissue_positions)[colnames(tissue_positions) == "cell"] <- temp$spatial_id
  colnames(tissue_positions)[colnames(tissue_positions) == "sample"] <- temp$sample
  colnames(tissue_positions)[colnames(tissue_positions) == "x"] <- temp$x
  colnames(tissue_positions)[colnames(tissue_positions) == "y"] <- temp$y

  merge_by <- if (temp$sample %in% colnames(data) && temp$sample %in% colnames(tissue_positions)) {
    c(temp$sample, temp$spatial_id)
  } else {
    temp$spatial_id
  }
  data <- merge(
    data,
    tissue_positions,
    by = merge_by,
    all.x = FALSE,
    all.y = FALSE,
    sort = FALSE
  )
  if (nrow(data) == 0L) {
    rlang::abort("No receptor IDs matched finite spatial coordinates in `seurat_object`.")
  }

  # A spatial_id can contain duplicate rows for the same sequence. Collapse those
  # first so the spatial graph only has one observation per spatial_id/sequence.
  cli::cli_progress_step("Ensuring one row per spatial_id/sequence...")
  group_cols <- c(intersect(temp$sample, colnames(data)), temp$spatial_id, temp$sequence)
  data <- data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(consensus_count_col), sum),
      dplyr::across(dplyr::all_of(c(temp$x, temp$y)), mean),
      dplyr::across(-dplyr::all_of(c(consensus_count_col, temp$x, temp$y)), dplyr::first),
      .groups = "drop"
    )
  cli::cli_progress_done()

  cli::cli_progress_step("Estimating distance threshold from full tissue positions...")

  # Estimate the spacing between neighboring spots from the full tissue grid.
  # The first density minimum separates immediate spot-neighbors from longer
  # distances across gaps, image boundaries, or unobserved tissue.
  all_coords <- as.matrix(tissue_positions[, c(temp$x, temp$y)])
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
  if (length(unique_distances) == 0L) {
    rlang::abort("Could not find positive neighbor distances in the spatial coordinate table.")
  } else if (length(unique_distances) < 3L) {
    threshold <- max(unique_distances)
    cli::cli_alert_info(
      "Neighbor distances are nearly uniform; using {round(threshold, 4)} as the distance threshold."
    )
  } else {
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
  }

  if (is.null(threshold) || is.na(threshold)) {
    stop("Could not determine a valid distance threshold from the density estimate.")
  }

  cli::cli_progress_done()
  cli::cli_alert_info("Distance threshold: {round(threshold, 4)}")

  cli::cli_progress_step("Filtering to spots with observed sequences...")

  # Build the graph only on spots with receptor observations, but keep the
  # threshold learned from all tissue positions so boundary spots are not
  # over-connected.
  spatial_keys <- unique(data[, merge_by, drop = FALSE])
  tissue_key <- .sirup_key_vector(.sirup_collapse_key_table(tissue_positions, temp))
  observed_key <- .sirup_key_vector(.sirup_collapse_key_table(spatial_keys, temp))
  tissue_data <- tissue_positions[tissue_key %in% observed_key, ]
  coords <- as.matrix(tissue_data[, c(temp$x, temp$y)])
  spatial_ids_ordered <- .sirup_key_vector(.sirup_collapse_key_table(tissue_data, temp))
  n <- nrow(tissue_data)
  n_pretty <- scales::comma(n)
  n_all_pretty <- scales::comma(n_all)
  cli::cli_progress_done()
  cli::cli_alert_info("Sequence-positive spots: {n_pretty} (of {n_all_pretty} total)")

  if (n < 2L) {
    result <- data[, setdiff(colnames(data), unlist(temp, use.names = FALSE)), drop = FALSE]
    cli::cli_alert_success("Done! Fewer than two sequence-positive spots were found, so no neighboring spots were collapsed.")
    return(result)
  }

  cli::cli_progress_step("Building adjacency matrix over spots with observed sequences...")

  k_observed <- min(k_search, n)

  # Keep undirected edges only between nearby observed spots. Connected
  # components inside this graph are candidate bleedover regions.
  knn_full <- RANN::nn2(
    data = coords,
    query = coords,
    k = k_observed,
    treetype = "kd",
    searchtype = "standard"
  )

  row_idx <- rep(seq_len(n), times = k_observed)
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
    dimnames = list(spatial_ids_ordered, spatial_ids_ordered)
  )

  cli::cli_progress_done()

  spatial_id_to_idx <- stats::setNames(seq_len(n), spatial_ids_ordered)
  unique_sequences <- unique(data[[temp$sequence]])
  n_seqs <- length(unique_sequences)
  component_labels <- vector("list", n_seqs)

  cli::cli_progress_bar(
    name = "Finding connected components",
    total = n_seqs,
    format = "{cli::pb_name} {cli::pb_bar} {cli::pb_current}/{cli::pb_total} | ETA: {cli::pb_eta}"
  )

  for (si in seq_len(n_seqs)) {
    seq_i <- unique_sequences[[si]]
    seq_rows <- data[data[[temp$sequence]] == seq_i, , drop = FALSE]
    seq_spatial_ids <- .sirup_key_vector(.sirup_collapse_key_table(seq_rows, temp))
    seq_idx <- spatial_id_to_idx[seq_spatial_ids]

    # For each exact sequence, connected components identify spatially
    # contiguous spots carrying that sequence.
    if (length(seq_idx) == 1L) {
      component_labels[[si]] <- data.frame(
        .sirup_spatial_id = seq_spatial_ids,
        .sirup_sequence = seq_i,
        component = 1L
      )
    } else {
      sub_adj <- adj_sparse[seq_idx, seq_idx, drop = FALSE]
      g <- igraph::graph_from_adjacency_matrix(sub_adj, mode = "undirected", diag = FALSE)
      components <- igraph::components(g)

      component_labels[[si]] <- data.frame(
        .sirup_spatial_id = seq_spatial_ids,
        .sirup_sequence = seq_i,
        component = components$membership
      )
    }

    cli::cli_progress_update()
  }

  cli::cli_progress_done()

  component_df <- do.call(rbind, component_labels)

  cli::cli_progress_step("Collapsing components to representative spots...")

  coords_df <- data.frame(spatial_ids_ordered, cx = coords[, 1], cy = coords[, 2])
  colnames(coords_df)[1] <- ".sirup_collapse_key"

  colnames(component_df)[colnames(component_df) == ".sirup_spatial_id"] <- ".sirup_collapse_key"
  colnames(component_df)[colnames(component_df) == ".sirup_sequence"] <- temp$sequence
  component_df <- merge(component_df, coords_df, by = ".sirup_collapse_key")
  data$.sirup_collapse_key <- .sirup_key_vector(.sirup_collapse_key_table(data, temp))
  component_df <- merge(
    component_df,
    data[, c(".sirup_collapse_key", temp$sequence, consensus_count_col)],
    by = c(".sirup_collapse_key", temp$sequence)
  )

  representative_rows <- lapply(
    split(component_df, list(component_df[[temp$sequence]], component_df$component), drop = TRUE),
    function(grp) {
      if (nrow(grp) == 1L) {
        return(grp)
      }

      # Prefer the most abundant consensus. If counts tie, use the point
      # closest to the component centroid; remaining ties are random.
      max_count <- max(grp[[consensus_count_col]])
      grp <- grp[grp[[consensus_count_col]] == max_count, , drop = FALSE]
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

  keep_keys <- paste(representative_df$.sirup_collapse_key, representative_df[[temp$sequence]], sep = "___")
  all_keys <- paste(data$.sirup_collapse_key, data[[temp$sequence]], sep = "___")
  result <- data[all_keys %in% keep_keys, ]
  result <- result[, setdiff(colnames(result), c(unlist(temp, use.names = FALSE), ".sirup_collapse_key")), drop = FALSE]
  n_result_pretty <- scales::comma(nrow(result))
  n_data_pretty <- scales::comma(nrow(data))
  n_removed_pretty <- scales::comma(nrow(data) - nrow(result))

  cli::cli_progress_done()
  cli::cli_alert_success(
    "Done! Retained {n_result_pretty} rows from {n_data_pretty} ({n_removed_pretty} bleedover sequences removed)"
  )

  result
}

# Reserve internal column names before temporarily joining receptor rows to
# Seurat coordinates.
.collapse_identical_neighbors_temp_cols <- function(data) {
  temp <- c(
    spatial_id = ".sirup_collapse_barcode",
    sample = ".sirup_collapse_sample",
    sequence = ".sirup_collapse_sequence",
    row = ".sirup_collapse_row",
    x = ".sirup_collapse_x",
    y = ".sirup_collapse_y"
  )
  overlapping <- intersect(unname(temp), colnames(data))
  if (length(overlapping) > 0L) {
    rlang::abort(paste0(
      "Temporary collapse column",
      if (length(overlapping) == 1L) " " else "s ",
      paste0("'", overlapping, "'", collapse = ", "),
      " already present in `receptors`."
    ))
  }
  as.list(temp)
}

.sirup_collapse_key_table <- function(dat, temp) {
  out <- data.frame(.sirup_spatial_id = as.character(dat[[temp$spatial_id]]), stringsAsFactors = FALSE)
  if (temp$sample %in% colnames(dat)) {
    out$.sirup_sample <- as.character(dat[[temp$sample]])
  }
  out
}

# Build the expected pRESTO/IgBLAST output filenames from the input FASTA stem
# and the clustering/consensus thresholds.
.presto_paths <- function(input_fasta, identity, consensus_freq, output_dir = NULL,
                          cluster_field = "barcode", consensus_min_count = 1) {
  if (is.null(output_dir)) {
    output_dir <- dirname(input_fasta)
  }
  stem <- tools::file_path_sans_ext(basename(input_fasta))
  cluster <- file.path(
    output_dir,
    sprintf("%s_ident_%s_%s_cluster_cluster-pass.fasta", stem, identity, cluster_field)
  )
  reheader <- sub("\\.fasta$", "_reheader.fasta", cluster)
  aligned <- sub("\\.fasta$", "_align-pass.fasta", reheader)
  consensus <- file.path(
    output_dir,
    sprintf(
      "%s_ident_%s_%s_cluster_freq_%s_n_%s_cluster-pass_reheader_align-pass_consensus-pass.fasta",
      stem, identity, cluster_field, consensus_freq, consensus_min_count
    )
  )
  fmt7 <- sub("\\.fasta$", "_igblast.fmt7", consensus)
  db_pass <- sub("_igblast\\.fmt7$", "_igblast_db-pass.tsv", fmt7)
  list(
    input = input_fasta,
    cluster = cluster,
    reheader = reheader,
    aligned = aligned,
    consensus = consensus,
    fmt7 = fmt7,
    db_pass = db_pass
  )
}

# Build the six pRESTO/Change-O commands for the consensus pipeline as
# exec/args pairs for system2. Pure function: it does not run anything.
# `germline_files` is the expanded IMGT germline file list for MakeDb.py -r
# (globs must expand in R because the runner quotes every argument).
.sirup_presto_commands <- function(paths,
                                   receptor,
                                   organism,
                                   identity,
                                   consensus_freq,
                                   germline_files,
                                   igblast_base,
                                   igblast_exec,
                                   muscle_exec,
                                   nproc,
                                   cluster_field,
                                   delimiters,
                                   parse_fields,
                                   parse_key,
                                   parse_act,
                                   consensus_min_count,
                                   consensus_copy_fields,
                                   consensus_copy_action,
                                   consensus_max_error = 0.1,
                                   consensus_max_gap = 0.3,
                                   align_log,
                                   consensus_log) {
  delim_args <- c("--delim", delimiters)
  align_log_args <- if (is.null(align_log)) NULL else c("--log", path.expand(align_log))
  consensus_log_args <- if (is.null(consensus_log)) NULL else c("--log", path.expand(consensus_log))
  consensus_copy_args <- NULL
  if (!is.null(consensus_copy_fields) && length(consensus_copy_fields) > 0L) {
    consensus_copy_args <- c("--cf", consensus_copy_fields, "--act", consensus_copy_action)
  }

  list(
    # Cluster reads per spatial ID, then preserve the spatial_id-cluster
    # grouping in the FASTA headers used by alignment and consensus calling.
    list(exec = "ClusterSets.py", out = paths$cluster, args = c(
      "set",
      "-s", path.expand(paths$input),
      "-f", cluster_field,
      "-k", "cluster",
      "--ident", as.character(identity),
      "--cluster", "vsearch",
      delim_args,
      "--nproc", as.character(nproc),
      "-o", path.expand(paths$cluster)
    )),
    list(exec = "ParseHeaders.py", out = paths$reheader, args = c(
      "merge",
      "-s", path.expand(paths$cluster),
      "-f", parse_fields,
      "-k", parse_key,
      delim_args,
      "--act", parse_act
    )),
    list(exec = "AlignSets.py", out = paths$aligned, args = c(
      "muscle",
      "-s", path.expand(paths$reheader),
      align_log_args,
      delim_args,
      "--nproc", as.character(nproc),
      "--bf", parse_key,
      "--exec", path.expand(muscle_exec)
    )),
    # Build consensus sequences from aligned spatial_id-cluster groups, assign
    # V(D)J genes with IgBLAST, then make an AIRR-style database table.
    list(exec = "BuildConsensus.py", out = paths$consensus, args = c(
      "-s", path.expand(paths$aligned),
      "-o", path.expand(paths$consensus),
      consensus_log_args,
      "--nproc", as.character(nproc),
      "-n", as.character(consensus_min_count),
      "-q", "20",
      "--freq", as.character(consensus_freq),
      "--bf", parse_key,
      "--maxerror", as.character(consensus_max_error),
      "--maxgap", as.character(consensus_max_gap),
      "--fasta",
      delim_args,
      consensus_copy_args
    )),
    list(exec = "AssignGenes.py", out = paths$fmt7, args = c(
      "igblast",
      "-s", path.expand(paths$consensus),
      "-b", path.expand(igblast_base),
      "--organism", organism,
      "--loci", receptor,
      "--format", "blast",
      "--exec", path.expand(igblast_exec),
      "--nproc", as.character(nproc)
    )),
    list(exec = "MakeDb.py", out = paths$db_pass, args = c(
      "igblast",
      "-i", path.expand(paths$fmt7),
      "-s", path.expand(paths$consensus),
      "-r", germline_files,
      "--extended"
    ))
  )
}

# Return the longest gap run in each aligned receptor sequence. Vectorised: the
# spatial tables run to hundreds of thousands of rows, and the scalar form had
# to be driven through vapply() by every caller. Both `.` and `-` count as gap
# characters, so a deletion is measured as well as IMGT padding.
.sirup_max_gap_length <- function(seq) {
  seq <- as.character(seq)
  out <- integer(length(seq))
  keep <- !is.na(seq)
  out[!keep] <- NA_integer_
  if (!any(keep)) {
    return(out)
  }
  runs <- regmatches(seq[keep], gregexpr("[.-]+", seq[keep]))
  out[keep] <- vapply(runs, function(x) if (length(x) == 0L) 0L else max(nchar(x)), integer(1))
  out
}

# Warn when a filter the caller explicitly asked for cannot run because its
# column is absent. Only explicit requests warn: the defaults would otherwise
# fire on every pre-consensus table, where `consensus_count` cannot exist yet.
# Skipping silently has already voided a filter in a real analysis: MakeDb.py
# only writes `consensus_count` when the pRESTO field-value delimiter is `=`, so
# a pipeline run with `:` delimiters made `min_consensus_count` a no-op.
.sirup_warn_missing_filter_col <- function(arg, column) {
  rlang::warn(paste0(
    "`", arg, "` was requested but column `", column,
    "` is absent. The filter did not run."
  ))
}

# Return the length of the leading gap run in an aligned receptor sequence.
.sirup_opening_gap_length <- function(seq) {
  seq <- as.character(seq)
  out <- integer(length(seq))
  keep <- !is.na(seq)
  out[!keep] <- NA_integer_
  if (!any(keep)) {
    return(out)
  }
  m <- regexpr("^[.-]+", seq[keep])
  lens <- attr(m, "match.length")
  lens[m == -1L] <- 0L
  out[keep] <- as.integer(lens)
  out
}
