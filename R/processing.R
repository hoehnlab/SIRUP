#' Filter AIRR sequences
#'
#' `filterAirr` filters sequences from an AIRR table. ]
#'
#' @details
#' The function keeps the rows that pass all of these filters. It skips a filter
#' if its column is missing, unless you set the argument yourself. In that case,
#' it gives a warning.
#' \itemize{
#'   \item the locus (`locus_val`);
#'   \item the V length, J length, and junction length (`min_v_length`,
#'     `min_j_length`, `min_junction_length`);
#'   \item productive rearrangements (`productive_only`);
#'   \item complete V(D)J rearrangements (`complete_vdj`);
#'   \item the consensus count (`min_consensus_count`);
#'   \item the gap lengths in `sequence_alignment` (`max_opening_gap`,
#'     `max_gap_length`).
#' }
#' The function prints the number of rows that it kept.
#'
#' @param airr_df AIRR-style data frame. It needs the columns
#'   `v_sequence_start`, `v_sequence_end`, `j_sequence_start`, `j_sequence_end`,
#'   and `junction_length`.
#' @param locus_val Loci to keep. If `NULL` (default), the function does not
#'   filter by locus.
#' @param min_v_length Smallest V length to keep (default 100).
#' @param min_j_length Smallest J length to keep (default 10).
#' @param min_junction_length Smallest junction length to keep (default 1).
#' @param productive_only Keep only productive rearrangements (default `TRUE`)?
#' @param complete_vdj Keep only rows with a complete V(D)J (default `FALSE`)?
#'   IgBLAST writes this column as `"T"` or `"F"`. The function skips the filter
#'   if the column does not exist.
#' @param min_consensus_count Smallest `consensus_count` to keep (default 2).
#'   Use `NULL` to skip the filter.
#' @param max_opening_gap Longest gap at the start of `sequence_alignment` to
#'   keep. If `NULL` (default), the function does not filter by this gap.
#' @param max_gap_length Longest run of gap characters in `sequence_alignment`
#'   to keep. If `NULL` (default), the function does not filter by this gap. A
#'   real somatic hypermutation indel has a length that is a multiple of three.
#'   It is rarely longer than 12 nucleotides. A longer gap points to a chimeric
#'   or truncated consensus. The function measures the gap in the IMGT-numbered
#'   alignment. In this alignment, the padding depends on the locus: IGHV has
#'   padding of 6, 9, or 12 nucleotides, and TRBV has padding of 18 or 21. A
#'   limit that suits IGH rejects every TRB sequence, so set the limit for each
#'   locus.
#' @param sequence_alignment Name of the alignment column for the gap filters
#'   (default `"sequence_alignment"`).
#' @param id Name of the column that holds the read or consensus IDs (default
#'   `"sequence_id"`).
#' @param spatial_id Name of a spatial ID column that exists in `airr_df`.
#' @param add_spatial_sequence_id Add a column that combines the sequence ID and
#'   the spatial ID (default `FALSE`)? You must give `spatial_id` if you use
#'   `TRUE`.
#' @param spatial_sequence_id Name of the combined column (default
#'   `"spatial_sequence_id"`).
#' @param spatial_id_key Field name for the spatial ID in the combined column
#'   (default `"spatial_id"`). A value looks like
#'   `sequence_id=<id>|<key>=<spatial ID>`.
#'
#' @return The filtered data frame. It has the added columns `v_sequence_length`
#'   and `j_sequence_length`, and the combined ID column if you asked for it.
#' @seealso [addSequenceQualityMetrics()], [plotLengthQc()],
#'   [processSpatialSample()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' receptors <- getReceptors(x)
#'
#' # The thresholds of the SIRUP paper.
#' filterAirr(
#'   receptors,
#'   locus_val = "IGH",
#'   min_v_length = 250,
#'   min_j_length = 30,
#'   min_junction_length = 10
#' )
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

#' Process one spatial sample from reads to an AIRR table
#'
#' `processSpatialSample`processes demultiplexed reads for one sample and
#' one receptor group (`ig` or `tr`). It outputs an AIRR table of filtered
#' consensus sequences, with the V(D)J genes assigned and IMGT-gaps added.
#'
#' @details
#' The function calls external programs: IgBLAST, pRESTO, vsearch, and MUSCLE.
#' It prints each command before it runs it. The function runs all commands with
#' `system2()` and no shell. A non-zero exit status stops the function, and the
#' error message shows the command that failed. The final output is the
#' `db_pass` AIRR TSV, which is in the returned list of paths.
#'
#' The `input` can be one of these:
#' \itemize{
#'   \item Raw reads: a FASTA or FASTQ file, optionally gzipped, from RRAD (or RAD).
#'     Export the rrad FASTA with `to_fasta = TRUE` and
#'     `header_format = "presto"`, so that the barcode and UMI tags stay in the
#'     IgBLAST identifier. The function first annotates the reads with
#'     [runIgblast()] (which calls `airrblast::assignGenes()`). This step needs
#'     `igblast_base` and `igblast_refs`.
#'   \item An AIRR table: a `.tsv` path or a data frame. The function skips the
#'     IgBLAST step. A data frame needs `output_dir`.
#'   \item A FASTA that is ready for pRESTO (`input_type = "spatial_fasta"`).
#'     The function skips the IgBLAST and filter steps, and starts with the
#'     consensus workflow.
#' }
#'
#' For raw reads and AIRR tables, the function first does these steps:
#' \enumerate{
#'   \item It splits the rrad tags `|BARCODE=...|UMI=...` (or the legacy colon
#'     tags) out of `sequence_id` with [parseRadSequenceIds()].
#'   \item It filters sequences with [filterAirr()]. The pipeline
#'     defaults are `complete_vdj = TRUE`, `spatial_id = "barcode"`, and a
#'     combined header ID with the key `cluster_field`. Change them with
#'     `filter_args`.
#'   \item It writes the raw and the filtered AIRR tables to `output_dir`, so
#'     that you can inspect them.
#'   \item It writes the filtered sequences to a spatial FASTA file. The headers
#'     of the file hold the spatial ID.
#' }
#'
#' The spatial FASTA then goes through a number of Immcantation functions:
#' \enumerate{
#'   \item `ClusterSets.py set` groups the reads within each spatial ID. It uses
#'     `cluster_field`, the identity threshold `identity`, and vsearch.
#'   \item `ParseHeaders.py merge` writes a new header field, `parse_key`.
#'     The field combines the spatial ID and the cluster ID (`parse_fields`)
#'     with the action `parse_act`. `delimiters` sets the delimiters of the
#'     FASTA header.
#'   \item `AlignSets.py muscle` aligns the reads within each `parse_key` group.
#'     `muscle_exec` is the path of the MUSCLE program
#'     (<https://github.com/rcedgar/muscle>). The function can write a log to
#'     `align_log`.
#'   \item `BuildConsensus.py` makes a consensus sequence for each
#'     `parse_key` group. `consensus_min_count` is the smallest number of
#'     sequences that give a valid consensus. `consensus_freq` is the fraction
#'     of occurrences below which the function assigns an ambiguous character.
#'     The function can write a log to `consensus_log`. To copy more header
#'     fields to the consensus file, give their names in
#'     `consensus_copy_fields` and the way to copy them in
#'     `consensus_copy_action` (`majority`, `first`, `last`, or `all`).
#'   \item `airrblast::assignGenes()` assigns the V(D)J genes with IgBLAST. It
#'     replaces `AssignGenes.py`. `receptor` is `ig` (BCR) or `tr` (TCR).
#'     `organism` is the species (human, mouse, or rhesus_monkey).
#'     `igblast_base` is the IgBLAST installation, which must hold
#'     `bin/igblastn`. `igblast_refs` is the folder with the IgBLAST `database/`
#'     and `fasta/` folders.
#'   \item `airrblast::addGaps()` adds the IMGT-gapped sequences from the
#'     germlines in `imgt_dir`. SIRUP then reads the pRESTO header fields into
#'     columns, as `MakeDb.py --extended` did. It keeps the sequences that have
#'     a V call, a J call, a productive flag, an aligned sequence, and a
#'     junction. The result is the `db_pass` table.
#' }
#'
#' @param input Raw reads (FASTA or FASTQ, optionally gzipped), the path of an
#'   AIRR TSV, an AIRR data frame, or a FASTA that is ready for pRESTO.
#' @param output_dir Directory for all outputs. The default is the directory of
#'   the input file. You must give it if `input` is a data frame.
#' @param receptor `"ig"` for BCR loci or `"tr"` for TCR loci.
#' @param organism Organism of the sample (default `"human"`).
#' @param input_type How to treat `input`: `"auto"` (default), `"raw_reads"`,
#'   `"airr"`, or `"spatial_fasta"`. With `"auto"`, the function finds the type
#'   from the class of the object and the file extension. A FASTA or FASTQ path
#'   becomes `"raw_reads"`. Use `"spatial_fasta"` to send a FASTA that is ready
#'   for pRESTO straight to the consensus workflow. For a new rrad FASTA, use
#'   `"auto"` or `"raw_reads"`, because `"spatial_fasta"` skips the first
#'   annotation of V(D)J genes and the filtering of the reads.
#' @param filter_args Named list that changes arguments of [filterAirr()]. The
#'   pipeline defaults are `complete_vdj = TRUE`, `min_consensus_count = NULL`,
#'   `spatial_id = "barcode"`, `add_spatial_sequence_id = TRUE`, and
#'   `spatial_id_key = cluster_field`. A `NULL` value turns off an optional
#'   filter (`locus_val`, `min_consensus_count`, `max_opening_gap`, or
#'   `spatial_id`). The other arguments do not accept `NULL`.
#' @param invalid_spatial_ids What to do with filtered rows whose spatial
#'   barcode is missing or holds a delimiter: `"abort"` (default) stops with a
#'   count, and `"drop"` removes the rows with a warning.
#' @param fasta_sequence_col Name of the AIRR column to write to the spatial
#'   FASTA (default `"sequence_alignment"`, the aligned V(D)J region). This
#'   keeps extra read sequence out of the clustering and consensus steps. Use
#'   `"sequence"` to write the full read.
#' @param sample_name Stem of the intermediate file names. The default is the
#'   input file name without the extension, or `"sample"` for a data frame.
#' @param identity Identity threshold for the clusters, passed to
#'   `ClusterSets.py` (default 0.975).
#' @param consensus_freq Frequency threshold for the consensus, passed to
#'   `BuildConsensus.py` (default 0.51).
#' @param germline_glob Pattern of the germline files, such as
#'   `<imgt_dir>/human/vdj/imgt_human_*.fasta`. The function uses it only to
#'   find `imgt_dir` if `imgt_dir` is `NULL`.
#' @param igblast_base IgBLAST installation folder (default
#'   `"~/software/igblast"`). The function passes it to
#'   `airrblast::assignGenes()` as `igblast`. The folder must hold
#'   `bin/igblastn`.
#' @param muscle_exec Path of the MUSCLE program, passed to `AlignSets.py
#'   --exec` (default `"~/software/muscle3.8.31_i86linux64"`).
#' @param nproc Number of processors for the IgBLAST step on raw reads, for
#'   pRESTO, and for the gene assignment (default 8).
#' @param cluster_field FASTA header field for the clustering within a spot
#'   (default `"barcode"`, the parsed rrad `BARCODE` field). The pipeline writes
#'   the headers with this field name. Change it for other spatial ID fields,
#'   such as `barcode_8um`.
#' @param delimiters Character vector of the pRESTO delimiters, passed after
#'   `--delim` (default `c("|", "=", ",")`).
#' @param parse_fields Header fields for `ParseHeaders.py merge -f` (default
#'   `c(cluster_field, "cluster")`).
#' @param parse_key Name of the header field that `ParseHeaders.py merge -k`
#'   makes (default `paste0(cluster_field, "_cluster")`).
#' @param parse_act Action for `ParseHeaders.py merge --act` (default `"set"`).
#'   With `"set"`, the function joins the spatial ID and the cluster ID with a
#'   comma, so that you can recover the spatial ID from the consensus header.
#'   With `"cat"`, it joins them with no separator, so that you cannot recover
#'   the spatial ID of a barcode that has a variable length. Both give the same
#'   groups, so use `"set"`.
#' @param consensus_min_count Smallest number of sequences for a consensus,
#'   passed to `BuildConsensus.py -n` (default 1). Use at least 2 for data that
#'   do not have single-cell resolution. This value decides which groups the
#'   function writes. It does not change the consensus of a group that passes.
#'   So you can run with `1` and filter to a higher value later with
#'   `filterAirr(min_consensus_count = )`. `BuildConsensus.py -q` is fixed at 20
#'   and has no effect here, because the pipeline input is FASTA and has no base
#'   qualities.
#' @param umi_support Group the reads by spatial barcode and UMI before the
#'   barcode step (default `FALSE`)? Reads of one molecule that differ by a
#'   sequencing error end up in different clusters at `identity = 1.0`, and the
#'   function then drops them as singletons. This stage first joins them at
#'   `umi_identity`, so that they support each other. Each barcode-UMI cluster
#'   gives one consensus. The consensus carries `cc_bc_umi` (the number of reads
#'   in the cluster, counted before the productive subset) and
#'   `n_productive_reads` (the number of reads that the consensus uses). The
#'   function discards a cluster that has no such reads. The barcode step then
#'   runs with `BuildConsensus -n 1` and sums both fields, and `min_umi_support`
#'   replaces the read minimum. With `umi_weight = TRUE`, the function writes
#'   each consensus once for each productive read, so that one read is one vote
#'   in the barcode consensus. A short read from one molecule then cannot trim
#'   the consensus of many full-length reads. For this reason, `consensus_count`
#'   counts reads and not molecules. The option needs an AIRR table with the
#'   columns `umi` and `productive`.
#' @param umi_weight Write each molecule consensus once for each productive read
#'   (default `TRUE`)? The function uses it only with `umi_support`. If `FALSE`,
#'   it writes each molecule once, so that one UMI cluster is one vote.
#' @param umi_field Name of the column that holds the UMI (default `"umi"`). The
#'   function uses it only with `umi_support`.
#' @param umi_identity Identity threshold for the barcode-UMI clusters (default
#'   0.95). It is lower than `identity` on purpose, so that reads of one
#'   molecule with errors form one group.
#' @param min_umi_support Smallest summed `cc_bc_umi` that a barcode consensus
#'   needs (default 2). The function applies it after the gene assignment. It
#'   writes a filtered copy next to `db_pass`, so that you can change the
#'   threshold without a new run of the pipeline.
#' @param consensus_copy_fields Header fields for `BuildConsensus.py --cf`. With
#'   `umi_support`, the function adds `CC_BC_UMI` and `N_PRODUCTIVE_READS` at
#'   the start of the list and sums them.
#' @param consensus_copy_action Action for `BuildConsensus.py --act`, if you
#'   give `consensus_copy_fields` (default `"majority"`).
#' @param consensus_max_error Largest error rate, passed to
#'   `BuildConsensus.py --maxerror` (default 0.1).
#' @param consensus_max_gap Largest gap fraction, passed to
#'   `BuildConsensus.py --maxgap` (default 0.3).
#' @param igblast_refs Folder with the IgBLAST `database/` and `fasta/` folders
#'   (the `refs` argument of `airrblast::assignGenes()`). If `NULL`, the
#'   function uses `igblast_base`.
#' @param imgt_dir IMGT folder with the gapped germlines, at
#'   `<imgt_dir>/<organism>/vdj` (the `gapdb` argument of
#'   `airrblast::addGaps()`). If `NULL`, the function finds it from
#'   `germline_glob`, which must point into a `<organism>/vdj` folder.
#' @param skip_existing Skip a pipeline command if its output already exists and
#'   is not empty (default `FALSE`)? Use it to restart an interrupted run. You
#'   can also use it to reuse one clustering and alignment pass for several
#'   `consensus_freq` values, because those steps depend only on `identity`.
#' @param align_log Path of a log file for `AlignSets.py --log`.
#' @param consensus_log Path of a log file for `BuildConsensus.py --log`.
#'
#' @return A named list, invisibly. It holds `raw_airr`, `filtered_airr`,
#'   `spatial_fasta`, `n_raw`, and `n_filtered`, which are `NULL` for the steps
#'   that did not run. It also holds the output paths of pRESTO and IgBLAST.
#'   `db_pass` is the final AIRR TSV. With `umi_support`, `db_pass_umi` is the
#'   same table, filtered to `min_umi_support`.
#' @seealso [processSpatialSamples()], [runIgblast()], [filterAirr()],
#'   [parseRadSequenceIds()]
#'
#' @examples
#' \dontrun{
#' # This example needs IgBLAST, pRESTO, vsearch, and MUSCLE. See
#' # vignette("processing", package = "SIRUP").
#' result <- processSpatialSample(
#'   "sample1_reads.fasta.gz",
#'   output_dir = "sample1_output",
#'   receptor = "ig",
#'   organism = "human",
#'   igblast_base = "~/software/igblast",
#'   imgt_dir = "~/share/germlines/imgt",
#'   muscle_exec = "~/software/muscle",
#'   nproc = 8
#' )
#' bcrs <- airr::read_rearrangement(result$db_pass)
#' }
#' @export
processSpatialSample <- function(input,
                                 output_dir = NULL,
                                 receptor = c("ig", "tr"),
                                 organism = "human",
                                 input_type = c("auto", "raw_reads", "airr", "spatial_fasta"),
                                 filter_args = list(),
                                 invalid_spatial_ids = c("abort", "drop"),
                                 fasta_sequence_col = "sequence_alignment",
                                 sample_name = NULL,
                                 identity = 0.975,
                                 consensus_freq = 0.51,
                                 germline_glob = NULL,
                                 igblast_refs = NULL,
                                 imgt_dir = NULL,
                                 igblast_base = "~/software/igblast",
                                 muscle_exec = "~/software/muscle3.8.31_i86linux64",
                                 nproc = 8,
                                 cluster_field = "barcode",
                                 delimiters = c("|", "=", ","),
                                 parse_fields = c(cluster_field, "cluster"),
                                 parse_key = paste0(cluster_field, "_cluster"),
                                 parse_act = "set",
                                 consensus_min_count = 1,
                                 umi_support = FALSE,
                                 umi_weight = TRUE,
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
  if (is.null(imgt_dir)) {
    # addGaps() reads <imgt_dir>/<organism>/vdj, so a glob such as
    # "/x/imgt/human/vdj/imgt_human_*.fasta" gives "/x/imgt".
    if (is.null(germline_glob)) {
      rlang::abort("Give `imgt_dir` or `germline_glob`.")
    }
    vdj <- dirname(path.expand(germline_glob))
    if (!identical(basename(vdj), "vdj") || !identical(basename(dirname(vdj)), organism)) {
      rlang::abort(paste0(
        "Cannot infer `imgt_dir` from `germline_glob` (", germline_glob, "). ",
        "The glob must point into `<imgt_dir>/", organism, "/vdj`. ",
        "Set `imgt_dir` instead."
      ))
    }
    imgt_dir <- dirname(dirname(vdj))
  }
  if (is.null(igblast_refs)) {
    igblast_refs <- igblast_base
  }
  # assignGenes() always runs <igblast>/bin/igblastn.
  igblast_exec <- file.path(igblast_base, "bin", "igblastn")

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
        igblast_base = igblast_base,
        igblast_refs = igblast_refs,
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
        umi_weight = umi_weight,
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

  if (skip_existing && file.exists(paths$db_pass) && file.info(paths$db_pass)$size > 0) {
    message("Skipping airrblast: output exists (", paths$db_pass, ").")
  } else {
    # Gene assignment. assignGenes() runs IgBLAST and replaces AssignGenes.py.
    # addGaps() adds the IMGT-gapped sequences. With the header and filter steps
    # below, it replaces MakeDb.py.
    locus <- c(ig = "Ig", tr = "TR")[[receptor]]
    assigned <- airrblast::assignGenes(
      file = path.expand(paths$consensus),
      igblast = path.expand(igblast_base),
      refs = path.expand(igblast_refs),
      organism = organism,
      outfile = path.expand(paths$igblast_tsv),
      nproc = nproc,
      locus = locus,
      return = TRUE,
      verbose = FALSE
    )
    if (nrow(assigned) == 0L) {
      rlang::abort(paste0("IgBLAST returned no rows for ", paths$consensus, "."))
    }
    db <- as.data.frame(
      airrblast::addGaps(
        assigned,
        gapdb = path.expand(imgt_dir),
        organism = organism,
        locus = locus,
        nproc = nproc
      ),
      stringsAsFactors = FALSE
    )

    # Read the pRESTO header fields, as MakeDb.py --extended did. The text before
    # the first `|` becomes `sequence_id`. Each `KEY=value` field after it becomes
    # a lower-case column and replaces any column of the same name that IgBLAST
    # wrote. CONSCOUNT and DUPCOUNT become `consensus_count` and `duplicate_count`,
    # as in Change-O. This uses the pRESTO default delimiters.
    header <- as.character(db$sequence_id)
    if (any(grepl("|", header, fixed = TRUE))) {
      parts <- strsplit(header, "|", fixed = TRUE)
      fields <- unlist(lapply(parts, function(x) x[-1L]), use.names = FALSE)
      if (any(!grepl("=", fields, fixed = TRUE))) {
        # MakeDb.py gives up on the whole file here and keeps the headers as they are.
        rlang::warn("Some consensus headers have fields without `=`. Headers are left unparsed.")
      } else {
        db$sequence_id <- vapply(parts, function(x) if (length(x)) x[[1L]] else "", character(1))
        header_row <- rep(seq_along(parts), lengths(parts) - 1L)
        header_key <- tolower(sub("=.*$", "", fields))
        header_value <- sub("^[^=]*=([^=]*).*$", "\\1", fields)
        header_key[header_key == "conscount"] <- "consensus_count"
        header_key[header_key == "dupcount"] <- "duplicate_count"
        for (k in unique(header_key)) {
          hit <- header_key == k
          col <- if (k %in% names(db)) as.character(db[[k]]) else rep(NA_character_, nrow(db))
          col[header_row[hit]] <- header_value[hit]
          db[[k]] <- col
        }
      }
    }

    # Keep the sequences that have a V call, a J call, a productive flag, an
    # aligned sequence, and a junction. This is MakeDb.py's check without the
    # test that the junction lies inside the gapped sequence.
    keep <- !is.na(db$productive)
    for (col in c("v_call", "j_call", "sequence_alignment", "junction")) {
      keep <- keep & !is.na(db[[col]]) & nzchar(as.character(db[[col]])) & as.character(db[[col]]) != "None"
    }
    message("airrblast: ", sum(keep), " of ", nrow(db), " sequences passed -> ", paths$db_pass)
    db <- db[keep, , drop = FALSE]

    # MakeDb.py wrote logical columns as T/F.
    for (col in names(db)[vapply(db, is.logical, logical(1))]) {
      db[[col]] <- ifelse(is.na(db[[col]]), NA_character_, ifelse(db[[col]], "T", "F"))
    }
    data.table::fwrite(db, paths$db_pass, sep = "\t", quote = FALSE, na = "")
  }

  db_pass_umi <- NULL
  if (isTRUE(umi_support) && file.exists(paths$db_pass)) {
    db <- .sirup_read_airr(paths$db_pass)
    sup_col <- intersect(c("cc_bc_umi", "CC_BC_UMI"), colnames(db))
    if (length(sup_col) == 0L) {
      rlang::warn("The gene-assignment table has no `cc_bc_umi` column, so the `min_umi_support` floor was not applied.")
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

#' Process several spatial samples
#'
#' `processSpatialSamples` runs [processSpatialSample()] on each input in turn.
#'
#' @param inputs List or character vector of inputs, one for each sample. Each
#'   input goes to [processSpatialSample()]. It can be a reads file, the path of
#'   an AIRR TSV, or an AIRR data frame.
#' @param output_dirs Vector of output directories, one for each input. If
#'   `NULL`, each sample uses the default of [processSpatialSample()].
#' @param sample_names Vector of sample names, one for each input. If `NULL`,
#'   the function uses the file name of each input, or `sample_<i>` for a data
#'   frame.
#' @param ... Other arguments for [processSpatialSample()].
#'
#' @return A list with one list of output paths for each sample, invisibly.
#' @seealso [processSpatialSample()]
#'
#' @examples
#' \dontrun{
#' # This example needs IgBLAST, pRESTO, vsearch, and MUSCLE.
#' results <- processSpatialSamples(
#'   c("sample1_reads.fasta.gz", "sample2_reads.fasta.gz"),
#'   output_dirs = c("sample1_output", "sample2_output"),
#'   receptor = "ig",
#'   imgt_dir = "~/share/germlines/imgt"
#' )
#' }
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

#' Keep the most abundant consensus sequence of each spatial ID and locus
#'
#' `collapseTopConsensusPerSpatialId` keeps the row with the highest consensus
#' count for each spatial ID and locus.
#'
#' @details
#' The function adds the column `prop_total_consensus`. It then adds the quality
#' columns from [addSequenceQualityMetrics()]: `non_N_nucleotide_prop`,
#' `max_gap_length`, and `opening_gap_length`. If several rows tie for the top
#' consensus count, the function prints how many groups have a tie.
#'
#' @param airr_df AIRR-style data frame.
#' @param spatial_id Name of the spatial ID column for the grouping (default
#'   `"spatial_id_2um"`).
#' @param locus Name of the locus column for the grouping (default `"locus"`).
#' @param consensus Name of the consensus count column for the ranking (default
#'   `"consensus_count"`).
#' @param sequence Name of the sequence column for the proportion of called
#'   bases (default `"sequence"`).
#' @param alignment Name of the alignment column for the gap measures (default
#'   `"sequence_alignment"`).
#' @param ties Keep `"all"` rows that tie for the top count (default), or only
#'   the `"first"` of them.
#'
#' @return The collapsed data frame. It has the column `prop_total_consensus`
#'   (the count of the row divided by the total count of its group) and the
#'   quality columns.
#' @seealso [addSequenceQualityMetrics()], [filterAirr()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' collapseTopConsensusPerSpatialId(getReceptors(x), spatial_id = "spatial_id")
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

#' Add quality measures to an AIRR table
#'
#' `addSequenceQualityMetrics` adds columns that can help indicate the quality of each
#' consensus sequence. It keeps all rows.
#'
#' @details
#' Use this function to look at the distributions of the measures before you
#' choose the thresholds for [filterAirr()].
#' [collapseTopConsensusPerSpatialId()] also adds these columns, but it removes
#' rows.
#'
#' The function adds these columns, if the source column exists:
#' \itemize{
#'   \item `non_N_nucleotide_prop`: the share of called bases (`A`, `C`, `G`,
#'     `T`) in `sequence`. A low value means an ambiguous consensus. The
#'     function takes this value from the ungapped sequence. The alignment
#'     would put IMGT gap characters in the denominator.
#'   \item `max_gap_length`: the longest run of gap characters in `alignment`.
#'   \item `opening_gap_length`: the gap run at the start of `alignment`.
#' }
#'
#' @param airr_df AIRR-style data frame.
#' @param sequence Name of the sequence column for the share of called bases
#'   (default `"sequence"`).
#' @param alignment Name of the alignment column for the gap measures (default
#'   `"sequence_alignment"`).
#'
#' @return `airr_df` with the quality columns added.
#' @seealso [filterAirr()], [collapseTopConsensusPerSpatialId()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' metrics <- addSequenceQualityMetrics(getReceptors(x))
#' head(metrics[, c("non_N_nucleotide_prop", "max_gap_length", "opening_gap_length")])
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

#' Count the sequences and spots of a receptor table
#'
#' `spatialReceptorSummary` returns the basic counts of a receptor table.
#'
#' @param airr_df AIRR-style data frame, or a SIRUP object.
#' @param spatial_id Name of the spatial ID column to count (default
#'   `"spatial_id_2um"`). For a SIRUP object, if you omit it, the function uses
#'   the stored spatial ID field. If the object has a sample field, the function
#'   counts each spatial ID once for each sample.
#'
#' @return A one-row tibble with the columns `n_sequences`, `n_spatial_ids`,
#'   `n_loci`, and `n_productive`. `n_loci` and `n_productive` are `NA` if the
#'   table has no `locus` or `productive` column.
#' @seealso [filterAirr()], [spatialCloneDominance()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' spatialReceptorSummary(x)
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

#' Collapse neighboring spots with identical receptor sequences
#'
#' `collapseIdenticalNeighbors` finds groups of neighboring spots that share a
#' receptor sequence. It keeps one representative spot for each group.
#'
#' @details
#' For each unique receptor sequence, the function finds the connected groups of
#' neighboring spots that hold that sequence. It keeps the spot with the highest
#' consensus count in each group. This hopefully reduces bleedover, but in a very 
#' coarse manner.
#'
#' @param object SIRUP object, or an AIRR-style receptor data frame.
#' @param seurat_object Seurat object. Use it only if `object` is a receptor
#'   data frame.
#' @param spatial_id Name of the receptor column that holds the Seurat spot
#'   names. For a SIRUP object, if you omit it, the function uses the stored
#'   spatial ID field. For a receptor table, if you omit it, the function looks
#'   for a standard spatial ID column.
#' @param sequence Name of the column that holds the receptor sequences (default
#'   `"sequence"`).
#' @param consensus_count Name of the column that holds the consensus counts
#'   (default `"consensus_count"`). The function uses it to choose the
#'   representative spot of a group.
#' @param images Names of the Seurat images to use for the coordinates.
#' @param image_scale Image scale for `SeuratObject::GetTissueCoordinates()`.
#' @param swap_coords Swap the x and y coordinates after the function reads the
#'   Seurat coordinates?
#' @param num_neighbors_non_boundary Number of nearby spots that makes a spot a
#'   non-boundary spot (default 6).
#' @param validate Validate the updated SIRUP object (default `TRUE`)? The
#'   function uses it only if `object` is a SIRUP object.
#'
#' @return For a SIRUP object, the updated SIRUP object with a new receptor
#'   table. For a receptor table, the collapsed table. Both keep only the
#'   representative spot of each group of identical sequences.
#' @seealso [filterAirr()], [updateReceptors()]
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
    spatial_samples <- unique(stats::na.omit(as.character(tissue_positions$sample)))
    receptor_samples <- as.character(data[[sample_col]])
    samples <- intersect(
      unique(stats::na.omit(receptor_samples)),
      spatial_samples
    )
    unmatched_samples <- data[is.na(receptor_samples) | !receptor_samples %in% spatial_samples, , drop = FALSE]
    if (length(samples) == 0L) {
      return(unmatched_samples)
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
    return(dplyr::bind_rows(c(out, list(unmatched_samples))))
  }

  original_cols <- colnames(data)
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
  receptor_key <- .sirup_key_vector(.sirup_collapse_key_table(data, temp))
  tissue_key <- .sirup_key_vector(.sirup_collapse_key_table(tissue_positions, temp))
  matched <- receptor_key %in% tissue_key
  unmatched <- data[!matched, original_cols, drop = FALSE]
  data <- data[matched, , drop = FALSE]
  data <- merge(
    data,
    tissue_positions,
    by = merge_by,
    all.x = FALSE,
    all.y = FALSE,
    sort = FALSE
  )
  if (nrow(data) == 0L) {
    return(unmatched)
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
    result <- dplyr::bind_rows(result, unmatched)
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
  result <- dplyr::bind_rows(result, unmatched)
  n_result_pretty <- scales::comma(nrow(result))
  n_data_pretty <- scales::comma(nrow(data) + nrow(unmatched))
  n_removed_pretty <- scales::comma(nrow(data) + nrow(unmatched) - nrow(result))

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
  igblast_tsv <- sub("\\.fasta$", "_igblast.tsv", consensus)
  db_pass <- sub("\\.fasta$", "_igblast_db-pass.tsv", consensus)
  list(
    input = input_fasta,
    cluster = cluster,
    reheader = reheader,
    aligned = aligned,
    consensus = consensus,
    igblast_tsv = igblast_tsv,
    db_pass = db_pass
  )
}

# Build the pRESTO commands for the consensus pipeline as
# exec/args pairs for system2. Pure function: it does not run anything.
.sirup_presto_commands <- function(paths,
                                   receptor,
                                   organism,
                                   identity,
                                   consensus_freq,
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

  commands <- list(
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
    # Build consensus sequences from aligned spatial_id-cluster groups.
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
    ))
  )

  # Gene assignment is not a shell command: processSpatialSample() runs
  # airrblast::assignGenes() and addGaps() after these.
  commands
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
.datatable.aware = TRUE

#' Write sequences to a FASTQ or FASTA file
#'
#' `writeFastqas` writes a table of sequences to a FASTQ or FASTA file.
#'
#' @param df Data frame or data.table with a `seq` column. For FASTQ, it also
#'   needs a `qual` column. The `header_col` column holds the record headers.
#' @param fn Path of the output file.
#' @param type File type: `"fq"` (default), `"fastq"`, `"fa"`, or `"fasta"`.
#' @param header_col Name of the column that holds the record headers (default
#'   `"id"`). Use `"full_id"` to write the full header.
#' @param append Add the records to an existing file (default `FALSE`)?
#' @param nthreads Number of threads for `data.table::fwrite()` (default 1).
#'
#' @return `fn`, invisibly.
#' @seealso [processSpatialSample()]
#'
#' @examples
#' reads <- data.frame(id = c("read1", "read2"), seq = c("ACGTACGT", "TTGACCGA"))
#' path <- tempfile(fileext = ".fasta")
#' writeFastqas(reads, path, type = "fasta")
#' readLines(path)
#' @export
writeFastqas <- function(df, fn, type = "fq", header_col = "id",
                          append = FALSE, nthreads = 1) {
  type <- .fastqas_normalize_type(type)
  if (!header_col %in% names(df)) {
    stop("header_col '", header_col, "' not found in df. ",
         "Use header_col = 'id' or 'full_id'.")
  }
  if (!"seq" %in% names(df)) stop("`df` must have a `seq` column.")
  prefix  <- if (type == "fq") "@" else ">"
  headers <- paste0(prefix, sub("^[@>]", "", df[[header_col]]))
  if (type == "fq") {
    if (!"qual" %in% names(df)) stop("FASTQ output requires a `qual` column.")
    out <- tibble::tibble(header = headers,
                          seq    = df$seq,
                          plus   = "+",
                          qual   = df$qual)
  } else {
    out <- tibble::tibble(header = headers,
                          seq    = df$seq)
  }
  flat <- as.character(t(as.matrix(out)))
  data.table::fwrite(
    data.table::data.table(line = flat),
    file = fn, compress = "auto", col.names = FALSE,
    quote = FALSE, append = append, nThread = nthreads
  )
  invisible(fn)
}


# Shared FASTQ/FASTA header contract for writeFastqas:
#   - `id`      always present; short id (everything before first space),
#               with no leading `@` (fastq) or `>` (fasta).
#   - `full_id` present only when `full_id = TRUE`; full header (name +
#               comment) with no leading `@`/`>`.
#   - `seq`     always present.
#   - `qual`    present for fastq.
#   - Column order: id, full_id (if requested), seq, qual (if fastq).
#
# `type` accepts "fq" / "fastq" / "fa" / "fasta".

.fastqas_normalize_type <- function(type) {
  type <- match.arg(type, c("fq", "fastq", "fa", "fasta"))
  if (type %in% c("fa", "fasta")) "fa" else "fq"
}
# Helpers and exported functions for the raw-read IgBLAST step of the
# spatial processing pipeline. Prepared reads come from the rrad package as
# FASTA/FASTQ files, optionally gzipped, with headers like
# `id|BARCODE=YYY|UMI=XXX`.

# Check that each external tool exists. `execs` is a named character vector;
# names are used in the error message. A tool passes if its expanded path is
# an existing file or if it resolves on the PATH.
.sirup_check_tools <- function(execs) {
  found <- vapply(execs, function(x) {
    file.exists(path.expand(x)) || nzchar(Sys.which(x))
  }, logical(1))
  if (any(!found)) {
    missing <- execs[!found]
    rlang::abort(paste0(
      "External tool", if (length(missing) == 1L) "" else "s", " not found: ",
      paste0(names(missing), " (", missing, ")", collapse = ", "),
      ". Install the tool or pass its full path."
    ))
  }
  invisible(TRUE)
}

# Run one external command with system2. Every argument is shell-quoted here,
# so callers pass plain strings. `log` redirects stdout and stderr to one
# file. `stdin`/`stdout` take file paths for redirection. Aborts on a
# non-zero exit status and includes the log tail when a log exists.
.sirup_run_command <- function(exec, args, log = NULL, stdin = "", stdout = "") {
  exec_path <- path.expand(exec)
  if (!file.exists(exec_path)) {
    exec_path <- exec
  }
  args <- as.character(args)
  shown <- paste(c(exec, args), collapse = " ")
  message("$ ", shown)
  if (!is.null(log)) {
    stdout <- log
    stderr <- log
  } else {
    stderr <- ""
  }
  status <- system2(exec_path, shQuote(args), stdout = stdout, stderr = stderr, stdin = stdin)
  if (!identical(status, 0L)) {
    detail <- ""
    if (!is.null(log) && file.exists(log)) {
      tail_lines <- utils::tail(readLines(log, warn = FALSE), 10L)
      detail <- paste0("\nLog tail (", log, "):\n", paste(tail_lines, collapse = "\n"))
    }
    rlang::abort(paste0("Command failed with status ", status, ": ", shown, detail))
  }
  invisible(status)
}

# Build a regex character class that matches any pRESTO header delimiter or
# whitespace. Values that contain one of these corrupt FASTA header parsing,
# so validation must use the same `delimiters` the pipeline passes to pRESTO.
.sirup_delims_regex_class <- function(delimiters) {
  # POSIX classes take no backslash escapes: `]` is literal only directly
  # after the opening `[`, and `-` only at the end. `^` cannot come first
  # because `[:space:]` always precedes the delimiter characters.
  delimiters <- unique(delimiters)
  lead <- if ("]" %in% delimiters) "]" else ""
  trail <- if ("-" %in% delimiters) "-" else ""
  middle <- setdiff(delimiters, c("]", "-"))
  paste0("[", lead, "[:space:]", paste(middle, collapse = ""), trail, "]")
}

# Classify a pipeline input. Returns "airr_df" for a data frame, "airr_tsv"
# for a tabular file path, or "raw_reads" for a FASTA/FASTQ path (optionally
# gzipped). Aborts for anything else.
.sirup_detect_input_type <- function(input) {
  if (is.data.frame(input)) {
    return("airr_df")
  }
  if (!is.character(input) || length(input) != 1L || is.na(input)) {
    rlang::abort("`input` must be a data frame or a single file path.")
  }
  base <- sub("\\.gz$", "", input, ignore.case = TRUE)
  ext <- tolower(tools::file_ext(base))
  if (ext %in% c("tsv", "tab", "txt")) {
    return("airr_tsv")
  }
  if (ext %in% c("fa", "fasta", "fna", "fq", "fastq")) {
    return("raw_reads")
  }
  rlang::abort(paste0(
    "Cannot detect the input type from extension '.", ext, "'. ",
    "Accepted: .tsv/.tab/.txt for AIRR tables, ",
    ".fa/.fasta/.fna/.fq/.fastq for reads (each optionally .gz)."
  ))
}

# Read an AIRR TSV into a data frame. fread handles gzipped files natively.
.sirup_read_airr <- function(path) {
  if (!file.exists(path.expand(path))) {
    rlang::abort(paste0("AIRR file not found: ", path))
  }
  out <- data.table::fread(path.expand(path), sep = "\t", header = TRUE, data.table = FALSE)
  out
}

# Normalize raw reads to a plain FASTA that igblastn can read. igblastn
# cannot read gzip and only accepts FASTA, so this streams on disk with
# gunzip and awk; it never loads sequences into R. A plain FASTA passes
# through unchanged. FASTQ quality strings are dropped.
.sirup_prepare_query_fasta <- function(input, output_dir, sample_name) {
  input <- path.expand(input)
  if (!file.exists(input)) {
    rlang::abort(paste0("Input reads file not found: ", input))
  }
  is_gz <- grepl("\\.gz$", input, ignore.case = TRUE)
  base <- sub("\\.gz$", "", input, ignore.case = TRUE)
  ext <- tolower(tools::file_ext(base))
  is_fastq <- ext %in% c("fq", "fastq")
  if (!is_gz && !is_fastq) {
    return(input)
  }
  out_fasta <- file.path(output_dir, paste0(sample_name, "_input.fasta"))
  if (!is_fastq) {
    # Gzipped FASTA: stream-decompress only.
    .sirup_run_command("gunzip", c("-c", input), stdout = out_fasta)
    return(out_fasta)
  }
  fastq <- input
  tmp_fastq <- NULL
  if (is_gz) {
    tmp_fastq <- file.path(output_dir, paste0(sample_name, "_input.fastq"))
    .sirup_run_command("gunzip", c("-c", input), stdout = tmp_fastq)
    fastq <- tmp_fastq
  }
  # FASTQ to FASTA in one streaming pass; qualities are dropped.
  .sirup_run_command(
    "awk",
    "NR % 4 == 1 { sub(/^@/, \">\"); print } NR % 4 == 2 { print }",
    stdin = fastq,
    stdout = out_fasta
  )
  if (!is.null(tmp_fastq)) {
    unlink(tmp_fastq)
  }
  out_fasta
}

#' Annotate raw reads with IgBLAST
#'
#' `runIgblast` runs `airrblast::assignGenes()` on raw reads. It writes an AIRR
#' table.
#'
#' @details
#' The input can be FASTA or FASTQ, optionally gzipped. The function first
#' writes gzipped and FASTQ inputs as a plain FASTA next to `output_tsv`. It
#' drops the FASTQ quality strings, because `igblastn` accepts only FASTA. Each
#' call runs one receptor group: `ig` or `tr`.
#'
#' `assignGenes()` does not search the C-region database, so the tables of raw
#' reads have no `c_call` from this step.
#'
#' @param input Path of a FASTA or FASTQ file, optionally gzipped. A plain FASTA
#'   file must end in `.fasta` or `.fa`.
#' @param output_tsv Path of the AIRR TSV file. It must end in `.tsv`.
#' @param receptor `"ig"` for BCR or `"tr"` for TCR.
#' @param organism `"human"`, `"mouse"`, or `"rhesus_monkey"`.
#' @param igblast_base IgBLAST installation folder. It must hold `bin/igblastn`,
#'   `internal_data`, and `optional_file`.
#' @param igblast_refs Folder with the IgBLAST `database/` and `fasta/` folders
#'   (the `refs` argument of `airrblast::assignGenes()`). If `NULL`, the
#'   function uses `igblast_base`.
#' @param num_threads Number of threads, passed to `assignGenes()` as `nproc`.
#'
#' @return `output_tsv`, invisibly.
#' @seealso [processSpatialSample()], [parseRadSequenceIds()]
#'
#' @examples
#' \dontrun{
#' # This example needs IgBLAST and the IgBLAST reference files.
#' runIgblast(
#'   "sample1_reads.fasta.gz",
#'   output_tsv = "sample1_igblast.tsv",
#'   receptor = "ig",
#'   organism = "human",
#'   igblast_base = "~/software/igblast"
#' )
#' }
#' @export
runIgblast <- function(input,
                       output_tsv,
                       receptor = c("ig", "tr"),
                       organism = "human",
                       igblast_base = "~/software/igblast",
                       igblast_refs = NULL,
                       num_threads = 8) {
  receptor <- match.arg(receptor)
  if (is.null(igblast_refs)) {
    igblast_refs <- igblast_base
  }
  output_tsv <- path.expand(output_tsv)
  output_dir <- dirname(output_tsv)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  sample_name <- tools::file_path_sans_ext(basename(output_tsv))
  query_fasta <- .sirup_prepare_query_fasta(input, output_dir, sample_name)
  airrblast::assignGenes(
    file = query_fasta,
    igblast = path.expand(igblast_base),
    refs = path.expand(igblast_refs),
    organism = organism,
    outfile = output_tsv,
    nproc = num_threads,
    locus = c(ig = "Ig", tr = "TR")[[receptor]],
    return = FALSE,
    verbose = FALSE
  )
  invisible(output_tsv)
}

#' Split the rrad tags out of sequence IDs
#'
#' `parseRadSequenceIds` reads the barcode and the UMI from the sequence IDs of
#' an AIRR table. It adds them as columns.
#'
#' @details
#' The function reads IDs of the form `id|BARCODE=YYY|UMI=XXX`.
#'
#' If an ID has no tag, the function keeps the values in an existing `barcode`
#' or `umi` column. A conflict between values, or a repeated tag, gives an
#' error. A missing UMI is allowed. For other header formats, use
#' [extractSpatialIds()].
#'
#' The metadata must be in the ID, before any white space, so that it stays in
#' the IgBLAST output. For unformatted rrad output with SAM-style comments,
#' reformat the reads first with `header_format = "presto"`.
#'
#' @param airr_df AIRR-style data frame.
#' @param id Name of the column that holds the read IDs (default
#'   `"sequence_id"`).
#' @param strip Replace the ID column with the read ID (default `TRUE`)? The
#'   read ID is the part before the first `|` or white space.
#'
#' @return The data frame with the columns `umi` and `barcode` added.
#' @seealso [extractSpatialIds()], [processSpatialSample()]
#'
#' @examples
#' reads <- data.frame(
#'   sequence_id = c("read1|BARCODE=ACGTACGT|UMI=TTGCA", "read2|BARCODE=GGGTTTAA|UMI=CCATG")
#' )
#' parseRadSequenceIds(reads)
#' @export
parseRadSequenceIds <- function(airr_df, id = "sequence_id", strip = TRUE) {
  .check_required_cols(airr_df, id, "airr_df")
  out <- .as_data_frame(airr_df)
  ids <- as.character(out[[id]])
  for (field in c("barcode", "umi")) {
    # Tags end at the next `|` or at whitespace: a FASTA description after
    # the identifier is not part of the tag value.
    rx <- stringr::regex(
      paste0("\\|", field, "[=:]([^|[:space:]]*)"),
      ignore_case = TRUE
    )
    repeated <- which(stringr::str_count(ids, rx) > 1L)
    if (length(repeated)) {
      rlang::abort(paste0("Duplicate `", field, "` tag at row ", repeated[[1L]], "."))
    }
    value <- stringr::str_match(ids, rx)[, 2]
    present <- !is.na(value) & nzchar(value)
    invalid <- which(present & grepl("[=,[:space:]]", value))
    if (length(invalid)) {
      rlang::abort(paste0(
        "The `", field, "` tag contains a delimiter or whitespace at row ",
        invalid[[1L]], "."
      ))
    }
    existing <- if (field %in% colnames(out)) {
      as.character(out[[field]])
    } else {
      rep(NA_character_, nrow(out))
    }
    conflict <- which(present & !is.na(existing) & nzchar(existing) & value != existing)
    if (length(conflict)) {
      rlang::abort(paste0(
        "Conflicting `", field, "` values in the sequence ID and AIRR column at row ",
        conflict[[1L]], "."
      ))
    }
    existing[present] <- value[present]
    out[[field]] <- existing
  }
  if (strip) {
    out[[id]] <- sub("[|[:space:]].*$", "", ids)
  }
  out
}
# UMI-aware pre-consensus stage for processSpatialSample(umi_support = TRUE).
#
# The default pipeline clusters reads by spatial barcode only. Two reads of one
# molecule that differ by a sequencing error therefore cluster apart at
# `identity = 1.0`, and each is dropped as a singleton by the
# `BuildConsensus -n` floor. This stage groups reads by barcode AND UMI first,
# at a relaxed identity, so those reads support each other before the barcode
# step runs.
#
# Per (barcode, UMI) cluster it records:
#   cc_bc_umi           reads in the cluster, counted BEFORE the productive
#                       subset, so support is not deflated by error-bearing
#                       reads that the consensus would have corrected
#   n_productive_reads  productive reads, which is what the consensus is built
#                       from. A cluster with none is discarded.
#
# Both travel into the barcode stage as header fields and are summed there with
# `BuildConsensus --act sum`, which is what the `min_umi_support` floor uses.

# Read a FASTA into a data frame of id and seq. Only used for the small
# intermediate files this stage writes, not for raw read data.
.sirup_read_fasta <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0L) {
    return(data.frame(id = character(), seq = character(), stringsAsFactors = FALSE))
  }
  is_header <- startsWith(lines, ">")
  ids <- sub("^>", "", lines[is_header])
  grp <- cumsum(is_header)
  seqs <- vapply(
    split(lines[!is_header], grp[!is_header]),
    function(x) paste(x, collapse = ""),
    character(1)
  )
  seqs <- seqs[as.character(seq_along(ids))]
  data.frame(id = ids, seq = unname(ifelse(is.na(seqs), "", seqs)),
             stringsAsFactors = FALSE)
}

# Pull one `field=value` annotation out of pRESTO-style FASTA headers. The match
# allows start-of-string as well as `|`, because the leading token is the
# sequence ID itself (written here as `sequence_id=<id>`) rather than a
# delimited annotation. pRESTO upper-cases annotation keys on output, so the
# match is case-insensitive.
.sirup_header_field <- function(ids, field) {
  m <- stringr::str_match(ids, paste0("(?i)(?:^|[|])", field, "=([^|]*)"))
  m[, 2]
}

# Weight each molecule consensus by its read support for the barcode stage.
# BuildConsensus gives every record one vote, and its `--maxgap`, `--freq` and
# `--maxerror` calls count records. A molecule built from 7 reads would count as
# much as a molecule built from 1 read, so one short read could trim the
# consensus of several full-length molecules. Write each consensus
# `n_productive_reads` times instead, so one vote is one read, as before the UMI
# stage. Only the first copy carries `cc_bc_umi` and `n_productive_reads`; the
# other copies carry 0. The `--act sum` totals are then unchanged.
.sirup_replicate_by_support <- function(grp, seq, barcode, cc_bc_umi,
                                        n_productive_reads, c_call = NULL,
                                        cluster_field,
                                        copies = n_productive_reads) {
  reps <- pmax(as.integer(rep_len(copies, length(grp))), 1L)
  i <- rep(seq_along(grp), reps)
  copy <- sequence(reps)
  first <- copy == 1L
  id <- paste0(
    "sequence_id=", grp[i], ifelse(first, "", paste0("_r", copy)),
    "|", cluster_field, "=", barcode[i],
    "|cc_bc_umi=", ifelse(first, cc_bc_umi[i], 0L),
    "|n_productive_reads=", ifelse(first, n_productive_reads[i], 0L),
    if (!is.null(c_call)) paste0("|c_call=", c_call[i]) else ""
  )
  data.frame(id = id, seq = seq[i], stringsAsFactors = FALSE)
}

# Build one consensus per (barcode, UMI) cluster and return the FASTA that the
# barcode stage should consume, plus the per-consensus support counts.
.sirup_umi_consensus <- function(df,
                                 output_dir,
                                 sample_name,
                                 receptor,
                                 cluster_field,
                                 umi_field,
                                 umi_identity,
                                 umi_weight = TRUE,
                                 consensus_freq,
                                 consensus_min_count,
                                 consensus_max_error,
                                 consensus_max_gap,
                                 delimiters,
                                 muscle_exec,
                                 nproc,
                                 fasta_sequence_col,
                                 id_col,
                                 c_call_col = "c_call",
                                 align_log = NULL,
                                 consensus_log = NULL) {
  .check_required_cols(df, c(cluster_field, umi_field, id_col, fasta_sequence_col), "airr_df")
  if (!"productive" %in% colnames(df)) {
    rlang::abort("`umi_support = TRUE` needs a `productive` column to pick the reads to build each UMI consensus from.")
  }

  barcode <- as.character(df[[cluster_field]])
  umi <- as.character(df[[umi_field]])
  bad_umi <- is.na(umi) | !nzchar(umi)
  if (all(bad_umi)) {
    rlang::abort(paste0("`", umi_field, "` is empty for every row, so there is nothing to group on."))
  }
  if (any(bad_umi)) {
    rlang::warn(paste0(sum(bad_umi), " of ", length(umi), " rows have no `",
                       umi_field, "`. Dropping them."))
    df <- df[!bad_umi, , drop = FALSE]
    barcode <- barcode[!bad_umi]
    umi <- umi[!bad_umi]
  }
  delim_class <- .sirup_delims_regex_class(delimiters)
  if (any(grepl(delim_class, umi))) {
    rlang::abort(paste0("`", umi_field, "` contains a pRESTO delimiter."))
  }

  read_id <- as.character(df[[id_col]])
  bc_umi <- paste0(barcode, "_", umi)

  # --- 1. cluster within barcode + UMI at the relaxed identity ------------
  stage1_fasta <- file.path(output_dir, sprintf("%s_%s_umi_reads.fasta", sample_name, receptor))
  writeFastqas(
    data.frame(id = paste0("sequence_id=", read_id, "|bc_umi=", bc_umi),
               seq = df[[fasta_sequence_col]], stringsAsFactors = FALSE),
    fn = stage1_fasta, type = "fa", header_col = "id"
  )
  message("UMI stage: ", nrow(df), " reads in ", length(unique(bc_umi)), " barcode-UMI pairs")

  cluster_fasta <- sub("\\.fasta$", sprintf("_ident_%s_cluster-pass.fasta", umi_identity), stage1_fasta)
  .sirup_run_command("ClusterSets.py", c(
    "set", "-s", path.expand(stage1_fasta), "-f", "bc_umi", "-k", "cluster",
    "--ident", as.character(umi_identity), "--cluster", "vsearch",
    "--delim", delimiters, "--nproc", as.character(nproc),
    "-o", path.expand(cluster_fasta)
  ))

  # --- 2. count support, then keep the productive reads -------------------
  clustered <- .sirup_read_fasta(cluster_fasta)
  key <- data.frame(
    read_id = .sirup_header_field(clustered$id, "sequence_id"),
    bc_umi = .sirup_header_field(clustered$id, "bc_umi"),
    cluster = .sirup_header_field(clustered$id, "cluster"),
    stringsAsFactors = FALSE
  )
  key$grp <- paste0(key$bc_umi, "_", key$cluster)
  idx <- match(key$read_id, read_id)
  if (anyNA(idx)) {
    rlang::abort("ClusterSets output holds read IDs that are absent from the input table.")
  }
  productive <- df$productive[idx]
  productive <- !is.na(productive) & productive %in% c(TRUE, "TRUE", "T")

  cc_bc_umi <- as.integer(table(key$grp)[key$grp])
  n_prod <- as.integer(tapply(productive, key$grp, sum)[key$grp])
  keep <- productive & n_prod > 0L
  n_groups <- length(unique(key$grp))
  message("UMI stage: ", n_groups, " barcode-UMI clusters; ",
          length(unique(key$grp[n_prod == 0L])), " have no productive read and are dropped")
  if (!any(keep)) {
    rlang::abort("No barcode-UMI cluster holds a productive read.")
  }

  align_in <- file.path(output_dir, sprintf("%s_%s_umi_productive.fasta", sample_name, receptor))
  writeFastqas(
    data.frame(id = paste0("sequence_id=", key$read_id[keep], "|bc_umi_cluster=", key$grp[keep]),
               seq = clustered$seq[keep], stringsAsFactors = FALSE),
    fn = align_in, type = "fa", header_col = "id"
  )

  # --- 3. align and build one consensus per barcode-UMI cluster -----------
  # -n 1 on purpose: a cluster of one read still carries its support forward,
  # and the barcode stage decides whether the summed support is enough.
  aligned <- sub("\\.fasta$", "_align-pass.fasta", align_in)
  .sirup_run_command("AlignSets.py", c(
    "muscle", "-s", path.expand(align_in),
    if (is.null(align_log)) NULL else c("--log", path.expand(align_log)),
    "--delim", delimiters, "--nproc", as.character(nproc),
    "--bf", "bc_umi_cluster", "--exec", path.expand(muscle_exec)
  ))

  umi_consensus <- file.path(output_dir, sprintf("%s_%s_umi_consensus.fasta", sample_name, receptor))
  .sirup_run_command("BuildConsensus.py", c(
    "-s", path.expand(aligned), "-o", path.expand(umi_consensus),
    if (is.null(consensus_log)) NULL else c("--log", path.expand(consensus_log)),
    "--nproc", as.character(nproc), "-n", "1", "-q", "20",
    "--freq", as.character(consensus_freq), "--bf", "bc_umi_cluster",
    "--maxerror", as.character(consensus_max_error),
    "--maxgap", as.character(consensus_max_gap),
    "--fasta", "--delim", delimiters
  ))

  # --- 4. hand the barcode stage one record per molecule ------------------
  cons <- .sirup_read_fasta(umi_consensus)
  cons$grp <- sub("[|].*$", "", cons$id)
  per_grp <- data.frame(
    grp = key$grp, cc_bc_umi = cc_bc_umi, n_productive_reads = n_prod,
    barcode = barcode[idx], stringsAsFactors = FALSE
  )
  per_grp <- per_grp[!duplicated(per_grp$grp), , drop = FALSE]
  if (c_call_col %in% colnames(df)) {
    calls <- as.character(df[[c_call_col]])[idx]
    calls[is.na(calls) | !nzchar(calls)] <- "NA"
    # Majority call among the reads of the cluster, ties broken alphabetically,
    # matching how BuildConsensus resolves a copy field.
    tab <- as.data.frame(table(key$grp, calls), stringsAsFactors = FALSE)
    names(tab) <- c("grp", "call", "n")
    tab <- tab[tab$n > 0, , drop = FALSE]
    tab <- tab[order(tab$grp, -tab$n, tab$call), , drop = FALSE]
    tab <- tab[!duplicated(tab$grp), , drop = FALSE]
    per_grp$c_call <- tab$call[match(per_grp$grp, tab$grp)]
  }
  m <- match(cons$grp, per_grp$grp)
  if (anyNA(m)) {
    rlang::abort("The UMI consensus FASTA holds groups that were never counted.")
  }
  out <- per_grp[m, , drop = FALSE]
  records <- .sirup_replicate_by_support(
    grp = cons$grp, seq = cons$seq, barcode = out$barcode,
    cc_bc_umi = out$cc_bc_umi, n_productive_reads = out$n_productive_reads,
    c_call = out$c_call, cluster_field = cluster_field,
    copies = if (isTRUE(umi_weight)) out$n_productive_reads else 1L
  )
  spatial_fasta <- file.path(output_dir, sprintf("%s_%s_spatial.fasta", sample_name, receptor))
  writeFastqas(records, fn = spatial_fasta, type = "fa", header_col = "id")

  message("UMI stage: wrote ", nrow(cons), " molecule consensuses (", nrow(records),
          " records after weighting by support) to ", spatial_fasta)
  message("UMI stage: reads per molecule, median ", stats::median(out$cc_bc_umi),
          ", ", round(100 * mean(out$cc_bc_umi >= 2), 1), "% with 2 or more")
  list(
    spatial_fasta = spatial_fasta,
    umi_reads_fasta = stage1_fasta,
    umi_cluster_fasta = cluster_fasta,
    umi_consensus_fasta = umi_consensus,
    n_reads = nrow(df),
    n_groups = n_groups,
    n_molecules = nrow(cons)
  )
}

# ---- Spatial ID and tissue position helpers ----------------------------------

#' Read a spatial ID mapping file
#'
#' `readSpatialIdMapping` reads a table that maps spatial IDs to other values,
#' such as tissue positions. The file can be CSV, TSV, or parquet.
#'
#' @details
#' The function chooses the format from the file extension. A parquet file needs
#' the `arrow` package.
#'
#' @param path Path of the mapping file.
#'
#' @return A data frame.
#' @seealso [addTissuePositions()]
#'
#' @examples
#' positions <- data.frame(
#'   barcode = c("s1", "s2"),
#'   pxl_col_in_fullres = c(100, 200),
#'   pxl_row_in_fullres = c(50, 60)
#' )
#' path <- tempfile(fileext = ".csv")
#' write.csv(positions, path, row.names = FALSE)
#' readSpatialIdMapping(path)
#' @export
readSpatialIdMapping <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "parquet") {
    .require_namespace("arrow")
    return(.as_data_frame(arrow::read_parquet(path)))
  }
  if (ext %in% c("tsv", "tab")) {
    return(.as_data_frame(data.table::fread(path, sep = "\t")))
  }
  .as_data_frame(data.table::fread(path))
}

#' Extract spatial IDs from sequence IDs
#'
#' `extractSpatialIds` copies the spatial ID out of the sequence ID with a
#' regular expression. It adds the spatial ID as a column.
#'
#' @param airr_df AIRR-style data frame.
#' @param id Name of the column that holds the read or consensus IDs (default
#'   `"sequence_id"`).
#' @param spatial_id Name of the new spatial ID column (default `"spatial_id"`).
#' @param regex Regular expression that matches the spatial ID. If it has a
#'   capture group, the function uses the first group. Otherwise, it uses the
#'   whole match.
#' @param prefix Text to add before the match (default `""`).
#' @param suffix Text to add after the match (default `""`).
#'
#' @return `airr_df` with the `spatial_id` column added. A row with no match
#'   gets `NA`.
#' @seealso [parseRadSequenceIds()], [addTissuePositions()]
#'
#' @examples
#' reads <- data.frame(
#'   sequence_id = c("consensus1,ACGTACGTACGTACGT", "consensus2,TTTTGGGGCCCCAAAA")
#' )
#'
#' # Take the 16-nucleotide barcode and add the suffix of the Seurat spot names.
#' extractSpatialIds(reads, regex = ",([ACGT]{16})$", suffix = "-1")
#' @export
extractSpatialIds <- function(airr_df,
                                    id = "sequence_id",
                                    spatial_id = "spatial_id",
                                    regex,
                                    prefix = "",
                                    suffix = "") {
  .check_required_cols(airr_df, id, "airr_df")
  out <- .as_data_frame(airr_df)
  ids <- as.character(out[[id]])

  # str_match makes a matrix with the full match in the first column and capture groups in subsequent columns
  match_matrix <- stringr::str_match(ids, regex)
  if (ncol(match_matrix) > 1) {
    extracted_spatial_id <- match_matrix[, 2]
  } else {
    extracted_spatial_id <- match_matrix[, 1]
  }
  out[[spatial_id]] <- ifelse(is.na(extracted_spatial_id), NA_character_, paste0(prefix, extracted_spatial_id, suffix))
  out
}

#' Add tissue coordinates to an AIRR table
#'
#' `addTissuePositions` joins the tissue positions of the spots to an AIRR
#' table. It adds the columns `x` and `y`.
#'
#' @param airr_df AIRR-style data frame with a spatial ID column.
#' @param positions Data frame of tissue positions, or the path of a CSV or
#'   parquet file.
#' @param spatial_id Name of the column in `airr_df` to match (default
#'   `"spatial_id_8um"`).
#' @param position_spatial_id Name of the spatial ID column in `positions`
#'   (default `"barcode"`).
#' @param x Name of the column in `positions` to copy to `x` (default
#'   `"pxl_col_in_fullres"`).
#' @param y Name of the column in `positions` to copy to `y` (default
#'   `"pxl_row_in_fullres"`).
#' @param keep_position_cols Keep all columns of `positions` (default `TRUE`)?
#'   If `FALSE`, the result has the columns of `airr_df`, `x`, and `y`.
#'
#' @return A data frame with the columns `x` and `y` added.
#' @seealso [readSpatialIdMapping()], [extractSpatialIds()]
#'
#' @examples
#' reads <- data.frame(sequence_id = c("a", "b"), spatial_id_8um = c("s1", "s2"))
#' positions <- data.frame(
#'   barcode = c("s1", "s2"),
#'   pxl_col_in_fullres = c(100, 200),
#'   pxl_row_in_fullres = c(50, 60)
#' )
#' addTissuePositions(reads, positions)
#' @export
addTissuePositions <- function(airr_df,
                                 positions,
                                 spatial_id = "spatial_id_8um",
                                 position_spatial_id = "barcode",
                                 x = "pxl_col_in_fullres",
                                 y = "pxl_row_in_fullres",
                                 keep_position_cols = TRUE) {
  .check_required_cols(airr_df, spatial_id, "airr_df")
  if (is.character(positions) && length(positions) == 1) {
    positions <- readSpatialIdMapping(positions)
  }
  .check_required_cols(positions, c(position_spatial_id, x, y), "positions")

  out <- dplyr::left_join(
    .as_data_frame(airr_df),
    .as_data_frame(positions),
    by = stats::setNames(position_spatial_id, spatial_id)
  )
  out$x <- out[[x]]
  out$y <- out[[y]]

  if (!keep_position_cols) {
    keep <- union(colnames(airr_df), c("x", "y"))
    out <- out[, keep, drop = FALSE]
  }
  out
}
