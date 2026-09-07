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
#' @param complete_vdj Keep only rows whose IgBLAST `complete_vdj` field is
#'   true. When enabled, the field is required so a malformed first-pass table
#'   cannot silently bypass this filter.
#' @param min_consensus_count Optional minimum consensus_count.
#' @param max_opening_gap Optional maximum opening gap in `sequence_alignment`.
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

  if (productive_only && "productive" %in% colnames(out)) {
    out <- out[out$productive %in% c(TRUE, "TRUE", "T"), , drop = FALSE]
  }
  if (complete_vdj) {
    .check_required_cols(out, "complete_vdj", "airr_df")
    out <- out[out$complete_vdj %in% c(TRUE, "TRUE", "T"), , drop = FALSE]
  }
  if (!is.null(min_consensus_count) && "consensus_count" %in% colnames(out)) {
    out <- out[!is.na(out$consensus_count) & out$consensus_count >= min_consensus_count, , drop = FALSE]
  }
  if (!is.null(max_opening_gap)) {
    .check_required_cols(out, sequence_alignment, "airr_df")
    gaps <- vapply(out[[sequence_alignment]], .sirap_opening_gap_length, integer(1))
    out <- out[gaps <= max_opening_gap, , drop = FALSE]
  }

  if (!is.null(spatial_id)) {
    .check_required_cols(out, spatial_id, "airr_df")
  }

  if (add_spatial_sequence_id) {
    if (is.null(spatial_id)) {
      rlang::abort("`spatial_id` is required when `add_spatial_sequence_id` is TRUE.")
    }
    .check_required_cols(out, c(id, spatial_id), "airr_df")
    out[[spatial_sequence_id]] <- .sirap_format_spatial_headers(
      out[[id]], out[[spatial_id]], barcode_key = spatial_id_key
    )
  }

  message(
    "filterAirr kept ", nrow(out), " of ", n_input,
    " rows (dropped ", n_input - nrow(out), ")."
  )
  out
}

.sirap_process_input_path <- function(input) {
  if (inherits(input, "rad_reformat_result")) {
    candidates <- c(input$artifacts$fasta, input$artifacts$fastx,
                    input$artifacts$fastq)
  } else if (inherits(input, c("rad_demux_result", "sirap_rad_result"))) {
    candidates <- input$artifacts$fastq
  } else {
    return(input)
  }
  candidates <- candidates[
    is.character(candidates) & !is.na(candidates) & nzchar(candidates)
  ]
  if (!length(candidates)) {
    rlang::abort("The RAD result does not contain an aggregate sequence artifact.")
  }
  candidates[[1L]]
}

.sirap_process_input_type <- function(input, input_type) {
  input_type <- match.arg(
    input_type, c("auto", "rad", "rad_fasta", "airr", "spatial_fasta")
  )
  if (input_type != "auto") {
    if (is.data.frame(input) && input_type != "airr") {
      rlang::abort("A data-frame input requires `input_type = \"airr\"`.")
    }
    return(input_type)
  }
  if (inherits(input, "rad_reformat_result")) {
    fasta <- input$artifacts$fasta
    if (is.character(fasta) && length(fasta) == 1L && !is.na(fasta) &&
        nzchar(fasta)) {
      return("rad_fasta")
    }
    return("rad")
  }
  if (inherits(input, c("rad_demux_result", "sirap_rad_result"))) {
    return("rad")
  }
  if (is.data.frame(input)) {
    return("airr")
  }
  if (!is.character(input) || length(input) != 1L || is.na(input) ||
      !nzchar(input)) {
    rlang::abort("`input_fasta` must be a RAD result, data frame, or one file path.")
  }
  base <- sub("\\.gz$", "", input, ignore.case = TRUE)
  ext <- tolower(tools::file_ext(base))
  if (ext %in% c("fq", "fastq")) {
    return("rad")
  }
  if (ext %in% c("tsv", "tab", "txt")) {
    return("airr")
  }
  if (ext %in% c("fa", "fasta", "fna")) {
    # The safe default is the complete two-pass workflow. Callers resuming
    # from a FASTA that has already passed read-level IgBLAST must opt into the
    # legacy bypass with input_type = "spatial_fasta".
    return("rad_fasta")
  }
  rlang::abort(paste0("Cannot infer an input type from: ", input, "."))
}

.sirap_process_sample_name <- function(input, input_path, sample_name) {
  if (!is.null(sample_name)) {
    if (!is.character(sample_name) || length(sample_name) != 1L ||
        is.na(sample_name) || !nzchar(sample_name) ||
        grepl("[/\\\\]", sample_name)) {
      rlang::abort("`sample_name` must be one non-empty filename stem.")
    }
    return(sample_name)
  }
  if (is.data.frame(input)) {
    return("sample")
  }
  stem <- sub("\\.gz$", "", basename(input_path), ignore.case = TRUE)
  tools::file_path_sans_ext(stem)
}

.sirap_inspect_rad_fasta <- function(path) {
  magic <- readBin(path, what = "raw", n = 2L)
  compressed <- length(magic) == 2L &&
    identical(as.integer(magic), c(31L, 139L))
  connection <- if (compressed) gzfile(path, open = "rt") else {
    file(path, open = "rt")
  }
  on.exit(close(connection), add = TRUE)
  header <- character()
  repeat {
    line <- readLines(connection, n = 1L, warn = FALSE)
    if (!length(line)) break
    if (nzchar(line)) {
      header <- line
      break
    }
  }
  if (!length(header) || !startsWith(header, ">")) {
    rlang::abort("A RAD FASTA input must begin with a FASTA header.")
  }
  tagged <- grepl(
    "(?:^|[\\t ])(?:CB|BC|UB):[A-Za-z]:[^\\t ]*",
    header,
    perl = TRUE
  )
  list(
    compressed = compressed,
    header_format = if (tagged) "presto" else "preserve"
  )
}

.sirap_prepare_rad_fasta <- function(path, output_dir, sample_name, threads) {
  contract <- .sirap_inspect_rad_fasta(path)
  if (!contract$compressed && identical(contract$header_format, "preserve")) {
    return(list(path = path, result = NULL))
  }
  output <- file.path(output_dir, paste0(sample_name, "_rad_plain.fasta"))
  result <- rad_reformat(
    fastq = path,
    output_fastq = output,
    threads = threads,
    to_fasta = TRUE,
    header_format = contract$header_format
  )
  list(path = result$artifacts$fasta, result = result)
}

#' Process RAD reads through two IgBLAST passes and spatial consensus
#'
#' The full workflow is intentionally two-pass. RAD first demultiplexes the
#' reads and [rad_reformat()] streams the tagged FASTQ to a plain FASTA.
#' [runIgblast()] then assigns V/D/J calls to every read. SIRAP restores the
#' RAD barcode/UMI fields, filters the read-level AIRR table, and writes a
#' canonical pRESTO FASTA. pRESTO builds spatial consensuses, after which
#' `AssignGenes.py igblast` and `MakeDb.py` perform the second V/D/J annotation
#' on the consensus sequences.
#'
#' A `rad_demux_result` can be passed directly. Bare FASTA paths default to
#' `input_type = "rad_fasta"` so the read-level pass is never silently skipped.
#' Gzip-compressed RAD FASTA is streamed to a plain staged FASTA automatically,
#' and SAM-style RAD tags are normalized to the canonical pRESTO header first.
#' Existing callers resuming from an already annotated, pRESTO-ready FASTA must
#' explicitly use `input_type = "spatial_fasta"` to bypass the first pass.
#'
#' @param input_fasta A RAD demultiplex/reformat result, RAD-tagged FASTQ, RAD
#'   FASTA, read-level AIRR table/path, or pRESTO-ready spatial FASTA.
#' @param output_dir Directory for all generated artifacts.
#' @param receptor Either `"ig"` or `"tr"`.
#' @param organism IgBLAST organism name.
#' @param input_type One of `"auto"`, `"rad"`, `"rad_fasta"`, `"airr"`, or
#'   `"spatial_fasta"`.
#' @param igblast_db_dir Optional legacy directory containing first-pass
#'   `IG_V`, `IG_D`, and `IG_J` database prefixes. Normally leave this `NULL`
#'   and use the conventional databases beneath `igblast_base` for both passes.
#' @param aux_file Optional first-pass IgBLAST auxiliary-data file. In this
#'   two-pass wrapper it must resolve to the conventional
#'   `<igblast_base>/optional_file/<organism>_gl.aux`, because `AssignGenes.py`
#'   uses that same file for pass 2. [runIgblast()] alone permits other files.
#' @param germline_db_v,germline_db_d,germline_db_j,c_region_db Optional
#'   database overrides used for both IgBLAST passes.
#' @param filter_args Named [filterAirr()] overrides.
#' @param fasta_sequence_col Read-level AIRR sequence column sent to pRESTO.
#' @param sample_name Output filename stem.
#' @param identity pRESTO clustering identity.
#' @param consensus_freq pRESTO consensus frequency.
#' @param germline_glob IMGT germline files for the second-pass `MakeDb.py`.
#'   These must be the same reference sequences used to build the selected
#'   IgBLAST V/D/J databases.
#' @param igblast_base IgBLAST installation shared by the direct read-level
#'   pass and consensus-level `AssignGenes.py` pass. Its conventional
#'   `database`, `optional_file`, `internal_data`, and `bin` directories are
#'   used for pass 1.
#' @param igblast_exec IgBLAST executable used by both passes.
#' @param muscle_exec MUSCLE executable used by pRESTO.
#' @param nproc Total CPU budget for each external stage. The read-level
#'   [runIgblast()] pass shares this budget across as many as four concurrent,
#'   base-balanced IgBLAST processes; small inputs stay single-process.
#' @param igblast_max_processes Maximum concurrent processes for the direct
#'   read-level IgBLAST pass. Set to one to force monolithic execution.
#' @param cluster_field pRESTO spatial grouping field.
#' @param delimiters pRESTO field, key/value, and list delimiters.
#' @param parse_fields Fields merged by `ParseHeaders.py`.
#' @param parse_key Name of the merged spatial-cluster field.
#' @param parse_act `ParseHeaders.py` merge action.
#' @param consensus_min_count Minimum reads per consensus.
#' @param consensus_copy_fields Fields copied into consensus headers. `NULL`
#'   copies `cluster_field`, keeping the spatial barcode in the final AIRR file.
#' @param consensus_copy_action Consensus field-copy action.
#' @param align_log Optional alignment log.
#' @param consensus_log Optional consensus log.
#' @return Invisibly, a `sirap_spatial_process_result` with explicit read-level
#'   and consensus-level stage artifacts.
#' @export
processSpatialSample <- function(input_fasta,
                                 output_dir = NULL,
                                 receptor = c("ig", "tr"),
                                 organism = "human",
                                 input_type = c("auto", "rad", "rad_fasta", "airr", "spatial_fasta"),
                                 igblast_db_dir = NULL,
                                 aux_file = NULL,
                                 germline_db_v = NULL,
                                 germline_db_d = NULL,
                                 germline_db_j = NULL,
                                 c_region_db = NULL,
                                 filter_args = list(),
                                 fasta_sequence_col = "sequence_alignment",
                                 sample_name = NULL,
                                 identity = 0.975,
                                 consensus_freq = 0.51,
                                 germline_glob,
                                 igblast_base = "~/software/igblast",
                                 igblast_exec = NULL,
                                 muscle_exec = "~/software/muscle3.8.31_i86linux64",
                                 nproc = 8L,
                                 igblast_max_processes = 4L,
                                 cluster_field = "BARCODE",
                                 delimiters = c("|", "=", ","),
                                 parse_fields = c(cluster_field, "CLUSTER"),
                                 parse_key = paste0(cluster_field, "_CLUSTER"),
                                 parse_act = "cat",
                                 consensus_min_count = 1L,
                                 consensus_copy_fields = NULL,
                                 consensus_copy_action = "majority",
                                 align_log = NULL,
  consensus_log = NULL) {
  receptor <- match.arg(receptor)
  organism <- .sirap_igblast_scalar_character(organism, "organism")
  input_type <- .sirap_process_input_type(input_fasta, input_type)
  input_path <- .sirap_process_input_path(input_fasta)
  if (!is.data.frame(input_fasta)) {
    if (!is.character(input_path) || length(input_path) != 1L ||
        !file.exists(path.expand(input_path)) || dir.exists(path.expand(input_path))) {
      rlang::abort(paste0("Input file does not exist: ", input_path, "."))
    }
    input_path <- normalizePath(path.expand(input_path), mustWork = TRUE)
  }
  if (is.data.frame(input_fasta) && is.null(output_dir)) {
    rlang::abort("`output_dir` is required for an AIRR data-frame input.")
  }
  if (is.null(output_dir)) {
    output_dir <- dirname(input_path)
  }
  output_dir <- path.expand(output_dir)
  if (!dir.exists(output_dir) &&
      !dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)) {
    rlang::abort(paste0("Could not create `output_dir`: ", output_dir, "."))
  }
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  sample_name <- .sirap_process_sample_name(
    input_fasta, input_path, sample_name
  )
  invalid_filter_names <- length(filter_args) &&
    (is.null(names(filter_args)) || anyNA(names(filter_args)) ||
       any(!nzchar(names(filter_args))))
  if (!is.list(filter_args) || invalid_filter_names) {
    rlang::abort("`filter_args` must be a named list.")
  }
  if (!is.character(cluster_field) || length(cluster_field) != 1L ||
      is.na(cluster_field) || !grepl("^[A-Za-z][A-Za-z0-9_]*$", cluster_field)) {
    rlang::abort("`cluster_field` must be one valid pRESTO field name.")
  }
  cluster_field <- toupper(cluster_field)
  parse_fields <- toupper(parse_fields)
  parse_key <- toupper(parse_key)
  if (!is.null(consensus_copy_fields)) {
    consensus_copy_fields <- toupper(consensus_copy_fields)
  }
  igblast_base <- .sirap_igblast_base(igblast_base)
  if (is.null(igblast_base)) {
    rlang::abort(
      "`igblast_base` is required for the consensus-level IgBLAST pass."
    )
  }
  if (!is.null(aux_file)) {
    aux_file <- .sirap_igblast_existing_file(aux_file, "aux_file")
    conventional_aux <- file.path(
      igblast_base, "optional_file", paste0(organism, "_gl.aux")
    )
    if (!file.exists(conventional_aux) || !identical(
      aux_file,
      normalizePath(conventional_aux, mustWork = TRUE)
    )) {
      rlang::abort(paste0(
        "A custom `aux_file` cannot be mirrored by the consensus-level ",
        "AssignGenes.py pass. Put the desired file at ", conventional_aux,
        " or omit `aux_file`."
      ))
    }
  }
  if (is.null(igblast_exec)) {
    igblast_exec <- file.path(igblast_base, "bin", "igblastn")
  }

  reformat_result <- NULL
  read_igblast <- NULL
  rad_fasta <- NULL
  read_airr_path <- NULL
  annotated_airr_path <- NULL
  filtered_airr_path <- NULL
  n_read_airr <- NULL
  n_filtered <- NULL
  pass2_databases <- NULL

  if (input_type == "spatial_fasta") {
    spatial_fasta <- input_path
  } else {
    if (input_type == "rad") {
      rad_fasta <- file.path(output_dir, paste0(sample_name, "_rad.fasta"))
      reformat_result <- rad_reformat(
        fastq = input_path,
        output_fastq = rad_fasta,
        threads = nproc,
        to_fasta = TRUE,
        header_format = "presto"
      )
      rad_fasta <- reformat_result$artifacts$fasta
    } else if (input_type == "rad_fasta") {
      prepared_fasta <- .sirap_prepare_rad_fasta(
        input_path, output_dir, sample_name, nproc
      )
      rad_fasta <- prepared_fasta$path
      reformat_result <- prepared_fasta$result
    }

    if (input_type %in% c("rad", "rad_fasta")) {
      read_airr_path <- file.path(
        output_dir, paste0(sample_name, "_read_igblast.tsv")
      )
      read_igblast <- runIgblast(
        input_fasta = rad_fasta,
        output_tsv = read_airr_path,
        receptor = receptor,
        organism = organism,
        igblast_base = igblast_base,
        igblast_db_dir = igblast_db_dir,
        igblast_exec = igblast_exec,
        aux_file = aux_file,
        germline_db_v = germline_db_v,
        germline_db_d = germline_db_d,
        germline_db_j = germline_db_j,
        c_region_db = c_region_db,
        num_threads = nproc,
        max_processes = igblast_max_processes,
        stage = "read"
      )
      read_airr_path <- read_igblast$artifacts$airr
      pass2_databases <- read_igblast$config$databases
      airr <- .sirap_read_airr(read_airr_path)
    } else if (is.data.frame(input_fasta)) {
      airr <- .as_data_frame(input_fasta)
    } else {
      read_airr_path <- input_path
      airr <- .sirap_read_airr(read_airr_path)
    }

    n_read_airr <- nrow(airr)
    airr <- parseRadSequenceIds(airr, require_barcode = FALSE)
    annotated_airr_path <- file.path(
      output_dir, paste0(sample_name, "_read_annotated.tsv")
    )
    data.table::fwrite(airr, annotated_airr_path, sep = "\t")
    filter_defaults <- list(
      complete_vdj = TRUE,
      min_consensus_count = NULL,
      spatial_id = "barcode",
      add_spatial_sequence_id = FALSE
    )
    selected_filter_args <- utils::modifyList(filter_defaults, filter_args)
    spatial_col <- selected_filter_args$spatial_id
    if (is.null(spatial_col) || !is.character(spatial_col) ||
        length(spatial_col) != 1L || is.na(spatial_col) || !nzchar(spatial_col)) {
      rlang::abort("`filter_args$spatial_id` must name the spatial barcode column.")
    }
    .check_required_cols(airr, spatial_col, "read-level AIRR table")
    missing_spatial <- is.na(airr[[spatial_col]]) | !nzchar(airr[[spatial_col]])
    if (any(missing_spatial)) {
      rlang::abort(paste0(
        "Read-level AIRR metadata is missing `", spatial_col, "` at row ",
        which(missing_spatial)[[1L]], "."
      ))
    }
    filtered <- do.call(
      filterAirr, c(list(airr_df = airr), selected_filter_args)
    )
    n_filtered <- nrow(filtered)
    if (!n_filtered) {
      rlang::abort("No read-level sequences passed `filterAirr()`.")
    }
    filtered_airr_path <- file.path(
      output_dir, paste0(sample_name, "_read_filtered.tsv")
    )
    data.table::fwrite(filtered, filtered_airr_path, sep = "\t")
    spatial_fasta <- file.path(
      output_dir, paste0(sample_name, "_spatial.fasta")
    )
    spatial_fasta <- .sirap_write_spatial_fasta(
      filtered,
      spatial_fasta,
      id = "sequence_id",
      barcode = spatial_col,
      umi = "umi",
      sequence = fasta_sequence_col,
      barcode_key = cluster_field
    )
  }

  if (missing(germline_glob) || !is.character(germline_glob) ||
      !length(germline_glob)) {
    rlang::abort("`germline_glob` is required for consensus MakeDb processing.")
  }
  germline_files <- unique(unlist(lapply(
    path.expand(germline_glob), Sys.glob
  ), use.names = FALSE))
  if (!length(germline_files)) {
    rlang::abort(paste0(
      "`germline_glob` matched no files: ",
      paste(germline_glob, collapse = ", "), "."
    ))
  }
  paths <- .presto_paths(spatial_fasta, identity, consensus_freq, output_dir)
  if (is.null(consensus_copy_fields)) {
    consensus_copy_fields <- cluster_field
  }
  commands <- .sirap_presto_commands(
    paths = paths,
    receptor = receptor,
    organism = organism,
    identity = identity,
    consensus_freq = consensus_freq,
    germline_files = germline_files,
    igblast_base = igblast_base,
    igblast_exec = igblast_exec,
    igblast_dbs = pass2_databases,
    muscle_exec = muscle_exec,
    nproc = nproc,
    cluster_field = cluster_field,
    delimiters = delimiters,
    parse_fields = parse_fields,
    parse_key = parse_key,
    parse_act = parse_act,
    consensus_min_count = consensus_min_count,
    consensus_copy_fields = consensus_copy_fields,
    consensus_copy_action = consensus_copy_action,
    align_log = align_log,
    consensus_log = consensus_log
  )
  tools <- unique(vapply(commands, function(command) command$exec, character(1)))
  .sirap_check_tools(c(
    stats::setNames(tools, basename(tools)),
    muscle = muscle_exec,
    igblastn = igblast_exec
  ))
  for (command in commands) {
    .sirap_run_stage_command(command)
  }

  artifacts <- list(
    rad_fasta = rad_fasta,
    read_igblast_airr = read_airr_path,
    annotated_read_airr = annotated_airr_path,
    filtered_read_airr = filtered_airr_path,
    spatial_fasta = spatial_fasta,
    clustered_fasta = paths$cluster,
    aligned_fasta = paths$aligned,
    consensus_fasta = paths$consensus,
    consensus_igblast = paths$fmt7,
    final_airr = paths$db_pass
  )
  result <- structure(
    c(
      list(
        status = "success",
        stages = list(
          rad_reformat = reformat_result,
          read_igblast = read_igblast,
          spatial_prep = list(
            n_read_airr = n_read_airr,
            n_filtered = n_filtered
          ),
          consensus = list(commands = commands[seq_len(4L)]),
          consensus_igblast = list(
            command = commands[[5L]],
            make_db = commands[[6L]]
          )
        ),
        artifacts = artifacts,
        # Compatibility aliases for callers of the original FASTA-only API.
        raw_airr = read_airr_path,
        filtered_airr = filtered_airr_path,
        spatial_fasta = spatial_fasta,
        n_raw = n_read_airr,
        n_filtered = n_filtered
      ),
      paths
    ),
    class = c("sirap_spatial_process_result", "list")
  )
  invisible(result)
}

#' Process multiple RAD or spatial samples
#'
#' @param input_fastas Inputs accepted by [processSpatialSample()].
#' @param output_dirs Optional output directory per input.
#' @param sample_names Optional filename stem per input.
#' @param ... Passed to [processSpatialSample()].
#' @return A named list of `sirap_spatial_process_result` objects, invisibly.
#' @export
processSpatialSamples <- function(input_fastas,
                                  output_dirs = NULL,
                                  sample_names = NULL,
                                  ...) {
  one_result <- inherits(input_fastas, c(
    "rad_demux_result", "rad_reformat_result", "sirap_rad_result"
  ))
  inputs <- if (one_result) {
    list(input_fastas)
  } else if (is.list(input_fastas) && !is.data.frame(input_fastas)) {
    input_fastas
  } else {
    as.list(input_fastas)
  }
  if (is.null(output_dirs)) {
    output_dirs <- rep(NA_character_, length(inputs))
  }
  if (length(output_dirs) != length(inputs)) {
    rlang::abort("`output_dirs` must be NULL or match the number of inputs.")
  }
  if (!is.null(sample_names) && length(sample_names) != length(inputs)) {
    rlang::abort("`sample_names` must be NULL or match the number of inputs.")
  }
  result_names <- vapply(seq_along(inputs), function(index) {
    if (!is.null(sample_names)) {
      sample_names[[index]]
    } else {
      path <- .sirap_process_input_path(inputs[[index]])
      if (is.character(path) && length(path) == 1L) basename(path) else {
        paste0("sample_", index)
      }
    }
  }, character(1))
  results <- vector("list", length(inputs))
  names(results) <- result_names
  for (index in seq_along(inputs)) {
    destination <- output_dirs[[index]]
    if (is.na(destination)) destination <- NULL
    results[[index]] <- processSpatialSample(
      input_fasta = inputs[[index]],
      output_dir = destination,
      sample_name = if (is.null(sample_names)) NULL else sample_names[[index]],
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
  if (sequence %in% colnames(out)) {
    out$non_N_nucleotide_prop <- stringr::str_count(out[[sequence]], "[ACGT]") / nchar(out[[sequence]])
  }
  if (alignment %in% colnames(out)) {
    out$max_gap_length <- vapply(out[[alignment]], .sirap_max_gap_length, integer(1))
    out$opening_gap_length <- vapply(out[[alignment]], .sirap_opening_gap_length, integer(1))
  }
  out
}

#' Summarize a spatial receptor table
#'
#' @param airr_df AIRR-style data frame or a SIRAP object.
#' @param spatial_id Spatial ID column to count. When omitted for a SIRAP object,
#'   the spatial ID field recorded in its metadata is used.
#' @return One-row tibble of basic counts.
#' @export
spatialReceptorSummary <- function(airr_df, spatial_id = "spatial_id_2um") {
  spatial_id_missing <- missing(spatial_id)
  inputs <- .sirap_receptor_input(
    airr_df,
    if (isSIRAP(airr_df) && spatial_id_missing) NULL else spatial_id
  )
  airr_df <- inputs$receptors
  spatial_id <- inputs$spatial_id_col %||% "spatial_id_2um"
  .check_required_cols(airr_df, spatial_id, "airr_df")
  tibble::tibble(
    n_sequences = nrow(airr_df),
    n_spatial_ids = dplyr::n_distinct(airr_df[[spatial_id]], na.rm = TRUE),
    n_loci = if ("locus" %in% colnames(airr_df)) dplyr::n_distinct(airr_df$locus, na.rm = TRUE) else NA_integer_,
    n_productive = if ("productive" %in% colnames(airr_df)) sum(airr_df$productive %in% c(TRUE, "TRUE", "T"), na.rm = TRUE) else NA_integer_
  )
}

#' Collapse spatially contiguous spots with identical receptor sequences
#'
#' For each unique receptor sequence, finds connected components of neighboring
#' spots sharing that sequence, then collapses each component to a single
#' representative spot. A SIRAP object is the preferred input: receptor IDs
#' are read from the SIRAP field mapping and spatial coordinates are read from
#' the stored Seurat object. If you are not using a SIRAP object, pass an
#' AIRR-style receptor table as `object` and the matching Seurat object as
#' `seurat_object`.
#'
#' When `object` is a SIRAP object, the returned value is the updated SIRAP
#' object with its receptor table replaced. When `object` is a receptor table,
#' the returned value is the collapsed receptor table.
#'
#' @param object A SIRAP object, or an AIRR-style receptor data frame.
#' @param seurat_object Optional Seurat object used only when `object` is a
#'   receptor data frame.
#' @param spatial_id Column in the receptor table matching Seurat spot/cell
#'   spatial IDs. When omitted for a SIRAP object, the recorded SIRAP spatial ID field
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
#' @param validate Whether to validate the updated SIRAP object before
#'   returning it. Only used when `object` is a SIRAP object.
#' @return A SIRAP object or receptor data frame with bleedover rows removed,
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
  is_sirap <- isSIRAP(object)
  if (is_sirap) {
    if (!is.null(seurat_object)) {
      rlang::abort("`seurat_object` should not be supplied when `object` is a SIRAP object.")
    }
    validateSIRAP(object)
    receptor_fields <- .sirap_receptor_fields(object)
    receptors <- object$receptors
    spatial_id <- spatial_id %||% receptor_fields$spatial_id
    seurat_object <- object$spatial
  } else {
    if (is.null(seurat_object)) {
      rlang::abort("`seurat_object` is required when `object` is a receptor table.")
    }
    receptors <- .as_data_frame(object)
    spatial_id <- .guess_spatial_id_col(receptors, spatial_id)
  }
  options <- .sirap_spatial_options(
    object = if (is_sirap) object else seurat_object,
    image_scale = image_scale,
    swap_coords = swap_coords
  )
  image_scale <- options$image_scale
  swap_coords <- options$swap_coords

  coords <- .sirap_spatial_coordinate_table(
    seurat_object,
    images = images,
    image_scale = image_scale,
    swap_coords = swap_coords
  )
  collapsed <- .collapse_identical_neighbors_table(
    receptors = receptors,
    tissue_positions = coords,
    spatial_id_col = spatial_id,
    sequence_col = sequence,
    consensus_count_col = consensus_count,
    num_neighbors_non_boundary = num_neighbors_non_boundary
  )

  if (!is_sirap) {
    return(collapsed)
  }

  object$receptors <- tibble::as_tibble(collapsed)
  if (isTRUE(validate)) {
    validateSIRAP(object)
  }
  object
}

# Collapse identical neighboring receptors after SIRAP/Seurat inputs have been
# resolved to a receptor table and a coordinate table.
.collapse_identical_neighbors_table <- function(receptors,
                                                tissue_positions,
                                                spatial_id_col,
                                                sequence_col,
                                                consensus_count_col,
                                                num_neighbors_non_boundary) {
  data <- .as_data_frame(receptors)
  .check_required_cols(data, c(spatial_id_col, sequence_col, consensus_count_col), "receptors")
  .check_required_cols(tissue_positions, c("cell", "x", "y"), "spatial coordinates")

  temp <- .collapse_identical_neighbors_temp_cols(data)
  data[[temp$spatial_id]] <- as.character(data[[spatial_id_col]])
  data[[temp$sequence]] <- as.character(data[[sequence_col]])
  data[[temp$row]] <- seq_len(nrow(data))

  tissue_positions <- .as_data_frame(tissue_positions)
  tissue_positions <- tissue_positions[is.finite(tissue_positions$x) & is.finite(tissue_positions$y), , drop = FALSE]
  tissue_positions <- tissue_positions[!duplicated(tissue_positions$cell), c("cell", "x", "y"), drop = FALSE]
  colnames(tissue_positions) <- c(temp$spatial_id, temp$x, temp$y)

  data <- merge(
    data,
    tissue_positions,
    by = temp$spatial_id,
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
  data <- data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(temp$spatial_id, temp$sequence)))) %>%
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
  spatial_ids <- unique(data[[temp$spatial_id]])
  tissue_data <- tissue_positions[tissue_positions[[temp$spatial_id]] %in% spatial_ids, ]
  coords <- as.matrix(tissue_data[, c(temp$x, temp$y)])
  spatial_ids_ordered <- tissue_data[[temp$spatial_id]]
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
    seq_spatial_ids <- data[[temp$spatial_id]][data[[temp$sequence]] == seq_i]
    seq_idx <- spatial_id_to_idx[seq_spatial_ids]

    # For each exact sequence, connected components identify spatially
    # contiguous spots carrying that sequence.
    if (length(seq_idx) == 1L) {
      component_labels[[si]] <- data.frame(
        .sirap_spatial_id = seq_spatial_ids,
        .sirap_sequence = seq_i,
        component = 1L
      )
    } else {
      sub_adj <- adj_sparse[seq_idx, seq_idx, drop = FALSE]
      g <- igraph::graph_from_adjacency_matrix(sub_adj, mode = "undirected", diag = FALSE)
      components <- igraph::components(g)

      component_labels[[si]] <- data.frame(
        .sirap_spatial_id = seq_spatial_ids,
        .sirap_sequence = seq_i,
        component = components$membership
      )
    }

    cli::cli_progress_update()
  }

  cli::cli_progress_done()

  component_df <- do.call(rbind, component_labels)

  cli::cli_progress_step("Collapsing components to representative spots...")

  coords_df <- data.frame(spatial_ids_ordered, cx = coords[, 1], cy = coords[, 2])
  colnames(coords_df)[1] <- temp$spatial_id

  colnames(component_df)[colnames(component_df) == ".sirap_spatial_id"] <- temp$spatial_id
  colnames(component_df)[colnames(component_df) == ".sirap_sequence"] <- temp$sequence
  component_df <- merge(component_df, coords_df, by = temp$spatial_id)
  component_df <- merge(
    component_df,
    data[, c(temp$spatial_id, temp$sequence, consensus_count_col)],
    by = c(temp$spatial_id, temp$sequence)
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

  keep_keys <- paste(representative_df[[temp$spatial_id]], representative_df[[temp$sequence]], sep = "___")
  all_keys <- paste(data[[temp$spatial_id]], data[[temp$sequence]], sep = "___")
  result <- data[all_keys %in% keep_keys, ]
  result <- result[, setdiff(colnames(result), unlist(temp, use.names = FALSE)), drop = FALSE]
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
    spatial_id = ".sirap_collapse_barcode",
    sequence = ".sirap_collapse_sequence",
    row = ".sirap_collapse_row",
    x = ".sirap_collapse_x",
    y = ".sirap_collapse_y"
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

# Build the expected pRESTO/IgBLAST output filenames from the input FASTA stem
# and the clustering/consensus thresholds.
.presto_paths <- function(input_fasta, identity, consensus_freq, output_dir = NULL) {
  if (is.null(output_dir)) {
    output_dir <- dirname(input_fasta)
  }
  stem <- tools::file_path_sans_ext(basename(input_fasta))
  cluster <- file.path(
    output_dir,
    sprintf("%s_ident_%s_barcode_cluster_cluster-pass.fasta", stem, identity)
  )
  reheader <- sub("\\.fasta$", "_reheader.fasta", cluster)
  aligned <- sub("\\.fasta$", "_align-pass.fasta", reheader)
  consensus <- file.path(
    output_dir,
    sprintf("%s_ident_%s_barcode_cluster_freq_%s_cluster-pass_reheader_align-pass_consensus-pass.fasta", stem, identity, consensus_freq)
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

# Build pRESTO/Change-O commands as executable/argument pairs. Keeping the
# arguments structured avoids shell parsing and preserves paths with spaces.
.sirap_presto_commands <- function(paths,
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
                                   align_log,
                                   consensus_log,
                                   igblast_dbs = NULL) {
  if (!is.character(delimiters) || length(delimiters) != 3L ||
      anyNA(delimiters) || any(nchar(delimiters, type = "bytes") != 1L)) {
    rlang::abort("`delimiters` must contain exactly three one-byte values.")
  }
  delim_args <- c("--delim", delimiters)
  align_log_args <- if (is.null(align_log)) {
    NULL
  } else {
    c("--log", path.expand(align_log))
  }
  consensus_log_args <- if (is.null(consensus_log)) {
    NULL
  } else {
    c("--log", path.expand(consensus_log))
  }
  consensus_copy_args <- if (length(consensus_copy_fields)) {
    c("--cf", consensus_copy_fields, "--act", consensus_copy_action)
  } else {
    NULL
  }
  assign_db_args <- character()
  if (!is.null(igblast_dbs)) {
    if (!is.list(igblast_dbs) ||
        !all(c("v", "d", "j") %in% names(igblast_dbs))) {
      rlang::abort("`igblast_dbs` must contain V, D, and J database prefixes.")
    }
    for (segment in c("v", "d", "j", "c")) {
      value <- igblast_dbs[[segment]]
      if (!is.null(value)) {
        option <- c(v = "--vdb", d = "--ddb", j = "--jdb", c = "--cdb")[[segment]]
        assign_db_args <- c(assign_db_args, option, path.expand(value))
      }
    }
  }

  list(
    list(stage = "cluster", exec = "ClusterSets.py", output = paths$cluster,
         args = c(
      "set",
      "-s", path.expand(paths$input),
      "-f", cluster_field,
      "-k", "CLUSTER",
      "--ident", as.character(identity),
      "--cluster", "vsearch",
      delim_args,
      "--nproc", as.character(nproc),
      "-o", path.expand(paths$cluster)
    )),
    list(stage = "merge_headers", exec = "ParseHeaders.py",
         output = paths$reheader, args = c(
      "merge",
      "-s", path.expand(paths$cluster),
      "-o", path.expand(paths$reheader),
      "-f", parse_fields,
      "-k", parse_key,
      delim_args,
      "--act", parse_act
    )),
    list(stage = "align", exec = "AlignSets.py", output = paths$aligned,
         args = c(
      "muscle",
      "-s", path.expand(paths$reheader),
      "-o", path.expand(paths$aligned),
      align_log_args,
      delim_args,
      "--nproc", as.character(nproc),
      "--bf", parse_key,
      "--exec", path.expand(muscle_exec)
    )),
    list(stage = "consensus", exec = "BuildConsensus.py",
         output = paths$consensus, args = c(
      "-s", path.expand(paths$aligned),
      "-o", path.expand(paths$consensus),
      consensus_log_args,
      "--nproc", as.character(nproc),
      "-n", as.character(consensus_min_count),
      "-q", "20",
      "--freq", as.character(consensus_freq),
      "--bf", parse_key,
      "--maxerror", "0.1",
      "--maxgap", "0.3",
      "--fasta",
      delim_args,
      consensus_copy_args
    )),
    list(stage = "consensus_igblast", exec = "AssignGenes.py",
         output = paths$fmt7, args = c(
      "igblast",
      "-s", path.expand(paths$consensus),
      "-o", path.expand(paths$fmt7),
      "-b", path.expand(igblast_base),
      "--organism", organism,
      "--loci", receptor,
      "--format", "blast",
      "--exec", path.expand(igblast_exec),
      "--nproc", as.character(nproc),
      assign_db_args
    )),
    list(stage = "make_db", exec = "MakeDb.py", output = paths$db_pass,
         args = c(
      "igblast",
      "-i", path.expand(paths$fmt7),
      "-s", path.expand(paths$consensus),
      "-r", germline_files,
      "-o", path.expand(paths$db_pass),
      "--extended"
    ))
  )
}

# Run one external spatial-processing stage into a same-directory temporary
# artifact. A successful exit without a non-empty output is still a failure;
# only a validated artifact is atomically published for the next stage.
.sirap_run_stage_command <- function(command) {
  required <- c("stage", "exec", "args", "output")
  if (!is.list(command) || !all(required %in% names(command))) {
    rlang::abort("Malformed external-stage command specification.")
  }
  output <- path.expand(command$output)
  parent <- dirname(output)
  if (!dir.exists(parent)) {
    rlang::abort(paste0(
      "Output directory for stage `", command$stage,
      "` does not exist: ", parent, "."
    ))
  }
  output_positions <- which(command$args == "-o")
  if (length(output_positions) != 1L ||
      output_positions[[1L]] == length(command$args)) {
    rlang::abort(paste0(
      "Stage `", command$stage, "` must define exactly one `-o` artifact."
    ))
  }

  extension <- tools::file_ext(output)
  extension <- if (nzchar(extension)) paste0(".", extension) else ""
  staged <- tempfile(
    pattern = paste0(".sirap-", command$stage, "-"),
    tmpdir = parent,
    fileext = extension
  )
  on.exit(unlink(staged), add = TRUE)
  execution_args <- command$args
  execution_args[output_positions[[1L]] + 1L] <- staged
  .sirap_run_command(command$exec, execution_args)

  if (!file.exists(staged) || isTRUE(file.info(staged)$isdir) ||
      is.na(file.info(staged)$size) || file.info(staged)$size <= 0) {
    rlang::abort(
      paste0(
        "External stage `", command$stage,
        "` exited successfully but produced no non-empty output artifact."
      ),
      class = "sirap_igblast_command_error"
    )
  }
  .sirap_publish_igblast_artifact(
    staged, output, paste0(command$stage, " stage")
  )
}

# Return the longest internal gap run in an aligned receptor sequence.
.sirap_max_gap_length <- function(seq) {
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

# Return the length of the leading gap run in an aligned receptor sequence.
.sirap_opening_gap_length <- function(seq) {
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
