#' Locate the rad CLI binary
#'
#' Resolves the path to the `rad` command-line binary. Resolution order:
#' explicit `path` argument, then `getOption("SIRAP.rad_path")`, then the
#' `SIRAP_RAD_PATH` environment variable, then `Sys.which("rad")`.
#'
#' @param path Optional explicit path. If supplied and exists, returned
#'   normalized.
#' @return Absolute path to the rad binary.
#' @export
rad_path <- function(path = NULL) {
  if (!is.null(path)) {
    if (!nzchar(path) || !file.exists(path)) {
      rlang::abort(paste0("rad binary not found at: ", path))
    }
    return(normalizePath(path))
  }
  opt <- getOption("SIRAP.rad_path")
  if (!is.null(opt) && nzchar(opt) && file.exists(opt)) {
    return(normalizePath(opt))
  }
  env <- Sys.getenv("SIRAP_RAD_PATH", unset = "")
  if (nzchar(env) && file.exists(env)) {
    return(normalizePath(env))
  }
  on_path <- unname(Sys.which("rad"))
  if (nzchar(on_path)) {
    return(on_path)
  }
  rlang::abort(
    paste0(
      "Could not locate the 'rad' CLI binary. Install it from ",
      "https://github.com/indianewok/rad, then either: ",
      "(1) place it on your PATH, ",
      "(2) set options(SIRAP.rad_path = '/path/to/rad'), ",
      "(3) set SIRAP_RAD_PATH=/path/to/rad in your environment, or ",
      "(4) pass rad = '/path/to/rad' explicitly."
    )
  )
}

.rad_flag <- function(flag, value) {
  if (is.null(value)) return(character())
  if (is.logical(value)) {
    if (length(value) != 1L) {
      rlang::abort(paste0("flag '", flag, "' must be length 1"))
    }
    if (isTRUE(value)) return(flag)
    return(character())
  }
  if (length(value) != 1L) {
    rlang::abort(paste0("flag '", flag, "' must be a length-1 scalar"))
  }
  c(flag, as.character(value))
}

.rad_call <- function(args,
                      rad = rad_path(),
                      capture = FALSE,
                      echo = TRUE,
                      check = TRUE) {
  args <- as.character(args)
  if (echo && !capture) {
    cat("$ ", rad, " ", paste(shQuote(args), collapse = " "), "\n", sep = "")
  }
  if (capture) {
    out <- suppressWarnings(
      system2(rad, args = args, stdout = TRUE, stderr = TRUE)
    )
    status <- attr(out, "status")
    if (is.null(status)) status <- 0L
    res <- list(
      status = as.integer(status),
      stdout = paste(out, collapse = "\n"),
      stderr = NA_character_,
      command = c(rad, args)
    )
  } else {
    status <- system2(rad, args = args)
    res <- list(
      status = as.integer(status),
      stdout = NULL,
      stderr = NULL,
      command = c(rad, args)
    )
  }
  if (check && res$status != 0L) {
    rlang::abort(
      paste0("rad exited with status ", res$status, ". Command: ",
             paste(shQuote(res$command), collapse = " "))
    )
  }
  invisible(res)
}

#' Prepare read layout and/or position map (`rad prep`)
#'
#' Wraps `rad prep`. At least one of `read_layout` or `position_map` must be
#' `TRUE`. `position_map = TRUE` requires both `fastq` and `output`.
#'
#' @param layout Layout key (e.g. `five_prime`, `three_prime`, `splitseq`,
#'   `visium`) or path to a layout CSV.
#' @param read_layout Logical; emit the parsed read-layout structure.
#' @param position_map Logical; compute a position map with misalignment
#'   statistics.
#' @param fastq Path to FASTQ/FASTQ.gz input. Required when
#'   `position_map = TRUE`.
#' @param output Output base path; rad appends `_layout.csv` and
#'   `_position_map.csv`.
#' @param max_reads Max reads for misalignment sampling. `NULL` leaves the
#'   rad default (50000).
#' @param threads Number of threads.
#' @param verbose Logical; verbose logging.
#' @param max_verbose Logical; maximum verbosity (implies `verbose = TRUE`).
#' @param ... Additional raw CLI tokens appended verbatim.
#' @param capture Logical; capture stdout/stderr into the return list instead
#'   of streaming to the console.
#' @param rad Path to the rad binary. Defaults to [rad_path()].
#' @return Invisible list with `status`, `stdout`, `stderr`, `command`.
#' @seealso [rad_demux()], [rad_reformat()]
#' @export
rad_prep <- function(layout,
                     read_layout = TRUE,
                     position_map = FALSE,
                     fastq = NULL,
                     output = NULL,
                     max_reads = NULL,
                     threads = 1L,
                     verbose = FALSE,
                     max_verbose = FALSE,
                     ...,
                     capture = FALSE,
                     rad = rad_path()) {
  if (!isTRUE(read_layout) && !isTRUE(position_map)) {
    rlang::abort("Set read_layout = TRUE and/or position_map = TRUE.")
  }
  if (isTRUE(position_map)) {
    if (is.null(fastq))  rlang::abort("position_map = TRUE requires `fastq`.")
    if (is.null(output)) rlang::abort("position_map = TRUE requires `output`.")
  }
  args <- c(
    "prep",
    .rad_flag("--layout", layout),
    if (isTRUE(read_layout))  "--read-layout"  else character(),
    if (isTRUE(position_map)) "--position-map" else character(),
    .rad_flag("--fastq",     fastq),
    .rad_flag("--output",    output),
    .rad_flag("--max-reads", max_reads),
    .rad_flag("--threads",   threads),
    if (isTRUE(verbose))     "--verbose"     else character(),
    if (isTRUE(max_verbose)) "--max-verbose" else character(),
    as.character(c(...))
  )
  .rad_call(args, rad = rad, capture = capture)
}

#' Run the full demultiplexing pipeline (`rad demux`)
#'
#' Wraps `rad demux`. Required: `layout` and `fastq`.
#'
#' @param layout Layout key or path to a layout CSV.
#' @param fastq Path to input FASTQ/FASTQ.gz.
#' @param kit Use this kit's default whitelist.
#' @param global_whitelist Path to a custom global whitelist CSV.
#' @param custom_whitelist Path to a custom whitelist CSV.
#' @param bc_correction_mode `"offensive"` (default) or `"defensive"`.
#' @param joint_bc_mode `"default"` (permissive) or `"strict"`.
#' @param whitelist_mutation Mutation space for the whitelist. rad default 2.
#' @param generated_mutation Mutation space for generated barcodes. rad
#'   default 2.
#' @param max_reads Maximum reads (default: all).
#' @param chunk_size Chunk size for processing. rad default 5000.
#' @param output Filename prefix. rad default `"output"`.
#' @param out_dir Output directory. Defaults to the working directory.
#' @param log_file Path to a log file; rad tees stdout/stderr there.
#' @param write_dbg Logical; write debug files (`.sig`, `.csv`, metrics).
#' @param bc_split Logical; split reads into per-barcode FASTQs (see note in
#'   the rad README about `rad reformat --split-bc` being the current way).
#' @param threads Number of threads.
#' @param verbose,max_verbose Verbosity controls.
#' @param ... Additional raw CLI tokens.
#' @param capture,rad See [rad_prep()].
#' @return Invisible list with `status`, `stdout`, `stderr`, `command`.
#' @export
rad_demux <- function(layout,
                      fastq,
                      kit = NULL,
                      global_whitelist = NULL,
                      custom_whitelist = NULL,
                      bc_correction_mode = NULL,
                      joint_bc_mode = NULL,
                      whitelist_mutation = NULL,
                      generated_mutation = NULL,
                      max_reads = NULL,
                      chunk_size = NULL,
                      output = NULL,
                      out_dir = NULL,
                      log_file = NULL,
                      write_dbg = FALSE,
                      bc_split = FALSE,
                      threads = 1L,
                      verbose = FALSE,
                      max_verbose = FALSE,
                      ...,
                      capture = FALSE,
                      rad = rad_path()) {
  args <- c(
    "demux",
    .rad_flag("--layout", layout),
    .rad_flag("--fastq",  fastq),
    .rad_flag("--kit", kit),
    .rad_flag("--global-whitelist",     global_whitelist),
    .rad_flag("--custom-whitelist",     custom_whitelist),
    .rad_flag("--bc-correction-mode",   bc_correction_mode),
    .rad_flag("--joint-bc-mode",        joint_bc_mode),
    .rad_flag("--whitelist-mutation",   whitelist_mutation),
    .rad_flag("--generated-mutation",   generated_mutation),
    .rad_flag("--max-reads",            max_reads),
    .rad_flag("--chunk-size",           chunk_size),
    .rad_flag("--output",               output),
    .rad_flag("--dir",                  out_dir),
    .rad_flag("--log-file",             log_file),
    if (isTRUE(write_dbg)) "--write-dbg" else character(),
    if (isTRUE(bc_split))  "--bc-split"  else character(),
    .rad_flag("--threads", threads),
    if (isTRUE(verbose))     "--verbose"     else character(),
    if (isTRUE(max_verbose)) "--max-verbose" else character(),
    as.character(c(...))
  )
  .rad_call(args, rad = rad, capture = capture)
}

#' Reformat fastx outputs (`rad reformat`)
#'
#' Wraps `rad reformat`. Required: `fastq`. `out_dir` is required if
#' `split_bc = TRUE`.
#'
#' @param fastq Input FASTQ/FASTA (`.fq` / `.fa` / `.gz`).
#' @param out_dir Output directory for per-barcode FASTQs (required when
#'   `split_bc = TRUE`).
#' @param split_bc Logical; split reads into per-barcode `.fq.gz` files via
#'   the `CB:Z` tag.
#' @param reformat_header Logical; collapse headers to `QNAME_CB_UB`.
#' @param reformat_delim Single character delimiter for collapsed headers
#'   (use `" "` for space).
#' @param coordinate Either `FALSE` (don't rewrite), `TRUE` (use default
#'   mode `vizHD-v1`), or a string MODE.
#' @param bin_size Target bin size in microns for coordinate mode. rad
#'   default 2.
#' @param threads Number of worker threads. rad default 2.
#' @param verbose Logical; verbose logging.
#' @param ... Additional raw CLI tokens.
#' @param capture,rad See [rad_prep()].
#' @return Invisible list with `status`, `stdout`, `stderr`, `command`.
#' @export
rad_reformat <- function(fastq,
                         out_dir = NULL,
                         split_bc = FALSE,
                         reformat_header = FALSE,
                         reformat_delim = NULL,
                         coordinate = FALSE,
                         bin_size = NULL,
                         threads = NULL,
                         verbose = FALSE,
                         ...,
                         capture = FALSE,
                         rad = rad_path()) {
  if (isTRUE(split_bc) && is.null(out_dir)) {
    rlang::abort("split_bc = TRUE requires `out_dir`.")
  }
  coord_arg <- if (isFALSE(coordinate) || is.null(coordinate)) {
    character()
  } else if (isTRUE(coordinate)) {
    "--coordinate"
  } else {
    paste0("--coordinate=", coordinate)
  }
  args <- c(
    "reformat",
    .rad_flag("--fastq",  fastq),
    .rad_flag("--outdir", out_dir),
    if (isTRUE(split_bc))        "--split-bc"        else character(),
    if (isTRUE(reformat_header)) "--reformat-header" else character(),
    .rad_flag("--reformat-delim", reformat_delim),
    coord_arg,
    .rad_flag("--bin-size", bin_size),
    .rad_flag("--threads",  threads),
    if (isTRUE(verbose)) "--verbose" else character(),
    as.character(c(...))
  )
  .rad_call(args, rad = rad, capture = capture)
}

#' List registered layouts and whitelist kits (`rad list`)
#'
#' Wraps `rad list`. If both flags are `FALSE`, rad lists everything.
#'
#' @param layouts Logical; list read-layout keys and resolved paths.
#' @param whitelists Logical; list whitelist kit keys and resolved paths.
#' @param ... Additional raw CLI tokens.
#' @param capture,rad See [rad_prep()]. Defaults to `capture = TRUE` because
#'   the output is structured table-like text most callers want to keep.
#' @return Invisible list with `status`, `stdout` (the listing), `stderr`,
#'   `command`.
#' @export
rad_list <- function(layouts = FALSE,
                     whitelists = FALSE,
                     ...,
                     capture = TRUE,
                     rad = rad_path()) {
  args <- c(
    "list",
    if (isTRUE(layouts))    "--layouts"    else character(),
    if (isTRUE(whitelists)) "--whitelists" else character(),
    as.character(c(...))
  )
  .rad_call(args, rad = rad, capture = capture)
}

#' Add or remove layout / whitelist mappings (`rad modify`)
#'
#' Wraps `rad modify`. Exactly one of `layout_key` / `whitelist_kit` must be
#' set, and exactly one of `add` / `remove`. Note from the rad README:
#' `rad_config set/rm` (this command) is currently process-local — changes
#' do not persist across separate invocations.
#'
#' @param layout_key Read-layout key to add/remove.
#' @param whitelist_kit Whitelist kit key to add/remove.
#' @param add Path to add/update the mapping to. Mutually exclusive with
#'   `remove`.
#' @param remove Logical; remove the mapping. Mutually exclusive with `add`.
#' @param ... Additional raw CLI tokens.
#' @param capture,rad See [rad_prep()].
#' @return Invisible list with `status`, `stdout`, `stderr`, `command`.
#' @export
rad_modify <- function(layout_key = NULL,
                       whitelist_kit = NULL,
                       add = NULL,
                       remove = FALSE,
                       ...,
                       capture = FALSE,
                       rad = rad_path()) {
  if (is.null(layout_key) == is.null(whitelist_kit)) {
    rlang::abort("Provide exactly one of `layout_key` or `whitelist_kit`.")
  }
  if (is.null(add) == isFALSE(remove)) {
    rlang::abort("Provide exactly one of `add = <path>` or `remove = TRUE`.")
  }
  args <- c(
    "modify",
    .rad_flag("--layout",    layout_key),
    .rad_flag("--whitelist", whitelist_kit),
    .rad_flag("--add",       add),
    if (isTRUE(remove)) "--remove" else character(),
    as.character(c(...))
  )
  .rad_call(args, rad = rad, capture = capture)
}

#' Scan a FASTQ for barcodes and generate/validate whitelists (`rad scan-wl`)
#'
#' Wraps `rad scan-wl`. Required: `adapter_seq`. Provide either `input`
#' (single run) or `batch_csv` (batched runs).
#'
#' @param adapter_seq Adapter/primer sequence to anchor the barcode search.
#' @param input Single-run FASTQ.
#' @param batch_csv CSV with two columns: FASTQ path, output prefix.
#' @param output_prefix Output prefix for generated `.txt` and `.csv` files.
#' @param barcode_length Bases to extract.
#' @param left_margin,right_margin Extra bases on each side.
#' @param max_reads Max reads to process. `NULL` = all.
#' @param max_error Maximum edit-distance ratio. rad default 0.1.
#' @param whitelist Optional whitelist kit key or path.
#' @param threads Threads. rad default auto.
#' @param chunk_size rad default 10000.
#' @param bc1_whitelist,bc2_whitelist Two-part-barcode mode (BC1+BC2).
#'   Setting `bc1_whitelist` activates this mode.
#' @param umi_length UMI length between primer and BC1. rad default 9.
#' @param offset_min,offset_max Min/max extra bases after UMI for BC1
#'   search. rad defaults 0 and 3.
#' @param verbose Logical; verbose logging.
#' @param ... Additional raw CLI tokens.
#' @param capture,rad See [rad_prep()].
#' @return Invisible list with `status`, `stdout`, `stderr`, `command`.
#' @export
rad_scan_wl <- function(adapter_seq,
                        input = NULL,
                        batch_csv = NULL,
                        output_prefix = NULL,
                        barcode_length = NULL,
                        left_margin = NULL,
                        right_margin = NULL,
                        max_reads = NULL,
                        max_error = NULL,
                        whitelist = NULL,
                        threads = NULL,
                        chunk_size = NULL,
                        bc1_whitelist = NULL,
                        bc2_whitelist = NULL,
                        umi_length = NULL,
                        offset_min = NULL,
                        offset_max = NULL,
                        verbose = FALSE,
                        ...,
                        capture = FALSE,
                        rad = rad_path()) {
  if (is.null(input) && is.null(batch_csv)) {
    rlang::abort("Provide either `input` (single run) or `batch_csv` (batched).")
  }
  args <- c(
    "scan-wl",
    .rad_flag("--adapter_seq",   adapter_seq),
    .rad_flag("--input",         input),
    .rad_flag("--batch-csv",     batch_csv),
    .rad_flag("--output-prefix", output_prefix),
    .rad_flag("--barcode-length", barcode_length),
    .rad_flag("--left-margin",   left_margin),
    .rad_flag("--right-margin",  right_margin),
    .rad_flag("--max-reads",     max_reads),
    .rad_flag("--max-error",     max_error),
    .rad_flag("--whitelist",     whitelist),
    .rad_flag("--threads",       threads),
    .rad_flag("--chunk-size",    chunk_size),
    .rad_flag("--bc1-whitelist", bc1_whitelist),
    .rad_flag("--bc2-whitelist", bc2_whitelist),
    .rad_flag("--umi-length",    umi_length),
    .rad_flag("--offset-min",    offset_min),
    .rad_flag("--offset-max",    offset_max),
    if (isTRUE(verbose)) "--verbose" else character(),
    as.character(c(...))
  )
  .rad_call(args, rad = rad, capture = capture)
}

#' Print the rad CLI help for a subcommand
#'
#' Convenience wrapper around `rad help <topic>`.
#'
#' @param topic One of `"prep"`, `"demux"`, `"reformat"`, `"list"`,
#'   `"modify"`, `"scan-wl"`, or `NULL` for top-level help.
#' @param rad Path to the rad binary.
#' @return Invisible list (see [rad_prep()]).
#' @export
rad_help <- function(topic = NULL, rad = rad_path()) {
  args <- if (is.null(topic)) c("--help") else c("help", topic)
  .rad_call(args, rad = rad, capture = TRUE, echo = FALSE, check = FALSE)
}
