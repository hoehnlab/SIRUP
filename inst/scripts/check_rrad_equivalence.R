#!/usr/bin/env Rscript
# A real external-tool comparison. Sourcing this file defines helpers only.

equivalence_compare <- function(expected, actual) {
  expected <- as.data.frame(expected, stringsAsFactors = FALSE)
  actual <- as.data.frame(actual, stringsAsFactors = FALSE)
  if (!nrow(expected) || !nrow(actual)) {
    return(list(pass = FALSE, detail = "Empty checkpoint: equivalence is not established."))
  }
  if (!identical(sort(names(expected)), sort(names(actual))) ||
      anyDuplicated(names(expected)) || anyDuplicated(names(actual))) {
    return(list(pass = FALSE, detail = "Column names differ or are duplicated."))
  }
  canonical <- function(x) {
    x <- x[, sort(names(x)), drop = FALSE]
    x[] <- lapply(x, as.character)
    x <- x[do.call(order, c(unname(x), list(na.last = TRUE, method = "radix"))), , drop = FALSE]
    rownames(x) <- NULL
    x
  }
  expected <- canonical(expected)
  actual <- canonical(actual)
  pass <- identical(expected, actual)
  list(pass = pass, detail = if (pass) "All columns and row multiplicities match; row order ignored."
       else paste(utils::head(all.equal(expected, actual), 3L), collapse = "; "))
}

equivalence_read_fasta <- function(path, limit = Inf) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  ids <- sequences <- character()
  id <- NULL
  sequence <- character()
  add <- function() {
    if (is.null(id)) return()
    if (!length(sequence)) stop("Empty FASTA record: ", id)
    ids <<- c(ids, id)
    sequences <<- c(sequences, paste0(sequence, collapse = ""))
  }
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) {
      add()
      break
    }
    if (startsWith(line, ">")) {
      add()
      if (length(ids) >= limit) break
      id <- substring(line, 2L)
      if (!nzchar(id)) stop("Empty FASTA identifier.")
      sequence <- character()
    } else {
      if (is.null(id)) stop("Input is not FASTA: sequence appears before a header.")
      if (!nzchar(line) || grepl("[[:space:]]", line)) stop("Blank or whitespace-containing FASTA sequence.")
      sequence <- c(sequence, line)
    }
  }
  data.frame(sequence_id = ids, sequence = sequences, stringsAsFactors = FALSE)
}

equivalence_write_fasta <- function(records, path) {
  writeLines(as.vector(rbind(paste0(">", records$sequence_id), records$sequence)), path)
}

# Deliberately independent of parseRadSequenceIds(), restricted to rrad's
# canonical pRESTO export. No identifiers or biological columns are discarded.
equivalence_parse_tags <- function(x) {
  tags <- strsplit(x$sequence_id, "|", fixed = TRUE)
  fields <- lapply(tags, function(parts) {
    kv <- strsplit(parts[-1L], "=", fixed = TRUE)
    if (any(lengths(kv) != 2L)) stop("Expected rrad pRESTO tags: |BARCODE=...|UMI=...")
    keys <- vapply(kv, function(z) toupper(z[[1L]]), character(1))
    if (anyDuplicated(keys)) stop("Duplicate FASTA header tag.")
    stats::setNames(vapply(kv, `[[`, character(1), 2L), keys)
  })
  for (key in c("BARCODE", "UMI")) {
    values <- vapply(fields, function(z) if (key %in% names(z)) z[[key]] else NA_character_, character(1))
    present <- !is.na(values)
    if (any(present & (!nzchar(values) | grepl("[|=,[:space:]]", values)))) stop("Invalid ", key, " tag.")
    x[[tolower(key)]] <- values
  }
  if (anyNA(x$barcode)) stop("Every selected FASTA read needs a |BARCODE=... tag.")
  x$sequence_id <- vapply(tags, `[[`, character(1), 1L)
  if (any(!nzchar(x$sequence_id)) || any(grepl("[[:space:]]", x$sequence_id))) stop("Invalid read identifier.")
  x
}

equivalence_help <- function() {
  cat(paste(
    "Compare 500 rrad-prepared reads through legacy and current SIRAP workflows.",
    "",
    "Rscript inst/scripts/check_rrad_equivalence.R \\",
    "  --input reads.fasta.gz --output-dir NEW_DIRECTORY \\",
    "  --igblast-base /path/igblast --igblast-db-dir /path/raw_databases \\",
    "  --germline-glob '/path/imgt/*.fasta' --muscle /path/muscle \\",
    "  --receptor ig --organism human --threads 1",
    "",
    "Optional: --sirap-repo CHECKOUT (load working source with pkgload)",
    "  --baseline-repo CHECKOUT (default: --sirap-repo or working directory)",
    "  --baseline-ref REF (default: bb1b912ebc97f9daadd0c338d0d30031b1d1e1b3)",
    "  --igblast-exec PATH --aux-file PATH",
    "  --germline-db-v PATH --germline-db-d PATH --germline-db-j PATH",
    "  --c-region-db PATH (IG only; defaults to IG_C in --igblast-db-dir)",
    "",
    "Requires real IgBLAST, pRESTO, Change-O, vsearch and MUSCLE 3 plus databases.",
    "No downloads. No mocks. Outputs must not already exist. Nonzero exit on failure.",
    sep = "\n"
  ), "\n")
}

equivalence_options <- function(args) {
  if (!length(args) || "--help" %in% args) {
    equivalence_help()
    return(NULL)
  }
  allowed <- c("input", "output-dir", "sirap-repo", "baseline-repo", "baseline-ref",
               "igblast-base", "igblast-db-dir", "germline-glob", "muscle", "receptor",
               "organism", "threads", "aux-file", "igblast-exec", "germline-db-v",
               "germline-db-d", "germline-db-j", "c-region-db")
  if (length(args) %% 2L) stop("Options require --name value pairs; see --help.")
  flags <- args[seq(1L, length(args), 2L)]
  if (any(!startsWith(flags, "--"))) stop("Option names must start with --; see --help.")
  keys <- sub("^--", "", flags)
  if (any(!keys %in% allowed) || anyDuplicated(keys)) stop("Unknown or repeated option; see --help.")
  opt <- as.list(stats::setNames(args[seq(2L, length(args), 2L)], keys))
  defaults <- list(receptor = "ig", organism = "human", threads = "1",
                   `baseline-ref` = "bb1b912ebc97f9daadd0c338d0d30031b1d1e1b3")
  opt <- utils::modifyList(defaults, opt)
  for (key in c("input", "output-dir", "igblast-base", "germline-glob", "muscle")) {
    if (is.null(opt[[key]]) || !nzchar(opt[[key]])) stop("Required option: --", key)
  }
  if (!opt$receptor %in% c("ig", "tr")) stop("--receptor must be ig or tr.")
  threads <- suppressWarnings(as.integer(opt$threads))
  if (!grepl("^[1-9][0-9]*$", opt$threads) || is.na(threads)) {
    stop("--threads must be a positive integer no greater than ", .Machine$integer.max, ".")
  }
  if (!grepl("^[[:alnum:]_]+$", opt$organism)) stop("Invalid --organism.")
  opt
}

equivalence_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  opt <- equivalence_options(args)
  if (is.null(opt)) return(invisible(NULL))
  input <- normalizePath(path.expand(opt$input), mustWork = TRUE)
  out <- path.expand(opt[["output-dir"]])
  if (file.exists(out) || dir.exists(out)) stop("Output directory must be new: ", out)
  if (!is.null(opt[["sirap-repo"]])) pkgload::load_all(opt[["sirap-repo"]], quiet = TRUE, export_all = FALSE)
  ns <- asNamespace("SIRAP")
  if (!"input" %in% names(formals(ns$processSpatialSample))) stop("Load the reviewed SIRAP version with --sirap-repo.")
  baseline_repo <- opt[["baseline-repo"]]
  if (is.null(baseline_repo)) baseline_repo <- if (is.null(opt[["sirap-repo"]])) getwd() else opt[["sirap-repo"]]
  git <- function(args) {
    result <- system2("git", shQuote(c("-C", baseline_repo, args)), stdout = TRUE, stderr = TRUE)
    if (!is.null(attr(result, "status"))) stop("Cannot read baseline Git revision: ", paste(result, collapse = "\n"))
    result
  }
  baseline_sha <- git(c("rev-parse", paste0(opt[["baseline-ref"]], "^{commit}")))
  legacy_source <- git(c("show", paste0(baseline_sha, ":R/processing.R")))
  legacy <- new.env(parent = ns)
  eval(parse(text = legacy_source), envir = legacy)
  if (!"input_fasta" %in% names(formals(legacy$processSpatialSample))) stop("Baseline must expose the legacy input_fasta workflow.")
  igblast_base <- normalizePath(path.expand(opt[["igblast-base"]]), mustWork = TRUE)
  igblast_exec <- opt[["igblast-exec"]]
  if (is.null(igblast_exec)) igblast_exec <- file.path(igblast_base, "bin", "igblastn")
  resolve_tool <- function(x) {
    expanded <- path.expand(x)
    found <- if (file.exists(expanded)) normalizePath(expanded) else unname(Sys.which(x))
    if (!nzchar(found) || file.access(found, 1L) != 0L) stop("Required real executable unavailable: ", x)
    found
  }
  tools <- c(igblastn = igblast_exec, muscle = opt$muscle,
             ClusterSets.py = "ClusterSets.py", ParseHeaders.py = "ParseHeaders.py",
             AlignSets.py = "AlignSets.py", BuildConsensus.py = "BuildConsensus.py",
             AssignGenes.py = "AssignGenes.py", MakeDb.py = "MakeDb.py", vsearch = "vsearch")
  tool_paths <- vapply(tools, resolve_tool, character(1))
  germlines <- Sys.glob(path.expand(opt[["germline-glob"]]))
  if (!length(germlines)) stop("--germline-glob matched no files.")
  germlines <- normalizePath(germlines, mustWork = TRUE)
  # Database-path resolution is shared configuration for both arms, so use
  # the package's own resolver instead of a copy that could drift from it.
  dbs <- ns$.sirap_igblast_db_paths(
    igblast_db_dir = opt[["igblast-db-dir"]],
    receptor = opt$receptor,
    germline_db_v = opt[["germline-db-v"]],
    germline_db_d = opt[["germline-db-d"]],
    germline_db_j = opt[["germline-db-j"]],
    c_region_db = opt[["c-region-db"]]
  )
  aux <- if (is.null(opt[["aux-file"]])) NULL else normalizePath(path.expand(opt[["aux-file"]]), mustWork = TRUE)
  dir.create(out, recursive = TRUE)
  out <- normalizePath(out)
  input_info <- file.info(input)
  manifest <- list(started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), options = opt,
                   baseline_sha = baseline_sha, input = input,
                   input_bytes = input_info$size,
                   input_modified_utc = format(input_info$mtime, tz = "UTC", usetz = TRUE),
                   R = R.version.string,
                   SIRAP_version = as.character(utils::packageVersion("SIRAP")),
                   rrad_version = if (requireNamespace("rrad", quietly = TRUE)) as.character(utils::packageVersion("rrad")) else NA_character_,
                   tools = tool_paths, database_paths = dbs,
                   germline_md5 = tools::md5sum(germlines), status = "running")
  if (!is.null(opt[["sirap-repo"]])) {
    manifest$current_git <- system2("git", shQuote(c("-C", opt[["sirap-repo"]], "rev-parse", "HEAD")), stdout = TRUE)
    manifest$current_worktree <- system2("git", shQuote(c("-C", opt[["sirap-repo"]], "status", "--porcelain")), stdout = TRUE)
    manifest$current_source_md5 <- tools::md5sum(file.path(opt[["sirap-repo"]],
      c("R/processing.R", "R/igblast.R", "R/r_commands.R", "inst/scripts/check_rrad_equivalence.R")))
  }
  manifest$tool_versions <- lapply(names(tool_paths), function(key) {
    flag <- if (key %in% c("igblastn", "muscle")) "-version" else "--version"
    suppressWarnings(system2(tool_paths[[key]], flag, stdout = TRUE, stderr = TRUE))
  })
  names(manifest$tool_versions) <- names(tool_paths)
  save_manifest <- function() {
    saveRDS(manifest, file.path(out, "manifest.rds"))
    capture.output(str(manifest), file = file.path(out, "manifest.txt"))
  }
  save_manifest()
  stage <- "select_reads"
  result <- tryCatch({
    reads <- equivalence_read_fasta(input, 500L)
    if (nrow(reads) != 500L) stop("At least 500 complete FASTA records are required; found ", nrow(reads), ".")
    metadata <- equivalence_parse_tags(reads)
    if (anyDuplicated(metadata$sequence_id)) stop("Selected reads have duplicate bare sequence identifiers.")
    subset <- file.path(out, "reads_500.fasta")
    equivalence_write_fasta(reads, subset)
    manifest$subset_md5 <- unname(tools::md5sum(subset))
    manifest$selected_reads <- nrow(reads)
    writeLines(reads$sequence_id, file.path(out, "selected_read_ids.txt"))
    for (arm in c("legacy", "current")) dir.create(file.path(out, arm))
    legacy_raw <- file.path(out, "legacy", "sample_raw.tsv")
    raw_args <- c("-query", subset, "-ig_seqtype", if (opt$receptor == "ig") "Ig" else "TCR",
                  "-organism", opt$organism, "-germline_db_V", dbs$v,
                  "-germline_db_D", dbs$d, "-germline_db_J", dbs$j)
    if (!is.null(dbs$c)) raw_args <- c(raw_args, "-c_region_db", dbs$c)
    if (!is.null(aux)) raw_args <- c(raw_args, "-auxiliary_data", aux)
    raw_args <- c(raw_args, "-outfmt", "19", "-num_threads", opt$threads, "-out", legacy_raw)
    manifest$legacy_raw_command <- c(tool_paths[["igblastn"]], raw_args)
    save_manifest()
    stage <- "legacy_raw_airr"
    # Note: only the final db_pass checkpoint is a genuine old-vs-new
    # comparison. The raw-IgBLAST, tag-parse, and filter stages of this arm
    # are script-side reconstructions that mirror the current pipeline's
    # defaults, because R/igblast.R did not exist at the baseline commit.
    message("Running legacy-workflow arm on 500 reads.")
    log <- file.path(out, "legacy", "raw_igblast.log")
    status <- system2(tool_paths[["igblastn"]], shQuote(raw_args), stdout = log, stderr = log)
    if (status != 0L) stop("Legacy raw IgBLAST failed; see ", log)
    old_airr <- data.table::fread(legacy_raw, data.table = FALSE)
    if (nrow(old_airr) != 500L) stop("Legacy raw AIRR must contain 500 reads; found ", nrow(old_airr), ".")
    stage <- "legacy_filtered_airr"
    old_airr <- equivalence_parse_tags(old_airr)
    if ("complete_vdj" %in% names(old_airr)) old_airr <- old_airr[old_airr$complete_vdj %in% c(TRUE, "TRUE", "T"), , drop = FALSE]
    old_filtered <- legacy$filterAirr(old_airr, min_consensus_count = NULL, spatial_id = "barcode",
                                      add_spatial_sequence_id = TRUE, spatial_id_key = "barcode")
    if (!nrow(old_filtered)) stop("Legacy filtered AIRR is empty; choose a receptor-appropriate subset.")
    legacy_filtered <- file.path(out, "legacy", "sample_filtered.tsv")
    data.table::fwrite(old_filtered, legacy_filtered, sep = "\t")
    legacy_fasta <- file.path(out, "legacy", "sample_spatial.fasta")
    equivalence_write_fasta(data.frame(sequence_id = old_filtered$spatial_sequence_id,
                                      sequence = old_filtered$sequence_alignment), legacy_fasta)
    common <- list(receptor = opt$receptor, organism = opt$organism,
                   igblast_base = igblast_base, igblast_exec = tool_paths[["igblastn"]],
                   muscle_exec = tool_paths[["muscle"]], nproc = as.integer(opt$threads),
                   cluster_field = "barcode", consensus_min_count = 1)
    # The old workflow expects a shell glob. Quote each resolved file explicitly
    # so paths containing spaces are safe and both arms see identical references.
    stage <- "legacy_consensus_pipeline"
    legacy_paths <- do.call(legacy$processSpatialSample, c(list(input_fasta = legacy_fasta,
      output_dir = file.path(out, "legacy"), germline_glob = paste(shQuote(germlines), collapse = " ")), common))
    message("Running current processSpatialSample on the identical 500 reads.")
    stage <- "current_pipeline"
    current_paths <- do.call(ns$processSpatialSample, c(list(input = subset,
      output_dir = file.path(out, "current"), sample_name = "sample",
      germline_glob = opt[["germline-glob"]], aux_file = aux,
      germline_db_v = dbs$v, germline_db_d = dbs$d, germline_db_j = dbs$j, c_region_db = dbs$c,
      filter_args = list(min_consensus_count = NULL)), common))
    checkpoints <- list(raw_airr = c(legacy_raw, current_paths$raw_airr),
                        filtered_airr = c(legacy_filtered, current_paths$filtered_airr),
                        spatial_fasta = c(legacy_fasta, current_paths$spatial_fasta),
                        db_pass = c(legacy_paths$db_pass, current_paths$db_pass))
    stage <- "compare_checkpoints"
    results <- lapply(names(checkpoints), function(key) {
      paths <- checkpoints[[key]]
      if (length(paths) != 2L || any(!file.exists(paths))) stop("Missing ", key, " output.")
      tables <- lapply(paths, function(path) {
        if (key == "spatial_fasta") equivalence_read_fasta(path)
        else data.table::fread(path, sep = "\t", data.table = FALSE, colClasses = "character", na.strings = NULL)
      })
      checked <- equivalence_compare(tables[[1L]], tables[[2L]])
      data.frame(checkpoint = key, pass = checked$pass, legacy_rows = nrow(tables[[1L]]),
                 current_rows = nrow(tables[[2L]]), detail = checked$detail)
    })
    results <- do.call(rbind, results)
    data.table::fwrite(results, file.path(out, "equivalence.tsv"), sep = "\t")
    manifest$checkpoints <- checkpoints
    print(results[, 1:4], row.names = FALSE)
    if (!all(results$pass)) stop("Equivalence failed; inspect equivalence.tsv and both output directories.")
    manifest$status <- "passed"
    invisible(results)
  }, error = function(e) {
    manifest$status <<- "failed"
    manifest$error <<- conditionMessage(e)
    manifest$failed_stage <<- stage
    if (!file.exists(file.path(out, "equivalence.tsv"))) {
      data.table::fwrite(data.frame(checkpoint = stage, pass = FALSE,
        legacy_rows = NA_integer_, current_rows = NA_integer_, detail = conditionMessage(e)),
        file.path(out, "equivalence.tsv"), sep = "\t")
    }
    NULL
  })
  manifest$finished_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  save_manifest()
  if (manifest$status != "passed") stop(manifest$error, call. = FALSE)
  message("500-read equivalence passed. Evidence: ", out)
  invisible(result)
}

if (sys.nframe() == 0L) equivalence_main()
