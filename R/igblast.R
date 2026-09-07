# Direct IgBLAST support ----------------------------------------------------

.sirap_igblast_abort <- function(message, class = "sirap_igblast_error") {
  rlang::abort(message, class = class)
}

.sirap_igblast_scalar_character <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    .sirap_igblast_abort(
      paste0("`", name, "` must be one non-empty character value.")
    )
  }
  enc2utf8(x)
}

.sirap_igblast_positive_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x != trunc(x) || x < 1 || x > .Machine$integer.max) {
    .sirap_igblast_abort(
      paste0("`", name, "` must be one positive whole number.")
    )
  }
  as.integer(x)
}

.sirap_igblast_existing_file <- function(path, name) {
  path <- path.expand(.sirap_igblast_scalar_character(path, name))
  if (!file.exists(path) || isTRUE(file.info(path)$isdir)) {
    .sirap_igblast_abort(
      paste0("`", name, "` does not identify an existing file: ", path, ".")
    )
  }
  normalizePath(path, mustWork = TRUE)
}

.sirap_igblast_output_path <- function(path, name, extensions = NULL) {
  path <- path.expand(.sirap_igblast_scalar_character(path, name))
  if (dir.exists(path)) {
    .sirap_igblast_abort(
      paste0("`", name, "` identifies a directory, not a file: ", path, ".")
    )
  }
  if (!is.null(extensions)) {
    extension <- tolower(tools::file_ext(path))
    if (!(extension %in% extensions)) {
      .sirap_igblast_abort(paste0(
        "`", name, "` must end in ",
        paste0(".", extensions, collapse = ", "), "."
      ))
    }
  }
  parent <- dirname(path)
  if (file.exists(parent) && !dir.exists(parent)) {
    .sirap_igblast_abort(
      paste0("The parent of `", name, "` is not a directory: ", parent, ".")
    )
  }
  if (!dir.exists(parent) &&
      !dir.create(parent, recursive = TRUE, showWarnings = FALSE)) {
    .sirap_igblast_abort(
      paste0("Could not create the directory for `", name, "`: ", parent, ".")
    )
  }
  file.path(normalizePath(parent, mustWork = TRUE), basename(path))
}

.sirap_resolve_executable <- function(exec, name) {
  exec <- path.expand(.sirap_igblast_scalar_character(exec, name))
  resolved <- if (file.exists(exec)) exec else unname(Sys.which(exec))
  if (!length(resolved) || !nzchar(resolved) || !file.exists(resolved) ||
      isTRUE(file.info(resolved)$isdir) || file.access(resolved, mode = 1L) != 0L) {
    .sirap_igblast_abort(paste0(
      "External tool not found or not executable: ", name, " (", exec,
      "). Install it, add it to PATH, or pass its full path."
    ), class = "sirap_igblast_tool_error")
  }
  normalizePath(resolved, mustWork = TRUE)
}

# Resolve and validate external commands. The resolved paths are returned so
# the exact executable used can be recorded in result provenance.
.sirap_check_tools <- function(execs) {
  if (!is.character(execs) || !length(execs) || anyNA(execs)) {
    .sirap_igblast_abort("`execs` must be a non-empty character vector.")
  }
  labels <- names(execs)
  if (is.null(labels)) {
    labels <- execs
  }
  empty_labels <- is.na(labels) | !nzchar(labels)
  labels[empty_labels] <- execs[empty_labels]
  resolved <- vapply(seq_along(execs), function(index) {
    .sirap_resolve_executable(execs[[index]], labels[[index]])
  }, character(1))
  names(resolved) <- labels
  resolved
}

.sirap_tail_lines <- function(path, n = 10L) {
  if (!file.exists(path) || n < 1L) {
    return(character())
  }
  connection <- file(path, open = "rt", encoding = "UTF-8")
  on.exit(close(connection), add = TRUE)
  result <- character()
  repeat {
    lines <- readLines(connection, n = 1000L, warn = FALSE)
    if (!length(lines)) {
      break
    }
    result <- utils::tail(c(result, lines), n)
  }
  result
}

.sirap_command_env <- function(env) {
  if (is.null(env) || !length(env)) {
    return(character())
  }
  if (!is.character(env) || anyNA(env) || is.null(names(env)) ||
      anyNA(names(env)) || any(!nzchar(names(env))) ||
      anyDuplicated(names(env)) ||
      any(!grepl("^[A-Za-z_][A-Za-z0-9_]*$", names(env)))) {
    .sirap_igblast_abort(paste0(
      "`env` must be NULL or a uniquely named character vector with valid ",
      "environment-variable names."
    ))
  }
  paste0(
    names(env), "=",
    vapply(unname(env), shQuote, character(1), type = "sh")
  )
}

# Execute a command without invoking an intermediate shell. Arguments are
# individually quoted for system2(), and both output streams are captured in
# one log file when requested.
.sirap_run_command <- function(exec, args, log = NULL, env = NULL) {
  exec <- .sirap_igblast_scalar_character(exec, "exec")
  if (!is.character(args) || anyNA(args)) {
    .sirap_igblast_abort("`args` must be a character vector without missing values.")
  }
  command_env <- .sirap_command_env(env)
  if (!is.null(log)) {
    log <- .sirap_igblast_scalar_character(log, "log")
  }

  status <- tryCatch(
    system2(
      command = exec,
      args = vapply(args, shQuote, character(1)),
      stdout = if (is.null(log)) "" else log,
      stderr = if (is.null(log)) "" else log,
      env = command_env
    ),
    error = function(error) {
      .sirap_igblast_abort(
        paste0(
          "Could not run external command ", exec, ": ",
          conditionMessage(error)
        ),
        class = "sirap_igblast_command_error"
      )
    }
  )

  if (!identical(as.integer(status), 0L)) {
    detail <- ""
    if (!is.null(log) && file.exists(log)) {
      log_tail <- .sirap_tail_lines(log, 10L)
      if (length(log_tail)) {
        detail <- paste0("\nLog tail (", log, "):\n", paste(log_tail, collapse = "\n"))
      }
    }
    .sirap_igblast_abort(
      paste0(
        "External command failed with status ", status, ": ", exec, ".",
        detail
      ),
      class = "sirap_igblast_command_error"
    )
  }
  invisible(as.integer(status))
}

# Stream over a plain FASTA so malformed or accidentally compressed/FASTQ
# input is rejected before a potentially expensive IgBLAST run.
.sirap_validate_igblast_fasta <- function(path) {
  path <- .sirap_igblast_existing_file(path, "input_fasta")
  if (isTRUE(file.info(path)$size == 0)) {
    .sirap_igblast_abort("`input_fasta` is empty.")
  }
  magic <- readBin(path, what = "raw", n = 2L)
  if (length(magic) == 2L && identical(as.integer(magic), c(31L, 139L))) {
    .sirap_igblast_abort(
      "`input_fasta` is gzip-compressed. Convert it to a plain FASTA before calling `runIgblast()`."
    )
  }

  connection <- file(path, open = "rt", encoding = "UTF-8")
  on.exit(close(connection), add = TRUE)
  query_count <- 0
  base_count <- 0
  line_number <- 0
  have_header <- FALSE
  have_sequence <- FALSE

  repeat {
    lines <- readLines(connection, n = 10000L, warn = FALSE)
    if (!length(lines)) {
      break
    }
    for (line in lines) {
      line_number <- line_number + 1
      if (!nzchar(line)) {
        next
      }
      if (startsWith(line, ">")) {
        if (have_header && !have_sequence) {
          .sirap_igblast_abort(paste0(
            "FASTA record ", query_count, " has no sequence before line ",
            line_number, "."
          ))
        }
        header <- substring(line, 2L)
        identifier <- sub("[[:space:]].*$", "", header)
        if (!nzchar(identifier)) {
          .sirap_igblast_abort(
            paste0("FASTA header on line ", line_number, " has no identifier.")
          )
        }
        query_count <- query_count + 1
        have_header <- TRUE
        have_sequence <- FALSE
        next
      }
      if (!have_header) {
        .sirap_igblast_abort(
          paste0("`input_fasta` does not begin with a FASTA header (line ", line_number, ").")
        )
      }
      if (grepl("[[:space:]]", line) ||
          !grepl("^[ACGTURYSWKMBDHVNXacgturyswkmbdhvnx.*-]+$", line)) {
        .sirap_igblast_abort(paste0(
          "Invalid nucleotide sequence on FASTA line ", line_number, "."
        ))
      }
      have_sequence <- TRUE
      base_count <- base_count + nchar(line, type = "bytes")
    }
  }

  if (!have_header) {
    .sirap_igblast_abort("`input_fasta` contains no FASTA records.")
  }
  if (!have_sequence) {
    .sirap_igblast_abort(
      paste0("FASTA record ", query_count, " has no sequence at end of file.")
    )
  }
  list(path = path, queries = query_count, bases = base_count)
}

.sirap_igblast_chunk_bases <- function(chunk_bases) {
  if (is.null(chunk_bases)) {
    return(NULL)
  }
  .sirap_igblast_positive_integer(chunk_bases, "chunk_bases")
}

# Choose a conservative outer-process count. Automatic planning requires at
# least two records per process on average, avoiding parallel startup for the
# tiny FASTA files that are common in examples and resumed workflows.
.sirap_igblast_process_plan <- function(fasta,
                                        num_threads,
                                        max_processes,
                                        chunk_bases = NULL) {
  num_threads <- .sirap_igblast_positive_integer(num_threads, "num_threads")
  max_processes <- .sirap_igblast_positive_integer(
    max_processes, "max_processes"
  )
  chunk_bases <- .sirap_igblast_chunk_bases(chunk_bases)

  process_limit <- min(num_threads, max_processes)
  minimum_parallel_bases <- 100000
  force_single_process <- max_processes == 1L ||
    (is.null(chunk_bases) && fasta$bases < minimum_parallel_bases)
  if (force_single_process) {
    desired_chunks <- 1L
    target_bases <- fasta$bases
    max_chunks <- 1L
  } else if (is.null(chunk_bases)) {
    record_limit <- max(1, floor(fasta$queries / 2))
    desired_chunks <- as.integer(min(process_limit, record_limit))
    target_bases <- max(1, ceiling(fasta$bases / desired_chunks))
    max_chunks <- desired_chunks
  } else {
    desired_chunks <- max(1, ceiling(fasta$bases / chunk_bases))
    target_bases <- chunk_bases
    max_chunks <- min(desired_chunks, fasta$queries)
  }

  list(
    thread_budget = num_threads,
    max_processes = max_processes,
    desired_chunks = desired_chunks,
    target_bases = target_bases,
    max_chunks = max_chunks,
    minimum_parallel_bases = minimum_parallel_bases,
    requested_chunk_bases = chunk_bases
  )
}

# Split a validated FASTA without materializing it. Rotation occurs only at a
# header, so wrapped sequences and very long individual records stay intact.
.sirap_split_igblast_fasta <- function(fasta, directory, plan) {
  if (plan$desired_chunks <= 1L) {
    return(list(list(
      index = 1L,
      path = fasta$path,
      queries = fasta$queries,
      bases = fasta$bases,
      temporary = FALSE
    )))
  }
  if (!dir.exists(directory) &&
      !dir.create(directory, recursive = TRUE, showWarnings = FALSE)) {
    .sirap_igblast_abort(paste0(
      "Could not create the IgBLAST chunk directory: ", directory, "."
    ))
  }

  input <- file(fasta$path, open = "rt", encoding = "UTF-8")
  on.exit(close(input), add = TRUE)
  chunk_index <- 1L
  chunk_paths <- file.path(directory, sprintf("query-%05d.fasta", chunk_index))
  output <- file(chunk_paths[[chunk_index]], open = "wt", encoding = "UTF-8")
  output_open <- TRUE
  on.exit(if (output_open) close(output), add = TRUE)
  chunk_queries <- 0
  chunk_base_counts <- 0
  have_record <- FALSE
  queries_seen <- 0
  bases_seen <- 0

  repeat {
    lines <- readLines(input, n = 10000L, warn = FALSE)
    if (!length(lines)) {
      break
    }
    segment_start <- 1L
    for (line_index in seq_along(lines)) {
      line <- lines[[line_index]]
      if (startsWith(line, ">")) {
        remaining_chunks <- plan$max_chunks - chunk_index
        remaining_records <- fasta$queries - queries_seen
        reached_boundary <- bases_seen >= ceiling(
          fasta$bases * chunk_index / plan$max_chunks
        )
        must_leave_one_record_per_chunk <- remaining_chunks > 0 &&
          remaining_records <= remaining_chunks
        rotate <- have_record &&
          chunk_index < plan$max_chunks &&
          (reached_boundary || must_leave_one_record_per_chunk)
        if (rotate) {
          if (line_index > segment_start) {
            writeLines(
              lines[segment_start:(line_index - 1L)], output,
              sep = "\n", useBytes = TRUE
            )
          }
          close(output)
          output_open <- FALSE
          chunk_index <- chunk_index + 1L
          chunk_paths[[chunk_index]] <- file.path(
            directory, sprintf("query-%05d.fasta", chunk_index)
          )
          chunk_queries[[chunk_index]] <- 0
          chunk_base_counts[[chunk_index]] <- 0
          output <- file(
            chunk_paths[[chunk_index]], open = "wt", encoding = "UTF-8"
          )
          output_open <- TRUE
          segment_start <- line_index
        }
        chunk_queries[[chunk_index]] <-
          chunk_queries[[chunk_index]] + 1
        queries_seen <- queries_seen + 1
        have_record <- TRUE
      } else if (nzchar(line)) {
        chunk_base_counts[[chunk_index]] <-
          chunk_base_counts[[chunk_index]] + nchar(line, type = "bytes")
        bases_seen <- bases_seen + nchar(line, type = "bytes")
      }
    }
    if (segment_start <= length(lines)) {
      writeLines(
        lines[segment_start:length(lines)], output,
        sep = "\n", useBytes = TRUE
      )
    }
  }
  close(output)
  output_open <- FALSE

  lapply(seq_along(chunk_paths), function(index) {
    list(
      index = as.integer(index),
      path = normalizePath(chunk_paths[[index]], mustWork = TRUE),
      queries = chunk_queries[[index]],
      bases = chunk_base_counts[[index]],
      temporary = TRUE
    )
  })
}

.sirap_paths_equivalent <- function(left, right) {
  if (!file.exists(left) || !file.exists(right)) {
    return(FALSE)
  }
  identical(
    normalizePath(left, mustWork = TRUE),
    normalizePath(right, mustWork = TRUE)
  )
}

.sirap_publish_igblast_artifact <- function(staged, destination, label) {
  if (!file.rename(staged, destination)) {
    .sirap_igblast_abort(paste0(
      "Could not atomically publish the IgBLAST ", label, " artifact to: ",
      destination, "."
    ))
  }
  normalizePath(destination, mustWork = TRUE)
}

.sirap_igblast_db_prefix_exists <- function(prefix) {
  prefix <- path.expand(prefix)
  database_file_ok <- function(path) {
    if (!file.exists(path)) {
      return(FALSE)
    }
    info <- file.info(path)
    !isTRUE(info$isdir) && !is.na(info$size) && info$size > 0 &&
      file.access(path, mode = 4L) == 0L
  }
  if (database_file_ok(paste0(prefix, ".nal"))) {
    return(TRUE)
  }
  core <- paste0(prefix, c(".nhr", ".nin", ".nsq"))
  if (all(vapply(core, database_file_ok, logical(1)))) {
    return(TRUE)
  }
  parent <- dirname(prefix)
  stem <- basename(prefix)
  files <- if (dir.exists(parent)) {
    list.files(parent, all.files = TRUE, full.names = TRUE, no.. = TRUE)
  } else {
    character()
  }
  filenames <- basename(files)
  volume_files <- files[
    startsWith(filenames, paste0(stem, ".")) &
      grepl("\\.(nhr|nin|nsq)$", filenames)
  ]
  if (!length(volume_files)) {
    return(FALSE)
  }
  volume_prefixes <- unique(sub("\\.(nhr|nin|nsq)$", "", volume_files))
  numeric_volume <- vapply(volume_prefixes, function(volume_prefix) {
    suffix <- substring(volume_prefix, nchar(prefix) + 2L)
    nzchar(suffix) && grepl("^[0-9]+$", suffix)
  }, logical(1))
  volume_prefixes <- volume_prefixes[numeric_volume]
  length(volume_prefixes) > 0L && all(vapply(volume_prefixes, function(volume_prefix) {
    files <- paste0(volume_prefix, c(".nhr", ".nin", ".nsq"))
    all(vapply(files, database_file_ok, logical(1)))
  }, logical(1)))
}

.sirap_validate_igblast_db_prefix <- function(prefix, segment) {
  prefix <- path.expand(.sirap_igblast_scalar_character(prefix, segment))
  if (grepl("\\.(nal|nhr|nin|nsq|ndb|not|ntf|nto)$", prefix, ignore.case = TRUE)) {
    .sirap_igblast_abort(paste0(
      "`", segment, "` must be a BLAST database prefix, without a database-file extension: ",
      prefix, "."
    ))
  }
  parent <- dirname(prefix)
  if (!dir.exists(parent)) {
    .sirap_igblast_abort(
      paste0("The directory for `", segment, "` does not exist: ", parent, ".")
    )
  }
  prefix <- file.path(normalizePath(parent, mustWork = TRUE), basename(prefix))
  if (!.sirap_igblast_db_prefix_exists(prefix)) {
    .sirap_igblast_abort(paste0(
      "No nucleotide BLAST database was found for `", segment, "` prefix ",
      prefix, ". Expected .nhr/.nin/.nsq files (or a .nal alias)."
    ))
  }
  prefix
}

.sirap_igblast_base <- function(igblast_base) {
  if (is.null(igblast_base)) {
    return(NULL)
  }
  igblast_base <- path.expand(
    .sirap_igblast_scalar_character(igblast_base, "igblast_base")
  )
  if (!dir.exists(igblast_base)) {
    .sirap_igblast_abort(paste0(
      "`igblast_base` does not identify an IgBLAST installation directory: ",
      igblast_base, "."
    ))
  }
  normalizePath(igblast_base, mustWork = TRUE)
}

# Resolve IgBLAST database prefixes. V, D, and J are required. The legacy
# `igblast_db_dir` convention uses IG_V/IG_D/IG_J names. An installation root
# uses Change-O's imgt_<organism>_<receptor>_<segment> names beneath database/.
# An explicit C-region prefix is used; otherwise the matching conventional C
# prefix is included only when it exists.
.sirap_igblast_db_paths <- function(igblast_db_dir = NULL,
                                    receptor = c("ig", "tr"),
                                    germline_db_v = NULL,
                                    germline_db_d = NULL,
                                    germline_db_j = NULL,
                                    c_region_db = NULL,
                                    validate = TRUE,
                                    igblast_base = NULL,
                                    organism = "human") {
  receptor <- match.arg(receptor)
  organism <- .sirap_igblast_scalar_character(organism, "organism")
  igblast_base <- .sirap_igblast_base(igblast_base)
  if (!is.logical(validate) || length(validate) != 1L || is.na(validate)) {
    .sirap_igblast_abort("`validate` must be TRUE or FALSE.")
  }
  if (!is.null(igblast_db_dir)) {
    igblast_db_dir <- path.expand(
      .sirap_igblast_scalar_character(igblast_db_dir, "igblast_db_dir")
    )
    if (!dir.exists(igblast_db_dir)) {
      .sirap_igblast_abort(
        paste0("`igblast_db_dir` does not identify a directory: ", igblast_db_dir, ".")
      )
    }
    igblast_db_dir <- normalizePath(igblast_db_dir, mustWork = TRUE)
  }

  legacy_prefix <- toupper(receptor)
  from_legacy_dir <- function(segment) {
    if (is.null(igblast_db_dir)) {
      NULL
    } else {
      file.path(igblast_db_dir, paste0(legacy_prefix, "_", toupper(segment)))
    }
  }
  from_installation <- function(segment) {
    if (is.null(igblast_base)) {
      NULL
    } else {
      file.path(
        igblast_base,
        "database",
        paste0("imgt_", organism, "_", receptor, "_", tolower(segment))
      )
    }
  }
  choose <- function(override, segment) {
    if (!is.null(override)) {
      override
    } else if (!is.null(igblast_db_dir)) {
      from_legacy_dir(segment)
    } else {
      from_installation(segment)
    }
  }

  dbs <- list(
    v = choose(germline_db_v, "V"),
    d = choose(germline_db_d, "D"),
    j = choose(germline_db_j, "J"),
    c = NULL
  )
  missing <- names(dbs)[names(dbs) %in% c("v", "d", "j") &
                         vapply(dbs, is.null, logical(1))]
  if (length(missing)) {
    .sirap_igblast_abort(paste0(
      "Cannot resolve germline database prefixes for ",
      paste(toupper(missing), collapse = ", "),
      ". Provide `igblast_base`, `igblast_db_dir`, or each per-segment prefix."
    ))
  }

  if (!is.null(c_region_db)) {
    dbs$c <- c_region_db
  } else {
    candidate <- if (!is.null(igblast_db_dir)) {
      from_legacy_dir("C")
    } else {
      from_installation("C")
    }
    if (!is.null(candidate) && .sirap_igblast_db_prefix_exists(candidate)) {
      dbs$c <- candidate
    }
  }

  if (validate) {
    dbs$v <- .sirap_validate_igblast_db_prefix(dbs$v, "germline_db_v")
    dbs$d <- .sirap_validate_igblast_db_prefix(dbs$d, "germline_db_d")
    dbs$j <- .sirap_validate_igblast_db_prefix(dbs$j, "germline_db_j")
    if (!is.null(dbs$c)) {
      dbs$c <- .sirap_validate_igblast_db_prefix(dbs$c, "c_region_db")
    }
  } else {
    dbs <- lapply(dbs, function(value) {
      if (is.null(value)) NULL else path.expand(value)
    })
  }
  dbs
}

.sirap_igblast_aux_file <- function(aux_file, igblast_base, organism) {
  if (is.null(aux_file)) {
    if (is.null(igblast_base)) {
      return(NULL)
    }
    aux_file <- file.path(
      igblast_base, "optional_file", paste0(organism, "_gl.aux")
    )
  }
  .sirap_igblast_existing_file(aux_file, "aux_file")
}

.sirap_igblast_extra_args <- function(extra_args) {
  if (is.null(extra_args)) {
    return(character())
  }
  if (!is.character(extra_args) || anyNA(extra_args) || any(!nzchar(extra_args))) {
    .sirap_igblast_abort(
      "`extra_args` must be NULL or a character vector of non-empty arguments."
    )
  }
  protected <- c(
    "-query", "-out", "-outfmt", "-ig_seqtype", "-organism",
    "-germline_db_V", "-germline_db_D", "-germline_db_J",
    "-c_region_db", "-auxiliary_data", "-domain_system", "-num_threads"
  )
  option_names <- sub("=.*$", "", extra_args)
  conflicting <- unique(option_names[option_names %in% protected])
  if (length(conflicting)) {
    .sirap_igblast_abort(paste0(
      "`extra_args` cannot override managed IgBLAST options: ",
      paste(conflicting, collapse = ", "), "."
    ))
  }
  extra_args
}

# Build a complete, deterministic igblastn command for AIRR output.
.sirap_igblast_args <- function(input_fasta,
                                output_tsv,
                                receptor = c("ig", "tr"),
                                organism = "human",
                                dbs,
                                aux_file = NULL,
                                num_threads = 8L,
                                extra_args = NULL) {
  receptor <- match.arg(receptor)
  organism <- .sirap_igblast_scalar_character(organism, "organism")
  num_threads <- .sirap_igblast_positive_integer(num_threads, "num_threads")
  extra_args <- .sirap_igblast_extra_args(extra_args)
  required_dbs <- c("v", "d", "j")
  if (!is.list(dbs) || !all(required_dbs %in% names(dbs)) ||
      any(vapply(dbs[required_dbs], is.null, logical(1)))) {
    .sirap_igblast_abort("`dbs` must contain non-NULL V, D, and J prefixes.")
  }

  args <- c(
    "-query", input_fasta,
    "-ig_seqtype", if (receptor == "ig") "Ig" else "TCR",
    "-organism", organism,
    "-germline_db_V", dbs$v,
    "-germline_db_D", dbs$d,
    "-germline_db_J", dbs$j
  )
  if (!is.null(dbs$c)) {
    args <- c(args, "-c_region_db", dbs$c)
  }
  if (!is.null(aux_file)) {
    args <- c(args, "-auxiliary_data", aux_file)
  }
  c(
    args,
    "-domain_system", "imgt",
    "-outfmt", "19",
    "-num_threads", as.character(num_threads),
    "-out", output_tsv,
    extra_args
  )
}

.sirap_inspect_airr_output <- function(path) {
  path <- .sirap_igblast_existing_file(path, "IgBLAST AIRR output")
  if (isTRUE(file.info(path)$size == 0)) {
    .sirap_igblast_abort("IgBLAST created an empty AIRR output file.")
  }
  connection <- file(path, open = "rt", encoding = "UTF-8")
  on.exit(close(connection), add = TRUE)
  header_line <- readLines(connection, n = 1L, warn = FALSE)
  if (!length(header_line)) {
    .sirap_igblast_abort("IgBLAST created an empty AIRR output file.")
  }
  fields <- strsplit(sub("^\ufeff", "", header_line), "\t", fixed = TRUE)[[1L]]
  required <- c("sequence_id", "v_call", "d_call", "j_call")
  missing <- setdiff(required, fields)
  if (length(missing)) {
    .sirap_igblast_abort(paste0(
      "IgBLAST output is missing required AIRR fields: ",
      paste(missing, collapse = ", "), "."
    ))
  }
  rows <- 0
  annotated_rows <- 0
  call_columns <- match(c("v_call", "d_call", "j_call"), fields)
  sequence_id_column <- match("sequence_id", fields)
  repeat {
    lines <- readLines(connection, n = 10000L, warn = FALSE)
    if (!length(lines)) {
      break
    }
    lines <- lines[nzchar(lines)]
    if (length(lines)) {
      # Appending a sentinel preserves trailing empty TSV fields that base
      # strsplit() would otherwise discard for a no-hit AIRR row.
      values <- lapply(strsplit(paste0(lines, "."), "\t", fixed = TRUE),
                       function(row) {
        final <- length(row)
        row[[final]] <- substr(row[[final]], 1L, nchar(row[[final]]) - 1L)
        row
      })
      widths <- lengths(values)
      malformed <- which(widths != length(fields))
      if (length(malformed)) {
        .sirap_igblast_abort(paste0(
          "IgBLAST AIRR data row ", rows + malformed[[1L]], " has ",
          widths[[malformed[[1L]]]], " fields; the header has ",
          length(fields), "."
        ))
      }
      empty_id <- which(!vapply(values, function(row) {
        nzchar(trimws(row[[sequence_id_column]]))
      }, logical(1)))
      if (length(empty_id)) {
        .sirap_igblast_abort(paste0(
          "IgBLAST AIRR data row ", rows + empty_id[[1L]],
          " has an empty sequence_id."
        ))
      }
      annotated_rows <- annotated_rows + sum(vapply(values, function(row) {
        calls <- row[call_columns]
        any(!is.na(calls) & nzchar(trimws(calls)))
      }, logical(1)))
      rows <- rows + length(lines)
    }
  }
  list(
    path = path,
    header = sub("^\ufeff", "", header_line),
    fields = fields,
    rows = rows,
    annotated_rows = annotated_rows,
    bytes = unname(file.info(path)$size)
  )
}

.sirap_merge_igblast_airr <- function(results, destination) {
  expected_header <- results[[1L]]$airr$header
  output <- file(destination, open = "wt", encoding = "UTF-8")
  on.exit(close(output), add = TRUE)
  writeLines(expected_header, output, sep = "\n", useBytes = TRUE)

  for (result in results) {
    if (!identical(result$airr$header, expected_header)) {
      .sirap_igblast_abort(paste0(
        "IgBLAST chunk ", result$index,
        " returned an AIRR schema that differs from chunk 1."
      ))
    }
    input <- file(result$airr$path, open = "rt", encoding = "UTF-8")
    readLines(input, n = 1L, warn = FALSE)
    repeat {
      lines <- readLines(input, n = 10000L, warn = FALSE)
      if (!length(lines)) {
        break
      }
      lines <- lines[nzchar(lines)]
      if (length(lines)) {
        writeLines(lines, output, sep = "\n", useBytes = TRUE)
      }
    }
    close(input)
  }
  invisible(destination)
}

.sirap_merge_igblast_logs <- function(results, destination) {
  if (length(results) == 1L) {
    if (!file.copy(results[[1L]]$log, destination, overwrite = TRUE)) {
      .sirap_igblast_abort("Could not stage the IgBLAST log artifact.")
    }
    return(invisible(destination))
  }

  output <- file(destination, open = "wt", encoding = "UTF-8")
  on.exit(close(output), add = TRUE)
  total <- length(results)
  for (result in results) {
    writeLines(paste0(
      "===== IgBLAST chunk ", result$index, "/", total,
      "; queries=", result$queries,
      "; bases=", result$bases,
      "; threads=", result$threads, " ====="
    ), output)
    input <- file(result$log, open = "rt", encoding = "UTF-8")
    repeat {
      lines <- readLines(input, n = 1000L, warn = FALSE)
      if (!length(lines)) {
        break
      }
      writeLines(lines, output, sep = "\n", useBytes = TRUE)
    }
    close(input)
  }
  invisible(destination)
}

.sirap_run_igblast_job <- function(job, executable, command_env) {
  tryCatch({
    if (!file.create(job$log)) {
      .sirap_igblast_abort(paste0(
        "Could not create the log for IgBLAST chunk ", job$index, "."
      ))
    }
    .sirap_run_command(
      executable, job$args, log = job$log, env = command_env
    )
    airr <- .sirap_inspect_airr_output(job$output)
    if (airr$rows > job$queries) {
      .sirap_igblast_abort(paste0(
        "IgBLAST chunk ", job$index, " returned ", airr$rows,
        " AIRR rows for ", job$queries, " FASTA queries. A direct AIRR ",
        "pass cannot contain more data rows than input queries."
      ))
    }
    list(
      ok = TRUE,
      index = job$index,
      queries = job$queries,
      bases = job$bases,
      threads = job$threads,
      airr = airr,
      log = job$log
    )
  }, error = function(error) {
    list(
      ok = FALSE,
      index = job$index,
      message = conditionMessage(error),
      classes = class(error)
    )
  })
}

.sirap_run_igblast_jobs <- function(jobs, executable, command_env,
                                    processes) {
  run_one <- function(job) {
    .sirap_run_igblast_job(job, executable, command_env)
  }
  results <- if (processes <= 1L || length(jobs) <= 1L) {
    lapply(jobs, run_one)
  } else {
    parallel::mclapply(
      jobs,
      run_one,
      mc.cores = processes,
      mc.preschedule = FALSE,
      mc.set.seed = FALSE,
      mc.cleanup = TRUE
    )
  }

  unusable <- which(vapply(results, function(result) {
    !is.list(result) || is.null(result$ok)
  }, logical(1)))
  if (length(unusable)) {
    .sirap_igblast_abort(
      paste0("IgBLAST worker ", unusable[[1L]], " ended unexpectedly."),
      class = "sirap_igblast_command_error"
    )
  }
  failures <- which(!vapply(results, `[[`, logical(1), "ok"))
  if (length(failures)) {
    failure <- results[[failures[[1L]]]]
    error_class <- if ("sirap_igblast_command_error" %in% failure$classes) {
      "sirap_igblast_command_error"
    } else {
      "sirap_igblast_error"
    }
    .sirap_igblast_abort(
      paste0(
        "IgBLAST chunk ", failure$index, "/", length(jobs),
        " failed: ", failure$message
      ),
      class = error_class
    )
  }
  results
}

# Read a completed AIRR TSV when an in-memory table is wanted by the caller.
.sirap_read_airr <- function(path) {
  path <- .sirap_igblast_existing_file(path, "path")
  data.table::fread(path, sep = "\t", header = TRUE, data.table = FALSE)
}

#' Run one direct IgBLAST annotation pass
#'
#' `runIgblast()` runs `igblastn` directly on a plain FASTA file and requests
#' AIRR rearrangement output (`-outfmt 19`). It is used for the first,
#' read-level V(D)J assignment after RAD demultiplexing and can also label an
#' explicit consensus-level pass. FASTQ and gzip conversion deliberately live
#' outside this function; use
#' `rad_reformat(..., to_fasta = TRUE, header_format = "presto")` and pass its
#' `artifacts$fasta` output to prepare RAD reads for this step.
#'
#' V, D, and J database prefixes can be given individually or resolved from
#' a Change-O-style `igblast_base`. The default database names beneath its
#' `database/` directory are
#' `imgt_<organism>_<receptor>_{v,d,j,c}`, and the default auxiliary file is
#' `optional_file/<organism>_gl.aux`. `IGDATA` is set to `igblast_base` only
#' for the IgBLAST subprocess. The older `igblast_db_dir` convention using
#' `IG_V`/`IG_D`/`IG_J` or `TR_V`/`TR_D`/`TR_J` remains supported. Explicit
#' database and auxiliary-file arguments override either convention. Database
#' arguments are BLAST prefixes, not individual `.nhr`, `.nin`, or `.nsq`
#' files.
#'
#' The FASTA is streamed into complete-record chunks and, when useful, those
#' chunks are handled by several concurrent IgBLAST processes. `num_threads`
#' is a total CPU budget shared across at most `max_processes`; it is never
#' multiplied by the process count. Automatic chunks are balanced by sequence
#' bases, not record count. Set `max_processes = 1` for the original monolithic
#' execution. `chunk_bases` can request smaller chunks for benchmarking or
#' unusually heterogeneous inputs; `NULL` creates one chunk per active
#' process. With automatic chunking, inputs below 100,000 total bases remain
#' single-process. Supplying `chunk_bases` opts into chunking below that
#' threshold. Chunk outputs are validated and merged in original input order.
#'
#' IgBLAST AIRR output can contain a blank-call row for a query with no
#' reportable rearrangement. Therefore a read-stage result reports
#' `query_count`, `airr_rows`, `annotated_queries`, and `unannotated_queries`
#' separately. Only rows with a non-empty V, D, or J call count as annotated.
#' More AIRR rows than FASTA queries is treated as invalid output.
#'
#' @param input_fasta Existing uncompressed FASTA query file. FASTQ and
#'   gzip-compressed input are rejected.
#' @param output_tsv Destination for IgBLAST AIRR TSV output, ending in
#'   `.tsv`, `.tab`, or `.txt`.
#' @param receptor Receptor group: `"ig"` or `"tr"`.
#' @param organism Organism name passed to `igblastn -organism`.
#' @param igblast_base Optional IgBLAST installation/IGDATA root. This is the
#'   same directory passed to `AssignGenes.py igblast -b`. It supplies the
#'   conventional databases, auxiliary file, and subprocess `IGDATA` value.
#' @param igblast_db_dir Optional legacy directory containing `IG_V`/`IG_D`/
#'   `IG_J` or `TR_V`/`TR_D`/`TR_J` prefixes.
#' @param igblast_exec Executable name on `PATH` or full path to `igblastn`.
#'   When omitted with `igblast_base`, `<igblast_base>/bin/igblastn` is used
#'   when present.
#' @param aux_file Optional existing IgBLAST auxiliary-data file. This
#'   overrides the file derived from `igblast_base`.
#' @param germline_db_v,germline_db_d,germline_db_j Optional explicit BLAST
#'   database prefixes. These override paths derived from `igblast_db_dir`.
#' @param c_region_db Optional explicit IG constant-region database prefix.
#' @param num_threads Positive whole-number total CPU budget. Concurrent
#'   IgBLAST processes share this budget.
#' @param max_processes Maximum number of concurrent IgBLAST processes. The
#'   default is four; fewer are used for small inputs or thread budgets.
#' @param chunk_bases Optional positive whole-number target bases per FASTA
#'   chunk. `NULL` creates one approximately base-balanced chunk per active
#'   process. Records are never split.
#' @param extra_args Optional additional IgBLAST arguments. Options managed by
#'   this function, including `-outfmt`, cannot be overridden.
#' @param stage Provenance label for this pass: `"read"` (the default) or
#'   `"consensus"`.
#'
#' @return A `sirap_igblast_result` list containing status and stage,
#'   `query_fasta`, `airr`, and `log` artifacts, the resolved configuration,
#'   command provenance, streaming record counts, and timing.
#' @export
runIgblast <- function(input_fasta,
                       output_tsv,
                       receptor = c("ig", "tr"),
                       organism = "human",
                       igblast_base = NULL,
                       igblast_db_dir = NULL,
                       igblast_exec = Sys.getenv("IGBLASTN", unset = "igblastn"),
                       aux_file = NULL,
                       germline_db_v = NULL,
                       germline_db_d = NULL,
                       germline_db_j = NULL,
                       c_region_db = NULL,
                       num_threads = 8L,
                       max_processes = 4L,
                       chunk_bases = NULL,
                       extra_args = NULL,
                       stage = c("read", "consensus")) {
  use_default_executable <- missing(igblast_exec) || is.null(igblast_exec)
  started <- Sys.time()
  receptor <- match.arg(receptor)
  stage <- match.arg(stage)
  organism <- .sirap_igblast_scalar_character(organism, "organism")
  igblast_base <- .sirap_igblast_base(igblast_base)
  num_threads <- .sirap_igblast_positive_integer(num_threads, "num_threads")
  max_processes <- .sirap_igblast_positive_integer(
    max_processes, "max_processes"
  )
  chunk_bases <- .sirap_igblast_chunk_bases(chunk_bases)
  extra_args <- .sirap_igblast_extra_args(extra_args)

  if (use_default_executable) {
    installation_executable <- if (is.null(igblast_base)) {
      NULL
    } else {
      file.path(igblast_base, "bin", "igblastn")
    }
    igblast_exec <- if (!is.null(installation_executable) &&
                        file.exists(installation_executable)) {
      installation_executable
    } else {
      Sys.getenv("IGBLASTN", unset = "igblastn")
    }
  }

  fasta <- .sirap_validate_igblast_fasta(input_fasta)
  output_tsv <- .sirap_igblast_output_path(
    output_tsv, "output_tsv", c("tsv", "tab", "txt")
  )
  log <- .sirap_igblast_output_path(
    paste0(tools::file_path_sans_ext(output_tsv), ".log"), "log"
  )
  if (.sirap_paths_equivalent(output_tsv, fasta$path) ||
      .sirap_paths_equivalent(log, fasta$path)) {
    .sirap_igblast_abort(
      "IgBLAST output and log paths must be different from `input_fasta`."
    )
  }
  if (.sirap_paths_equivalent(output_tsv, log)) {
    .sirap_igblast_abort(
      "IgBLAST output and log paths must identify different files."
    )
  }

  staged_airr <- tempfile(
    pattern = ".sirap-igblast-airr-",
    tmpdir = dirname(output_tsv),
    fileext = ".tsv"
  )
  staged_log <- tempfile(
    pattern = ".sirap-igblast-log-",
    tmpdir = dirname(log),
    fileext = ".log"
  )
  on.exit(unlink(c(staged_airr, staged_log)), add = TRUE)
  chunk_directory <- tempfile(pattern = "sirap-igblast-chunks-")
  on.exit(
    if (dir.exists(chunk_directory)) {
      unlink(chunk_directory, recursive = TRUE, force = TRUE)
    },
    add = TRUE
  )

  executable <- unname(.sirap_check_tools(c(igblastn = igblast_exec))[[1L]])
  dbs <- .sirap_igblast_db_paths(
    igblast_db_dir = igblast_db_dir,
    receptor = receptor,
    germline_db_v = germline_db_v,
    germline_db_d = germline_db_d,
    germline_db_j = germline_db_j,
    c_region_db = c_region_db,
    validate = TRUE,
    igblast_base = igblast_base,
    organism = organism
  )
  aux_file <- .sirap_igblast_aux_file(aux_file, igblast_base, organism)
  command_env <- if (is.null(igblast_base)) {
    NULL
  } else {
    c(IGDATA = igblast_base)
  }

  plan <- .sirap_igblast_process_plan(
    fasta = fasta,
    num_threads = num_threads,
    max_processes = max_processes,
    chunk_bases = chunk_bases
  )
  chunks <- .sirap_split_igblast_fasta(fasta, chunk_directory, plan)
  processes <- as.integer(min(
    max_processes, num_threads, length(chunks)
  ))
  threads_per_process <- if (processes == 1L) {
    num_threads
  } else {
    as.integer(floor(num_threads / processes))
  }

  args <- .sirap_igblast_args(
    input_fasta = fasta$path,
    output_tsv = output_tsv,
    receptor = receptor,
    organism = organism,
    dbs = dbs,
    aux_file = aux_file,
    num_threads = num_threads,
    extra_args = extra_args
  )
  jobs <- lapply(seq_along(chunks), function(index) {
    chunk <- chunks[[index]]
    chunk_output <- if (length(chunks) == 1L) {
      staged_airr
    } else {
      file.path(chunk_directory, sprintf("airr-%05d.tsv", index))
    }
    chunk_log <- if (length(chunks) == 1L) {
      staged_log
    } else {
      file.path(chunk_directory, sprintf("igblast-%05d.log", index))
    }
    list(
      index = as.integer(index),
      input = chunk$path,
      output = chunk_output,
      log = chunk_log,
      queries = chunk$queries,
      bases = chunk$bases,
      threads = threads_per_process,
      args = .sirap_igblast_args(
        input_fasta = chunk$path,
        output_tsv = chunk_output,
        receptor = receptor,
        organism = organism,
        dbs = dbs,
        aux_file = aux_file,
        num_threads = threads_per_process,
        extra_args = extra_args
      )
    )
  })
  results <- .sirap_run_igblast_jobs(
    jobs = jobs,
    executable = executable,
    command_env = command_env,
    processes = processes
  )

  if (length(results) > 1L) {
    .sirap_merge_igblast_airr(results, staged_airr)
    .sirap_merge_igblast_logs(results, staged_log)
  }
  airr <- .sirap_inspect_airr_output(staged_airr)
  chunk_rows <- sum(vapply(results, function(result) result$airr$rows, numeric(1)))
  chunk_annotations <- sum(vapply(
    results, function(result) result$airr$annotated_rows, numeric(1)
  ))
  if (airr$rows != chunk_rows || airr$annotated_rows != chunk_annotations) {
    .sirap_igblast_abort(
      "Merged IgBLAST AIRR counts do not match the validated chunk outputs."
    )
  }
  if (airr$rows > fasta$queries) {
    .sirap_igblast_abort(paste0(
      "IgBLAST returned ", airr$rows, " AIRR rows for ", fasta$queries,
      " FASTA queries. A direct AIRR pass cannot contain more data rows ",
      "than input queries."
    ))
  }
  annotated_queries <- airr$annotated_rows
  unannotated_queries <- fasta$queries - annotated_queries
  log <- .sirap_publish_igblast_artifact(staged_log, log, "log")
  output_tsv <- .sirap_publish_igblast_artifact(
    staged_airr, output_tsv, "AIRR"
  )
  completed <- Sys.time()

  structure(
    list(
      status = "success",
      stage = stage,
      artifacts = list(
        query_fasta = fasta$path,
        airr = output_tsv,
        log = log
      ),
      config = list(
        stage = stage,
        receptor = receptor,
        organism = organism,
        igblast_base = igblast_base,
        igblast_db_dir = if (is.null(igblast_db_dir)) NULL else {
          normalizePath(path.expand(igblast_db_dir), mustWork = TRUE)
        },
        databases = dbs,
        auxiliary_data = aux_file,
        executable = executable,
        threads = num_threads,
        thread_budget = num_threads,
        max_processes = max_processes,
        processes = processes,
        threads_per_process = threads_per_process,
        minimum_parallel_bases = plan$minimum_parallel_bases,
        requested_chunk_bases = chunk_bases,
        extra_args = extra_args
      ),
      command = list(
        executable = executable,
        args = args,
        execution_args = if (length(jobs) == 1L) jobs[[1L]]$args else NULL,
        executions = lapply(jobs, function(job) {
          list(
            index = job$index,
            queries = job$queries,
            bases = job$bases,
            threads = job$threads,
            args = job$args
          )
        }),
        env = command_env
      ),
      stats = list(
        query_count = fasta$queries,
        queries = fasta$queries,
        query_bases = fasta$bases,
        airr_rows = airr$rows,
        annotated_queries = annotated_queries,
        unannotated_queries = unannotated_queries,
        annotation_rate = annotated_queries / fasta$queries,
        chunks = length(chunks),
        chunk_queries = vapply(chunks, `[[`, numeric(1), "queries"),
        chunk_bases = vapply(chunks, `[[`, numeric(1), "bases"),
        airr_bytes = airr$bytes,
        log_bytes = unname(file.info(log)$size)
      ),
      timing = list(
        started = started,
        completed = completed,
        elapsed_seconds = as.numeric(difftime(completed, started, units = "secs"))
      )
    ),
    class = c("sirap_igblast_result", "list")
  )
}
