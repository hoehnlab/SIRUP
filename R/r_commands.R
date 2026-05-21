.datatable.aware = TRUE

# Shared FASTQ/FASTA header contract for read_fastqas / write_fastqas /
# stream_fastqas:
#   - `id`      always present; short id (everything before first space),
#               with no leading `@` (fastq) or `>` (fasta).
#   - `full_id` present only when `full_id = TRUE`; full header (name +
#               comment) with no leading `@`/`>`.
#   - `seq`     always present.
#   - `qual`    present for fastq.
#   - Column order: id, full_id (if requested), seq, qual (if fastq).
#
# `type` accepts "fq" / "fastq" / "fa" / "fasta" in all three functions.

.fastqas_normalize_type <- function(type) {
  type <- match.arg(type, c("fq", "fastq", "fa", "fasta"))
  if (type %in% c("fa", "fasta")) "fa" else "fq"
}

#' Write a data.table of sequences to FASTQ/FASTA
#'
#' @param df data.table or data.frame with at least `seq` (and `qual` for
#'   fastq). The header column is chosen by `header_col`.
#' @param fn output path.
#' @param type one of `"fq"`, `"fastq"`, `"fa"`, `"fasta"`.
#' @param header_col column name to use as the record header. Defaults to
#'   `"id"`. Pass `"full_id"` to write the full header instead.
#' @param append append to an existing file?
#' @param nthreads `data.table::fwrite` thread count.
#' @return Invisibly returns `fn`.
#' @export
write_fastqas <- function(df, fn, type = "fq", header_col = "id",
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

#' Read FASTQ/FASTA into a data.table
#'
#' Returns a data.table matching the shared contract documented above:
#' columns `id`, optionally `full_id`, then `seq`, then `qual` (fastq only).
#' Headers are stripped of their leading `@` / `>`. `id` is the first
#' whitespace-delimited token of the header; `full_id` is the entire
#' header.
#'
#' @param fn input path, or a character vector / list of paths.
#' @param type one of `"fq"`, `"fastq"`, `"fa"`, `"fasta"`.
#' @param full_id include a `full_id` column with the full header.
#' @param ... passed through to `data.table::fread`.
#' @return data.table.
#' @export
read_fastqas <- function(fn, type = "fq", full_id = FALSE, ...) {
  type <- .fastqas_normalize_type(type)
  include_full <- isTRUE(full_id)

  rfq_single <- function(path) {
    ext <- tools::file_ext(path)
    sysname <- Sys.info()[["sysname"]]
    cat_cmd <- if (ext == "gz") {
      switch(sysname,
             Linux  = "zcat",
             Darwin = "gunzip -c",
             stop("Unsupported OS: ", sysname))
    } else {
      "cat"
    }

    if (type == "fq") {
      res <- data.table::fread(
        cmd = glue::glue("{cat_cmd} {path} | paste - - - - | cut -f1,2,4"),
        col.names = c("hdr", "seq", "qual"),
        sep = "\t", header = FALSE, quote = "",
        data.table = TRUE, ...
      )
      res[, hdr := sub("^@", "", hdr)]
      res[, id := data.table::tstrsplit(hdr, " ", fixed = TRUE)[[1]]]
      if (include_full) {
        data.table::setnames(res, "hdr", "full_id")
        data.table::setcolorder(res, c("id", "full_id", "seq", "qual"))
      } else {
        res[, hdr := NULL]
        data.table::setcolorder(res, c("id", "seq", "qual"))
      }
      return(res[])
    }

    awk_cmd <- 'awk \'/^>/{if (NR>1) printf("\\n"); printf("%s\\t",$0); next} {printf("%s",$0);} END {printf("\\n");}\''
    res <- data.table::fread(
      cmd = glue::glue("{cat_cmd} {path} | {awk_cmd}"),
      col.names = c("hdr", "seq"),
      sep = "\t", data.table = TRUE, ...
    )
    res[, hdr := sub("^>", "", hdr)]
    res[, id := data.table::tstrsplit(hdr, " ", fixed = TRUE)[[1]]]
    if (include_full) {
      data.table::setnames(res, "hdr", "full_id")
      data.table::setcolorder(res, c("id", "full_id", "seq"))
    } else {
      res[, hdr := NULL]
      data.table::setcolorder(res, c("id", "seq"))
    }
    res[]
  }

  if (is.list(fn) || length(fn) > 1) {
    out <- pbapply::pblapply(fn, rfq_single)
    data.table::rbindlist(out)
  } else {
    rfq_single(fn)
  }
}

#' @export
read_paf<-function(path, ...){
  df<-data.table::fread(input = path, 
    fill = Inf, sep = "\t", strip.white = FALSE, header = FALSE, na.strings = "")
  df[df == '']<-NA
  ids<-c("id","ql","qs","qe","qr","sn","sl","ss","se","mm","mt","mq")
  data.table::setnames(df, c(ids,
    df[,ncol(df):ncol(df),] %>% {
      which(!is.na(.))} %>% 
      df[.,13:ncol(df),] %>% 
      {substr(unlist(.), start = 1, stop = 2)} %>% 
      unique))
  for (col in (length(ids) + 1):ncol(df)) {
    tag_prefixes<-substr(df[[col]], start = 1, stop = 2)
    expected_tag<-colnames(df)[col]
    mismatch_indices<-which(tag_prefixes != expected_tag & !is.na(tag_prefixes))
    if (length(mismatch_indices) > 0) {
      # Shifting columns to the right for mismatched rows
      df[mismatch_indices, (col + 1):ncol(df) := .SD, .SDcols = col:(ncol(df) - 1)]
      df[mismatch_indices, (col) := NA]  # Set current column to NA for mismatched rows
    }
    df[, (expected_tag) := gsub(pattern = paste0(expected_tag,":.:"), replacement = "", x = get(expected_tag))]
  }
  for (col in names(df)) {
    numeric_conversion <- suppressWarnings(as.numeric(df[[col]]))
    if (!any(is.na(numeric_conversion)) || all(is.na(numeric_conversion) == is.na(df[[col]]))) {
      df[[col]] <- numeric_conversion
    }
  }
  return(df)
}

#' @export
aggc <- function() {
  gc()
  sysname <- Sys.info()[["sysname"]]
  if (sysname == "Linux") {
    mallinfo::malloc.trim()
  } 
  if (sysname == "Darwin") {  
    # macOS
    system("purge", ignore.stdout = TRUE, ignore.stderr = TRUE)
  }
}

#' @export
substring_odds<-function(subsequence_length, sequence_length){
  return((mapply(sequence_length, subsequence_length, FUN = 
                   function(seq_len, sub_seq){
                     1-(1-(0.25)^sub_seq)^(seq_len-sub_seq+1)
  }, USE.NAMES = FALSE, SIMPLIFY = TRUE)))
}

# ---- R helper to build a sparse matrix ----
#' @export
sparse_onehot_matrix <- function(cigars, aln_starts, aln_ends,
                                 region_start, region_end,
                                 rownames = NULL,
                                 logical_matrix = FALSE,
                                 bin = 1L,
                                 agg = c("any", "sum")) {
  stopifnot(requireNamespace("Matrix", quietly = TRUE))
  agg <- match.arg(agg)
  if (!is.integer(bin)) bin <- as.integer(bin)
  
  trip <- cigar_onehot_triplets(cigars, aln_starts, aln_ends, region_start, region_end)
  i  <- trip$i
  j  <- trip$j
  nr <- trip$nrow
  nc <- trip$ncol
  
  if (length(i) == 0L) {
    # Empty matrix fast-path
    nc_eff <- if (bin > 1L) ceiling(nc / bin) else nc
    M <- Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(),
                              dims = c(nr, nc_eff), repr = "C")
    if (!is.null(rownames)) rownames(M) <- rownames
    if (logical_matrix) M <- (M != 0)
    attr(M, "bin") <- bin
    attr(M, "region_start") <- region_start
    attr(M, "region_end") <- region_end
    attr(M, "bin_edges0") <- seq.int(from = region_start,
                                     to   = region_start + nc_eff*bin,
                                     by   = bin)[1:nc_eff]
    return(M)
  }
  
  if (bin > 1L) {
    # Map original 1-based columns j to 1-based bin columns jb
    jb     <- ((as.integer(j) - 1L) %/% bin) + 1L
    nc_bin <- as.integer((nc + bin - 1L) %/% bin)
    
    if (agg == "any") {
      # Any coverage in a bin -> 1
      x <- rep.int(1L, length(jb))
      M <- Matrix::sparseMatrix(i = i, j = jb, x = x,
                                dims = c(nr, nc_bin), repr = "C")
      if (logical_matrix) {
        M <- (M != 0)  # lgCMatrix
      } else {
        # Binarize in case duplicates summed > 1
        M@x[] <- 1L
      }
    } else { # agg == "sum" = number of covered bases per bin
      # Each triplet represents one base; summing duplicates yields base count per bin.
      x <- rep.int(1L, length(jb))
      M <- Matrix::sparseMatrix(i = i, j = jb, x = x,
                                dims = c(nr, nc_bin), repr = "C")
      if (logical_matrix) {
        # For logical output, turn any nonzero sum to TRUE
        M <- (M != 0)
      }
    }
    
    if (!is.null(rownames)) rownames(M) <- rownames
    attr(M, "bin") <- bin
    attr(M, "region_start") <- region_start
    attr(M, "region_end") <- region_end
    # 0-based half-open starts for each bin on the reference
    attr(M, "bin_edges0") <- seq.int(from = region_start,
                                     to   = region_start + nc_bin*bin,
                                     by   = bin)[1:nc_bin]
    return(M)
  }
  
  # No binning (bin == 1): original behavior
  x <- rep.int(1L, length(i))
  M <- Matrix::sparseMatrix(i = i, j = j, x = x,
                            dims = c(nr, nc), repr = "C")
  if (!is.null(rownames)) rownames(M) <- rownames
  if (logical_matrix) {
    M <- (M != 0)
  } else {
    M@x[] <- 1L
  }
  attr(M, "bin") <- 1L
  attr(M, "region_start") <- region_start
  attr(M, "region_end") <- region_end
  attr(M, "bin_edges0") <- seq.int(from = region_start,
                                   to   = region_start + nc,
                                   by   = 1L)[1:nc]
  M
}
