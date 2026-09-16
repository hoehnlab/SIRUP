.datatable.aware = TRUE

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
substringOdds<-function(subsequence_length, sequence_length){
  return((mapply(sequence_length, subsequence_length, FUN = 
                   function(seq_len, sub_seq){
                     1-(1-(0.25)^sub_seq)^(seq_len-sub_seq+1)
  }, USE.NAMES = FALSE, SIMPLIFY = TRUE)))
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
