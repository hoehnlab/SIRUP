.datatable.aware = TRUE

#' @export
write_fastqas<-function(df, fn, type = "fq", append = FALSE, nthreads = 1){
  if(type == "fq"){
    fq_list<-df %>%
      dplyr::mutate(dummy = "+") %>%
      dplyr::select(id, seq, dummy, qual) %>%
      as.matrix() %>%
      t() %>%
      as.character() %>%
      tibble::as_tibble()
    data.table::fwrite(fq_list, file = fn, compress = "auto", col.names = FALSE, quote = FALSE, append = append, nThread = nthreads)
  }
  if(type == "fa"){
    if(any(grepl("^>", df$id)) == FALSE){
      df$id<-paste0(">",df$id)
      fa_list<-df %>%
        dplyr::select(id, seq) %>%
        as.matrix() %>%
        t() %>%
        as.character() %>%
        tibble::as_tibble()
      data.table::fwrite(fa_list, file = fn, compress = "auto", col.names = FALSE, quote = FALSE, append = append, nThread = nthreads)
    } else {
    fa_list<-df %>%
      dplyr::select(id, seq) %>%
      as.matrix() %>%
      t() %>%
      as.character() %>%
      tibble::as_tibble()
    data.table::fwrite(fa_list, file = fn, compress = "auto", col.names = FALSE, quote = FALSE, append = append, nThread = nthreads)
    }
  }
}

#' @export
read_fastqas<-function(fn, type = "fq", full_id = FALSE, ...){
  rfq_single<-function(fn, type, full_id = FALSE, ...){
    ext<-tools::file_ext(fn)
    sysname<-Sys.info()[["sysname"]]
    if (ext == "gz"){
      if(sysname == "Linux"){
        cat_cmd<-"zcat"
      }else if(sysname =="Darwin"){
        cat_cmd<-"gunzip -c"
      }else{
        stop("Unsupported OS")
      }
    } else {
      cat_cmd<-"cat"
    }
    if (type == "fq"){
      res<-data.table::fread(
        cmd = glue::glue("{cat_cmd} {fn} | paste - - - - | cut -f1,2,4"),
        col.names = c("id", "seq", "qual"),
        sep = "\t", header = FALSE, quote = "", data.table = TRUE, ...
      )
      res<-res[, .(id = data.table::tstrsplit(id, " ", fixed = TRUE)[[1]], qual, seq)]
      if (full_id) {
        res[, full_id := id]
      }
      return(res)
    }
    if(type == "fa"){
      awk_cmd<-'awk \'/^>/{if (NR>1) printf("\\n"); printf("%s\\t",$0); next} {printf("%s",$0);} END {printf("\\n");}\''
      res<-data.table::fread(
        cmd = glue::glue("{cat_cmd} {fn} | {awk_cmd}"),
        col.names = c("id", "seq"),
        sep = "\t", data.table = TRUE, ...
      )
      return(res)
    }
  }
  if (is.list(fn) || length(fn) > 1) {
    out<-pbapply::pblapply(fn, FUN = function(X) {
      rfq_single(fn = X, type = type, full_id = full_id, ...)
    }, ...)
    dt<-data.table::rbindlist(out)
    return(dt)
  } else {
    return(rfq_single(fn = fn, type = type, full_id = full_id, ...))
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
