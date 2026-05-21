.datatable.aware = TRUE

#' @export
whitelist_importer<-function(whitelist_path, save = NULL, convert = FALSE){
  cat("Importing a little bit of your whitelist to figure out its type...\n")
  whitelist<-data.table::fread(input = whitelist_path, header = FALSE, data.table = TRUE, nrows = 10, nThread = 1)
  if(ncol(whitelist) == 1){
    data.table::setnames(whitelist, new = "whitelist_bcs")
  } else {
    cat("More than one column detected! Assuming the first one is the one with the barcodes.\n")
    data.table::setnames(whitelist, old = "V1", new = "whitelist_bcs")
  }
  if(class(whitelist$whitelist_bcs) == "character" || 
      !grepl(pattern = "bitlist", x = whitelist_path, fixed = TRUE)){
    cat("Character type detected! Importing the whitelist.\n")
    if(ncol(whitelist) == 1){
      whitelist<-data.table::fread(
        input = whitelist_path,
        header = FALSE,
        data.table = TRUE,
        nThread = 1,
        col.names = "whitelist_bcs"
        )
    } else {
      whitelist<-data.table::fread(input = whitelist_path, header = FALSE, data.table = TRUE, nThread = 1)
      whitelist<-whitelist[,1]
      data.table::setnames(whitelist, "whitelist_bcs")
    }
    if(convert != FALSE){
      cat("Converting the whitelist into a more efficient format!\n")
      whitelist$whitelist_bcs<-sequence_to_bits(whitelist$whitelist_bcs)
    }
    whitelist<-data.table::setkey(whitelist, "whitelist_bcs")
    gc()
    if(!is.null(save)){
      if(save == TRUE){
        new_path<-paste0(dirname(whitelist_path),"/",tools::file_path_sans_ext(whitelist_path)
          %>% basename(.),
          "_bitlist.csv.gz")
        cat(paste0("Now saving the integer whitelist to ",new_path,"!\n"))
        data.table::fwrite(x = whitelist, file = new_path, compress = "auto", nThread = 1, showProgress = TRUE, col.names = FALSE)
        gc()
      }
    }
    return(whitelist)
  }
  else{
    if(class(whitelist$whitelist_bcs) == "integer"){
      cat("Integer64-type whitelist detected! Now importing...\n")
      whitelist<-data.table::fread(input = whitelist_path, header = FALSE, data.table = TRUE, nThread = 1,
        col.names = "whitelist_bcs", key = "whitelist_bcs", integer64 = "integer64")
      gc()
      cat(paste0("Your whitelist is ", nrow(whitelist), " barcodes long!\n"))
      print(head(whitelist))
      return(whitelist)
    }
  }
}

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
prep_read_layout<-function(read_layout_form, return_env = .GlobalEnv) {
  cat("Importing read layout and figuring out its order!\n")
  tidytable::setDTthreads(threads = 1)
  read_layout<-data.table::fread(file = read_layout_form, header = TRUE, fill = TRUE, na.strings = "", skip = 1)
  read_layout$expected_length[which(!is.na(read_layout$seq))]<-stringr::str_length(read_layout$seq[which(!is.na(read_layout$seq))])
  
  read_layout[(class %in% c("poly_a", "poly_t")) & is.na(expected_length), expected_length := 12]
  read_layout[class == "poly_a", seq := paste0("A{", expected_length, ",}+")]
  read_layout[class == "poly_t", seq := paste0("T{", expected_length, ",}+")]
  
  new_row_top<-data.table::data.table(id = "seq_start", seq = NA, expected_length = 0,
    type = "static", class = "start", direction = "forward", class_id = "seq_start")
  new_row_bottom<-data.table::data.table(id = "seq_stop", seq = NA, expected_length = 0, 
    type = "static", class = "stop", direction = "forward", class_id = "seq_stop")
  
  read_layout<-data.table::rbindlist(list(new_row_top, read_layout, new_row_bottom), fill = TRUE)
  
  if(any(read_layout$class != "forw_primer")){
    start_index<-which(read_layout$class == "start")
    if(read_layout$type[start_index+1] == "static"|is.na(read_layout$class[start_index+1])){
      read_layout[start_index+1, class := "forw_primer"]
    }
  }
  
  if(any(read_layout$class != "rev_primer")){
    stop_index<-which(read_layout$class == "stop")
    if(read_layout$type[stop_index-1] == "static"|is.na(read_layout$class[stop_index-1])){
      read_layout[stop_index-1, class := "rev_primer"]
    }
  }
  
  if (any(duplicated(read_layout$id))) {
    read_layout$id[which(duplicated(read_layout$id) | duplicated(read_layout$id, fromLast = TRUE))]<-
      which(duplicated(read_layout$id) | duplicated(read_layout$id, fromLast = TRUE)) %>%
      read_layout$id[.] %>% paste0(., "_", seq(length(.)))
  }
  if (any(duplicated(read_layout$class))) {
    read_layout$class[which(duplicated(read_layout$class) | duplicated(read_layout$class, fromLast = TRUE))] %>% unique %>%
      sapply(., FUN = function(x){
        out<-which(read_layout$class == x) %>%
          read_layout[., class_id := paste0(class, "_", seq(.))]
        return(out)
      }, simplify = FALSE, USE.NAMES = FALSE) %>% invisible(.)
    read_layout$class_id[which(is.na(read_layout$class_id))]<-read_layout$class[which(is.na(read_layout$class_id))]
  } else {
    read_layout$class_id <- ifelse(test = (read_layout$class == "start" | read_layout$class == "stop"), 
      yes = read_layout$class_id, 
      no = read_layout$class)
  }
  
  read_layout<-read_layout %>%
    .[, direction := "forward"] %>%
    data.table::copy(.) %>%
    .[!class %in% c("start", "stop")] %>%
    .[, seq := ifelse(class == "poly_t", "A{12,}+", ifelse(class == "poly_a", "T{12,}+", sapply(seq, revcomp)))] %>%
    .[, id := ifelse(class %in% c("poly_a", "poly_t"), ifelse(class == "poly_a", "poly_t", "poly_a"), paste0("rc_", id))] %>%
    .[, class_id := ifelse(class_id == "poly_a", "poly_t", ifelse(class_id == "poly_t", "poly_a", class_id))] %>%
    .[, class := ifelse(class == "poly_a", "poly_t", ifelse(class == "poly_t", "poly_a", class))] %>%
    .[, direction := "reverse"] %>%
    .[order(seq(nrow(.)), decreasing = TRUE), ] %>%
    {data.table::rbindlist(list(data.table::data.table(id = "seq_start", seq = NA, expected_length = 0, type = "static",
      class = "start", direction = "reverse", class_id = "seq_start"),
      ., data.table::data.table(id = "seq_stop", seq = NA, expected_length = 0, type = "static",
        class = "stop", direction = "reverse", class_id = "seq_stop")),
      fill = TRUE)} %>%
    {rbind(read_layout, .)} %>%
    {.[, class_id := ifelse(test = {.$direction == "reverse" &!
        (.$class_id == "poly_a" | .$class_id == "poly_t")},# | .$class_id == "seq_start" | .$class_id == "seq_stop")},
      yes = paste0("rc_", .$class_id),
      no = .$class_id)]} %>%
    {.[, "order" := seq(1:nrow(.))]} %>%
    {.[, class := ifelse(test = {.$class == "poly_a" | .$class == "poly_t"},
      yes = "poly_tail",
      no = .$class)]}
  
  adapters <- read_layout$seq[-which(is.na(read_layout$seq))]
  names(adapters) <- read_layout$class_id[-which(is.na(read_layout$seq))]
  return(list2env(x = list(read_layout = read_layout, adapters = adapters), envir = return_env))
}

#' @export
stat_collector<-function(df, read_layout, mode = "stats"){
  forward_adapters<-read_layout[type %in% c("static") & direction == "forward" & class != "poly_tail" & class != "start" & class != "stop", class_id]
  reverse_adapters<-read_layout[type %in% c("static") & direction == "reverse" & class != "poly_tail" & class != "start" & class != "stop", class_id]
  df<-tidytable::as_tidytable(df)
  #changed best_edit_distance from == 0 to <= 1. 
  #not the best fix, as it leaves room for perfect concatenates with no errors or something weird like that
  #but it also adds for some more sequences that might otherwise be discarded.
  forward_zero_ids<-df %>%
    tidytable::filter(query_id %in% forward_adapters) %>%
    tidytable::group_by(id) %>%
    tidytable::filter(all(best_edit_distance <= 1)) %>%
    tidytable::pull(id)
  reverse_zero_ids<-df %>%
    tidytable::filter(query_id %in% reverse_adapters) %>%
    tidytable::group_by(id) %>%
    tidytable::filter(all(best_edit_distance <= 1)) %>%
    tidytable::pull(id)
  if(mode == "stats"){
    stats_reverse<-df %>%
      tidytable::filter(id %in% forward_zero_ids, query_id %in% reverse_adapters) %>%
      tidytable::group_by(query_id) %>%
      tidytable::summarise(
        misal_threshold = mean(best_edit_distance),
        misal_sd = sd(best_edit_distance))
    stats_forward<-df %>%
      tidytable::filter(id %in% reverse_zero_ids, query_id %in% forward_adapters) %>%
      tidytable::group_by(query_id) %>%
      tidytable::summarise(
        misal_threshold = mean(best_edit_distance),
        misal_sd = sd(best_edit_distance))
    misalignment_thresholds<-rbind(stats_forward, stats_reverse)
    return(misalignment_thresholds)
  }
  if(mode == "graphics"){
    stats_reverse<-df %>%
      tidytable::filter(id %in% forward_zero_ids, query_id %in% reverse_adapters) %>%
      tidytable::group_by(query_id)
    stats_forward<-df %>%
      tidytable::filter(id %in% reverse_zero_ids, query_id %in% forward_adapters) %>%
      tidytable::group_by(query_id)
    total_stats<-data.table::rbindlist(list(stats_forward, stats_reverse))
    return(total_stats)
  }
}

#' @export
prep_external_path<-function(external_path_form, create_output_dir = TRUE, ...){
  print("Checking executable paths!")
  path_layout<-data.table::fread(file = external_path_form, header = TRUE, fill = TRUE, na.strings = "", skip = 1,
    key = "path_type", data.table = TRUE)
  if(dir.exists(path_layout[path_type == "output_dir", actual_path]) == FALSE & create_output_dir == TRUE){
    dir.create(path = path_layout[path_type == "output_dir", actual_path])
  }
  whitelist<-whitelist_importer(whitelist_path = path_layout[path_type == "whitelist", actual_path], ...)
  return(list2env(x = list(path_layout = path_layout, whitelist = whitelist), envir = .GlobalEnv))
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
process_barcodes<-function(barcode_path, read_layout) {
  barcodes<-data.table::fread(barcode_path, col.names = c("seq", "count","concatenate_count"))
  data.table::setkey(read_layout, "class_id")
  
  dir_path<-paste0(dirname(barcode_path),"/")
  barcode_id<-basename(tools::file_path_sans_ext(barcode_path, compression = TRUE))
  
  expected_length<-read_layout[barcode_id, expected_length]
  barcodes[, filtered := NA_character_]
  barcodes[(stringr::str_length(seq) != expected_length), filtered := "length_filter"]
  rescue_filter<-ceiling(expected_length * 0.5)
  patterns<-sprintf("%s{%d,}", c("A", "C", "T", "G"), rescue_filter + 1)
  poly_runs<-data.table::rbindlist(lapply(patterns, function(pattern) {
    barcodes[is.na(filtered) & grepl(pattern, seq), .(seq, count, filtered = "poly_run_filtered")]
  }))
  if (nrow(poly_runs) > 0){
    barcodes[poly_runs, on = .(seq), filtered := i.filtered]
  }
  static_seqs<-read_layout[type == "static" & !is.na(expected_length) & expected_length > 0 & class != "poly_tail", seq]
  static_qgrams<-colnames(stringdist::qgrams(static_seqs, q = expected_length, useNames = TRUE))
  static_qgrams_mutations<-static_qgrams %>% 
    sequence_to_bits(.) %>% 
    generate_recursive_quaternary_mutations(., mutation_rounds =  2, sequence_length = 16) %>%
    bits_to_sequence(.)
  static_qgrams<-append(static_qgrams, static_qgrams_mutations)
  
  barcodes[is.na(filtered) & seq %in% static_qgrams, filtered := "adapter_segment_filtered"]
  
  barcodes[, `:=`(
    int64_seq = bit64::as.integer64(NA),
    int64_rcseq = bit64::as.integer64(NA)
  )]
  candidate_idx <- which(is.na(barcodes$filtered) & stringr::str_length(barcodes$seq) == expected_length)
  if (length(candidate_idx) > 0) {
    seq_bits <- bit64::as.integer64(sequence_to_bits(barcodes[candidate_idx, seq]))
    rc_bits <- bit64::as.integer64(sequence_to_bits(revcomp(barcodes[candidate_idx, seq])))
    detected <- rad_whitelist_detect_cpp(
      sequences = barcodes[candidate_idx, seq],
      barcode_class_id = barcode_id,
      read_layout = read_layout,
      nthreads = 1L,
      verbose = FALSE
    )

    barcodes[candidate_idx, `:=`(
      int64_seq = seq_bits,
      int64_rcseq = rc_bits,
      filtered = detected
    )]
  }
  data.table::fwrite(x = barcodes, file = barcode_path, col.names = TRUE)
}

#' @export
correct_barcodes<-function(input_file,
                           verbose = FALSE,
                           nthreads = 1,
                           mutation_rounds = 3,
                           bc_correction_mode = c("offensive", "defensive"),
                           chunk_size = 50000,
                           read_layout = read_layout,
                           input_fastq_path = NULL,
                           filtered_fastq_path = NULL,
                           unfiltered_fastq_path = NULL
                           ){
  bc_correction_mode <- match.arg(bc_correction_mode)
  barcode<-data.table::fread(
    input_file, key = "filtered")
  data.table::setkey(read_layout, "class_id")
  compressed<-grepl(
    pattern = ".gz", 
    x = input_file, fixed = TRUE)
  barcode_id<-basename(tools::file_path_sans_ext(input_file, compression = compressed))
  
  dir_path<-paste0(dirname(input_file),"/")
  if (is.null(input_fastq_path)) {
    input_fastq_path<-paste0(dirname(dir_path),"/demuxed_reads.fastq.gz")
  }
  if (is.null(filtered_fastq_path)) {
    filtered_fastq_path<-paste0(dirname(dir_path),"/filtered_reads.fastq.gz")
  }
  if (is.null(unfiltered_fastq_path)) {
    unfiltered_fastq_path<-paste0(dirname(dir_path),"/unfiltered_reads.fastq.gz")
  }
  barcode_mutation_data_path<-paste0(dir_path, barcode_id, "_mutation_counts.csv.gz")
  
  invalid_barcodes <- unique(barcode[filtered == "barcode_to_correct", seq])
  invalid_barcodes <- invalid_barcodes[!is.na(invalid_barcodes) & nzchar(invalid_barcodes)]

  if (length(invalid_barcodes) == 0L) {
    correction_map <- data.table::data.table(
      incorrect_barcodes = character(),
      corrected_barcodes = character(),
      resolved_status = logical(),
      mutation_type = character()
    )
    data.table::fwrite(correction_map, file = barcode_mutation_data_path)
    if (!is.null(input_fastq_path) && !is.null(filtered_fastq_path) && file.exists(input_fastq_path)) {
      file.copy(input_fastq_path, filtered_fastq_path, overwrite = TRUE)
    }
    return(invisible(correction_map))
  }

  correction_map <- rad_barcode_correct_cpp(
    sequences = invalid_barcodes,
    barcode_class_id = barcode_id,
    read_layout = read_layout,
    mode = bc_correction_mode,
    mut_dist = mutation_rounds,
    nthreads = nthreads,
    verbose = verbose
  )
  data.table::fwrite(correction_map, file = barcode_mutation_data_path)
  correction_map<-data.table::fread(input = barcode_mutation_data_path, na.strings = "")
  
  correction_map<-correction_map[which(correction_map$resolved_status == TRUE),]
  if (nrow(correction_map) == 0L) {
    if (!is.null(input_fastq_path) && !is.null(filtered_fastq_path) && file.exists(input_fastq_path)) {
      file.copy(input_fastq_path, filtered_fastq_path, overwrite = TRUE)
    }
    return(invisible(correction_map))
  }
  
  #this guy really runs slow, but i think it's chunk size dependent
  fastq_correction(
    input_file = input_fastq_path,
    output_file = filtered_fastq_path,
    incorrect_bcs = correction_map$incorrect_barcodes,
    correct_bcs = correction_map$corrected_barcodes,
    barcode_id = barcode_id,
    nthreads = nthreads,
    chunk_size = chunk_size,
    downsample = -1
  )
}

#' @export
process_sig<-function(file_path,
                      output_prefix,
                      chunk_size = 50000,
                      max_sequences = -1,
                      compress = TRUE,
                      output_type = "summary",
                      nthreads = 1){
  tabulate_sigs(file_path = file_path, 
                output_prefix = output_prefix, 
                chunk_size = chunk_size,
                max_sequences = max_sequences, 
                compress = compress, 
                nthreads = nthreads,
                output_type = output_type)
  summary_results<-data.table::fread(file = paste0(output_prefix,"/summary_table.csv",  
                                                   ifelse(test = compress, yes = ".gz", no = "")),
                                     header = TRUE)
  valid_reads<-sum(summary_results$unique_id_count[which(summary_results$direction == "F"|summary_results$direction == "R")])
  total_reads<-sum(summary_results$unique_id_count)
  summary_plot<-ggpubr::ggbarplot(data = summary_results, 
                                  x = "direction", 
                                  y = "unique_id_count",
                                  color = "concatenate",
                                  fill = "concatenate", 
                                  title = "SigString Summary", 
                                  legend = "right", 
                                  legend.title = "Concatenate Status", 
                                  xlab = "Read (or Filter) Type",
                                  ylab = "Number of Reads", 
                                  order = c("F", 
                                            "R", 
                                            "missing_read_or_barcode_undecided", 
                                            "invalid_barcode_length_undecided",
                                            "invalid_read_length_undecided"), 
                                  subtitle = paste0(
                                    valid_reads, 
                                    " valid reads out of a total ", 
                                    total_reads, " total reads."), 
                                  x.text.angle = 45, palette = "npj")
  ggplot2::ggsave(filename = paste0(output_prefix,"summary_barplot.pdf"), 
                  plot = summary_plot, 
                  device = "pdf", 
                  width = 8, 
                  height = 11, 
                  units = "in", 
                  dpi = "retina")
  results<-data.table::data.table(unique_id_count = paste0(valid_reads, " valid reads captured from ", total_reads, " total reads."))
  summary_results<-data.table::rbindlist(l = list(summary_results, results), fill = TRUE)
  data.table::fwrite(summary_results, paste0(output_prefix,"/summary_table.csv", ifelse(test = compress, yes = ".gz", no = "")))
}

#' @export
rad_run<-function(
    fastq_file_or_directory_path, 
    read_layout_path, 
    output_directory_path, 
    compress = TRUE,
    generate_sigstring_diagnostics = FALSE,
    sigstring_diagnostic_verbose = FALSE,
    misalignment_threshold_path = NULL,
    tabulated_sigstring_count = -1,
    chunk_size = 50000,
    mutation_rounds = 3,
    bc_correction_mode = c("offensive", "defensive"),
    nthreads = 1){
  bc_correction_mode <- match.arg(bc_correction_mode)
  start_time = Sys.time()
  if(!dir.exists(output_directory_path)){
    dir.create(output_directory_path)
  }
  prep_read_layout(read_layout_form = read_layout_path)
  if(!is.null(misalignment_threshold_path)){
    misalignment_threshold<-data.table::fread(misalignment_threshold_path)
    out <- list(
      misalignment_threshold = misalignment_threshold,
      masked_misalignment_threshold = misalignment_threshold,
      masked_adapters = adapters,
      adapters = adapters
    )
  } else {
    out<-generate_misalignment_info(
      read_input_path = fastq_file_or_directory_path,
      read_layout_path = read_layout_path,
      nthreads = nthreads,
      max_sequences = 500000,
      chunk_size = chunk_size
    )
    data.table::fwrite(out$misalignment_threshold, file = paste0(output_directory_path,"/misalignment_threshold.csv"))
    data.table::fwrite(out$masked_misalignment_threshold, file = paste0(output_directory_path,"/masked_misalignment_threshold.csv"))
  }
  
  sigstream(input_path = fastq_file_or_directory_path, 
            nthreads = nthreads, 
            compress = compress, 
            write_fastq = TRUE, 
            min_length = 100,
            chunkSize = chunk_size,
            outputPath = paste0(output_directory_path, "/demuxed_reads.fastq"),
            adapters = out$masked_adapters,
            misalignment_threshold = out$masked_misalignment_threshold,
            read_layout = read_layout,
            sigstringsFilePath = paste0(output_directory_path, "/sigstrings.txt")
  )
  gc()
  if(generate_sigstring_diagnostics){
    if(sigstring_diagnostic_verbose){
      output_type<-"diagnostic"
    } else {
      output_type<-"summary"
    }
    process_sig(
      file_path = paste0(output_directory_path, "/sigstrings.txt", ifelse(test = compress, yes = ".gz", no = "")),
      output_prefix = output_directory_path,
      max_sequences = tabulated_sigstring_count,
      chunk_size = chunk_size,
      compress = compress,
      nthreads = nthreads,
      output_type = output_type
    )
    gc()
  } else {
    cat("Generate sigstring diagnostics set to FALSE--skipping sigstring tabulation.\n")
  }
  dir.create(path = paste0(output_directory_path,"/variable_seqs/"))
  cat("Tabulating barcodes...\n")
  extractable_elements<-read_layout[type == "variable" & class != "read" & direction == "forward", class_id]
  tabulate_variable_sequences(
    input_files = paste0(output_directory_path, "/demuxed_reads.fastq", 
                         ifelse(test = compress, yes = ".gz", no = "")), 
    identity_elements = c(extractable_elements),
    output_prefix = paste0(output_directory_path,"/variable_seqs/"), 
    compress = TRUE)
  barcode_files<-list.files(path = paste0(output_directory_path,"/variable_seqs/"), full.names = TRUE, pattern = "barcode")
  #need to figure out how to loop over multiple barcodes and rewrite the output file
  out<-lapply(X = barcode_files, FUN = function(X){
    process_barcodes(barcode_path = X, read_layout = read_layout)
    #think I need to do it here, because the input file gets chosen internally in correct_barcodes
    correct_barcodes(
      input_file = X,
      read_layout = read_layout,
      mutation_rounds = mutation_rounds,
      bc_correction_mode = bc_correction_mode,
      nthreads = nthreads,
      chunk_size = chunk_size
    )
  })
  end_time = Sys.time()
  cat(paste0("Finished running! Total run time: ", difftime(end_time, start_time, units='mins'), " minutes.\n"))
}

#' @export
generate_synthetic_reads<-function(
  read_layout_path,
  num_cells,
  read_length,
  num_reads = c(1,50)){
  # Prep layout and set key for class_id
  prep_read_layout(read_layout_path, return_env = parent.frame())
  data.table::setkey(read_layout, "class_id")
  read_layout <- read_layout[direction == "forward" & (expected_length > 0 | class == "read"), ]
  data.table::setkey(read_layout, "class_id")
  
  # Initialize layout as a data.table with num_cells rows
  layout <- data.table(cell_id = as.character(1:num_cells))
  
  # Pre-generate the number of reads per cell
  num_reads_per_cell <- sample(seq(num_reads[1], num_reads[2]), num_cells, replace = TRUE)
  
  # Now we generate the static elements (primers, poly_t, etc.)
  for (unique_class in unique(read_layout$class)) {
    # Handle static elements (e.g., primers, poly_t)
    if (unique_class %in% c("forw_primer", "poly_tail","tso", "rev_primer")) {
      print(paste("Dealing with static element:", unique_class))
      for (static_id in unique(read_layout[class == unique_class, class_id])) {
        static_seq<-read_layout[class_id == static_id, seq]
        if(unique_class == "poly_tail"){
          print(static_seq)
          min_reps<-as.numeric(DescTools::StrExtractBetween(x = static_seq, left = "\\{", right = ","))
          print(min_reps)
          poly_seq<-ifelse(
            grepl("A", static_seq, fixed = TRUE), "A",
            ifelse(grepl("C", static_seq, fixed = TRUE), "C",
              ifelse(grepl("T", static_seq, fixed = TRUE), "T",
                ifelse(grepl("G", static_seq, fixed = TRUE), "G",
                  stop("Invalid base detected")))))
          static_seq<-paste0(rep(x = poly_seq, times = ceiling(min_reps*1.5)), collapse = "")
          print(static_seq)
        }
        layout[, paste0(static_id) := static_seq]
      }
    }
    
    # Handle barcodes (per cell)
    if (unique_class == "barcode") {
      print("Dealing with barcode")
      for (barcode_id in unique(read_layout[class == "barcode", class_id])) {
        layout[, paste0(barcode_id) := {
          whitelist_path <- switch(
            read_layout[barcode_id, whitelist],
            "10x_3v1" = system.file(package = "rad", "extdata", "737K-august-2016_bitlist.csv.gz"),
            "10x_3v3" = system.file(package = "rad", "extdata", "3M-february-2018-3v3.txt_bitlist.csv.gz"),
            "10x_3v4" = system.file(package = "rad", "extdata", "3M-3pgex-may-2023.txt_bitlist.csv.gz"),
            "10x_5v1" = system.file(package = "rad", "extdata", "737K-august-2016_bitlist.csv.gz"),
            "10x_5v3" = system.file(package = "rad", "extdata", "3M-5pgex-jan-2023.txt_bitlist.csv.gz"),
            NA_character_
          )
          whitelist <- whitelist_importer(whitelist_path = whitelist_path)
          barcode_sample <- sample(whitelist$whitelist_bcs, size = num_cells, replace = FALSE)
          barcode_length <- read_layout[barcode_id, expected_length]
          bits_to_sequence(barcode_sample, sequence_length = barcode_length)
        }]
      }
      print("done with barcode")
    }
  }
  # Now clone the rows for each barcode according to num_reads_per_cell
  layout<-layout[rep(seq_len(num_cells), num_reads_per_cell)]
  # Now generate UMIs and reads for each read (post-cloning of rows)
  for (unique_class in unique(read_layout$class)) {
    
    # Handle UMIs
    if (unique_class == "umi") {
      print("Dealing with UMI")
      for (umi_id in unique(read_layout[class == "umi", class_id])) {
        umi_length <- as.numeric(read_layout[class_id == umi_id, expected_length])
        if (is.na(umi_length) | umi_length <= 0) stop(paste("Invalid UMI length for", umi_id))
        
        # Generate a unique UMI for each read within each cell
        layout[, paste0(umi_id) := unlist(lapply(1:num_cells, function(cell) {
          replicate(num_reads_per_cell[cell], paste0(sample(c("A", "T", "C", "G"), umi_length, replace = TRUE), collapse = ""))
        }))]
      }
    }
    
    # Handle reads
    if (unique_class == "read") {
      print("Dealing with reads")
      for (read_id in unique(read_layout[class == "read", class_id])) {
        # Generate a unique read for each read within each cell
        layout[, paste0(read_id) := unlist(lapply(1:num_cells, function(cell) {
          replicate(num_reads_per_cell[cell], paste0(sample(c("A", "T", "C", "G"), read_length, replace = TRUE), collapse = ""))
        }))]
      }
    }
  }
  
  #   # Combine multiple barcode and UMI columns if needed
  layout[, combined_barcode := do.call(paste, c(.SD, sep = "_")), .SDcols = grep("barcode", names(layout), value = TRUE)]
  layout[, combined_umi := do.call(paste, c(.SD, sep = "_")), .SDcols = grep("umi", names(layout), value = TRUE)]
  # Create cell_id by combining barcode, UMI, and read number
  layout[, cell_id := paste0(
    combined_barcode, "_",
    combined_umi, "_",
    sequence(.N)  # Reset numbering within each group
  ), by = .(combined_barcode)]  # Group by combined_barcode
  
  # Reshape and return the final read layout
  ordered_class_ids<-read_layout[order(order), class_id]
  reshaped_layout<-tidytable::unite(layout, final_read, all_of(ordered_class_ids), sep = "")
  reshaped_layout <- reshaped_layout[, .(id = cell_id, seq = final_read)]
  return(reshaped_layout)
}


#' @export
blaze_data_processor<-function(
  path_to_matched_fastq, 
  path_to_putative_bcs
  ){
  putative_bcs<-data.table::fread(path_to_putative_bcs, na.strings = "", select = 2) #2 is the putative_bc column
  putative_bcs<-putative_bcs[, .N, by = putative_bc]
  data.table::setnames(putative_bcs, "N", "original")
  final_counts<-extract_blaze_barcode(path_to_matched_fastq)
  final_counts<-as.data.table(final_counts)
  merged_barcodes<-final_counts[putative_bcs, on = c("barcode" = "putative_bc"), nomatch = 0L]
  merged_barcodes$corrected<-merged_barcodes$total-merged_barcodes$original
  return(merged_barcodes)
}

#' @export
generate_misalignment_info<-function(read_input_path,
                                     read_layout_path,
                                     nthreads = 8,
                                     max_sequences = 500000,
                                     chunk_size = 50000){
  prep_read_layout(read_layout_form = read_layout_path)
  misalignment_data<-n_masked_misalignment_cigar_string(
    input_path = read_input_path, 
    adapters = adapters, 
    nthreads = nthreads, 
    max_sequences = max_sequences, 
    chunk_size = chunk_size)
  data.table::setDT(misalignment_data)
  misalignment_threshold<-stat_collector_new(df = misalignment_data, read_layout = read_layout, mode = 'stats')
  data.table::setDT(misalignment_threshold, key = 'query_id')
  misalignment_threshold[, lower_misal_bounds := round(misal_threshold - misal_sd)]
  out<-sapply(misalignment_threshold[,query_id], FUN = function(adapter){
    seq<-misalignment_data[query_id == adapter & best_edit_distance < 
                             misalignment_threshold[adapter, lower_misal_bounds], aligned_sequence]
    full_consensus<-cmv_consensusMatrix(seq)
    seq<-full_consensus %>% 
      {.[,c(1:stringr::str_length(adapters[adapter]))]} %>%
      {prop.table(., margin = 2) * 100} %>%
      # Apply function to get consensus with 'N' for non-dominant bases
      {apply(., MARGIN = 2, FUN = function(X) {
        if (max(X) > 50) {
          rownames(full_consensus)[which.max(X)]
        } else {
          "N"
        }
      }, simplify = TRUE)} %>%
      paste0(., collapse = '')
    return(list(seq = seq, full_consensus = full_consensus))
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  masked_adapters<-sapply(out, FUN = function(x){x$seq}, USE.NAMES = TRUE, simplify = TRUE)
  masked_adapters<-append(masked_adapters,adapters[-which(names(adapters) %in% names(masked_adapters))])
  
  masked_misalignment_threshold<-n_masked_misalignment_cigar_string(
    input_path = read_input_path, 
    adapters = masked_adapters, 
    nthreads = nthreads, 
    max_sequences = max_sequences, 
    chunk_size = chunk_size) %>% {stat_collector_new(df = ., read_layout = read_layout, mode = 'stats')}
  
  return(
      list(
        misalignment_threshold = misalignment_threshold,
        masked_misalignment_threshold = masked_misalignment_threshold,
        masked_adapters = masked_adapters,
        adapters = adapters
      )
  )
}

#' @export
stat_collector_new<-function(df, read_layout, mode = "stats") {
  data.table::setDT(df)
  data.table::setDT(read_layout)
  adapters <- read_layout[type == "static" & 
      direction %in% c("forward", "reverse") & 
      !class %in% c("poly_tail", "start", "stop"),
    .(class_id, direction)]
  
  forward_adapters <- adapters[direction == "forward", class_id]
  reverse_adapters <- adapters[direction == "reverse", class_id]
  
  # Find IDs with good alignments (edit distance <= 1) more efficiently
  forward_zero_ids <- df[query_id %in% forward_adapters, 
    if(all(best_edit_distance == 0)) .SD[1], 
    by = id][, id]
  
  reverse_zero_ids <- df[query_id %in% reverse_adapters,
    if(all(best_edit_distance == 0)) .SD[1],
    by = id][, id]
  
  if(mode == "stats") {
    # Calculate misalignment stats (wrong direction)
    stats_reverse <- df[id %in% forward_zero_ids & 
        query_id %in% reverse_adapters,
      .(misal_threshold = mean(best_edit_distance),
        misal_sd = sd(best_edit_distance),
        misal_mean_start = mean(best_start_pos/sequence_length)*100,
        misal_mean_stop = mean(best_stop_pos/sequence_length)*100,
        misal_sd_start = sd(best_start_pos/sequence_length)*100,
        misal_sd_stop = sd(best_stop_pos/sequence_length)*100),
      by = query_id]
    
    stats_forward <- df[id %in% reverse_zero_ids & 
        query_id %in% forward_adapters,
      .(misal_threshold = mean(best_edit_distance),
        misal_sd = sd(best_edit_distance),
        misal_mean_start = mean(best_start_pos/sequence_length)*100,
        misal_mean_stop = mean(best_stop_pos/sequence_length)*100,
        misal_sd_start = sd(best_start_pos/sequence_length)*100,
        misal_sd_stop = sd(best_stop_pos/sequence_length)*100),
      by = query_id]
    
    # Calculate correct alignment stats (right direction)
    correct_stats_forward <- df[id %in% forward_zero_ids & 
        query_id %in% forward_adapters,
      .(align_mean_start = mean(best_start_pos/sequence_length)*100,
        align_mean_stop = mean(best_stop_pos/sequence_length)*100,
        align_sd_start = sd(best_start_pos/sequence_length)*100,
        align_sd_stop = sd(best_stop_pos/sequence_length)*100),
      by = query_id]
    
    correct_stats_reverse <- df[id %in% reverse_zero_ids & 
        query_id %in% reverse_adapters,
      .(align_mean_start = mean(best_start_pos/sequence_length)*100,
        align_mean_stop = mean(best_stop_pos/sequence_length)*100,
        align_sd_start = sd(best_start_pos/sequence_length)*100,
        align_sd_stop = sd(best_stop_pos/sequence_length)*100),
      by = query_id]
    
    # Combine misalignment results
    misalignment_thresholds <- data.table::rbindlist(list(
      stats_forward[, direction := "forward"],
      stats_reverse[, direction := "reverse"]
    ))
    
    # Add correct alignment stats
    correct_alignments <- data.table::rbindlist(list(
      correct_stats_forward[, direction := "forward"],
      correct_stats_reverse[, direction := "reverse"]
    ))
    
    # Merge misalignment and correct alignment stats
    results <- merge(
      misalignment_thresholds,
      correct_alignments,
      by = c("query_id", "direction"),
      all = TRUE
    )
    
    # Add metadata
    results <- merge(
      results,
      read_layout[, .(class_id, class, type)],
      by.x = "query_id",
      by.y = "class_id",
      all.x = TRUE
    )
    
    return(results)
    
  } else if(mode == "graphics") {
    out<-data.table::rbindlist(list(
      # Misaligned sequences
      df[id %in% reverse_zero_ids & 
          query_id %in% forward_adapters, 
        .SD, 
        .SDcols = c("id", "query_id", "best_edit_distance", 
          "best_start_pos", "best_stop_pos", 
          "sequence_length",'cigar_matches')
      ][, `:=`(direction = "forward", alignment_type = "misaligned")],
      
      df[id %in% forward_zero_ids & 
          query_id %in% reverse_adapters,
        .SD,
        .SDcols = c("id", "query_id", "best_edit_distance", 
          "best_start_pos", "best_stop_pos", 
          "sequence_length",'cigar_matches')
      ][, `:=`(direction = "reverse", alignment_type = "misaligned")],
      
      # Correctly aligned sequences
      df[id %in% forward_zero_ids & 
          query_id %in% forward_adapters,
        .SD,
        .SDcols = c("id", "query_id", "best_edit_distance", 
          "best_start_pos", "best_stop_pos", 
          "sequence_length",'cigar_matches')
      ][, `:=`(direction = "forward", alignment_type = "correct")],
      
      df[id %in% reverse_zero_ids & 
          query_id %in% reverse_adapters,
        .SD,
        .SDcols = c("id", "query_id", "best_edit_distance", 
          "best_start_pos", "best_stop_pos", 
          "sequence_length",'cigar_matches')
      ][, `:=`(direction = "reverse", alignment_type = "correct")]
    ))
    out[,'best_start_perc' := (best_start_pos/sequence_length)*100]
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
