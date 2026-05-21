# Rewrite FASTQ ids using a truth CSV (sctagger synthetic data utility).

#' Rewrite FASTQ IDs using a truth CSV
#' @param input_fq input FASTQ path
#' @param truth_csv CSV containing true IDs (assumes second column)
#' @param output_fq output FASTQ path
#' @export
generate_sctagger_ids <- function(input_fq, truth_csv, output_fq){
  fqs <- SIRAP::read_fastqas(fn = input_fq)
  truth <- data.table::fread(truth_csv)
  fqs[, id := truth$V2]
  fqs[, id := paste0(id, "-", seq_len(.N)), by = id]
  SIRAP::write_fastqas(df = fqs, fn = output_fq)
}
