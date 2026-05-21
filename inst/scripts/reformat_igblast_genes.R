# Normalize igblast V gene calls by collapsing alleles and OR annotations.

#' Reformat igblast gene calls
#' @param v_calls character vector (e.g., nolim_receptor_data$bcr_v_call)
#' @param mc_cores number of cores for pbmcmapply
#' @return character vector of reformatted gene calls
#' @export
reformat_igblast_genes <- function(v_calls, mc_cores = 4) {
  val <- v_calls %>%
    gsub("\\*..", replacement = "", x = .) %>%
    gsub("\\|OR(\\d+)-\\d+", "-OR\\1", .) %>%
    gsub("D-", "-", .) %>%
    strsplit(",") %>% lapply(unique, simplify = TRUE, USE.NAMES = FALSE)

  pbmcmapply(X = val, mc.cores = mc_cores, FUN = function(X){
    if(any(is.na(X))){
      return(NA)
    }
    if(length(X) <= 1){
      return(X)
    }
    root <- gsub(pattern = "(?<=-).+$", x = X, replacement = "", perl = TRUE)
    if((length(X) > 1) & any(duplicated(root))){
      dup_index <- duplicated(root) | duplicated(root, fromLast = TRUE)
      if(all(dup_index)){
        root <- unique(root)
        stripped_val <- gsub(pattern = root, x = X, replacement = "", perl = TRUE)
        stripped_val <- paste0(stripped_val, collapse = "|")
      } else {
        root <- unique(root[dup_index])
        stripped_val <- gsub(pattern = root, x = X[dup_index], replacement = "", perl = TRUE)
      }
      stripped_val <- paste0(stripped_val, collapse = "|")
      collapsed <- paste0(root, stripped_val)
      if(all(dup_index)){
        return(collapsed)
      }
      if(length(X) > 1 & !all(dup_index)){
        out <- paste(c(collapsed, X[!dup_index]), collapse = ":")
        return(out)
      }
    } else if(length(X) > 1){
      out <- paste0(X, collapse = ":")
      return(out)
    }
  })
}
