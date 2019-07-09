readKallistoh5 <- function(fpath, ...) {
  if (!requireNamespace("rhdf5", quietly=TRUE)) {
    stop("reading kallisto results from hdf5 files requires Bioconductor package `rhdf5`")
  }
  counts <- rhdf5::h5read(fpath, "est_counts")
  ids <- rhdf5::h5read(fpath, "aux/ids")
  efflens <- rhdf5::h5read(fpath, "aux/eff_lengths")
  lens <- rhdf5::h5read(fpath, "aux/lengths")
  # as suggested by https://support.bioconductor.org/p/96958/#101090
  ids <- as.character(ids)
  
  stopifnot(length(counts) == length(ids)) 
  stopifnot(length(efflens) == length(ids))
  
  result <- data.frame(target_id = ids,
                       length=lens,
                       eff_length = efflens,
                       est_counts = counts,
                       stringsAsFactors = FALSE)
  normfac <- with(result, (1e6)/sum(est_counts/eff_length))
  result$tpm <- with(result, normfac*(est_counts/eff_length))
  return(result)
}

transInh5 <- function(fpath){
  
}
