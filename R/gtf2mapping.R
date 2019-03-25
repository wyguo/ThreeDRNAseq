# gtfile <- "C:/Users/WG42010/Box Sync/Clock AS Project/Researchers folders/Runxuan/RTD2_19April2016/atRTD2_19April2016.gtf"
#' Generate transcript-gene mapping from transcriptome gtf file
#' @param gtfile the name of the transcriptome gtf file which the data are to be read from.
#' @return a matrix with first column of transcript names and second column of gene ids.
#' @export
#' 
gtf2mapping <- function(gtfile){
  x <- read.delim(gtfile, header=FALSE,sep = '\t',colClasses = c(rep("NULL", 8), rep("character", 1)))
  x <- unique(as.vector(t(x)))
  x <- base::strsplit(x = x,split = ';')
  x <- unique(x)
  message('Generate transcript-gene mapping...')
  pb <- txtProgressBar(min = 0, max = length(x), style = 3)
  mapping <- lapply(1:length(x),function(i){
    Sys.sleep(0)
    setTxtProgressBar(pb, i)
    y <- x[[i]]
    idx <- c(grep(pattern = 'transcript_id',x = y),grep(pattern = 'gene_id',x = y))
    gsub(pattern = 'transcript_id|gene_id|[[:space:]]',replacement = '',y[idx])
  })
  close(pb)
  mapping <- data.frame(do.call(rbind,mapping),stringsAsFactors = F)
  colnames(mapping) <- c('TXNAME','GENEIND')
  mapping <- mapping[!duplicated(mapping),]
  return(mapping)
}

