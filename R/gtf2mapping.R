# gtfile <- "C:/Users/WG42010/Box Sync/Clock AS Project/Researchers folders/Runxuan/RTD2_19April2016/atRTD2_19April2016.gtf"
#' Generate transcript-gene mapping from transcriptome gtf file
#' @param gtfile the name of the transcriptome gtf file which the data are to be read from.
#' @details Transcript-gene association mapping information will be extracted from the last column of the gtf file. Transcript names
#' are from the tag "transcript_id" and gene IDs are from the tag "gene_id".
#' @return a matrix with first column of transcript names and second column of gene ids.
#' @export
#' 
gtf2mapping <- function(gtfile){
  x <- read.delim(gtfile, header=FALSE,sep = '\t',colClasses = c(rep("NULL", 8), rep("character", 1)))
  x <- unique(as.vector(t(x)))
  start.time <- Sys.time()
  x <- base::strsplit(x = x,split = ';')
  x <- unique(x)
  message('Generate transcript-gene mapping. It may take a while ...')
  mapping <- lapply(1:length(x),function(i){
    y <- x[[i]]
    idx <- c(grep(pattern = 'transcript_id',x = y),grep(pattern = 'gene_id',x = y))
    if(length(idx)!=2){
      NA
    } else {
      gsub(pattern = 'transcript_id|gene_id|[[:space:]]',replacement = '',y[idx])
    }
  })
  mapping <- mapping[!is.na(mapping)]
  mapping <- data.frame(do.call(rbind,mapping),stringsAsFactors = F)
  colnames(mapping) <- c('TXNAME','GENEIND')
  mapping <- mapping[!duplicated(mapping),]
  return(mapping)
}

#' Generate transcript-gene mapping from transcriptome fasta (.fa) file
#' @param pafile the name of the transcriptome .fa file which the data are to be read from.
#' @details Transcript-gene association mapping information will be extracted from the description line of a sequence. 
#' Transcript names are the sequence names after tag ">" or "transcript" while gene IDs are from the sequence names after tag "gene".
#' @return a matrix with first column of transcript names and second column of gene ids.
#' @export
#' 
# pafile <- "N:/scratch/ImmediateEarlyRNAseq/AtRTDv2_QUASI_13April2017.fa"
fa2mapping <- function(pafile){
  Lines <- readLines(pafile)
  idx <- which(substr(Lines, 1L, 1L) == ">")
  x <- Lines[idx]
  mapping <- lapply(x, function(i){
    ann <- unlist(strsplit(i,'[[:space:]]|[|]|[=]|[:]|[>]'))
    tg <- c(ifelse("^transcript$" %in% ann,ann[grep(pattern = 'transcript',x = ann)+1],ann[2]),
      ann[grep(pattern = '^gene$',x = ann)+1])
    if(length(tg)!=2)
      NA else tg
  })
  mapping <- mapping[!is.na(mapping)]
  mapping <- data.frame(do.call(rbind,mapping),row.names = NULL)
  colnames(mapping) <- c('TXNAME','GENEID')
  mapping
}


