#' Split genome into user defined bins.
#' 
#' This function splid genome into a equaly-sized user defined genomic intervals.
#'
#' @param bsgenome A reference genome to get lengths of genomic sequences (eg. GRCh38).
#' @param fai A FASTA index to get lengths of genomic sequences.
#' @param chromosomes A user defined set of chromosomes for binning (eg. 'chr1')
#' @param binsize A size of the genomic bin to split genome into.
#' @param stepsize A size of the genomic interval to move each bin. For non-overlapping bins use the same size as binsize.
#' @importFrom GenomeInfoDb seqlengths seqlevels
#' @importFrom utils read.table 
#' @importFrom GenomicRanges GRangesList GRanges
#' @importFrom IRanges IRanges
#' @return A \code{\link{GRanges-class}} object with resized original set of ranges. 
#' @author David Porubsky
#' @export
makeBins <- function(bsgenome=NULL, fai=NULL, chromosomes=NULL, binsize=100000, stepsize=binsize/2) {
  
  if (!is.null(bsgenome)) {
    chr.lengths <- GenomeInfoDb::seqlengths(bsgenome)
  } else if (!is.null(fai)) {
    ## Get contigs/scaffolds names and sizes from fasta index
    fai.tab <- utils::read.table(fai)
    fai.tab <- fai.tab[order(fai.tab$V2, decreasing = TRUE),]
    chr.lengths <- fai.tab$V2
    names(chr.lengths) <- fai.tab$V1
    #chr.lengths <- chr.lengths[names(chr.lengths) %in% chromosomes]
  } else {
    warning("Please submit chromosome lengths in a form of BSgenome object or fasta index (.fai)!!!")
  } 
  
  ## Keep only chromosomes larger then the binsize
  chr.lengths <- chr.lengths[chr.lengths > binsize]
  ## Keep only user defined chromosomes
  chroms.in.data <- names(chr.lengths)
  if (!is.null(chromosomes)) {
    chromosomes <- chromosomes[chromosomes %in% chroms.in.data]
  } else {
    chromosomes <- chroms.in.data
  }  
  
  bins <- GenomicRanges::GRangesList()
  GenomeInfoDb::seqlevels(bins) <- chromosomes
  for (i in seq_along(chr.lengths)) {
    chr.len <- chr.lengths[i]
    
    bin.starts <- seq(from = 1, to = chr.len-binsize, by = stepsize)
    bin.ends <- seq(from = binsize, to = chr.len, by = stepsize)
    
    chr.bins <- GenomicRanges::GRanges(seqnames=names(chr.len), ranges=IRanges::IRanges(start=bin.starts, end=bin.ends))
    bins[[i]] <- chr.bins
  }
  bins <- unlist(bins, use.names = FALSE)
  GenomeInfoDb::seqlengths(bins) <- chr.lengths[chromosomes]
  return(bins)
}
