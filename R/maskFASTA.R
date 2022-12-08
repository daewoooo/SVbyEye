#' Mask FASTA sequence at defined regions.
#' 
#' @param fasta.file A path to a FASTA file to be masked.
#' @param mask.ranges A \code{\link{IRanges-class}} object of coordinates to be masked in input FASTA.
#' @param invert If set \code{TRUE} ranges defined in 'mask.ranges' will be kept while the rest of the FASTA will be masked.
#' @param mask.character A single character to be used for masking (Default : 'N').
#' @inheritParams paf2FASTA
#' @importFrom Rsamtools scanFa
#' @importFrom Biostrings replaceAt
#' @importFrom GenomicRanges setdiff
#' @importFrom IRanges IRanges ranges reduce
#' @importFrom S4Vectors sapply
#' @author David Porubsky
#' @export
#'
maskFasta <- function(fasta.file, mask.ranges=NULL, invert=FALSE, mask.character='N', fasta.save=NULL) {
  if (file.exists(fasta.file)) {
    ## Read in FASTA sequence
    gr.seq <- Rsamtools::scanFa(file = fasta.file, as = "DNAStringSet")
  } else {
    stop("Input FASTA: '", fasta.file, "' does not exists !!!")
  } 
  
  if (!is.null(mask.ranges)) {
    ## Make sure input ranges are IRanges and nonoverlapping
    mask.ranges <- IRanges::reduce(IRanges::ranges(mask.ranges))
    ## Get ranges to mask
    seq.range <- IRanges::IRanges(start=1, end=width(gr.seq))
    if (invert) {
      mask.at.ranges <- GenomicRanges::setdiff(seq.range, mask.ranges)
    } else {
      mask.at.ranges <- mask.ranges
    } 
    ## Create mask
    mask.n <- S4Vectors::sapply(width(mask.at.ranges), function(x) paste(rep(mask.character, times=x), collapse = ''))
    ## Mask FASTA
    masked.seq <- Biostrings::replaceAt(gr.seq, at = mask.at.ranges, value = mask.n)
  } else {
    masked.seq <- gr.seq
  }
  
  ## Write final FASTA
  if (is.character(fasta.save)) {
    ## Remove comment character from sequence names
    names(masked.seq) <- gsub(names(masked.seq), pattern = '#', replacement = '_')
    Biostrings::writeXStringSet(x = masked.seq, filepath = fasta.save, format = 'fasta')
  } else {
    warning("Please speficify 'fasta.save' if you want to export FASTA into a file!!!")
  } 
  ## Return masked FASTA
  return(masked.seq)
}
