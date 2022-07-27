#' Read PAF alignments from a file and apply a set of post-processing steps.
#' 
#' This function takes an PAF output file (e.g. minimap2 output) and loads the file and applies 
#' desired CIGAR operations and filtering steps.
#'
#' @param min.align.len Filter our PAF alignemnt smaller than this size.
#' @param min.selfaln.dist Keep alignment pairs with this or larger distance from each other [Applied only for FASTA self-alignments].
#' @param collapse.overlaps Set to \code{TRUE} to merge overlapping pair of alignments with the same relative orientation.
#' @param break.paf.aln Set to \code{TRUE} in order to split CIGAR string at insertions and deletions.
#' @param bin.paf.aln Set to \code{TRUE} in order to bin each alignment into bins defined in 'binsize' parameter.
#' @param binsize A user defined binsize (in bp) to split each PAF alignment into.
#' @inheritParams readPaf
#' @inheritParams breakPafAlignment
#' @return A \code{list} of \code{tibble} objects storing matched ('M') alignments as well as structurally variable ('SV') regions.
#' @author David Porubsky
#' @export
crunchPaf <- function(paf.file=NULL, min.align.len=1000, min.selfaln.dist=10000, min.deletion.size=1000, min.insertion.size=1000, break.paf.aln=TRUE, bin.paf.aln=FALSE, binsize=10000) {
  ## Check user input ##
  if (!nchar(paf.file) > 0 | !file.exists(paf.file)) {
    stop("Submitted file in 'aln.coords' does not exists !!!")
  }
  
  ## Data input ##
  ## Read in paf file (e.g. minimap2 output)
  paf.data <- readPaf(paf.file = paf.file, restrict.paf.tags = 'cg')
  
  ## Data filtering ##
  ## Filter by alignment length
  if (min.align.len > 0) {
    paf.data <- paf.data[paf.data$aln.len >= min.align.len,]
  }
  ## When processing self-alignments apply separate set of filters
  if (all(paf.data$q.name == paf.data$t.name)) {
    ## Remove diagonals (self-alignments) with the same start and end position
    paf.data <- paf.data[!(paf.data$q.start == paf.data$t.start & paf.data$q.end == paf.data$t.end),]
    ## Due to the minimap2 self-alignment redundancy keep only alignments where query start is smaller than the target start
    paf.data <- paf.data[paf.data$q.start < paf.data$t.start,]
    ## Filter by alignment distance [self-alignments only]
    aln.dist <- pmin(paf.data$t.start, paf.data$t.end) - pmax(paf.data$q.start, paf.data$q.end)
    if (min.selfaln.dist > 0) {
      paf.data <- paf.data[aln.dist >= min.selfaln.dist,]
    }
  }
  
  ## Data post-processing ##
  ## Make sure cigar string ('cg') column is present
  if ('cg' %in% colnames(paf.data)) {
    ## Break paf alignments
    if (break.paf.aln) {
      paf.data <- breakPaf(paf.table = paf.data, min.deletion.size = min.deletion.size, min.insertion.size = min.insertion.size, report.sv = TRUE)
      paf.data.sv <- paf.data$SVs
      paf.data <- paf.data$M
    } else {
      paf.data.sv <- NULL
    }
    ## Bin paf alignments
    if (bin.paf.aln) {
      paf.data <- pafToBins(paf.table = paf.data, binsize = binsize)
    }  
  }  
  
  ## Return paf alignments
  return(list('M' = paf.data, 'SVs' = paf.data.sv))
}