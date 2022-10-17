#' Function to lift local sequence coordinates to the alignment of a given sequence to the reference.
#'
#' @param gr A \code{\link{GRanges-class}} object containing single or multiple ranges in local sequence coordinates.
#' @param paf.file A path to a PAF file containing alignments of a local sequence to the reference sequence.
#' @importFrom GenomicRanges GRanges findOverlaps mcols
#' @importFrom GenomicAlignments GAlignments mapFromAlignments
#' @importFrom S4Vectors queryHits subjectHits
#' @return A \code{\link{GRanges-class}} object with resized original set of ranges. 
#' @author David Porubsky
#'
liftRangesToAlignment_deprecated <- function(gr=NULL, paf.file=NULL) {
  ## Read PAF alignments to reference
  paf.aln <- readPaf(paf.file = paf.file, restrict.paf.tags = 'cg')
  ## Make sure defined ranges are present in the PAF alignment
  paf.query.gr <- GenomicRanges::GRanges(seqnames=paf.aln$q.name, ranges = IRanges::IRanges(start=paf.aln$q.start, end=paf.aln$q.end))
  hits <- suppressWarnings( GenomicRanges::findOverlaps(gr, paf.query.gr) )
  if (length(hits) == 0) {
    message('None of the user defined ranges are present in the PAF query coordinates!!!')
    ## Return zero positioned ranges if there is no overlap
    #gr.lifted <- GenomicRanges::GRanges(seqnames = GenomeInfoDb::seqnames(gr), ranges = IRanges::IRanges(start=1, end=0))
    gr.lifted <- GenomicRanges::GRanges(seqnames = rep('Failed', length(gr)), ranges = IRanges::IRanges(start=1, end=0))
    return(gr.lifted)
  } else {
    gr.lift <- gr[S4Vectors::queryHits(hits)]
    names(gr.lift) <- paste0('aln', S4Vectors::subjectHits(hits))
    
    ## Define PAF alignment object
    alignment <- GenomicAlignments::GAlignments(seqnames = paf.aln$t.name, 
                                                pos = as.integer(paf.aln$t.start), 
                                                cigar = paf.aln$cg, 
                                                strand=GenomicRanges::strand(paf.aln$strand), 
                                                names = paste0('aln', 1:nrow(paf.aln)))
    ## Lift ranges to PAF alignment
    gr.lifted <- GenomicAlignments::mapFromAlignments(x = gr.lift, alignments = alignment)
    names(gr.lifted) <- NULL
    ## Add index corresponding to original input ranges
    gr.lifted$idx <- S4Vectors::queryHits(hits)[gr.lifted$xHits]
    ## Add empty ranges to coordinates that failed to lift
    if (!all(1:length(gr) %in% gr.lifted$idx)) {
      failed.idx <- (1:length(gr))[-gr.lifted$idx]
      decoy.gr <- GenomicRanges::GRanges(seqnames = rep('Failed', length(failed.idx)), 
                                         ranges = IRanges::IRanges(start=1, end=0),
                                         xHits = 1L,
                                         alignmentsHits = 1L,
                                         idx = failed.idx)
      gr.lifted <- suppressWarnings( c(gr.lifted, decoy.gr) )
      gr.lifted <- gr.lifted[order(gr.lifted$idx)]
    }  
    ## Keep extra metacolumns if present
    if (ncol(GenomicRanges::mcols(gr)) > 0) {
      GenomicRanges::mcols(gr.lifted) <- GenomicRanges::mcols(gr)
    } else {
      gr.lifted <- gr.lifted[,0]
    } 
    
    ## Return lifted coordinates
    return(gr.lifted)
  }  
}
