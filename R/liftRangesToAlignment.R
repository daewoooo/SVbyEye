#' Function to lift coordinates to the alignment in PAF format.
#'
#' @param gr A \code{\link{GRanges-class}} object containing single or multiple ranges in query or target sequence coordinates.
#' @param paf.file A path to a PAF file containing alignments of a query sequence to a target sequence.
#' @param direction One of the possible, lift ranges from query to target 'query2target' or vice versa 'target2query'.
#' @importFrom GenomicRanges GRanges findOverlaps mcols
#' @importFrom GenomicAlignments GAlignments mapFromAlignments
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom methods is
#' @return A \code{\link{GRanges-class}} object with resized original set of ranges. 
#' @author David Porubsky
#' @export
#' @examples 
#'## Define range(s) to lift
#'roi.gr <- as('chr17:46645907-46697277', 'GRanges')
#'## Get PAF alignment to lift to
#'paf.file <- system.file("extdata", "test_lift1.paf", package="SVbyEye")
#'## Lift target range to query coordinates
#'liftRangesToAlignment(gr = roi.gr, paf.file = paf.file, direction = 'target2query')
#'
liftRangesToAlignment <- function(gr=NULL, paf.file=NULL, direction='query2target') {
  ## Check user input
  stopifnot(methods::is(gr, "GRanges"), file.exists(paf.file), methods::is(direction, 'character'))
  
  ## Read PAF alignments to reference
  paf.aln <- readPaf(paf.file = paf.file, restrict.paf.tags = 'cg')
  ## Make sure defined ranges are present in the PAF alignment
  if (direction == 'query2target') {
    paf.query.gr <- GenomicRanges::GRanges(seqnames=paf.aln$q.name, ranges = IRanges::IRanges(start=paf.aln$q.start, end=paf.aln$q.end))
    hits <- suppressWarnings( GenomicRanges::findOverlaps(gr, paf.query.gr) )
  } else if (direction == 'target2query') {
    paf.target.gr <- GenomicRanges::GRanges(seqnames=paf.aln$t.name, ranges = IRanges::IRanges(start=paf.aln$t.start, end=paf.aln$t.end))
    hits <- suppressWarnings( GenomicRanges::findOverlaps(gr, paf.target.gr) )
  } else {
    stop("Parameter 'direction' can take only values: 'query2target' or 'target2query'")
  } 
  
  if (length(hits) == 0) {
    message('None of the user defined ranges are present in the PAF query coordinates!!!')
    ## Return zero positioned ranges if there is no overlap
    gr.lifted <- GenomicRanges::GRanges(seqnames = rep('Failed', length(gr)), ranges = IRanges::IRanges(start=1, end=0))
    return(gr.lifted)
  } else {
    ## Lift ranges to PAF alignment
    if (direction == 'query2target') {
      ## Set names to ranges to be lifted
      gr.lift <- gr[S4Vectors::queryHits(hits)]
      names(gr.lift) <- paste0('aln', S4Vectors::subjectHits(hits))
      ## Define PAF alignment object
      paf.aln <- paf.aln[unique(S4Vectors::subjectHits(hits)),]
      alignment <- GenomicAlignments::GAlignments(seqnames = paf.aln$t.name, 
                                                  pos = as.integer(paf.aln$t.start) + 1L, 
                                                  cigar = paf.aln$cg, 
                                                  strand=GenomicRanges::strand(paf.aln$strand), 
                                                  names = paste0('aln', 1:nrow(paf.aln)))
      gr.lifted <- GenomicAlignments::mapFromAlignments(x = gr.lift, alignments = alignment)
      names(gr.lifted) <- NULL
    } else {
      gr.lift <- gr
      ## Define PAF alignment object
      paf.aln <- paf.aln[unique(S4Vectors::subjectHits(hits)),]
      alignment <- GenomicAlignments::GAlignments(seqnames = paf.aln$t.name, 
                                                  pos = as.integer(paf.aln$t.start) + 1L, 
                                                  cigar = paf.aln$cg, 
                                                  strand=GenomicRanges::strand(paf.aln$strand), 
                                                  names = paf.aln$q.name)
      
      gr.lifted <- GenomicAlignments::mapToAlignments(x = gr.lift, alignments = alignment)
      names(gr.lifted) <- NULL
      ## Flip coordinates in case of reverse alignment
      if (any(S4Vectors::runValue(GenomicRanges::strand(alignment)) == '-')) {
        #start.reverse <- (paf.aln$q.end[gr.lifted$alignmentsHits] - end(gr.lifted)) + 1
        #end.reverse <- (paf.aln$q.end[gr.lifted$alignmentsHits] - start(gr.lifted)) + 1
        #ranges(gr.lifted) <- IRanges::IRanges(start=start.reverse, end=end.reverse)
        ## Get alignments to flip
        aln.id <- which(S4Vectors::runValue(GenomicRanges::strand(alignment)) == '-')
        bounds <- IRanges::IRanges(start = 1L, end = paf.aln$q.end[aln.id])
        mask <- gr.lifted$alignmentsHits == aln.id
        ## Flip reverse alignments
        IRanges::ranges(gr.lifted[mask]) <- IRanges::reflect(x = IRanges::ranges(gr.lifted[mask]), bounds = bounds)
      }
    }  
    
    ## Add index corresponding to original input ranges
    #gr.lifted$idx <- S4Vectors::queryHits(hits)[gr.lifted$xHits]
    ## Add empty ranges to coordinates that failed to lift
    if (!all(1:length(gr) %in% gr.lifted$xHits)) {
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
      GenomicRanges::mcols(gr.lifted) <- c(GenomicRanges::mcols(gr.lifted), GenomicRanges::mcols(gr))
    } 
    
    ## Return lifted coordinates
    return(gr.lifted)
  }  
}
