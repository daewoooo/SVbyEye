#' This function will take in a \code{\link{GRanges-class}} or \code{\link{GRangesList-class}} object of alignments of a 
#' single contig to a reference and collapses alignment gaps based user specified maximum allowed gap.
#'
#' @param ranges A \code{\link{GRanges-class}} or \code{\link{GRangesList-class}} object of regions of a single or multiple contigs aligned to a reference.
#' @param max.gap A maximum length of a gap within a single contig alignments to be collapsed.
#' @return A \code{\link{GRanges-class}} object.
#' @author David Porubsky
#' @export
#' 
collapseGaps <- function(ranges, max.gap=100000) {
  ## Helper function definition
  fillGaps <- function(gr, max.gap=100000) {
    gr <- GenomeInfoDb::keepSeqlevels(gr, value = as.character(unique(GenomeInfoDb::seqnames(gr))), pruning.mode = 'coarse')
    gr <- GenomicRanges::sort(gr)
    gap.gr <- GenomicRanges::gaps(gr, start = min(start(gr)))
    gap.gr <- gap.gr[width(gap.gr) <= max.gap]
    if (length(gap.gr) > 0) {
      red.gr <- GenomicRanges::reduce(c(gr[,0], gap.gr))
      GenomicRanges::mcols(red.gr) <- GenomicRanges::mcols(gr)[length(gr),]
    } else {
      red.gr <- GenomicRanges::reduce(gr)
      GenomicRanges::mcols(red.gr) <- GenomicRanges::mcols(gr)[length(gr),]
    }
    return(red.gr)
  }
  
  if (class(ranges) == "CompressedGRangesList") {
    ## Process only contigs with split alignments
    to.fill <- which(lengths(ranges) > 1)
    red.ranges <-  suppressWarnings( S4Vectors::endoapply(ranges[to.fill], function(gr) fillGaps(gr=gr, max.gap=max.gap)) )
    red.ranges <- unlist(red.ranges, use.names = FALSE)
    ## Add ranges with no split alignment
    red.ranges <- c(unlist(ranges[-to.fill], use.names = FALSE), red.ranges)
  } else if (class(ranges) == "GRanges") {
    red.ranges <- fillGaps(gr=ranges, max.gap = max.gap)
  } else {
    stop("Only objects of class 'GRanges' or 'GRangesList' are allowed as input for parameter 'ranges' !!!")
  }
  return(red.ranges)
}

#' This function will take in a \code{\link{GRanges-class}} or \code{\link{GRangesList-class}} object of genomic ranges of a 
#' single a contig aligned to a reference and reports all alignment gaps.
#'
#' @param id.col A column number from the original \code{\link{GRanges-class}} object to be reported as an unique ID.
#' @inheritParams collapseGaps
#' @return A \code{\link{GRanges-class}} object.
#' @author David Porubsky
#' @export
#' 
reportGaps <- function(ranges, id.col=NULL) {
  ## Helper function definitions
  getGaps <- function(gr=NULL) {
    gap.gr <- GenomicRanges::gaps(gr, start = min(start(gr)))
    gap.gr <- gap.gr[GenomicRanges::strand(gap.gr) == '*']
    return(gap.gr)
  }  
  
  processGaps <- function(gr, id.col=NULL) {
    ## Make sure only seqlevels present in the submitted ranges are kept
    gr <- GenomeInfoDb::keepSeqlevels(gr, value = as.character(unique(GenomeInfoDb::seqnames(gr))), pruning.mode = 'coarse')
    ## Sort ranges by position
    GenomicRanges::strand(gr) <- '*'
    gr <- GenomicRanges::sort(gr)
    ## Make sure rows are not named
    names(gr) <- NULL
    ## Keep only ranges with the same seqnames
    if (length(GenomeInfoDb::seqlevels(gr)) > 1) {
      max.seqname <- GenomeInfoDb::seqlevels(gr)[which.max(S4Vectors::runLength(GenomeInfoDb::seqnames(gr)))]
      gr <- gr[GenomeInfoDb::seqnames(gr) == max.seqname]
      warning("Multiple 'seqlevels' present in submitted ranges, keeping only ranges for: ", max.seqname)
    }
    ## TODO for alignment landing on different chromosomes report gap == 0 and both alignments
    if (length(gr) > 1) {
      ## Calculate gaps
      #gap.gr <- GenomicRanges::gaps(gr, start = min(start(gr)))
      #gap.gr <- gap.gr[strand(gap.gr) == '*']
      gap.gr <- getGaps(gr)
      
      ## Add ID column from the original gr object if defined
      if (!is.null(id.col)) {
        if (id.col > 0 & ncol(GenomicRanges::mcols(gr)) >= id.col) {
          GenomicRanges::mcols(gap.gr) <- rep(unique(GenomicRanges::mcols(gr)[id.col]), length(gap.gr))
        } else {
          warning("User defined 'id.col' number is larger the then total number of columns in input 'gr', skipping adding id column ...")
        }
      }
            
      ## Report alignment on each side of the gap
      if (length(gap.gr) > 0) {
        ## Report upstream and downstream target ranges
        ## Upstream
        up.idx <- IRanges::follow(gap.gr, gr)
        up.idx.keep <- !is.na(up.idx)
        gap.gr$up.gr <- GenomicRanges::GRanges(seqnames=GenomeInfoDb::seqnames(gap.gr), 
                                               ranges=IRanges::IRanges(start=GenomicRanges::start(gap.gr), end=GenomicRanges::start(gap.gr)))
        gap.gr$up.gr[up.idx.keep] <- gr[up.idx[up.idx.keep]][,0]
        ## Downstream
        down.idx <- IRanges::precede(gap.gr, gr)
        down.idx.keep <- !is.na(down.idx)
        gap.gr$down.gr <- GenomicRanges::GRanges(seqnames=GenomeInfoDb::seqnames(gap.gr), 
                                                 ranges=IRanges::IRanges(start=GenomicRanges::end(gap.gr), end=GenomicRanges::end(gap.gr)))
        gap.gr$down.gr[down.idx.keep] <- gr[down.idx[down.idx.keep]][,0]
        
        ## ## Report upstream and downstream query ranges and gaps if defined
        if ('query.gr' %in% names(mcols(gr)) & class(gr$query.gr) == 'GRanges') {
          ## Upstream
          gap.gr$query.up.gr <- GenomicRanges::GRanges(seqnames=GenomeInfoDb::seqnames(gap.gr), 
                                                       ranges=IRanges::IRanges(start=GenomicRanges::start(gap.gr), end=GenomicRanges::start(gap.gr)))
          suppressWarnings( gap.gr$query.up.gr[up.idx.keep] <- gr$query.gr[up.idx[up.idx.keep]][,0] )
          ## Downstream
          gap.gr$query.down.gr <- GenomicRanges::GRanges(seqnames=GenomeInfoDb::seqnames(gap.gr), 
                                                         ranges=IRanges(start=GenomicRanges::end(gap.gr), end=GenomicRanges::end(gap.gr)))
          suppressWarnings( gap.gr$query.down.gr[down.idx.keep] <- gr$query.gr[down.idx[down.idx.keep]][,0] )
          ## Calculate gaps
          gap.grl <- split(gap.gr, 1:length(gap.gr))
          query.gap.grl <- S4Vectors::endoapply(gap.grl, function(x) getGaps(c(x$query.up.gr, x$query.down.gr)))
          query.gap.gr <- unlist(query.gap.grl, use.names = FALSE)
          ## Initialize gap ranges
          gap.gr$query.gap.gr <- GenomicRanges::GRanges(seqnames=rep('unknown', length(gap.gr)), 
                                                        ranges=IRanges::IRanges(start=1, end=1))
          if (length(query.gap.gr) > 0) {
            gap.idx <- which(lengths(query.gap.grl) > 0)
            gap.gr$query.gap.gr[gap.idx] <- query.gap.gr
          }
        }
        return(gap.gr)
      } else {
        dummy.gr <- GRanges(seqnames = 'dummy', ranges = IRanges(start=1, end=1), ID='dummy')
        dummy.gr$up.gr <- GRanges(seqnames = 'dummy', ranges = IRanges(start=1, end=1))
        dummy.gr$down.gr <- GRanges(seqnames = 'dummy', ranges = IRanges(start=1, end=1))
        dummy.gr$query.up.gr <- GRanges(seqnames = 'dummy', ranges = IRanges(start=1, end=1))
        dummy.gr$query.down.gr <- GRanges(seqnames = 'dummy', ranges = IRanges(start=1, end=1))
        dummy.gr$query.gap.gr <- GRanges(seqnames = 'dummy', ranges = IRanges(start=1, end=1))
        return(dummy.gr)
      }
    } else {
      dummy.gr <- GRanges(seqnames = 'dummy', ranges = IRanges(start=1, end=1), ID='dummy')
      dummy.gr$up.gr <- GRanges(seqnames = 'dummy', ranges = IRanges(start=1, end=1))
      dummy.gr$down.gr <- GRanges(seqnames = 'dummy', ranges = IRanges(start=1, end=1))
      dummy.gr$query.up.gr <- GRanges(seqnames = 'dummy', ranges = IRanges(start=1, end=1))
      dummy.gr$query.down.gr <- GRanges(seqnames = 'dummy', ranges = IRanges(start=1, end=1))
      dummy.gr$query.gap.gr <- GRanges(seqnames = 'dummy', ranges = IRanges(start=1, end=1))
      return(dummy.gr)
      #return(GRanges())
    }
  }  
  
  if (class(ranges) == "CompressedGRangesList") {
    ## Keep only contigs with split alignments
    ptm <- startTimedMessage("Reporting gaps")
    
    ranges <- ranges[lengths(ranges) > 1]
    gaps <-  suppressWarnings( S4Vectors::endoapply(ranges, function(gr) processGaps(gr=gr, id.col=id.col)) )
    gaps <- unlist(gaps, use.names = FALSE)
    ## Remove empty ranges
    gaps <- gaps[seqnames(gaps) != 'dummy']
    #gaps <-  suppressWarnings( S4Vectors::lapply(ranges, function(gr) processGaps(gr=gr, id.col=id.col)) )
    #do.call(c, gaps)
    
    stopTimedMessage(ptm)
  } else if (class(ranges) == "GRanges") {
    ptm <- startTimedMessage("Reporting gaps")
    
    gaps <- processGaps(gr=ranges, id.col=id.col)
    
    stopTimedMessage(ptm)
  } else {
    stop("Only objescts of class 'GRanges' or 'GRangesList' are allowed as input for parameter 'ranges' !!!")
  }
  return(gaps)
}


#' This function will take in a \code{\link{GRanges-class}} or \code{\link{GRangesList-class}} object of genomic ranges of a 
#' single a contig aligned to a reference and reports coverage of regions that overlaps each other. 
#'
#' @inheritParams collapseGaps
#' @inheritParams reportGaps
#' @return A \code{\link{GRanges-class}} object.
#' @author David Porubsky
#' @export
#' 
reportContigCoverage <- function(ranges, id.col=NULL) {
  ## Helper function definitions
  processContigCoverage <- function(gr=NULL, id.col=NULL) {
    ## Get disjoin ranges
    gr.disjoint <- GenomicRanges::disjoin(gr, ignore.strand=TRUE)
    ## Add ID column from the original gr object if defined
    if (!is.null(id.col)) {
      if (id.col > 0 & ncol(GenomicRanges::mcols(gr)) >= id.col) {
        GenomicRanges::mcols(gr.disjoint) <- rep(unique(GenomicRanges::mcols(gr)[id.col]), length(gr.disjoint))
      } else {
        warning("User defined 'id.col' number is larger the then total number of columns in input 'gr', skipping adding id column ...")
      }
    }
    ## Report coverage of each disjoint range
    gr.disjoint$cov <- IRanges::countOverlaps(gr.disjoint, gr)
    ## Export final ranges with coverage
    return(gr.disjoint)
  }
  
  if (class(ranges) == "CompressedGRangesList") {
    covs <-  suppressWarnings( S4Vectors::endoapply(ranges, function(gr) processContigCoverage(gr=gr, id.col=id.col)) )
    covs <- unlist(covs, use.names = FALSE)
  } else if (class(ranges) == "GRanges") {
    covs <- processContigCoverage(gr=ranges, id.col=id.col)
  } else {
    stop("Only objects of class 'GRanges' or 'GRangesList' are allowed as input for parameter 'ranges' !!!")
  }
  return(covs)
}

#' This function will take in a \code{\link{GRanges-class}} object of genomic ranges and enumerate the number of overlapping bases
#' with other set of user defined 'query' genomic ranges.
#'
#' @param ranges A \code{\link{GRanges-class}} object of genomic regions to get the number of overlapping bases with 'query.ranges'.
#' @param query.ranges A \code{\link{GRanges-class}} object of genomic regions for which one want to report overlaps. 
#' @param index A user defined name of a column where the number of overlapping bases between 'ranges' and 'query.ranges' will be reported.
#' @return A \code{\link{GRanges-class}} object.
#' @author David Porubsky
#' @export
#'
reportOverlapBases <- function(ranges=NULL, query.ranges=NULL, index=NULL) {
  if (!is.null(query.ranges) & !is.null(ranges)) {
    if (class(ranges) == "GRanges" & class(query.ranges) == "GRanges") {
      ## Get disjoint ranges between query and user.defined set of ranges
      disj.gr <- suppressWarnings( GenomicRanges::disjoin(c(ranges[,0], query.ranges[,0])) )
      ## Get disjoint ranges that are overlapping with query.ranges
      disj.query.gr <- IRanges::subsetByOverlaps(disj.gr, query.ranges)
      ## Get disjoint ranges overlapping with ranges of interest
      disj.query.roi.gr <- IRanges::subsetByOverlaps(disj.query.gr, ranges)
      ## Split by regions of interest
      hits <- IRanges::findOverlaps(ranges, disj.query.roi.gr)
      disj.query.roi.grl <- split(disj.query.roi.gr[S4Vectors::subjectHits(hits)], S4Vectors::queryHits(hits))
      query.bases <- sapply(disj.query.roi.grl, function(gr) sum(width(reduce(gr))))
      ## Add overlapping bases counts
      if (!is.null(index) & nchar(index) > 0) {
        new.col.idx <- ncol(GenomicRanges::mcols(ranges)) + 1
        GenomicRanges::mcols(ranges)[new.col.idx] <- 0
        colnames(GenomicRanges::mcols(ranges))[new.col.idx] <- index
        GenomicRanges::mcols(ranges)[new.col.idx][unique(S4Vectors::queryHits(hits)),] <- query.bases
      } else {
        ranges$query.bases <- 0
        ranges$query.bases[unique(S4Vectors::queryHits(hits))] <- query.bases
      } 
    }
  }
  return(ranges)
}  
      
