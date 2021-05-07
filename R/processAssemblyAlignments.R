#' This function will take in a \code{\link{GRanges-class}} object of alignments of a single contig to a reference
#' and collapses gaps in alignment based user specified maximum allowed gap.
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
      mcols(red.gr) <- mcols(gr)[length(gr),]
    } else {
      red.gr <- GenomicRanges::reduce(gr)
      mcols(red.gr) <- mcols(gr)[length(gr),]
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
    stop("Only objescts of class 'GRanges' or 'GRangesList' are allowed as input for parameter 'ranges' !!!")
  }
  return(red.ranges)
}

#' This function will take in a \code{\link{GRanges-class}} object of alignments of a single contig to a reference
#' and reports all alignment gaps.
#'
#' @param id.col A column number from the original \code{\link{GRanges-class}} object to be reported as na unique ID.
#' @inheritParams collapseGaps
#' @return A \code{\link{GRanges-class}} object.
#' @author David Porubsky
#' @export
#' 
reportGaps <- function(ranges, id.col=NULL) {
  ## Helper function definition
  getGaps <- function(gr, id.col=1) {
    gr <- GenomeInfoDb::keepSeqlevels(gr, value = as.character(unique(GenomeInfoDb::seqnames(gr))), pruning.mode = 'coarse')
    strand(gr) <- '*'
    gr <- GenomicRanges::sort(gr)
    ## Make sure rows are not named
    names(gr) <- NULL
    ## Keep only ranges with the same seqnames
    #max.seqname <- seqlevels(gr)[which.max(runLength(seqnames(gr)))]
    #gr <- gr[seqnames(gr) == max.seqname]
    ## TODO for alignment landing on different chromosomes report gap == 0 and both alignments
    if (length(gr) > 1) {
      ## Calculate gaps
      gap.gr <- GenomicRanges::gaps(gr, start = min(start(gr)))
      gap.gr <- gap.gr[strand(gap.gr) == '*']
      ## Add ID column from the original gr object if defined
      if (!is.null(id.col)) {
        if (id.col > 0 & ncol(mcols(gr)) >= id.col) {
          mcols(gap.gr) <- rep(unique(mcols(gr)[id.col]), length(gap.gr))
        } else {
          warning("User defined 'id.col' number is larger the then total number of columns in input 'gr', skipping adding id column ...")
        }
      }  
      if (length(gap.gr) > 0) {
        ## Report upstream and downstream ranges
        up.idx <- follow(gap.gr, gr)
        toKeep <- !is.na(up.idx)
        gap.gr$up.gr <- GRanges(seqnames=seqnames(gap.gr), ranges=IRanges(start=start(gap.gr), end=start(gap.gr)))
        gap.gr$up.gr[toKeep] <- gr[up.idx[toKeep]][,0]
        down.idx <- precede(gap.gr, gr)
        toKeep <- !is.na(down.idx)
        gap.gr$down.gr <- GRanges(seqnames=seqnames(gap.gr), ranges=IRanges(start=end(gap.gr), end=end(gap.gr)))
        gap.gr$down.gr[toKeep] <- gr[down.idx[toKeep]][,0]
      }  
    }
    ## Export final gap ranges
    return(gap.gr)
  }  
  
  if (class(ranges) == "CompressedGRangesList") {
    ## Keep only contigs with split alignments
    ranges <- ranges[lengths(ranges) > 1]
    gaps <-  suppressWarnings( S4Vectors::endoapply(ranges, function(gr) getGaps(gr=gr, id.col=id.col)) )
    gaps <- unlist(gaps, use.names = FALSE)
  } else if (class(ranges) == "GRanges") {
    gaps <- getGaps(gr=ranges, id.col=id.col)
  } else {
    stop("Only objescts of class 'GRanges' or 'GRangesList' are allowed as input for parameter 'ranges' !!!")
  }
  return(gaps)
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
      
