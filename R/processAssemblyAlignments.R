#' Load an de novo assembly aligned to the reference in BED file into a \code{\link{GRanges-class}} object.
#' 
#' This function will take a BED file of contig alignments to a reference genome and converts them 
#' into a \code{\link{GRanges-class}} object. This function also aims to collapse single unique contigs
#' that are aligned in multiple pieces after the alignment to the reference.
#'
#' @param bedfile A BED file of contig alignments to a reference genome.
#' @param index A unique identifier to be added as an 'ID' field. 
#' @param min.mapq A minimum mapping quality of alignments reported in submitted BED file.
#' @param min.align A minimum length of an aligned sequence to a reference genome.
#' @param min.ctg.size A minimum length a final contig after gaps are collapsed.
#' @param max.gap A maximum length of a gap within a single contig alignments to be collapsed.
#' @param report.ctg.ends Set to \code{TRUE} if alignment ends of each contig to the reference should reported as well.
#' @return A \code{\link{GRanges-class}} object.
#' @author David Porubsky
#' 
bed2ranges <- function(bedfile=NULL, index=NULL, min.mapq=10, min.align=10000, min.ctg.size=500000, max.gap=100000, report.ctg.ends=FALSE) {
  if (file.exists(bedfile)) {
    message("Loading BED file: ", bedfile)
    bed.df <- utils::read.table(file = bedfile, header = FALSE, stringsAsFactors = FALSE)
    ## Add column names
    colnames(bed.df) <- c('seqnames','start','end','ctg','mapq','strand')[1:ncol(bed.df)]
  } else {
    stop(paste0("BED file ", bedfile, " doesn't exists !!!"))
  }  
  ## Add index if defined
  if (!is.null(index) & is.character(index)) {
    bed.df$ID <- index
  }
  ## Filter by mapping quality
  if (min.mapq > 0) {
    bed.df <- bed.df[bed.df$mapq >= min.mapq,]
  }
  if (nrow(bed.df) == 0) {
    stop("None of the BED alignments reach user defined mapping quality (min.mapq) !!!")
  }
  ## Convert data.frame to GRanges object
  bed.gr <- GenomicRanges::makeGRangesFromDataFrame(bed.df, keep.extra.columns = TRUE)
  ## Ignore strand
  GenomicRanges::strand(bed.gr) <- '*'
  ## Filter out small alignments
  if (min.align > 0) {
    bed.gr <- bed.gr[width(bed.gr) >= min.align]
  }
  ## Split ranges per contig
  bed.grl <- GenomicRanges::split(bed.gr, bed.gr$ctg)
  ## Make sure data frame is sorted
  #bed.gr <- GenomicRanges::sort(bed.gr, ignore.strand=TRUE)
  ## Collapse merge splitted continuous alignments of the same contig
  #bed.gr <- primatR::collapseBins(bed.gr, id.field = 1)
  #bed.grl <- endoapply(bed.grl, function(gr) primatR::collapseBins(gr, id.field = 1))
  if (max.gap > 0) {
    message("Filling gaps of max size: ", max.gap, 'bp')
    bed.grl <- suppressWarnings( S4Vectors::endoapply(bed.grl, function(gr) fillGaps(gr=gr, max.gap=max.gap)) )
    bed.gr <- unlist(bed.grl, use.names = FALSE)
  }  
  ## Filter out small contigs after contig concatenation
  if (min.ctg.size > 0) {
    message("Keeping collapsed alignments of min size: ", min.ctg.size, 'bp')
    bed.gr <- bed.gr[width(bed.gr) >= min.ctg.size]
  }
  ## Report end positions of each contig
  if (report.ctg.ends == TRUE) {
    message("Reporting alignment ends for each contig")
    bed.grl <- GenomicRanges::split(bed.gr, bed.gr$ctg)
    ctg.ends.grl <- S4Vectors::endoapply(bed.grl, range)
    ctg.ends <- sapply(ctg.ends.grl, function(x) paste(x, collapse = '; '))
    bed.gr$ctg.ends <- rep(ctg.ends, lengths(bed.grl))
  }
  return(bed.gr)
}


#' This function will takes in a \code{\link{GRanges-class}} object of alignments of a single contig to a reference
#' and collapses gaps in alignment based user specified maximum allowed gap.
#'
#' @param gr A \code{\link{GRanges-class}} object of regions of a single contigs aligned to a reference.
#' @param max.gap A maximum length of a gap within a single contig alignments to be collapsed.
#' @return A \code{\link{GRanges-class}} object.
#' @author David Porubsky
#' 
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

#' This function will takes in a \code{\link{GRanges-class}} object of alignments of a single contig to a reference
#' and reports all alignment gaps.
#'
#' @param gr A \code{\link{GRanges-class}} object of regions of a single contigs aligned to a reference.
#' @param id.col A column number from the original \code{\link{GRanges-class}} object to be reported as na unique ID.
#' @return A \code{\link{GRanges-class}} object.
#' @author David Porubsky
#' 
reportGaps <- function(gr, id.col=1) {
  gr <- GenomeInfoDb::keepSeqlevels(gr, value = as.character(unique(GenomeInfoDb::seqnames(gr))), pruning.mode = 'coarse')
  gr <- GenomicRanges::sort(gr)
  gap.gr <- GenomicRanges::gaps(gr, start = min(start(gr)))
  gap.gr <- gap.gr[strand(gap.gr) == '*']
  
  if (id.col > 0 & ncol(mcols(gr)) >= id.col) {
    mcols(gap.gr) <- rep(unique(mcols(gr)[id.col]), length(gap.gr))
  } else {
    warning("User defined 'id.col' number is larger the then total number of columns in input 'gr', skipping adding id column ...")
  }
  return(gap.gr)
}
