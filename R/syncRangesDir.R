#' Make sure that directionality of set of ranges reflect desired majority strand
#'
#' @param ranges A \code{\link{GRanges}} object.
#' @param majority.strand A desired majority strand directionality to be reported.
#' @param strand.only If set to \code{TRUE} simple character string of strand directionality ('+', '-') is reported.
#' @importFrom GenomicRanges strand width
#' @author David Porubsky
#' @export
#'
syncRangesDir <- function(ranges, majority.strand = '+', strand.only = FALSE) {
  ## Define majority and minority strand
  if (majority.strand == '+') {
    minority.strand = '-'
  } else if (majority.strand == '-') {
    minority.strand = '+'
  } else {
    stop("Parameter 'majority.strand' can only be defined as '+' or '-' !!!")
  }
  ## Check if the input object is GRangesList or single GRanges object
  if (grepl(class(ranges), pattern = 'GRangesList')) {
    for (i in seq_along(ranges)) {
      gr <- ranges[[i]]
      ## Flip ranges directionality to make sure majority strand covers the most bases
      if (sum(GenomicRanges::width(gr[GenomicRanges::strand(gr) == majority.strand])) > sum(GenomicRanges::width(gr[GenomicRanges::strand(gr) == minority.strand]))) {
        ranges[[i]] <- gr
      } else {
        gr.new <- gr
        GenomicRanges::strand(gr.new)[GenomicRanges::strand(gr) == majority.strand] <- minority.strand
        GenomicRanges::strand(gr.new)[GenomicRanges::strand(gr) == minority.strand] <- majority.strand
        ranges[[i]] <- gr.new
      }
    }
    if (strand.only) {
      return(as.character(GenomicRanges::strand(unlist(ranges, use.names = FALSE))))
    } else {
      return(ranges)
    }
  } else if (is(ranges, 'GRanges')) {
    gr <- ranges
    ## Flip directionality based to make sure majority strand covers the most bases
    if (sum(GenomicRanges::width(gr[GenomicRanges::strand(gr) == majority.strand])) > sum(GenomicRanges::width(gr[GenomicRanges::strand(gr) == minority.strand]))) {
      if (strand.only) {
        return(as.character(GenomicRanges::strand(gr)))
      } else {
        return(gr)
      }
    } else {
      gr.new <- gr
      GenomicRanges::strand(gr.new)[GenomicRanges::strand(gr) == majority.strand] <- minority.strand
      GenomicRanges::strand(gr.new)[GenomicRanges::strand(gr) == minority.strand] <- majority.strand
      if (strand.only) {
        return(as.character(GenomicRanges::strand(gr.new)))
      } else {
        return(gr.new)
      }
    }
  } else {
    stop("Only valid 'GRanges' or 'GRangesList' object can be processed !!!")
  }
}
