#' Synchronize orientation of genomic ranges given the desired majority direction.
#'
#' This function takes a set of ranges in \code{\link{GRanges-class}} or \code{\link{GRangesList-class}} object
#' and return the same ranges in same class object such that total length of plus (direct, '+') and minus (minus, '-') .
#' oriented ranges meet the user defined `majority.strand` orientation.
#'
#' @param ranges A \code{\link{GRanges}} object.
#' @param majority.strand A desired majority strand directionality to be reported.
#' @param strand.only If set to \code{TRUE} simple character string of strand directionality ('+', '-') is reported.
#' @importFrom methods is
#' @importFrom GenomicRanges strand width
#' @return A \code{\link{GRanges-class}} or \code{\link{GRangesList-class}} object.
#' @author David Porubsky
#' @export
#' @examples
#' ## Define test genomic ranges to synchronize directionality ##
#' test.gr <- GenomicRanges::GRanges(
#'     seqnames = "test",
#'     ranges = IRanges::IRanges(
#'         start = c(1, 200, 300),
#'         end = c(200, 300, 500)
#'     ), strand = c("-", "+", "-")
#' )
#' syncRangesDir(ranges = test.gr, majority.strand = "+")
#'
syncRangesDir <- function(ranges, majority.strand = "+", strand.only = FALSE) {
    ## Define majority and minority strand
    if (majority.strand == "+") {
        minority.strand <- "-"
    } else if (majority.strand == "-") {
        minority.strand <- "+"
    } else {
        stop("Parameter 'majority.strand' can only be defined as '+' or '-' !!!")
    }
    ## Check if the input object is GRangesList or single GRanges object
    if (methods::is(ranges, "GRangesList")) {
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
    } else if (methods::is(ranges, "GRanges")) {
        gr <- ranges
        ## Flip directionality to make sure majority strand covers the most bases
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
