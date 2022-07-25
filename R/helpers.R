#' Function to add bezier control points for horizontal layout
#'
#' @param data A \code{data.frame} containing x and y coordinates.
#' @param strength The proportion to move the control point along the y-axis towards the other end of the bezier curve.
#' @return A \code{vector} of rescaled coordinate values.
#' @author David Porubsky
#' @export
add.control.points <- function(data=NULL, strength=0.5) {
  start <- data[c(TRUE, FALSE), ]
  end <- data[c(FALSE, TRUE), ]
  y_diff <- (end$y - start$y) * strength
  mid1 <- start
  mid1$y <- mid1$y + y_diff
  mid2 <- end
  mid2$y <- mid2$y - y_diff
  rbind(start, mid1, mid2, end)
}

#' Function to convert between different coordinate scales
#'
#' This function takes as input query coordinates and a query range and convert them into 
#' coordinate scale defined by a target range.
#'
#' @param x A \code{vector} of coordinate values to be rescaled.
#' @param q.range  A \code{vector} containing min and max value for a query range.
#' @param t.range A \code{vector} containing min and max value for a target range.
#' @return A \code{vector} of rescaled coordinate values.
#' @author David Porubsky
#' @export
q2t <- function(x, q.range, t.range) {
  if (is.numeric(x) && is.numeric(q.range) && is.numeric(t.range)) {
    coord.factor <- (t.range[2] - t.range[1]) / (q.range[2] - q.range[1])
    return(t.range[1] + (x - q.range[1]) * coord.factor)
  }
} 

#' Mirror/reflect genomic ranges given the sequence length.
#'
#' This function flips set of ranges in a mirrored fashion.
#'
#' @param gr A \code{\link{GRanges-class}} object with one or more ranges to be mirrored/reflected given the sequence length.
#' @param seqlength  A \code{numeric} value containing the sequence length from where the submitted ranges originate from.
#' @return A \code{\link{GRanges-class}} object with mirrored coordinates.
#' @author David Porubsky
#' @export
mirrorRanges <- function(gr, seqlength=NULL) {
  if (!is.null(seqlength)) {
    gr.len <- seqlength
  } else if (!is.na(seqlengths(gr))) {
    gr.len <- seqlengths(gr)
  } else {
    stop('No seglength provided!!!')
  }
  if (!all(end(gr) <= gr.len)) {
    stop("One or all submitted ranges are outside of defined seqlength!!!")
  }
  starts <- gr.len - end(gr)
  ends <- (gr.len - start(gr)) 
  #starts <- gr.len - cumsum(width(gr))
  #ends <- starts + width(gr)
  new.gr <- GenomicRanges::GRanges(seqnames = seqnames(gr), ranges = IRanges(start=starts, end=ends), strand = strand(gr))
  suppressWarnings( seqlengths(new.gr) <- gr.len )
  return(new.gr)
}