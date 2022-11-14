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


#' Add color scheme based on certain values and defined breaks.
#'
#' This function takes a \code{data.frame} or \code{tibble} object and given the defined value field split this field into 
#' chunks based on defined breaks. Each chunk has assigned unique color on a gradient scale.
#'
#' @param data.table A \code{data.frame} or \code{tibble} object to be processed.
#' @param value.field Either a column index or a column name present in submitted 'data.table'  
#' @param breaks User defined breaks in defined 'value.field' in order to split these values into chunks.
#' @return A \code{\link{list}} containing original 'data.table' with extra column adding 'col.levels' based on defined 'breaks'.
#' List also contains element 'color' with a gradient color assigned to each 'col.level'.
#' @importFrom wesanderson wes_palette
#' @importFrom dplyr pull
#' @author David Porubsky
#' @export
getColorScheme <- function(data.table=NULL, value.field=NULL, breaks=c(90, 95, 96, 97, 98, 99, 99.5, 99.9)) {
  ## Check user input
  if (is.null(data.table)) {
    stop("No data submitted, please define 'data.table' parameter !!!")
  }
  if (is.numeric(value.field)) {
    if (ncol(data.table) < value.field) {
      stop("Defined value field index is larger than number of columns in submitted 'data.table' !!!")
    }
  } else if (is.character(value.field)) {
    if (!value.field %in% colnames(data.table)) {
      stop("Defined field name does not match any column name in submitted 'data.table' !!!")
    }
  } else {
    stop("No 'value.field' defined, please define 'value.field' either as column index or column name !!!")
  }
  
  ## Define break ranges
  levels <- c(paste0('<', breaks[1]), 
              paste(breaks[-length(breaks)], breaks[-1], sep = ':'), 
              paste0('>', breaks[length(breaks)]))
  ## Get break intervals
  vals <- data.table %>% dplyr::pull(eval(value.field))
  ids <- findInterval(vals, vec = breaks) + 1
  data.table$col.levels <- levels[ids]
  colors <- wesanderson::wes_palette(name = "Zissou1", n = length(levels), type = 'continuous')
  colors <- setNames(as.list(colors), levels)
  ## Return color scheme
  return(list(data=data.table, colors=colors))
}
