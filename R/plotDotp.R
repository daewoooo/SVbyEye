#' Plot simple dotplot of two sequences.
#' 
#' This function takes alignment coordinates from 'nucmer' or 'minimap2' aligner and 
#' construct a simple dot plot.
#'
#' @param coords.file A path to a PAF or NUCMER specific file containing alignments to be loaded.
#' @param format Is either 'mm2' or 'nucmer' for minimap2 or nucmer specific input 'coords.file', respectively.
#' @param shape A shape used to plot aligned sequences: Either 'segment' or 'point'.
#' @param keep.longest.aln Set to \code{TRUE} if a target sequence with a single longest alignment should be kept.
#' @param genome.coords Set to \code{TRUE} if the target sequence should be reported in genomic coordinates.
#' @inheritParams filterPaf
#' @inheritParams plotSelf
#' @return A \code{ggplot} object.
#' @import ggplot2
#' @importFrom scales comma
#' @author David Porubsky
#' @export
#' 
simpledotplot <- function(aln.coords=NULL, format='nucmer', shape='segment', min.align.len=1000, keep.longest.aln=FALSE, highlight.pos=NULL, highlight.region=NULL, title=NULL, genome.coord=FALSE) {
  
  ## Helper function
  remapCoord <- function(x = NULL, new.range = NULL) {
    offset <- (min(new.range) + x[1]) - 1
    dist <- cumsum(diff(x))
    new.x <- c(offset, (offset + dist))
    return(new.x)
  } 
  
  ## Check if the submitted file exists
  if (!nchar(coords.file) > 0 | !file.exists(coords.file)) {
    stop("Submitted file in 'coords.file' does not exists !!!")
  }
  
  ## Data input ##
  ################
  if (format == 'nucmer') {
    ## Read in coordinates from nucmer output
    coords <- utils::read.table(coords.file, skip=5, stringsAsFactors = FALSE, comment.char = '&')
    coords.df <- data.frame(s1.start=coords$V1,
                            s1.end=coords$V2,
                            s2.start=coords$V4,
                            s2.end=coords$V5,
                            s1.width=abs(coords$V1 - coords$V2),
                            s2.width=abs(coords$V4 - coords$V5),
                            s1.id=coords$V12,
                            s2.id=coords$V13, 
                            stringsAsFactors = FALSE)
  } else if (format == 'mm2') {
    paf.data <- readPaf(paf.file = coords.file, include.paf.tags = FALSE)
    coords.df <- data.frame(s1.start=paf.data$q.start,
                            s1.end=paf.data$q.end,
                            s2.start=paf.data$t.start,
                            s2.end=paf.data$t.end,
                            s1.width=paf.data$aln.len,
                            s2.width=paf.data$aln.len,
                            s1.id=paf.data$q.name,
                            s2.id=paf.data$t.name,
                            stringsAsFactors = FALSE)
    
    ## Flip start and end for reverse oriented alignments
    coords.df[paf.data$strand == '-',] <- transform(coords.df[paf.data$strand == '-',], 's2.start' = s2.end, 's2.end' = s2.start)
  }
  
  ## Get max position (for x-axis plotting)
  max.pos <- max(c(coords.df$s1.start, coords.df$s1.end, coords.df$s2.start, coords.df$s2.end))
  
  ## Translate sequence coordinates to genome-wide coordinates [only works for nucmer]
  if (format == 'mm2') {genome.coord <- FALSE}
  if (genome.coord) {
    s2.region <- unique(coords.df$s2.id)
    s2.region <- strsplit(s2.region, ":")[[1]][2]
    s2.range <- as.numeric( strsplit(s2.region, "-")[[1]] )
    coords.df$s2.start <- remapCoord(x = coords.df$s2.start, new.range = s2.range)
    coords.df$s2.end <- remapCoord(x = coords.df$s2.end, new.range = s2.range)
  }
  
  ## Data filtering ##
  ####################
  ## Filter by alignment length
  if (min.align.len > 0) {
    coords.df <- coords.df[coords.df$s1.width >= min.align.len & coords.df$s2.width >= min.align.len,]
  }
  ## Remove duplicated ranges
  xmin <- apply(coords.df[,c('s1.start', 's1.end', 's2.start', 's2.end')], 1, min)
  xmax <- apply(coords.df[,c('s1.start', 's1.end', 's2.start', 's2.end')], 1, max)
  coords.df <- coords.df[!duplicated(paste(xmin, xmax, sep = '_')),]
  
  ## Keep only the target sequence with the longest alignment
  if (keep.longest.aln) {
    s2.aln.len <- sapply(split(coords.df$s2.width, coords.df$s2.id), sum)
    to.keep <- names(s2.aln.len)[which.max(s2.aln.len)]
    coords.df <- coords.df[coords.df$s2.id == to.keep,]
  }
  
  ## Data transformations ##
  ##########################
  ## Add alignment directionality
  coords.df$dir <- 'rev'
  forw.mask <- (coords.df$s1.start < coords.df$s1.end) & (coords.df$s2.start < coords.df$s2.end)
  coords.df$dir[forw.mask] <- 'forw'
  
  ## Prepare data for plotting ##
  ###############################
  ## Order alignments by size
  coords.df$aln.size <- coords.df$s1.width + coords.df$s2.width
  coords.df <- coords.df[order(coords.df$aln.size, decreasing = TRUE),]
  
  ## Coerce shape segment in minimap alignments are loaded
  if (format == 'mm2' & shape == 'point') {
    shape <- 'segment'
    message("     Coercing shape == 'segm' for alignments in 'mm2' format!")
  }
  
  if (shape == 'segment') {
    ## Plot alignments
    plt <- ggplot2::ggplot(coords.df, ggplot2::aes(x=s1.start,xend=s1.end,y=s2.start,yend=s2.end, color=dir)) + 
      ggplot2::geom_segment()
  } else if (shape == 'point') {
    coords.df$s1.midpoint <- coords.df$s1.start + ((coords.df$s1.end - coords.df$s1.start)/2)
    coords.df$s2.midpoint <- coords.df$s2.start + ((coords.df$s2.end - coords.df$s2.start)/2)
    plt <- ggplot2::ggplot(coords.df, ggplot2::aes(x=s1.midpoint, y=s2.midpoint, color=dir)) + 
      ggplot2::geom_point()
  } else {
    warning("Paremeter shape can only take values 'segment' or 'point' !!!")
    plt <- ggplot2::ggplot()
  }  
  ## If there are multiple alignment to the target sequence use facets to divide the dotplot
  n.seq <- length(unique(coords.df$s2.id))
  if (n.seq > 1) {
    plt <- plt +
      ggplot2::facet_grid(s2.id ~ ., space='free', scales='free') +
      ggplot2::xlab(unique(coords.df$s1.id)) +
      ggplot2::ylab(paste(unique(coords.df$s2.id), collapse = ';')) +
      ggplot2::theme_bw() + 
      ggplot2::scale_color_manual(values = c('chartreuse4', 'darkgoldenrod2')) +
      ggplot2::scale_x_continuous(labels = scales::comma) +
      ggplot2::scale_y_continuous(labels = scales::comma)
  } else {
    plt <- plt +  
      ggplot2::xlab(unique(coords.df$s1.id)) +
      ggplot2::ylab(unique(coords.df$s2.id)) +
      ggplot2::theme_bw() + 
      #theme(aspect.ratio=1) + #force x and y axis to have the same proportion
      ggplot2::scale_color_manual(values = c('chartreuse4', 'darkgoldenrod2')) +
      ggplot2::scale_x_continuous(labels = scales::comma) +
      ggplot2::scale_y_continuous(labels = scales::comma)
  }
  
  ## Highlight user defined positions 
  if (!is.null(highlight.pos) & is.numeric(highlight.pos)) {
    highlight.pos <- highlight.pos[highlight.pos > 0 & highlight.pos <= max.pos]
    
    if (length(highlight.pos) > 0) {
      plt <- plt + ggplot2::geom_vline(xintercept = highlight.pos)
    }  
  }
  
  ## Highlight user defined region
  if (!is.null(highlight.region)) {
    highlight.region <- highlight.region[highlight.region$xmin > 0 & highlight.region$xmax <= max.pos,]
    
    if (nrow(highlight.region) > 0) {
      plt <- plt + 
        ggplot2::geom_rect(data = highlight.region, ggplot2::aes(xmin=xmin, xmax=xmax, ymin=0, ymax=Inf), color='red', alpha=0.25)
    }  
  }
  
  ## Add title if defined
  if (!is.null(title)) {
    if (nchar(title) > 0) {
      plt <- plt + ggplot2::ggtitle(title)
    }
  }
  ## Return final plot
  return(plt)
}
