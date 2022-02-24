#' Make a horizontal self-alignment dot plot.
#' 
#' This function takes self-alignment coordinates from 'nucmer' or 'minimap2' aligner and 
#' construct a horizontal dot plot.
#'
#'
#' @param aln.coords A path to a file containing self-alignment coordinates.
#' @param format Define format of the input alignment coordinates as either 'nucmer' or minimap2 'mm2', [default: 'nucmer']
#' @param min.align.dist Keep alignment pairs with this or larger distance from each other.
#' @param collapse.overlaps Set to \code{TRUE} to merge overlapping pair of alignments with the same relative orientation.
#' @param highlight.pos A single or a set of positions to be highlighted as vertical solid lines.
#' @param highlight.region A pair of positions to be highlighted as vertical range.
#' @param title A title to be added to the final dotplot.
#' @param return Define desired output, either alignment plot as 'plot' or self-alignments in stored \code{\link{GRanges-class}} object as 'selfaln', [default: 'plot']
#' @return A \code{ggplot} object.
#' @importFrom scales comma
#' @author David Porubsky
#' @export
selfdotplot <- function(aln.coords=NULL, format='nucmer', min.align.len=1000, min.align.dist=1000, collapse.overlaps=TRUE, highlight.pos=NULL, highlight.region=NULL, title=NULL, return='plot') {
  ## Check if the submitted file exists
  if (!nchar(aln.coords) > 0 | !file.exists(aln.coords)) {
    stop("Submitted file in 'aln.coords' does not exists !!!")
  }
  
  ## Data input ##
  ################
  if (format == 'nucmer') {
    ## Read in coordinates from nucmer output
    coords <- utils::read.table(aln.coords, skip=5, stringsAsFactors = FALSE)
    coords.df <- data.frame(s1.start=coords$V1,
                            s1.end=coords$V2,
                            s2.start=coords$V4,
                            s2.end=coords$V5,
                            s1.width=abs(coords$V1 - coords$V2),
                            s2.width=abs(coords$V4 - coords$V5),
                            s1.id=coords$V12,
                            s2.id=coords$V12, 
                            stringsAsFactors = FALSE)
  } else if (format == 'mm2') {
    # ## Read in coordinates from minimap2 output
    # coords <- utils::read.table(aln.coords, stringsAsFactors = FALSE, comment.char = '&', fill = TRUE)
    # coords.df <- data.frame(s1.start=coords$V3,
    #                         s1.end=coords$V4,
    #                         s2.start=coords$V8,
    #                         s2.end=coords$V9,
    #                         s1.width=coords$V11,
    #                         s2.width=coords$V11,
    #                         s1.id=coords$V1,
    #                         s2.id=coords$V2, 
    #                         stringsAsFactors = FALSE)
    paf.data <- readPaf(paf.file = aln.coords, include.paf.tags = FALSE)
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
  ## Get distance between alignments
  coords.df <- transform(coords.df, dist = abs(pmin(s2.start, s2.end) - pmax(s1.start, s1.end)))
  
  ## Get max position (for x-axis plotting)
  max.pos <- max(c(coords.df$s1.start, coords.df$s1.end, coords.df$s2.start, coords.df$s2.end))
  
  ## Data filtering ##
  ####################
  ## Remove self-alignments
  ## Remove alignments with the same start and end position
  coords.df <- coords.df[!(coords.df$s1.start == coords.df$s2.start & coords.df$s1.end == coords.df$s2.end),]
  ## Filter by alignment length
  if (min.align.len > 0) {
    coords.df <- coords.df[coords.df$s1.width >= min.align.len & coords.df$s2.width >= min.align.len,]
  }
  ## Filter by alignment distance
  if (min.align.dist > 0) {
    coords.df <- coords.df[coords.df$dist >= min.align.dist,]
  }
  ## Remove duplicated ranges
  #paste0(coords.df$s1.start, coords.df$s1.end) %in% paste0(coords.df$s2.start, coords.df$s2.end)
  xmin <- apply(coords.df[,c('s1.start', 's1.end', 's2.start', 's2.end')], 1, min)
  xmax <- apply(coords.df[,c('s1.start', 's1.end', 's2.start', 's2.end')], 1, max)
  coords.df <- coords.df[!duplicated(paste(xmin, xmax, sep = '_')),]
  
  ## Data transformations ##
  ##########################
  ## Add alignment directionality
  coords.df$dir <- 'rev'
  forw.mask <- (coords.df$s1.start < coords.df$s1.end) & (coords.df$s2.start < coords.df$s2.end)
  coords.df$dir[forw.mask] <- 'forw'
  
  ## Make sure s1 coords are always smaller than s2 coords
  toFlip <- which(pmin(coords.df$s1.start, coords.df$s1.end) > pmax(coords.df$s2.start, coords.df$s2.end))
  coords.df[toFlip,] <- transform(coords.df[toFlip,], 's1.start' = s2.start, 's1.end' = s2.end, 's2.start' = s1.start, 's2.end' = s1.end)
  
  ## Collapse overlapping alignments with the same directionality ##
  ##################################################################
  ## Convert alignments to Genomic ranges
  s1.gr <- GenomicRanges::GRanges(seqnames = 's1', ranges = IRanges::IRanges(start=pmin(coords.df$s1.start,coords.df$s1.end), end=pmax(coords.df$s1.start,coords.df$s1.end)), dir=coords.df$dir)
  s2.gr <- GenomicRanges::GRanges(seqnames = 's2', ranges = IRanges::IRanges(start=pmin(coords.df$s2.start,coords.df$s2.end), end=pmax(coords.df$s2.start,coords.df$s2.end)), dir=coords.df$dir)
  ## Order by S1 coords
  ord <- GenomicRanges::order(s1.gr)
  s1.gr <- s1.gr[ord]
  s2.gr <- s2.gr[ord]
  
  if (collapse.overlaps) {
    ## Get self-alignments
    hits1 <- IRanges::findOverlaps(s1.gr, drop.self=TRUE)
    hits2 <- IRanges::findOverlaps(s2.gr, drop.self=TRUE)
    if (length(hits1) > 1 & length(hits2) > 1) {
      ## Keep unique overlap pairs
      mask1 <- !duplicated(paste0(pmin(queryHits(hits1), subjectHits(hits1)), pmax(queryHits(hits1), subjectHits(hits1))))
      hits1 <- hits1[mask1]
      mask2 <- !duplicated(paste0(pmin(queryHits(hits2), subjectHits(hits2)), pmax(queryHits(hits2), subjectHits(hits2))))
      hits2 <- hits2[mask2]
      ## Keep pairs overlapping each other
      mask1 <- paste0(queryHits(hits1), subjectHits(hits1)) %in% paste0(queryHits(hits2), subjectHits(hits2))
      mask2 <- paste0(queryHits(hits2), subjectHits(hits2)) %in% paste0(queryHits(hits1), subjectHits(hits1))
      hits1 <- hits1[mask1]
      hits2 <- hits2[mask2]
      ## Keep pairs having the same dir
      mask <- s1.gr$dir[queryHits(hits1)] == s1.gr$dir[subjectHits(hits1)] & s2.gr$dir[queryHits(hits2)] == s2.gr$dir[subjectHits(hits2)]
      hits1 <- hits1[mask]
      hits2 <- hits2[mask]
      
      ## Get groups alignments overlapping each other
      groups <- list()
      #for (i in seq_along(hits1)) {
      for (i in order(subjectHits(hits1))) {
        pair <- hits1[i]
        pair <- c(queryHits(pair), subjectHits(pair))
        if (length(groups) == 0) {
          groups[[length(groups) + 1]] <- pair
          group.id <- 1
        } else {
          if (any(pair %in% groups[[group.id]])) {
            groups[[group.id]] <- unique(c(groups[[group.id]], pair))
          } else {
            group.id <- group.id + 1
            groups[[group.id]] <- pair
          }
        }
      }
      ## Get alignment groups
      grp <- unlist(groups)
      names(grp) <- rep(1:length(groups), lengths(groups))
      ## Assign alignment groups to GRanges
      s1.gr$group <- 0
      s2.gr$group <- 0
      s1.gr$group[grp] <- names(grp)
      s2.gr$group[grp] <- names(grp)
      ## Collapse ranges from the same group
      ## alignment1
      s1.collapse.grl <- GenomicRanges::split(s1.gr[s1.gr$group > 0], s1.gr$group[s1.gr$group > 0])
      s1.collapse.gr <- unlist(S4Vectors::endoapply(s1.collapse.grl, range))
      s1.collapse.gr$dir <- sapply(s1.collapse.grl, function(x) unique(x$dir))
      s1.collapse.gr$group <- unique(s1.gr$group[s1.gr$group > 0])
      names(s1.collapse.gr) <- NULL
      ## alignment2
      s2.collapse.grl <- split(s2.gr[s1.gr$group > 0], s2.gr$group[s2.gr$group > 0])
      s2.collapse.gr <- unlist(endoapply(s2.collapse.grl, range))
      s2.collapse.gr$dir <- sapply(s2.collapse.grl, function(x) unique(x$dir))
      s2.collapse.gr$group <- unique(s2.gr$group[s2.gr$group > 0])
      names(s2.collapse.gr) <- NULL
      ## Replace collapsed ranges
      s1.gr <- c(s1.gr[s1.gr$group == 0], s1.collapse.gr)
      s2.gr <- c(s2.gr[s2.gr$group == 0], s2.collapse.gr) 
    }
  }  
  ## Prepare object of self-alignments for export
  self.gr <- s1.gr[,0]
  self.gr$s2 <- s2.gr[,0]
  strand(self.gr) <- '+'
  strand(self.gr$s2) <- ifelse(s1.gr$dir == 'forw', '+', '-')
  
  ## Prepare data for plotting ##
  ###############################
  coords.df <- data.frame(s1.start = start(s1.gr), 
                          s1.end = end(s1.gr),
                          s2.start = start(s2.gr),
                          s2.end = end(s2.gr),
                          s1.width = width(s1.gr),
                          s2.width = width(s2.gr),
                          dir = s1.gr$dir)
  ## Flip start and end for reverse oriented alignments
  coords.df[coords.df$dir == 'rev',] <- transform(coords.df[coords.df$dir == 'rev',], 's2.start' = s2.end, 's2.end' = s2.start)
  ## Order alignments by size
  coords.df$aln.size <- coords.df$s1.width + coords.df$s2.width
  coords.df <- coords.df[order(coords.df$aln.size, decreasing = TRUE),]
  
  plt.df <- coords.df
  ## Get y-axis coordinates
  # s1.sums <- cumsum(plt.df$s1.width)
  # s2.sums <- cumsum(plt.df$s2.width)
  # plt.df$y1 <- c(1, s1.sums[-length(s1.sums)])
  # plt.df$y1end <- s1.sums
  # plt.df$y2 <- c(1, s2.sums[-length(s2.sums)])
  # plt.df$y2end <- s2.sums
  plt.df$y1 <- 0
  plt.df$y1end <- 0
  plt.df$y2 <- 0
  plt.df$y2end <- 0
  for (i in 1:nrow(plt.df)) {
    row.df <- plt.df[i,]
    if (i == 1) {
      row.df$y1 <- 1
      row.df$y2 <- 1
      row.df$y1end <- row.df$s1.width
      row.df$y2end <- row.df$s2.width
      offset <- max(row.df$y1end, row.df$y2end)
    } else {
      row.df$y1 <- offset
      row.df$y2 <- offset
      row.df$y1end <- offset + row.df$s1.width
      row.df$y2end <- offset + row.df$s2.width
      offset <- max(row.df$y1end, row.df$y2end)
    }
    plt.df[i,] <- row.df
  }
  
  ## Get polygon coordinates
  plt.dir.df <- plt.df[plt.df$dir == 'forw',]
  plt.rev.df <- plt.df[plt.df$dir == 'rev',]
  if (nrow(plt.dir.df) > 0) {
    poly.dir.df <- data.frame(x=c(rbind(plt.dir.df$s1.start, plt.dir.df$s2.start, plt.dir.df$s2.end, plt.dir.df$s1.end)),
                              #y=c(rbind(plt.dir.df$y1, plt.dir.df$y2, plt.dir.df$y1end, plt.dir.df$y2end)),
                              y=c(rbind(plt.dir.df$y1, plt.dir.df$y2, plt.dir.df$y2end, plt.dir.df$y1end)),
                              group=rep(1:nrow(plt.dir.df), each=4),
                              direction=rep(plt.dir.df$dir, each=4))
  } else {
    poly.dir.df <- data.frame(x=c(rbind(NaN, NaN, NaN, NaN)),
                              y=c(rbind(NaN, NaN, NaN, NaN)),
                              group=rep(1, each=4),
                              direction=rep('forw', each=4))
  }  
  
  if (nrow(plt.rev.df) > 0) {
    poly.rev.df <- data.frame(x=c(rbind(plt.rev.df$s1.start, plt.rev.df$s1.end, plt.rev.df$s2.end, plt.rev.df$s2.start)),
                              y=c(rbind(plt.rev.df$y1, plt.rev.df$y1end, plt.rev.df$y2end, plt.rev.df$y2)),
                              group=rep(1:nrow(plt.rev.df), each=4),
                              direction=rep(plt.rev.df$dir, each=4))
  } else {
    poly.rev.df <- data.frame(x=c(rbind(NaN, NaN, NaN, NaN)),
                              y=c(rbind(NaN, NaN, NaN, NaN)),
                              group=rep(1, each=4),
                              direction=rep('rev', each=4))
  }  
  
  ## Make dotplot ##
  ##################
  y.limit <- max(c(plt.df$y1end, plt.df$y2end))
  ## Plot alignment pairs
  plt <- ggplot(plt.df) +
    geom_segment(aes(x=s1.start, xend=s1.end, y=y1, yend=y1end)) +
    geom_segment(aes(x=s2.start, xend=s2.end, y=y2, yend=y2end)) +
    geom_polygon(data=poly.dir.df, aes(x=x, y=y, group=group, fill=direction), alpha=0.25, inherit.aes=FALSE) +
    geom_polygon(data=poly.rev.df, aes(x=x, y=y, group=group, fill=direction), alpha=0.25, inherit.aes=FALSE) +
    scale_x_continuous(limits = c(0, max.pos), labels = scales::comma, expand = c(0, 0)) +
    scale_y_continuous(limits = c(-1, y.limit), expand = c(0.1, 0.1)) +
    ylab('Self-alignments') +
    xlab('Contig position (bp)') +
    scale_fill_manual(values = c('forw'='chartreuse4', 'rev'='darkgoldenrod2')) +
    coord_fixed(ratio = 1) +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  
  ## Add alignment dotted lines
  plt <- plt + 
    geom_linerange(aes(x=s1.start, ymin=0, ymax=y1), linetype='dotted') +
    geom_linerange(aes(x=s1.end, ymin=0, ymax=y1end), linetype='dotted') +
    geom_linerange(aes(x=s2.start, ymin=0, ymax=y2), linetype='dotted') +
    geom_linerange(aes(x=s2.end, ymin=0, ymax=y2end), linetype='dotted')
  
  ## Add alignment arrows
  arrow.df <- data.frame('xmin' = c(rbind(coords.df$s1.start, coords.df$s2.start)),
                         'xmax' = c(rbind(coords.df$s1.end, coords.df$s2.end)),
                         'dir' = c(rbind('forw', coords.df$dir))) # to make sure s1 is always forward
  arrow.df$direction <- ifelse(arrow.df$dir == 'forw', 1, -1)         
  ## Make sure start is always smaller than end of the alignment
  arrow.df[,c('xmin', 'xmax')] <- t(apply(arrow.df[,c('xmin', 'xmax')], 1, sort))
  
  plt <- plt + geom_gene_arrow(data=arrow.df, aes(xmin=xmin, xmax=xmax, y=0, forward=direction, fill=dir))
  
  ## Highlight user defined positions 
  if (!is.null(highlight.pos) & is.numeric(highlight.pos)) {
    highlight.pos <- highlight.pos[highlight.pos > 0 & highlight.pos <= max.pos]

    if (length(highlight.pos) > 0) {
      plt <- plt + geom_vline(xintercept = highlight.pos)
    }  
  }
  
  ## Highlight user defined region
  if (!is.null(highlight.region)) {
    highlight.region <- highlight.region[highlight.region$xmin > 0 & highlight.region$xmax <= max.pos,]
    
    if (nrow(highlight.region) > 0) {
      plt <- plt + geom_rect(data = highlight.region, aes(xmin=xmin, xmax=xmax, ymin=0, ymax=Inf), color='red', alpha=0.25)
    }  
  }
  
  ## Add title if defined
  if (!is.null(title)) {
    if (nchar(title) > 0) {
      plt <- plt + ggtitle(title)
    }
  }
  
  if (return == 'plot') {
    ## Return final plot
    return(plt)
  } else if (return == 'selfaln') {
    ## Return final self-alignments
    return(self.gr)    
  } else {
    warning("Please define as 'plot' or 'selfaln', returing plot by default !!!")
    return(plt)
  }
}


#' Plot simple dotplot of two sequences.
#' 
#' This function takes alignment coordinates from 'nucmer' or 'minimap2' aligner and 
#' construct a simple dot plot.
#'
#' @param shape A shape used to plot aligned sequences: Either 'segm' or 'point'.
#' @param genome.coords Set to \code{TRUE} if the target sequence should be reported in genomic coordinates.
#' @inheritParams selfdotplot
#' @return A \code{ggplot} object.
#' @importFrom scales comma
#' @author David Porubsky
#' @export
#' 
simpledotplot <- function(aln.coords=NULL, format='nucmer', min.align.len=1000, highlight.pos=NULL, highlight.region=NULL, shape='segm', title=NULL, genome.coord=FALSE) {
  
  ## Helper function
  remapCoord <- function(x = NULL, new.range = NULL) {
    offset <- (min(new.range) + x[1]) - 1
    dist <- cumsum(diff(x))
    new.x <- c(offset, (offset + dist))
    return(new.x)
  } 
  
  ## Check if the submitted file exists
  if (!nchar(aln.coords) > 0 | !file.exists(aln.coords)) {
    stop("Submitted file in 'aln.coords' does not exists !!!")
  }
  
  ## Data input ##
  ################
  if (format == 'nucmer') {
    ## Read in coordinates from nucmer output
    coords <- utils::read.table(aln.coords, skip=5, stringsAsFactors = FALSE)
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
    paf.data <- readPaf(paf.file = aln.coords, include.paf.tags = FALSE)
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
  
  
  ## Translate sequence coordinates to genome-wide coordinates
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
  if (format == 'mm2') {
    shape <- 'segm'
    message("     Coercing shape == 'segm' for alignments in 'mm2' format!")
  }
  
  if (shape == 'segm') {
    ## Plot alignments
    plt <- ggplot2::ggplot(coords.df, aes(x=s1.start,xend=s1.end,y=s2.start,yend=s2.end, color=dir)) + 
      geom_segment()
  } else if (shape == 'point') {
    coords.df$s1.midpoint <- coords.df$s1.start + ((coords.df$s1.end - coords.df$s1.start)/2)
    coords.df$s2.midpoint <- coords.df$s2.start + ((coords.df$s2.end - coords.df$s2.start)/2)
    plt <- ggplot2::ggplot(coords.df, aes(x=s1.midpoint, y=s2.midpoint, color=dir)) + 
      geom_point()
  }  
  plt <- plt +  
    xlab(unique(coords.df$s1.id)) +
    ylab(unique(coords.df$s2.id)) +
    theme_bw() + 
    theme(aspect.ratio=1) + #force x and y axis to have the same proportion
    scale_color_manual(values = c('chartreuse4', 'darkgoldenrod2')) +
    scale_x_continuous(labels = scales::comma) +
    scale_y_continuous(labels = scales::comma)
  
  ## Highlight user defined positions 
  if (!is.null(highlight.pos) & is.numeric(highlight.pos)) {
    highlight.pos <- highlight.pos[highlight.pos > 0 & highlight.pos <= max.pos]
    
    if (length(highlight.pos) > 0) {
      plt <- plt + geom_vline(xintercept = highlight.pos)
    }  
  }
  
  ## Highlight user defined region
  if (!is.null(highlight.region)) {
    highlight.region <- highlight.region[highlight.region$xmin > 0 & highlight.region$xmax <= max.pos,]
    
    if (nrow(highlight.region) > 0) {
      plt <- plt + geom_rect(data = highlight.region, aes(xmin=xmin, xmax=xmax, ymin=0, ymax=Inf), color='red', alpha=0.25)
    }  
  }
  
  ## Add title if defined
  if (!is.null(title)) {
    if (nchar(title) > 0) {
      plt <- plt + ggtitle(title)
    }
  }
  ## Return final plot
  return(plt)
}
