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
#' @inheritParams paf2coords
#' @return A \code{ggplot} object.
#' @author David Porubsky
#' @export
selfdotplot <- function(aln.coords=NULL, format='nucmer', min.align.len=1000, min.align.dist=1000, collapse.overlaps=TRUE) {
  ## Check if the submitted file exists
  if (!nchar(aln.coords) > 0 | !file.exists(aln.coords)) {
    stop("Submitted file in 'aln.coords' does not exists !!!")
  }
  
  ## Data input ##
  ################
  if (format == 'nucmer') {
    ## Read in coordinates from nucmer output
    coords <- read.table(aln.coords, skip=5, stringsAsFactors = FALSE)
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
    ## Read in coordinates from minimap2 output
  }
  ## Get distance between alignments
  coords.df <- transform(coords.df, dist = abs(pmin(s2.start, s2.end) - pmax(s1.start, s1.end)))
  
  ## Get max position
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
  ## Add alignment directionality
  coords.df$dir <- 'rev'
  forw.mask <- (coords.df$s1.start < coords.df$s1.end) & (coords.df$s2.start < coords.df$s2.end)
  coords.df$dir[forw.mask] <- 'forw'
  
  ## Collapse overlapping alignments with the same directionality ##
  ##################################################################
  ## Convert alignments to Genomic ranges
  s1.gr <- GenomicRanges::GRanges(seqnames = 's1', ranges = IRanges::IRanges(start=pmin(coords.df$s1.start,coords.df$s1.end), end=pmax(coords.df$s1.start,coords.df$s1.end)), dir=coords.df$dir)
  s2.gr <- GenomicRanges::GRanges(seqnames = 's2', ranges = IRanges::IRanges(start=pmin(coords.df$s2.start,coords.df$s2.end), end=pmax(coords.df$s2.start,coords.df$s2.end)), dir=coords.df$dir)
  
  if (collapse.overlaps) {
    ## Get self-alignments
    hits1 <- IRanges::findOverlaps(s1.gr, drop.self=TRUE)
    hits2 <- IRanges::findOverlaps(s2.gr, drop.self=TRUE)
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
                          s2.width =width(s2.gr),
                          dir = s1.gr$dir)
  ## Flip start and end for reverse oriented alignments
  coords.df[coords.df$dir == 'rev',] <- transform(coords.df[coords.df$dir == 'rev',], 's2.start' = s2.end, 's2.end' = s2.start)
  ## Order alignments by size
  coords.df$aln.size <- coords.df$s1.width + coords.df$s2.width
  coords.df <- coords.df[order(coords.df$aln.size, decreasing = TRUE),]
  
  plt.df <- coords.df
  ## Get y-axis coordinates
  s1.sums <- cumsum(plt.df$s1.width)
  s2.sums <- cumsum(plt.df$s2.width)
  plt.df$y1 <- c(1, s1.sums[-length(s1.sums)])
  plt.df$y1end <- s1.sums
  plt.df$y2 <- c(1, s2.sums[-length(s2.sums)])
  plt.df$y2end <- s2.sums
  ## Get polygon coordinates
  plt.dir.df <- plt.df[plt.df$dir == 'forw',]
  plt.rev.df <- plt.df[plt.df$dir == 'rev',]
  poly.dir.df <- data.frame(x=c(rbind(plt.dir.df$s1.start, plt.dir.df$s2.start, plt.dir.df$s2.end, plt.dir.df$s1.end)),
                            y=c(rbind(plt.dir.df$y1, plt.dir.df$y2, plt.dir.df$y1end, plt.dir.df$y2end)),
                            group=rep(1:nrow(plt.dir.df), each=4),
                            direction=rep(plt.dir.df$dir, each=4))
  poly.rev.df <- data.frame(x=c(rbind(plt.rev.df$s1.start, plt.rev.df$s1.end, plt.rev.df$s2.end, plt.rev.df$s2.start)),
                            y=c(rbind(plt.rev.df$y1, plt.rev.df$y1end, plt.rev.df$y2end, plt.rev.df$y2)),
                            group=rep(1:nrow(plt.rev.df), each=4),
                            direction=rep(plt.rev.df$dir, each=4))
  
  ## Make dotplot ##
  ##################
  y.limit <- max(c(plt.df$y1end, plt.df$y2end))
  ## Plot alignment pairs
  plt <- ggplot(plt.df) +
    geom_segment(aes(x=s1.start, xend=s1.end, y=y1, yend=y1end)) +
    geom_segment(aes(x=s2.start, xend=s2.end, y=y2, yend=y2end)) +
    geom_polygon(data=poly.dir.df, aes(x=x, y=y, group=group, fill=direction), alpha=0.25, inherit.aes=FALSE) +
    geom_polygon(data=poly.rev.df, aes(x=x, y=y, group=group, fill=direction), alpha=0.25, inherit.aes=FALSE) +
    scale_x_continuous(limits = c(0, max.pos), labels = comma, expand = c(0, 0)) +
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
  
  ## Return final plot
  return(plt)
  ## Return final self-alignments
  #return(self.gr)
}
