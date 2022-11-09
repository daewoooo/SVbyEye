#' Make a horizontal self-alignment dot plot.
#' 
#' This function takes self-alignment coordinates from 'nucmer' or 'minimap2' aligner and 
#' construct a horizontal dot plot.
#'
#' @param paf.table A \code{data.frame} or \code{tibble} containing a single or multiple PAF record(s) with 12 mandatory columns.
#' @param shape A shape used to plot aligned sequences: Either 'segment' or 'arc'.
#' @param sort.by Order PAF alignments by relative left-most coordinates ('position') or by the alignment length ('length').
#' @param color.by Color PAF alignments by relative orientation ('direction') or by their sequence identity ('identity').
#' @param highlight.pos A single or a set of positions to be highlighted as vertical solid lines.
#' @param highlight.region A pair of positions to be highlighted as vertical range.
#' @param title A title to be added to the final dotplot.
#' @return A \code{ggplot} object.
#' @import ggplot2
#' @importFrom scales comma
#' @importFrom wesanderson wes_palette
#' @importFrom dplyr group_by mutate arrange
#' @author David Porubsky
#' @export

selfdotplot <- function(paf.table=NULL, shape='segment', sort.by='position', color.by='direction', highlight.pos=NULL, highlight.region=NULL, title=NULL) {
  if (nrow(paf.table$M) > 0) {
    ## Get self-alignments
    paf.data <- paf.table$M
    ## Define alignment id if not defined
    if (!'aln.id' %in% colnames(paf.data)) {
      paf.data$aln.id <- 1:nrow(paf.data)
    }
    ## Get relative orientation
    direction <- ifelse(paf.data$strand == '+', 'forw', 'rev')
    ## Prepare data for plotting
    coords.df <- data.frame(s1.start=paf.data$q.start,
                            s1.end=paf.data$q.end,
                            s2.start=paf.data$t.start,
                            s2.end=paf.data$t.end,
                            s1.width=(paf.data$q.end - paf.data$q.start) + 1,
                            s2.width=(paf.data$t.end - paf.data$t.start) + 1,
                            s1.id=paf.data$q.name,
                            s2.id=paf.data$t.name,
                            aln.len=paf.data$aln.len,
                            dir=direction,
                            identity=(paf.data$n.match / paf.data$aln.len) * 100,
                            aln.id = paf.data$aln.id,
                            stringsAsFactors = FALSE)
    
    ## Get max position (for x-axis plotting)
    max.pos <- unique(paf.data$q.len)
    ## Flip start and end for reverse oriented alignments
    coords.df[coords.df$dir == 'rev',] <- transform(coords.df[coords.df$dir == 'rev',], 's2.start' = s2.end, 's2.end' = s2.start)
    if (sort.by == 'position') {
      ## Order alignments by query position
      ord.aln.id <- coords.df %>% 
        dplyr::group_by(aln.id) %>% 
        dplyr::summarise(start.id = min(s1.start)) %>%
        arrange(start.id)
      coords.df$aln.id <- factor(coords.df$aln.id, levels = ord.aln.id$aln.id)  
      
      coords.df <- coords.df %>% 
        dplyr::group_by(aln.id) %>%
        dplyr::arrange(s1.start, .by_group = TRUE)
      
    } else if (sort.by == 'length') {
      ## Order alignments by alignment length
      coords.df <- coords.df %>% 
        dplyr::group_by(aln.id) %>% 
        dplyr::mutate(length.id = sum(aln.len)) %>% 
        dplyr::arrange(desc(length.id), s1.start)
    } else {
      ## Order by alignment id
      coords.df <- coords.df %>% dplyr::arrange(aln.id, s1.start)
    }
    
    ## Make dotplot ##
    ##################
    plt.df <- coords.df
    ## Color by
    if (color.by == 'direction') {
      col.by <- 'direction'
      colors <- c('forw'='chartreuse4', 'rev'='darkgoldenrod2')
    } else if (color.by == 'identity') {
      ## TODO use a 'getColorScheme' function from plotMiro wrapper (see helpers.R) !!!
      col.by <- 'identity'
      identity.breaks <- c(90, 95, 96, 97, 98, 99, 99.5, 99.9)
      identity.levels <- c('<90', '90:95', '95:96', '96:97', '97:98', '98:99', '99:99.5', '99.5:99.9', '>99.9')
      ids <- findInterval(plt.df$identity, vec = identity.breaks) + 1
      plt.df$identity <- identity.levels[ids]
      colors <- wesanderson::wes_palette(name = "Zissou1", n = length(identity.levels), type = 'continuous')
      colors <- setNames(as.list(colors), identity.levels)
    } else {
      col.by <- 'direction'
    }  
    
    if (shape == 'segment') {
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
                                  direction=rep(plt.dir.df$dir, each=4),
                                  identity=rep(plt.dir.df$identity, each=4))
      } else {
        poly.dir.df <- data.frame(x=c(rbind(NaN, NaN, NaN, NaN)),
                                  y=c(rbind(NaN, NaN, NaN, NaN)),
                                  group=rep(1, each=4),
                                  direction=rep('forw', each=4),
                                  identity=rep('<80', each=4))
      }  
      
      if (nrow(plt.rev.df) > 0) {
        poly.rev.df <- data.frame(x=c(rbind(plt.rev.df$s1.start, plt.rev.df$s1.end, plt.rev.df$s2.end, plt.rev.df$s2.start)),
                                  y=c(rbind(plt.rev.df$y1, plt.rev.df$y1end, plt.rev.df$y2end, plt.rev.df$y2)),
                                  group=rep(1:nrow(plt.rev.df), each=4),
                                  direction=rep(plt.rev.df$dir, each=4),
                                  identity=rep(plt.rev.df$identity, each=4))
      } else {
        poly.rev.df <- data.frame(x=c(rbind(NaN, NaN, NaN, NaN)),
                                  y=c(rbind(NaN, NaN, NaN, NaN)),
                                  group=rep(1, each=4),
                                  direction=rep('rev', each=4),
                                  identity=rep('<80', each=4))
      }  
      
      ## Make segment dotplot
      y.limit <- max(c(plt.df$y1end, plt.df$y2end))
      ## Plot alignment pairs
      plt <- ggplot(plt.df) +
        geom_segment(aes(x=s1.start, xend=s1.end, y=y1, yend=y1end)) +
        geom_segment(aes(x=s2.start, xend=s2.end, y=y2, yend=y2end)) +
        geom_polygon(data=poly.dir.df, aes_string(x='x', y='y', group='group', fill=eval(col.by)), alpha=0.5, inherit.aes=FALSE) +
        geom_polygon(data=poly.rev.df, aes_string(x='x', y='y', group='group', fill=eval(col.by)), alpha=0.5, inherit.aes=FALSE) +
        scale_x_continuous(labels = comma, expand = c(0, 0)) +
        scale_y_continuous(limits = c(-1, y.limit), expand = c(0.1, 0.1)) +
        coord_cartesian(xlim = c(0, max.pos)) +
        ylab('Self-alignments') +
        xlab('Contig position (bp)') +
        scale_fill_manual(values = colors, drop=FALSE) +
        #coord_fixed(ratio = 1) +
        theme_minimal() +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
    } else if (shape == 'arc') {
      ## Make Arc plot
      x <- c(rbind(plt.df$s1.start, plt.df$s2.start, plt.df$s1.end, plt.df$s2.end))
      group <- rep(1:nrow(plt.df), each=4)
      seq.id <- c(rbind('s1', 's2', 's1', 's2'))
      direction <- rep(plt.df$dir, each=4)
      identity <- rep(plt.df$identity, each=4)
      
      plt.df <- data.frame(x=x,
                           y=0,
                           group=group,
                           seq.id=seq.id,
                           direction=direction,
                           identity=identity)
      
      plt <- ggplot(plt.df) +
        geom_wide_arc(aes_string(x='x', y='y', group='group', fill=eval(col.by)), color='gray', size=0.1, alpha=0.5) +
        scale_x_continuous(labels = comma, expand = c(0, 0)) +
        coord_cartesian(xlim = c(0, max.pos)) +
        ylab('Self-alignments') +
        xlab('Contig position (bp)') +
        scale_fill_manual(values = colors, drop=FALSE) +
        theme_minimal() +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
    } else {
      warning("Paremeter shape can only take values 'segment' or 'arc' !!!")
      plt <- ggplot() +
        scale_fill_manual(values = c('forw'='chartreuse4', 'rev'='darkgoldenrod2')) +
        ylab('Self-alignments') +
        xlab('Contig position (bp)') +
        theme_minimal() +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
    }
    
    ## Add alignment arrows
    # arrow.df <- data.frame('xmin' = c(rbind(coords.df$s1.start, coords.df$s2.start)),
    #                        'xmax' = c(rbind(coords.df$s1.end, coords.df$s2.end)),
    #                        'dir' = c(rbind('forw', coords.df$dir))) # to make sure s1 is always forward
    # arrow.df$direction <- ifelse(arrow.df$dir == 'forw', 1, -1)         
    # ## Make sure start is always smaller than end of the alignment
    # arrow.df[,c('xmin', 'xmax')] <- t(apply(arrow.df[,c('xmin', 'xmax')], 1, sort))
    # 
    # plt <- plt + geom_gene_arrow(data=arrow.df, aes(xmin=xmin, xmax=xmax, y=0, forward=direction, fill=dir))
    
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
  } else {
    message("No alignments in submitted 'paf.table', returning empty plot!!!")
    plt <- ggplot()
  }  

  ## Return final plot
  return(plt)
}

## DEPRACATED
# selfdotplot <- function(selfaln.gr=NULL, shape='segment', highlight.pos=NULL, highlight.region=NULL, title=NULL) {
#   ## Get self-alignments
#   self.gr <- selfaln.gr$SelfAlnPairs
#   s1.gr <- self.gr[,0]
#   s2.gr <- self.gr$s2[,0]
#   ## Get relative orientation
#   direction <- ifelse(strand(s1.gr) == strand(s2.gr), 'forw', 'rev')
#   ## Prepare data for plotting ##
#   ###############################
#   coords.df <- data.frame(s1.start = start(s1.gr), 
#                           s1.end = end(s1.gr),
#                           s2.start = start(s2.gr),
#                           s2.end = end(s2.gr),
#                           s1.width = width(s1.gr),
#                           s2.width = width(s2.gr),
#                           dir = direction)
#   ## Get max position (for x-axis plotting)
#   max.pos <- max(c(coords.df$s1.start, coords.df$s1.end, coords.df$s2.start, coords.df$s2.end))
#   ## Flip start and end for reverse oriented alignments
#   coords.df[coords.df$dir == 'rev',] <- transform(coords.df[coords.df$dir == 'rev',], 's2.start' = s2.end, 's2.end' = s2.start)
#   ## Order alignments by size
#   coords.df$aln.size <- coords.df$s1.width + coords.df$s2.width
#   coords.df <- coords.df[order(coords.df$aln.size, decreasing = TRUE),]
#   
#   ## Make dotplot ##
#   ##################
#   plt.df <- coords.df
#   
#   if (shape == 'segment') {
#     plt.df$y1 <- 0
#     plt.df$y1end <- 0
#     plt.df$y2 <- 0
#     plt.df$y2end <- 0
#     for (i in 1:nrow(plt.df)) {
#       row.df <- plt.df[i,]
#       if (i == 1) {
#         row.df$y1 <- 1
#         row.df$y2 <- 1
#         row.df$y1end <- row.df$s1.width
#         row.df$y2end <- row.df$s2.width
#         offset <- max(row.df$y1end, row.df$y2end)
#       } else {
#         row.df$y1 <- offset
#         row.df$y2 <- offset
#         row.df$y1end <- offset + row.df$s1.width
#         row.df$y2end <- offset + row.df$s2.width
#         offset <- max(row.df$y1end, row.df$y2end)
#       }
#       plt.df[i,] <- row.df
#     }
#     
#     ## Get polygon coordinates
#     plt.dir.df <- plt.df[plt.df$dir == 'forw',]
#     plt.rev.df <- plt.df[plt.df$dir == 'rev',]
#     if (nrow(plt.dir.df) > 0) {
#       poly.dir.df <- data.frame(x=c(rbind(plt.dir.df$s1.start, plt.dir.df$s2.start, plt.dir.df$s2.end, plt.dir.df$s1.end)),
#                                 #y=c(rbind(plt.dir.df$y1, plt.dir.df$y2, plt.dir.df$y1end, plt.dir.df$y2end)),
#                                 y=c(rbind(plt.dir.df$y1, plt.dir.df$y2, plt.dir.df$y2end, plt.dir.df$y1end)),
#                                 group=rep(1:nrow(plt.dir.df), each=4),
#                                 direction=rep(plt.dir.df$dir, each=4))
#     } else {
#       poly.dir.df <- data.frame(x=c(rbind(NaN, NaN, NaN, NaN)),
#                                 y=c(rbind(NaN, NaN, NaN, NaN)),
#                                 group=rep(1, each=4),
#                                 direction=rep('forw', each=4))
#     }  
#     
#     if (nrow(plt.rev.df) > 0) {
#       poly.rev.df <- data.frame(x=c(rbind(plt.rev.df$s1.start, plt.rev.df$s1.end, plt.rev.df$s2.end, plt.rev.df$s2.start)),
#                                 y=c(rbind(plt.rev.df$y1, plt.rev.df$y1end, plt.rev.df$y2end, plt.rev.df$y2)),
#                                 group=rep(1:nrow(plt.rev.df), each=4),
#                                 direction=rep(plt.rev.df$dir, each=4))
#     } else {
#       poly.rev.df <- data.frame(x=c(rbind(NaN, NaN, NaN, NaN)),
#                                 y=c(rbind(NaN, NaN, NaN, NaN)),
#                                 group=rep(1, each=4),
#                                 direction=rep('rev', each=4))
#     }  
#     
#     ## Make segment dotplot
#     y.limit <- max(c(plt.df$y1end, plt.df$y2end))
#     ## Plot alignment pairs
#     plt <- ggplot(plt.df) +
#       geom_segment(aes(x=s1.start, xend=s1.end, y=y1, yend=y1end)) +
#       geom_segment(aes(x=s2.start, xend=s2.end, y=y2, yend=y2end)) +
#       geom_polygon(data=poly.dir.df, aes(x=x, y=y, group=group, fill=direction), alpha=0.25, inherit.aes=FALSE) +
#       geom_polygon(data=poly.rev.df, aes(x=x, y=y, group=group, fill=direction), alpha=0.25, inherit.aes=FALSE) +
#       scale_x_continuous(labels = scales::comma, expand = c(0, 0)) +
#       scale_y_continuous(limits = c(-1, y.limit), expand = c(0.1, 0.1)) +
#       coord_cartesian(xlim = c(0, max.pos)) +
#       ylab('Self-alignments') +
#       xlab('Contig position (bp)') +
#       scale_fill_manual(values = c('forw'='chartreuse4', 'rev'='darkgoldenrod2')) +
#       #coord_fixed(ratio = 1) +
#       theme_minimal() +
#       theme(axis.text.y = element_blank(),
#             axis.ticks.y = element_blank())
#     
#     ## Add alignment dotted lines
#     # plt <- plt + 
#     #   geom_linerange(aes(x=s1.start, ymin=0, ymax=y1), linetype='dotted') +
#     #   geom_linerange(aes(x=s1.end, ymin=0, ymax=y1end), linetype='dotted') +
#     #   geom_linerange(aes(x=s2.start, ymin=0, ymax=y2), linetype='dotted') +
#     #   geom_linerange(aes(x=s2.end, ymin=0, ymax=y2end), linetype='dotted')
#   } else if (shape == 'arc') {
#     ## Make Arc plot
#     x <- c(rbind(plt.df$s1.start, plt.df$s2.start, plt.df$s1.end, plt.df$s2.end))
#     group <- rep(1:nrow(plt.df), each=4)
#     seq.id <- c(rbind('s1', 's2', 's1', 's2'))
#     direction <- rep(plt.df$dir, each=4)
#     
#     plt.df <- data.frame(x=x,
#                          y=0,
#                          group=group,
#                          seq.id=seq.id,
#                          direction=direction)
#     
#     plt <- ggplot(plt.df) +
#       #geom_wide_arc(aes(x=x, y=y, group=group, fill=direction), color='gray', alpha=0.25) +
#       geom_wide_arc(aes(x=x, y=y, group=group, fill=direction, color=direction), size=0.1, alpha=0.25) +
#       scale_x_continuous(labels = comma, expand = c(0, 0)) +
#       coord_cartesian(xlim = c(0, max.pos)) +
#       ylab('Self-alignments') +
#       xlab('Contig position (bp)') +
#       scale_fill_manual(values = c('forw'='chartreuse4', 'rev'='darkgoldenrod2')) +
#       scale_color_manual(values = c('forw'='chartreuse4', 'rev'='darkgoldenrod2')) +
#       theme_minimal() +
#       theme(axis.text.y = element_blank(),
#             axis.ticks.y = element_blank())
#   } else {
#     warning("Paremeter shape can only take values 'segment' or 'arc' !!!")
#     plt <- ggplot() +
#       scale_fill_manual(values = c('forw'='chartreuse4', 'rev'='darkgoldenrod2')) +
#       ylab('Self-alignments') +
#       xlab('Contig position (bp)') +
#       theme_minimal() +
#       theme(axis.text.y = element_blank(),
#             axis.ticks.y = element_blank())
#   }
#   
#   ## Add alignment arrows
#   arrow.df <- data.frame('xmin' = c(rbind(coords.df$s1.start, coords.df$s2.start)),
#                          'xmax' = c(rbind(coords.df$s1.end, coords.df$s2.end)),
#                          'dir' = c(rbind('forw', coords.df$dir))) # to make sure s1 is always forward
#   arrow.df$direction <- ifelse(arrow.df$dir == 'forw', 1, -1)         
#   ## Make sure start is always smaller than end of the alignment
#   arrow.df[,c('xmin', 'xmax')] <- t(apply(arrow.df[,c('xmin', 'xmax')], 1, sort))
#   
#   plt <- plt + geom_gene_arrow(data=arrow.df, aes(xmin=xmin, xmax=xmax, y=0, forward=direction, fill=dir))
#   
#   ## Highlight user defined positions 
#   if (!is.null(highlight.pos) & is.numeric(highlight.pos)) {
#     highlight.pos <- highlight.pos[highlight.pos > 0 & highlight.pos <= max.pos]
#     
#     if (length(highlight.pos) > 0) {
#       plt <- plt + geom_vline(xintercept = highlight.pos)
#     }  
#   }
#   
#   ## Highlight user defined region
#   if (!is.null(highlight.region)) {
#     highlight.region <- highlight.region[highlight.region$xmin > 0 & highlight.region$xmax <= max.pos,]
#     
#     if (nrow(highlight.region) > 0) {
#       plt <- plt + geom_rect(data = highlight.region, aes(xmin=xmin, xmax=xmax, ymin=0, ymax=Inf), color='red', alpha=0.25)
#     }  
#   }
#   
#   ## Add title if defined
#   if (!is.null(title)) {
#     if (nchar(title) > 0) {
#       plt <- plt + ggtitle(title)
#     }
#   }
#   
#   ## Return final plot
#   return(plt)
# }


#' Plot simple dotplot of two sequences.
#' 
#' This function takes alignment coordinates from 'nucmer' or 'minimap2' aligner and 
#' construct a simple dot plot.
#'
#' @param shape A shape used to plot aligned sequences: Either 'segm' or 'point'.
#' @param genome.coords Set to \code{TRUE} if the target sequence should be reported in genomic coordinates.
#' @param keep.longest.aln Set to \code{TRUE} if a target sequence with a single longest alignment should be kept.
#' @inheritParams selfdotplot
#' @return A \code{ggplot} object.
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
  if (!nchar(aln.coords) > 0 | !file.exists(aln.coords)) {
    stop("Submitted file in 'aln.coords' does not exists !!!")
  }
  
  ## Data input ##
  ################
  if (format == 'nucmer') {
    ## Read in coordinates from nucmer output
    coords <- utils::read.table(aln.coords, skip=5, stringsAsFactors = FALSE, comment.char = '&')
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
    plt <- ggplot2::ggplot(coords.df, aes(x=s1.start,xend=s1.end,y=s2.start,yend=s2.end, color=dir)) + 
      geom_segment()
  } else if (shape == 'point') {
    coords.df$s1.midpoint <- coords.df$s1.start + ((coords.df$s1.end - coords.df$s1.start)/2)
    coords.df$s2.midpoint <- coords.df$s2.start + ((coords.df$s2.end - coords.df$s2.start)/2)
    plt <- ggplot2::ggplot(coords.df, aes(x=s1.midpoint, y=s2.midpoint, color=dir)) + 
      geom_point()
  } else {
    warning("Paremeter shape can only take values 'segment' or 'point' !!!")
    plt <- ggplot()
  }  
  ## If there are multiple alignment to the target sequence use facets to divide the dotplot
  n.seq <- length(unique(coords.df$s2.id))
  if (n.seq > 1) {
    plt <- plt +
      facet_grid(s2.id ~ ., space='free', scale='free') +
      xlab(unique(coords.df$s1.id)) +
      ylab(paste(unique(coords.df$s2.id), collapse = ';')) +
      theme_bw() + 
      scale_color_manual(values = c('chartreuse4', 'darkgoldenrod2')) +
      scale_x_continuous(labels = scales::comma) +
      scale_y_continuous(labels = scales::comma)
  } else {
    plt <- plt +  
      xlab(unique(coords.df$s1.id)) +
      ylab(unique(coords.df$s2.id)) +
      theme_bw() + 
      #theme(aspect.ratio=1) + #force x and y axis to have the same proportion
      scale_color_manual(values = c('chartreuse4', 'darkgoldenrod2')) +
      scale_x_continuous(labels = scales::comma) +
      scale_y_continuous(labels = scales::comma)
  }
  
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
