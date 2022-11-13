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
#' @examples
#'## Get PAF to plot
#'paf.file <- system.file("extdata", "test2.paf", package="SVbyEye")
#'## Read in PAF 
#'paf.table <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = 'cg')
#'## Make a plot
#'## Color by alignment directionality
#'selfdotplot(paf.table = paf.table, color.by = 'direction')
#'
plotSelf <- function(paf.table=NULL, min.deletion.size=NULL, min.insertion.size=NULL, highlight.sv=NULL, binsize=NULL, shape='segment', sort.by='position', color.by='direction', highlight.pos=NULL, highlight.region=NULL, title=NULL) {
  ## Check user input
  ## Make sure submitted paf.table has at least 12 mandatory fields
  if (ncol(paf.table) >= 12) {
    paf <- paf.table
    paf$direction.flip <- FALSE
  } else {
    stop('Submitted PAF alignments do not contain a minimum of 12 mandatory fields, see PAF file format definition !!!')
  }
  
  ## Break PAF at insertion/deletions defined in cigar string
  if (!is.null(min.deletion.size) | !is.null(min.insertion.size)) {
    if (min.deletion.size > 0 | min.insertion.size > 0) {
      paf.l <- breakPaf(paf.table = paf.table, min.deletion.size = min.deletion.size, min.insertion.size = min.insertion.size, collapse.mismatches = TRUE, report.sv = TRUE)
      paf <- paf.l$M
      paf.svs <- paf.l$SVs
    }  
  } else {
    paf$aln.id <- 1:nrow(paf)
    paf.svs <- NULL
  }
  ## Store PAF alignments for later addition of 'geom_gene_arrow' 
  paf.copy <- paf
  paf.copy$ID <- 'M'
  
  ## Bin PAF alignments 
  if (!is.null(binsize)) {
    if (binsize > 0) {
      if (binsize < 10) {
        binsize <- 10
        warning('Minimum allowed bin size is 10, forced binsize=10!!!')
      }  
      paf <- pafToBins(paf.table = paf, binsize = binsize)
      ## If the PAF alignments are binned only color.by = 'fraction.matches' is allowed
      color.by <- 'identity'
    } 
  }
  ## Mark alignments ranges by 'M' (matches)
  paf$ID <- 'M'
  
  ## Add SVs to the alignment table
  if (!is.null(paf.svs)) {
    if (nrow(paf.svs) > 0) {
      paf.svs$ID <- 'INS'
      paf.svs$ID[grep(paf.svs$cg, pattern = 'D', ignore.case = TRUE)] <- 'DEL'
      paf <- dplyr::bind_rows(paf, paf.svs)
    }  
  }  
  
  ## Prepare data for plotting
  paf <- paf[paf$ID == 'M',]
  coords <- data.frame(s1.start=paf$q.start,
                       s1.end=paf$q.end,
                       s2.start=paf$t.start,
                       s2.end=paf$t.end,
                       s1.width=(paf$q.end - paf$q.start) + 1,
                       s2.width=(paf$t.end - paf$t.start) + 1,
                       s1.id=paf$q.name,
                       s2.id=paf$t.name,
                       dir=paf$strand,
                       n.match=paf$n.match,
                       aln.len=paf$aln.len,
                       identity=(paf$n.match / paf$aln.len) * 100,
                       aln.id = paf$aln.id,
                       stringsAsFactors = FALSE)
    
    ## Get max position (for x-axis plotting)
    max.pos <- unique(paf$q.len)
    ## Flip start and end for reverse oriented alignments
    coords[coords$dir == '-',] <- transform(coords[coords$dir == '-',], 's2.start' = s2.end, 's2.end' = s2.start)
    ## Sort alignments
    if (sort.by == 'position') {
      ## Order alignments by query position
      ord.aln.id <- coords %>% 
        dplyr::group_by(aln.id) %>% 
        dplyr::summarise(start.id = min(s1.start)) %>%
        dplyr::arrange(start.id)
      coords$aln.id <- factor(coords$aln.id, levels = ord.aln.id$aln.id)  
      
      coords <- coords %>% 
        dplyr::group_by(aln.id) %>%
        dplyr::arrange(s1.start, .by_group = TRUE)
      
    } else if (sort.by == 'length') {
      ## Order alignments by alignment length
      coords <- coords %>% 
        dplyr::group_by(aln.id) %>% 
        dplyr::mutate(length.id = sum(aln.len)) %>% 
        dplyr::arrange(desc(length.id), s1.start)
    } else {
      ## Order by alignment id
      coords <- coords %>% dplyr::arrange(aln.id, s1.start)
    }
    
    
    ## Make self-dotplot ##
    #######################
    ## Color by
    if (color.by == 'direction') {
      #colors <- c('forw'='chartreuse4', 'rev'='darkgoldenrod2')
      colors <- c('-' = 'cornflowerblue', '+' = 'forestgreen')
    } else if (color.by == 'identity') {
      coords$identity[is.nan(coords$identity) | is.na(coords$identity)] <- 0
      ## Define color scheme
      coords.l <- getColorScheme(data.table = coords, value.field = 'identity', breaks=c(90, 95, 99, 99.5, 99.6, 99.7, 99.8, 99.9))
      coords <- coords.l$data
      colors <- coords.l$colors
    } else {
      col.by <- 'direction'
    }  
    
    if (shape == 'segment') {
      coords$y1 <- 0
      coords$y1end <- 0
      coords$y2 <- 0
      coords$y2end <- 0
      for (i in 1:nrow(coords)) {
        row.df <- coords[i,]
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
        coords[i,] <- row.df
      }
      
      ## Get polygon coordinates
      plt.dir.df <- coords[coords$dir == '+',]
      plt.rev.df <- coords[coords$dir == '-',]
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
      y.limit <- max(c(coords$y1end, coords$y2end))
      ## Plot alignment pairs
      plt <- ggplot(coords) +
        geom_segment(aes(x=s1.start, xend=s1.end, y=y1, yend=y1end)) +
        geom_segment(aes(x=s2.start, xend=s2.end, y=y2, yend=y2end)) +
        geom_polygon(data=poly.dir.df, aes_string(x='x', y='y', group='group', fill=eval(color.by)), alpha=0.5, inherit.aes=FALSE) +
        geom_polygon(data=poly.rev.df, aes_string(x='x', y='y', group='group', fill=eval(color.by)), alpha=0.5, inherit.aes=FALSE) +
        scale_x_continuous(labels = scales::comma, expand = c(0, 0)) +
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
      x <- c(rbind(coords$s1.start, coords$s2.start, coords$s1.end, coords$s2.end))
      group <- rep(1:nrow(coords), each=4)
      seq.id <- c(rbind('s1', 's2', 's1', 's2'))
      direction <- rep(coords$dir, each=4)
      identity <- rep(coords$identity, each=4)
      
      plt.df <- data.frame(x=x,
                           y=0,
                           group=group,
                           seq.id=seq.id,
                           direction=direction,
                           identity=identity)
      
      plt <- ggplot(plt.df) +
        geom_wide_arc(aes_string(x='x', y='y', group='group', fill=eval(color.by)), color='gray', size=0.1, alpha=0.5) +
        scale_x_continuous(labels = scales::comma, expand = c(0, 0)) +
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
        scale_fill_manual(values = c('-' = 'cornflowerblue', '+' = 'forestgreen')) +
        ylab('Self-alignments') +
        xlab('Contig position (bp)') +
        theme_minimal() +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
    }
    
    ## Add indels
    # if (!is.null(highlight.sv)) {
    #   if (nrow(coords[coords$ID != 'M',]) > 0) { 
    #     ## Add SVs to the plot
    #     if (highlight.sv == 'outline') {
    #       plt <- plt + ggnewscale::new_scale_color() +
    #         geom_miropeats(data=coords[coords$ID != 'M',], aes(x, y, group = group, color=ID), fill=NA, alpha=0.5, inherit.aes = FALSE) +
    #         scale_color_manual(values = c('DEL' = 'firebrick3', 'INS' = 'dodgerblue3'), name='SV class')
    #     } else if (highlight.sv == 'fill') {
    #       plt <- plt + ggnewscale::new_scale_fill() +
    #         geom_miropeats(data=coords[coords$ID != 'M',], aes(x, y, group = group, fill=ID), alpha=0.5, inherit.aes = FALSE) +
    #         scale_fill_manual(values = c('DEL' = 'firebrick3', 'INS' = 'dodgerblue3'), name='SV class')
    #     } else {
    #       warning("Parameter 'highlight.sv' can only take values 'outline' or 'fill', see function documentation!!!")
    #     } 
    #   } else {
    #     warning("There are no SVs to highlight. Make sure parameters 'min.deletion.size' and 'min.insertion.size' are set or decrease their values!!!")
    #   }  
    # }
    
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
  
  ## Return final plot
  return(plt)
}