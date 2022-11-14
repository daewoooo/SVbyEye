#' Plot PAF file as miropeat alignments 
#' 
#' This function takes PAF output file from minimap2 alignments, and visualize the alignments
#' in miropeat style. 
#'
#' @param sd.annot A \code{\link[GenomicRanges]{GRanges}} object containing segmental duplication coordinates.
#' @inheritParams paf2coords
#' @return A \code{list} of miropeat style plots.
#' @importFrom scales comma
#' @importFrom wesanderson wes_palette
#' @importFrom gggenes geom_gene_arrow
#' @importFrom ggnewscale new_scale_fill new_scale_color
#' @author David Porubsky
plotMiro_deprecated <- function(paf.file = paf.file, min.mapq = 10, min.align.len = 100, min.align.n = 1, min.deletion.size=NULL, min.insertion.size=NULL, sd.annot = NULL, drop.self.align = FALSE) {
  ## Load PAF file
  #coords.data <- paf2coords(paf.file = paf.file, min.mapq = min.mapq, min.align.len = min.align.len, min.align.n = min.align.n, min.deletion.size=min.deletion.size, min.insertion.size=min.insertion.size, drop.self.align = drop.self.align)  
  ## Process data per alignment
  coords.data.l <- split(coords.data, coords.data$align.id)
  plots <- list()
  for (i in seq_along(coords.data.l)) {
    coords <- coords.data.l[[i]]
    target.seqname <- unique(coords$seq.name[coords$seq.id == 'target'])
    ## Get y-axis labels
    q.range <- range(coords$seq.pos[coords$seq.id == 'query'])
    t.range <- range(coords$seq.pos[coords$seq.id == 'target'])
    ## Adjust target ranges given the size difference with respect query ranges
    range.offset <- diff(q.range) - diff(t.range)
    t.range[2] <- t.range[2] + range.offset ## Make a start position as offset and change only end position
    ## Get x-axis labels
    q.labels <- pretty(q.range)
    t.labels <- pretty(t.range)
    ## Covert query to target coordinates
    q.breaks <- SVbyEye::q2t(x = q.labels, q.range = q.range, t.range = t.range)
    t.breaks <- t.labels
    
    ## Get x-axis labels
    seq.labels <- c(unique(coords$seq.name[coords$seq.id == 'query']), 
                    unique(coords$seq.name[coords$seq.id == 'target']))
    
    coords$frac.match <- coords$n.match / coords$aln.len
    
    plt <- ggplot2::ggplot(coords) +
      geom_miropeats(aes(x, y, group = group, fill=frac.match), alpha=0.5) +
      scale_y_continuous(breaks = c(1, 2), labels = seq.labels) +
      scale_x_continuous(breaks = q.breaks, labels = scales::comma(q.labels),
                         sec.axis = sec_axis(trans = y ~ ., breaks = t.breaks, labels = scales::comma(t.labels)), expand = c(0,0)) +
      scale_fill_gradient(low = 'gray', high = 'red') +
      xlab('Genomic position (bp)') +
      ylab('') +
      theme_bw()
    
    ## Add arrows
    start <- coords$x[c(T, T, F, F)]
    end <- coords$x[c(F, F, T, T)]
    y <- coords$y[c(T, T, F, F)]
    group <- coords$group[c(T, T, F, F)]
    plt.df <- data.frame(start=start, end=end, y=y, group=group)
    plt.df$direction <- ifelse(plt.df$start < plt.df$end, '+', '-')
    
    plt <- plt + ggnewscale::new_scale_fill() + ggnewscale::new_scale_color() +
      gggenes::geom_gene_arrow(data=plt.df, aes(xmin = start, xmax = end, y = y, color= direction, fill = direction), arrowhead_height = unit(3, 'mm')) +
      scale_fill_manual(values = c('cornflowerblue',' forestgreen')) +
      scale_color_manual(values = c('cornflowerblue',' forestgreen')) +
      theme_bw()
    
    ## Add SD arrowheads
    if (!is.null(sd.annot) & class(sd.annot) == 'GRanges') {
      ## Restrict a target region 
      sd.annot <- sd.annot[start(sd.annot) > t.range[1] & end(sd.annot) < t.range[2] & GenomeInfoDb::seqnames(sd.annot) == target.seqname]
      if (length(sd.annot) > 0) {
        sd.annot.df <- as.data.frame(sd.annot)
        if ('fracMatch' %in% colnames(sd.annot.df)) {
          ## Define SD colors
          sd.categ <- findInterval(sd.annot.df$fracMatch, vec = c(0.95, 0.98, 0.99))
          sd.categ <- dplyr::recode(sd.categ, '0' = '<95%', '1' = '95-98%', '2' = '98-99%', '3'='>=99%')
          sd.categ <- factor(sd.categ, levels=c('<95%', '95-98%', '98-99%', '>=99%'))
          sd.annot.df$sd.categ <- sd.categ
          ## Get colors
          pal <- wesanderson::wes_palette("Zissou1", 4, type = "continuous")
        } else {
          sd.annot.df$sd.categ <- 'unknown'
          pal <- 'black'
        }  
        
        plt <- plt + new_scale_fill() + new_scale_color() +
          geom_arrowhead(data=sd.annot.df, aes(xmin=start, xmax=end, y=2.05, color=sd.categ, fill=sd.categ)) +
          scale_fill_manual(values = pal) +
          scale_color_manual(values = pal)
      }  
    }
    ## Save plot
    plots[[i]] <- plt
  }
  return(plots)
}


#' Plot all-versus-all alignments stored in PAf format.
#' 
#' This function takes PAF output file from minimap2 reporting all-versus-all alignments of multiple FASTA sequences 
#' and visualize the alignments in a miropeat style. 
#'
#' @param seqnames.order A user defined order sequence names to be plotted from top to the bottom.
#' @param outline.miro Set to \code{TRUE} if separate miropeat segments should be outlined.
#' @inheritParams readPaf
#' @inheritParams paf2coords
#' @inheritParams syncRangesDir
#' @return A \code{list} of miropeat style plots.
#' @importFrom scales comma
#' @author David Porubsky
#' @export
plotAVA_deprecated <- function(paf.file = NULL, seqnames.order=NULL, seqnames.grep=NULL, target.region=NULL, min.align.n=1, min.mapq=0, min.align.len=1000, drop.self.align=TRUE, majority.strand=NULL, outline.miro=FALSE) {
  ## Check user input
  if (!file.exists(paf.file)) {
    stop("Submitted file ", basename(paf.file), " doesn't exists !!!")
  }
  ## Load PAF file ##
  paf <- readPaf(paf.file = paf.file, include.paf.tags = TRUE)
  
  ## Rename sequences if named vector defined by a user
  if (any(names(seqnames.order) %in% paf$q.name)) {
    paf$q.name <- dplyr::recode(paf$q.name, !!!(seqnames.order))
    paf$t.name <- dplyr::recode(paf$t.name, !!!(seqnames.order))
  }
  ## Desired sequence order
  seq.ids <- unique(paf$q.name)
  if (is.character(seqnames.order)) {
    if (all(seq.ids %in% seqnames.order)) {
      seq.ord <- seqnames.order
    } else {
      seq.ord <- c(seqnames.order, setdiff(seq.ids, seqnames.order))
    }
  } else {
    seq.ord <- NULL
  } 
  
  ## Filter alignments by target region
  if (!is.null(target.region)) {
    if (grepl(target.region, pattern = '\\w+\\d+:\\d+-\\d+')) {
      target.region.gr <- as(target.region, 'GRanges')
    } else if (class(target.region.gr) == 'GRanges') {
      target.region.gr <- target.region
    } else {
      message("Parameter 'target.region' can either be 'GRanges' object or character string 'chr#:start-end'!!!")
    }
    ## Subset PAF by ranges overlaps
    target.gr <- GenomicRanges::makeGRangesFromDataFrame(paf, seqnames.field = 't.name', start.field = 't.start', end.field = 't.end')
    hits <- GenomicRanges::findOverlaps(target.gr, target.region.gr)
    paf <- paf[S4Vectors::queryHits(hits),]
  }
  
  ## Get unique alignment ID
  paf$seq.pair <- paste0(paf$q.name, '___', paf$t.name)
  
  ## Filter PAF alignments ##
  ###########################
  ## Get number of alignments per sequence pair
  if (min.align.n > 0) {
    paf <- paf %>% dplyr::group_by(seq.pair) %>% dplyr::mutate(align.n = dplyr::n())
    paf <- paf[paf$align.n >= min.align.n,]
  }
  ## Keep only specific sequence/region name
  if (!is.null(seqnames.grep)) {
    paf <- paf[grep(paf$q.name, pattern = seqname.grep),]
    paf <- paf[grep(paf$t.name, pattern = seqname.grep),]
  }
  ## Filter by mapping quality
  if (min.mapq > 0 & is.numeric(paf$mapq)) {
    paf <- paf[paf$mapq >= min.mapq,]
  }
  ## Filter by alignment length
  if (min.align.len > 0 & is.numeric(paf$aln.len)) {
    paf <- paf[paf$aln.len >= min.align.len,]
  }
  ## Filter out self-alignments
  if (drop.self.align) {
    paf <- paf[!(paf$q.name == paf$t.name),]
  }
  
  ## Check if there are any alignments left after filtering to plot
  if (nrow(paf) == 0) {
    stop("All alignments were filtered out, Please check if the input 'paf.file' contains correct all-vs-all alignments!")
  }
  
  ## Sync by majority strand directionality
  if (!is.null(majority.strand)) {
    ## Define majority and minority strand
    if (majority.strand == '+') {
      minority.strand = '-'
    } else if (majority.strand == '-') {
      minority.strand = '+'
    } else {
      stop("Parameter 'majority.strand' can only be defined as '+' or '-' !!!")
    }
    
    paf.l <- split(paf, paf$q.name)
    for (i in seq_along(paf.l)) {
      paf.ctg <- paf.l[[i]]
      ## Flip directionality based to make sure majority strand covers the most bases
      if (sum(paf.ctg$aln.len[paf.ctg$strand == majority.strand]) > sum(paf.ctg$aln.len[paf.ctg$strand == minority.strand])) {
        paf.l[[i]] <- paf.ctg
      } else {
        paf.ctg.new <- paf.ctg
        ## Flip alignment orientation
        paf.ctg.new$strand[paf.ctg$strand == majority.strand] <- minority.strand
        paf.ctg.new$strand[paf.ctg$strand == minority.strand] <- majority.strand
        ## Flip query coordinates
        paf.ctg.new$q.end <- paf.ctg$q.len - paf.ctg$q.start
        paf.ctg.new$q.start <- paf.ctg$q.len - paf.ctg$q.end
        
        paf.l[[i]] <- paf.ctg.new
      }
    }
    paf <- do.call(rbind, paf.l)
  }
  
  ## Flip start-end if strand == '-'
  paf[paf$strand == '-', c('t.start','t.end')] <- rev(paf[paf$strand == '-', c('t.start','t.end')])
  
  ## Get number of matching bases in plus an minus orientation
  #paf %>% dplyr::group_by(seq.pair, strand) %>% dplyr::summarise(match = sum(n.match)) %>% dplyr::mutate(diff = match[strand == '+'] - match[strand == '-'])
  
  if (is.null(seq.ord)) {
    ## Order alignments based on number of mismatches
    paf.ord <- paf %>% dplyr::group_by(seq.pair) %>% dplyr::summarise(n.nm = sum(NM)) %>% dplyr::arrange(n.nm)
    paf.ord.pairs <- unlist(strsplit(paf.ord$seq.pair, '___'))
    ## Assign level to seq.names ordered by number of matching bases in plus an minus orientation
    seq.names <- paf.ord.pairs[!duplicated(paf.ord.pairs)]
  } else {
    ## Assign user defined assembly order
    seq.names <- seq.ord
  } 
  
  ## Assign level to seq.names
  seq.ids <- length(seq.names):1
  names(seq.ids) <- seq.names
  paf$y1 <- seq.ids[match(paf$q.name, names(seq.ids))]
  paf$y2 <- seq.ids[match(paf$t.name, names(seq.ids))]
  ## Keep subsequent comparisons only
  paf <- paf[abs(paf$y2 - paf$y1) == 1,]
  
  ## Flip query and target for alignments where query comes first
  flipQT <- which(paf$y1 > paf$y2)
  paf[flipQT,] <- transform(paf[flipQT,], 
                            q.name = t.name, q.start = t.start, q.end = t.end,
                            t.name = q.name, t.start = q.start, t.end = q.end,
                            y1 = y2, y2 = y1)
  paf$seq.pair[flipQT] <- paste0(paf$q.name[flipQT], '___', paf$t.name[flipQT])
  
  ## Vectorize data transformation ##
  x <- c(rbind(paf$q.start, paf$t.start, paf$q.end, paf$t.end))
  y <- c(rbind(paf$y1, paf$y2, paf$y1, paf$y2))
  group <- rep(1:nrow(paf), each=4)
  seq.name <- c(rbind(paf$q.name, paf$t.name, paf$q.name, paf$t.name))
  seq.pos <- c(rbind(paf$q.start, paf$t.start, paf$q.end, paf$t.end))
  seq.id <- c(rbind('query', 'target', 'query', 'target'))
  n.match <- rep(paf$n.match, each=4)
  aln.len <- rep(paf$aln.len, each=4)
  mapq <- rep(paf$mapq, each=4)
  align.id <- rep(paf$seq.pair, each=4)
  direction <- rep(paf$strand, each=4)
  
  coords <- data.frame(x=x, 
                       y=y, 
                       group=group, 
                       seq.pos=seq.pos,
                       direction=direction,
                       seq.name=seq.name, 
                       seq.id=seq.id,
                       n.match=n.match,
                       aln.len=aln.len,
                       mapq=mapq,
                       align.id=align.id,
                       stringsAsFactors = FALSE)
  
  ## Get x-axis labels
  y.labels <- unique(coords$seq.name)
  y.breaks <- coords$y[match(y.labels, coords$seq.name)]
  
  coords$frac.match <- coords$n.match / coords$aln.len
  
  if (outline.miro) {
    plt <- ggplot2::ggplot(coords) +
      geom_miropeats(aes(x, y, group = group, fill=direction), alpha=0.25, color='gray', size=0.5) +
      geom_hline(yintercept = y.breaks, size=1)
  } else {
    plt <- ggplot2::ggplot(coords) +
      geom_miropeats(aes(x, y, group = group, fill=direction), alpha=0.25) +
      geom_hline(yintercept = y.breaks, size=1)
  }
  plt <- plt +
    scale_fill_manual(values = c('+' = 'chartreuse4', '-' = 'darkgoldenrod2')) +
    scale_x_continuous(expand = c(0,0), labels=comma) +
    scale_y_continuous(breaks = y.breaks, labels = y.labels) +
    xlab('Genomic position (bp)') +
    ylab('') +
    theme_minimal()
  
  return(plt)
}  

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
#'
selfdotplot_deprecated <- function(paf.table=NULL, shape='segment', sort.by='position', color.by='direction', highlight.pos=NULL, highlight.region=NULL, title=NULL) {
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

#' Read self-alignments as \code{\link{GRanges-class}} object.
#' 
#' This function reads FASTA self-alignments generated by 'nucmer' or 'minimap2' aligner.
#' This function considers proximal (leftmost) duplication as a query(s1)
#' and corresponding distal (rightmost) duplication as a target (s2)
#'
#' @param aln.coords A path to a file containing self-alignment coordinates.
#' @param format Define format of the input alignment coordinates as either 'nucmer' or minimap2 'mm2', [default: 'nucmer']
#' @param min.align.len ...
#' @param min.align.dist Keep alignment pairs with this or larger distance from each other.
#' @param collapse.overlaps Set to \code{TRUE} to merge overlapping pair of alignments with the same relative orientation.
#' @param break.paf.aln Set to \code{TRUE} in order to split CIGAR string at insertions and deletions.
#' @inheritParams breakPafAlignment
#' @return A \code{list} of \code{\link{GRanges-class}} objects.
#' @importFrom scales comma
#' @importFrom S4Vectors subjectHits queryHits endoapply
#' @import GenomicRanges
#' @author David Porubsky
#' @export
readSelfAlignments <- function(aln.coords=NULL, format='nucmer', min.align.len=1000, min.align.dist=1000, collapse.overlaps=TRUE, break.paf.aln=TRUE, min.deletion.size = 1000, min.insertion.size = 1000) {
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
    ## Get max genomic position as sequence length
    seq.len <- max(c(coords.df$s1.start, coords.df$s1.end, coords.df$s2.start, coords.df$s2.end))
  } else if (format == 'mm2') {
    ## Read in coordinates from minimap2 output
    paf.data <- readPaf(paf.file = aln.coords, restrict.paf.tags = 'cg')
    ## Get sequence length
    seq.len <- unique(paf.data$q.len)
    ## Due to the minimap2 self-alignment redundancy keep only alignments where query start is smaller than the target start
    paf.data <- paf.data[paf.data$q.start < paf.data$t.start,]
    ## Break paf alignments
    if (break.paf.aln) {
      paf.data <- breakPaf(paf.table = paf.data, min.deletion.size = min.deletion.size, min.insertion.size = min.insertion.size, report.sv = TRUE)
      paf.data.sv <- paf.data$SVs
      paf.data <- paf.data$M
    } else {
      paf.data.sv <- NULL
    } 
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
  
  ## Data transformations ##
  ##########################
  ## Add alignment directionality
  coords.df$dir <- 'rev'
  forw.mask <- (coords.df$s1.start < coords.df$s1.end) & (coords.df$s2.start < coords.df$s2.end)
  coords.df$dir[forw.mask] <- 'forw'
  
  ## Make sure s1 coordinates are always smaller than s2 coordinates
  toFlip <- which(pmin(coords.df$s1.start, coords.df$s1.end) > pmax(coords.df$s2.start, coords.df$s2.end))
  coords.df[toFlip,] <- transform(coords.df[toFlip,], 's1.start' = s2.start, 's1.end' = s2.end, 's2.start' = s1.start, 's2.end' = s1.end)
  
  ## Get distance between alignments
  coords.df <- transform(coords.df, dist = pmin(s2.start, s2.end) - pmax(s1.start, s1.end))
  
  ## Data filtering ##
  ####################
  ## Keep diagonals (self-alignments) for export
  # diagonals.df <- coords.df[coords.df$s1.start == coords.df$s2.start & coords.df$s1.end == coords.df$s2.end,]
  # if (nrow(diagonals.df) > 0) {
  #   diagonals.gr <- GenomicRanges::makeGRangesFromDataFrame(diagonals.df, seqnames.field = 's1.id', start.field = 's1.start', end.field = 's1.end')
  #   names(diagonals.gr) <- NULL
  # } else {
  #   diagonals.gr <- NULL
  # }  
  ## Remove diagonals (self-alignments) with the same start and end position
  coords.df <- coords.df[!(coords.df$s1.start == coords.df$s2.start & coords.df$s1.end == coords.df$s2.end),]
  ## Filter by alignment length
  if (min.align.len > 0) {
    coords.df <- coords.df[coords.df$s1.width >= min.align.len & coords.df$s2.width >= min.align.len,]
  }
  ## Remove duplicated ranges
  #paste0(coords.df$s1.start, coords.df$s1.end) %in% paste0(coords.df$s2.start, coords.df$s2.end)
  xmin <- apply(coords.df[,c('s1.start', 's1.end', 's2.start', 's2.end')], 1, min)
  xmax <- apply(coords.df[,c('s1.start', 's1.end', 's2.start', 's2.end')], 1, max)
  coords.df <- coords.df[!duplicated(paste(xmin, xmax, sep = '_')),]
  ## Filter by alignment distance
  if (min.align.dist > 0) {
    coords.df <- coords.df[coords.df$dist >= min.align.dist,]
  }
  
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
    ## Keep unique overlap pairs
    mask1 <- !duplicated(paste0(pmin(S4Vectors::queryHits(hits1), S4Vectors::subjectHits(hits1)), pmax(S4Vectors::queryHits(hits1), S4Vectors::subjectHits(hits1))))
    hits1 <- hits1[mask1]
    mask2 <- !duplicated(paste0(pmin(S4Vectors::queryHits(hits2), S4Vectors::subjectHits(hits2)), pmax(S4Vectors::queryHits(hits2), S4Vectors::subjectHits(hits2))))
    hits2 <- hits2[mask2]
    ## Keep pairs overlapping each other
    mask1 <- paste0(S4Vectors::queryHits(hits1), S4Vectors::subjectHits(hits1)) %in% paste0(S4Vectors::queryHits(hits2), S4Vectors::subjectHits(hits2))
    mask2 <- paste0(S4Vectors::queryHits(hits2), S4Vectors::subjectHits(hits2)) %in% paste0(S4Vectors::queryHits(hits1), S4Vectors::subjectHits(hits1))
    hits1 <- hits1[mask1]
    hits2 <- hits2[mask2]
    ## Keep pairs having the same dir
    mask <- s1.gr$dir[S4Vectors::queryHits(hits1)] == s1.gr$dir[S4Vectors::subjectHits(hits1)] & s2.gr$dir[S4Vectors::queryHits(hits2)] == s2.gr$dir[S4Vectors::subjectHits(hits2)]
    hits1 <- hits1[mask]
    hits2 <- hits2[mask]
    
    #if (length(hits1) > 1 & length(hits2) > 1) {
    while (length(hits1) > 1 & length(hits2) > 1) {
      ## Get groups alignments overlapping each other
      groups <- list()
      #for (i in seq_along(hits1)) {
      for (i in order(S4Vectors::subjectHits(hits1))) {
        pair <- hits1[i]
        pair <- c(S4Vectors::queryHits(pair), S4Vectors::subjectHits(pair))
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
      names(grp) <- rep(1:length(groups), times=lengths(groups))
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
      s2.collapse.gr <- unlist(S4Vectors::endoapply(s2.collapse.grl, range))
      s2.collapse.gr$dir <- sapply(s2.collapse.grl, function(x) unique(x$dir))
      s2.collapse.gr$group <- unique(s2.gr$group[s2.gr$group > 0])
      names(s2.collapse.gr) <- NULL
      ## Replace collapsed ranges
      s1.gr <- c(s1.gr[s1.gr$group == 0], s1.collapse.gr)
      s2.gr <- c(s2.gr[s2.gr$group == 0], s2.collapse.gr)
      
      ## Recalculate self-alignments
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
      ## Keep pairs having the same dir
      mask <- s1.gr$dir[queryHits(hits1)] == s1.gr$dir[subjectHits(hits1)] & s2.gr$dir[queryHits(hits2)] == s2.gr$dir[subjectHits(hits2)]
      hits1 <- hits1[mask]
      hits2 <- hits2[mask]
    }
  }  
  ## Remove remaining self-overlapping ranges
  # if (min.align.dist > 0) {
  #   s2.copy.gr <- s2.gr
  #   seqlevels(s2.copy.gr) <- 's1'
  #   mask <- which(IRanges::distance(s1.gr, s2.copy.gr) == 0)
  #   if (length(mask) > 0) {
  #     s1.gr <- s1.gr[-mask]
  #     s2.gr <- s2.gr[-mask]
  #   }    
  # }
  
  ## Prepare object of self-alignments for export
  if (length(s1.gr) > 0 & length(s2.gr) > 0) {
    self.gr <- GenomicRanges::GRanges(seqnames=unique(coords.df$s1.id), ranges=ranges(s1.gr[,0]))
    #self.gr$s2 <- s2.gr[,0]
    self.gr$s2 <- GenomicRanges::GRanges(seqnames=unique(coords.df$s1.id), ranges=ranges(s2.gr[,0]))
    strand(self.gr) <- '+'
    strand(self.gr$s2) <- ifelse(s1.gr$dir == 'forw', '+', '-')
    ## Add sequence length
    seqlengths(self.gr) <- seq.len
    seqlengths(self.gr$s2) <- seq.len
  } else {
    self.gr <- GenomicRanges::GRanges()
  } 
  ## Prepare reported SVs for export
  if (format == 'mm2') {
    if (nrow(paf.data.sv) > 0) {
      sv.gr <- GenomicRanges::GRanges(seqnames=paf.data.sv$q.name, ranges=IRanges(start=paf.data.sv$q.start, end=paf.data.sv$q.end), strand = strand('*'))
      sv.gr$s2 <- GenomicRanges::GRanges(seqnames=paf.data.sv$t.name, ranges=IRanges(start=paf.data.sv$t.start, end=paf.data.sv$t.end), strand = strand('*'))
      sv.gr$sv.type <- gsub(paf.data.sv$cg, pattern = '\\d+', replacement = '')
      sv.gr$sv.size <- gsub(paf.data.sv$cg, pattern = '[A-Z,=]', replacement = '', ignore.case = TRUE)
      ## Add sequence length
      seqlengths(sv.gr) <- seq.len
    } else {
      sv.gr <- GenomicRanges::GRanges()
    } 
  } else {
    sv.gr <- GenomicRanges::GRanges()
  }
  
  ## Return self-alignments
  return(list('SelfAlnPairs' = self.gr, 'SVs' = sv.gr))
} 