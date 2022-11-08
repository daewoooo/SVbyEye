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
plotAVA <- function(paf.file = NULL, seqnames.order=NULL, seqnames.grep=NULL, target.region=NULL, min.align.n=1, min.mapq=0, min.align.len=1000, drop.self.align=TRUE, majority.strand=NULL, outline.miro=FALSE) {
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

