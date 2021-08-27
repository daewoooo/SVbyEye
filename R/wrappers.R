#' Plot PAF file as miropeat alignments
#' 
#' This function takes PAF output file from minimap2 alignemts, and visualize the alignments
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
#' @export
plotMiro <- function(paf.file = paf.file, min.mapq = 10, min.align.len = 100, min.align.n = 1, sd.annot = NULL) {
  ## Load PAF file
  coords.data <- paf2coords(paf.file = paf.file, min.mapq = 10, min.align.len = 100, min.align.n = 1)  
  ## Process data per alignment
  coords.data.l <- split(coords.data, coords.data$align.id)
  plots <- list()
  for (i in seq_along(coords.data.l)) {
    coords <- coords.data.l[[i]]
    target.seqname <- unique(coords$seq.name[coords$seq.id == 'target'])
    ## Get y-axis labels
    q.range <- range(coords$seq.pos[coords$seq.id == 'query'])
    t.range <- range(coords$seq.pos[coords$seq.id == 'target'])
    q.labels <- pretty(q.range)
    t.labels <- pretty(t.range)
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
