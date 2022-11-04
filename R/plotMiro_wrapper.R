#' Visualize PAF alignments.
#' 
#' This function takes PAF output file from minimap2 alignments, and visualize the alignments
#' in miropeat style. 
#'
#' @param color.by 
#' @param flip.alignment Set to \code{TRUE} if query PAF alignments should be flipped.
#' @inheritParams readPaf
#' @inheritParams filterPaf
#' @inheritParams breakPaf
#' @inheritParams flipPaf
#' @return A \code{ggplot2} object or \code{list} of  \code{ggplot2} objects.
#' @importFrom scales comma
#' @importFrom wesanderson wes_palette
#' @importFrom gggenes geom_gene_arrow
#' @importFrom ggnewscale new_scale_fill new_scale_color
#' @author David Porubsky
#' @export
#' @examples
#'## Get PAF to plot
#'paf.file <- system.file("extdata", "test1.paf", package="SVbyEye")
#'## Make a plot
#'## Color by alignment directionality
#'plotMiro(paf.file = paf.file, color.by = 'direction')
#'## Color by fraction of matched bases in each alignment
#'plotMiro(paf.file = paf.file, color.by = 'fraction.matches')
#'
plotMiro <- function(paf.file = paf.file, min.mapq = 10, min.align.len = 100, min.align.n = 1, target.region = NULL, query.region = NULL, min.deletion.size=NULL, min.insertion.size=NULL, drop.self.align = FALSE, majority.strand = '+', color.by = 'direction', flip.alignment = FALSE) {
  ## Check user input
  ## Make sure submitted files exists
  if (!file.exists(paf.file)) {
    stop(paste0("Submitted 'paf.file' ", paf.file, " does not exists !!!"))
  }
  
  ## Load PAF file
  paf <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = 'cg') 
  ## Filter PAF file
  paf <- filterPaf(paf.table = paf, min.mapq = min.mapq, min.align.len = min.align.len, min.align.n = min.align.n, target.region = target.region, query.region = query.region, drop.self.align = drop.self.align)
  ## Break PAF at insertion/deletions defined in cigar string
  paf.l <- breakPaf(paf.table = paf, min.deletion.size = min.deletion.size, min.insertion.size = min.insertion.size, collapse.mismatches = TRUE, report.sv = TRUE)
  paf <- paf.l$M
  paf.svs <- paf.l$SVs
  ## Flip PAF alignments given the user defined majority orientation
  paf <- flipPaf(paf.table = paf, majority.strand = majority.strand, force=flip.alignment)
  ## Convert PAF alignments to plotting coordinates
  paf <- paf2coords(paf.table = paf)
  
  ## Process data per alignment
  coords.data.l <- split(paf, paf$align.id)
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
    ## Make sure axis labels are always positive numbers
    q.labels <- abs(q.labels)
    t.labels <- abs(t.labels)
    
    ## Get x-axis labels
    seq.labels <- c(unique(coords$seq.name[coords$seq.id == 'query']), 
                    unique(coords$seq.name[coords$seq.id == 'target']))
    
    ## Color alignments by variable
    if (color.by == 'direction') {
      plt <- ggplot2::ggplot(coords) +
        geom_miropeats(aes(x, y, group = group, fill=direction), alpha=0.5) +
        scale_fill_manual(values = c('-' = 'cornflowerblue', '+' = 'forestgreen'))
    } else if (color.by == 'fraction.matches') {
      coords$frac.match <- coords$n.match / coords$aln.len
      plt <- ggplot2::ggplot(coords) +
        geom_miropeats(aes(x, y, group = group, fill=frac.match), alpha=0.5) +
        scale_fill_gradient(low = 'gray', high = 'red', name='Fraction\nmatches')
    } else if (color.by == 'mapq') {
      plt <- ggplot2::ggplot(coords) +
        geom_miropeats(aes(x, y, group = group, fill=mapq), alpha=0.5) +
        scale_fill_gradient(low = 'gray', high = 'red')
    } else {
      plt <- ggplot2::ggplot(coords) +
        geom_miropeats(aes(x, y, group = group), alpha=0.5, fill='gray')
    }
    
    ## Add x and y scales
    plt <- plt +
      scale_y_continuous(breaks = c(1, 2), labels = seq.labels) +
      scale_x_continuous(breaks = q.breaks, labels = scales::comma(q.labels),
                         sec.axis = sec_axis(trans = y ~ ., breaks = t.breaks, labels = scales::comma(t.labels)), expand = c(0,0)) +
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
      gggenes::geom_gene_arrow(data=plt.df, ggplot2::aes(xmin = start, xmax = end, y = y, color= direction, fill = direction), arrowhead_height = unit(3, 'mm')) +
      scale_fill_manual(values = c('-' = 'cornflowerblue', '+' = 'forestgreen')) +
      scale_color_manual(values = c('-' = 'cornflowerblue', '+' = 'forestgreen')) +
      theme_bw()
    
    ## Store plot
    plots[[i]] <- plt
  }
  if (length(plots) == 1) {
    return(plots[[1]])
  } else {
    return(plots)
  }  
}


#' Add annotation ranges to miropeat style plot.
#'
#' This function takes a \code{ggplot2} object generated using \code{\link{plotMiro}} function and adds extra annotation on top of query 
#' or target coordinates. These ranges are specified in 'annot.gr' object and are visualized either as arrowheads or rectangles.
#'
#' @param ggplot.obj A \code{ggplot2} object generated using \code{\link{plotMiro}} function.
#' @param annot.gr A \code{\link{GRanges-class}} object with a set of ranges to be added as extra annotation.
#' @param shape A user defined shape ranges in 'annot.gr' are visualized, either 'arrowhead' or 'rectangle'.
#' @param fill.by A name of an extra field present in 'annot.gr' to be used to define color scheme.   
#' @param coordinate.space A coordinate space ranges in 'annot.gr' are reported, either 'target' or 'query'.
#' @return A \code{list} of miropeat style plots.
#' @importFrom scales comma
#' @importFrom wesanderson wes_palette
#' @importFrom gggenes geom_gene_arrow
#' @importFrom ggnewscale new_scale_fill new_scale_color
#' @author David Porubsky
#' @export
#' @examples
#'## Get PAF to plot
#'paf.file <- system.file("extdata", "test1.paf", package="SVbyEye")
#'## Make a plot
#'plt <- plotMiro(paf.file = paf.file)
#'## Load target annotation file
#'target.annot <- system.file("extdata", "test1_target_annot.txt", package="SVbyEye")
#'target.annot.df <- read.table(target.annot, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
#'target.annot.gr <- makeGRangesFromDataFrame(target.annot.df)
#'## Add target annotation as arrowhead
#'plt <- add_annotation(ggplot.obj = plt, annot.gr = target.annot.gr, coordinate.space='target')
#'## Load query annotation file
#'query.annot <- system.file("extdata", "test1_query_annot.txt", package="SVbyEye")
#'query.annot.df <- read.table(query.annot, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
#'query.annot.gr <- makeGRangesFromDataFrame(query.annot.df)
#'## Add query annotation as rectangle
#'add_annotation(ggplot.obj = plt, annot.gr = query.annot.gr, shape='rectangle', coordinate.space='query')
#'## Lift target annotation to query and plot
#'lifted.annot.gr <- liftRangesToAlignment(gr = target.annot.gr, paf.file = paf.file, direction = 'target2query')
#'add_annotation(ggplot.obj = plt, annot.gr = lifted.annot.gr, shape='rectangle', coordinate.space='query')
#'
add_annotation <- function(ggplot.obj=NULL, annot.gr=NULL, shape='arrowhead', fill.by=NULL, coordinate.space='target') {
  ## Get plotted data
  gg.data <- ggplot.obj$data
  target.id <- unique(gg.data$seq.name[gg.data$seq.id == 'target'])
  query.id <- unique(gg.data$seq.name[gg.data$seq.id == 'query'])
  ## Get query and target coordinate ranges
  t.range <- range(gg.data$seq.pos[gg.data$seq.id == 'target'])
  q.range <- range(gg.data$seq.pos[gg.data$seq.id == 'query'])
  
  if (coordinate.space == 'query') {
    ## Restrict to plotted query ranges
    annot.gr <- annot.gr[start(annot.gr) >= q.range[1] & end(annot.gr) <= q.range[2]]
    if (length(annot.gr) > 0) {
      ## Covert query to target coordinates
      new.start <- SVbyEye::q2t(x = start(annot.gr), q.range = q.range, t.range = t.range)
      new.end <- SVbyEye::q2t(x = end(annot.gr), q.range = q.range, t.range = t.range)
      new.ranges <- IRanges::IRanges(start = new.start, end = new.end)
      ranges(annot.gr) <- new.ranges
    }  
  }  
    
  ## Get x and y-axis limits
  xlim <- layer_scales(ggplot.obj)$x$range$range
  ylim <- layer_scales(ggplot.obj)$y$range$range
  
  ## Get annotation track offset
  if (coordinate.space == 'target') {
    y.offset <- max(ylim) + 0.05
  } else if (coordinate.space == 'query') {
    y.offset <- min(ylim) - 0.05
  } else {
    stop("Please specify 'coordinate.space' as either 'target' or 'query' !!!")
  }
  
  
  if (!is.null(annot.gr) & class(annot.gr) == 'GRanges') {
    ## Restrict to plot x-limits
    annot.gr <- annot.gr[start(annot.gr) >= xlim[1] & end(annot.gr) <= xlim[2]]
    if (length(annot.gr) > 0) {
      ## Convert to data frame object
      annot.df <- as.data.frame(annot.gr)
      if (!is.null(fill.by)) {
        if (fill.by %in% colnames(annot.df)) {
          ## Get color scales
          if (is.numeric(annot.df[,eval(fill.by)])) {
            pal <- wesanderson::wes_palette("Zissou1", 100, type = "continuous")
            col.scale <- 'gradient'
          } else {
            n.uniq <- length(unique(annot.df[,eval(fill.by)]))
            if (length(n.uniq) <= 20) {
              pal <- wesanderson::wes_palette("Zissou1", n.uniq, type = "continuous")
              col.scale <- 'discrete'
            } else {
              warning("More than 20 color levels, legend won't be reported!!!")
              col.scale <- 'discreteNoLegend'
            }  
          }
        }  
      } else {
        pal <- 'black'
        col.scale <- 'discreteNoLegend'
      }  
      
      if (shape == 'arrowhead') {
        plt <- ggplot.obj + ggnewscale::new_scale_fill() + ggnewscale::new_scale_color() +
          geom_arrowhead(data=annot.df, aes(xmin=start, xmax=end, y=y.offset, color=eval(fill.by), fill=eval(fill.by)))
      } else if (shape == 'rectangle') {
        plt <- ggplot.obj + ggnewscale::new_scale_fill() + ggnewscale::new_scale_color() +
          #geom_rect(data=annot.df, aes(xmin=start, xmax=end, ymin=y.offset, ymax=y.offset + 0.01, color=eval(fill.by), fill=eval(fill.by)))
          geom_roundrect(data=annot.df, aes(xmin=start, xmax=end, y=y.offset + 0.01, color=eval(fill.by), fill=eval(fill.by)), radius = grid::unit(0, "mm"))
      }  
      
      if (col.scale == 'gradient') {
        plt <- plt + scale_fill_gradientn(colours = pal) + scale_color_gradientn(colours = pal)
      } else if (col.scale == 'discrete') {
        plt <- plt + scale_fill_manual(values = pal) + scale_color_manual(values = pal)
      } else if (col.scale == 'discreteNoLegend') {
        plt <- plt + scale_fill_manual(values = pal, guide='none') + scale_color_manual(values = pal, guide='none')
      }
      return(plt)
    } else {
      warning("None of the ranges reported in 'annot.gr' falls into the x-axis limits !!!")
      return(ggplot.obj)
    } 
  } else {
    return(ggplot.obj)
    warning("Make sure that 'annot.gr' is 'GRanges' class object !!!")
  }
}

