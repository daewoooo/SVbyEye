#' Visualize PAF alignments.
#' 
#' This function takes PAF output file from minimap2 alignments, and visualize the alignments
#' in miropeat style. 
#'
#' @param highlight.sv Visualize alignment embedded structural variation either as an outlined ('outline') or filled ('fill') miropeats. 
#' @param color.by Color alignments either by directionality ('direction') or fraction of matched base pairs ('identity').
#' @param outline.alignments Set to \code{TRUE} if boundaries of each alignment should be highlighted by gray outline.
#' @inheritParams breakPaf
#' @inheritParams pafAlignmentToBins
#' @inheritParams paf2coords
#' @return A \code{ggplot2} object.
#' @importFrom scales comma
#' @importFrom wesanderson wes_palette
#' @importFrom gggenes geom_gene_arrow
#' @importFrom ggnewscale new_scale_fill new_scale_color
#' @author David Porubsky
#' @export
#' @examples
#'## Get PAF to plot
#'paf.file <- system.file("extdata", "test1.paf", package="SVbyEye")
#'## Read in PAF 
#'paf.table <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = 'cg')
#'## Optional steps include PAF filtering and flipping query coordinates (see filterPaf and flipPaf function documentation)
#'## Make a plot
#'## Color by alignment directionality
#'plotMiro(paf.table = paf.table, color.by = 'direction')
#'## Color by fraction of matched bases in each alignment
#'plotMiro(paf.table = paf.table, color.by = 'identity')
#'## Outline PAF alignments
#'plotMiro(paf.table = paf.table, outline.alignments = TRUE)
#'## Offset target PAF alignments
#'plotMiro(paf.table = paf.table, offset.alignments = TRUE)
#'## Highlight structural variants
#'paf.file <- system.file("extdata", "test3.paf", package="SVbyEye")
#'paf.table <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = 'cg')
#'plotMiro(paf.table = paf.table, min.deletion.size=50, highlight.sv='outline')
#'## Bin PAF alignments into user defined bin and color them by sequence identity (% of matched bases)
#'plotMiro(paf.table = paf.table, binsize=10000)
#'
#plotMiro <- function(paf.table, min.mapq = 10, min.align.len = 100, min.align.n = 1, target.region = NULL, query.region = NULL, drop.self.align = FALSE, min.deletion.size=NULL, min.insertion.size=NULL, highlight.sv=NULL, majority.strand = '+', color.by = 'direction', flip.alignment = FALSE) {
plotMiro <- function(paf.table, min.deletion.size=NULL, min.insertion.size=NULL, highlight.sv=NULL, binsize=NULL, color.by='direction', outline.alignments=FALSE, offset.alignments=FALSE) {
  ## Check user input
  ## Make sure submitted paf.table has at least 12 mandatory fields
  if (ncol(paf.table) >= 12) {
    paf <- paf.table
    paf$direction.flip <- FALSE
  } else {
    stop('Submitted PAF alignments do not contain a minimum of 12 mandatory fields, see PAF file format definition !!!')
  }
  ## Make sure submitted files exists
  # if (!file.exists(paf.file)) {
  #   stop(paste0("Submitted 'paf.file' ", paf.file, " does not exists !!!"))
  # }
  
  ## Load PAF file
  #paf <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = 'cg') 
  ## Filter PAF file
  #paf <- filterPaf(paf.table = paf, min.mapq = min.mapq, min.align.len = min.align.len, min.align.n = min.align.n, target.region = target.region, query.region = query.region, drop.self.align = drop.self.align)
  
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
    if (!is.null(highlight.sv)) {
      highlight.sv <- NULL
      warning("Please specify 'min.deletion.size' and 'min.insertion.size' in order to make parameter 'highlight.sv' to work !!!")
    }
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
  
  ## Flip PAF alignments given the user defined majority orientation
  #paf <- flipPaf(paf.table = paf, majority.strand = majority.strand, force=flip.alignment)
  
  ## Convert PAF alignments to plotting coordinates
  coords <- paf2coords(paf.table = paf, offset.alignments = offset.alignments)
  
  ## Prepare data for plotting
  target.seqname <- unique(coords$seq.name[coords$seq.id == 'target'])
  ## Get y-axis labels
  q.range <- range(coords$seq.pos[coords$seq.id == 'query'])
  t.range <- range(coords$seq.pos[coords$seq.id == 'target'])
  ## Adjust target ranges given the size difference with respect to query ranges
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
    plt <- ggplot2::ggplot(coords[coords$ID == 'M',]) +
      geom_miropeats(aes(x, y, group = group, fill=direction), alpha=0.5) +
      scale_fill_manual(values = c('-' = 'cornflowerblue', '+' = 'forestgreen'), name='Alignment\ndirection')
  } else if (color.by == 'identity') {
    coords$identity <- (coords$n.match / coords$aln.len) * 100
    coords$identity[is.nan(coords$identity) | is.na(coords$identity)] <- 0
    ## Define color scheme
    coords.l <- getColorScheme(data.table = coords, value.field = 'identity', breaks=c(90, 95, 99, 99.5, 99.6, 99.7, 99.8, 99.9))
    coords <- coords.l$data
    colors <- coords.l$colors
    
    plt <- ggplot2::ggplot(coords[coords$ID == 'M',]) +
      geom_miropeats(aes(x, y, group = group, fill=col.levels), alpha=0.5) +
      scale_fill_manual(values = colors, drop=FALSE, name='Identity')
  } else if (color.by == 'mapq') {
    plt <- ggplot2::ggplot(coords[coords$ID == 'M',]) +
      geom_miropeats(aes(x, y, group = group, fill=mapq), alpha=0.5) +
      scale_fill_gradient(low = 'gray', high = 'red')
  } else {
    plt <- ggplot2::ggplot(coords[coords$ID == 'M',]) +
      geom_miropeats(aes(x, y, group = group), alpha=0.5, fill='gray')
  }
  
  ## Add alignment outlines 
  if (outline.alignments) {
    plt <- plt + geom_miropeats(data=coords[coords$ID == 'M',], aes(x, y, group = group),  fill=NA, color='gray', size=0.25)
  }  
  
  ## Add indels
  if (!is.null(highlight.sv)) {
    if (nrow(coords[coords$ID != 'M',]) > 0) { 
      ## Add SVs to the plot
      if (highlight.sv == 'outline') {
        plt <- plt + ggnewscale::new_scale_color() +
          geom_miropeats(data=coords[coords$ID != 'M',], aes(x, y, group = group, color=ID), fill=NA, alpha=0.5, inherit.aes = FALSE) +
          scale_color_manual(values = c('DEL' = 'firebrick3', 'INS' = 'dodgerblue3'), name='SV class')
      } else if (highlight.sv == 'fill') {
        plt <- plt + ggnewscale::new_scale_fill() +
          geom_miropeats(data=coords[coords$ID != 'M',], aes(x, y, group = group, fill=ID), alpha=0.5, inherit.aes = FALSE) +
          scale_fill_manual(values = c('DEL' = 'firebrick3', 'INS' = 'dodgerblue3'), name='SV class')
      } else {
        warning("Parameter 'highlight.sv' can only take values 'outline' or 'fill', see function documentation!!!")
      } 
    } else {
      warning("There are no SVs to highlight. Try to decrease 'min.deletion.size' and 'min.insertion.size' values!!!")
    }  
  }
  
  ## Add x and y scales
  suppressMessages(
    plt <- plt +
      scale_y_continuous(breaks = c(1, 2), labels = seq.labels) +
      scale_x_continuous(breaks = q.breaks, labels = scales::comma(q.labels),
                         sec.axis = sec_axis(trans = y ~ ., breaks = t.breaks, labels = scales::comma(t.labels)), expand = c(0,0)) +
      xlab('Genomic position (bp)') +
      ylab('')
  )
  
  ## Add arrows to mark start and end of each alignment
  ## Always used unbinned version of PAF alignments
  coords.arrow <- paf2coords(paf.table = paf.copy, offset.alignments = offset.alignments)
  start <- coords.arrow$x[c(T, T, F, F)]
  end <- coords.arrow$x[c(F, F, T, T)]
  y <- coords.arrow$y[c(T, T, F, F)]
  group <- coords.arrow$group[c(T, T, F, F)]
  plt.df <- data.frame(start=start, end=end, y=y, group=group)
  plt.df$direction <- ifelse(plt.df$start < plt.df$end, '+', '-')
  
  plt <- plt + ggnewscale::new_scale_fill() + ggnewscale::new_scale_color() +
    gggenes::geom_gene_arrow(data=plt.df, ggplot2::aes(xmin = start, xmax = end, y = y, color= direction, fill = direction), arrowhead_height = unit(3, 'mm')) +
    scale_fill_manual(values = c('-' = 'cornflowerblue', '+' = 'forestgreen'), name='Alignment\ndirection') +
    scale_color_manual(values = c('-' = 'cornflowerblue', '+' = 'forestgreen'), name='Alignment\ndirection')
  
  ## Set the theme
  theme_miro <- theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(), 
                      axis.line.x = element_line(size = 1),
                      axis.ticks.x = element_line(size=1),
                      axis.ticks.length.x = unit(2, 'mm'))
  plt <- plt + theme_miro
  
  ## Return final plot
  return(plt)
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
#' @param color.palette User defined discrete color palette as named character vector (elements = colors, names = discrete levels).  
#' @param coordinate.space A coordinate space ranges in 'annot.gr' are reported, either 'target' or 'query'.
#' @return A \code{ggplot2} object.
#' @import ggplot2
#' @importFrom scales comma
#' @importFrom wesanderson wes_palette
#' @importFrom gggenes geom_gene_arrow
#' @importFrom ggnewscale new_scale_fill new_scale_color
#' @author David Porubsky
#' @export
#' @examples
#'## Get PAF to plot
#'paf.file <- system.file("extdata", "test1.paf", package="SVbyEye")
#'## Read in PAF
#'paf.table <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = 'cg')
#'## Make a plot
#'plt <- plotMiro(paf.table = paf.table)
#'## Load target annotation file
#'target.annot <- system.file("extdata", "test1_target_annot.txt", package="SVbyEye")
#'target.annot.df <- read.table(target.annot, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
#'target.annot.gr <- makeGRangesFromDataFrame(target.annot.df)
#'## Add target annotation as arrowhead
#'plt <- add_annotation(ggplot.obj = plt, annot.gr = target.annot.gr, coordinate.space = 'target')
#'## Load query annotation file
#'query.annot <- system.file("extdata", "test1_query_annot.txt", package="SVbyEye")
#'query.annot.df <- read.table(query.annot, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
#'query.annot.gr <- makeGRangesFromDataFrame(query.annot.df)
#'## Add query annotation as rectangle
#'add_annotation(ggplot.obj = plt, annot.gr = query.annot.gr, shape = 'rectangle', coordinate.space = 'query')
#'## Lift target annotation to query and plot
#'lifted.annot.gr <- liftRangesToAlignment(gr = target.annot.gr, paf.file = paf.file, direction = 'target2query')
#'add_annotation(ggplot.obj = plt, annot.gr = lifted.annot.gr, shape = 'rectangle', coordinate.space = 'query')
#'## Add segmental duplication annotation
#'plt <- plotMiro(paf.table = paf.table)
#'sd.annot <- system.file("extdata", "test1.sd.annot.RData", package="SVbyEye")
#'sd.annot.gr <- get(load(sd.annot))
#'## Create a custom discrete levels
#'sd.categ <- findInterval(sd.annot.gr$fracMatch, vec = c(0.95, 0.98, 0.99))
#'sd.categ <- dplyr::recode(sd.categ, '0' = '<95%', '1' = '95-98%', '2' = '98-99%', '3'='>=99%')
#'sd.categ <- factor(sd.categ, levels=c('<95%', '95-98%', '98-99%', '>=99%'))
#'sd.annot.gr$sd.categ <- sd.categ
## Create a custom color palette
#'color.palette <- c('<95%' = 'gray72', '95-98%' = 'gray47', '98-99%' = '#cccc00', '>=99%' = '#ff6700')
#'## Add annotation to the plot
#'add_annotation(ggplot.obj = plt, annot.gr = sd.annot.gr, fill.by = 'sd.categ', color.palette = color.palette, coordinate.space = 'target')
#'
add_annotation <- function(ggplot.obj=NULL, annot.gr=NULL, shape='arrowhead', fill.by=NULL, color.palette=NULL, coordinate.space='target') {
  ## Get plotted data
  gg.data <- ggplot.obj$data
  target.id <- unique(gg.data$seq.name[gg.data$seq.id == 'target'])
  query.id <- unique(gg.data$seq.name[gg.data$seq.id == 'query'])
  ## Get query and target coordinate ranges
  t.range <- range(gg.data$seq.pos[gg.data$seq.id == 'target'])
  q.range <- range(gg.data$seq.pos[gg.data$seq.id == 'query'])
  ## Adjust target ranges given the size difference with respect to query ranges
  range.offset <- diff(q.range) - diff(t.range)
  t.range[2] <- t.range[2] + range.offset ## Make a start position as offset and change only end position
  
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
  xlim <- ggplot2::layer_scales(ggplot.obj)$x$range$range
  ylim <- ggplot2::layer_scales(ggplot.obj)$y$range$range
  
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
          ## Get color scales ##
          ## Define continuous color scale
          if (is.numeric(annot.df[,eval(fill.by)])) {
            pal <- wesanderson::wes_palette("Zissou1", 100, type = "continuous")
            col.scale <- 'gradient'
          ## Define discrete color scale
          } else {
            dicrete.levels <- unique(annot.df[,eval(fill.by)])
            n.uniq <- length(dicrete.levels)
            ## Get user define discrete color palette
            if (all(dicrete.levels %in% names(color.palette))) {
              pal <- color.palette
              col.scale <- 'discrete'
            ## Create default discrete color palette  
            } else {
              if (n.uniq <= 20) {
                pal <- wesanderson::wes_palette("Zissou1", n.uniq, type = "continuous")
                col.scale <- 'discrete'
              } else {
                warning("More than 20 color levels, legend won't be reported!!!")
                col.scale <- 'discreteNoLegend'
              }  
            } 
          }
        } else {
          stop("User defined 'fill.by' value is not a valid field in 'annot.gr' !!!")
        } 
      } else {
        pal <- 'black'
        col.scale <- 'discreteNoLegend'
      }  
      
      if (shape == 'arrowhead') {
        ## Make sure field 'strand' is defined
        if (!'strand' %in% colnames(annot.df)) {
          annot.df$strand <- '*'
        }
        plt <- ggplot.obj + ggnewscale::new_scale_fill() + ggnewscale::new_scale_color() +
          geom_arrowhead(data=annot.df, aes_string(xmin='start', xmax='end', y='y.offset', color=eval(fill.by), fill=eval(fill.by), strand='strand'))
      } else if (shape == 'rectangle') {
        plt <- ggplot.obj + ggnewscale::new_scale_fill() + ggnewscale::new_scale_color() +
          #geom_rect(data=annot.df, aes(xmin=start, xmax=end, ymin=y.offset, ymax=y.offset + 0.01, color=eval(fill.by), fill=eval(fill.by)))
          geom_roundrect(data=annot.df, aes_string(xmin='start', xmax='end', y=y.offset + 0.01, color=eval(fill.by), fill=eval(fill.by)), radius = grid::unit(0, "mm"))
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


#' Flip query annotation ranges 
#' 
#' This function takes loaded PAF alignments using \code{\link{readPaf}} function and postprocessed using 
#' \code{\link{flipPaf}} function. In case the PAF alignments were flipped ranges defined in 'query.annot.gr'
#' will be flipped accordingly to match query coordinates defined in 'paf.table'.
#' 
#' @param query.annot.gr A \code{\link{GRanges-class}} object with a set of ranges in query coordinates. See function
#' \code{\link{liftRangesToAlignment}} if coordinates need to be lifted from target space.
#' @inheritParams breakPafAlignment
#' @return A \code{\link{GRanges-class}} object.
#' @author David Porubsky
#' @export
#' @examples 
#'## Get PAF to process
#'paf.file <- system.file("extdata", "test1.paf", package="SVbyEye")
#'## Read in PAF
#'paf.table <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = 'cg')
#'## Flip PAF alignments
#'paf.table <- flipPaf(paf.table = paf.table, force=TRUE)
#'## Load query annotation file
#'query.annot <- system.file("extdata", "test1_query_annot.txt", package="SVbyEye")
#'query.annot.df <- read.table(query.annot, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
#'query.annot.gr <- makeGRangesFromDataFrame(query.annot.df)
#'## Synchronize orientation of query annotation file with flipped PAF alignments
#'flipQueryAnnotation(paf.table = paf.table, query.annot.gr=query.annot.gr)
#'
flipQueryAnnotation <- function(paf.table, query.annot.gr=NULL) {
  ## Check user input ##
  ## Check if submitted query annotation is GRanges object
  if (!is.null(query.annot.gr)) {
    if (class(query.annot.gr) != 'GRanges') {
      stop("Parameter 'query.annot.gr' has to be of 'GenomicRanges' object !!!")
    }
  }
  ## Make sure submitted paf.table has at least 12 mandatory fields
  if (ncol(paf.table) >= 12) {
    paf <- paf.table
    paf$direction.flip <- FALSE
  } else {
    stop('Submitted PAF alignments do not contain a minimum of 12 mandatory fields, see PAF file format definition !!!')
  }
  
  ## Process per query name
  paf.table.l <- split(paf.table, paf.table$q.name)
  query.annot.grl <- split(query.annot.gr, GenomeInfoDb::seqnames(query.annot.gr)) 
  ## Keep only non-empty list elements
  query.annot.grl <- query.annot.grl[lengths(query.annot.grl) > 0]
  
  flipped.annot <- list()
  for (q.name in names(query.annot.grl)) {
    annot.gr <- query.annot.grl[[q.name]]
    paf.subtable <- paf.table.l[[q.name]]
    if (nrow(paf.subtable) > 0) {
      ## Flip query annotation if query.flip is TRUE
      if (unique(paf.subtable$query.flip)) {
        #annot.gr.flipped <- mirrorRanges(gr = annot.gr, seqlength = unique(paf.subtable$q.len))
        annot.gr.flipped <- annot.gr
        ranges(annot.gr.flipped) <- IRanges::reflect(x = ranges(annot.gr), bounds = IRanges::IRanges(start = 1L, end = unique(paf.subtable$q.len)))
        ## Flip ranges strand orientation
        strand(annot.gr.flipped) <- dplyr::recode(as.character(strand(annot.gr.flipped)), '+' = '-', '-' = '+')
        flipped.annot[[length(flipped.annot) + 1]] <- annot.gr.flipped 
      } else {
        flipped.annot[[length(flipped.annot) + 1]] <- annot.gr 
      }
    }
  } 
  flipped.annot.gr <- do.call(c, flipped.annot)
  return(flipped.annot.gr)
}

