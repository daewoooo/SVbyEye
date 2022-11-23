#' Add annotation ranges to miropeat style plot.
#'
#' This function takes a \code{ggplot2} object generated using \code{\link{plotMiro}} function and adds extra annotation on top of query 
#' or target coordinates. These ranges are specified in 'annot.gr' object and are visualized either as arrowheads or rectangles.
#'
#' @param ggplot.obj A \code{ggplot2} object generated using \code{\link{plotMiro}} function.
#' @param annot.gr A \code{\link{GRanges-class}} object with a set of ranges to be added as extra annotation.
#' @param shape A user defined shape ranges in 'annot.gr' are visualized, either 'arrowhead' or 'rectangle'.
#' @param fill.by A name of an extra field present in 'annot.gr' to be used to define color scheme.
#' @param color.palette A discrete color palette defined as named character vector (elements = colors, names = discrete levels).  
#' @param coordinate.space A coordinate space ranges in 'annot.gr' are reported, either 'target', 'query' or 'self'.
#' @param offset.annotation Set to \code{TRUE} if subsequent annotation ranges should be offsetted below and above the midline.
#' @param annotation.label A \code{character} string to be used as a label to added annotation track.
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
#'plt <- addAnnotation(ggplot.obj = plt, annot.gr = target.annot.gr, coordinate.space = 'target')
#'## Load query annotation file
#'query.annot <- system.file("extdata", "test1_query_annot.txt", package="SVbyEye")
#'query.annot.df <- read.table(query.annot, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
#'query.annot.gr <- makeGRangesFromDataFrame(query.annot.df)
#'## Add query annotation as rectangle
#'addAnnotation(ggplot.obj = plt, annot.gr = query.annot.gr, shape = 'rectangle', coordinate.space = 'query')
#'## Lift target annotation to query and plot
#'lifted.annot.gr <- liftRangesToAlignment(gr = target.annot.gr, paf.file = paf.file, direction = 'target2query')
#'addAnnotation(ggplot.obj = plt, annot.gr = lifted.annot.gr, shape = 'rectangle', coordinate.space = 'query')
#'## Add segmental duplication annotation
#'plt <- plotMiro(paf.table = paf.table)
#'sd.annot <- system.file("extdata", "test1.sd.annot.RData", package="SVbyEye")
#'sd.annot.gr <- get(load(sd.annot))
#'## Create a custom discrete levels
#'sd.categ <- findInterval(sd.annot.gr$fracMatch, vec = c(0.95, 0.98, 0.99))
#'sd.categ <- dplyr::recode(sd.categ, '0' = '<95%', '1' = '95-98%', '2' = '98-99%', '3'='>=99%')
#'sd.categ <- factor(sd.categ, levels=c('<95%', '95-98%', '98-99%', '>=99%'))
#'sd.annot.gr$sd.categ <- sd.categ
#'## Create a custom color palette
#'color.palette <- c('<95%' = 'gray72', '95-98%' = 'gray47', '98-99%' = '#cccc00', '>=99%' = '#ff6700')
#'## Add annotation to the plot
#'addAnnotation(ggplot.obj = plt, annot.gr = sd.annot.gr, fill.by = 'sd.categ', color.palette = color.palette, coordinate.space = 'target')
#'## Offset annotation ranges
#'addAnnotation(ggplot.obj = plt, annot.gr = sd.annot.gr, fill.by = 'sd.categ', color.palette = color.palette, coordinate.space = 'target', offset.annotation = TRUE)
#'## Add label to the added annotation
#'addAnnotation(ggplot.obj = plt, annot.gr = sd.annot.gr, fill.by = 'sd.categ', color.palette = color.palette, coordinate.space = 'target', annotation.label = 'SD')
#'
addAnnotation <- function(ggplot.obj=NULL, annot.gr=NULL, shape='arrowhead', fill.by=NULL, color.palette=NULL, coordinate.space='target', new.annotation.level=TRUE, offset.annotation=FALSE, annotation.label=NULL) {
  ## Get plotted data
  gg.data <- ggplot.obj$data
  target.id <- unique(gg.data$seq.name[gg.data$seq.id == 'target'])
  query.id <- unique(gg.data$seq.name[gg.data$seq.id == 'query'])
  
  ## Get x and y-axis limits
  #xlim <- ggplot2::layer_scales(ggplot.obj)$x$range$range
  ## For x-axis range also consider user defined cartesian coordinates for target region
  xlim <- range(c(gg.data$seq.pos[gg.data$seq.id == 'target'], ggplot.obj$coordinates$limits$x))
  ylim <- ggplot2::layer_scales(ggplot.obj)$y$range$range
  ylabels <- ggplot2::layer_scales(ggplot.obj)$y$labels
  
  ## Define the offset value to be 5% of the y axis range
  if (new.annotation.level) {
    offset <- diff(ylim) * 0.05
  } else {
    offset <- 0
  }  
  
  ## Get query and target coordinate ranges
  if (coordinate.space == 'query' | coordinate.space == 'target') {
    t.range <- range(gg.data$seq.pos[gg.data$seq.id == 'target'])
    q.range <- range(gg.data$seq.pos[gg.data$seq.id == 'query'])
    ## Adjust target ranges given the size difference with respect to query ranges
    range.offset <- diff(q.range) - diff(t.range)
    t.range[2] <- t.range[2] + range.offset ## Make a start position as offset and change only end position
  }  
  
  ## Translate query coordinates to target coordinates
  if (coordinate.space == 'query') {
    if (length(annot.gr) > 0) {
      ## Covert query to target coordinates
      new.start <- SVbyEye::q2t(x = start(annot.gr), q.range = q.range, t.range = t.range)
      new.end <- SVbyEye::q2t(x = end(annot.gr), q.range = q.range, t.range = t.range)
      new.ranges <- IRanges::IRanges(start = new.start, end = new.end)
      suppressWarnings( ranges(annot.gr) <- new.ranges )
    }  
  }  
  
  ## Get annotation track offset
  if (coordinate.space == 'target') {
    y.offset <- max(ylim) + offset
  } else if (coordinate.space == 'query') {
    y.offset <- min(ylim) + -offset
  } else if (coordinate.space == 'self') {
    y.offset <- min(ylim) + -offset
  } else {
    stop("Please specify 'coordinate.space' as either 'target', 'query' or 'self' !!!")
  }
  
  if (!is.null(annot.gr) & class(annot.gr) == 'GRanges') {
    ## Restrict to plot x-limits
    annot.gr <- annot.gr[start(annot.gr) >= xlim[1] & end(annot.gr) <= xlim[2]]
    if (length(annot.gr) > 0) {
      ## Offset overlapping annotation ranges up&down based on start position
      if (offset.annotation) {
        annot.gr <- GenomicRanges::sort(annot.gr, ignore.strand=TRUE)
        if (coordinate.space == 'target') {
          y.offset <- rep(y.offset + c(0, offset), times=ceiling(length(annot.gr) / 2))[1:length(annot.gr)]
        } else if (coordinate.space == 'query') {
          y.offset <- rep(y.offset - c(0, offset), times=ceiling(length(annot.gr) / 2))[1:length(annot.gr)]
        }  
      }
      ## Convert to data frame object
      annot.df <- as.data.frame(annot.gr)
      ## Add offset column
      annot.df$y.offset <- y.offset
      
      ## Define color scale ##
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
              if (length(pal) >= 20) {
                col.scale <- 'discreteNoLegend'
              } else {
                col.scale <- 'discrete'
              }  
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
      
      ## Add annotation to the ggplot obkect ##
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
      
      ## Add y-label to annotation track if defined
      if (!is.null(annotation.label)) {
        if (nchar(annotation.label) > 0) {
          annot.yrange <- range(annot.df$y.offset)
          annot.break <- annot.yrange[1] + (diff(annot.yrange) / 2)
          y.breaks <- c(ylim, annot.break)
          y.labels <- c(ylabels, annotation.label)
          suppressMessages(
            plt <- plt + scale_y_continuous(breaks = y.breaks, labels = y.labels)
          )  
        }
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
        suppressWarnings(
          ranges(annot.gr.flipped) <- IRanges::reflect(x = ranges(annot.gr), bounds = IRanges::IRanges(start = 1L, end = unique(paf.subtable$q.len)))
        )
        ## Flip ranges strand orientation
        strand(annot.gr.flipped) <- dplyr::recode(as.character(strand(annot.gr.flipped)), '+' = '-', '-' = '+')
        flipped.annot[[length(flipped.annot) + 1]] <- annot.gr.flipped 
      } else {
        flipped.annot[[length(flipped.annot) + 1]] <- annot.gr 
      }
    }
  } 
  flipped.annot.gr <- suppressWarnings( do.call(c, flipped.annot) )
  return(flipped.annot.gr)
}