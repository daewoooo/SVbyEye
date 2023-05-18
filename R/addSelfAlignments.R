#' Add PAF self-alignments to a SVbyEye miropeat style plot.
#'
#' This function takes a \code{ggplot2} object generated using \code{\link{plotMiro}} function and adds PAF self-alignments
#' stored in the `paf.table` to the plot.
#'
#' @inheritParams addAnnotation
#' @inheritParams breakPaf
#' @inheritParams plotMiro
#' @return A \code{ggplot2} object.
#' @import ggplot2
#' @importFrom ggnewscale new_scale_fill new_scale_color
#' @author David Porubsky
#' @export
#' @examples
#' ## Get PAF to process ##
#' paf.file <- system.file("extdata", "test1.paf", package = "SVbyEye")
#' ## Read in PAF
#' paf.table <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = "cg")
#' ## Make a plot
#' plt <- plotMiro(paf.table = paf.table)
#' ## Create custom PAF alignments
#' paf.annot <- tibble::tibble(q.name = 'target.region',
#'                             q.len = 0,
#'                             q.start = c(19000000, 19200000),
#'                             q.end = c(19100000, 19220000),
#'                             strand = c('+', '-'),
#'                             t.name = 'target.region',
#'                             t.len = 0,
#'                             t.start = c(19200000, 19300000),
#'                             t.end = c(19250000, 19330000),
#'                             n.match = 0,
#'                             aln.len = 0,
#'                             mapq = 0)
#' ## Add self-alignments to the plot
#' addSelfAlignments(ggplot.obj = plt, paf.table = paf.annot, coordinate.space = 'target')
#'
addSelfAlignments <- function(ggplot.obj = NULL, paf.table = NULL, color.by = "direction", color.palette = NULL, coordinate.space = "target", annotation.level = 0.05, annotation.label = NULL) {
    ## Check user input ##
    ## Make sure submitted paf.table has at least 12 mandatory fields
    if (ncol(paf.table) >= 12) {
        paf <- paf.table
    } else {
        stop("Submitted PAF alignments do not contain a minimum of 12 mandatory fields, see PAF file format definition !!!")
    }

  ## Get plotted data
  gg.data <- ggplot.obj$data

  ## Get x-axis limits (expected to be always continuous)
  ## For x-axis range also consider user defined cartesian coordinates for the target region
  if (coordinate.space == 'target') {
    xlim <- range(c(gg.data$seq.pos[gg.data$seq.id == "target"], ggplot.obj$coordinates$limits$x))
  } else {
    #xlim <- ggplot.obj$coordinates$limits$x
    xlim <- ggplot2::layer_scales(ggplot.obj)$x$range$range
  }

  ## Get x-axis limits
  if ("ScaleContinuous" %in% class(ggplot2::layer_scales(ggplot.obj)$y)) { ## To finish!!!
    ylim <- ggplot2::layer_scales(ggplot.obj)$y$range$range
    ylabels <- ggplot2::layer_scales(ggplot.obj)$y$labels
    ybreaks <- ggplot2::layer_scales(ggplot.obj)$y$breaks
    if (length(ylabels) == 0) {ylabels <- ''}
    if (length(ybreaks) == 0) {ybreaks <- 0}
    #ylabels.ord <- ggplot2::layer_scales(ggplot.obj)$y$breaks
    ylabels <- ylabels[order(ybreaks)]
    ybreaks <- sort(ybreaks) ## [Testing]
  } else {
    stop("'addSelfAnnotation' function works only for ggplot2 objects with continuous scale for both x and y-axis !!!")
  }

  ## Define the offset value to be the user defined fraction of the y-axis range [default: 0.05]
  offset <- diff(ylim) * annotation.level

  ## Get query and target coordinate ranges
  if (coordinate.space == "query" | coordinate.space == "target") {
    t.range <- range(gg.data$seq.pos[gg.data$seq.id == "target"])
    q.range <- range(gg.data$seq.pos[gg.data$seq.id == "query"])
    ## Adjust target ranges given the size difference with respect to query ranges
    range.offset <- diff(q.range) - diff(t.range)
    t.range[2] <- t.range[2] + range.offset ## Make a start position as offset and change only end position
  }

  ## Translate query coordinates to target coordinates
  if (coordinate.space == "query") {
    ## Covert query to target coordinates
    paf$q.start <- SVbyEye::q2t(x = paf$q.start, q.range = q.range, t.range = t.range)
    paf$q.end <- SVbyEye::q2t(x = paf$q.end, q.range = q.range, t.range = t.range)
    paf$t.start <- SVbyEye::q2t(x = paf$t.start, q.range = q.range, t.range = t.range)
    paf$t.end <- SVbyEye::q2t(x = paf$t.end, q.range = q.range, t.range = t.range)
  }

  ## Get annotation track offset
  if (coordinate.space == "target") {
    y.offset <- max(ylim) + offset
  } else if (coordinate.space == "query") {
    y.offset <- min(ylim) + -offset
  } else {
    stop("Please specify 'coordinate.space' as either 'target' or 'query' !!!")
  }

  ## Restrict selfalignments to x-limits
  mask <- paf$q.start >= xlim[1] & paf$t.start >= xlim[1] & paf$q.end <= xlim[2] & paf$t.end <= xlim[2]
  paf <- paf[mask,]
  if (nrow(paf) > 0) {
    ## Prepare data for plotting ##
    ## Flip start and end for reverse oriented alignments
    paf[paf$strand == "-", ] <- transform(paf[paf$strand == "-", ],
                                             "t.start" = paf[paf$strand == "-", ]$t.end,
                                             "t.end" = paf[paf$strand == "-", ]$t.start
    )

    ## Define color palette for alignment directionality
    if (!is.null(color.palette)) {
      if (all(c("+", "-") %in% names(color.palette))) {
        pal <- color.palette
      } else {
        pal <- c("-" = "cornflowerblue", "+" = "forestgreen")
        warning("User defined 'color.palette' does not contain both '+' and '-' directions, using default values instead!!!")
      }
    } else {
      pal <- c("-" = "cornflowerblue", "+" = "forestgreen")
    }

    ## Define discrete color palettes
    if (color.by == "direction") {
      paf$col.levels <- factor(paf$strand, levels = c("+", "-"))
    } else if (color.by == "identity") {
      paf$identity <- (paf$n.match / paf$aln.len) * 100
      paf$identity[is.nan(paf$identity) | is.na(paf$identity)] <- 0
      ## Define color scheme
      colors.l <- getColorScheme(data.table = paf, value.field = "identity", breaks = c(90, 95, 99, 99.5, 99.6, 99.7, 99.8, 99.9))
      paf <- colors.l$data
      #paf$identity <- paf$col.levels
      pal <- colors.l$colors
    } else if (color.by %in% colnames(paf)) {
      ## Define color scheme
      colors.l <- getColorScheme(data.table = paf, value.field = color.by)
      paf <- colors.l$data
      pal <- colors.l$colors
    } else {
      color.by <- "direction"
    }

    ## Prepare data for plotting
    x <- c(rbind(paf$q.start, paf$t.start, paf$q.end, paf$t.end))
    group <- rep(seq_len(nrow(paf)), each = 4)
    col.levels <- rep(paf$col.levels, each = 4)

    plt.df <- data.frame(
      x = x,
      group = group,
      #seq.id = seq.id,
      col.levels
    )

    ## Make a plot
    if (coordinate.space == "target") {
      plt <- ggplot.obj + ggnewscale::new_scale_fill() +
        geom_wide_arc(data = plt.df,
                      ggplot2::aes(x = .data$x, group = .data$group, fill = .data$col.levels),
                      alpha = 0.5, y.offset = y.offset, max.width = 0.5)
    } else if (coordinate.space == "query") {
      plt <- ggplot.obj + ggnewscale::new_scale_fill() +
        geom_wide_arc(data = plt.df,
                      ggplot2::aes(x = .data$x, group = .data$group, fill = .data$col.levels),
                      alpha = 0.5, y.offset = y.offset, max.width = 0.5, y.reverse = TRUE)
    }

    ## Add color scale
    if (color.by == "direction") {
      plt <- plt + scale_fill_manual(values = pal, name = 'Self-alignment\ndirection', drop = FALSE)
    } else if (color.by == "identity") {
      plt <- plt + scale_fill_manual(values = pal, name = 'Identity', drop = FALSE)
    }

    ## Add y-label to annotation track if defined
    if (!is.null(annotation.label)) {
      if (nchar(annotation.label) > 0) {
        annot.break <- y.offset
        if (coordinate.space == "target") {
          y.breaks <- c(ybreaks, annot.break)
          y.labels <- c(ylabels, annotation.label)
        } else {
          y.breaks <- c(annot.break, ybreaks)
          y.labels <- c(annotation.label, ylabels)
        }
        suppressMessages(
          plt <- plt + ggplot2::scale_y_continuous(breaks = y.breaks, labels = y.labels) +
            ggplot2::theme(axis.text.y = ggplot2::element_text(),
                           axis.ticks.y = ggplot2::element_line())
        )
      }
    }

    ## Return updated plot
    return(plt)
  } else {
    ## Return original plot
    return(ggplot.obj)
  }
}
