#' Add PAF alignments to a SVbyEye miropeat style plot.
#'
#' This function takes a \code{ggplot2} object generated using \code{\link{plotMiro}} function and adds extra PAF alignments to it
#' stored in the `paf.table`. This function can also be used to highlight already present alignment or to add other features such as
#' alignment gaps reported by the \code{\link{paf2gaps}} function.
#'
#' @param color.by A name of a column present in submitted `paf.table`.to be used to define color scheme.
#' @param fill.by A name of a column present in submitted `paf.table`.to be used to define fill scheme.
#' @param color.palette A discrete color palette defined as named character vector (elements = colors, names = discrete levels).
#' @param fill.palette Same as parameter `color.palette`.
#' @param linetype One of the following types ['solid', 'dashed', 'dotted', 'dotdash', 'longdash', 'twodash'].
#' @param coord.strict If set to \code{TRUE} PAF alignments from the `paf.table` will be subsetted in order to be within
#' x-axis limits of the `ggplot.obj`.
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
#' ## Report gaps between PAF alignments
#' paf.gaps <- paf2gaps(paf.table = paf.table)
#' ## Add annotation field for plotting
#' paf.gaps$SV <- "INS"
#' paf.gaps$SV[grep(paf.gaps$cg, pattern = "D", ignore.case = TRUE)] <- "DEL"
#' ## Define color palette for plotting
#' color.palette <- c("DEL" = "firebrick3", "INS" = "dodgerblue3")
#' ## Add gap position to the plot as extra PAF alignmetns
#' addAlignments(
#'     ggplot.obj = plt, paf.table = paf.gaps,
#'     color.by = "SV", color.palette = color.palette, linetype = "dashed"
#' )
#'
addAlignments <- function(ggplot.obj = NULL, paf.table = NULL, color.by = NULL, color.palette = NULL, fill.by = NULL, fill.palette = NULL, outline.alignments = FALSE, linetype = "solid", coord.strict = TRUE) {
    ## Check user input ##
    ## Make sure submitted paf.table has at least 12 mandatory fields
    if (ncol(paf.table) >= 12) {
        paf <- paf.table
    } else {
        stop("Submitted PAF alignments do not contain a minimum of 12 mandatory fields, see PAF file format definition !!!")
    }

    ## Get plotted data
    gg.data <- ggplot.obj$data
    target.id <- unique(gg.data$seq.name[gg.data$seq.id == "target"])
    query.id <- unique(gg.data$seq.name[gg.data$seq.id == "query"])

    ## Get x and y-axis limits
    # xlim <- ggplot2::layer_scales(ggplot.obj)$x$range$range
    ## For x-axis range also consider user defined cartesian coordinates for target region
    xlim <- range(c(gg.data$seq.pos[gg.data$seq.id == "target"], ggplot.obj$coordinates$limits$x))
    ylim <- ggplot2::layer_scales(ggplot.obj)$y$range$range
    ylabels <- ggplot2::layer_scales(ggplot.obj)$y$labels

    ## Get query and target coordinate ranges
    t.range <- range(gg.data$seq.pos[gg.data$seq.id == "target"])
    q.range <- range(gg.data$seq.pos[gg.data$seq.id == "query"])
    ## Adjust target ranges given the size difference with respect to query ranges
    range.offset <- diff(q.range) - diff(t.range)
    t.range[2] <- t.range[2] + range.offset ## Make a start position as offset and change only end position

    ## Make sure that submitted PAF alignments are within ggplot.obj coordinates
    if (coord.strict) {
        q.filt <- paf$q.start >= q.range[1] & paf$q.end <= q.range[2]
        t.filt <- paf$t.start >= t.range[1] & paf$t.end <= t.range[2]
        paf.table <- paf.table[q.filt & t.filt, ]
    }

    ## Convert PAF alignments to plotting coordinates
    if (color.by %in% colnames(paf)) {
        coords <- paf2coords(paf.table = paf, add.col = color.by)
    } else {
        coords <- paf2coords(paf.table = paf)
    }

    ## Define linetype given the user input
    if (!is.null(linetype)) {
        allowed.linetypes <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
        if (!linetype %in% allowed.linetypes) {
            w.msg <- c(
                "Parameter linetype can only take values [",
                paste(shQuote(allowed.linetypes), collapse = ", "),
                "], using default 'solid' linetype !!!"
            )
            warning(paste0(w.msg))
            linetype <- "solid"
        }
    }

    ## Define color or fill of the added alignments
    if (!is.null(color.by)) {
        if (!color.by %in% colnames(paf)) {
            color.by <- NULL
        } else {
            ## If color.by is defined it has precedence over fill.by
            fill.by <- NULL
        }
    }
    if (!is.null(fill.by)) {
        if (!fill.by %in% colnames(paf)) {
            fill.by <- NULL
        }
    }

    ## Add alignments to the plot
    if (!is.null(color.by)) {
        ggplot.obj <- ggplot.obj + ggnewscale::new_scale_color() +
            geom_miropeats(
                data = coords, ggplot2::aes(x = .data$x, y = .data$y, group = .data$group, color = .data[[color.by]]),
                fill = NA, linetype = linetype, alpha = 0.5, inherit.aes = FALSE
            ) +
            ggplot2::scale_color_manual(values = color.palette, drop = FALSE)
    } else if (!is.null(fill.by)) {
        ggplot.obj <- ggplot.obj + ggnewscale::new_scale_fill() +
            geom_miropeats(
                data = coords, ggplot2::aes(x = .data$x, y = .data$y, group = .data$group, fill = .data[[fill.by]]),
                alpha = 0.5, inherit.aes = FALSE
            ) +
            ggplot2::scale_fill_manual(values = fill.palette, drop = FALSE)
        ## Add alignment outlines
        if (outline.alignments) {
            ggplot.obj <- ggplot.obj +
                geom_miropeats(data = coords, ggplot2::aes(x = .data$x, y = .data$y, group = .data$group), fill = NA, color = "gray", linewidth = 0.25)
        }
    } else {
        ggplot.obj <- ggplot.obj +
            geom_miropeats(
                data = coords, ggplot2::aes(x = .data$x, y = .data$y, group = .data$group),
                fill = NA, color = "black", linetype = linetype, alpha = 0.5, inherit.aes = FALSE
            )
    }
    ## Return updated ggplot object
    return(ggplot.obj)
}
