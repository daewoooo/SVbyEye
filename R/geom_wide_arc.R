#' A 'ggplot2' geom to draw arcs between genomic alignments.
#'
#' `geom_wide_arc()` draws wide polygons between two sets of start and end coordinates.
#'
#' This geom is intended to draws wide arc polygons between self-alignments defined in PAF format.
#' Such alignments can be directly visualized using a wrapper function \code{\link{plotSelf}}.
#' Input data is a \code{data.frame} object that contains required aesthetics `x` and `y` coordinates as well as
#' `group` field that is required in order to determine which coordinates represent a single alignment.
#'
#' @section Aesthetics:
#' `geom_wide_arc()` require or can take the following aesthetics
#' (required aesthetics are in bold):
#'
#' - **x**
#' - **group**
#' - y (if not defined it will set to zero)
#' - color
#' - linewidth
#' - linetype
#' - alpha
#' - fill
#' - size
#'
#' @inheritParams ggplot2::geom_polygon
#' @inheritParams ggplot2::stat_identity
#'
#' @param mapping,data,stat,position,na.rm,show.legend,inherit.aes,... As is standard for ggplot2.
#' @param n The number of points to create for each alignment polygon (Default : `100`).
#' @param max.width The maximum width of the arc in y-coordinate units.
#' @param y.offset The y-axis coordinate from where the arc will start.
#' @param y.reverse Set to `TRUE` if the arc orientation should be flipped upside down.
#'
#' @seealso [plotSelf()]
#' @return Plotting coordinates
#' @author David Porubsky
#'
#' @name geom_wide_arc
#' @rdname geom_wide_arc
#'
#' @examples
#' ## Create example data.frame to plot ##
#' ## Each link between sequence region1 (seq1) and sequence region2 (seq2)
#' ## is expected to have 4 x-coordinates (start.seq1, start.seq2, end.seq1, end.seq2)
#' plt.df <- data.frame(
#'     x = c(100, 500, 200, 1000),
#'     group = 1
#' )
#' ## Make a plot
#' ggplot2::ggplot(plt.df) +
#'     geom_wide_arc(ggplot2::aes(x = x, group = group), alpha = 0.5)
#'
NULL

#' A ggproto class definition to extend the functionality of ggplot2.
#' @importFrom ggforce StatBezier
#' @format NULL
#' @usage NULL
#' @return Plotting coordinates
#' @export
StatWideArc <- ggplot2::ggproto("StatWideArc", ggplot2::Stat,
    setup_data = function(data, params) {
        if (any(table(data$group) != 4)) {
            stop("Each group must consist of 4 points")
        }
        data
    },
    compute_panel = function(data, scales, n = 100, max.width = NULL, y.offset = 0, y.reverse = FALSE) {
        if (!"y" %in% colnames(data)) {
            data$y <- 0
        }
        # data <- data[order(data$group, data$y), ]
        data <- data[order(data$group), ]
        coords1 <- data[c(TRUE, FALSE, FALSE, TRUE), ]
        coords2 <- data[c(FALSE, TRUE, TRUE, FALSE), ]
        coords1 <- add.widearc.controls.points(coords1)
        coords2 <- add.widearc.controls.points(coords2)
        ## Scale y height to user defined max.width
        if (!is.null(max.width)) {
            if (max.width > 0) {
                max.y <- max(c(coords1$y, coords2$y))
                coords1$y <- q2t(coords1$y, q.range = c(0, max.y), t.range = c(0, max.width))
                coords2$y <- q2t(coords2$y, q.range = c(0, max.y), t.range = c(0, max.width))
            }
        }
        ## Reverse arc upside down
        if (y.reverse) {
            coords1$y <- coords1$y * -1
            coords2$y <- coords2$y * -1
        }
        ## Add offset to y axis coordinates
        coords1$y <- coords1$y + y.offset
        coords2$y <- coords2$y + y.offset
        ## Calculate bezier curve
        coords1 <- ggforce::StatBezier$compute_panel(coords1, scales, n)
        coords2 <- ggforce::StatBezier$compute_panel(coords2, scales, n)
        arc.coords <- rbind(coords1, coords2)
    },
    required_aes = c("x", "group"),
    extra_params = c("na.rm", "n", "max.width", "y.offset")
)

#' @rdname geom_wide_arc
#' @export
geom_wide_arc <- function(mapping = NULL, data = NULL, geom = "polygon",
                          stat = "wide_arc", position = "identity",
                          n = 100, max.width = NULL, y.offset = 0,
                          y.reverse = FALSE, na.rm = FALSE,
                          show.legend = NA, inherit.aes = TRUE, ...) {
    ggplot2::layer(
        data = data, mapping = mapping, stat = stat, geom = geom,
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(
            na.rm = na.rm, n = n, max.width = max.width,
            y.offset = y.offset, y.reverse = y.reverse, ...
        )
    )
}

## Helper function definition
#' @format NULL
#' @usage NULL
add.widearc.controls.points <- function(data) {
    start <- data[c(TRUE, FALSE), ]
    end <- data[c(FALSE, TRUE), ]
    y_height <- sqrt(abs(end$x - start$x))
    mid1 <- start
    mid1$y <- y_height
    mid2 <- end
    mid2$y <- y_height
    rbind(start, mid1, mid2, end)
}
