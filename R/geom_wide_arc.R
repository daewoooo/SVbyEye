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
#' - **y**
#' - **group**
#' - color
#' - linewidth
#' - linetype
#' - alpha
#' - fill
#' - size
#' 
#' @param mapping,data,stat,position,na.rm,show.legend,inherit.aes, etc... As standard for ggplot2.
#' @seealso [plotSelf()]
#' @author David Porubsky
#' @export
#' @examples 
#'## Create example data.frame to plot
#'plt.df <- data.frame(x=c(100, 500, 200, 1000),
#'                     y=0,
#'                     group=1)
#'## Make a plot 
#'ggplot2::ggplot(plt.df) + 
#'  geom_wide_arc(ggplot2::aes(x=x, y=y, group=group), alpha=0.5)
#'
geom_wide_arc <- function(mapping = NULL, data = NULL, geom = 'polygon',
                           stat = 'wide_arc', position = 'identity',
                           n = 100, na.rm = FALSE,
                           show.legend = NA, inherit.aes = TRUE, ...) {
  ggplot2::layer(
    data = data, mapping = mapping, stat = stat, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, n = n, ...)
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

#' @importFrom ggforce StatBezier
#' @format NULL
#' @usage NULL
#' @export
StatWideArc <- ggplot2::ggproto('StatWideArc', ggplot2::Stat,
                                setup_data = function(data, params) {
                                  if (any(table(data$group) != 4)) {
                                    stop('Each group must consist of 4 points')
                                  }
                                  data
                                },
                                compute_panel = function(data, scales, n = 100) {
                                  data <- data[order(data$group, data$y), ]
                                  coords1 <- data[c(TRUE, FALSE, FALSE, TRUE), ]
                                  coords2 <- data[c(FALSE, TRUE, TRUE, FALSE), ]
                                  coords1 <- add.widearc.controls.points(coords1)
                                  coords2 <- add.widearc.controls.points(coords2)
                                  coords1 <- ggforce::StatBezier$compute_panel(coords1, scales, n)
                                  coords2 <- ggforce::StatBezier$compute_panel(coords2, scales, n)
                                  arc.coords <- rbind(coords1, coords2)
                                },
                                required_aes = c('x', 'y', 'group'),
                                extra_params = c('na.rm', 'n')
)