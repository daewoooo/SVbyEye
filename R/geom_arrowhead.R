#' A 'ggplot2' geom to draw genomic ranges as arrowheads.
#'
#' `geom_arrowhead()` draws ranges defined by `xmin` and `xmax` as triangular polygon.
#' draws genomic ranges as arrowheads, allowing to draw for instance segmental
#' duplication maps.
#'
#' This geom draws triangular polygons as arrowheads between defined start and end coordinates. 
#' Intended application of this geom is to visualize genomic coordinates defined by start and end
#' position.
#' 
#' @section Aesthetics:
#' `geom_roundrect()` require or can take the following aesthetics
#' (required aesthetics are in bold):
#'
#' - **xmin**
#' - **xmax**
#' - **y**
#' - color
#' - linewidth
#' - linetype
#' - alpha
#' - fill
#' - size
#' 
#' @param mapping,data,stat,position,na.rm,show.legend,inherit.aes, etc... As standard for ggplot2.
#' @param arrowhead_height A `grid::unit()` object providing the height of the arrowhead.  [Default: 3 mm].
#' @author David Porubsky
#' @export
#' @examples 
#'## Create example data.frame to plot
#'plt.df <- data.frame(xmin=c(10, 100, 200),
#'                     xmax=c(100, 190, 400)
#'                     )
#'## Plot rectangles with rounded edges
#'ggplot2::ggplot(plt.df) + 
#'  geom_arrowhead(ggplot2::aes(xmin=xmin, xmax=xmax, y=1))
#'
geom_arrowhead <- function(
  mapping = NULL,
  data = NULL,
  stat = "identity",
  position = "identity",
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE,
  arrowhead_height = grid::unit(3, "mm"),
  ...
) {
  ggplot2::layer(
    geom = GeomArrowhead, mapping = mapping, data = data, stat = stat,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      arrowhead_height = arrowhead_height,
      ...
    )
  )
}

#' GeomArrowhead
#' @noRd
GeomArrowhead <- ggplot2::ggproto("GeomArrowhead", ggplot2::Geom,
                                  required_aes = c("xmin", "xmax", "y"),
                                  default_aes = ggplot2::aes(
                                    strand = TRUE,
                                    alpha = 1,
                                    colour = "black",
                                    fill = "white",
                                    linetype = 1,
                                    size = 0.3
                                  ),
                                  draw_key = function(data, params, size) {
                                    grid::rectGrob(
                                      width = grid::unit(1, "npc") - grid::unit(1, "mm"),
                                      height = grid::unit(1, "npc") - grid::unit(1, "mm"),
                                      gp = grid::gpar(
                                        col = data$colour,
                                        fill = ggplot2::alpha(data$fill, data$alpha),
                                        lty = data$linetype,
                                        lwd = data$size * ggplot2::.pt
                                      )
                                    )
                                  },
                                  draw_panel = function(
                                    data,
                                    panel_scales,
                                    coord,
                                    arrowhead_height
                                  ) {
                                    
                                    data <- coord$transform(data, panel_scales)
                                    #str(data)
                                    
                                    gt <- grid::gTree(
                                      data = data,
                                      cl = "arrowheadtree",
                                      arrowhead_height = arrowhead_height
                                    )
                                    gt$name <- grid::grobName(gt, "geom_arrowhead")
                                    gt
                                  }
)

#' @importFrom grid makeContent
#' @export
makeContent.arrowheadtree <- function(x) {
  
  data <- x$data
  
  # Prepare grob for genomic range
  grobs <- lapply(1:nrow(data), function(i) {
    
    range <- data[i, ]
    
    # Reverse non-forward genes
    if (range$strand == '-') {
      range[, c("xmin", "xmax")] <- range[, c("xmax", "xmin")]
    }
    
    # Set arrowhead heights. Divide by 2 for convenience to calculate polygon coordinates
    arrowhead_height <- as.numeric(grid::convertHeight(x$arrowhead_height, "native")) / 2
    
    # Create polygon grob
    pg <- grid::polygonGrob(
      x = c(
        range$xmin,
        range$xmax,
        range$xmin
      ),
      y = c(
        range$y - arrowhead_height,
        range$y,
        range$y + arrowhead_height
      ),
      gp = grid::gpar(
        fill = ggplot2::alpha(range$fill, range$alpha),
        col = ggplot2::alpha(range$colour, range$alpha),
        lty = range$linetype,
        lwd = range$size * ggplot2::.pt
      )
    )
    
    # Return the polygon grob
    pg
  })
  
  class(grobs) <- "gList"
  grid::setChildren(x, grobs)
}
