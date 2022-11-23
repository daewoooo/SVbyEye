#' A 'ggplot2' geom to draw genomic ranges as round rectangles.
#'
#' `geom_roundrect()` draws ranges defined by `xmin` and `xmax` coordinates with rounded edges.
#'
#' This geom draws rectangle with round or sharp edges between defined start and end coordinates. 
#' Intended application of this geom is to visualize genomic coordinates defined by start and end
#' position. Rounded edges will help to observe boundaries between closely positioned genomic ranges.
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
#' @param rect_height A `grid::unit()` object providing the height of the rectangle.  [Default: 3 mm].
#' @param radius A `grid::unit()` object providing required curvature of rectangle edges. [Default: 1 mm].
#' @author David Porubsky
#' @export
#' @examples 
#'## Create example data.frame to plot
#'plt.df <- data.frame(xmin=c(10, 100, 200),
#'                     xmax=c(100, 190, 400)
#'                     )
#'## Plot rectangles with rounded edges
#'ggplot2::ggplot(plt.df) + 
#'  geom_roundrect(ggplot2::aes(xmin=xmin, xmax=xmax, y=1))
#'## Plot rectangles without rounded edges
#'ggplot2::ggplot(plt.df) + 
#'  geom_roundrect(ggplot2::aes(xmin=xmin, xmax=xmax, y=1), radius=grid::unit(0, "mm"))
#'
geom_roundrect <- function(
  mapping = NULL,
  data = NULL,
  stat = "identity",
  position = "identity",
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE,
  rect_height = grid::unit(3, "mm"),
  radius = grid::unit(1, "mm"),
  ...
) {
  ggplot2::layer(
    geom = GeomRoundRect, mapping = mapping, data = data, stat = stat,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      rect_height = rect_height,
      radius = radius,
      ...
    )
  )
}

#' GeomRoundRect
#' @noRd
GeomRoundRect <- ggplot2::ggproto("GeomRoundRect", ggplot2::Geom,
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
                                    rect_height,
                                    radius
                                  ) {
                                    
                                    data <- coord$transform(data, panel_scales)
                                    #str(data)
                                    #print(data)
                                    
                                    gt <- grid::gTree(
                                      data = data,
                                      cl = "roundrecttree",
                                      rect_height = rect_height,
                                      radius = radius
                                    )
                                    gt$name <- grid::grobName(gt, "geom_roundrect")
                                    gt
                                  }
)

#' @importFrom grid makeContent
#' @export
makeContent.roundrecttree <- function(x) {
  
  data <- x$data
  
  # Prepare grob for genomic range
  grobs <- lapply(1:nrow(data), function(i) {
    
    range <- data[i, ]
    
    # Set arrowhead heights. Divide by 2 for convenience to calculate polygon coordinates
    rect_height <- as.numeric(grid::convertHeight(x$rect_height, "native")) / 2
    radius <- x$radius
    
    # Create polygon grob
    pg <- grid::roundrectGrob(
      range$xmin, range$y + rect_height,
      width = (range$xmax - range$xmin),
      height = rect_height * 2,
      r = radius,
      default.units = "native",
      just = c("left", "top"),
      gp = grid::gpar(
        fill = ggplot2::alpha(range$fill, range$alpha),
        col = ggplot2::alpha(range$colour, range$alpha),
        lty = range$linetype,
        lwd = range$size * ggplot2::.pt,
        lineend = "butt"
      )
    )
    # Return the roundrect grob
    pg
  })
  
  class(grobs) <- "gList"
  grid::setChildren(x, grobs)
}
