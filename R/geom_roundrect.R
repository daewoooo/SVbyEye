#' A 'ggplot2' geom to draw genomic annotations as round rectangles
#'
#' `geom_roundrect()` draws genomic ranges as round rectangles, allowing to draw for instance ...
#'
#' This geom draws ...
#'
#' @export
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

# grid::roundrectGrob(
#   1, 2,
#   width = 1,
#   height = 1 * 2,
#   r = 1
# )  

# gr <- GRanges(seqnames=paste0('chr', c(1:3)), ranges = IRanges(start=c(10, 100, 200), end=c(100, 190, 400)))
# df <- as.data.frame(gr)
# 
# ggplot(df) +
#   geom_roundrect(aes(xmin=start, xmax=end, y=1, fill=seqnames)) +
#   facet_grid(seqnames ~ .)
# 
# ggplot(df) +
#   geom_roundrect(aes(xmin=start, xmax=end, y=1, fill=seqnames), rect_height = unit(5, 'mm'), radius = unit(2, 'mm'))
# 
# ggplot(df) +
#   geom_roundrect(aes(xmin=start, xmax=end, y=seqnames, fill=seqnames))
