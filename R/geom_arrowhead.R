#' A 'ggplot2' geom to draw genomic annotations as arrowheads
#'
#' `geom_arrowhead()` draws genomic ranges as arrowheads, allowing to draw for instance segmental
#' duplication maps.
#'
#' This geom draws genes as arrows along a horizontal line representing the
#' molecule. The start and end locations of the gene are expressed with the
#' `xmin` and `xmax` aesthetics, while the molecule can be specified with the
#' `y` aesthetic. Optionally, an additional `forward` aesthetic can be used to
#' reverse the orientation of some or all genes from that implied by `xmin` and
#' `xmax`.
#'
#' Unless the plot is faceted with a free x scale, all the molecules will share
#' a common x axis. This means that if the locations are very different across
#' different molecules, the genes might appear very small and squished together
#' with a lot of unnecessary empty space. To get around this, either facet the
#' plot with `scales = "free_x"`, or normalise the gene locations if their
#' exact locations are not important.
#'
#' See `make_alignment_dummies()` for a method to align genes between molecules.
#'
#' @section Aesthetics:
#'
#' - xmin,xmax (start and end of the gene; will be used to determine gene
#' orientation)
#' - y (molecule)
#' - forward (if any value that is not TRUE, or coercible to TRUE, the gene
#' arrow will be drawn in the opposite direction to that determined by `xmin`
#' and `xmax`)
#' - alpha
#' - colour
#' - fill
#' - linetype
#' - size
#'
#' @param mapping,data,stat,position,na.rm,show.legend,inherit.aes,... As
#' standard for ggplot2.
#' @param arrowhead_height `grid::unit()` object giving the height of the
#' arrowhead.  Defaults to 3 mm.
#'
#' @examples
#'
#' ggplot2::ggplot(example_genes, ggplot2::aes(xmin = start, xmax = end,
#'                                             y = molecule, fill = gene)) +
#' geom_arrowhead() +
#' ggplot2::facet_wrap(~ molecule, scales = "free")
#'
#' @export
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
