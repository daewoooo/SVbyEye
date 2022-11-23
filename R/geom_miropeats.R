#' A 'ggplot2' geom to draw genomic alignments in a miropeats style.
#' 
#' `geom_miropeats()` draws miropeat style polygons between query to target genomic alignments.
#' 
#' This geom draws polygons between query to target alignments defined in PAF format.
#' Such alignments are first loaded using \code{\link{readPaf}} function into a \code{tibble} object.
#' Then loaded alignments are converted into the plotting coordinates using \code{\link{paf2coords}} function.
#' molecule. Resulting \code{data.frame} object contains required aesthetics `x` and `y` coordinates as well as 
#' `group` field that is required in order to determine which coordinates represent a single alignment.
#' 
#' @section Aesthetics:
#' `geom_miropeats()` require or can take the following aesthetics
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
#' @seealso [readPaf()], [paf2coords()]
#' @author David Porubsky
#' @export
#' @examples 
#'## Get PAF to plot
#'paf.file <- system.file("extdata", "test1.paf", package="SVbyEye")
#'## Read in PAF 
#'paf.table <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = 'cg')
#'## Convert PAF alignments to plotting coordinates
#'coords <- paf2coords(paf.table = paf.table)
#'## Make a plot 
#'ggplot2::ggplot(coords) +
#'  geom_miropeats(ggplot2::aes(x, y, group = group, fill = direction), alpha = 0.5)
#'
geom_miropeats <- function(mapping = NULL, data = NULL, geom = 'polygon',
                               stat = 'miropeats', position = 'identity',
                               n = 100, na.rm = FALSE, strength = 0.5,
                               show.legend = NA, inherit.aes = TRUE, ...) {
  ggplot2::layer(
    data = data, mapping = mapping, stat = stat, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, n = n, strength = strength, ...)
  )
}


#' @importFrom ggforce StatBezier
#' @format NULL
#' @usage NULL
#' @export
StatMiropeats <- ggplot2::ggproto('StatMiropeats', ggplot2::Stat,
                                  setup_data = function(data, params) {
                                    if (any(table(data$group) != 4)) {
                                      stop('Each group must consist of 4 points')
                                    }
                                    data
                                  },
                                  compute_panel = function(data, scales, strength = 0.5, n = 100) {
                                    data <- data[order(data$group, data$y), ]
                                    coords1 <- data[c(TRUE, FALSE, TRUE, FALSE), ]
                                    coords2 <- data[c(FALSE, TRUE, FALSE, TRUE), ]
                                    coords1 <- add.control.points(coords1, strength)
                                    coords1 <- ggforce::StatBezier$compute_panel(coords1, scales, n)
                                    coords2 <- add.control.points(coords2[rev(seq_len(nrow(coords2))), ], strength)
                                    coords2 <- ggforce::StatBezier$compute_panel(coords2, scales, n)
                                    diagonals <- rbind(coords1, coords2)
                                    
                                    # data <- data[order(data$group, data$x, data$y), ]
                                    # lower <- data[c(TRUE, FALSE, TRUE, FALSE), ]
                                    # upper <- data[c(FALSE, TRUE, FALSE, TRUE), ]
                                    # lower <- add_controls(lower, strength)
                                    # upper <- add_controls(upper[rev(seq_len(nrow(upper))), ], strength)
                                    # lower <- StatBezier$compute_panel(lower, scales, n)
                                    # upper <- StatBezier$compute_panel(upper, scales, n)
                                    # diagonals <- rbind(lower, upper)
                                    # diagonals$index <- NULL
                                    # diagonals[order(diagonals$group), ]
                                  },
                                  required_aes = c('x', 'y', 'group'),
                                  extra_params = c('na.rm', 'n', 'strength')
)