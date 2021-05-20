#' @rdname ggforce-extensions
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
# @rdname stat_miropeats
# @export
# stat_miropeats <- function(mapping = NULL, data = NULL, geom = 'polygon',
#                                position = 'identity', n = 100, strength = 0.5,
#                                na.rm = FALSE, show.legend = NA,
#                                inherit.aes = TRUE, ...) {
#   ggplot2::layer(
#     stat = StatMiropeats, data = data, mapping = mapping, geom = geom,
#     position = position, show.legend = show.legend, inherit.aes = inherit.aes,
#     params = list(na.rm = na.rm, n = n, strength = strength, ...)
#   )
# }

#' @rdname geom_miropeats
#' @export
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
