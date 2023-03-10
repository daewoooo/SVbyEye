#' Plot all-versus-all alignments stored in PAf format.
#'
#' This function takes PAF output file from minimap2 reporting all-versus-all alignments of multiple FASTA sequences
#' and visualize the alignments in a miropeat style.
#'
#' @param seqnames.order A user defined order sequence names to be plotted from top to the bottom.
#' @inheritParams readPaf
#' @inheritParams plotMiro
#' @return A \code{list} of miropeat style plots.
#' @import ggplot2
#' @importFrom grid unit
#' @importFrom scales comma
#' @importFrom dplyr bind_rows group_by summarise arrange
#' @importFrom GenomicAlignments explodeCigarOpLengths
#' @author David Porubsky
#' @export
#' @examples
#' ## Get PAF to plot
#' paf.file <- system.file("extdata", "test_ava.paf", package = "SVbyEye")
#' ## Read in PAF
#' paf.table <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = "cg")
#' ## Make a plot
#' ## Color by alignment directionality
#' plotAVA(paf.table = paf.table, color.by = "direction")
#' ## Color by fraction of matched bases in each alignment
#' plotAVA(paf.table = paf.table, color.by = "identity")
#' ## Define custom sample order
#' seqnames.order <- c("HG00438_2", "HG01358_2", "HG02630_2", "HG03453_2")
#' plotAVA(paf.table = paf.table, color.by = "direction", seqnames.order = seqnames.order)
#' ## Only samples present in custom sample order are being plotted
#' seqnames.order <- c("HG00438_2", "HG01358_2", "HG03453_2")
#' plotAVA(paf.table = paf.table, color.by = "direction", seqnames.order = seqnames.order)
#' ## Outline PAF alignments
#' plotAVA(paf.table = paf.table, outline.alignments = TRUE)
#' ## Highlight structural variants
#' plotAVA(
#'     paf.table = paf.table, min.deletion.size = 1000, min.insertion.size = 1000,
#'     highlight.sv = "outline"
#' )
#' ## Bin PAF alignments into user defined bin and color them by sequence identity (% of matched bases)
#' plotAVA(paf.table = paf.table, binsize = 10000)
#' ## Add annotation to self-alignments ##
#' plt <- plotAVA(paf.table = paf.table, color.by = 'direction')
#' annot.file <- system.file("extdata", "test_annot_ava.RData", package="SVbyEye")
#' annot.gr <- get(load(annot.file))
#' addAnnotation(ggplot.obj = plt, annot.gr = annot.gr, coordinate.space = 'self', y.label.id = 'ID')
#'
plotAVA <- function(paf.table, seqnames.order = NULL, min.deletion.size = NULL, min.insertion.size = NULL, highlight.sv = NULL, binsize = NULL, color.by = "direction", color.palette = NULL, outline.alignments = FALSE) {
    ## Check user input
    ## Make sure submitted paf.table has at least 12 mandatory fields
    if (ncol(paf.table) >= 12) {
        paf <- paf.table
        ## Add PAF alignment IDs if it doesn't exists
        if (!"aln.id" %in% colnames(paf)) {
            paf$aln.id <- seq_len(nrow(paf))
        }
    } else {
        stop("Submitted PAF alignments do not contain a minimum of 12 mandatory fields, see PAF file format definition !!!")
    }

  ## Get desired sequence order ##
  ## Based on user input
  seq.ids <- unique(paf$q.name)
  if (is.character(seqnames.order)) {
    if (all(seq.ids %in% seqnames.order)) {
      seq.ord <- seqnames.order
    } else {
      #seq.ord <- c(seqnames.order, setdiff(seq.ids, seqnames.order))
      message("Not all query and target IDs present in user defined 'seqnames.order', subsetting !!!")
      ## Subset PAF to only those samples defined in seqnames.order
      paf <- paf[paf$q.name %in% seqnames.order & paf$t.name %in% seqnames.order,]
      seq.ord <- seqnames.order
    }
  } else {
    seq.ord <- NULL
  }

    ## Break PAF at insertion/deletions defined in cigar string
    if (!is.null(min.deletion.size) | !is.null(min.insertion.size)) {
        paf.l <- breakPaf(paf.table = paf, min.deletion.size = min.deletion.size, min.insertion.size = min.insertion.size, collapse.mismatches = TRUE, report.sv = TRUE)
        paf <- paf.l$M
        paf.svs <- paf.l$SVs
    } else {
        # paf$aln.id <- 1:nrow(paf)
        paf.svs <- NULL
        if (!is.null(highlight.sv)) {
            highlight.sv <- NULL
            warning("Please specify 'min.deletion.size' and 'min.insertion.size' in order to make parameter 'highlight.sv' to work !!!")
        }
    }
    ## Store PAF alignments for later addition of 'geom_gene_arrow'
    paf.copy <- paf
    paf.copy$ID <- "M"

    ## Bin PAF alignments
    if (!is.null(binsize)) {
        if (binsize > 0) {
            if (binsize < 10) {
                binsize <- 10
                warning("Minimum allowed bin size is 10, forced binsize=10!!!")
            }
            paf <- pafToBins(paf.table = paf, binsize = binsize)
            ## If the PAF alignments are binned only color.by = 'fraction.matches' is allowed
            color.by <- "identity"
        }
    }
    ## Mark alignments ranges by 'M' (matches)
    paf$ID <- "M"

    ## Add SVs to the alignment table
    if (!is.null(paf.svs)) {
        if (nrow(paf.svs) > 0) {
            paf.svs$ID <- "INS"
            paf.svs$ID[grep(paf.svs$cg, pattern = "D", ignore.case = TRUE)] <- "DEL"
            paf <- dplyr::bind_rows(paf, paf.svs)
        }
    }

    ## If user sequence order not defined get one based on number of mismatches between all pairs of alignments ##
    ## Get unique alignment ID
    paf$seq.pair <- paste0(paf$q.name, "__", paf$t.name)
    if (is.null(seq.ord)) {
      ## Order alignments based on the number of mismatches
      # paf$NM <- S4Vectors::sapply(paf$cg, function(cg) sum(GenomicAlignments::explodeCigarOpLengths(cigar = cg, ops = c('X'))[[1]]), USE.NAMES = FALSE)
      paf$NM <- vapply(paf$cg, function(cg) sum(GenomicAlignments::explodeCigarOpLengths(cigar = cg, ops = c("X"))[[1]]), USE.NAMES = FALSE, FUN.VALUE = numeric(1))
      paf.ord <- paf %>%
        dplyr::group_by(seq.pair) %>%
        dplyr::summarise(n.nm = sum(.data$NM)) %>%
        dplyr::arrange(.data$n.nm)
      # paf.ord <- paf %>% dplyr::group_by(seq.pair) %>% dplyr::summarise(identity =  sum(n.match) / sum(aln.len)) %>% dplyr::arrange(identity)
      paf.ord.pairs <- unlist(strsplit(paf.ord$seq.pair, "__"))
      ## Assign level to seq.names ordered by number of matching bases in plus an minus orientation
      seq.names <- paf.ord.pairs[!duplicated(paf.ord.pairs)]
    } else {
      ## Assign user defined assembly order
      seq.names <- seq.ord
    }

    ## Order PAF
    # desired.pairs <- paste0(seq.names[-length(seq.names)], '__', seq.names[-1])
    # paf <- paf %>% dplyr::arrange(match(paf$q.name, seq.names))
    # paf <- paf[match(paf$q.name, seqnames.order),]
    # paf <- paf[paf$seq.pair %in% desired.pairs,]

    ## Flip start-end if strand == '-'
    paf[paf$strand == "-", c("t.start", "t.end")] <- rev(paf[paf$strand == "-", c("t.start", "t.end")])

    ## Assign level to seq.names
    seq.ids <- length(seq.names):1
    names(seq.ids) <- seq.names
    paf$y1 <- seq.ids[match(paf$q.name, names(seq.ids))]
    paf$y2 <- seq.ids[match(paf$t.name, names(seq.ids))]
    ## Keep subsequent comparisons only
    paf <- paf[abs(paf$y2 - paf$y1) == 1, ]

    ## Flip query and target for alignments where query comes first
    flipQT <- which(paf$y1 > paf$y2)
    # paf[flipQT,] <- transform(paf[flipQT,],
    #                           'q.name' = t.name, 'q.start' = t.start, 'q.end' = t.end,
    #                           't.name' = q.name, 't.start' = q.start, 't.end' = q.end,
    #                           'y1' = y2, 'y2' = y1)
    if (length(flipQT) > 0) {
        paf.sub <- paf[flipQT, ]
        paf[flipQT, ] <- transform(paf.sub,
            "q.name" = paf.sub$t.name, "q.start" = paf.sub$t.start, "q.end" = paf.sub$t.end,
            "t.name" = paf.sub$q.name, "t.start" = paf.sub$q.start, "t.end" = paf.sub$q.end,
            "y1" = paf.sub$y2, "y2" = paf.sub$y1
        )
    }
    paf$seq.pair[flipQT] <- paste0(paf$q.name[flipQT], "___", paf$t.name[flipQT])

    ## Translate paf alignments to plotting coordinates ##
    x <- c(rbind(paf$q.start, paf$t.start, paf$q.end, paf$t.end))
    y <- c(rbind(paf$y1, paf$y2, paf$y1, paf$y2))
    group <- rep(seq_len(nrow(paf)), each = 4)
    seq.name <- c(rbind(paf$q.name, paf$t.name, paf$q.name, paf$t.name))
    seq.pos <- c(rbind(paf$q.start, paf$t.start, paf$q.end, paf$t.end))
    seq.id <- c(rbind("query", "target", "query", "target"))
    n.match <- rep(paf$n.match, each = 4)
    aln.len <- rep(paf$aln.len, each = 4)
    mapq <- rep(paf$mapq, each = 4)
    # aln.id <- rep(paf$aln.id, each=4)
    ID <- rep(paf$ID, each = 4)
    seq.pair <- rep(paf$seq.pair, each = 4)
    direction <- rep(paf$strand, each = 4)

    coords <- data.frame(
        x = x,
        y = y,
        group = group,
        seq.pos = seq.pos,
        direction = direction,
        seq.name = seq.name,
        seq.id = seq.id,
        n.match = n.match,
        aln.len = aln.len,
        mapq = mapq,
        # aln.id=aln.id,
        ID = ID,
        # direction.flip=direction.flip,
        seq.pair = seq.pair,
        stringsAsFactors = FALSE
    )

    ## Get x-axis labels
    y.labels <- unique(coords$seq.name)
    y.breaks <- coords$y[match(y.labels, coords$seq.name)]

    ## Define color palette
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

    ## Plot alignments and color by direction or identity
    if (color.by == "direction") {
        plt <- ggplot2::ggplot(coords[coords$ID == "M", ]) +
            geom_miropeats(ggplot2::aes(x, y, group = group, fill = direction), alpha = 0.5) +
            ggplot2::scale_fill_manual(values = pal, name = "Alignment\ndirection")
    } else if (color.by == "identity") {
        coords$identity <- (coords$n.match / coords$aln.len) * 100
        coords$identity[is.nan(coords$identity) | is.na(coords$identity)] <- 0
        ## Define color scheme
        coords.l <- getColorScheme(data.table = coords, value.field = "identity", breaks = c(90, 95, 99, 99.5, 99.6, 99.7, 99.8, 99.9))
        coords <- coords.l$data
        colors <- coords.l$colors

        plt <- ggplot2::ggplot(coords[coords$ID == "M", ]) +
            geom_miropeats(ggplot2::aes(x, y, group = group, fill = .data$col.levels), alpha = 0.5) +
            ggplot2::scale_fill_manual(values = colors, drop = FALSE, name = "Identity")
    } else if (color.by == "mapq") {
        plt <- ggplot2::ggplot(coords[coords$ID == "M", ]) +
            geom_miropeats(ggplot2::aes(x, y, group = group, fill = mapq), alpha = 0.5) +
            ggplot2::scale_fill_gradient(low = "gray", high = "red")
    } else {
        plt <- ggplot2::ggplot(coords[coords$ID == "M", ]) +
            geom_miropeats(ggplot2::aes(x, y, group = group), alpha = 0.5, fill = "gray")
    }

    ## Add alignment outlines
    if (outline.alignments) {
        plt <- plt + geom_miropeats(data = coords[coords$ID == "M", ], ggplot2::aes(x, y, group = group), fill = NA, color = "gray", linewidth = 0.25)
    }

    ## Add indels
    if (!is.null(highlight.sv)) {
        if (nrow(coords[coords$ID != "M", ]) > 0) {
            ## Add SVs to the plot
            if (highlight.sv == "outline") {
                plt <- plt + ggnewscale::new_scale_color() +
                    geom_miropeats(data = coords[coords$ID != "M", ], ggplot2::aes(x, y, group = group, color = ID), fill = NA, alpha = 0.5, inherit.aes = FALSE) +
                    ggplot2::scale_color_manual(values = c("DEL" = "firebrick3", "INS" = "dodgerblue3"), name = "SV class")
            } else if (highlight.sv == "fill") {
                plt <- plt + ggnewscale::new_scale_fill() +
                    geom_miropeats(data = coords[coords$ID != "M", ], ggplot2::aes(x, y, group = group, fill = ID), alpha = 0.5, inherit.aes = FALSE) +
                    ggplot2::scale_fill_manual(values = c("DEL" = "firebrick3", "INS" = "dodgerblue3"), name = "SV class")
            } else {
                warning("Parameter 'highlight.sv' can only take values 'outline' or 'fill', see function documentation!!!")
            }
        } else {
            message("There are no SVs to highlight. Try to decrease 'min.deletion.size' and 'min.insertion.size' values!!!")
        }
    }

    ## Add x and y scales
    suppressMessages(
        plt <- plt +
            ggplot2::scale_x_continuous(expand = c(0, 0), labels = scales::comma) +
            ggplot2::scale_y_continuous(breaks = y.breaks, labels = y.labels) +
            ggplot2::xlab("Genomic position (bp)") +
            ggplot2::ylab("")
    )

    ## Add arrows to mark start and end of each alignment
    # start <- coords$x[c(T, T, F, F)]
    # end <- coords$x[c(F, F, T, T)]
    # y <- coords$y[c(T, T, F, F)]
    # group <- coords$group[c(T, T, F, F)]
    # plt.df <- data.frame(start=start, end=end, y=y, group=group)
    # plt.df$direction <- ifelse(plt.df$start < plt.df$end, '+', '-')
    #
    # plt <- plt + ggnewscale::new_scale_fill() + ggnewscale::new_scale_color() +
    #   gggenes::geom_gene_arrow(data=plt.df, ggplot2::aes(xmin = start, xmax = end, y = y, color= direction, fill = direction), arrowhead_height = unit(3, 'mm')) +
    #   scale_fill_manual(values = c('-' = 'cornflowerblue', '+' = 'forestgreen'), name='Alignment\ndirection') +
    #   scale_color_manual(values = c('-' = 'cornflowerblue', '+' = 'forestgreen'), name='Alignment\ndirection')

    ## Add sequence length lines
    seq.lines <- data.frame(y.breaks = y.breaks, y.labels = y.labels)
    seq.lengths <- dplyr::bind_cols("seq.id" = c(paf.copy$q.name, paf.copy$t.name), "seq.len" = c(paf.copy$q.len, paf.copy$t.len))
    seq.lengths <- seq.lengths[!duplicated(seq.lengths$seq.id), ]
    seq.lines$seq.len <- seq.lengths$seq.len[match(seq.lines$y.labels, seq.lengths$seq.id)]
    plt <- plt + ggplot2::geom_segment(data = seq.lines, ggplot2::aes(x = 1, xend = .data$seq.len, y = y.breaks, yend = y.breaks), linewidth = 1)

    ## Set the theme and scales
    theme_ava <- ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        axis.line.x = ggplot2::element_line(linewidth = 1),
        axis.ticks.x = ggplot2::element_line(linewidth = 1),
        axis.ticks.length.x = grid::unit(2, "mm")
    )
    plt <- plt + theme_ava

    ## Return final plot
    return(plt)
}



# paf2coords_test <- function(paf.table, offset.alignments=FALSE) {
#   ## Check user input ##
#   ## Make sure PAF has at least 12 mandatory fields
#   if (ncol(paf.table) >= 12) {
#     paf <- paf.table
#   } else {
#     stop('Submitted PAF alignments do not contain a minimum of 12 mandatory fields, see PAF file format definition !!!')
#   }
#
#   ## Flip start-end if strand == '-'
#   paf[paf$strand == '-', c('t.start','t.end')] <- rev(paf[paf$strand == '-', c('t.start','t.end')])
#   #paf[paf$strand == '-', c('q.start','q.end')] <- rev(paf[paf$strand == '-', c('q.start','q.end')])
#
#   ## Get subsequent seqnames order
#   seq.names <- unique(c(paf$q.name, paf$t.name))
#   #seq.names <- unique(c(rbind(paf$q.name, paf$t.name)))
#
#   ## Get unique alignment ID
#   if (!'seq.pair' %in% colnames(paf)) {
#     paf$seq.pair <- paste0(paf$q.name, '__', paf$t.name)
#   }
#
#   ## Sync scales between alignments [per region id]
#   paf.l <- split(paf, paf$seq.pair)
#   for (i in seq_along(paf.l)) {
#     paf.sub <- paf.l[[i]]
#     q.range <- range(c(paf.sub$q.start, paf.sub$q.end))
#     t.range <- range(c(paf.sub$t.start, paf.sub$t.end))
#     ## Adjust target ranges given the size difference with respect query ranges
#     range.offset <- diff(q.range) - diff(t.range)
#     t.range[2] <- t.range[2] + range.offset ## Make a start position as offset and change only end position
#     ## Covert query to target coordinates
#     paf.sub$q.start.trans <- q2t(x = paf.sub$q.start, q.range = q.range, t.range = t.range)
#     paf.sub$q.end.trans <- q2t(x = paf.sub$q.end, q.range = q.range, t.range = t.range)
#     # q.range <- range(c(paf$q.start, paf$q.end))
#     # t.range <- range(c(paf$t.start, paf$t.end))
#     # paf$q.start.trans <- q2t(x = paf$q.start, q.range = q.range, t.range = t.range)
#     # paf$q.end.trans <- q2t(x = paf$q.end, q.range = q.range, t.range = t.range)
#     paf.l[[i]] <- paf.sub
#   }
#   paf <- do.call(rbind, paf.l)
#
#   ## Assign level to seq.names
#   seq.ids <- length(seq.names):1
#   names(seq.ids) <- seq.names
#   paf$y1 <- seq.ids[match(paf$q.name, names(seq.ids))]
#   paf$y2 <- seq.ids[match(paf$t.name, names(seq.ids))]
#   ## Keep subsequent comparisons only
#   paf <- paf[abs(paf$y2 - paf$y1) == 1,]
#
#   ## Flip query and target for alignments where query comes first
#   flipQT <- which(paf$y1 > paf$y2)
#   paf[flipQT,] <- transform(paf[flipQT,],
#                             q.name = t.name, q.start = t.start, q.end = t.end,
#                             t.name = q.name, t.start = q.start, t.end = q.end,
#                             y1 = y2, y2 = y1)
#   paf$seq.pair[flipQT] <- paste0(paf$q.name[flipQT], '___', paf$t.name[flipQT])
#
#   ## Vectorize data transformation ##
#   #x <- c(rbind(paf$q.start.trans, paf$t.start, paf$q.end.trans, paf$t.end))
#   #y <- rep(c(1,2,1,2), times=nrow(paf))
#   x <- c(rbind(paf$q.start, paf$t.start, paf$q.end, paf$t.end))
#   y <- c(rbind(paf$y1, paf$y2, paf$y1, paf$y2))
#   ## Offset target alignments
#   if (offset.alignments) {
#     offset <- rep(c(0,0,0,0,0,0.05,0,0.05), times=ceiling(nrow(paf) / 2))[1:length(y)]
#     y <- y + offset
#   }
#   group <- rep(1:nrow(paf), each=4)
#   seq.name <- c(rbind(paf$q.name, paf$t.name, paf$q.name, paf$t.name))
#   seq.pos <- c(rbind(paf$q.start, paf$t.start, paf$q.end, paf$t.end))
#   seq.id <- c(rbind('query', 'target', 'query', 'target'))
#   n.match <- rep(paf$n.match, each=4)
#   aln.len <- rep(paf$aln.len, each=4)
#   mapq <- rep(paf$mapq, each=4)
#   aln.id <- rep(paf$aln.id, each=4)
#   ID <- rep(paf$ID, each=4)
#   seq.pair <- rep(paf$seq.pair, each=4)
#   direction <- rep(paf$strand, each=4)
#
#   ## Create final data coordinate data frame
#   coords <- data.frame(x=x,
#                        y=y,
#                        group=group,
#                        seq.pos=seq.pos,
#                        direction=direction,
#                        seq.name=seq.name,
#                        seq.id=seq.id,
#                        n.match=n.match,
#                        aln.len=aln.len,
#                        mapq=mapq,
#                        aln.id=aln.id,
#                        ID=ID,
#                        #direction.flip=direction.flip,
#                        seq.pair=seq.pair,
#                        stringsAsFactors = FALSE)
#
#   return(coords)
# }
#
