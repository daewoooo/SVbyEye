#' Visualize PAF alignments.
#'
#' This function takes PAF output file from minimap2 alignments, and visualize the alignments
#' in a miropeat style.
#'
#' @param highlight.sv Visualize alignment embedded structural variation either as an outlined ('outline') or filled ('fill') miropeats.
#' @param color.by Color alignments either by directionality ('direction'), fraction of matched base pairs ('identity'),
#' or a custom column name present in submitted `paf.table`.
#' @param color.palette A discrete color palette defined as named character vector (elements = colors, names = discrete levels)
#' to color alignment directionality, `[default: color.palette <- c('-' = 'cornflowerblue', '+' = 'forestgreen')]`.
#' @param outline.alignments Set to \code{TRUE} if boundaries of each alignment should be highlighted by gray outline.
#' @param genomic.scale Report genomic coordinates in base pairs ('bp') kilobase pairs ('kbp') or megabase pairs ('Mbp') `[default: 'bp']`.
#' @inheritParams breakPaf
#' @inheritParams pafAlignmentToBins
#' @inheritParams paf2coords
#' @inheritParams cutPafAlignments
#' @return A \code{ggplot2} object.
#' @import ggplot2
#' @importFrom grid unit
#' @importFrom scales comma
#' @importFrom wesanderson wes_palette
#' @importFrom gggenes geom_gene_arrow
#' @importFrom GenomicRanges start end
#' @importFrom ggnewscale new_scale_fill new_scale_color
#' @author David Porubsky
#' @export
#' @examples
#' ## Get PAF to plot ##
#' paf.file <- system.file("extdata", "test1.paf", package = "SVbyEye")
#' ## Read in PAF
#' paf.table <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = "cg")
#' ## Optional steps include PAF filtering and flipping query coordinates
#' ## (see filterPaf and flipPaf function documentation)
#' ## Make a plot ##
#' ## Color by alignment directionality
#' plotMiro(paf.table = paf.table, color.by = "direction")
#' ## Color by fraction of matched bases in each alignment
#' plotMiro(paf.table = paf.table, color.by = "identity")
#' ## Use custom color palette to color alignment directionality
#' plotMiro(paf.table = paf.table, color.palette = c("+" = "azure3", "-" = "yellow3"))
#' ## Change genomic scale to megabase pairs
#' plotMiro(paf.table = paf.table, genomic.scale = "Mbp")
#' ## Outline PAF alignments
#' plotMiro(paf.table = paf.table, outline.alignments = TRUE)
#' ## Offset target PAF alignments
#' plotMiro(paf.table = paf.table, offset.alignments = TRUE)
#' ## Bin PAF alignments into user defined bin and color them by sequence identity (% of matched bases)
#' plotMiro(paf.table = paf.table, binsize = 10000)
#' ## Highlight structural variants
#' paf.file <- system.file("extdata", "test3.paf", package = "SVbyEye")
#' paf.table <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = "cg")
#' plotMiro(paf.table = paf.table, min.deletion.size = 50, highlight.sv = "outline")
#'
plotMiro <- function(paf.table, min.deletion.size = NULL, min.insertion.size = NULL, highlight.sv = NULL, binsize = NULL, color.by = "direction", color.palette = NULL, outline.alignments = FALSE, offset.alignments = FALSE, target.region = NULL, genomic.scale = "bp") {
    ## Check user input ##
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
    ## Make sure submitted paf.table contains single query and target sequence id
    if (length(unique(paf.table$q.name)) > 1) {
        stop("Currently, function [plotMiro] does not support visualization of more than one query sequence id !!!")
    }
    if (length(unique(paf.table$t.name)) > 1) {
        stop("Currently, function [plotMiro] does not support visualization of more than one target sequence id !!!")
    }

    ## Subset/cut PAF alignments to user defined target region
    if (!is.null(target.region)) {
        paf <- cutPafAlignments(paf.table = paf, target.region = target.region)
    }

    ## Break PAF at insertion/deletions defined in cigar string
    if (!is.null(min.deletion.size) | !is.null(min.insertion.size)) {
        paf.l <- breakPaf(paf.table = paf.table, min.deletion.size = min.deletion.size, min.insertion.size = min.insertion.size, collapse.mismatches = TRUE, report.sv = TRUE)
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

    ## Convert PAF alignments to plotting coordinates
    if (color.by %in% colnames(paf)) {
        coords <- paf2coords(paf.table = paf, offset.alignments = offset.alignments, add.col = color.by)
    } else {
        coords <- paf2coords(paf.table = paf, offset.alignments = offset.alignments)
    }

    ## Prepare data for plotting
    target.seqname <- unique(coords$seq.name[coords$seq.id == "target"])
    ## Get y-axis labels
    q.range <- range(coords$seq.pos[coords$seq.id == "query"])
    t.range <- range(coords$seq.pos[coords$seq.id == "target"])
    ## Adjust target ranges given the size difference with respect to query ranges
    range.offset <- diff(q.range) - diff(t.range)
    t.range[2] <- t.range[2] + range.offset ## Make a start position as offset and change only end position
    ## Get x-axis labels
    q.labels <- pretty(q.range)
    t.labels <- pretty(t.range)
    ## Covert query to target coordinates
    q.breaks <- SVbyEye::q2t(x = q.labels, q.range = q.range, t.range = t.range)
    t.breaks <- t.labels

    ## Convert axis labels to desired genomic scale
    if (genomic.scale == "kbp") {
        q.labels <- round(abs(q.labels) / 1000, digits = 3)
        t.labels <- round(abs(t.labels) / 1000, digits = 3)
    } else if (genomic.scale == "Mbp") {
        q.labels <- round(abs(q.labels) / 1000000, digits = 1)
        t.labels <- round(abs(t.labels) / 1000000, digits = 1)
    } else {
        ## Make sure axis labels are always positive numbers
        q.labels <- abs(q.labels)
        t.labels <- abs(t.labels)
    }

    ## Get y-axis labels
    seq.labels <- c(
        unique(coords$seq.name[coords$seq.id == "query"]),
        unique(coords$seq.name[coords$seq.id == "target"])
    )

    ## Define color palette for alignment directionality
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

    ## Plot alignments and color by a user defined variable
    if (color.by == "direction") {
        plt <- ggplot2::ggplot(coords[coords$ID == "M", ]) +
            geom_miropeats(ggplot2::aes(x = .data$x, y = .data$y, group = .data$group, fill = .data$direction), alpha = 0.5) +
            ggplot2::scale_fill_manual(values = pal, name = "Alignment\ndirection")
    } else if (color.by == "identity") {
        coords$identity <- (coords$n.match / coords$aln.len) * 100
        coords$identity[is.nan(coords$identity) | is.na(coords$identity)] <- 0
        ## Define color scheme
        coords.l <- getColorScheme(data.table = coords, value.field = "identity", breaks = c(90, 95, 99, 99.5, 99.6, 99.7, 99.8, 99.9))
        coords <- coords.l$data
        colors <- coords.l$colors

        plt <- ggplot2::ggplot(coords[coords$ID == "M", ]) +
            geom_miropeats(ggplot2::aes(x = .data$x, y = .data$y, group = .data$group, fill = .data$col.levels), alpha = 0.5) +
            ggplot2::scale_fill_manual(values = colors, drop = FALSE, name = "Identity")
    } else if (color.by %in% colnames(paf)) {
        ## Define color scheme
        coords.l <- getColorScheme(data.table = coords, value.field = color.by)
        coords <- coords.l$data
        colors <- coords.l$colors

        plt <- ggplot2::ggplot(coords[coords$ID == "M", ]) +
            geom_miropeats(ggplot2::aes(x = .data$x, y = .data$y, group = .data$group, fill = .data$col.levels), alpha = 0.5) +
            ggplot2::scale_fill_manual(values = colors, drop = FALSE, name = eval(color.by))
    } else {
        plt <- ggplot2::ggplot(coords[coords$ID == "M", ]) +
            geom_miropeats(ggplot2::aes(x = .data$x, y = .data$y, group = .data$group), alpha = 0.5, fill = "gray")
    }

    ## Add alignment outlines
    if (outline.alignments) {
        plt <- plt +
            geom_miropeats(data = coords[coords$ID == "M", ], ggplot2::aes(x = .data$x, y = .data$y, group = .data$group), fill = NA, color = "gray", linewidth = 0.25)
    }

    ## Add indels
    if (!is.null(highlight.sv)) {
        if (nrow(coords[coords$ID != "M", ]) > 0) {
            ## Add SVs to the plot
            if (highlight.sv == "outline") {
                plt <- plt + ggnewscale::new_scale_color() +
                    geom_miropeats(data = coords[coords$ID != "M", ], ggplot2::aes(x = .data$x, y = .data$y, group = .data$group, color = .data$ID), fill = NA, alpha = 0.5, inherit.aes = FALSE) +
                    ggplot2::scale_color_manual(values = c("DEL" = "firebrick3", "INS" = "dodgerblue3"), name = "SV class")
            } else if (highlight.sv == "fill") {
                plt <- plt + ggnewscale::new_scale_fill() +
                    geom_miropeats(data = coords[coords$ID != "M", ], ggplot2::aes(x = .data$x, y = .data$y, group = .data$group, fill = .data$ID), alpha = 0.5, inherit.aes = FALSE) +
                    ggplot2::scale_fill_manual(values = c("DEL" = "firebrick3", "INS" = "dodgerblue3"), name = "SV class")
            } else {
                warning("Parameter 'highlight.sv' can only take values 'outline' or 'fill', see function documentation!!!")
            }
        } else {
            warning("There are no SVs to highlight. Try to decrease 'min.deletion.size' and 'min.insertion.size' values!!!")
        }
    }

    ## Add custom x and y scales
    suppressMessages(
        plt <- plt +
            ggplot2::scale_y_continuous(breaks = c(1, 2), labels = seq.labels) +
            ggplot2::scale_x_continuous(
                breaks = q.breaks, labels = scales::comma(q.labels),
                sec.axis = ggplot2::sec_axis(trans = y ~ ., breaks = t.breaks, labels = scales::comma(t.labels)), expand = c(0, 0)
            ) +
            # xlab('Genomic position (bp)') +
            ggplot2::xlab(paste0("Genomic position (", genomic.scale, ")")) +
            ggplot2::ylab("")
    )

    ## Set plot boundaries based on user defined target region
    if (!is.null(target.region)) {
        if (is.character(target.region)) {
            target.region.gr <- as(target.region, "GRanges")
            plt <- plt + ggplot2::coord_cartesian(xlim = c(GenomicRanges::start(target.region.gr), GenomicRanges::end(target.region.gr)))
        } else if (is(target.region.gr, "GRanges")) {
            target.region.gr <- target.region
            plt <- plt + ggplot2::coord_cartesian(xlim = c(GenomicRanges::start(target.region.gr), GenomicRanges::end(target.region.gr)))
        } else {
            warning("User defined 'target.region' has be either character string 'chr:start-end' or GenomicRanges object, skipping!!!")
        }
    }

    ## Add arrows to mark start and end of each alignment ##
    ## Use an un-binned version of PAF alignments (group by bin ID if present)
    if ("bin.id" %in% colnames(paf.copy)) {
        paf.copy <- collapsePaf(paf.table = paf.copy, collapse.by = "bin.id")
    }
    ## Convert to plotting coordinates
    coords.arrow <- paf2coords(paf.table = paf.copy, offset.alignments = offset.alignments)
    start <- coords.arrow$x[c(T, T, F, F)]
    end <- coords.arrow$x[c(F, F, T, T)]
    y <- coords.arrow$y[c(T, T, F, F)]
    group <- coords.arrow$group[c(T, T, F, F)]
    plt.df <- data.frame(start = start, end = end, y = y, group = group)
    plt.df$direction <- ifelse(plt.df$start < plt.df$end, "+", "-")
    ## Add arrows to the plot
    plt <- plt + ggnewscale::new_scale_fill() + ggnewscale::new_scale_color() +
        gggenes::geom_gene_arrow(data = plt.df, ggplot2::aes(xmin = .data$start, xmax = .data$end, y = .data$y, color = .data$direction, fill = .data$direction), arrowhead_height = grid::unit(3, "mm")) +
        ggplot2::scale_fill_manual(values = pal, name = "Alignment\ndirection") +
        ggplot2::scale_color_manual(values = pal, name = "Alignment\ndirection")

    ## Set the theme
    theme_miro <- ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        axis.line.x = ggplot2::element_line(linewidth = 1),
        axis.ticks.x = ggplot2::element_line(linewidth = 1),
        axis.ticks.length.x = grid::unit(2, "mm")
    )
    plt <- plt + theme_miro

    ## Return final plot
    return(plt)
}
