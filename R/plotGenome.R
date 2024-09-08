#' Visualize genome-wide PAF alignments.
#'
#' This function takes genome-wide alignments of de novo assembly to the reference genome PAF format
#' and visualize the alignments with respect to all reference chromosomes in a single plot
#'
#' @param chromosomes User defined chromosomes (target sequence IDs) to be plotted.
#' @param chromosome.bar.width A width of chromosome/contig bars over each alignment `[Default : 2mm]`.
#' @param min.query.aligned.bp Only contigs with this minimum number of aligned base pairs will be plotted.
#' @param target.centromeres A \code{\link{GRanges-class}} object containing positions of centromeres in target chromosomes.
#' @inheritParams breakPaf
#' @inheritParams pafAlignmentToBins
#' @inheritParams plotMiro
#' @return A \code{ggplot2} object.
#' @import ggplot2
#' @importFrom grid unit
#' @importFrom methods is
#' @importFrom scales comma
#' @importFrom randomcoloR randomColor
#' @importFrom ggnewscale new_scale_fill new_scale_color
#' @importFrom dplyr group_by filter reframe bind_rows
#' @importFrom GenomeInfoDb keepSeqlevels
#' @author David Porubsky
#' @export
#' @examples
#' ## Get PAF to plot ##
#' paf.file <- system.file("extdata", "PTR_test.paf", package = "SVbyEye")
#' ## Read in PAF
#' paf.table <- readPaf(paf.file = paf.file)
#' ## Make a plot ##
#' ## Color by alignment directionality
#' plotGenome(paf.table = paf.table, chromosome.bar.width = grid::unit(2, 'mm'),
#'            min.query.aligned.bp = 5000000)
#'
plotGenome <- function(paf.table, min.deletion.size = NULL, min.insertion.size = NULL, highlight.sv = NULL, binsize = NULL, color.by = "direction", color.palette = NULL, perc.identity.breaks = c(90, 95, 99, 99.5, 99.6, 99.7, 99.8, 99.9), chromosomes = NULL, chromosome.bar.width = grid::unit(2, "mm"), min.query.aligned.bp = NULL, genomic.scale = "bp", target.centromeres = NULL) {
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

  ## Subset to desired target sequences chromosomes
  targets <- unique(paf.table$t.name)
  if (!is.null(chromosomes)) {
    if (any(chromosomes %in% targets)) {
      paf <- paf[paf$t.name %in% chromosomes,]
    } else {
      warning("None of the sequence names defined in 'chromosomes' is present in submitted PAF file!!!")
      chromosomes <- NULL
    }
  }

  ## Keep only query sequences whose total alignment size is >= user defined threshold
  if (!is.null(min.query.aligned.bp)) {
    query.aln.summary <- paf %>% dplyr::group_by(.data$t.name, .data$q.name) %>% dplyr::reframe(total.bp = sum(.data$aln.len))
    ## Get query sequences to keep
    query2keep <- query.aln.summary %>% dplyr::group_by(.data$t.name) %>% dplyr::filter(.data$total.bp >= min.query.aligned.bp)
    if (nrow(query2keep) > 0) {
      paf.l <- split(paf, paf$t.name)
      for (i in seq_along(paf.l)) {
        sub.paf <- paf.l[[i]]
        paf.l[[i]] <- sub.paf[sub.paf$q.name %in% query2keep$q.name[query2keep$t.name == names(paf.l[i])],]
      }
      paf <- dplyr::bind_rows(paf.l)
    } else {
      stop("None of the query sequences meet the 'min.query.aligned.bp' threshold!!!")
    }
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

  ## Process by target sequence to get plotting coordinates
  paf.l <- split(paf, paf$t.name)
  max.group <- 0
  coords.l <- list()
  for (i in seq_along(paf.l)) {
    sub.paf <- paf.l[[i]]
    ## Sync alignment orientation
    sub.paf <- suppressMessages( flipPaf(paf.table = sub.paf, majority.strand = '+') )
    ## Convert PAF alignments to plotting coordinates
    coords <- paf2coords(paf.table = sub.paf)
    ## Make sure unique groups
    coords$group <- coords$group + max.group
    max.group <- max(coords$group)
    ## Make sure query and target sequences appear on the same level
    coords$y[coords$seq.id == 'target'] <- 2
    coords$y[coords$seq.id == 'query'] <- 1
    ## Add target ID
    coords$target.id <- names(paf.l[i])
    ## Store
    coords.l[[i]] <- coords
  }
  coords <- dplyr::bind_rows(coords.l)

  ## Order target.id
  if (requireNamespace("gtools", quietly = TRUE)) {
    coords$target.id <- factor(coords$target.id, levels = gtools::mixedsort(unique(coords$target.id)))
  } else {
    print("Install R package 'gtools' to have chromosome order alphanumerically!!!")
    coords$target.id <- factor(coords$target.id, levels = unique(coords$target.id))
  }

  ## Get axis breaks and labels
  t.range <- range(coords$seq.pos[coords$seq.id == "target"])
  t.labels <- pretty(t.range)
  t.breaks <- t.labels

  ## Convert axis labels to desired genomic scale
  if (genomic.scale == "kbp") {
    t.labels <- round(abs(t.labels) / 1000, digits = 3)
  } else if (genomic.scale == "Mbp") {
    t.labels <- round(abs(t.labels) / 1000000, digits = 1)
  } else {
    ## Make sure axis labels are always positive numbers
    t.labels <- abs(t.labels)
  }

  ## Set default direction color
  if (!is.null(color.palette)) {
    if (all(c("+", "-") %in% names(color.palette))) {
      if (is.list(color.palette)) {
        pal <- unlist(color.palette, use.names = TRUE)
      } else {
        pal <-  color.palette
      }
    } else {
      pal <- c("-" = "cornflowerblue", "+" = "forestgreen")
      #warning("User defined 'color.palette' does not contain both '+' and '-' directions, using default values instead!!!")
    }
  } else {
    pal <- c("-" = "cornflowerblue", "+" = "forestgreen")
  }

  ## Plot alignments and color by a user defined variable
  if (color.by == "direction") {
    ## Make a plot
    plt <- ggplot2::ggplot(coords[coords$ID == "M",]) +
      geom_miropeats(ggplot2::aes(x = .data$x, y = .data$y, group = .data$group, fill = .data$direction), alpha = 0.5) +
      ggplot2::scale_fill_manual(values = pal, name = "Alignment\ndirection") +
      ggplot2::facet_grid(target.id ~ ., switch='y')
  } else if (color.by == "identity") {
    coords$identity <- (coords$n.match / coords$aln.len) * 100
    coords$identity[is.nan(coords$identity) | is.na(coords$identity)] <- 0
    ## Define color scheme
    coords.l <- getColorScheme(data.table = coords, value.field = "identity", breaks = perc.identity.breaks)
    coords <- coords.l$data
    colors <- coords.l$colors
    ## Make a plot
    plt <- ggplot2::ggplot(coords[coords$ID == "M",]) +
      geom_miropeats(ggplot2::aes(x = .data$x, y = .data$y, group = .data$group, fill = .data$col.levels), alpha = 0.5) +
      ggplot2::scale_fill_manual(values = colors, drop = FALSE, name = "Identity") +
      ggplot2::facet_grid(target.id ~ ., switch='y')
  } else if (color.by %in% colnames(paf)) {
    col.levels <- unique(coords[,eval(color.by), drop = TRUE])
    if (!all(col.levels %in% names(color.palette))) {
      ## Define color random scheme
      coords.l <- getColorScheme(data.table = coords, value.field = color.by)
      colors <- coords.l$colors
    } else {
      colors <- color.palette
    }
    ## Make a plot
    plt <- ggplot2::ggplot(coords[coords$ID == "M",]) +
      geom_miropeats(ggplot2::aes(x = .data$x, y = .data$y, group = .data$group, fill = .data[[color.by]]), alpha = 0.5) +
      ggplot2::scale_fill_manual(values = colors, drop = FALSE, name = eval(color.by)) +
      ggplot2::facet_grid(target.id ~ .)
  } else {
    plt <- ggplot2::ggplot(coords[coords$ID == "M",]) +
      geom_miropeats(ggplot2::aes(x = .data$x, y = .data$y, group = .data$group), alpha = 0.5, fill = "gray") +
      ggplot2::facet_grid(target.id ~ ., switch='y')
  }

  ## Add indels
  if (!is.null(highlight.sv)) {
    if (nrow(coords[coords$ID != "M",]) > 0) {
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
      message("There are no SVs to highlight. Try to decrease 'min.deletion.size' and 'min.insertion.size' values!!!")
    }
  }

  ## Add custom x and y scales
  suppressMessages(
    plt <- plt +
      ggplot2::scale_y_continuous(expand = c(0.05, 0.05)) +
      ggplot2::scale_x_continuous(breaks = t.breaks, labels = scales::comma(t.labels), expand = c(0, 0)) +
      ggplot2::xlab(paste0("Genomic position (", genomic.scale, ")")) +
      ggplot2::ylab("")
  )

  ## Add chromosome and contig bars ##
  coords.m <- coords[coords$ID == "M",]
  coords.m.l <- split(coords.m, coords.m$target.id)
  ranges.l <- list()
  for (i in seq_along(coords.m.l)) {
    sub.coords <- coords.m.l[[i]]
    ## Get query ranges
    q.ranges <- sub.coords %>% dplyr::filter(seq.id == 'query') %>%
      dplyr::group_by(seq.name) %>% dplyr::reframe(pos = range(.data$x))
    start <- q.ranges$pos[c(TRUE, FALSE)]
    end <- q.ranges$pos[c(FALSE, TRUE)]
    y <- 1
    q.ranges <- data.frame(start = start, end = end, y = y)
    ## Get target ranges
    t.ranges <- sub.coords %>% dplyr::filter(seq.id == 'target') %>%
      dplyr::group_by(seq.name) %>% dplyr::reframe(pos = range(.data$x))
    start <- t.ranges$pos[c(TRUE, FALSE)]
    end <- t.ranges$pos[c(FALSE, TRUE)]
    y <- 2
    t.ranges <- data.frame(start = start, end = end, y = y)
    ## Merge
    ranges <- dplyr::bind_rows(t.ranges, q.ranges)
    ranges$target.id <- names(coords.m.l[i])
    ranges.l[[i]] <- ranges
  }
  plt.data <- dplyr::bind_rows(ranges.l)
  plt.data$target.id <- factor(plt.data$target.id, levels = levels(coords$target.id))
  plt.data$width <- (plt.data$end - plt.data$start)
  ## Get unique color per target id
  n.uniq <- length(unique(plt.data$target.id))
  target.col <- randomcoloR::randomColor(count = n.uniq)
  target.pal <- setNames(target.col, unique(plt.data$target.id))
  plt <- plt + ggnewscale::new_scale_fill() + ggnewscale::new_scale_color() +
    geom_roundrect(data = plt.data[plt.data$y == 1 & plt.data$width >= 10000000,], ggplot2::aes(xmin = .data$start, xmax = .data$end, y = .data$y), fill='gray80', rect_height = chromosome.bar.width, radius = chromosome.bar.width / 2) +
    geom_roundrect(data = plt.data[plt.data$y == 1 & (plt.data$width >= 1000000 & plt.data$width < 10000000),], ggplot2::aes(xmin = .data$start, xmax = .data$end, y = .data$y), fill='gray80', rect_height = chromosome.bar.width, radius = chromosome.bar.width / 4) +
    geom_roundrect(data = plt.data[plt.data$y == 1 & (plt.data$width >= 100000 & plt.data$width < 1000000),], ggplot2::aes(xmin = .data$start, xmax = .data$end, y = .data$y), fill='gray80', rect_height = chromosome.bar.width, radius = chromosome.bar.width / 8) +
    geom_roundrect(data = plt.data[plt.data$y == 1 & plt.data$width < 100000,], ggplot2::aes(xmin = .data$start, xmax = .data$end, y = .data$y), fill='gray80', rect_height = chromosome.bar.width, radius = chromosome.bar.width * 0) +
    geom_roundrect(data = plt.data[plt.data$y == 2 & plt.data$width >= 10000000,], ggplot2::aes(xmin = .data$start, xmax = .data$end, y = .data$y, fill = .data$target.id), rect_height = chromosome.bar.width, radius = chromosome.bar.width / 2) +
    geom_roundrect(data = plt.data[plt.data$y == 2 & (plt.data$width >= 1000000 & plt.data$width < 10000000),], ggplot2::aes(xmin = .data$start, xmax = .data$end, y = .data$y, fill = .data$target.id), rect_height = chromosome.bar.width, radius = chromosome.bar.width / 4) +
    geom_roundrect(data = plt.data[plt.data$y == 2 & (plt.data$width >= 100000 & plt.data$width < 1000000),], ggplot2::aes(xmin = .data$start, xmax = .data$end, y = .data$y, fill = .data$target.id), rect_height = chromosome.bar.width, radius = chromosome.bar.width / 8) +
    geom_roundrect(data = plt.data[plt.data$y == 2 & plt.data$width < 100000,], ggplot2::aes(xmin = .data$start, xmax = .data$end, y = .data$y, fill = .data$target.id), rect_height = chromosome.bar.width, radius = chromosome.bar.width * 0) +
    scale_fill_manual(values = target.pal, name='Target\nsequence')

  ## Add centromere positions
  if (!is.null(target.centromeres)) {
    if (methods::is(target.centromeres, 'GRanges')) {
      if (!is.null(chromosomes)) {
        cent.gr <- GenomeInfoDb::keepSeqlevels(target.centromeres, value = chromosomes, pruning.mode = 'coarse')
      } else {
        cent.gr <- target.centromeres
      }
      cent.gr <- range(cent.gr)
    } else {
      message("Parameter 'target.centromeres' has to be valid GRanges object!!!")
    }
    ## Add centromeres to the plot
    cent.df <- as.data.frame(cent.gr)
    cent.df$target.id <- cent.df$seqnames
    plt <- plt + ggnewscale::new_scale_fill() +
      geom_roundrect(data = cent.df, ggplot2::aes(xmin = .data$start, xmax = .data$end, y = 2, fill = 'Centromere'), rect_height = chromosome.bar.width, radius = grid::unit(0, "mm"), inherit.aes = TRUE) +
      scale_fill_manual(values = c('Centromere' = 'black'), name='Annotation')
  }

  ## Set plot theme ##
  theme_genome <- ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.line.x = ggplot2::element_line(linewidth = 1),
    axis.ticks.x = ggplot2::element_line(linewidth = 1),
    axis.ticks.length.x = grid::unit(2, "mm"),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    strip.text.y.left = element_text(angle = 0)
  )
  plt <- plt + theme_genome

  ## Return final plot
  return(plt)
}
