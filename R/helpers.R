#' Function to add bezier control points for horizontal layout.
#'
#' @param data A \code{data.frame} containing x and y coordinates.
#' @param strength The proportion to move the control point along the y-axis towards the other end of the bezier curve.
#' @return A \code{vector} of rescaled coordinate values.
#' @author David Porubsky
#'
add.control.points <- function(data = NULL, strength = 0.5) {
    start <- data[c(TRUE, FALSE), ]
    end <- data[c(FALSE, TRUE), ]
    y_diff <- (end$y - start$y) * strength
    mid1 <- start
    mid1$y <- mid1$y + y_diff
    mid2 <- end
    mid2$y <- mid2$y - y_diff
    rbind(start, mid1, mid2, end)
}

#' Function to convert coordinates between different coordinate scales.
#'
#' This function takes as input query coordinates and a query range and convert them into
#' coordinate scale defined by a target range.
#'
#' @param x A \code{vector} of coordinate values to be rescaled.
#' @param q.range  A \code{vector} containing min and max value for a query range.
#' @param t.range A \code{vector} containing min and max value for a target range.
#' @return A \code{vector} of rescaled coordinate values.
#' @author David Porubsky
#' @export
#' @examples
#' ## Convert query coordinates (x) from query range (q.range) to target range (t.range)
#' q2t(x = c(100, 1000), q.range = c(1, 2000), t.range = c(1, 100))
#'
q2t <- function(x, q.range, t.range) {
    if (is.numeric(x) && is.numeric(q.range) && is.numeric(t.range)) {
        coord.factor <- (t.range[2] - t.range[1]) / (q.range[2] - q.range[1])
        return(t.range[1] + (x - q.range[1]) * coord.factor)
    }
}

#' Mirror/reflect genomic ranges given the sequence length.
#'
#' This function flips set of ranges in a mirrored fashion.
#'
#' @param gr A \code{\link{GRanges-class}} object with one or more ranges to be mirrored/reflected given the sequence length.
#' @param seqlength  A \code{numeric} value containing the sequence length from where the submitted ranges originate from.
#' @importFrom GenomeInfoDb seqlengths seqnames
#' @importFrom GenomicRanges GRanges strand start end
#' @importFrom IRanges IRanges
#' @return A \code{\link{GRanges-class}} object with mirrored coordinates.
#' @author David Porubsky
#'
mirrorRanges <- function(gr, seqlength = NULL) {
    if (!is.null(seqlength)) {
        gr.len <- seqlength
    } else if (!is.na(GenomeInfoDb::seqlengths(gr))) {
        gr.len <- GenomeInfoDb::seqlengths(gr)
    } else {
        stop("No seglength provided!!!")
    }
    if (!all(GenomicRanges::end(gr) <= gr.len)) {
        stop("One or all submitted ranges are outside of defined seqlength!!!")
    }
    starts <- gr.len - GenomicRanges::end(gr)
    ends <- gr.len - GenomicRanges::start(gr)
    # starts <- gr.len - cumsum(width(gr))
    # ends <- starts + width(gr)
    new.gr <- GenomicRanges::GRanges(
        seqnames = GenomeInfoDb::seqnames(gr),
        ranges = IRanges::IRanges(start = starts, end = ends),
        strand = GenomicRanges::strand(gr)
    )
    suppressWarnings(GenomeInfoDb::seqlengths(new.gr) <- gr.len)
    return(new.gr)
}


#' Add color scheme based on certain values and defined breaks.
#'
#' This function takes a \code{data.frame} or \code{tibble} object and given the defined value field split this field into
#' chunks based on defined breaks. Each chunk has assigned unique color on a gradient scale.
#'
#' @param data.table A \code{data.frame} or \code{tibble} object to be processed.
#' @param value.field Either a column index or a column name present in submitted 'data.table'
#' @param breaks User defined breaks in defined 'value.field' in order to split these values into chunks
#' (Default : `5` roughly equally sized chunks).
#' @return A \code{\link{list}} containing original 'data.table' with extra column adding 'col.levels' based on defined 'breaks'.
#' List also contains element 'color' with a gradient color assigned to each 'col.level'.
#' @importFrom wesanderson wes_palette
#' @importFrom dplyr pull
#' @importFrom stats setNames
#' @importFrom ggplot2 cut_number
#' @importFrom randomcoloR randomColor
#' @author David Porubsky
#' @export
#' @examples
#' ## Define a data table with column containig 100 random values
#' data.t <- data.frame("Identity" = runif(100, min = 0, max = 1) * 100)
#' ## By default value column is divided into 5 equally sized chunks
#' ## Each chunk has assigned unique gradient color
#' colorScheme.l <- getColorScheme(data.table = data.t, value.field = "Identity")
#' ## Function returns a list containing original data table with extra color levels column
#' colorScheme.l$data
#' ## as well as corresponding color scheme
#' colorScheme.l$colors
#' ## User can also define custom breaks used to divide value.field into chunks.
#' getColorScheme(data.table = data.t, value.field = "Identity", breaks = c(0.25, 0.5, 0.75))
#'
getColorScheme <- function(data.table = NULL, value.field = NULL, breaks = NULL) {
    ## Check user input
    if (is.null(data.table)) {
        stop("No data submitted, please define 'data.table' parameter !!!")
    }
    if (is.numeric(value.field)) {
        if (ncol(data.table) < value.field) {
            stop("Defined value field index is larger than number of columns in submitted 'data.table' !!!")
        }
    } else if (is.character(value.field)) {
        if (!value.field %in% colnames(data.table)) {
            stop("Defined field name does not match any column name in submitted 'data.table' !!!")
        }
    } else {
        stop("No 'value.field' defined, please define 'value.field' either as column index or column name !!!")
    }

    ## Define continuous color scale ##
    if (is.numeric(data.table[, eval(value.field), drop = TRUE])) {
        vals <- data.table %>% dplyr::pull(eval(value.field))
        ## Define break ranges ##
        if (!is.null(breaks)) {
            levels <- c(
                paste0("<", breaks[1]),
                paste(breaks[-length(breaks)], breaks[-1], sep = ":"),
                paste0(">", breaks[length(breaks)])
            )
            ## Get break intervals
            ids <- findInterval(vals, vec = breaks) + 1
            data.table$col.levels <- factor(levels[ids], levels = levels)
            colors <- wesanderson::wes_palette(name = "Zissou1", n = length(levels), type = "continuous")
            colors <- stats::setNames(as.list(colors), levels)
        } else {
            ## If breaks are not defined split data into 5 chunks
            data.table$col.levels <- ggplot2::cut_number(vals, n = 5)
            colors <- wesanderson::wes_palette(name = "Zissou1", n = 5, type = "continuous")
            colors <- stats::setNames(as.list(colors), levels(data.table$col.levels))
        }
    }

    ## Define discrete color scale ##
    if (is.character(data.table[, eval(value.field), drop = TRUE])) {
        discrete.levels <- unique(data.table[, eval(value.field), drop = TRUE])
        n.uniq <- length(discrete.levels)
        colors <- randomcoloR::randomColor(count = n.uniq)
        colors <- stats::setNames(as.list(colors), discrete.levels)
    }
    ## Return color scheme
    return(list(data = data.table, colors = colors))
}


#' Collapse PAF query and target ranges based on unique identifier.
#'
#' This function takes PAF alignments stored in using \code{tibble} object and collapse them
#' based on grouping variable defined in `collapse.by` that have to be a valid column name
#' present in `paf.table`.
#'
#' @param collapse.by A user defined column name present in `paf.table` to serve as grouping variable.
#' @inheritParams breakPaf
#' @return A \code{tibble} of collapsed PAF alignments.
#' @importFrom dplyr group_by across all_of summarise relocate last_col
#' @author David Porubsky
#' @export
#' @examples
#' ## Get PAF to plot ##
#' paf.file <- system.file("extdata", "test1.paf", package = "SVbyEye")
#' ## Read in PAF
#' paf.table <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = "cg")
#' ## Split PAF alignments into user defined bins
#' paf.table <- pafToBins(paf.table = paf.table, binsize = 1000)
#' ## Collapse PAF alignments by bin id
#' collapsePaf(paf.table = paf.table, collapse.by = "bin.id")
#'
collapsePaf <- function(paf.table, collapse.by = NULL) {
    ## Check user input ##
    ## Make sure PAF has at least 12 mandatory fields
    if (ncol(paf.table) >= 12) {
        paf <- paf.table
    } else {
        stop("Submitted PAF alignments do not contain a minimum of 12 mandatory fields, see PAF file format definition !!!")
    }
    ## Make sure parameter collapse.by is defined and present as column in paf.table
    if (!is.null(collapse.by)) {
        if (!(nchar(collapse.by) > 0 & collapse.by %in% colnames(paf.table))) {
            stop("User defined parameter 'collapse.by' is not a valid column in submitted 'paf.table' !!!")
        }
    }
    ## Take minimum (start) and maximum (end) position for query and target coordinates grouped by 'collapse.by' variable
    paf <- paf %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(collapse.by))) %>%
        dplyr::summarise(
            q.name = paste(unique(q.name), collapse = ";"),
            q.len = paste(unique(q.len), collapse = ";"),
            q.start = min(q.start),
            q.end = max(q.end),
            strand = paste(unique(strand), collapse = ";"),
            t.name = paste(unique(t.name), collapse = ";"),
            t.len = paste(unique(t.len), collapse = ";"),
            t.start = min(t.start),
            t.end = max(t.end),
            n.match = 0, ## Try to calculate this
            aln.len = 0, ## Try to calculate this
            mapq = paste(unique(mapq), collapse = ";")
        ) %>%
        dplyr::relocate(!!collapse.by, .after = dplyr::last_col())
    ## Return collapsed PAF
    return(paf)
}


#' Assign each genomic range a non-overlapping level value
#'
#' This function takes a set of genomic ranges in \code{\link{GRanges-class}} object and assign to each range
#' a unique level such that each range or group of ranges (if `annotation.group` is defined) resides on different levels.
#'
#' @param annot.gr A \code{\link{GRanges-class}} object with a set of ranges to report unique non-overlapping levels for.
#' @param offset A \code{numeric} value from which unique annotation levels start from.
#' @param annotation.group A name of an extra field present in 'annot.gr' to be used to report unique levels for a set of
#'  ranges that belong to the same group.
#' @param direction Set to `positive` or `negative` to make sure reported levels are either decreasing or increasing
#' from the defined 'offset' (default: `positive`).
#' @return A \code{numeric} vector of unique levels corresponding to submitted ranges via 'annot.gr'.
#' @importFrom GenomicRanges mcols disjointBins GRanges
#' @importFrom IRanges ranges
#' @author David Porubsky
#' @export
#' @examples
#' ## Define random set of ranges grouped by gene name
#' test.gr1 <- GenomicRanges::GRanges(
#'     seqnames = "target.region",
#'     ranges = IRanges::IRanges(
#'         start = c(19000000, 19030000, 19070000),
#'         end = c(19010000, 19050000, 19090000)
#'     )
#' )
#' test.gr1$ID <- "gene1"
#' test.gr2 <- GenomicRanges::shift(test.gr1, shift = 10000)
#' test.gr2$ID <- "gene2"
#' test.gr3 <- GenomicRanges::shift(test.gr1, shift = 90000)
#' test.gr3$ID <- "gene3"
#' test.gr <- c(test.gr1, test.gr2, test.gr3)
#' ## Obtain unique annotation levels grouped by ID
#' getAnnotationLevels(annot.gr = test.gr, offset = 1, annotation.group = "ID")
#'
getAnnotationLevels <- function(annot.gr = NULL, offset = NULL, annotation.group = NULL, direction = "positive") {
    if (!is.null(annotation.group)) {
        if (annotation.group %in% names(GenomicRanges::mcols(annot.gr))) {
            group.seqnames <- unlist(GenomicRanges::mcols(annot.gr[, annotation.group]))
            ranges.gr <- GenomicRanges::GRanges(seqnames = group.seqnames, ranges = IRanges::ranges(annot.gr))
            ranges.gr <- range(ranges.gr)
            ranges.gr$levels <- GenomicRanges::disjointBins(IRanges::ranges(ranges.gr))
            levels <- ranges.gr$levels[match(group.seqnames, as.character(GenomeInfoDb::seqnames(ranges.gr)))]
        }
    } else {
        levels <- GenomicRanges::disjointBins(annot.gr)
    }
    ## Define offsets
    n.levels <- length(unique(levels))
    if (direction == "positive") {
        annot.levels <- seq(from = offset, by = 0.05, length.out = n.levels)
    } else if (direction == "negative") {
        annot.levels <- seq(from = offset, by = 0.05, length.out = n.levels)
    } else {
        warning("Parameter 'direction' can only take values 'positive' or 'negative', using default value !!!")
        annot.levels <- seq(from = offset, by = 0.05, length.out = n.levels)
    }
    annot.levels <- annot.levels[levels]
    return(annot.levels)
}


#' Define shift of genomic coordinates in order to get continuous coordinate scale.
#'
#' This function takes a PAF alignments loaded by \code{\link{readPaf}} function and defines how much to shift genomic coordinates
#' in order to get continuous coordinates in cases when there more than one query or target sequences in the input PAF alignments.
#'
#' @inheritParams breakPaf
#' @importFrom dplyr group_by reframe
#' @return A \code{tibble} of PAF alignments reported with a continuous coordinates.
#' @author David Porubsky
#' @export
#' @examples
#' ## Get PAF to plot
#' paf.file <- system.file("extdata", "test_ava.paf", package = "SVbyEye")
#' ## Read in PAF
#' paf.table <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = "cg")
#' ## Filter PAF
#' paf.table <- paf.table[paf.table$t.name == "HG03453_2", ]
#' ## Define shift in genomic coordiantes to get continuous scale
#' paf2continuousScale(paf.table)
#'
paf2continuousScale <- function(paf.table) {
    ## Make sure PAF has at least 12 mandatory fields
    if (ncol(paf.table) >= 12) {
        paf <- paf.table
    } else {
        stop("Submitted PAF alignments do not contain a minimum of 12 mandatory fields, see PAF file format definition !!!")
    }

    ## Check if there more than one unique query sequences
    if (length(unique(paf$q.name)) > 1) {
        ## If yes make sure query coordinates are continuous
        q.limits <- paf %>%
            dplyr::arrange(t.start) %>% ## Make sure query coordinates are sorted by target coordinates
            dplyr::mutate(q.name = factor(q.name, levels = unique(q.name))) %>%
            dplyr::group_by(q.name) %>%
            dplyr::reframe(range = range(c(q.start, q.end)))
        q.gaps <- diff(q.limits$range) * -1
        q.shift <- c(0, q.gaps[seq(length(q.gaps)) %% 2 == 0] + 1)
        q.shift <- cumsum(q.shift)
        names(q.shift) <- unique(q.limits$q.name)
        # paf$q.genomic.start <- paf$q.start
        # paf$q.genomic.end <- paf$q.end
        # paf$q.start <- paf$q.start + q.shift[paf$q.name]
        # paf$q.end <- paf$q.end + q.shift[paf$q.name]
        paf$q.shift <- q.shift[paf$q.name]
    } else {
        paf$q.shift <- 0
        # paf$q.genomic.start <- paf$q.start
        # paf$q.genomic.end <- paf$q.end
    }

    ## Check if there more than one unique target sequences
    if (length(unique(paf$t.name)) > 1) {
        ## If yes make sure target coordinates are continuous
        t.limits <- paf %>%
            dplyr::arrange(q.start) %>% ## Make sure query coordinates are sorted by query coordinates
            dplyr::mutate(t.name = factor(t.name, levels = unique(t.name))) %>%
            dplyr::group_by(t.name) %>%
            dplyr::reframe(range = range(c(t.start, t.end)))
        t.gaps <- diff(t.limits$range) * -1
        t.shift <- c(0, t.gaps[seq(length(t.gaps)) %% 2 == 0] + 1)
        t.shift <- cumsum(t.shift)
        names(t.shift) <- unique(t.limits$t.name)
        # paf$t.genomic.start <- paf$t.start
        # paf$t.genomic.end <- paf$t.end
        # paf$t.start <- paf$t.start + t.shift[paf$t.name]
        # paf$t.end <- paf$t.end + t.shift[paf$t.name]
        paf$t.shift <- t.shift[paf$t.name]
    } else {
        paf$t.shift <- 0
        # paf$t.genomic.start <- paf$t.start
        # paf$t.genomic.end <- paf$t.end
    }
    return(paf)
}
