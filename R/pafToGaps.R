#' Report gaps between PAF alignments.
#'
#' This function takes loaded PAF alignments using \code{\link{readPaf}} function and then reports all gaps
#' in subsequent target and query alignments and classify them as either deletions ('D') or insertions ('I').
#'
#' @param min.gap.diff A user defined minimum gap size difference in target and query coordinates (and vice versa)
#' to be retained.
#' @inheritParams breakPaf
#' @return A \code{tibble} of PAF alignment gaps.
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom dplyr bind_rows
#' @author David Porubsky
#' @export
#' @examples
#' ## Get PAF to process ##
#' paf.file <- system.file("extdata", "test1.paf", package = "SVbyEye")
#' ## Read in PAF
#' paf.table <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = "cg")
#' ## Report gaps between PAF alignments
#' paf2gaps(paf.table = paf.table)
#'
paf2gaps <- function(paf.table, min.gap.diff = NULL) {
    ## Check user input ##
    ## Make sure submitted paf.table has at least 12 mandatory fields
    if (ncol(paf.table) >= 12) {
        paf <- paf.table
    } else {
        stop("Submitted PAF alignments do not contain a minimum of 12 mandatory fields, see PAF file format definition !!!")
    }

    ## Process by unique query-target pairs ##
    seq.pair <- paste0(paf$q.name, "__", paf$t.name)
    paf.l <- split(paf, seq.pair)

    all.gaps <- list()
    for (i in seq_along(paf.l)) {
        sub.paf <- paf.l[[i]]
        ## Convert paf.table to Genomic ranges object
        t.ranges <- GenomicRanges::GRanges(
            seqnames = sub.paf$t.name,
            ranges = IRanges::IRanges(start = sub.paf$t.start, end = sub.paf$t.end)
        )
        q.ranges <- GenomicRanges::GRanges(
            seqnames = sub.paf$q.name,
            ranges = IRanges::IRanges(start = sub.paf$q.start, end = sub.paf$q.end)
        )

        ## Get gaps in respect to target ranges
        t.gaps <- reportGaps(t.ranges = t.ranges, q.ranges = q.ranges)
        ## Get gaps in respect to query ranges
        q.gaps <- reportGaps(t.ranges = t.ranges, q.ranges = q.ranges, invert = TRUE)
        ## Concatenate target and query gaps
        t.gaps.gr <- c(t.gaps[, 0], q.gaps$aln.gap.gr[, 0])
        q.gaps.gr <- c(t.gaps$aln.gap.gr[, 0], q.gaps[, 0])

        ## Report gaps in PAF format ##
        paf.gaps <- tibble::tibble(
            q.name = unique(sub.paf$q.name),
            q.len = unique(sub.paf$q.len),
            q.start = GenomicRanges::start(q.gaps.gr),
            q.end = GenomicRanges::end(q.gaps.gr),
            strand = "*",
            t.name = unique(sub.paf$t.name),
            t.len = unique(sub.paf$t.len),
            t.start = GenomicRanges::start(t.gaps.gr),
            t.end = GenomicRanges::end(t.gaps.gr),
            n.match = 0,
            aln.len = 0,
            mapq = 0
        )

        ## Define CIGAR string ##
        ## Get size difference between target and corresponding alignment gap width
        size.diff <- abs(GenomicRanges::width(t.gaps.gr) - GenomicRanges::width(q.gaps.gr))
        ## Mark the gap size difference as either delection ('D') or insertion ('I')
        cg.label <- ifelse(GenomicRanges::width(t.gaps.gr) > GenomicRanges::width(q.gaps.gr), "D", "I")
        paf.gaps$cg <- paste0(size.diff, cg.label)
        ## Store gaps per unique sequence pair
        all.gaps[[i]] <- paf.gaps
    }
    all.gaps <- dplyr::bind_rows(all.gaps)

    ## Filter by gap size difference
    if (!is.null(min.gap.diff)) {
        if (min.gap.diff > 0) {
            gap.diff.len <- as.numeric(gsub(all.gaps$cg, pattern = "D|I", replacement = ""))
            all.gaps <- all.gaps[gap.diff.len >= min.gap.diff, ]
        }
    }
    ## Return PAF gaps
    if (nrow(all.gaps) > 0) {
        return(all.gaps)
    } else {
        NULL
    }
}


#' Report gaps between a set of target and query ranges.
#'
#' This function takes set of target and query ranges and reports all gaps in target ranges (Default behavior)
#' and their corresponding query ranges. This behavior is reversed using `invert` parameter.
#' Note that target and query ranges are expected to be 1-to-1 alignments as reported in a PAF format.
#'
#' @param t.ranges A \code{\link{GRanges-class}} object of target ranges where gaps are searched for by default.
#' @param q.ranges A \code{\link{GRanges-class}} object of query ranges.
#' @param invert If set to \code{TRUE} `q.ranges` are used to search for gaps instead of `t.ranges`.
#' @return A \code{\link{GRanges-class}} object reporting gaps as well as upstream (up) and downstream (down)
#' alignments around each gap.
#' @import GenomicRanges
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges IRanges follow precede
#' @importFrom tibble tibble
#' @importFrom dplyr select all_of
#' @author David Porubsky
#' @export
#' @examples
#' ## Define target and query ranges to process
#' target.ranges <- GenomicRanges::GRanges(
#'     seqnames = "target",
#'     IRanges::IRanges(
#'         start = c(1, 111, 150),
#'         end = c(100, 150, 200)
#'     )
#' )
#' query.ranges <- GenomicRanges::GRanges(
#'     seqnames = "query",
#'     IRanges::IRanges(
#'         start = c(1000, 1100, 1301),
#'         end = c(1100, 1200, 1500)
#'     )
#' )
#' ## Report gaps within target ranges
#' reportGaps(t.ranges = target.ranges, q.ranges = query.ranges)
#' ## Report gaps within query ranges
#' reportGaps(t.ranges = target.ranges, q.ranges = query.ranges, invert = TRUE)
#'
reportGaps <- function(t.ranges = NULL, q.ranges = NULL, invert = FALSE) {
    ## Helper function ##
    getGaps <- function(gr = NULL) {
        gap.gr <- GenomicRanges::gaps(gr, start = min(GenomicRanges::start(gr)))
        gap.gr <- gap.gr[GenomicRanges::strand(gap.gr) == "*"]
        return(gap.gr)
    }

    ## Make sure that both t.ranges and q.ranges are both GRanges objects
    if (!(is(t.ranges, "GRanges") & is(q.ranges, "GRanges"))) {
        stop("Both parameters 't.ranges' and 'q.ranges' have to be valid GRanges class objects !!!")
    }
    ## Make sure that both t.ranges and q.ranges are of the same length
    if (length(t.ranges) != length(q.ranges)) {
        stop("Both parameters 't.ranges' and 'q.ranges' have to be of the same length !!!")
    }

    ## Construct object for gap analysis
    if (invert) {
        gr <- q.ranges
        gr$aln.gr <- t.ranges
    } else {
        gr <- t.ranges
        gr$aln.gr <- q.ranges
    }

    ## Sort ranges by position
    GenomicRanges::strand(gr) <- "*"
    gr <- GenomicRanges::sort(gr)
    ## Make sure rows are not named
    names(gr) <- NULL

    if (length(gr) > 1) {
        ## Obtain gap ranges ##
        gap.gr <- getGaps(gr)

        ## Get alignments on each side of the gap ##
        if (length(gap.gr) > 0) {
            ## Report upstream and downstream target ranges
            ## Upstream
            up.idx <- IRanges::follow(gap.gr, gr)
            up.idx.keep <- !is.na(up.idx)
            gap.gr$up.gr <- GenomicRanges::GRanges(
                seqnames = GenomeInfoDb::seqnames(gap.gr),
                ranges = IRanges::IRanges(
                    start = GenomicRanges::start(gap.gr),
                    end = GenomicRanges::start(gap.gr)
                )
            )
            gap.gr$up.gr[up.idx.keep] <- gr[up.idx[up.idx.keep]][, 0]
            ## Downstream
            down.idx <- IRanges::precede(gap.gr, gr)
            down.idx.keep <- !is.na(down.idx)
            gap.gr$down.gr <- GenomicRanges::GRanges(
                seqnames = GenomeInfoDb::seqnames(gap.gr),
                ranges = IRanges::IRanges(
                    start = GenomicRanges::end(gap.gr),
                    end = GenomicRanges::end(gap.gr)
                )
            )
            gap.gr$down.gr[down.idx.keep] <- gr[down.idx[down.idx.keep]][, 0]

            ## Report upstream and downstream ranges and gaps for corresponding alignments (aln.gr)
            ## Upstream
            gap.gr$aln.up.gr <- GenomicRanges::GRanges(
                seqnames = GenomeInfoDb::seqnames(gap.gr),
                ranges = IRanges::IRanges(
                    start = GenomicRanges::start(gap.gr),
                    end = GenomicRanges::start(gap.gr)
                )
            )
            suppressWarnings(gap.gr$aln.up.gr[up.idx.keep] <- gr$aln.gr[up.idx[up.idx.keep]][, 0])
            ## Downstream
            gap.gr$aln.down.gr <- GenomicRanges::GRanges(
                seqnames = GenomeInfoDb::seqnames(gap.gr),
                ranges = IRanges::IRanges(
                    start = GenomicRanges::end(gap.gr),
                    end = GenomicRanges::end(gap.gr)
                )
            )
            suppressWarnings(gap.gr$aln.down.gr[down.idx.keep] <- gr$aln.gr[down.idx[down.idx.keep]][, 0])

            ## Define gaps in query coordinates
            starts <- pmin(GenomicRanges::end(gap.gr$aln.up.gr), GenomicRanges::start(gap.gr$aln.down.gr))
            ends <- pmax(GenomicRanges::end(gap.gr$aln.up.gr), GenomicRanges::start(gap.gr$aln.down.gr))
            aln.gap.gr <- GenomicRanges::GRanges(
                seqnames = GenomeInfoDb::seqnames(gap.gr$aln.up.gr),
                ranges = IRanges::IRanges(start = starts, end = ends)
            )
            ## Set ranges with equal start and end to width zero
            aln.gap.gr[GenomicRanges::start(aln.gap.gr) == GenomicRanges::end(aln.gap.gr)] <-
                GenomicRanges::resize(aln.gap.gr[GenomicRanges::start(aln.gap.gr) == GenomicRanges::end(aln.gap.gr)], width = 0)

            ## Return gaps
            gap.gr$aln.gap.gr <- aln.gap.gr
            return(gap.gr)
        } else {
            return(NULL)
        }
    } else {
        return(NULL)
    }
}
