#' Subset PAF alignments at desired genomic range.
#'
#' This function takes loaded PAF alignments using \code{\link{readPaf}} function and then subsets
#' as well as cuts PAF alignments at desired target coordinates.
#'
#' @param target.region A user defined target region either as character string ('chr:start-end') or as
#' a \code{\link{GRanges-class}} object containing a single genomic region to which PAF alignments will be
#' narrowed down.
#' @inheritParams breakPaf
#' @return A \code{tibble} of filtered PAF alignments.
#' @importFrom dplyr group_by mutate n
#' @importFrom utils read.table
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps resize start end
#' @importFrom IRanges subsetByOverlaps
#' @importFrom methods as
#' @author David Porubsky
#' @export
#' @examples
#' ## Get PAF to plot
#' paf.file <- system.file("extdata", "test1.paf", package = "SVbyEye")
#' ## Read in PAF
#' paf.table <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = "cg")
#' ## Filter PAF alignments based on desired target region
#' subsetPafAlignments(paf.table = paf.table, target.region = "target.region:19050000-19200000")
#'
subsetPafAlignments <- function(paf.table, target.region = NULL) {
    ptm <- startTimedMessage("[subsetPafAlignments] Subsetting&cutting PAF alignments")
    ## Check user input ##
    ## Make sure PAF has at least 12 mandatory fields
    if (ncol(paf.table) >= 12) {
        paf <- paf.table
    } else {
        stop("Submitted PAF alignments do not contain a minimum of 12 mandatory fields, see PAF file format definition !!!")
    }

    ## Filter alignments by target region ##
    if (!is.null(target.region)) {
        if (is.character(target.region)) {
            target.region.gr <- methods::as(target.region, "GRanges")
        } else if (is(target.region, "GRanges")) {
            target.region.gr <- target.region
        } else {
            message("Parameter 'target.region' can either be 'GRanges' object or character string 'chr#:start-end'!!!")
        }
        ## Subset PAF by ranges overlaps
        target.gr <- GenomicRanges::makeGRangesFromDataFrame(paf, seqnames.field = "t.name", start.field = "t.start", end.field = "t.end")
        hits <- GenomicRanges::findOverlaps(target.gr, target.region.gr) ## TODO Take only overlaps, regions itself are not shrinked!!!
        if (length(hits) > 0) {
            paf <- paf[S4Vectors::queryHits(hits), ]
        } else {
            stop("None of the PAF ranges overlap user defined 'target.region', exiting ...")
        }
    }

    ## Check if any paf alignment overlap user defined target region
    if (!all(target.gr == IRanges::subsetByOverlaps(target.gr, target.region.gr, type = "within"))) {
        if (nrow(paf) > 1) {
            ## Narrow alignment at desired start position ##
            cut.start.idx <- which(paf$t.start < GenomicRanges::start(target.region.gr))
            if (length(cut.start.idx) != 0) {
                ## Create PAF alignment object
                paf.aln <- paf[cut.start.idx, ]
                # aln.len <- GenomicAlignments::cigarWidthAlongReferenceSpace(cigar = paf.aln$cg)
                alignment <- GenomicAlignments::GAlignments(seqnames = paf.aln$t.name, pos = as.integer(paf.aln$t.start), cigar = paf.aln$cg, strand = GenomicRanges::strand(paf.aln$strand), names = "target")
                ## Map start position to alignment
                target.start <- GenomicRanges::resize(target.region.gr, width = 1, fix = "start")
                new.aln.start <- GenomicAlignments::mapToAlignments(x = target.start, alignments = alignment)
                ## Cut cigar string
                new.start.cigar <- GenomicAlignments::cigarNarrow(cigar = paf.aln$cg, start = start(new.aln.start), end = width(alignment))[1]
                ## Get alignment length from regional cigars
                new.start.aln.len <- sum(GenomicAlignments::explodeCigarOpLengths(cigar = new.start.cigar)[[1]])
                ## Get matched bases from regional cigars
                new.start.n.match <- sum(GenomicAlignments::explodeCigarOpLengths(cigar = new.start.cigar, ops = c("=", "M"))[[1]])
                ## Cut query coordinates
                alignment <- GenomicAlignments::GAlignments(seqnames = paf.aln$t.name, pos = as.integer(paf.aln$q.start), cigar = paf.aln$cg, strand = GenomicRanges::strand(paf.aln$strand), names = "query")
                names(new.aln.start) <- "query"
                new.query.start <- GenomicAlignments::mapFromAlignments(x = new.aln.start, alignments = alignment)
                ## Flip query ranges when alignment is in minus orientation
                if (paf.aln$strand == "-") {
                    new.query.start <- mirrorRanges(gr = new.query.start, seqlength = paf.aln$q.len)
                }
                ## Update alignment start
                paf[cut.start.idx, ]$q.start <- GenomicRanges::start(new.query.start)
                paf[cut.start.idx, ]$t.start <- GenomicRanges::start(target.start)
                paf[cut.start.idx, ]$n.match <- new.start.n.match
                paf[cut.start.idx, ]$aln.len <- new.start.aln.len
                paf[cut.start.idx, ]$cg <- new.start.cigar
            }

            ## Narrow alignment at desired end position ##
            cut.end.idx <- which(paf$t.end > GenomicRanges::end(target.region.gr))
            if (length(cut.end.idx) != 0) {
                ## Create PAF alignment object
                paf.aln <- paf[cut.end.idx, ]
                alignment <- GenomicAlignments::GAlignments(seqnames = paf.aln$t.name, pos = as.integer(paf.aln$t.start), cigar = paf.aln$cg, strand = GenomicRanges::strand(paf.aln$strand), names = "target")
                ## Map start position to alignment
                target.end <- GenomicRanges::resize(target.region.gr, width = 1, fix = "end")
                new.aln.end <- GenomicAlignments::mapToAlignments(x = target.end, alignments = alignment)
                ## Cut cigar string
                new.end.cigar <- GenomicAlignments::cigarNarrow(cigar = paf.aln$cg, start = 1, end = start(new.aln.end))[1]
                ## Get alignment length from regional cigars
                new.end.aln.len <- sum(GenomicAlignments::explodeCigarOpLengths(cigar = new.end.cigar)[[1]])
                ## Get matched bases from regional cigars
                new.end.n.match <- sum(GenomicAlignments::explodeCigarOpLengths(cigar = new.end.cigar, ops = c("=", "M"))[[1]])
                ## Cut query coordinates
                alignment <- GenomicAlignments::GAlignments(seqnames = paf.aln$t.name, pos = as.integer(paf.aln$q.start), cigar = paf.aln$cg, strand = GenomicRanges::strand(paf.aln$strand), names = "query")
                names(new.aln.end) <- "query"
                new.query.end <- GenomicAlignments::mapFromAlignments(x = new.aln.end, alignments = alignment)
                ## Flip query ranges when alignment is in minus orientation
                if (paf.aln$strand == "-") {
                    new.query.end <- mirrorRanges(gr = new.query.end, seqlength = paf.aln$q.len)
                }
                ## Update alignment end
                paf[cut.end.idx, ]$q.end <- GenomicRanges::start(new.query.end)
                paf[cut.end.idx, ]$t.end <- GenomicRanges::start(target.end)
                paf[cut.end.idx, ]$n.match <- new.end.n.match
                paf[cut.end.idx, ]$aln.len <- new.end.aln.len
                paf[cut.end.idx, ]$cg <- new.end.cigar
            }
        } else {
            ## Create PAF alignment object
            paf.aln <- paf
            alignment <- GenomicAlignments::GAlignments(seqnames = paf.aln$t.name, pos = as.integer(paf.aln$t.start), cigar = paf.aln$cg, strand = GenomicRanges::strand(paf.aln$strand), names = "target")
            ## Map target positions to alignment
            new.aln.coords <- GenomicAlignments::mapToAlignments(x = target.region.gr, alignments = alignment)
            ## Cut cigar string
            new.cigar <- GenomicAlignments::cigarNarrow(cigar = paf.aln$cg, start = start(new.aln.coords), end = end(new.aln.coords))[1]
            ## Get alignment length from regional cigars
            new.aln.len <- sum(GenomicAlignments::explodeCigarOpLengths(cigar = new.cigar)[[1]])
            ## Get matched bases from regional cigars
            new.n.match <- sum(GenomicAlignments::explodeCigarOpLengths(cigar = new.cigar, ops = c("=", "M"))[[1]])
            ## Cut query coordinates
            alignment <- GenomicAlignments::GAlignments(seqnames = paf.aln$t.name, pos = as.integer(paf.aln$q.start), cigar = paf.aln$cg, strand = GenomicRanges::strand(paf.aln$strand), names = "query")
            names(new.aln.coords) <- "query"
            new.query.coords <- GenomicAlignments::mapFromAlignments(x = new.aln.coords, alignments = alignment)
            ## Flip query ranges when alignment is in minus orientation
            if (paf.aln$strand == "-") {
                new.query.coords <- mirrorRanges(gr = new.query.coords, seqlength = paf.aln$q.len)
            }

            ## Update PAF coordinates and cigar string ##
            paf$q.start <- GenomicRanges::start(new.query.coords)
            paf$q.end <- GenomicRanges::end(new.query.coords)
            paf$t.start <- GenomicRanges::start(target.region.gr)
            paf$t.end <- GenomicRanges::end(target.region.gr)
            paf$n.match <- new.n.match
            paf$aln.len <- new.aln.len
            paf$cg <- new.cigar
        }
    }
    stopTimedMessage(ptm)
    return(paf)
}
