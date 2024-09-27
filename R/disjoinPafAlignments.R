#' Disjoin overlapping PAF alignments.
#'
#' This function takes loaded PAF alignments using \code{\link{readPaf}} function and identify overlapping
#' genomic ranges either in target or query coordinates and then split PAF alignments at the positions of
#' these overlaps into a disjoined set of genomic ranges. Alternatively, user defined set of genomic ranges
#' in parameter `disjoin.gr` can be used to disjoin PAF alignments.
#'
#' @param min.overlap A minimum number of overlapping base pairs to disjoin PAF alignments.
#' @param coordinates PAF coordinates to disjoin, either 'target' or 'query'.
#' @param disjoin.gr A \code{\link{GRanges-class}} object containing genomic ranges where PAF alignments will be disjoined
#' (must be in defined `coordinates`).
#' @inheritParams breakPaf
#' @return A \code{tibble} of disjoined PAF alignments.
#' @import GenomicRanges
#' @importFrom GenomicAlignments GAlignments cigarNarrow explodeCigarOpLengths
#' @importFrom GenomeInfoDb seqnames
#' @importFrom S4Vectors subjectHits queryHits
#' @importFrom IRanges IRanges
#' @importFrom methods is
#' @importFrom tibble tibble
#' @author David Porubsky
#' @export
#' @examples
#' ## Get PAF to plot
#' paf.file <- system.file("extdata", "test1.paf", package = "SVbyEye")
#' ## Read in PAF
#' paf.table <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = "cg")
#' ## Disjoin PAF at target coordinates
#' disj.paf.table <- disjoinPafAlignments(paf.table = paf.table, coordinates = "target")
#' ## Disjoin PAF at user defined coordinates
#' disj.gr <- GenomicRanges::GRanges(
#'     seqnames = "query.region",
#'     ranges = IRanges::IRanges(start = 16300000, end = 16350000)
#' )
#' disj.paf.table <- disjoinPafAlignments(
#'     paf.table = paf.table,
#'     coordinates = "query", disjoin.gr = disj.gr
#' )
#'
disjoinPafAlignments <- function(paf.table, min.overlap = 1000, coordinates = "target", disjoin.gr = NULL) {
    ptm <- startTimedMessage(paste0("[disjoinPafAlignments] Disjoining ", coordinates, " alignments"))

    ## Check user input ##
    ## Make sure PAF has at least 12 mandatory fields
    if (ncol(paf.table) >= 12) {
        paf <- paf.table
    } else {
        stop("Submitted PAF alignments do not contain a minimum of 12 mandatory fields, see PAF file format definition !!!")
    }

    ## Process target ranges ##
    if (coordinates == "target") {
        target.gr <- GenomicRanges::makeGRangesFromDataFrame(paf, seqnames.field = "t.name", start.field = "t.start", end.field = "t.end")

        if (methods::is(disjoin.gr, "GRanges")) {
            hits <- GenomicRanges::findOverlaps(disjoin.gr, target.gr)
            if (length(hits) > 0) {
                todisj.gr <- disjoin.gr
            } else {
                message("No genomic ranges in 'disjoin.gr' overlaps target coordinates")
            }
        } else {
            ## Get regions covered more than once
            cov.gr <- as(GenomicRanges::coverage(target.gr), "GRanges")
            cov.gr <- cov.gr[cov.gr$score > 1]
            ## Filter by minimum overlapping regions
            todisj.gr <- cov.gr[GenomicRanges::width(cov.gr) >= min.overlap, 0]
        }

        ## Get disjoined set of query and target ranges
        if (length(todisj.gr) > 0) {
            ## Remove regions to disjoin from target
            subtr.gr <- suppressWarnings(unlist(GenomicRanges::subtract(target.gr, todisj.gr, ignore.strand = TRUE)))
            ## Get intersection between regions to disjoin and target
            inters.gr <- suppressWarnings(GenomicRanges::intersect(target.gr, todisj.gr, ignore.strand = TRUE))
            target.gr <- suppressWarnings(sort(c(subtr.gr, inters.gr), ignore.strand = TRUE))
            ## Lift target ranges to query coordinates
            query.gr <- suppressMessages(liftRangesToAlignment(paf.table = paf, gr = target.gr, direction = "target2query", report.cigar.str = TRUE))
            ## Expand target ranges to match mapped query ranges
            target.gr <- target.gr[query.gr$idx]
            target.gr$alignmentsHits <- query.gr$alignmentsHits
            GenomicRanges::strand(target.gr) <- GenomicRanges::strand(query.gr)
            cigars.region <- query.gr$cg
            processCIGAR <- TRUE
        } else {
            processCIGAR <- FALSE
        }
        ## Process query ranges ##
    } else if (coordinates == "query") {
        query.gr <- GenomicRanges::makeGRangesFromDataFrame(paf, seqnames.field = "q.name", start.field = "q.start", end.field = "q.end")

        if (methods::is(disjoin.gr, "GRanges")) {
            hits <- GenomicRanges::findOverlaps(disjoin.gr, query.gr)
            if (length(hits) > 0) {
                todisj.gr <- disjoin.gr
            } else {
                message("No genomic ranges in 'disjoin.gr' overlaps target coordinates")
            }
        } else {
            ## Get regions covered more than once
            cov.gr <- as(GenomicRanges::coverage(query.gr), "GRanges")
            cov.gr <- cov.gr[cov.gr$score > 1]
            ## Filter by minimum overlapping regions
            todisj.gr <- cov.gr[GenomicRanges::width(cov.gr) >= min.overlap, 0]
        }

        ## Get disjoined set of query and target ranges
        if (length(todisj.gr) > 0) {
            ## Remove regions to disjoin from target
            subtr.gr <- suppressWarnings(unlist(GenomicRanges::subtract(query.gr, todisj.gr, ignore.strand = TRUE)))
            ## Get intersection between regions to disjoin and target
            inters.gr <- suppressWarnings(GenomicRanges::intersect(query.gr, todisj.gr, ignore.strand = TRUE))
            query.gr <- suppressWarnings(sort(c(subtr.gr, inters.gr), ignore.strand = TRUE))
            ## Lift query ranges to target coordinates
            target.gr <- suppressMessages(liftRangesToAlignment(paf.table = paf, gr = query.gr, direction = "query2target", report.cigar.str = TRUE))
            ## Expand query ranges to match mapped target ranges
            query.gr <- query.gr[target.gr$idx]
            GenomicRanges::strand(query.gr) <- GenomicRanges::strand(target.gr)
            cigars.region <- target.gr$cg
            processCIGAR <- TRUE
        } else {
            processCIGAR <- FALSE
        }
    }

    if (processCIGAR) {
        #   ## Subset the cigar string per target region ##
        # ## Define PAF alignment object
        # alignments <- GenomicAlignments::GAlignments(
        #   seqnames = paf$t.name,
        #   pos = as.integer(paf$t.start) + 1L,
        #   #pos = as.integer(paf$t.start),
        #   cigar = paf$cg,
        #   strand = GenomicRanges::strand(paf$strand),
        #   names = paf$t.name
        # )
        # ## Get overlaps between target ranges and alignments
        # #hits <- findOverlaps(target.gr, alignments, minoverlap = min.overlap)
        # hits <- GenomicRanges::findOverlaps(target.gr, alignments)
        #   end(target.gr[which.max(end(target.gr))]) <- max(end(alignments))
        #   hits <- GenomicRanges::findOverlaps(target.gr, alignments, type = 'within')
        #
        #   ## Keep only longest overlaps for each target range
        #   intersect.gr <- GenomicRanges::pintersect(target.gr[S4Vectors::queryHits(hits)], as(alignments[S4Vectors::subjectHits(hits)], 'GRanges'))
        #   overlap.width <- GenomicRanges::width(intersect.gr)
        #   max.overlap <- sapply(split(overlap.width, S4Vectors::queryHits(hits)), max)
        #   hits <- hits[overlap.width %in% max.overlap]
        #   #hits <- hits[which(unique(overlap.width) %in% max.overlap)]
        #
        #   ## Get local alignment coordinates [in parallel fashion]
        #   #lifted.gr <- GenomicAlignments::pmapToAlignments(x = target.gr, alignments = alignments[S4Vectors::subjectHits(hits)])
        #   #names(lifted.gr) <- NULL
        #
        #   ## Set target range to local alignment coordinates
        #   suppressWarnings( target.local.gr <- GenomicRanges::shift(target.gr, shift = -(start(alignments)[S4Vectors::subjectHits(hits)] - 1)) )
        #   suppressWarnings( target.local.gr <- GenomicRanges::shift(target.gr, shift = -(start(alignments)[target.gr$alignmentsHits] - 1)) )
        #   start(target.local.gr[start(target.local.gr) == 0]) <- 1
        #   #target.local.gr <- GenomicRanges::trim(target.local.gr)
        #
        #   ## Define start and end position for CIGAR extraction
        #   #starts <- pmax(GenomicRanges::start(lifted.gr), 1)
        #   #width <- GenomicRanges::width(target.gr) - 1
        #   #ends <- pmin(GenomicRanges::end(lifted.gr), GenomicRanges::width(alignments)[S4Vectors::subjectHits(hits)])
        #   starts <- GenomicRanges::start(target.local.gr)
        #   width <- GenomicRanges::width(target.local.gr)
        #
        #   cigars <- alignments@cigar[S4Vectors::subjectHits(hits)]
        #   cigars.region <- vapply(seq_along(starts), function(i) {
        #     tryCatch(
        #       {
        #         GenomicAlignments::cigarNarrow(cigar = cigars[i], start = starts[i], width = width[i])[1]
        #       },
        #       error = function(e) {
        #         return("1=")
        #       }
        #     )
        #   }, FUN.VALUE = character(1))

        ## Get alignment length from regional cigars
        aln.len <- vapply(cigars.region, function(cg) sum(GenomicAlignments::explodeCigarOpLengths(cigar = cg)[[1]]), USE.NAMES = FALSE, FUN.VALUE = numeric(1))
        ## Get matched bases from regional cigars
        n.match <- vapply(cigars.region, function(cg) sum(GenomicAlignments::explodeCigarOpLengths(cigar = cg, ops = c("=", "M"))[[1]]), USE.NAMES = FALSE, FUN.VALUE = numeric(1))
        ## Get single-base mismatches from regional cigars
        n.mismatch <- vapply(cigars.region, function(cg) sum(GenomicAlignments::explodeCigarOpLengths(cigar = cg, ops = c("X"))[[1]]), USE.NAMES = FALSE, FUN.VALUE = numeric(1))

        ## Prepare modified PAF for export ##
        # return.paf <- tibble::tibble("q.name" = as.character(GenomeInfoDb::seqnames(query.gr))
        #                           ,"q.len" = paf$q.len[S4Vectors::subjectHits(hits)]
        #                           ,"q.start" = GenomicRanges::start(query.gr)
        #                           ,"q.end" = GenomicRanges::end(query.gr)
        #                           ,"strand" = as.character( GenomicRanges::strand(query.gr) )
        #                           ,"t.name" = as.character(GenomeInfoDb::seqnames(target.gr))
        #                           ,"t.len" = paf$t.len[S4Vectors::subjectHits(hits)]
        #                           ,"t.start" = GenomicRanges::start(target.gr)
        #                           ,"t.end" = GenomicRanges::end(target.gr)
        #                           ,"n.match" = n.match
        #                           ,"aln.len" = aln.len
        #                           ,"mapq" = paf$mapq[S4Vectors::subjectHits(hits)]
        #                           ,"cg" = cigars.region)

        return.paf <- tibble::tibble(
            "q.name" = as.character(GenomeInfoDb::seqnames(query.gr)),
            "q.len" = paf$q.len[target.gr$alignmentsHits],
            "q.start" = GenomicRanges::start(query.gr),
            "q.end" = GenomicRanges::end(query.gr),
            "strand" = as.character(GenomicRanges::strand(query.gr)),
            "t.name" = as.character(GenomeInfoDb::seqnames(target.gr)),
            "t.len" = paf$t.len[target.gr$alignmentsHits],
            "t.start" = GenomicRanges::start(target.gr),
            "t.end" = GenomicRanges::end(target.gr),
            "n.match" = n.match,
            "aln.len" = aln.len,
            "mapq" = paf$mapq[target.gr$alignmentsHits],
            "cg" = cigars.region
        )
    } else {
        return.paf <- paf
        message("    No genomic ranges to disjoin were found !!!")
    }

    stopTimedMessage(ptm)
    ## Return PAF
    return(return.paf)
}
