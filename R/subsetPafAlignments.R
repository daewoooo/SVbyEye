#' Subset PAF alignments at desired genomic range.
#'
#' This function takes loaded PAF alignments using \code{\link{readPaf}} function and then subsets
#' as well as cuts PAF alignments at desired target coordinates. This function can only be applied
#' to PAF alignments containing a single query and target sequence.
#'
#' @param target.region A user defined target region either as character string ('chr:start-end') or as
#' a \code{\link{GRanges-class}} object containing a single genomic region to which PAF alignments will be
#' narrowed down.
#' @inheritParams breakPaf
#' @return A \code{tibble} of subsetted PAF alignments.
#' @importFrom GenomicRanges makeGRangesFromDataFrame resize start end width strand
#' @importFrom GenomicAlignments cigarNarrow GAlignments explodeCigarOpLengths
#' @importFrom IRanges subsetByOverlaps
#' @importFrom tibble tibble
#' @importFrom methods as
#' @importFrom dplyr bind_rows
#' @author David Porubsky
#' @export
#' @examples
#' ## Get PAF to plot
#' paf.file <- system.file("extdata", "test1.paf", package = "SVbyEye")
#' ## Read in PAF
#' paf.table <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = "cg")
#' ## Subset PAF alignments based on desired target region
#' subsetPafAlignments(paf.table = paf.table, target.region = "target.region:19050000-19200000")
#'
subsetPafAlignments <- function(paf.table, target.region = NULL) {
    ptm <- startTimedMessage("[subsetPafAlignments] Subsetting&cutting PAF alignments")
    ## Check user input ##
    ## Make sure PAF has at least 12 mandatory fields
    if (ncol(paf.table) >= 12) {
        paf <- paf.table
    } else {
        stop("\nSubmitted PAF alignments do not contain a minimum of 12 mandatory fields, see PAF file format definition !!!")
    }
    ## Check if PAF contains expected cg field containing CIGAR string
    if (!"cg" %in% colnames(paf)) {
        stop(paste0(
            "\nExpected PAF field 'cg' containing CIGAR string is missing, ",
            "therefore alignments cannot be narrowed down exactly to the user defined 'target.region'!!!"
        ))
    }

    ## Filter alignments by target region ##
    if (!is.null(target.region)) {
        if (is.character(target.region)) {
            target.region.gr <- methods::as(target.region, "GRanges")
        } else if (is(target.region, "GRanges")) {
            target.region.gr <- target.region
        } else {
            message("\nParameter 'target.region' can either be 'GRanges' object or character string 'chr#:start-end'!!!")
        }
        ## Subset PAF by ranges overlaps
        target.gr <- GenomicRanges::makeGRangesFromDataFrame(paf, seqnames.field = "t.name", start.field = "t.start", end.field = "t.end")
        hits <- IRanges::findOverlaps(target.gr, target.region.gr) ## TODO Take only overlaps, regions itself are not shrinked!!!
    }

    ## Return NULL if there are no overlaps with the target region
    if (length(hits) > 0) {
        paf <- paf[S4Vectors::queryHits(hits), ]
        target.gr <- target.gr[S4Vectors::queryHits(hits), ]

        ## Get alignments embedded within target region
        hits <- IRanges::findOverlaps(target.gr, target.region.gr, type = "within", ignore.strand = TRUE)
        within.idx <- S4Vectors::queryHits(hits)
        ## Get alignments with target region being embedded within
        hits <- IRanges::findOverlaps(target.region.gr, target.gr, type = "within", ignore.strand = TRUE)
        embedded.idx <- S4Vectors::subjectHits(hits)
        ## Get alignments with target region being partially overlapping
        partial.idx <- setdiff(seq_along(target.gr), c(embedded.idx, within.idx))

        ## Get alignments embedded within target region
        if (length(within.idx) > 0) {
            paf.within.aln <- paf[within.idx, ]
        } else {
            paf.within.aln <- tibble::tibble()
        }

        ## Process alignments with target region being embedded within
        if (length(embedded.idx) > 0) {
            ## Create PAF alignment object
            paf.aln <- paf[embedded.idx, ]
            alignments <- GenomicAlignments::GAlignments(
                seqnames = paf.aln$t.name,
                pos = as.integer(paf.aln$t.start) + 1L,
                cigar = paf.aln$cg,
                strand = GenomicRanges::strand(paf.aln$strand),
                names = paf.aln$t.name
            )
            ## Get overlaps between target ranges and alignments
            hits <- IRanges::findOverlaps(target.region.gr, alignments)
            ## Lift target positions to query coordinates
            query.coords <- suppressMessages(liftRangesToAlignment(paf.table = paf.aln, gr = target.region.gr, direction = "target2query", report.cigar.str = TRUE))
            ## Set target range to local alignment coordinates
            suppressWarnings(target.local.gr <- GenomicRanges::shift(target.region.gr[S4Vectors::queryHits(hits)],
                shift = -(GenomicRanges::start(alignments)[S4Vectors::subjectHits(hits)] - 1)
            ))
            GenomicRanges::start(target.local.gr[GenomicRanges::start(target.local.gr) == 0]) <- 1

            ## Cut cigar string
            starts <- GenomicRanges::start(target.local.gr)
            width <- GenomicRanges::width(target.local.gr)
            cigars <- alignments@cigar[S4Vectors::subjectHits(hits)]
            new.cigar <- vapply(seq_along(starts), function(i) {
                tryCatch(
                    {
                        GenomicAlignments::cigarNarrow(cigar = cigars[i], start = starts[i], width = width[i])[1]
                    },
                    error = function(e) {
                        return("1=")
                    }
                )
            }, FUN.VALUE = character(1))
            ## Get alignment length from regional cigars
            new.aln.len <- sum(GenomicAlignments::explodeCigarOpLengths(cigar = new.cigar)[[1]])
            ## Get matched bases from regional cigars
            new.n.match <- sum(GenomicAlignments::explodeCigarOpLengths(cigar = new.cigar, ops = c("=", "M"))[[1]])

            ## Update PAF coordinates and cigar string ##
            paf.aln$q.start <- GenomicRanges::start(query.coords)
            paf.aln$q.end <- GenomicRanges::end(query.coords)
            paf.aln$t.start <- GenomicRanges::start(target.region.gr)
            paf.aln$t.end <- GenomicRanges::end(target.region.gr)
            paf.aln$n.match <- new.n.match
            paf.aln$aln.len <- new.aln.len
            paf.aln$cg <- new.cigar
            paf.embedded.aln <- paf.aln
        } else {
            paf.embedded.aln <- tibble::tibble()
        }

        ## Process alignments with target region being partially overlapping
        if (length(partial.idx) > 0) {
            paf.partial.aln <- paf[partial.idx, ]
            ## Narrow alignment at desired start position ##
            cut.start.idx <- which(paf.partial.aln$t.start < GenomicRanges::start(target.region.gr))
            if (length(cut.start.idx) != 0) {
                ## Create PAF alignment object
                paf.aln <- paf.partial.aln[cut.start.idx, ]
                alignments <- GenomicAlignments::GAlignments(
                    seqnames = paf.aln$t.name,
                    pos = as.integer(paf.aln$t.start) + 1L,
                    cigar = paf.aln$cg,
                    strand = GenomicRanges::strand(paf.aln$strand),
                    names = paf.aln$t.name
                )
                ## Get overlaps between target ranges and alignments
                hits <- GenomicRanges::findOverlaps(target.region.gr, alignments)
                ## Lift target positions to query coordinates
                target.start <- GenomicRanges::resize(target.region.gr, width = 1, fix = "start")
                query.start <- suppressMessages(liftRangesToAlignment(paf.table = paf.aln, gr = target.start, direction = "target2query", report.cigar.str = TRUE))
                ## Set target range to local alignment coordinates
                suppressWarnings(target.local.gr <- GenomicRanges::shift(target.start[S4Vectors::queryHits(hits)],
                    shift = -(start(alignments)[S4Vectors::subjectHits(hits)] - 1)
                ))
                GenomicRanges::start(target.local.gr[GenomicRanges::start(target.local.gr) == 0]) <- 1

                ## Cut cigar string
                starts <- GenomicRanges::start(target.local.gr)
                cigars <- alignments@cigar[S4Vectors::subjectHits(hits)]
                widths <- GenomicRanges::width(alignments[S4Vectors::subjectHits(hits)])
                new.cigar <- vapply(seq_along(starts), function(i) {
                    tryCatch(
                        {
                            GenomicAlignments::cigarNarrow(cigar = cigars[i], start = starts[i], end = widths[i])[1]
                        },
                        error = function(e) {
                            return("1=")
                        }
                    )
                }, FUN.VALUE = character(1))
                ## Get alignment length from regional cigars
                new.aln.len <- vapply(GenomicAlignments::explodeCigarOpLengths(cigar = new.cigar), FUN = sum, FUN.VALUE = numeric(1))
                ## Get matched bases from regional cigars
                new.n.match <- vapply(GenomicAlignments::explodeCigarOpLengths(cigar = new.cigar, ops = c("=", "M")), FUN = sum, FUN.VALUE = numeric(1))

                ## Update PAF coordinates and cigar string ##
                paf.aln$q.end[paf.aln$strand == "-"] <- GenomicRanges::start(query.start[GenomicRanges::strand(query.start) == "-"])
                paf.aln$q.start[paf.aln$strand == "+"] <- GenomicRanges::start(query.start[GenomicRanges::strand(query.start) == "+"])
                paf.aln$t.start <- GenomicRanges::start(target.region.gr)
                # paf.aln$t.end <- GenomicRanges::end(target.region.gr)
                paf.aln$n.match <- new.n.match
                paf.aln$aln.len <- new.aln.len
                paf.aln$cg <- new.cigar
                paf.partial.aln[cut.start.idx, ] <- paf.aln
            }

            ## Narrow alignment at desired end position ##
            cut.end.idx <- which(paf.partial.aln$t.end > GenomicRanges::end(target.region.gr))
            if (length(cut.end.idx) != 0) {
                ## Create PAF alignment object
                paf.aln <- paf.partial.aln[cut.end.idx, ]
                alignments <- GenomicAlignments::GAlignments(
                    seqnames = paf.aln$t.name,
                    pos = as.integer(paf.aln$t.start) + 1L,
                    cigar = paf.aln$cg,
                    strand = GenomicRanges::strand(paf.aln$strand),
                    names = paf.aln$t.name
                )
                ## Get overlaps between target ranges and alignments
                hits <- GenomicRanges::findOverlaps(target.region.gr, alignments)
                ## Lift target positions to query coordinates
                target.end <- GenomicRanges::resize(target.region.gr, width = 1, fix = "end")
                query.end <- suppressMessages(liftRangesToAlignment(paf.table = paf.aln, gr = target.end, direction = "target2query", report.cigar.str = TRUE))
                ## Set target range to local alignment coordinates
                suppressWarnings(target.local.gr <- GenomicRanges::shift(target.end[S4Vectors::queryHits(hits)],
                    shift = -(start(alignments)[S4Vectors::subjectHits(hits)] - 1)
                ))
                GenomicRanges::start(target.local.gr[GenomicRanges::start(target.local.gr) == 0]) <- 1

                ## Cut cigar string
                ends <- GenomicRanges::start(target.local.gr)
                cigars <- alignments@cigar[S4Vectors::subjectHits(hits)]
                new.cigar <- vapply(seq_along(ends), function(i) {
                    tryCatch(
                        {
                            GenomicAlignments::cigarNarrow(cigar = cigars[i], start = 1, end = ends[i])[1]
                        },
                        error = function(e) {
                            return("1=")
                        }
                    )
                }, FUN.VALUE = character(1))
                ## Get alignment length from regional cigars
                new.aln.len <- vapply(GenomicAlignments::explodeCigarOpLengths(cigar = new.cigar), FUN = sum, FUN.VALUE = numeric(1))
                ## Get matched bases from regional cigars
                new.n.match <- vapply(GenomicAlignments::explodeCigarOpLengths(cigar = new.cigar, ops = c("=", "M")), FUN = sum, FUN.VALUE = numeric(1))

                ## Update PAF coordinates and cigar string ##
                paf.aln$q.start[paf.aln$strand == "-"] <- GenomicRanges::start(query.end[GenomicRanges::strand(query.end) == "-"])
                paf.aln$q.end[paf.aln$strand == "+"] <- GenomicRanges::start(query.end[GenomicRanges::strand(query.end) == "+"])
                # paf.aln$t.start <- GenomicRanges::start(target.region.gr)
                paf.aln$t.end <- GenomicRanges::end(target.region.gr)
                paf.aln$n.match <- new.n.match
                paf.aln$aln.len <- new.aln.len
                paf.aln$cg <- new.cigar
                paf.partial.aln[cut.end.idx, ] <- paf.aln
            }
        } else {
            paf.partial.aln <- tibble::tibble()
        }
        ## Merge final PAF
        paf <- dplyr::bind_rows(paf.within.aln, paf.embedded.aln, paf.partial.aln)

        stopTimedMessage(ptm)
        return(paf)
    } else {
        warning("\nNone of the PAF ranges overlap user defined 'target.region', exiting ...")
        return(NULL)
    }
}


#' Subset PAF alignments at desired genomic range.
#'
#' This function takes loaded PAF alignments using \code{\link{readPaf}} function and then subsets
#' as well as cuts PAF alignments at desired target coordinates.
#'
#' @inheritParams breakPaf
#' @inheritParams subsetPafAlignments
#' @return A \code{tibble} of subsetted PAF alignments.
#' @importFrom GenomicRanges makeGRangesFromDataFrame strand
#' @importFrom IRanges findOverlaps
#' @importFrom GenomeInfoDb seqlevels seqnames
#' @importFrom methods as
#' @importFrom dplyr bind_rows
#' @author David Porubsky
#' @export
#' @examples
#' ## Get PAF to plot
#' paf.file <- system.file("extdata", "test_ava.paf", package = "SVbyEye")
#' ## Read in PAF
#' paf.table <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = "cg")
#' ## Subset multiple PAF alignments based on desired target region
#' target.region <- "HG01358_2:100000-200000"
#' subsetPaf(paf.table, target.region = target.region)
#'
subsetPaf <- function(paf.table, target.region = NULL) {
    ## Check user input ##
    ## Make sure PAF has at least 12 mandatory fields
    if (ncol(paf.table) >= 12) {
        paf <- paf.table
    } else {
        stop("\nSubmitted PAF alignments do not contain a minimum of 12 mandatory fields, see PAF file format definition !!!")
    }
    ## Check if PAF contains expected cg field containing CIGAR string
    if (!"cg" %in% colnames(paf)) {
        stop(paste0(
            "\nExpected PAF field 'cg' containing CIGAR string is missing, ",
            "therefore alignments cannot be narrowed down exactly to the user defined 'target.region'!!!"
        ))
    }

    ## Filter alignments by target region ##
    if (!is.null(target.region)) {
        if (is.character(target.region)) {
            target.region.gr <- methods::as(target.region, "GRanges")
        } else if (is(target.region, "GRanges")) {
            target.region.gr <- target.region
        } else {
            message("\nParameter 'target.region' can either be 'GRanges' object or character string 'chr#:start-end'!!!")
        }
        ## Make GRanges from all target alignment coordinates
        target.gr <- GenomicRanges::makeGRangesFromDataFrame(paf, seqnames.field = "t.name", start.field = "t.start", end.field = "t.end")
        ## Make sure all sequence levels are present
        seq.lvs <- unique(c(unique(paf$t.name), unique(paf$q.name)))
        GenomeInfoDb::seqlevels(target.gr) <- seq.lvs
        GenomeInfoDb::seqlevels(target.region.gr) <- seq.lvs
        ## Check if user defined target region overlaps any target sequence in submitted paf.table
        hits <- IRanges::findOverlaps(target.gr, target.region.gr)
        if (length(hits) > 0) {
            ## Get all possible target coordinates
            all.targets <- list()
            while (length(IRanges::findOverlaps(target.region.gr, target.gr)) > 0) {
                gr <- liftRangesToAlignment(paf.table = paf, gr = target.region.gr, direction = "target2query")
                all.targets[[length(all.targets) + 1]] <- target.region.gr[, 0]
                target.region.gr <- gr
            }
            all.targets[[length(all.targets) + 1]] <- target.region.gr[, 0]
            all.targets <- do.call(c, all.targets)
            all.targets <- all.targets[GenomeInfoDb::seqnames(all.targets) != "Failed"]
            ## IF not all sequence levels present in original PAF are reported in all.targets then
            ## make sure to include any remaining ranges transferable from query to target
            if (!all(seq.lvs %in% as.character(GenomeInfoDb::seqnames(all.targets)))) {
                gr <- liftRangesToAlignment(paf.table = paf, gr = all.targets, direction = "query2target")
                gr <- gr[GenomeInfoDb::seqnames(gr) != "Failed"]
                all.targets <- c(all.targets, gr[, 0])
                GenomicRanges::strand(all.targets) <- "*"
                all.targets <- unique(all.targets)
            }

            ## Remove self-alignments
            paf <- paf[paf$q.name != paf$t.name, ]

            ## Subset alignments for all target coordinates
            paf.l <- split(paf, paf$t.name)
            target.paf <- list()
            for (i in seq_along(paf.l)) {
                sub.paf <- paf.l[[i]]
                new.paf <- subsetPafAlignments(
                    paf.table = sub.paf,
                    target.region = all.targets[as.character(GenomeInfoDb::seqnames(all.targets)) == unique(sub.paf$t.name)]
                )
                target.paf[[i]] <- new.paf
            }
            subset.paf <- dplyr::bind_rows(target.paf)
            ## Reassign new coordinates
            subset.paf$q.len <- (subset.paf$q.end - subset.paf$q.start) + 1
            subset.paf$t.len <- (subset.paf$t.end - subset.paf$t.start) + 1
            subset.paf$q.start <- 1
            subset.paf$t.start <- 1
            subset.paf$q.end <- subset.paf$q.len
            subset.paf$t.end <- subset.paf$t.len
        }
    } else {
        message("\nNone of the PAF ranges overlap user defined 'target.region', exiting ...")
        return(NULL)
    }
    return(subset.paf)
}
