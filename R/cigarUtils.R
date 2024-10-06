#' Function to parse CIGAR string into a set interval ranges.
#'
#' @param cigar.str A character string containing alignment represented as a CIGAR string.
#' @param coordinate.space A used defined coordinate space given CIGAR should be parsed against, either 'reference' or 'query'.
#' @importFrom GenomicAlignments explodeCigarOps explodeCigarOpLengths cigarRangesAlongReferenceSpace cigarRangesAlongQuerySpace
#' @importFrom dplyr recode
#' @return A \code{\link{IRangesList}} object.
#' @author David Porubsky
#' @export
#' @examples
#' ## Get PAF to process
#' paf.file <- system.file("extdata", "test3.paf", package = "SVbyEye")
#' ## Read in PAF
#' paf.table <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = "cg")
#' ## Parse CIGAR in reference (target) space
#' parseCigarString(cigar.str = paf.table$cg)
#' ## Parse CIGAR in query space
#' parseCigarString(cigar.str = paf.table$cg, coordinate.space = "query")
#'
parseCigarString <- function(cigar.str = NULL, coordinate.space = "reference") {
    ## Parse CIGAR string ##
    cigar.id <- GenomicAlignments::explodeCigarOps(cigar.str)[[1]]
    cigar.len <- GenomicAlignments::explodeCigarOpLengths(cigar.str)[[1]]
    ## Get Reference space coords
    if (coordinate.space == "reference") {
        cigar.ranges <- GenomicAlignments::cigarRangesAlongReferenceSpace(cigar.str)[[1]]
        # start(cigar.ranges[cigar.id == 'I']) <- end(cigar.ranges[cigar.id == 'I'])
    } else if (coordinate.space == "query") {
        cigar.ranges <- GenomicAlignments::cigarRangesAlongQuerySpace(cigar.str, after.soft.clipping = TRUE)[[1]]
    } else {
        stop("Parameter 'coordinate.space' can only take values 'reference' or 'query' !!!")
    }
    ## Translate cigar symbols
    cigar.id.trans <- dplyr::recode(cigar.id, "=" = "match", "M" = "match", "I" = "insertion", "D" = "deletion", "X" = "mismatch")
    ## Export ranges
    irl <- split(cigar.ranges, cigar.id.trans)
    return(irl)
}


#' Function to load CIGAR string reported in PAF alignments into a set of genomic ranges.
#'
#' @param min.insertion.size A minimum size (in base pairs) of an insertion to be retained.
#' @param min.deletion.size A minimum size (in base pairs) of a deletion to be retained.
#' @param collapse.mismatches Set to \code{TRUE} if mismatches should be collapsed in order expand matched regions.
#' @inheritParams breakPaf
#' @inheritParams parseCigarString
#' @importFrom GenomicRanges GRanges width shift reduce strand
#' @importFrom IRanges ranges reflect
#' @return A \code{\link{GRanges}} object.
#' @author David Porubsky
#' @export
#' @examples
#' ## Get PAF to process
#' paf.file <- system.file("extdata", "test3.paf", package = "SVbyEye")
#' ## Parse CIGAR into a set of genomic ranges
#' paf.table <- readPaf(paf.file = paf.file)
#' cigar2ranges(paf.table = paf.table)
#'
cigar2ranges <- function(paf.table = NULL, coordinate.space = "reference", min.insertion.size = 50, min.deletion.size = 50, collapse.mismatches = TRUE) {
    ## Check user input ##
    ## Make sure PAF has at least 12 mandatory fields
    if (ncol(paf.table) >= 12) {
        paf <- paf.table
    } else {
        stop("Submitted PAF alignments do not contain a minimum of 12 mandatory fields, see PAF file format definition !!!")
    }

    ## Process alignments ##
    matches <- list()
    mismatches <- list()
    insertions <- list()
    deletions <- list()
    for (i in seq_len(nrow(paf.table))) {
        paf.aln <- paf.table[i, ]
        ## Parse CIGAR string ##
        cg.ranges <- parseCigarString(cigar.str = paf.aln$cg, coordinate.space = coordinate.space)
        ## Get sequence ID
        if (coordinate.space == "reference") {
            seqname <- paf.aln$t.name
        } else if (coordinate.space == "query") {
            seqname <- paf.aln$q.name
        } else {
            stop("Parameter 'coordinate.space' can only be either 'reference' or 'query'!!!")
        }
        ## Get cigar ranges and offset ranges based on alignment starting position
        if (length(cg.ranges$match) > 0) {
            match.gr <- GenomicRanges::GRanges(seqnames = seqname, ranges = cg.ranges$match)
        } else {
            match.gr <- GenomicRanges::GRanges()
        }
        if (length(cg.ranges$mismatch) > 0) {
            mismatch.gr <- GenomicRanges::GRanges(seqnames = seqname, ranges = cg.ranges$mismatch)
            mismatch.gr$size <- 1
        } else {
            mismatch.gr <- GenomicRanges::GRanges()
        }
        if (length(cg.ranges$deletion) > 0) {
            del.gr <- GenomicRanges::GRanges(seqnames = seqname, ranges = cg.ranges$deletion)
            ## In case of query coordinates get deletion size from the reference coordinates
            if (coordinate.space == "query") {
                ref.ranges <- parseCigarString(cigar.str = paf.aln$cg, coordinate.space = "reference")
                del.gr$size <- GenomicRanges::width(ref.ranges$deletion)
            } else {
                del.gr$size <- GenomicRanges::width(del.gr)
            }
            ## Filter deletions by size [in reference coordinates]
            del2reduce <- del.gr[del.gr$size < min.deletion.size]
            del.gr <- del.gr[del.gr$size >= min.deletion.size]
        } else {
            del2reduce <- GenomicRanges::GRanges()
            del.gr <- GenomicRanges::GRanges()
        }
        if (length(cg.ranges$insertion) > 0) {
            ins.gr <- GRanges(seqnames = seqname, ranges = cg.ranges$insertion)
            ## In case of reference coordinates get insertion size from the query coordinates
            if (coordinate.space == "reference") {
                qry.ranges <- parseCigarString(cigar.str = paf.aln$cg, coordinate.space = "query")
                ins.gr$size <- GenomicRanges::width(qry.ranges$insertion)
            } else {
                ins.gr$size <- GenomicRanges::width(ins.gr)
            }
            ## Filter insertions by size [in query coordinates]
            ins2reduce <- ins.gr[ins.gr$size < min.insertion.size]
            ins.gr <- ins.gr[ins.gr$size >= min.insertion.size]
        } else {
            ins2reduce <- GenomicRanges::GRanges()
            ins.gr <- GenomicRanges::GRanges()
        }

        ## Collapse simple mismatches
        if (collapse.mismatches) {
            match.gr <- GenomicRanges::reduce(c(match.gr, mismatch.gr))
        } else {
            match.gr <- GenomicRanges::reduce(match.gr)
        }
        ## Collapse filtered deletions
        if (!is.null(del2reduce) & length(del2reduce) > 0) {
            match.gr <- GenomicRanges::reduce(c(match.gr, del2reduce))
        }
        ## Collapse filtered insertions
        if (!is.null(ins2reduce) & length(ins2reduce) > 0) {
            match.gr <- GenomicRanges::reduce(c(match.gr, ins2reduce))
        }
        ## Set size of the final matched bases
        match.gr$size <- GenomicRanges::width(match.gr)

        ## Convert to query or target coordinates
        if (coordinate.space == "reference") {
            match.gr <- GenomicRanges::shift(match.gr, shift = paf.aln$t.start)
            mismatch.gr <- GenomicRanges::shift(mismatch.gr, shift = paf.aln$t.start)
            del.gr <- GenomicRanges::shift(del.gr, shift = paf.aln$t.start)
            ins.gr <- GenomicRanges::shift(ins.gr, shift = paf.aln$t.start)
        } else {
            ## Flip query coordinates in case of reverse alignment
            if (paf.aln$strand == '-') {
                bounds <- IRanges::IRanges(start = 0L, end = paf.aln$q.end)
                IRanges::ranges(match.gr) <- IRanges::reflect(x = IRanges::ranges(match.gr), bounds = bounds)
                IRanges::ranges(mismatch.gr) <- IRanges::reflect(x = IRanges::ranges(mismatch.gr), bounds = bounds)
                IRanges::ranges(del.gr) <- IRanges::reflect(x = IRanges::ranges(del.gr), bounds = bounds)
                IRanges::ranges(ins.gr) <- IRanges::reflect(x = IRanges::ranges(ins.gr), bounds = bounds)
            } else {
               match.gr <- GenomicRanges::shift(match.gr, shift = paf.aln$q.start)
               mismatch.gr <- GenomicRanges::shift(mismatch.gr, shift = paf.aln$q.start)
               del.gr <- GenomicRanges::shift(del.gr, shift = paf.aln$q.start)
               ins.gr <- GenomicRanges::shift(ins.gr, shift = paf.aln$q.start)
            }
        }

        ## Prepare data for export
        if (length(match.gr) > 0) {
            GenomicRanges::strand(match.gr) <- paf.aln$strand
            match.gr$aln.id <- paste0("aln", i)
            match.gr$cg <- "M"
        }
        if (length(mismatch.gr) > 0) {
            GenomicRanges::strand(mismatch.gr) <- paf.aln$strand
            mismatch.gr$aln.id <- paste0("aln", i)
            mismatch.gr$cg <- "X"
        }
        if (length(ins.gr) > 0) {
            GenomicRanges::strand(ins.gr) <- paf.aln$strand
            ins.gr$aln.id <- paste0("aln", i)
            ins.gr$cg <- "I"
        }
        if (length(del.gr) > 0) {
            GenomicRanges::strand(del.gr) <- paf.aln$strand
            del.gr$aln.id <- paste0("aln", i)
            del.gr$cg <- "D"
        }

        matches[[length(matches) + 1]] <- match.gr
        mismatches[[length(mismatches) + 1]] <- mismatch.gr
        insertions[[length(insertions) + 1]] <- ins.gr
        deletions[[length(deletions) + 1]] <- del.gr
    }
    matches.gr <- do.call(c, matches)
    mismatches.gr <- do.call(c, mismatches)
    insertions.gr <- do.call(c, insertions)
    deletions.gr <- do.call(c, deletions)

    ## Add level to each separate alignment
    n.aln <- length(unique(matches.gr$aln.id))
    matches.gr$level <- rep(rep(seq_len(2), times = n.aln)[seq_len(n.aln)], times = rle(matches.gr$aln.id)$lengths)
    insertions.gr$level <- matches.gr$level[match(insertions.gr$aln.id, matches.gr$aln.id)]
    deletions.gr$level <- matches.gr$level[match(deletions.gr$aln.id, matches.gr$aln.id)]
    mismatches.gr$level <- matches.gr$level[match(mismatches.gr$aln.id, matches.gr$aln.id)]

    ## Return
    final.gr <- c(matches.gr, mismatches.gr, insertions.gr, deletions.gr)
    return(final.gr)
}

# cigar2divergence <- function(cigar.str=NULL, coordinate.space = 'reference') {
#   ## Parse CIGAR string ##
#   cg.ranges <- parseCigarString(cigar.str = cigar.str, coordinate.space = coordinate.space)
# }
