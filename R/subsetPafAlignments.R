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
        stop("\nSubmitted PAF alignments do not contain a minimum of 12 mandatory fields, see PAF file format definition !!!")
    }
    ## Check if PAF contains expected cg field containing CIGAR string
    if (!'cg' %in% colnames(paf)) {
      stop(paste0("\nExpected PAF field 'cg' containing CIGAR string is missing, ",
                  "therefore alignments cannot be narrowed down exactly to the user defined 'target.region'!!!"))
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
        hits <- GenomicRanges::findOverlaps(target.gr, target.region.gr) ## TODO Take only overlaps, regions itself are not shrinked!!!
        if (length(hits) > 0) {
            paf <- paf[S4Vectors::queryHits(hits), ]
            target.gr <- target.gr[S4Vectors::queryHits(hits),]
        } else {
            message("\nNone of the PAF ranges overlap user defined 'target.region', exiting ...")
            return(NULL)
            break
        }
    }

    ## Get alignments with target region being embedded within
    hits <- IRanges::findOverlaps(target.region.gr, target.gr, type = "within")
    within.idx <- S4Vectors::subjectHits(hits)
    ## Get alignments with target region being partially overlapping
    partial.idx <- setdiff(seq_along(target.gr), within.idx)

    ## Process alignments with target region being embedded within
    if (length(within.idx) > 0) {
      ## Create PAF alignment object
      paf.aln <- paf[within.idx,]
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
      query.coords <- suppressMessages( liftRangesToAlignment(paf.table = paf.aln, gr = target.region.gr, direction = 'target2query') )
      ## Set target range to local alignment coordinates
      suppressWarnings( target.local.gr <- GenomicRanges::shift(target.region.gr[S4Vectors::queryHits(hits)],
                                                                shift = -(start(alignments)[S4Vectors::subjectHits(hits)] - 1)) )
      start(target.local.gr[start(target.local.gr) == 0]) <- 1

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
      paf.within.aln <- paf.aln
    } else {
      paf.within.aln <- tibble::tibble()
    }

    ## Process alignments with target region being partially overlapping
    if (length(partial.idx) > 0) {
      paf.partial.aln <- paf[partial.idx,]
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
        query.start <- suppressMessages( liftRangesToAlignment(paf.table = paf.aln, gr = target.start, direction = 'target2query') )
        ## Set target range to local alignment coordinates
        suppressWarnings( target.local.gr <- GenomicRanges::shift(target.start[S4Vectors::queryHits(hits)],
                                                                  shift = -(start(alignments)[S4Vectors::subjectHits(hits)] - 1)) )
        start(target.local.gr[start(target.local.gr) == 0]) <- 1

        ## Cut cigar string
        starts <- GenomicRanges::start(target.local.gr)
        cigars <- alignments@cigar[S4Vectors::subjectHits(hits)]
        widths <- width(alignments[S4Vectors::subjectHits(hits)])
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
        new.aln.len <- sapply(GenomicAlignments::explodeCigarOpLengths(cigar = new.cigar), sum)
        ## Get matched bases from regional cigars
        new.n.match <- sapply(GenomicAlignments::explodeCigarOpLengths(cigar = new.cigar, ops = c("=", "M")), sum)

        ## Update PAF coordinates and cigar string ##
        paf.aln$q.end[paf.aln$strand == '-'] <- start(query.start[strand(query.start) == '-'])
        paf.aln$q.start[paf.aln$strand == '+'] <- start(query.start[strand(query.start) == '+'])
        paf.aln$t.start <- GenomicRanges::start(target.region.gr)
        paf.aln$t.end <- GenomicRanges::end(target.region.gr)
        paf.aln$n.match <- new.n.match
        paf.aln$aln.len <- new.aln.len
        paf.aln$cg <- new.cigar
        paf.partial.aln[cut.start.idx,] <- paf.aln
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
        query.end <- suppressMessages( liftRangesToAlignment(paf.table = paf.aln, gr = target.end, direction = 'target2query') )
        ## Set target range to local alignment coordinates
        suppressWarnings( target.local.gr <- GenomicRanges::shift(target.end[S4Vectors::queryHits(hits)],
                                                                  shift = -(start(alignments)[S4Vectors::subjectHits(hits)] - 1)) )
        start(target.local.gr[start(target.local.gr) == 0]) <- 1

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
        new.aln.len <- sapply(GenomicAlignments::explodeCigarOpLengths(cigar = new.cigar), sum)
        ## Get matched bases from regional cigars
        new.n.match <- sapply(GenomicAlignments::explodeCigarOpLengths(cigar = new.cigar, ops = c("=", "M")), sum)

        ## Update PAF coordinates and cigar string ##
        paf.aln$q.start[paf.aln$strand == '-'] <- start(query.end[strand(query.end) == '-'])
        paf.aln$q.end[paf.aln$strand == '+'] <- start(query.end[strand(query.end) == '+'])
        paf.aln$t.start <- GenomicRanges::start(target.region.gr)
        paf.aln$t.end <- GenomicRanges::end(target.region.gr)
        paf.aln$n.match <- new.n.match
        paf.aln$aln.len <- new.aln.len
        paf.aln$cg <- new.cigar
        paf.partial.aln[cut.end.idx,] <- paf.aln
      }
    } else {
      paf.partial.aln <- tibble::tibble()
    }
    ## Merge final PAF
    paf <- dplyr::bind_rows(paf.within.aln, paf.partial.aln)

    stopTimedMessage(ptm)
    return(paf)
}
