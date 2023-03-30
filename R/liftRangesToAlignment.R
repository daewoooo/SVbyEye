#' Function to lift coordinates to the alignment in PAF format.
#'
#' @param gr A \code{\link{GRanges-class}} object containing single or multiple ranges in query or target sequence coordinates.
#' @param direction One of the possible, lift ranges from query to target 'query2target' or vice versa 'target2query'.
#' @inheritParams breakPaf
#' @importFrom GenomicRanges GRanges findOverlaps mcols shift
#' @importFrom GenomicAlignments GAlignments mapFromAlignments
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom methods is
#' @return A \code{\link{GRanges-class}} object with resized original set of ranges.
#' @author David Porubsky
#' @export
#' @examples
#' ## Define range(s) to lift
#' roi.gr <- as("chr17:46645907-46697277", "GRanges")
#' ## Get PAF alignment to lift to
#' paf.file <- system.file("extdata", "test_lift1.paf", package = "SVbyEye")
#' ## Read in PAF
#' paf.table <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = "cg")
#' ## Lift target range to query coordinates
#' liftRangesToAlignment(paf.table = paf.table, gr = roi.gr, direction = "target2query")
#'
liftRangesToAlignment <- function(paf.table, gr = NULL, direction = "query2target") {
    ptm <- startTimedMessage(paste0("[liftRangesToAlignment] Lifting coordinates: ", direction))

    ## Check user input ##
    stopifnot(methods::is(gr, "GRanges"), methods::is(direction, "character"))
    ## Make sure PAF has at least 12 mandatory fields
    if (ncol(paf.table) >= 12 & 'cg' %in% colnames(paf.table)) {
      paf <- paf.table
    } else {
      stop("Submitted PAF alignments do not either contain required 12 mandatory fields or a CIGAR string reported in 'cg' column !!!")
    }

    ## Make sure defined ranges are present in the PAF alignment
    if (direction == "query2target") {
        paf.query.gr <- GenomicRanges::GRanges(seqnames = paf$q.name, ranges = IRanges::IRanges(start = paf$q.start, end = paf$q.end))
        ## Optional step to get only overlapping ranges [To be tested]
        #gr <- GenomicRanges::pintersect(paf.query.gr, gr)
        #gr <- gr[gr$hit == TRUE]
        hits <- suppressWarnings(GenomicRanges::findOverlaps(gr, paf.query.gr, ignore.strand=TRUE))
        from <- 'query'
    } else if (direction == "target2query") {
        paf.target.gr <- GenomicRanges::GRanges(seqnames = paf$t.name, ranges = IRanges::IRanges(start = paf$t.start, end = paf$t.end))
        hits <- suppressWarnings(GenomicRanges::findOverlaps(gr, paf.target.gr, ignore.strand=TRUE))
        from <- 'target'
    } else {
        stop("Parameter 'direction' can take only values: 'query2target' or 'target2query'")
    }

    if (length(hits) == 0) {
        message(paste0("None of the user defined ranges are present in the PAF ", from, " coordinates!!!"))
        ## Return zero positioned ranges if there is no overlap
        gr.lifted <- GenomicRanges::GRanges(seqnames = rep("Failed", length(gr)), ranges = IRanges::IRanges(start = 1, end = 0))
        return(gr.lifted)
    } else {
        ## Lift ranges to PAF alignment
        if (direction == "query2target") {
          ## Map from alignment ##
          gr.lift <- gr[S4Vectors::queryHits(hits)]
          GenomicRanges::mcols(gr.lift)$idx <- S4Vectors::queryHits(hits)
          names(gr.lift) <- paste0("aln", S4Vectors::subjectHits(hits))
          ## Define Genomic alignments object
          paf.aln <- paf[unique(S4Vectors::subjectHits(hits)),]
          #paf.aln <- paf
          alignments <- GenomicAlignments::GAlignments(
            seqnames = paf.aln$t.name,
            pos = as.integer(paf.aln$t.start) + 1L,
            #pos = pos,
            cigar = paf.aln$cg,
            strand = GenomicRanges::strand(paf.aln$strand),
            names = paste0("aln", unique(S4Vectors::subjectHits(hits)))
          )
          ## Flip desired coordinates in case of minus alignments and adjust to alignment bounds
          if (any(as.character(GenomicRanges::strand(alignments)) == "-")) {
            ## Get alignments to flip
            aln.idx <- which(as.character(GenomicRanges::strand(alignments)) == "-")
            mask <- which(names(gr.lift) %in% names(alignments)[aln.idx])
            ## Flip reverse alignments
            #bounds <- IRanges::IRanges(start = 0L, end = paf.aln$q.end[aln.idx])
            bounds <- IRanges::IRanges(start = 0L, end = paf.aln$q.end[match(names(gr.lift[mask]), names(alignments))])
            suppressWarnings( IRanges::ranges(gr.lift[mask]) <- IRanges::reflect(x = IRanges::ranges(gr.lift[mask]), bounds = bounds) )
          }
          ## Adjust start position for plus alignments alignment
          if (any(as.character(GenomicRanges::strand(alignments)) == "+")) {
            ## Get alignments to flip
            aln.idx <- which(as.character(GenomicRanges::strand(alignments)) == "+")
            mask <- which(names(gr.lift) %in% names(alignments)[aln.idx])
            ## Adjust start position of plus alignments to a local alignment coordinates
            suppressWarnings(
              gr.lift[mask] <- GenomicRanges::shift(gr.lift[mask], shift = paf.aln$q.start[match(names(gr.lift[mask]), names(alignments))] * -1)
            )
            ## Remove negative ranges
            remove <- start(gr.lift) < 0
            gr.lift <- gr.lift[!remove]
          }
          ## Lift ranges
          gr.lifted <- GenomicAlignments::mapFromAlignments(x = gr.lift, alignments = alignments)
          names(gr.lifted) <- NULL
          ## Add index corresponding to original input ranges
          gr.lifted$idx <- gr.lift$idx[gr.lifted$xHits]

        } else {
          ## Map to alignment ##
          gr.lift <- gr
          ## Define PAF alignment object
          paf.aln <- paf[unique(S4Vectors::subjectHits(hits)),]
          #paf.aln <- paf
          alignments <- GenomicAlignments::GAlignments(
            seqnames = paf.aln$t.name,
            pos = as.integer(paf.aln$t.start) + 1L,
            cigar = paf.aln$cg,
            strand = GenomicRanges::strand(paf.aln$strand),
            names = paf.aln$q.name
          )
          ## Lift ranges
          gr.lifted <- GenomicAlignments::mapToAlignments(x = gr.lift, alignments = alignments)
          names(gr.lifted) <- NULL

          ## Flip coordinates in case of reverse alignment
          if (any(as.character(GenomicRanges::strand(alignments)) == "-")) {
            # start.reverse <- (paf.aln$q.end[gr.lifted$alignmentsHits] - end(gr.lifted)) + 1
            # end.reverse <- (paf.aln$q.end[gr.lifted$alignmentsHits] - start(gr.lifted)) + 1
            # ranges(gr.lifted) <- IRanges::IRanges(start=start.reverse, end=end.reverse)
            ## Get alignments to flip
            aln.idx <- which(as.character(GenomicRanges::strand(alignments)) == "-")
            mask <- gr.lifted$alignmentsHits %in% aln.idx
            ## Flip reverse alignments
            #bounds <- IRanges::IRanges(start = 0L, end = paf.aln$q.end[aln.idx])
            bounds <- IRanges::IRanges(start = 0L, end = paf.aln$q.end[gr.lifted$alignmentsHits])
            IRanges::ranges(gr.lifted[mask]) <- IRanges::reflect(x = IRanges::ranges(gr.lifted[mask]), bounds = bounds[mask])
          }
          ## Adjust start position for plus alignments alignment
          if (any(as.character(GenomicRanges::strand(alignments)) == "+")) {
            ## Get alignments to flip
            aln.idx <- which(as.character(GenomicRanges::strand(alignments)) == "+")
            mask <- gr.lifted$alignmentsHits %in% aln.idx
            ## Adjust start position of plus alignments to a local alignment coordinates
            #GenomicRanges::end(gr.lifted[mask]) <- abs( GenomicRanges::end(gr.lifted[mask]) + paf.aln$q.start[gr.lifted$alignmentsHits[mask]] )
            #GenomicRanges::start(gr.lifted[mask]) <- abs( GenomicRanges::start(gr.lifted[mask]) + paf.aln$q.start[gr.lifted$alignmentsHits[mask]] )
            gr.lifted[mask] <- GenomicRanges::shift(gr.lifted[mask], shift = paf.aln$q.start[gr.lifted$alignmentsHits[mask]])
          }
          ## Add index corresponding to original input ranges
          gr.lifted$idx <- gr.lifted$xHits
        }

        if (length(gr.lifted) > 0) {
          ## Add alignment strand information
          GenomicRanges::strand(gr.lifted) <- GenomicAlignments::strand(alignments)[gr.lifted$alignmentsHits]
        }

        ## Add empty ranges to coordinates that failed to lift
        if (!all(seq_along(gr) %in% gr.lifted$idx)) {
            failed.idx <- (seq_along(gr))[-gr.lifted$idx]
            decoy.gr <- GenomicRanges::GRanges(
                seqnames = rep("Failed", length(failed.idx)),
                ranges = IRanges::IRanges(start = 1, end = 0),
                xHits = 1L,
                alignmentsHits = 1L,
                idx = failed.idx
            )
            gr.lifted <- suppressWarnings(c(gr.lifted, decoy.gr))
            gr.lifted <- gr.lifted[order(gr.lifted$idx)]
        }
        ## Keep extra metacolumns if present
        if (ncol(GenomicRanges::mcols(gr)) > 0) {
            GenomicRanges::mcols(gr.lifted) <- c(GenomicRanges::mcols(gr.lifted), GenomicRanges::mcols(gr))
        }

        stopTimedMessage(ptm)
        ## Return lifted coordinates
        return(gr.lifted)
    }
}
