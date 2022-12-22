#' Function to break PAF alignment into matching bases between query and target sequence.
#' In addition, locations of inserted bases in query and target sequence can be reported as well.
#'
#' @param binsize A size of a bin in base pairs to split a PAF alignment into.
#' @inheritParams breakPafAlignment
#' @importFrom GenomicRanges GRanges shift width start end strand
#' @importFrom GenomicAlignments GAlignments mapToAlignments qwidth cigarNarrow explodeCigarOpLengths
#' @importFrom GenomeInfoDb seqlengths seqlevels seqnames
#' @importFrom dplyr tibble
#' @importFrom S4Vectors sapply
#' @return  A \code{tibble} object storing binned PAF alignments.
#' @author David Porubsky
#' @export
#' @examples
#'## Get PAF to bin ##
#'paf.file <- system.file("extdata", "test3.paf", package="SVbyEye")
#'## Read in PAF
#'paf.table <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = 'cg')
#'## Split a single PAF alignment into user defined bins
#'pafToBins(paf.table = paf.table[1,], binsize = 100)
#'
pafAlignmentToBins <- function(paf.aln=NULL, binsize=10000) {
  ## Get target and query ranges
  target.gr <- GenomicRanges::GRanges(seqnames = paf.aln$t.name, ranges=IRanges::IRanges(start = paf.aln$t.start, end = paf.aln$t.end))
  query.gr <- GenomicRanges::GRanges(seqnames = paf.aln$q.name, ranges=IRanges::IRanges(start = paf.aln$q.start, end = paf.aln$q.end))
  ## Create PAF alignment object
  alignment <- GenomicAlignments::GAlignments(seqnames = paf.aln$t.name, pos = 1L, cigar = paf.aln$cg, strand=GenomicRanges::strand(paf.aln$strand), names = 'target')
  ## Check if the alignments size is at least twice the size of desired bin size
  if ((GenomicAlignments::qwidth(alignment) / binsize) >= 2) {
    ## Get bins in target coordinates
    #target.len <- GenomicRanges::width(target.gr)
    #names(target.len) <- paf.aln$t.name
    aln.len <- GenomicAlignments::cigarWidthAlongReferenceSpace(cigar = paf.aln$cg)
    names(aln.len) <- paf.aln$t.name
    bins.gr <- GenomicRanges::tileGenome(seqlengths = aln.len, tilewidth = binsize, cut.last.tile.in.chrom = TRUE)
    ## Get bins in query coordinates
    query.bins.gr <- GenomicAlignments::mapToAlignments(bins.gr, alignments = alignment)
    GenomeInfoDb::seqlengths(query.bins.gr) <- GenomicAlignments::qwidth(alignment)
    GenomeInfoDb::seqlevels(query.bins.gr) <- paf.aln$q.name

    ## Convert to target coordinates
    target.bins.gr <- suppressWarnings( GenomicRanges::shift(bins.gr, shift = paf.aln$t.start - 1) )
    ## Convert to query coordinates
    if (paf.aln$strand == '-') {
      query.bins.gr <- mirrorRanges(gr = query.bins.gr)
      query.bins.gr <- suppressWarnings( GenomicRanges::shift(query.bins.gr, shift = paf.aln$q.start) )
    } else {
      query.bins.gr <- suppressWarnings( GenomicRanges::shift(query.bins.gr, shift = paf.aln$q.start) )
    }

    ## Subset the cigar string per target region
    starts <- as.numeric(GenomicRanges::start(bins.gr))
    ends <- as.numeric(GenomicRanges::end(bins.gr))
    #cigars.region <- S4Vectors::sapply(1:length(starts), function(i) GenomicAlignments::cigarNarrow(cigar = paf.aln$cg, start = starts[i], end = ends[i])[1])
    cigars.region <- S4Vectors::sapply(1:length(starts), function(i) tryCatch(
      {GenomicAlignments::cigarNarrow(cigar = paf.aln$cg, start = starts[i], end = ends[i])[1]}
      , error = function(e) {return('1=')}
        )
      )
    ## Get alignment length from regional cigars
    aln.len <- S4Vectors::sapply(cigars.region, function(cg) sum(GenomicAlignments::explodeCigarOpLengths(cigar = cg)[[1]]), USE.NAMES = FALSE)
    ## Get matched bases from regional cigars
    n.match <- S4Vectors::sapply(cigars.region, function(cg) sum(GenomicAlignments::explodeCigarOpLengths(cigar = cg, ops = c('=','M'))[[1]]), USE.NAMES = FALSE)
    ## Create binned paf alignment
    binned.paf.aln <- dplyr::tibble(q.name=as.character(GenomeInfoDb::seqnames(query.bins.gr)),
                                    q.len=paf.aln$q.len,
                                    q.start=as.numeric(GenomicRanges::start(query.bins.gr)),
                                    q.end=as.numeric(GenomicRanges::end(query.bins.gr)),
                                    strand=paf.aln$strand,
                                    t.name=as.character(GenomeInfoDb::seqnames(target.bins.gr)),
                                    t.len=paf.aln$t.len,
                                    t.start=as.numeric(GenomicRanges::start(target.bins.gr)),
                                    t.end=as.numeric(GenomicRanges::end(target.bins.gr)),
                                    n.match=n.match,
                                    aln.len=aln.len,
                                    mapq=paf.aln$mapq,
                                    cg=cigars.region)
    ## Remove alignments set to size 1 due to 'CIGAR is empty after narrowing' error message
    binned.paf.aln <- binned.paf.aln[binned.paf.aln$aln.len > 1,]
    ## Return binned paf alignments
    return(binned.paf.aln)
  } else {
    return(paf.aln)
  }
}


#' A wrapper function for \code{\link{pafAlignmentToBins}} binning multiple PAF alignments
#' into a set of binned alignments between query and target sequence.
#'
#' @inheritParams breakPaf
#' @inheritParams pafAlignmentToBins
#' @importFrom dplyr bind_rows
#' @return A \code{list} of \code{tibble} objects storing matched ('M') alignments as well as structurally variable ('SV') bases if 'report.sv' is TRUE.
#' @author David Porubsky
#' @export
#' @examples
#'## Get PAF to bin ##
#'paf.file <- system.file("extdata", "test1.paf", package="SVbyEye")
#'## Read in PAF
#'paf.table <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = 'cg')
#'## Split multiple PAF alignments into user defined bins
#'pafToBins(paf.table = paf.table, binsize = 10000)
#'
pafToBins <- function(paf.table=NULL, binsize=10000) {

  ptm <- startTimedMessage(paste0("[pafToBins] Binning PAF alignments, binsize=", binsize, "bp"))
  ## Split each PAF record into user defined bins
  binned <- list()
  for (i in 1:nrow(paf.table)) {
    paf.aln <- paf.table[i,]
    paf.aln.binned <- pafAlignmentToBins(paf.aln = paf.aln, binsize = binsize)
    paf.aln.binned$aln.id <- i
    binned[[i]] <- paf.aln.binned
  }

  stopTimedMessage(ptm)
  ## Return binned paf alignments
  return(dplyr::bind_rows(binned))
}

