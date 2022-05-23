#' Function to parse CIGAR string into a set interval ranges.
#'
#' @param cigar.str A character string containing alignment represented as a CIGAR string.
#' @param coordinate.space A used defined coordinate space given CIGAR should be parsed against, either 'reference' or 'query'. 
#' @importFrom GenomicAlignments explodeCigarOps explodeCigarOpLengths cigarRangesAlongReferenceSpace cigarRangesAlongQuerySpace
#' @importFrom dplyr recode
#' @importFrom GenomicRanges resize
#' @return A \code{\link{IRangesList}} object.
#' @author David Porubsky
#' @export
#'
parseCigarString <- function(cigar.str=NULL, coordinate.space = 'reference') {
  ## Parse CIGAR string ##
  cigar.id <- GenomicAlignments::explodeCigarOps(cigar.str)[[1]]
  cigar.len <- GenomicAlignments::explodeCigarOpLengths(cigar.str)[[1]]
  ## Get Reference space coords
  if (coordinate.space == 'reference') {
    cigar.ranges <- GenomicAlignments::cigarRangesAlongReferenceSpace(cigar.str)[[1]]
  } else if (coordinate.space == 'query') {
    cigar.ranges <- GenomicAlignments::cigarRangesAlongQuerySpace(cigar.str)[[1]]
  } else {
    stop("Parameter 'coordinate.space' can only take values 'reference' or 'query' !!!")
  }
  ## Add insertion sizes
  cigar.ranges[cigar.id == 'I'] <- GenomicRanges::resize(cigar.ranges[cigar.id == 'I'], width = cigar.len[cigar.id == 'I'], fix = 'center')
  ## Translate cigar symbols
  cigar.id.trans <- dplyr::recode(cigar.id, '=' = 'match', 'I' = 'insertion', 'D' = 'deletion', 'X' = 'mismatch')
  ## Export ranges
  irl <- split(cigar.ranges, cigar.id.trans)
  return(irl)
}


#' Function to load PAF alignments into a set of genomic ranges.
#'
#' @param paf.file A character string containing alignment represented as a CIGAR string.
#' @param min.insertion.size A minimum size (in basepairs) of an insertion to be retained.
#' @param min.deletion.size A minimum size (in basepairs) of a deletion to be retained.
#' @param collapse.mismatches Set to \code{TRUE} if mismatches should be collapsed in order expand matched regions.
#' @inheritParams readPaf
#' @inheritParams parseCigarString
#' @importFrom GenomicRanges GRanges shift reduce
#' @return A \code{\link{GRanges}} object.
#' @author David Porubsky
#' @export
#'
cigar2ranges <- function(paf.file=NULL, coordinate.space='reference', min.insertion.size=50, min.deletion.size=50, collapse.mismatches=TRUE) {
  ## Read in coordinates from minimap2 output in PAF format
  paf.data <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = 'cg')
  qname <- paste(unique(paf.data$q.name), collapse = ';')
  ## Order by target sequence
  paf.data <- paf.data[order(paf.data$t.start),]
  
  ## Process alignments ##
  matches <- list()
  mismatches <- list()
  insertions <- list()
  deletions <- list()
  for (i in 1:nrow(paf.data)) {
    paf.aln <- paf.data[i,]
    ## Parse CIGAR string ##
    cg.ranges <- parseCigarString(cigar.str = paf.aln$cg, coordinate.space = coordinate.space)
    ## Get cigar ranges and offset ranges based on alignment starting position
    if (length(cg.ranges$match) > 0) {
      match.gr <- GenomicRanges::GRanges(seqnames = paf.aln$q.name, ranges = cg.ranges$match)
      match.gr <- GenomicRanges::shift(match.gr, shift = paf.aln$t.start)
    } else {
      match.gr <- NULL
    }  
    if (length(cg.ranges$mismatch) > 0) {
      mismatch.gr <- GenomicRanges::GRanges(seqnames = paf.aln$q.name, ranges = cg.ranges$mismatch)
      mismatch.gr <- GenomicRanges::shift(mismatch.gr, shift = paf.aln$t.start)
    } else {
      mismatch.gr <- NULL
    }  
    if (length(cg.ranges$deletion) > 0) { 
      del.gr <- GenomicRanges::GRanges(seqnames = paf.aln$q.name, ranges = cg.ranges$deletion)
      del.gr <- GenomicRanges::shift(del.gr, shift = paf.aln$t.start)
      ## Filter deletions by size
      del2reduce <- del.gr[width(del.gr) < min.deletion.size]
      del.gr <- GenomicRanges::reduce(del.gr[width(del.gr) >= min.deletion.size])
    } else {
      del2reduce <- NULL
      del.gr <- NULL
    }  
    if (length(cg.ranges$insertion) > 0) {  
      ins.gr <- GRanges(seqnames = paf.aln$q.name, ranges = cg.ranges$insertion)
      ins.gr <- GenomicRanges::shift(ins.gr, shift = paf.aln$t.start)
      ## Filter insertions by size
      ins.gr <- GenomicRanges::reduce(ins.gr[width(ins.gr) >= min.insertion.size])
    } else {
      ins.gr <- NULL
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
    
    ## Prepare data for export
    if (length(match.gr) > 0) {
      strand(match.gr) <- paf.aln$strand
      match.gr$aln.id <- paste0('aln', i)
      match.gr$cg <- 'M'
    } else {
      match.gr <- NULL
    } 
    if (length(mismatch.gr) > 0) {
      strand(mismatch.gr) <- paf.aln$strand
      mismatch.gr$aln.id <- paste0('aln', i)
      mismatch.gr$cg <- 'X'
    } else {
      mismatch.gr <- NULL
    }
    if (length(ins.gr) > 0) {
      strand(ins.gr) <- paf.aln$strand
      ins.gr$aln.id <- paste0('aln', i)
      ins.gr$cg <- 'I'
    } else {
      ins.gr <- NULL
    }
    if (length(del.gr) > 0) {
      strand(del.gr) <- paf.aln$strand
      del.gr$aln.id <- paste0('aln', i)
      del.gr$cg <- 'D'
    } else {
      del.gr <- NULL
    }
    
    matches[[i]] <- match.gr
    mismatches[[i]] <- mismatch.gr
    insertions[[i]] <- ins.gr
    deletions[[i]] <- del.gr
  }
  matches.gr <- do.call(c, matches)
  mismatches.gr <- do.call(c, mismatches)
  insertions.gr <- do.call(c, insertions)
  deletions.gr <- do.call(c, deletions)
  
  ## Add level to each separate alignment
  n.aln <- length(unique(matches.gr$aln.id))
  matches.gr$level <- rep(rep(c(1:2), times=n.aln)[1:n.aln], times=rle(matches.gr$aln.id)$lengths)
  insertions.gr$level <- matches.gr$level[match(insertions.gr$aln.id, matches.gr$aln.id)]
  deletions.gr$level <- matches.gr$level[match(deletions.gr$aln.id, matches.gr$aln.id)]
  mismatches.gr$level <- matches.gr$level[match(mismatches.gr$aln.id, matches.gr$aln.id)]
  
  ## Return 
  final.gr <- c(matches.gr, mismatches.gr, insertions.gr, deletions.gr)
  return(final.gr)
}

cigar2divergence <- function(cigar.str=NULL, coordinate.space = 'reference') {
  ## Parse CIGAR string ##
  cg.ranges <- parseCigarString(cigar.str = cigar.str, coordinate.space = coordinate.space)
}
