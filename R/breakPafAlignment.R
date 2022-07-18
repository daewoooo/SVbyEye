#' Function to break PAF alignment into matching bases between query and target sequence.
#' In addition, locations of inserted bases in query and target sequence can be reported as well.
#'
#' @param paf.aln A \code{data.frame} or \code{tibble} containing a single PAF record with 12 mandatory columns.
#' @param report.sv Set to \code{TRUE} if to report also ranges of deleted and inserted bases.
#' @inheritParams cigar2ranges
#' @importFrom GenomicRanges GRanges shift reduce
#' @importFrom GenomicAlignments GAlignments mapToAlignments qwidth cigarNarrow explodeCigarOpLengths
#' @importFrom dplyr tibble bind_rows
#' @importFrom S4Vectors sapply
#' @return  A \code{list} of \code{tibble} objects storing matched ('M') alignments as well as structurally vairable ('SV') bases if 'report.sv' is TRUE.
#' @author David Porubsky
#' @export
#'
breakPafAlignment <- function(paf.aln=NULL, min.deletion.size=50, min.insertion.size=50, collapse.mismatches=TRUE, report.sv=TRUE) {
  ## Helper function ##
  ## This function flips set of ranges in a mirrored fashion
  mirrorRanges <- function(gr, seqlength=NULL) {
    if (!is.null(seqlength)) {
      gr.len <- seqlength
    } else if (!is.na(seqlengths(gr))) {
      gr.len <- seqlengths(gr)
    } else {
      stop('No seglength provided!!!')
    }
    if (!all(end(gr) <= gr.len)) {
      stop("One or all submitted ranges are outside of defined seqlength!!!")
    }
    starts <- gr.len - end(gr)
    ends <- (gr.len - start(gr)) 
    #starts <- gr.len - cumsum(width(gr))
    #ends <- starts + width(gr)
    new.gr <- GenomicRanges::GRanges(seqnames = seqnames(gr), ranges = IRanges(start=starts, end=ends), strand = strand(gr))
    suppressWarnings( seqlengths(new.gr) <- gr.len )
    return(new.gr)
  }
  
  ## Parse CIGAR string ##
  t.ranges <- parseCigarString(cigar.str = paf.aln$cg, coordinate.space = 'reference')
  q.ranges <- parseCigarString(cigar.str = paf.aln$cg, coordinate.space = 'query')
  ## Create paf alignment object
  alignment <- GenomicAlignments::GAlignments(seqnames = paf.aln$t.name, pos = 1L, cigar = paf.aln$cg, strand=GenomicRanges::strand(paf.aln$strand), names = 'target')
  #alignment <- GenomicAlignments::GAlignments(seqnames = paf.aln$t.name, pos = 1L, cigar = paf.aln$cg, strand=strand(paf.aln$strand), names = 'query')
  #test.gr <- GRanges(seqnames = 'target', ranges=IRanges(start=95514, width = width(match.gr)))
  ## Get cigar ranges and offset ranges based on alignment starting position
  if (length(t.ranges$match) > 0) {
    match.gr <- GenomicRanges::GRanges(seqnames = paf.aln$t.name, ranges = t.ranges$match)
  } else {
    match.gr <- GenomicRanges::GRanges()
  }  
  if (length(t.ranges$mismatch) > 0) {
    mismatch.gr <- GenomicRanges::GRanges(seqnames = paf.aln$t.name, ranges = t.ranges$mismatch)
  } else {
    mismatch.gr <- GenomicRanges::GRanges()
  }  
  if (length(t.ranges$deletion) > 0) {
    ## Get deletion position in reference space
    target.del.gr <- GenomicRanges::GRanges(seqnames = paf.aln$t.name, ranges = t.ranges$deletion)
    ## Get deletion ranges in query space
    query.del.gr <- GenomicRanges::GRanges(seqnames = paf.aln$q.name, ranges = q.ranges$del)
    ## Filter deletion by size
    del.mask <- width(target.del.gr) >= min.deletion.size
    del2reduce <- target.del.gr[!del.mask]
    target.del.gr <- target.del.gr[del.mask]
    query.del.gr <- query.del.gr[del.mask]
  } else {
    del2reduce <- GenomicRanges::GRanges()
    target.del.gr <- GenomicRanges::GRanges()
    query.del.gr <- GenomicRanges::GRanges()
  } 
    #del.gr <- GenomicRanges::GRanges(seqnames = paf.aln$t.name, ranges = cg.ranges$deletion)
    ### Filter deletions by size
    #del2reduce <- del.gr[GenomicRanges::width(del.gr) < min.deletion.size]
    #del.gr <- GenomicRanges::reduce(del.gr[GenomicRanges::width(del.gr) >= min.deletion.size])
  #} else {
  #  del2reduce <- NULL
  #  del.gr <- NULL
  #}  
  if (length(t.ranges$insertion) > 0) { 
    ## Get insertion position in reference space
    target.ins.gr <- GenomicRanges::GRanges(seqnames = paf.aln$t.name, ranges = t.ranges$insertion)
    ## Get insertion ranges in query space
    query.ins.gr <- GenomicRanges::GRanges(seqnames = paf.aln$q.name, ranges = q.ranges$insertion)
    ## Filter insertions by size
    ins.mask <- width(query.ins.gr) >= min.insertion.size
    target.ins.gr <- target.ins.gr[ins.mask]
    query.ins.gr <- query.ins.gr[ins.mask]
  } else {
    target.ins.gr <- GenomicRanges::GRanges()
    query.ins.gr <- GenomicRanges::GRanges()
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
  
  ## Break alignment [matched bases]
  query.match.gr <- GenomicAlignments::mapToAlignments(match.gr, alignments = alignment)
  seqlengths(query.match.gr) <- GenomicAlignments::qwidth(alignment)
  seqlevels(query.match.gr) <- paf.aln$q.name
  ## Break alignment [insertion bases]
  if (length(query.ins.gr) > 0) {
      ## Break target at insertion sites [should be a 0-bp break in target coords]
      match.gr <- disjoin(c(match.gr, target.ins.gr))
      ## Break query at insertion sites [should be a range in query coords]
      query.match.gr <- disjoin(c(query.match.gr, query.ins.gr))
      ## Get inserted range(s) to be removed from query
      hits <- findOverlaps(query.match.gr, query.ins.gr)
      remove.ins <- queryHits(hits)
  }
  
  ## Convert to target coordinates
  target.gr <- GenomicRanges::shift(match.gr, shift = paf.aln$t.start)
  
  ## Convert to query coordinates
  if (paf.aln$strand == '-') {
    query.match.gr <- mirrorRanges(gr = query.match.gr)
    query.gr <- suppressWarnings( GenomicRanges::shift(query.match.gr, shift = paf.aln$q.start) )
  } else {
    query.gr <- suppressWarnings( GenomicRanges::shift(query.match.gr, shift = paf.aln$q.start) )
  }  

  # starts <- max(end(query.match.gr)) - cumsum(width(query.match.gr))
  # ends <- starts + width(query.match.gr)
  #test.gr <- GRanges(seqnames = 'test', ranges = IRanges(start=7011, end=7012))
  # GenomicRanges::shift(test.gr, shift = paf.aln$q.start)
  
  # query.aln.widths <- width(query.match.gr)
  # if (paf.aln$strand == '-') {
  #   starts <- cumsum(c(paf.aln$q.end, -query.aln.widths))[-1]
  #   ends <- cumsum(c(paf.aln$q.end, -query.aln.widths))[-(length(query.aln.widths) + 1)]
  #   query.gr <- GenomicRanges::GRanges(seqnames = paf.aln$q.name, ranges = IRanges(start = starts, end = ends))
  # } else {
  #   starts <- cumsum(c(paf.aln$q.start, query.aln.widths))[-(length(query.aln.widths) + 1)]
  #   ends <- cumsum(c(paf.aln$q.start, query.aln.widths))[-1]
  #   query.gr <- GenomicRanges::GRanges(seqnames = paf.aln$q.name, ranges = IRanges(start = starts, end = ends))
  # }
  
  ## Separate inserted ranges
  if (length(query.ins.gr) > 0) {
    query.ins.gr <- query.gr[remove.ins]
    query.gr <- query.gr[-remove.ins]
  }
  
  ## Adjust ranges start to match initial PAF alignment
  #target.gr <- resize(target.gr, width(target.gr) + 1, fix = 'end')
  #query.gr <- resize(query.gr, width(query.gr) + 1, fix = 'end')
  
  ## Subset the cigar string per target region
  starts <- as.numeric(start(match.gr))
  ends <- as.numeric(end(match.gr))
  cigars.region <- S4Vectors::sapply(1:length(starts), function(i) GenomicAlignments::cigarNarrow(cigar = paf.aln$cg, start = starts[i], end = ends[i])[1])
  ## Get alignment length from regional cigars
  aln.len <- S4Vectors::sapply(cigars.region, function(cg) sum(GenomicAlignments::explodeCigarOpLengths(cigar = cg)[[1]]), USE.NAMES = FALSE)
  ## Get matched bases from regional cigars
  n.match <- S4Vectors::sapply(cigars.region, function(cg) sum(GenomicAlignments::explodeCigarOpLengths(cigar = cg, ops = c('=','M'))[[1]]), USE.NAMES = FALSE)
  ## Update paf alignment
  match.paf.aln <- dplyr::tibble( q.name=as.character(seqnames(query.gr)),
                                  q.len=paf.aln$q.len,
                                  q.start=as.numeric(start(query.gr)),
                                  q.end=as.numeric(end(query.gr)),
                                  strand=paf.aln$strand,
                                  t.name=as.character(seqnames(target.gr)),
                                  t.len=paf.aln$t.len,
                                  t.start=as.numeric(start(target.gr)),
                                  t.end=as.numeric(end(target.gr)),
                                  n.match=n.match,
                                  aln.len=aln.len,
                                  mapq=paf.aln$mapq,
                                  cg=cigars.region)
  
  
  ## Report SV regions
  sv.paf.aln <- NULL
  if (report.sv) {
    ## Report insertions
    if (length(query.ins.gr) > 0) {
      ## Convert to query or target coordinates
      target.ins.gr <- GenomicRanges::shift(target.ins.gr, shift = paf.aln$t.start)
      #query.ins.gr <- GenomicRanges::shift(query.ins.gr[,0], shift = paf.aln$q.start)
      ## Create insertion ID
      ins.cg <- paste0(width(query.ins.gr) - 1, 'I')
      ## Create paf alignment
      ins.paf.aln <- dplyr::tibble( q.name=as.character(seqnames(query.ins.gr)),
                                    q.len=paf.aln$q.len,
                                    q.start=as.numeric(start(query.ins.gr)),
                                    q.end=as.numeric(end(query.ins.gr)),
                                    strand=paf.aln$strand,
                                    t.name=as.character(seqnames(target.ins.gr)),
                                    t.len=paf.aln$t.len,
                                    t.start=as.numeric(start(target.ins.gr)),
                                    t.end=as.numeric(end(target.ins.gr)),
                                    n.match=0,
                                    aln.len=0,
                                    mapq=paf.aln$mapq,
                                    cg=ins.cg)
      sv.paf.aln[[length(sv.paf.aln) + 1]] <- ins.paf.aln
    } 
    ## Report deletions
    if (length(target.del.gr) > 0) {
      ## Convert to query or target coordinates
      target.del.gr <- GenomicRanges::shift(target.del.gr, shift = paf.aln$t.start)
      #end(query.del.gr) <- seqlengths(query.match.gr) - end(query.del.gr)
      #start(query.del.gr) <- seqlengths(query.match.gr) - start(query.del.gr)
      if (paf.aln$strand == '-') {
        query.del.gr <- mirrorRanges(gr = query.del.gr, seqlength = seqlengths(query.match.gr))
        suppressWarnings( query.del.gr <- GenomicRanges::shift(query.del.gr[,0], shift = paf.aln$q.start) )
      } else {  
        suppressWarnings( query.del.gr <- GenomicRanges::shift(query.del.gr, shift = paf.aln$q.start) )
      }  
      ## Create deletion ID
      del.cg <- paste0(width(target.del.gr) - width(query.del.gr), 'D')
      ## Create paf alignment
      del.paf.aln <- dplyr::tibble( q.name=as.character(seqnames(query.del.gr)),
                                    q.len=paf.aln$q.len,
                                    q.start=as.numeric(start(query.del.gr)),
                                    q.end=as.numeric(end(query.del.gr)),
                                    strand=paf.aln$strand,
                                    t.name=as.character(seqnames(target.del.gr)),
                                    t.len=paf.aln$t.len,
                                    t.start=as.numeric(start(target.del.gr)),
                                    t.end=as.numeric(end(target.del.gr)),
                                    n.match=0,
                                    aln.len=0,
                                    mapq=paf.aln$mapq,
                                    cg=del.cg)
      sv.paf.aln[[length(sv.paf.aln) + 1]] <- del.paf.aln
    }
    if (length(sv.paf.aln) > 0) {
      sv.paf.aln <- dplyr::bind_rows(sv.paf.aln)
    }  
  } else {
    sv.paf.aln <- NULL
  }  
  
  ## Return broken paf alignments
  if (report.sv) {
    return(list('M' = match.paf.aln, 'SVs' = sv.paf.aln))
  } else {
    return(list('M' = match.paf.aln, 'SVs' = NULL))
  }  
}



#' A wrapper function for \code{\link{breakPafAlignment}} expanding multiple PAF alignments 
#' into a set of matching bases between query and target sequence.
#'
#' @param paf.table A \code{data.frame} or \code{tibble} containing a single PAF record with 12 mandatory columns.
#' @inheritParams cigar2ranges
#' @inheritParams breakPafAlignment
#' @importFrom dplyr bind_rows
#' @return A \code{list} of \code{tibble} objects storing matched ('M') alignments as well as structurally vairable ('SV') bases if 'report.sv' is TRUE.
#' @author David Porubsky
#' @export
#'
breakPaf <- function(paf.table=NULL, min.deletion.size=50, min.insertion.size=50, collapse.mismatches=TRUE, report.sv=TRUE) {
  ## Extract matching bases for each PAF record
  matches <- list()
  svs <- list()
  for (i in 1:nrow(paf.table)) {
    paf.aln <- paf.table[i,]
    paf.aln.exp <- breakPafAlignment(paf.aln = paf.aln, 
                      min.deletion.size = min.deletion.size, 
                      min.insertion.size = min.insertion.size, 
                      collapse.mismatches = collapse.mismatches, 
                      report.sv = report.sv)
    paf.aln.exp$M$aln.id <- i
    matches[[i]] <- paf.aln.exp$M
    if (!is.null(paf.aln.exp$SVs)) {
      paf.aln.exp$SVs$aln.id <- i
      svs[[length(svs) + 1]] <- paf.aln.exp$SVs
    }  
  }
  ## Return broken paf alignments
  if (report.sv) {
    return(list('M' = dplyr::bind_rows(matches), 'SVs' = dplyr::bind_rows(svs)))
  } else {
    return(list('M' = dplyr::bind_rows(matches), 'SVs' = NULL))
  }  
}

