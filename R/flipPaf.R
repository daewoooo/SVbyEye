#' Flip PAF alignments
#' 
#' This function takes loaded PAF alignments using \code{\link{readPaf}} function and flips 
#' the alignments given the desired majority strand orientation (Either '+' or '-').
#' 
#' @param force Set to \code{TRUE} if query PAF alignments should be flipped.
#' @inheritParams breakPaf
#' @inheritParams syncRangesDir
#' @return A \code{tibble} of PAF alignments flipped based on desired majority strand orientation.
#' @importFrom tibble is_tibble as_tibble
#' @author David Porubsky
#' @export
flipPaf <- function(paf.table, majority.strand=NULL, force=FALSE) {
  ## Helper function
  flipPafAln <- function(paf.table=NULL) {
    paf.flip <- paf.table
    ## Flip alignment orientation
    paf.flip$strand[paf.table$strand == '+'] <- '-'
    paf.flip$strand[paf.table$strand == '-'] <- '+'
    ## Keep info on original query coordinates [To test]
    #paf.flip$q.end <- paf.table$q.start
    #paf.flip$q.start <- paf.table$q.end
    ## Flip query coordinates
    paf.flip$q.end <- paf.table$q.len - paf.table$q.start
    paf.flip$q.start <- paf.table$q.len - paf.table$q.end
    ## Add info if the PAF alignment was flipped
    paf.flip$query.flip <- TRUE
    ## Make sure q.start position is always smaller than q.end
    #q.start <- pmin(paf.flip$q.start, paf.flip$q.end)
    #q.end <- pmax(paf.flip$q.start, paf.flip$q.end)
    #paf.flip$q.start <- q.start
    #paf.flip$q.end <- q.end
    
    return(paf.flip)
  }
  
  ## Check user input ##
  ## Make sure submitted paf.table has at least 12 mandatory fields
  if (ncol(paf.table) >= 12) {
    paf <- paf.table
    paf$query.flip <- FALSE
  } else {
    stop('Submitted PAF alignments do not contain a minimum of 12 mandatory fields, see PAF file format definition !!!')
  }
  
  ptm <- startTimedMessage("[flipPaf] Flipping orientation of PAF alignments")
  
  ## Force alignment flip
  if (force) {
    paf.l <- split(paf, paf$q.name)
    for (i in seq_along(paf.l)) {
      paf.ctg <- paf.l[[i]]
      paf.l[[i]] <- flipPafAln(paf.table = paf.ctg)
    }
    paf <- do.call(rbind, paf.l)
  }
  
  ## Sync by majority strand directionality
  if (!is.null(majority.strand) & force == FALSE) {
    ## Define majority and minority strand
    if (majority.strand == '+') {
      minority.strand = '-'
    } else if (majority.strand == '-') {
      minority.strand = '+'
    } else {
      stop("Parameter 'majority.strand' can only be defined as '+' or '-' !!!")
    }
    
    ## Make sure alignment is defined before deciding on majority directionality
    if (!is.na(sum(paf$aln.len)))  {
      if (sum(paf$aln.len) > 0) {
        paf.l <- split(paf, paf$q.name)
        for (i in seq_along(paf.l)) {
          paf.ctg <- paf.l[[i]]
          ## Flip directionality based to make sure majority strand covers the most bases
          majority.bases <- sum(paf.ctg$aln.len[paf.ctg$strand == majority.strand])
          minority.bases <- sum(paf.ctg$aln.len[paf.ctg$strand == minority.strand])
          if (majority.bases < minority.bases) {
            paf.l[[i]] <- flipPafAln(paf.table = paf.ctg)
            
          } else {
            paf.l[[i]] <- paf.ctg
          }
        }
        paf <- do.call(rbind, paf.l)
      }  
    }  
  }
  
  stopTimedMessage(ptm)
  ## Export data object
  if (tibble::is_tibble(paf)) {
    return(paf)
  } else {
    return(tibble::as_tibble(paf))
  } 
}
