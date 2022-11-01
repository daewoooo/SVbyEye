#' Read and filter PAF input file
#' 
#' This function takes PAF output file from minimap2 alignments, loads the file and
#' perform user defined filtering of input alignments based on mapping quality and
#' alignment length.
#'
#' @inheritParams breakPafAlignment
#' @return A \code{data.frame} of PAF alignments reported as x and y coordinate values.
#' @author David Porubsky
#' @export
paf2coords <- function(paf.table) {
  ## Check user input ##
  ## Make sure PAF has at least 12 mandatory fields
  if (ncol(paf.table) >= 12) {
    paf <- paf.table
  } else {
    stop('Submitted PAF alignments do not contain a minimum of 12 mandatory fields, see PAF file format definition !!!')
  }
  
  ## Flip start-end if strand == '-'
  paf[paf$strand == '-', c('t.start','t.end')] <- rev(paf[paf$strand == '-', c('t.start','t.end')])
  #paf[paf$strand == '-', c('q.start','q.end')] <- rev(paf[paf$strand == '-', c('q.start','q.end')])
  
  ## Get unique alignment ID
  if (!'seq.pair' %in% colnames(paf)) {
    paf$seq.pair <- paste0(paf$q.name, '__', paf$t.name)
  }  
  
  ## Sync scales between alignments [per region id]
  paf.l <- split(paf, paf$seq.pair)
  for (i in seq_along(paf.l)) {
    paf.sub <- paf.l[[i]]
    q.range <- range(c(paf.sub$q.start, paf.sub$q.end))
    t.range <- range(c(paf.sub$t.start, paf.sub$t.end))
    ## Adjust target ranges given the size difference with respect query ranges
    range.offset <- diff(q.range) - diff(t.range)
    t.range[2] <- t.range[2] + range.offset ## Make a start position as offset and change only end position
    ## Covert query to target coordinates
    paf.sub$q.start.trans <- q2t(x = paf.sub$q.start, q.range = q.range, t.range = t.range)
    paf.sub$q.end.trans <- q2t(x = paf.sub$q.end, q.range = q.range, t.range = t.range)
    # q.range <- range(c(paf$q.start, paf$q.end))
    # t.range <- range(c(paf$t.start, paf$t.end))
    # paf$q.start.trans <- q2t(x = paf$q.start, q.range = q.range, t.range = t.range)
    # paf$q.end.trans <- q2t(x = paf$q.end, q.range = q.range, t.range = t.range)
    paf.l[[i]] <- paf.sub
  }  
  paf <- do.call(rbind, paf.l)
  
  ## Vectorize data transformation ##
  x <- c(rbind(paf$q.start.trans, paf$t.start, paf$q.end.trans, paf$t.end))
  y <- rep(c(1,2,1,2), times=nrow(paf))
  group <- rep(1:nrow(paf), each=4)
  seq.name <- c(rbind(paf$q.name, paf$t.name, paf$q.name, paf$t.name))
  seq.pos <- c(rbind(paf$q.start, paf$t.start, paf$q.end, paf$t.end))
  seq.id <- c(rbind('query', 'target', 'query', 'target'))
  n.match <- rep(paf$n.match, each=4)
  aln.len <- rep(paf$aln.len, each=4)
  mapq <- rep(paf$mapq, each=4)
  align.id <- rep(paf$seq.pair, each=4)
  direction <- rep(paf$strand, each=4)
  
  coords <- data.frame(x=x, 
                       y=y, 
                       group=group, 
                       seq.pos=seq.pos,
                       direction=direction,
                       seq.name=seq.name, 
                       seq.id=seq.id,
                       n.match=n.match,
                       aln.len=aln.len,
                       mapq=mapq,
                       align.id=align.id,
                       stringsAsFactors = FALSE)
  return(coords)
}
