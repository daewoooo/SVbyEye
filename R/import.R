#' Read and filter PAF input file
#' 
#' This function takes PAF output file from minimap2 alignemts, loads the file and
#' perform user defined filtering of input alignements based on mapping quality and
#' alignment length.
#'
#'
#' @param paf.file A \code{data.frame} containing x and y coordinates.
#' @param min.mapq Minimum mapping quality to retain for PAF alignment file.
#' @param min.align.len Minimum alignment length to retain for PAF alignment file.
#' @param min.align.n Minimum number of alignments between a pair of sequences/regions.
#' @param seqname.grep Retain only a specific sequence/region name with a given character string.
#' @return A \code{vector} of rescaled coordinate values.
#' @importFrom dplyr group_by mutate
#' @importFrom utils read.table
#' @author David Porubsky
#' @export
paf2coords <- function(paf.file, min.mapq=10, min.align.len=1000, min.align.n=1, seqname.grep=NULL) {
  paf <- utils::read.table(paf.file, stringsAsFactors = FALSE)
  ## Keep only first 12 columns
  paf <- paf[,c(1:12)]
  ## Add header
  header <- c('q.name', 'q.len', 'q.start', 'q.end', 'strand', 't.name', 't.len', 't.start', 't.end', 'n.match', 'aln.len', 'mapq') 
  colnames(paf) <- header
  ## Flip start-end if strand == '-'
  paf[paf$strand == '-', c('t.start','t.end')] <- rev(paf[paf$strand == '-', c('t.start','t.end')])
  
  ## Get number of alignments per sequence pair
  if (min.align.n > 0) {
    paf$seq.pair <- paste0(paf$q.name, '__', paf$t.name)
    paf <- paf %>% dplyr::group_by(seq.pair) %>% dplyr::mutate(align.n = n())
    paf <- paf[paf$align.n >= min.align.n,]
  }
  
  ## Filter PAF alignments ##
  ## Keep only specific sequence/region name
  if (!is.null(seqname.grep)) {
    paf <- paf[grep(paf$q.name, pattern = seqname.grep),]
    paf <- paf[grep(paf$t.name, pattern = seqname.grep),]
  }
  ## Filter by mapping quality
  if (min.mapq > 0) {
    paf <- paf[paf$mapq >= min.mapq,]
  }
  ## Filter by alignment length
  if (min.align.len > 0) {
    paf <- paf[paf$aln.len >= min.align.len,]
  }
  
  ## Sync scales between alignments [per region id]
  paf.l <- split(paf, paf$seq.pair)
  for (i in seq_along(paf.l)) {
    paf.sub <- paf.l[[i]]
    q.range <- range(c(paf.sub$q.start, paf.sub$q.end))
    t.range <- range(c(paf.sub$t.start, paf.sub$t.end))
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
  #direction <- rep(paf$strand, each=4)
  
  coords <- data.frame(x=x, 
                       y=y, 
                       group=group, 
                       seq.pos=seq.pos,
                       #direction=direction,
                       seq.name=seq.name, 
                       seq.id=seq.id,
                       n.match=n.match,
                       aln.len=aln.len,
                       mapq=mapq,
                       align.id=align.id,
                       stringsAsFactors = FALSE)
  return(coords)
}
