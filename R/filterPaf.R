#' Filter PAF alignments
#' 
#' This function takes loaded PAF alignments from \function{readPaf} and perform
#' user defined filtering of input alignments based on user defined region, mapping quality, and
#' alignment length (etc.).
#'
#' @param min.mapq Minimum mapping quality to retain.
#' @param min.align.len Minimum alignment length to retain.
#' @param min.align.n Minimum number of alignments between a pair of sequences/regions.
#' @param target.region A user defined target region to load in a character string ('chr:start-end') or as
#' a \code{\link{GRanges-class}} object containing a single genomic region.
#' @param query.region A user defined query region to load in a character string ('chr:start-end') or as
#' a \code{\link{GRanges-class}} object containing a single genomic region.
#' @param drop.self.align Remove alignments of a given sequence to itself.
#' @inheritParams breakPaf
#' @return A \code{tibble} of filtered PAF alignments.
#' @importFrom dplyr group_by mutate n
#' @importFrom utils read.table
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @importFrom S4Vectors queryHits
#' @importFrom tibble is_tibble as_tibble 
#' @author David Porubsky
#' @export
filterPaf <- function(paf.table, min.mapq=10, min.align.len=1000, min.align.n=1, target.region=NULL, query.region=NULL, drop.self.align=TRUE) {
  ## Check user input ##
  ## Make sure PAF has at least 12 mandatory fields
  if (ncol(paf.table) >= 12) {
    paf <- paf.table
  } else {
    stop('Submitted PAF alignments do not contain a minimum of 12 mandatory fields, see PAF file format definition !!!')
  }
  
  ## Filer alignments by target region
  if (!is.null(target.region)) {
    if (is.character(target.region)) {
      target.region.gr <- as(target.region, 'GRanges')
    } else if (class(target.region.gr) == 'GRanges') {
      target.region.gr <- target.region
    } else {
      message("Parameter 'target.region' can either be 'GRanges' object or character string 'chr#:start-end'!!!")
    }
    ## Subset PAF by ranges overlaps
    target.gr <- GenomicRanges::makeGRangesFromDataFrame(paf, seqnames.field = 't.name', start.field = 't.start', end.field = 't.end')
    hits <- GenomicRanges::findOverlaps(target.gr, target.region.gr) ## TODO Take only overlaps, regions itself are not shrinked!!!
    if (length(hits) > 0) {
      paf <- paf[S4Vectors::queryHits(hits),]
      #paf <- paf[paf$t.name %in% as.character(seqnames(target.region.gr)) & paf$t.start >= start(target.region.gr) & paf$t.end <= end(target.region.gr),]
    } else {
      warning("None of the PAF ranges overlap user defined 'target.region', skipping ...")
    }  
  }
  
  ## Filer alignments by query region
  if (!is.null(query.region)) {
    if (is.character(query.region)) {
      query.region.gr <- as(query.region, 'GRanges')
    } else if (class(query.region.gr) == 'GRanges') {
      query.region.gr <- query.region
    } else {
      message("Parameter 'query.region' can either be 'GRanges' object or character string 'chr#:start-end'!!!")
    }
    ## Subset PAF by ranges overlaps
    query.gr <- GenomicRanges::makeGRangesFromDataFrame(paf, seqnames.field = 't.name', start.field = 't.start', end.field = 't.end')
    hits <- GenomicRanges::findOverlaps(query.gr, query.region.gr) ## TODO Take only overlaps, regions itself are not shrinked!!!
    if (length(hits) > 0) {
      paf <- paf[S4Vectors::queryHits(hits),]
      #paf <- paf[paf$t.name %in% as.character(seqnames(target.region.gr)) & paf$t.start >= start(target.region.gr) & paf$t.end <= end(target.region.gr),]
    } else {
      warning("None of the PAF ranges overlap user defined 'query.region', skipping ...")
    }  
  }
  
  ## Filter by the number of alignments per sequence alignment pair
  ## Get unique alignment ID
  paf$seq.pair <- paste0(paf$q.name, '__', paf$t.name)
  if (min.align.n > 0) {
    paf <- paf %>% dplyr::group_by(seq.pair) %>% dplyr::mutate(align.n = dplyr::n())
    paf <- paf[paf$align.n >= min.align.n,]
  }
  
  ## Keep only specific sequence/region name
  # @param seqname.grep Retain only a specific sequence/region name with a given character string.
  # if (!is.null(seqnames.grep)) {
  #   paf <- paf[grep(paf$q.name, pattern = seqnames.grep),]
  #   paf <- paf[grep(paf$t.name, pattern = seqnames.grep),]
  # }
  
  ## Filter by mapping quality
  ## Make sure mapping quality is defined
  if (all(is.na(paf$mapq))) {
    paf$mapq[is.na(paf$mapq)] <- min.mapq
    if (min.mapq > 0 & is.numeric(paf$mapq)) {
      paf <- paf[paf$mapq >= min.mapq,]
    }
  }  
  
  ## Filter by alignment length
  if (all(is.na(paf$aln.len))) {
    paf$aln.len[is.na(paf$aln.len)] <- min.align.len
    if (min.align.len > 0 & is.numeric(paf$aln.len)) {
      paf <- paf[paf$aln.len >= min.align.len,]
    }
  }  
  
  ## Filter out self-alignments
  if (drop.self.align) {
    paf <- paf[!(paf$q.name == paf$t.name),]
  }
  
  ## Export data object
  if (tibble::is_tibble(paf)) {
    return(paf)
  } else {
    return(tibble::as_tibble(paf))
  }  
}