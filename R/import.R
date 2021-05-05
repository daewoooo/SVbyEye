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
  if (file.exists(paf.file)) {
    message("Loading PAF file: ", paf.file)
    paf <- utils::read.table(paf.file, stringsAsFactors = FALSE)
    ## Keep only first 12 columns
    paf <- paf[,c(1:12)]
    ## Add header
    header <- c('q.name', 'q.len', 'q.start', 'q.end', 'strand', 't.name', 't.len', 't.start', 't.end', 'n.match', 'aln.len', 'mapq') 
    colnames(paf) <- header
  } else {
    stop(paste0("PAF file ", bedfile, " doesn't exists !!!"))
  }  
  
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

#' Load a de novo assembly aligned to the reference in BED format into a \code{\link{GRanges-class}} object.
#' 
#' This function will take a BED file of contig alignments to a reference genome and converts them 
#' into a \code{\link{GRanges-class}} object. This function also allows to collapse single unique contigs
#' that are aligned in multiple pieces after the alignment to the reference.
#'
#' @param bed.file A BED file of contig alignments to a reference genome.
#' @param index A unique identifier to be added as an 'ID' field. 
#' @param min.mapq A minimum mapping quality of alignments reported in submitted BED file.
#' @param min.aln.width A minimum width of an aligned sequence to a reference genome.
#' @param max.gap A maximum length of a gap within a single contig alignments to be collapsed.
#' @param report.mapped.ends Set to \code{TRUE} if alignment ends of each contig to the reference should reported as well.
#' @return A \code{\link{GRanges-class}} object.
#' @import GenomicRanges
#' @import GenomeInfoDb
#' @author David Porubsky
#' @export
#' 
bed2ranges <- function(bed.file=NULL, index=NULL, min.mapq=10, min.aln.width=10000, max.gap=100000, report.mapped.ends=FALSE) {
  if (file.exists(bed.file)) {
    message("Loading BED file: ", bed.file)
    bed.df <- utils::read.table(file = bed.file, header = FALSE, stringsAsFactors = FALSE)
    ## Add column names
    colnames(bed.df) <- c('seqnames','start','end','ctg','mapq','strand')[1:ncol(bed.df)]
  } else {
    stop(paste0("BED file ", bed.file, " doesn't exists !!!"))
  }  
  ## Add index if defined
  if (!is.null(index) & is.character(index)) {
    bed.df$ID <- index
  }
  ## Filter by mapping quality
  if (min.mapq > 0) {
    bed.df <- bed.df[bed.df$mapq >= min.mapq,]
  }
  if (nrow(bed.df) == 0) {
    stop("None of the BED alignments reach user defined mapping quality (min.mapq) !!!")
  }
  ## Convert data.frame to GRanges object
  bed.gr <- GenomicRanges::makeGRangesFromDataFrame(bed.df, keep.extra.columns = TRUE)
  ## Ignore strand
  #GenomicRanges::strand(bed.gr) <- '*'
  ## Filter out small alignments
  if (min.aln.width > 0) {
    message("Keeping alignments of min width: ", min.aln.width, 'bp')
    bed.gr <- bed.gr[width(bed.gr) >= min.aln.width]
  }
  ## Keep only seqlevels that remained after data filtering
  if (length(bed.gr) > 0) {
    bed.gr <- GenomeInfoDb::keepSeqlevels(bed.gr, value = unique(as.character(GenomeInfoDb::seqnames(bed.gr))))
  } else {
    stop("None of the BED alignments reach user defined alignment size (min.aln.width !!!")
  }
  
  ## Split ranges per contig
  bed.grl <- GenomicRanges::split(bed.gr, bed.gr$ctg)
  ## Make sure data frame is sorted
  #bed.gr <- GenomicRanges::sort(bed.gr, ignore.strand=TRUE)
  ## Collapse merge splitted continuous alignments of the same contig
  #bed.gr <- primatR::collapseBins(bed.gr, id.field = 1)
  #bed.grl <- endoapply(bed.grl, function(gr) primatR::collapseBins(gr, id.field = 1))
  if (max.gap > 0) {
    message("Filling gaps of max size: ", max.gap, 'bp')
    bed.gr <- collapseGaps(ranges = bed.grl, max.gap = max.gap)
  }
  
  ## Report end positions of each contig
  if (report.mapped.ends) {
    message("Reporting mapped range for each contig")
    bed.grl <- GenomicRanges::split(bed.gr, bed.gr$ctg)
    
    to.collapse <- which(lengths(bed.grl) > 1)
    ctg.ends.grl <- S4Vectors::endoapply(bed.grl[to.collapse], range)
    
    simple.ends.gr <- unlist(bed.grl[-to.collapse], use.names = FALSE)
    simple.ends <- as.character(simple.ends.gr)
    names(simple.ends) <- unique(as.character(GenomeInfoDb::seqnames(simple.ends.gr)))
    
    split.ends <- as.character(unlist(ctg.ends.grl))
    split.ends.l <- split(split.ends, names(split.ends))
    split.ends <- sapply(split.ends.l, function(x) paste(x, collapse = ';'))
    
    ctg.ends <- c(simple.ends, split.ends)
    bed.gr$ctg.end.pos <- ctg.ends[match(bed.gr$ctg, names(ctg.ends))]
  }
  names(bed.gr) <- NULL
  return(bed.gr)
}


#' Load a de novo assembly aligned to the reference in PAF format into a \code{\link{GRanges-class}} object.
#' 
#' This function will take a PAF file of contig alignments to a reference genome and converts them 
#' into a \code{\link{GRanges-class}} object.
#'
#' @param paf.file A BED file of contig alignments to a reference genome.
#' @param min.ctg.size A minimum length a final contig after gaps are collapsed.
#' @param report.ctg.ends Set to \code{TRUE} if aligned position of each contig ends should be reported
#' @return A \code{\link{GRanges-class}} object.
#' @import GenomicRanges
#' @import GenomeInfoDb
#' @author David Porubsky
#' @export
#' 
paf2ranges <- function(paf.file=NULL, index=NULL, min.mapq=10, min.aln.width=10000, min.ctg.size=500000, report.ctg.ends=FALSE) {
  if (file.exists(paf.file)) {
    message("Loading PAF file: ", paf.file)
    paf <- utils::read.table(paf.file, stringsAsFactors = FALSE)
    ## Keep only first 12 columns
    paf <- paf[,c(1:12)]
    ## Add header
    header <- c('q.name', 'q.len', 'q.start', 'q.end', 'strand', 't.name', 't.len', 't.start', 't.end', 'n.match', 'aln.len', 'mapq') 
    colnames(paf) <- header
  } else {
    stop(paste0("PAF file ", bedfile, " doesn't exists !!!"))
  }  
  ## Filter by mapping quality
  if (min.mapq > 0) {
    paf <- paf[paf$mapq >= min.mapq,]
  }
  if (nrow(paf) == 0) {
    stop("None of the PAF alignments reach user defined mapping quality (min.mapq) !!!")
  }
  ## Convert data.frame to GRanges object
  paf.gr <- GenomicRanges::makeGRangesFromDataFrame(paf, keep.extra.columns = FALSE, seqnames.field = 'q.name', start.field = 'q.start', end.field = 'q.end')
  mcols(paf.gr) <- paf[,c('q.len', 't.len', 'n.match', 'aln.len', 'mapq')]
  names(paf.gr) <- NULL
  paf.gr$target.gr <- GenomicRanges::GRanges(seqnames=paf$t.name, ranges=IRanges::IRanges(start=paf$t.start, end=paf$t.end), 'q.name'= paf$q.name)
  ## Add index if defined
  if (!is.null(index) & is.character(index)) {
    paf.gr$ID <- index
  }
  ## Ignore strand
  #GenomicRanges::strand(paf.gr) <- '*'
  ## Filter out small alignments
  if (min.aln.width > 0) {
    message("Keeping alignments of min width: ", min.aln.width, 'bp')
    paf.gr <- paf.gr[width(paf.gr) >= min.aln.width]
  }
  if (length(paf.gr) == 0) {
    stop("None of the PAF alignments reach user defined alignment size (min.aln.width !!!")
  }
  ## Filter out small contig sizes
  if (min.ctg.size > 0) {
    message("Keeping contigs of min size: ", min.ctg.size, 'bp')
    paf.gr <- paf.gr[paf.gr$q.len >= min.ctg.size]
  }
  ## Keep only seqlevels that remained after data filtering
  if (length(paf.gr) > 0) {
    paf.gr <- GenomeInfoDb::keepSeqlevels(paf.gr, value = unique(as.character(GenomeInfoDb::seqnames(paf.gr))))
  } else {
    stop("None of the PAF alignments reach user defined contig size (min.ctg.size) !!!")
  }
  
  if (report.ctg.ends == TRUE) {
    message("Reporting end positions for each contig")
    paf.grl <- split(paf.gr, GenomeInfoDb::seqnames(paf.gr))
    to.collapse <- which(lengths(paf.grl) > 1)
    ctg.ends.grl <- S4Vectors::endoapply( paf.grl[to.collapse], function(x) range(sort(x, ignore.strand=TRUE)$target.gr) )
    
    simple.ends.gr <- unlist(paf.grl[-to.collapse], use.names = FALSE)
    simple.ends <- as.character(simple.ends.gr$target.gr)
    names(simple.ends) <- unique(as.character(GenomeInfoDb::seqnames(simple.ends.gr)))
    
    split.ends <- as.character(unlist(ctg.ends.grl))
    split.ends.l <- split(split.ends, names(split.ends))
    split.ends <- sapply(split.ends.l, function(x) paste(x, collapse = ';'))
    
    ctg.ends <- c(simple.ends, split.ends)
    paf.gr$ctg.end.pos <- ctg.ends[match(as.character(seqnames(paf.gr)), names(ctg.ends))]
  }
  return(paf.gr)
}  

