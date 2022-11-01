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
#' @param min.ctg.ends A minimum length of alignment to be considered when reporting contig end alignments.
#' @return A \code{\link{GRanges-class}} object.
#' @import GenomicRanges
#' @import GenomeInfoDb
#' @author David Porubsky
#' @export
#' 
paf2ranges <- function(paf.file=NULL, index=NULL, min.mapq=10, min.aln.width=10000, min.ctg.size=500000, report.ctg.ends=FALSE, min.ctg.ends=50000) {
  ## Get total processing time
  #ptm <- proc.time()
  
  if (file.exists(paf.file)) {
    ptm <- startTimedMessage("\nLoading PAF file: ", paf.file)
    paf <- utils::read.table(paf.file, stringsAsFactors = FALSE, comment.char = '&')
    ## Keep only first 12 columns
    paf <- paf[,c(1:12)]
    ## Add header
    header <- c('q.name', 'q.len', 'q.start', 'q.end', 'strand', 't.name', 't.len', 't.start', 't.end', 'n.match', 'aln.len', 'mapq') 
    colnames(paf) <- header
    stopTimedMessage(ptm)
  } else {
    stop(paste0("PAF file ", bedfile, " doesn't exists !!!"))
  }  
  ## Filter by mapping quality
  if (min.mapq > 0) {
    ptm <- startTimedMessage("    Keeping alignments of min.mapq: ", min.mapq)
    paf <- paf[paf$mapq >= min.mapq,]
    stopTimedMessage(ptm)
  }
  if (nrow(paf) == 0) {
    stop("None of the PAF alignments reach user defined mapping quality (min.mapq) !!!")
  }
  ## Convert data.frame to GRanges object
  #paf.gr <- GenomicRanges::makeGRangesFromDataFrame(paf, keep.extra.columns = FALSE, seqnames.field = 'q.name', start.field = 'q.start', end.field = 'q.end')
  paf.gr <- GenomicRanges::makeGRangesFromDataFrame(paf, keep.extra.columns = FALSE, seqnames.field = 't.name', start.field = 't.start', end.field = 't.end')
  mcols(paf.gr) <- paf[,c('q.len', 't.len', 'n.match', 'aln.len', 'mapq', 'q.name')]
  names(paf.gr) <- NULL
  #paf.gr$target.gr <- GenomicRanges::GRanges(seqnames=paf$t.name, ranges=IRanges::IRanges(start=paf$t.start, end=paf$t.end), 'q.name'= paf$q.name)
  paf.gr$query.gr <- GenomicRanges::GRanges(seqnames=paf$q.name, ranges=IRanges::IRanges(start=paf$q.start, end=paf$q.end), 't.name'= paf$t.name)
  ## Add index if defined
  if (!is.null(index) & is.character(index)) {
    paf.gr$ID <- index
  }
  ## Ignore strand
  #GenomicRanges::strand(paf.gr) <- '*'
  ## Filter out small alignments
  if (min.aln.width > 0) {
    ptm <- startTimedMessage("    Keeping alignments of min width: ", min.aln.width, 'bp')
    #paf.gr <- paf.gr[width(paf.gr) >= min.aln.width]
    paf.gr <- paf.gr[paf.gr$aln.len >= min.aln.width]
    stopTimedMessage(ptm)
  }
  if (length(paf.gr) == 0) {
    stop("None of the PAF alignments reach user defined alignment size (min.aln.width !!!")
  }
  ## Filter out small contig sizes
  if (min.ctg.size > 0) {
    ptm <- startTimedMessage("    Keeping contigs of min size: ", min.ctg.size, 'bp')
    paf.gr <- paf.gr[paf.gr$q.len >= min.ctg.size]
    stopTimedMessage(ptm)
  }
  ## Keep only seqlevels that remained after data filtering
  if (length(paf.gr) > 0) {
    paf.gr <- GenomeInfoDb::keepSeqlevels(paf.gr, value = unique(as.character(GenomeInfoDb::seqnames(paf.gr))))
    paf.gr$query.gr <- GenomeInfoDb::keepSeqlevels(paf.gr$query.gr, value = unique(as.character(GenomeInfoDb::seqnames(paf.gr$query.gr))))
  } else {
    stop("None of the PAF alignments reach user defined contig size (min.ctg.size) !!!")
  }
  
  if (report.ctg.ends == TRUE) {
    ptm <- startTimedMessage("    Reporting end positions of each contig")
    #paf.grl <- split(paf.gr, GenomeInfoDb::seqnames(paf.gr))
    paf.grl <- split(paf.gr, GenomeInfoDb::seqnames(paf.gr$query.gr))
    to.collapse <- which(lengths(paf.grl) > 1)
    #ctg.ends.grl <- S4Vectors::endoapply( paf.grl[to.collapse], function(x) range(sort(x, ignore.strand=TRUE)$target.gr) )
    
    ## Helper function
    gr2ranges <- function(gr, min.ctg.ends=0, allowed.ctg.length=0.05) {
      gr <- GenomeInfoDb::keepSeqlevels(gr, value = unique(as.character(seqnames(gr))))
      q.name <- unique(gr$q.name)
      if (length(q.name) > 1) {
        stop("Submitted GRanges object contains more than one 'q.name', not allowed !!!")
      }
      if (min.ctg.ends > 0) {
        gr <- gr[gr$aln.len >= min.ctg.ends]
      }
      if (length(gr) > 0) {
        ## Order by contig alignments
        gr <- gr[order(gr$query.gr)]
        if (length(gr) > 1) {
          ## Report target ranges corresponding to the contig ends
          gr.ends <- gr[c(1, length(gr))]
          ## Report genomic range for contig ends mapping
          gr.range <- range(gr.ends, ignore.strand=TRUE)
          ## Get best contiguous alignment ##
          ## Report target sequence with the highest number of continuously aligned bases
          target2keep <- names(which.max(sum(split(width(gr), seqnames(gr)))))
          gr.contig <- gr[seqnames(gr) == target2keep]
          gr.contig <- gr.contig[,0]
          strand(gr.contig) <- '*'
          gr.gaps <- gaps(gr.contig, start = min(start(gr.contig)))
          if (length(gr.gaps) > 0) {
            ## Remove gaps longer than query(contig) length
            q.len <- unique(gr$q.len)
            ## Remove gaps longer than total alignment length
            #aln.len <- sum(gr$aln.len)
            ## Adjust allowed contig length based on fraction allowed to be added to the contigs size
            if (!is.null(allowed.ctg.length)) {
              if (allowed.ctg.length > 0) {
                q.len <- q.len + (q.len * allowed.ctg.length)
                #aln.len <- aln.len + (aln.len * allowed.ctg.length)
              }
            }
            gr.gaps <- gr.gaps[width(gr.gaps) < q.len]
            #gr.gaps <- gr.gaps[width(gr.gaps) < aln.len]
            ## Keep only gaps that together with contig length are no longer than query(contig) length
            ctg.size <- sum(width(reduce(gr.contig)))
            gap.size <- sum(width(reduce(gr.gaps)))  
            while ((ctg.size + gap.size) > q.len & length(gr.gaps) > 0) {
            #while ((ctg.size + gap.size) > aln.len & length(gr.gaps) > 0) {  
              ## Remove longest gap
              gr.gaps <- gr.gaps[-which.max(width(gr.gaps))]
              gap.size <- sum(width(gr.gaps))
            }
          }  
          ## Allow only gaps that are no longer than smallest alignment
          #gr.gaps <- gr.gaps[width(gr.gaps) <= min(gr$aln.len)]
          ## Collapse contiguous ranges
          gr.contig <- reduce(c(gr.contig, gr.gaps))
          ## Add contig names and split alignment ids
          if (length(gr.contig) > 1) {
            gr.contig$q.name <- q.name
            aln <- paste0('p', 1:length(gr.contig))
            gr.contig$aln <- aln
          } else {
            gr.contig$q.name <- q.name
            gr.contig$aln <- 'single'
          }  
          #gr.contig <- gr.contig[which.max(width(gr.contig))]
          # gr.ranges <- range(gr[order(gr$query.gr)], ignore.strand=TRUE)
          # if (report.longest.range) {
          #   return(gr.ranges[which.max(width(gr.ranges))])
          # } else {
          #   return(gr.ranges)
          # }
          gr.contig$ends <- paste(as.character(gr.range), collapse = ';')
          #gr.range$longest.aln <- gr.contig
        } else {
          #gr.range <- gr[,0]
          #strand(gr.range) <- '*'
          #gr.range$longest.aln <- gr.range
          gr.contig <- gr[,0]
          strand(gr.contig) <- '*'
          gr.contig$q.name <- q.name
          gr.contig$aln <- 'single'
          gr.contig$ends <- paste(as.character(gr.contig), collapse = ';')
        }  
      } else {
        #gr.range <- GRanges()
        gr.contig <- GRanges()
      }  
      #return(gr.range)
      return(gr.contig)
    }
    
    #ctg.ends.grl <- S4Vectors::endoapply( paf.grl[to.collapse], function(x) range(x[order(x$query.gr)], ignore.strand=TRUE) )
    
    ## Report contigs alignment ranges only for alignments of a certain size (default: 50kb)
    if (min.ctg.ends > 0) {
      ctg.ends.grl <- S4Vectors::endoapply( paf.grl[to.collapse], 
                                            function(x) gr2ranges(gr = x, min.ctg.ends = min.ctg.ends) )  
    } else {
      ctg.ends.grl <- S4Vectors::endoapply( paf.grl[to.collapse], 
                                            function(x) gr2ranges(gr = x) )
    }
    
    simple.ends.gr <- unlist(paf.grl[-to.collapse], use.names = FALSE)
    #simple.ends <- as.character(simple.ends.gr$target.gr)
    simple.ends <- as.character(simple.ends.gr)
    #names(simple.ends) <- unique(as.character(GenomeInfoDb::seqnames(simple.ends.gr)))
    names(simple.ends) <- unique(as.character(GenomeInfoDb::seqnames(simple.ends.gr$query.gr)))
    
    simple.gr <- simple.ends.gr[,0]
    strand(simple.gr) <- '*'
    #simple.gr$longest.aln <- simple.gr
    simple.gr$q.name <- as.character(seqnames(simple.ends.gr$query.gr))
    simple.gr$aln <- 'single'
    simple.gr$ends <- as.character(simple.gr)
    names(simple.gr) <- NULL
    
    split.ends <- as.character(unlist(ctg.ends.grl))
    split.ends.l <- split(split.ends, names(split.ends))
    split.ends <- sapply(split.ends.l, function(x) paste(x, collapse = ';'))
    
    split.ends.gr <- unlist(ctg.ends.grl)
    #split.ends.gr$q.name <- names(split.ends.gr)
    names(split.ends.gr) <- NULL
    
    #ctg.ends <- c(simple.ends, split.ends)
    #paf.gr$ctg.end.pos <- ctg.ends[match(as.character(seqnames(paf.gr)), names(ctg.ends))]
    #paf.gr$ctg.end.pos <- ctg.ends[match(as.character(seqnames(paf.gr$query.gr)), names(ctg.ends))]
    
    ends.gr <- sort(c(simple.gr, split.ends.gr))
    ends.gr$ID <- unique(paf.gr$ID)
    
    stopTimedMessage(ptm)
    return(list('ctg.aln'=paf.gr, 'ctg.ends'=ends.gr))
  } else {
    return(list('ctg.aln'=paf.gr, 'ctg.ends'=NULL))
  }
  
  ## Report total processing time
  #time <- proc.time() - ptm
  #message("Total time: ", round(time[3],2), "s")
  
  #return(paf.gr)
}  

