#' Export FASTA sequences from a set of alignments reported in PAF formatted file.
#' 
#' @param paf.file paf.file A path to a PAF file containing alignments to be loaded.
#' @param bsgenome A \code{\link{BSgenome-class}} object of reference genome to get the genomic sequence from. 
#' @param asm.fasta An assembly FASTA file to extract DNA sequence determined by 'gr' parameter.
#' @param majority.strand A desired majority strand directionality to be reported.
#' @param fasta.save A path to a filename where to store final FASTA file.
#' @importFrom Rsamtools indexFa FaFile scanFa scanFaIndex
#' @importFrom BSgenome getSeq
#' @importFrom Biostrings writeXStringSet
#' @author David Porubsky
#' @export
#'
paf2FASTA <- function(paf.file, bsgenome=NULL, asm.fasta=NULL, majority.strand='+', fasta.save=NULL) {
  ## Load BSgenome object
  if (class(bsgenome) != 'BSgenome') {
    if (is.character(bsgenome)) {
      bsgenome <- tryCatch({
        suppressPackageStartupMessages(library(bsgenome, character.only=TRUE))
        bsgenome <- eval(parse(text=bsgenome)) ## replacing string by object
      }, error = function(e) {return(NULL)})
    }
  }
  ## Check if submitted fasta file is indexed
  if (!is.null(asm.fasta)) {	
    asm.fasta.idx <- paste0(asm.fasta, ".fai")
    if (!file.exists(asm.fasta.idx)) {
      fa.idx <- Rsamtools::indexFa(file = asm.fasta)
    }
  }
  
  ## Read in PAF file
  if (file.exists(paf.file)) {
    message("Loading PAF file: ", paf.file)
    paf <- tryCatch(
      #utils::read.table(paf.file, stringsAsFactors = FALSE, comment.char = '&', fill = TRUE), error = function(e) NULL
      readPaf(paf.file = paf.file, include.paf.tags = FALSE), error = function(e) NULL
      )
    # if (!is.null(paf)) {
    #   ## Keep only first 12 columns
    #   paf <- paf[,c(1:12)]
    #   ## Add header
    #   header <- c('q.name', 'q.len', 'q.start', 'q.end', 'strand', 't.name', 't.len', 't.start', 't.end', 'n.match', 'aln.len', 'mapq') 
    #   colnames(paf) <- header
    # }  
  } else {
    stop(paste0("PAF file ", paf.file, " doesn't exists !!!"))
  }  
  if (!is.null(paf)) {
    ## Convert query coordinates to GRanges
    paf.gr <- GenomicRanges::makeGRangesFromDataFrame(paf, seqnames.field = 'q.name', start.field = 'q.start', end.field = 'q.end')
    target.gr <- GenomicRanges::makeGRangesFromDataFrame(paf, seqnames.field = 't.name', start.field = 't.start', end.field = 't.end', strand='*')
    paf.gr$target.gr <- target.gr
    ## Make sure no genomic region starts with zero
    start(paf.gr) <- pmax(start(paf.gr), 1)
    
    ## Sync alignment directionality based on preferred majority strand
    paf.grl <- split(paf.gr, seqnames(paf.gr))
    for (i in seq_along(paf.grl)) {
      gr <- paf.grl[[i]]
      qy.red <- range(gr, ignore.strand=TRUE)
      tg.red <- range(gr$target.gr, ignore.strand=TRUE)
      if (majority.strand %in% c('+', '-')) {
        new.strand <- syncRangesDir(ranges = gr, majority.strand = majority.strand, strand.only = TRUE)
        if (all(new.strand != GenomicRanges::strand(gr))) {
          GenomicRanges::strand(gr) <- new.strand
          gr <- qy.red
          gr$target.gr <- tg.red
          gr$revcomp <- TRUE
          #reverseComp <- TRUE
        } else {
          #reverseComp <- FALSE
          gr <- qy.red
          gr$target.gr <- tg.red
          gr$revcomp <- FALSE
        }  
      } else {
        gr <- qy.red
        gr$target.gr <- tg.red
        gr$revcomp <- FALSE
        warning("Parameter 'majority.strand' can only takes values '+' or '-'!!!")
      }
      paf.grl[[i]] <- gr
    }
    paf.gr <- unlist(paf.grl, use.names = FALSE)
    
    ## Order regions by query position
    #paf.gr <- GenomicRanges::sort(paf.gr, ignore.strand=TRUE)
    ## Order regions by target position
    paf.gr <- paf.gr[GenomicRanges::order(paf.gr$target.gr)]
    ## Collapse consecutive alignments coming from the same contig/sequence
    #paf.gr$q.id <- as.character(GenomeInfoDb::seqnames(paf.gr))
    #paf.gr <- primatR::collapseBins(paf.gr, id.field = 3)
    ## Extract FASTA sequence
    if (!is.null(bsgenome)) {
      ## Extract FASTA from BSgenome object
      gr.seq <- BSgenome::getSeq(bsgenome, paf.gr)
      names(gr.seq) <- as.character(paf.gr)
    } else if (is.character(asm.fasta)) {
      ## Extract FASTA from user defined FASTA file
      fa.file <- open(Rsamtools::FaFile(asm.fasta))
      ## Remove sequences not present in the FASTA index
      fa.idx <- Rsamtools::scanFaIndex(fa.file)
      ## Make sure all sequence ranges are present in submitted FASTA sequences
      if(!all(as.character(GenomeInfoDb::seqnames(paf.gr)) %in% as.character(GenomeInfoDb::seqnames(fa.idx)))) {
        warning('Not all PAF ranges are present in the submitted FASTA file, subsetting !!!')
        paf.gr <- suppressWarnings( IRanges::subsetByOverlaps(paf.gr, fa.idx) )
        if (length(paf.gr) == 0) {
          stop('None of the PAF ranges present in the submitted FASTA file, likely wrong FASTA file submitted!!!')
        }
      } 
      ## Read in sequence for a given range(s)
      gr.seq <- Rsamtools::scanFa(file = fa.file, param = paf.gr, as = "DNAStringSet")
      ## Reverse complement if the strand was switched during setting the majority strand step
      #if (reverseComp) {
      gr.seq[which(paf.gr$revcomp == TRUE)] <- Biostrings::reverseComplement(gr.seq[which(paf.gr$revcomp == TRUE)])
      #}  
      #names(gr.seq) <- as.character(gr)
    } else {
      stop("Please set a 'bsgenome' or 'asm.fasta' parameter!!!")
    }
    
    ## Concatenate multiple sequence into a single FASTA
    if (length(paf.gr) > 1) {
      ## Concatenate all sequences into a single FASTA separated by 100 N's.
      delim <- paste(rep('N', 100), collapse = '')
      gr.seq.collapsed <- Biostrings::DNAStringSet(paste(gr.seq, collapse = delim))
      names(gr.seq.collapsed) <- paste(names(gr.seq), collapse = ';')
      gr.seq <- gr.seq.collapsed
    }  
    
    ## Write final FASTA
    if (is.character(fasta.save)) {
      ## Remove comment character from sequence names
      names(gr.seq) <- gsub(names(gr.seq), pattern = '#', replacement = '_')
      Biostrings::writeXStringSet(x = gr.seq, filepath = fasta.save, format = 'fasta')
    } else {
      warning("Please speficify 'fasta.save' if you want to export FASTA into a file!!!")
    }  
    return(gr.seq)
  } else {
    return(paf)
  }  
}
