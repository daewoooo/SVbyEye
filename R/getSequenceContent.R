#' Add FASTA sequence content to PAF alignments.
#'
#' This function takes a PAF table and for each alignment (rows) will report counts and frequencies of user defined
#' `sequence.pattern` (such as exact DNA pattern, e.g. 'GA') or `nucleotide.content` (such as sequence GC content).
#'
#' @inheritParams breakPaf
#' @inheritParams paf2FASTA
#' @inheritParams getNucleotideContent
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom dplyr bind_cols
#' @author David Porubsky
#' @export
#' @examples
#' ## Get PAF to process ##
#' paf.file <- system.file("extdata", "test_getFASTA.paf", package = "SVbyEye")
#' ## Read in PAF
#' paf.table <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = "cg")
#' ## Split PAF alignments into user defined bins
#' paf.table <- pafToBins(paf.table = paf.table, binsize = 10000)
#' ## Get FASTA to process
#' asm.fasta <- system.file("extdata", "test_getFASTA_query.fasta", package = "SVbyEye")
#' ## Add sequence and nucleotide content to submitted paf.table
#' paf2nucleotideContent(
#'     paf.table = paf.table, asm.fasta = asm.fasta,
#'     alignment.space = "query", sequence.pattern = "GA"
#' )
#'
paf2nucleotideContent <- function(paf.table = NULL, asm.fasta = NULL, alignment.space = NULL, sequence.pattern = NULL, nucleotide.content = NULL) {
    ## Check user input ##
    ## Make sure submitted paf.table has at least 12 mandatory fields
    if (ncol(paf.table) >= 12) {
        paf <- paf.table
    } else {
        stop("Submitted PAF alignments do not contain a minimum of 12 mandatory fields, see PAF file format definition !!!")
    }

    ## Get regions to extract FASTA sequence from
    if (alignment.space == "query") {
        paf.gr <- GenomicRanges::GRanges(
            seqnames = unique(paf.table$q.name),
            ranges = IRanges::IRanges(
                start = paf.table$q.start,
                end = paf.table$q.end
            ),
            index = paste(unique(paf.table$q.name), seq_len(nrow(paf.table)), sep = "_")
        )
    } else {
        paf.gr <- GenomicRanges::GRanges(
            seqnames = unique(paf.table$t.name),
            ranges = IRanges::IRanges(
                start = paf.table$t.start,
                end = paf.table$t.end
            ),
            index = paste(unique(paf.table$t.name), seq_len(nrow(paf.table)), sep = "_")
        )
    }
    ## Get FASTA sequences per PAF region
    fasta.seq <- regions2FASTA(gr = paf.gr, asm.fasta = asm.fasta, index.field = 1)
    ## Get nucleotide or pattern count and frequency
    nuc.content <- getNucleotideContent(fasta.seq = fasta.seq, sequence.pattern = sequence.pattern, nucleotide.content = nucleotide.content)
    ## Add defined FASTA nucleotide content to submitted paf.table
    paf.table <- dplyr::bind_cols(paf.table, nuc.content)
    return(paf.table)
}


#' Get sequence content from a single FASTA file.
#'
#' This function with takes a path to FASTA file containing a single sequence and report counts and frequencies of user defined
#' `sequence.pattern` (such as exact DNA pattern, e.g. 'GA') or `nucleotide.content` (such as sequence GC content) across the
#' whole sequence or in user defined bins.
#'
#' @param fasta.file A path to a valid FASTA file to be processed.
#' @param binsize A user defined number of base pairs to split the loaded FASTA sequence into
#' @inheritParams getNucleotideContent
#' @importFrom Biostrings readDNAStringSet Views
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges ranges
#' @author David Porubsky
#' @export
#' @examples
#' ## Get FASTA to process ##
#' asm.fasta <- system.file("extdata", "test_getFASTA_query.fasta", package = "SVbyEye")
#' ## Report sequence content for a given FASTA file in 5000bp long bins
#' fasta2nucleotideContent(fasta.file = asm.fasta, binsize = 5000, sequence.pattern = "AT")
#'
fasta2nucleotideContent <- function(fasta.file, binsize=NULL, sequence.pattern = NULL, nucleotide.content = NULL) {
  ptm <- startTimedMessage("[fasta2nucleotideContent] Calculating FASTA sequence content")

  ## Check user input
  if (file.exists(fasta.file)) {
    fa.seq <- Biostrings::readDNAStringSet(fasta.file)
    seq.len <- GenomeInfoDb::seqlengths(fa.seq)
  } else {
    stop('Please submit a path to a valid FASTA file !!!')
  }

  ## Bin FASTA sequence
  if (!is.null(binsize)) {
    if (binsize > 0) {
      views.obj <- Biostrings::Views(unlist(fa.seq), start = seq(from=1, to=(seq.len - binsize), by=binsize), width = binsize)
      views.gr <- GenomicRanges::GRanges(seqnames = names(fa.seq), ranges = IRanges::ranges(views.obj))
      seq.views <- as(views.obj, 'DNAStringSet')
    } else {
      views.gr <- GenomicRanges::GRanges(seqnames = names(fa.seq), ranges = IRanges::IRanges(start = 1, end = seq.len))
      seq.views <- fa.seq
    }
  } else {
    views.gr <- GenomicRanges::GRanges(seqnames = names(fa.seq), ranges = IRanges::IRanges(start = 1, end = seq.len))
    seq.views <- fa.seq
  }
  ## Get nucleotide or pattern count and frequency
  seq.content <- suppressMessages( getNucleotideContent(fasta.seq = seq.views,
                                                        sequence.pattern = sequence.pattern,
                                                        nucleotide.content = nucleotide.content) )
  mcols(views.gr) <- seq.content
  ## Return calcualted sequence content
  stopTimedMessage(ptm)
  return(views.gr)
}


#' Get sequence content from one or multiple FASTA sequences.
#'
#' This function takes a single or set of FASTA sequences and for each sequence will report counts and frequencies of user defined
#' `sequence.pattern` (such as exact DNA pattern, e.g. 'GA') or `nucleotide.content` (such as sequence GC content).
#'
#' @param fasta.seq A \code{\link{DNAStringSet-class}} object containing one or multiple FASTA sequences.
#' @param sequence.pattern A user defined DNA sequence pattern which occurrences will be counted in submitted `fasta.seq`.
#' (e.g. set `sequence.pattern` to 'GA' to obtain counts of all 'GA' occurrences per FASTA sequence).
#' @param nucleotide.content A user defined nucleotides which total content will be counted in submitted `fasta.seq`.
#' (e.g. to obtain 'GC' content set `nucleotide.content` to 'GC')
#' @importFrom Biostrings vcountPattern alphabetFrequency width
#' @importFrom dplyr bind_cols
#' @importFrom tibble tibble
#' @author David Porubsky
#' @export
#' @examples
#' ## Get PAF to process ##
#' paf.file <- system.file("extdata", "test_getFASTA.paf", package = "SVbyEye")
#' ## Read in PAF
#' paf.table <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = "cg")
#' ## Get FASTA to process
#' asm.fasta <- system.file("extdata", "test_getFASTA_query.fasta", package = "SVbyEye")
#' ## Read in FASTA sequence
#' fasta.seq <- paf2FASTA(paf.table = paf.table, asm.fasta = asm.fasta)
#' ## Report sequence and nucleotide content for a given FASTA sequence
#' getNucleotideContent(fasta.seq = fasta.seq, sequence.pattern = "AT", nucleotide.content = "GC")
#'
getNucleotideContent <- function(fasta.seq, sequence.pattern = NULL, nucleotide.content = NULL) {
    ptm <- startTimedMessage("[getNucleotideContent] Calculating FASTA sequence content")

    ## Check user input ##
    ## Check if submitted FASTA sequence is in accepted format
    if (!is(fasta.seq, "DNAStringSet")) {
        stop("Parameter 'fasta.seq' has to be 'DNAStringSet' class object !!!")
    }
    ## Check if submitted sequence.pattern nucleotide.content are in accepted format
    if (!is.null(sequence.pattern)) {
        if (is.character(sequence.pattern) & nchar(sequence.pattern) > 0) {
            if (all(grepl(strsplit(sequence.pattern, split = "")[[1]], pattern = "A|C|G|T", ignore.case = TRUE))) {

            } else {
                sequence.pattern <- NULL
                warning("User defined parameter 'sequence.pattern' is not a valid DNA sequence !!!")
            }
        }
    } else if (!is.null(nucleotide.content)) {
        if (is.character(nucleotide.content) & nchar(nucleotide.content) > 0) {
            if (all(grepl(strsplit(nucleotide.content, split = "")[[1]], pattern = "A|C|G|T", ignore.case = TRUE))) {

            } else {
                nucleotide.content <- NULL
                warning("User defined parameter 'nucleotide.content' is not a valid DNA sequence !!!")
            }
        }
    } else if (!is.null(sequence.pattern) & !is.null(nucleotide.content)) {
        stop("Parameters 'sequence.pattern' and/or 'nucleotide.content' are either not defined or contain invalid DNA sequence string !!!")
    }

    ## Get sequence content counts and frequencies ##
    fasta.content <- tibble::tibble(.rows = length(fasta.seq))
    if (!is.null(sequence.pattern)) {
        ## Get counts of a specific sequence pattern per FASTA sequence
        pattern.counts <- Biostrings::vcountPattern(pattern = sequence.pattern, subject = fasta.seq)
        ## Get pattern content per sequence length [Perhaps it will be better to find all locations]
        pattern.freq <- (pattern.counts * nchar(sequence.pattern)) / Biostrings::width(fasta.seq)

        # ## Get positions of a specific sequence pattern [needed for binned counts]
        # pattern.pos <- Biostrings::vmatchPattern(pattern = sequence.pattern, subject = fasta.seq)
        # ## Convert to GRanges object
        # pattern.pos.gr <- as(pattern.pos, 'GRanges')
        # seq.bins <- tileGenome(seqlengths = dna.len, tilewidth = bin.width, cut.last.tile.in.chrom = TRUE)
        # ## Count number of di-nucleotides in each bin
        # binned.counts <- countOverlaps(seq.bins, pattern.pos.gr)

        ## Report counts and frequencies
        seq.content <- suppressMessages(dplyr::bind_cols(as.numeric(pattern.counts), pattern.freq))
        colnames(seq.content) <- c(paste0(sequence.pattern, "_seq.count"), paste0(sequence.pattern, "_seq.freq"))
        fasta.content <- dplyr::bind_cols(fasta.content, seq.content)
    }

    if (!is.null(nucleotide.content)) {
        ## Get counts of a specific nucleotide content
        base.freq <- Biostrings::alphabetFrequency(fasta.seq)
        bases <- toupper(strsplit(nucleotide.content, split = "")[[1]])
        ## Get counts of a specific nucleotides per FASTA sequence
        bases.counts <- rowSums(base.freq[, bases, drop = FALSE])
        ## Get frequency of a specific nucleotides per FASTA sequence
        bases.freq <- bases.counts / rowSums(base.freq[, c('A', 'C', 'G', 'T'), drop = FALSE])

        ## Report counts and frequencies
        nuc.content <- suppressMessages(dplyr::bind_cols(as.numeric(bases.counts), bases.freq))
        colnames(nuc.content) <- c(paste0(nucleotide.content, "_nuc.count"), paste0(nucleotide.content, "_nuc.freq"))
        # if (!is.null(sequence.pattern)) {
        #  fasta.content <- dplyr::right_join(x = fasta.content, y = nuc.content, by = c('ID' = 'ID'))
        # } else {
        fasta.content <- dplyr::bind_cols(fasta.content, nuc.content)
        # }
    }

    stopTimedMessage(ptm)
    ## Return sequence content table
    return(fasta.content)
}
