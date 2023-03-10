#' Export FASTA sequences from a set of Genomic Ranges.
#'
#' This function takes a \code{\link{GRanges-class}} object and extracts a genomic sequence
#' from these regions either from an original range or a range expanded on each side by
#' defined number of bases.
#'
#' @param gr A \code{\link{GRanges-class}} object of genomic regions to extract genomic sequence from.
#' @param index.field A user defined column number to be used as a sequence ID for the exported FASTA.
#' @param expand Expand edges of each genomic region by this length (in bp).
#' @inheritParams paf2FASTA
#' @importFrom Rsamtools indexFa FaFile scanFa scanFaIndex
#' @importFrom BSgenome getSeq
#' @importFrom methods is
#' @importFrom Biostrings writeXStringSet
#' @importFrom GenomeInfoDb seqlengths seqnames seqlevels
#' @importFrom GenomicRanges start end mcols
#' @importFrom IRanges subsetByOverlaps
#' @author David Porubsky
#' @export
#' @examples
#' ## Get PAF to process ##
#' paf.file <- system.file("extdata", "test_getFASTA.paf", package = "SVbyEye")
#' ## Read in PAF
#' paf.table <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = "cg")
#' ## Get FASTA sequence from a defined genomic region ##
#' query.fasta <- system.file("extdata", "test_getFASTA_query.fasta", package = "SVbyEye")
#' ## Define genomic ranges to export FASTA sequence from
#' query.gr <- GenomicRanges::GRanges(
#'     seqnames = unique(paf.table$q.name),
#'     ranges = IRanges::IRanges(
#'         start = paf.table$q.start,
#'         end = paf.table$q.end
#'     )
#' )
#' ## Extract FASTA
#' regions2FASTA(gr = query.gr, asm.fasta = query.fasta)
#' ## Get FASTA sequence around breakpoints of alignment indels ##
#' ## Break PAF alignment at indels of 10 kbp and longer
#' paf.table <- breakPaf(paf.table = paf.table, min.deletion.size = 10000, min.insertion.size = 10000)
#' paf.svs <- paf.table$SVs
#' ## Define breakpoint positions of detected indels
#' sv.bps.gr <- GenomicRanges::GRanges(
#'     seqnames = unique(paf.svs$q.name),
#'     ranges = IRanges::IRanges(
#'         start = paf.svs$q.start,
#'         end = paf.svs$q.start
#'     ),
#'     bp.index = paste0(paf.svs$q.name, "-left_bp-", 1:nrow(paf.svs))
#' )
#' ## Extract FASTA
#' regions2FASTA(gr = sv.bps.gr, asm.fasta = query.fasta, expand = 10, index = 1)
#'
regions2FASTA <- function(gr, bsgenome = NULL, asm.fasta = NULL, index.field = NULL, expand = 0, fasta.save = NULL) {
   ptm <- startTimedMessage("[regions2FASTA] Exporting GRanges object to FASTA")

   ## Load BSgenome object
   if (!methods::is(bsgenome, "BSgenome")) {
       if (is.character(bsgenome)) {
          warning("Parameter 'bsgenome' has to be a valid bsgenome class object !!!")
       } else {
          bsgenome <- NULL
       }
    }
    ## Check if submitted fasta file is indexed
    if (!is.null(asm.fasta)) {
        asm.fasta.idx <- paste0(asm.fasta, ".fai")
        if (!file.exists(asm.fasta.idx)) {
            fa.idx <- Rsamtools::indexFa(file = asm.fasta)
        }
    }
    ## Check if user defined regions are present in either in bsgenome object or submitted asm.fasta file
    if (!is.null(bsgenome)) {
        ## Make sure all sequence ranges are present in defined BSgenome object
        if (!all(as.character(GenomeInfoDb::seqnames(gr)) %in% as.character(GenomeInfoDb::seqnames(bsgenome)))) {
            warning("Not all PAF ranges are present in the submitted FASTA file, subsetting !!!")
            gr <- suppressWarnings(IRanges::subsetByOverlaps(gr, as(seqinfo(bsgenome), "GRanges")))
            if (length(gr) == 0) {
                stop("None of the PAF ranges present in the submitted FASTA file, likely wrong FASTA file submitted!!!")
            }
        }
        ## Add sequence lengths if missing
        if (any(is.na(GenomeInfoDb::seqlengths(gr)))) {
            suppressWarnings(GenomeInfoDb::seqlengths(gr) <- GenomeInfoDb::seqlengths(bsgenome)[GenomeInfoDb::seqlevels(gr)])
        }
    } else if (is.character(asm.fasta)) {
        ## Extract FASTA from user defined FASTA file
        fa.file <- open(Rsamtools::FaFile(asm.fasta))
        ## Remove sequences not present in the FASTA index
        fa.idx <- Rsamtools::scanFaIndex(fa.file)
        ## Make sure all sequence ranges are present in submitted FASTA sequences
        if (!all(as.character(GenomeInfoDb::seqnames(gr)) %in% as.character(GenomeInfoDb::seqnames(fa.idx)))) {
            warning("Not all PAF ranges are present in the submitted FASTA file, subsetting !!!")
            gr <- suppressWarnings(IRanges::subsetByOverlaps(gr, fa.idx))
            if (length(gr) == 0) {
                stop("None of the PAF ranges present in the submitted FASTA file, likely wrong FASTA file submitted!!!")
            }
        }
        ## Add sequence lengths if missing
        if (any(is.na(GenomeInfoDb::seqlengths(gr)))) {
            suppressWarnings(GenomeInfoDb::seqlengths(gr) <- GenomeInfoDb::seqlengths(fa.idx)[GenomeInfoDb::seqlevels(gr)])
        }
    }
    ## Make sure no genomic region starts with zero
    GenomicRanges::start(gr) <- pmax(GenomicRanges::start(gr), 1)

    ## Expand regions of interest by certain size (downstream and upstream)
    if (expand > 0) {
        gr.seqLen <- GenomeInfoDb::seqlengths(gr)[as.character(seqnames(gr))]
        if (all(!is.na(gr.seqLen))) {
            start(gr) <- pmax(1, GenomicRanges::start(gr) - expand)
            end(gr) <- pmin(gr.seqLen, GenomicRanges::end(gr) + expand)
        } else {
            warning("Could not expand the exported regions, missing 'seqlengths' in submitted 'gr' object!")
        }
    }

    ## Extract FASTA sequence
    if (!is.null(bsgenome)) {
        ## Extract FASTA from BSgenome object
        gr.seq <- BSgenome::getSeq(bsgenome, gr)
        names(gr.seq) <- as.character(gr)
    } else if (is.character(asm.fasta)) {
        ## Read in sequence for a given range(s)
        gr.seq <- Rsamtools::scanFa(file = fa.file, param = gr, as = "DNAStringSet")
    } else {
        stop("Please set a 'bsgenome' or 'asm.fasta' parameter!!!")
    }

    ## Used user defined column to name FASTA sequences
    if (!is.null(index.field)) {
        if (index.field <= length(GenomicRanges::mcols(gr))) {
            names(gr.seq) <- as.character(GenomicRanges::mcols(gr)[[index.field]])
        }
    }

    ## Write final FASTA
    if (is.character(fasta.save)) {
        ## Remove comment character from sequence names
        names(gr.seq) <- gsub(names(gr.seq), pattern = "#", replacement = "_")
        Biostrings::writeXStringSet(x = gr.seq, filepath = fasta.save, format = "fasta")
    } # else {
    # warning("Please speficify 'fasta.save' if you want to export FASTA into a file!!!")
    # }
    stopTimedMessage(ptm)
    return(gr.seq)
}
