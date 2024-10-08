#' Export FASTA sequences from a set of alignments reported in PAF formatted file.
#'
#' @param alignment.space What alignment coordinates should be exported as FASTA, either 'query' or 'target' (Default : `query`).
#' @param order.by Order alignment either by `query` or `target` coordinates.
#' @param bsgenome A \pkg{\link[BSgenome]{BSgenome-class}} object of reference genome to get the genomic sequence from.
#' @param asm.fasta An assembly FASTA file to extract DNA sequence from defined PAF alignments.
#' @param revcomp If set to \code{TRUE} FASTA sequence will be reverse complemented regardless of value defined in `majority.strand`.
#' @param report.longest.aln If set to \code{TRUE} only the sequence with the most aligned bases will be reported in final FASTA file.
#' @param report.query.name A single query (contig) name/id to be reported as FASTA sequence.
#' @param concatenate.aln Set to \code{TRUE} if multiple aligned contigs should be concatenated by 100 N's in to a single FASTA sequence (Default : `TRUE`).
#' @param fasta.save A path to a filename where to store final FASTA file.
#' @param return Set to either 'fasta' or 'index' to return either FASTA in \code{\link{DNAStringSet-class}} object or region index in \code{\link{GRanges-class}} object is returned.
#' @inheritParams breakPaf
#' @inheritParams syncRangesDir
#' @importFrom Rsamtools indexFa FaFile scanFa scanFaIndex
#' @importFrom BSgenome getSeq
#' @import GenomicRanges
#' @importFrom Biostrings writeXStringSet reverseComplement DNAStringSet
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges subsetByOverlaps
#' @return A \code{\link{DNAStringSet-class}} object with exported sequence
#' @author David Porubsky
#' @export
#' @examples
#' ## Get PAF to process ##
#' paf.file <- system.file("extdata", "test4.paf", package = "SVbyEye")
#' ## Read in PAF
#' paf.table <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = "cg")
#' ## Get FASTA using query alignment coordinates ##
#' ## Define assembly FASTA to get the sequence from
#' asm.fasta <- system.file("extdata", "test4_query.fasta", package = "SVbyEye")
#' paf2FASTA(paf.table = paf.table, alignment.space = "query", asm.fasta = asm.fasta)
#' \donttest{
#' ## Get FASTA using target alignment coordinates ##
#' ## Define BSgenome object to get the sequence from
#' paf2FASTA(
#'     paf.table = paf.table, alignment.space = "target",
#'     bsgenome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
#' )
#' }
paf2FASTA <- function(paf.table, alignment.space = "query", order.by = "query", bsgenome = NULL, asm.fasta = NULL, majority.strand = NULL, revcomp = NULL, report.longest.aln = FALSE, report.query.name = NULL, concatenate.aln = TRUE, fasta.save = NULL, return = "fasta") {
    ptm <- startTimedMessage("[paf2FASTA] Exporting PAF alignments to FASTA file")
    ## Check user input ##
    ## Make sure submitted paf.table has at least 12 mandatory fields
    if (ncol(paf.table) >= 12) {
        paf <- paf.table
        paf$direction.flip <- FALSE
    } else {
        stop("Submitted PAF alignments do not contain a minimum of 12 mandatory fields, see PAF file format definition !!!")
    }
    ## Load BSgenome object
    if (!is(bsgenome, "BSgenome")) {
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

    if (!is.null(paf)) {
        ## Order alignments
        if (order.by == "target") {
            paf <- paf[order(paf$t.start), ]
        } else if (order.by == "query") {
            paf <- paf[order(paf$q.start), ]
        }

        ## Convert query or target coordinates to GRanges
        if (alignment.space == "query") {
            paf.gr <- GenomicRanges::makeGRangesFromDataFrame(paf, seqnames.field = "q.name", start.field = "q.start", end.field = "q.end")
        } else {
            paf.gr <- GenomicRanges::makeGRangesFromDataFrame(paf, seqnames.field = "t.name", start.field = "t.start", end.field = "t.end")
        }
        # target.gr <- GenomicRanges::makeGRangesFromDataFrame(paf, seqnames.field = 't.name', start.field = 't.start', end.field = 't.end', strand='*')
        # paf.gr$target.gr <- target.gr

        ## Make sure no genomic region starts with zero
        GenomicRanges::start(paf.gr) <- pmax(GenomicRanges::start(paf.gr), 1)

        ## If defined process only a certain query name(s)
        if (!is.null(report.query.name)) {
            if (report.query.name %in% as.character(GenomeInfoDb::seqnames(paf.gr))) {
                paf.gr <- GenomeInfoDb::keepSeqlevels(paf.gr, value = report.query.name, pruning.mode = "coarse")
                ## Make sure parameter 'report.longest.aln' is set to FALSE
                report.longest.aln <- FALSE
            } else {
                warning("User defined 'report.query.name', ", report.query.name, " doesn't exist in submitted paf file!!!")
            }
        }

        ## Sync alignment directionality based on preferred majority strand
        split.f <- factor(as.character(GenomeInfoDb::seqnames(paf.gr)), levels = unique(as.character(GenomeInfoDb::seqnames(paf.gr))))
        paf.grl <- GenomicRanges::split(paf.gr, split.f)
        for (i in seq_along(paf.grl)) {
            gr <- paf.grl[[i]]
            qy.red <- base::range(gr, ignore.strand = TRUE)
            # tg.red <- range(gr$target.gr, ignore.strand=TRUE)
            if (!is.null(majority.strand)) {
                if (majority.strand %in% c("+", "-")) {
                    new.strand <- syncRangesDir(ranges = gr, majority.strand = majority.strand, strand.only = TRUE)
                    if (all(new.strand != GenomicRanges::strand(gr))) {
                        GenomicRanges::strand(gr) <- new.strand
                        gr <- qy.red
                        # gr$target.gr <- tg.red
                        gr$revcomp <- TRUE
                    } else {
                        gr <- qy.red
                        # gr$target.gr <- tg.red
                        gr$revcomp <- FALSE
                    }
                } else {
                    gr <- qy.red
                    # gr$target.gr <- tg.red
                    gr$revcomp <- FALSE
                    warning("Parameter 'majority.strand' can only takes values '+' or '-'!!!")
                }
            } else {
                gr <- qy.red
                gr$revcomp <- FALSE
            }
            paf.grl[[i]] <- gr
        }
        paf.gr <- unlist(paf.grl, use.names = FALSE)
        ## Force reverse complement if 'revcomp' parameter is set to TRUE/FALSE
        if (!is.null(revcomp)) {
            paf.gr$revcomp <- revcomp
        }

        ## Order regions by position
        # paf.gr <- GenomicRanges::sort(paf.gr, ignore.strand = TRUE)
        ## Order regions by target position [Prefered see above]
        # paf.gr <- paf.gr[GenomicRanges::order(paf.gr$target.gr)]

        ## Collapse consecutive alignments coming from the same contig/sequence
        # paf.gr$q.id <- as.character(GenomeInfoDb::seqnames(paf.gr))
        # paf.gr <- primatR::collapseBins(paf.gr, id.field = 3)

        ## Extract FASTA sequence
        if (!is.null(bsgenome)) {
            ## Make sure all sequence ranges are present in defined BSgenome object
            if (!all(as.character(GenomeInfoDb::seqnames(paf.gr)) %in% as.character(GenomeInfoDb::seqnames(bsgenome)))) {
                warning("Not all PAF ranges are present in the submitted FASTA file, subsetting !!!")
                paf.gr <- suppressWarnings(IRanges::subsetByOverlaps(paf.gr, as(seqinfo(bsgenome), "GRanges")))
                if (length(paf.gr) == 0) {
                    stop("None of the PAF ranges present in the submitted FASTA file, likely wrong FASTA file submitted!!!")
                }
            }
            ## Extract FASTA from BSgenome object
            gr.seq <- BSgenome::getSeq(bsgenome, paf.gr)
            names(gr.seq) <- as.character(paf.gr)
        } else if (is.character(asm.fasta)) {
            ## Extract FASTA from user defined FASTA file
            fa.file <- open(Rsamtools::FaFile(asm.fasta))
            ## Remove sequences not present in the FASTA index
            fa.idx <- Rsamtools::scanFaIndex(fa.file)
            ## Make sure all sequence ranges are present in submitted FASTA sequences
            if (!all(as.character(GenomeInfoDb::seqnames(paf.gr)) %in% as.character(GenomeInfoDb::seqnames(fa.idx)))) {
                warning("Not all PAF ranges are present in the submitted FASTA file, subsetting !!!")
                paf.gr <- suppressWarnings(IRanges::subsetByOverlaps(paf.gr, fa.idx))
                if (length(paf.gr) == 0) {
                    stop("None of the PAF ranges present in the submitted FASTA file, likely wrong FASTA file submitted!!!")
                }
            }
            ## Read in sequence for a given range(s)
            gr.seq <- Rsamtools::scanFa(file = fa.file, param = paf.gr, as = "DNAStringSet")
            ## Reverse complement if the strand was switched during setting the majority strand step
            # if (reverseComp) {
            gr.seq[which(paf.gr$revcomp == TRUE)] <- Biostrings::reverseComplement(gr.seq[which(paf.gr$revcomp == TRUE)])
            # }
            # names(gr.seq) <- as.character(gr)
        } else {
            stop("Please set a 'bsgenome' or 'asm.fasta' parameter!!!")
        }

        ## If TRUE report only the contig with the longest alignment
        if (report.longest.aln) {
            gr.seq <- gr.seq[which.max(width(gr.seq))]
            index.gr <- paf.gr[which.max(width(paf.gr))]
        } else if (concatenate.aln) {
            ## Concatenate multiple sequences into a single FASTA
            if (length(paf.gr) > 1) {
                ## Concatenate all sequences into a single FASTA separated by 100 N's.
                delim <- paste(rep("N", 100), collapse = "")
                gr.seq.collapsed <- Biostrings::DNAStringSet(paste(gr.seq, collapse = delim))
                names(gr.seq.collapsed) <- paste(names(gr.seq), collapse = ";")
                gr.seq <- gr.seq.collapsed
            }
            index.gr <- paf.gr
        } else {
            index.gr <- paf.gr
            gr.seq <- gr.seq
        }

        ## Write final FASTA
        if (is.character(fasta.save)) {
            ## Remove comment character from sequence names
            names(gr.seq) <- gsub(names(gr.seq), pattern = "#", replacement = "_")
            Biostrings::writeXStringSet(x = gr.seq, filepath = fasta.save, format = "fasta")
        } # else {
        # warning("Please speficify 'fasta.save' if you want to export FASTA into a file!!!")
        # }
        ## Return extracted FASTA or index
        if (return == "fasta") {
            stopTimedMessage(ptm)
            return(gr.seq)
        } else if (return == "index") {
            stopTimedMessage(ptm)
            return(index.gr)
        } else {
            stopTimedMessage(ptm)
            return(NULL)
        }
    } else {
        stopTimedMessage(ptm)
        return(paf)
    }
}
