#' Mask FASTA sequence at defined regions.
#'
#' @param fasta.file A path to a FASTA file to be masked.
#' @param mask.ranges A \code{\link{IRanges-class}} object of coordinates to be masked in input FASTA.
#' @param invert If set \code{TRUE} ranges defined in 'mask.ranges' will be kept while the rest of the FASTA will be masked.
#' @param mask.character A single character to be used for masking (Default : 'N').
#' @inheritParams paf2FASTA
#' @importFrom Rsamtools scanFa
#' @importFrom Biostrings replaceAt width DNA_ALPHABET
#' @importFrom GenomicRanges setdiff
#' @importFrom IRanges IRanges ranges reduce width
#' @author David Porubsky
#' @export
#' @examples
#' ## FASTA sequence to mask
#' fasta.file <- system.file("extdata", "test4_query.fasta", package = "SVbyEye")
#' ## Define ranges to mask
#' mask.ranges <- IRanges::IRanges(start = c(11, 1981), end = c(20, 1990))
#' ## Mask FASTA sequence at defined ranges
#' maskFASTA(fasta.file = fasta.file, mask.ranges = mask.ranges)
#' ## Mask FASTA sequence except of defined ranges
#' maskFASTA(fasta.file = fasta.file, mask.ranges = mask.ranges, invert = TRUE)
#'
maskFASTA <- function(fasta.file, mask.ranges = NULL, invert = FALSE, mask.character = "N", fasta.save = NULL) {
  ptm <- startTimedMessage(paste0("[maskFASTA] Masking FASTA file, ", fasta.file))

    if (file.exists(fasta.file)) {
        ## Read in FASTA sequence
        seq <- Rsamtools::scanFa(file = fasta.file, as = "DNAStringSet")
    } else {
        stop("Input FASTA: '", fasta.file, "' does not exists !!!")
    }

    ## Make sure user defined mask.character is in expected DNA_ALPHABET
    if (!mask.character %in% Biostrings::DNA_ALPHABET) {
        mask.character <- "N"
        warning("Parameter 'mask.character' can only take values defined in Biostrings::DNA_ALPHABET, using default value !!!")
    }

    if (!is.null(mask.ranges)) {
        ## Make sure input ranges are IRanges and nonoverlapping
        mask.ranges <- IRanges::reduce(IRanges::ranges(mask.ranges))
        ## Get ranges to mask
        seq.range <- IRanges::IRanges(start = 1, end = Biostrings::width(seq))
        if (invert) {
            mask.at.ranges <- GenomicRanges::setdiff(seq.range, mask.ranges)
        } else {
            mask.at.ranges <- mask.ranges
        }
        ## Create mask
        mask.n <- vapply(IRanges::width(mask.at.ranges), function(x) paste(rep(mask.character, times = x), collapse = ""), FUN.VALUE = character(1))
        ## Mask FASTA
        masked.seq <- Biostrings::replaceAt(seq, at = mask.at.ranges, value = mask.n)
    } else {
        masked.seq <- seq
    }

    ## Write final FASTA
    if (is.character(fasta.save)) {
        ## Remove comment character from sequence names
        names(masked.seq) <- gsub(names(masked.seq), pattern = "#", replacement = "_")
        Biostrings::writeXStringSet(x = masked.seq, filepath = fasta.save, format = "fasta")
    } else {
        message("Please speficify 'fasta.save' if you want to export FASTA into a file!!!")
    }

  stopTimedMessage(ptm)
    ## Return masked FASTA
    return(masked.seq)
}
