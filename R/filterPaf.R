#' Filter PAF alignments.
#'
#' This function takes loaded PAF alignments using \code{\link{readPaf}} function and perform
#' user defined filtering of input alignments based on mapping quality, alignment length, and
#' minimum alignments between target and query.
#'
#' @param min.mapq Minimum mapping quality to retain.
#' @param min.align.len Minimum alignment length to retain.
#' @param min.align.n Minimum number of alignments between a pair of sequences/regions.
#' @param is.selfaln Set to \code{TRUE} if processing alignments of a given sequence to itself.
#' @param min.selfaln.dist Keep alignment pairs with this or larger distance from each other (Applied only to FASTA self-alignments).
#' @inheritParams breakPaf
#' @return A \code{tibble} of filtered PAF alignments.
#' @importFrom dplyr group_by mutate n
#' @importFrom utils read.table
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @importFrom S4Vectors queryHits
#' @importFrom tibble is_tibble as_tibble
#' @author David Porubsky
#' @export
#' @examples
#' ## Get PAF to plot
#' paf.file <- system.file("extdata", "test1.paf", package = "SVbyEye")
#' ## Read in PAF
#' paf.table <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = "cg")
#' ## Filter PAF alignments based on desired target region
#' filterPaf(paf.table = paf.table, min.align.len = 50000)
#'
filterPaf <- function(paf.table, min.mapq = 10, min.align.len = 1000, min.align.n = 1, is.selfaln = FALSE, min.selfaln.dist = 0) {
    ## Check user input ##
    ## Make sure PAF has at least 12 mandatory fields
    if (ncol(paf.table) >= 12) {
        paf <- paf.table
    } else {
        stop("Submitted PAF alignments do not contain a minimum of 12 mandatory fields, see PAF file format definition !!!")
    }

    ptm <- startTimedMessage("[filterPaf] Filtering PAF alignments")

    # ## Filer alignments by target region
    # if (!is.null(target.region)) {
    #   if (is.character(target.region)) {
    #     target.region.gr <- as(target.region, 'GRanges')
    #   } else if (class(target.region.gr) == 'GRanges') {
    #     target.region.gr <- target.region
    #   } else {
    #     message("Parameter 'target.region' can either be 'GRanges' object or character string 'chr#:start-end'!!!")
    #   }
    #   ## Subset PAF by ranges overlaps
    #   target.gr <- GenomicRanges::makeGRangesFromDataFrame(paf, seqnames.field = 't.name', start.field = 't.start', end.field = 't.end')
    #   hits <- GenomicRanges::findOverlaps(target.gr, target.region.gr) ## TODO Take only overlaps, regions itself are not shrinked!!!
    #   if (length(hits) > 0) {
    #     paf <- paf[S4Vectors::queryHits(hits),]
    #     ## Narrow down target aligments to desired target region
    #     #paf$t.start[paf$t.start < GenomicRanges::start(target.region.gr)] <- GenomicRanges::start(target.region.gr)
    #     #paf$t.end[paf$t.end > GenomicRanges::end(target.region.gr)] <- GenomicRanges::end(target.region.gr)
    #
    #     #paf <- paf[paf$t.name %in% as.character(seqnames(target.region.gr)) & paf$t.start >= start(target.region.gr) & paf$t.end <= end(target.region.gr),]
    #   } else {
    #     warning("None of the PAF ranges overlap user defined 'target.region', skipping ...")
    #   }
    # }
    #
    # ## Filer alignments by query region
    # if (!is.null(query.region)) {
    #   if (is.character(query.region)) {
    #     query.region.gr <- as(query.region, 'GRanges')
    #   } else if (class(query.region.gr) == 'GRanges') {
    #     query.region.gr <- query.region
    #   } else {
    #     message("Parameter 'query.region' can either be 'GRanges' object or character string 'chr#:start-end'!!!")
    #   }
    #   ## Subset PAF by ranges overlaps
    #   query.gr <- GenomicRanges::makeGRangesFromDataFrame(paf, seqnames.field = 'q.name', start.field = 'q.start', end.field = 'q.end')
    #   hits <- GenomicRanges::findOverlaps(query.gr, query.region.gr) ## TODO Take only overlaps, regions itself are not shrinked!!!
    #   if (length(hits) > 0) {
    #     paf <- paf[S4Vectors::queryHits(hits),]
    #     #paf <- paf[paf$t.name %in% as.character(seqnames(target.region.gr)) & paf$t.start >= start(target.region.gr) & paf$t.end <= end(target.region.gr),]
    #   } else {
    #     warning("None of the PAF ranges overlap user defined 'query.region', skipping ...")
    #   }
    # }

    ## Filter by the number of alignments per sequence alignment pair
    ## Get unique alignment ID
    if (min.align.n > 0) {
        paf$seq.pair <- paste0(paf$q.name, "__", paf$t.name)
        paf <- paf %>%
            dplyr::group_by(.data$seq.pair) %>%
            dplyr::mutate(align.n = dplyr::n())
        paf <- paf[paf$align.n >= min.align.n,]
        paf <- paf[,-c('seq.pair', 'align.n')]
    }

    ## Keep only specific sequence/region name
    # @param seqname.grep Retain only a specific sequence/region name with a given character string.
    # if (!is.null(seqnames.grep)) {
    #   paf <- paf[grep(paf$q.name, pattern = seqnames.grep),]
    #   paf <- paf[grep(paf$t.name, pattern = seqnames.grep),]
    # }

    ## Filter by mapping quality
    ## Make sure mapping quality is defined
    if (all(!is.na(paf$mapq))) {
        paf$mapq[is.na(paf$mapq)] <- min.mapq
        if (min.mapq > 0 & is.numeric(paf$mapq)) {
            paf <- paf[paf$mapq >= min.mapq, ]
        }
    }

    ## Filter by alignment length
    if (all(!is.na(paf$aln.len))) {
        paf$aln.len[is.na(paf$aln.len)] <- min.align.len
        if (min.align.len > 0 & is.numeric(paf$aln.len)) {
            paf <- paf[paf$aln.len >= min.align.len, ]
        }
    }

    ## When processing self-alignments apply separate set of filters
    if (is.selfaln) {
      ## Remove diagonals (self-alignments) with the same start and end position
      paf <- paf[!(paf$q.start == paf$t.start & paf$q.end == paf$t.end),]
      ## Due to the minimap2 self-alignment redundancy keep only alignments where query start is smaller than the target start
      paf <- paf[paf$q.start < paf$t.start,]
      ## Filter by alignment distance [self-alignments only]
      aln.dist <- pmin(paf$t.start, paf$t.end) - pmax(paf$q.start, paf$q.end)
      if (min.selfaln.dist > 0) {
        paf <- paf[aln.dist >= min.selfaln.dist,]
      }
    } else {
      ## Filter out self-alignments
      paf <- paf[!(paf$q.name == paf$t.name), ]
    }

    stopTimedMessage(ptm)
    ## Export data object
    if (tibble::is_tibble(paf)) {
        return(paf)
    } else {
        return(tibble::as_tibble(paf))
    }
}
