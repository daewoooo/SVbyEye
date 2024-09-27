#' Flip orientation of PAF alignments.
#'
#' This function takes loaded PAF alignments using \code{\link{readPaf}} function and flips
#' the orientation of PAF alignments given the desired 'majority.strand' orientation (Either '+' or '-').
#'
#' @param force Set to \code{TRUE} if query PAF alignments should be flipped.
#' @param flip.seqnames User defined query or target sequence IDs to be flipped in orientation.
#' @inheritParams breakPaf
#' @inheritParams syncRangesDir
#' @return A \code{tibble} of PAF alignments flipped based on desired majority strand orientation.
#' @importFrom tibble is_tibble as_tibble
#' @author David Porubsky
#' @export
#' @examples
#' ## Get PAF to process
#' paf.file <- system.file("extdata", "test1.paf", package = "SVbyEye")
#' ## Read in PAF
#' paf.table <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = "cg")
#' ## Force to flip PAF alignments
#' flipPaf(paf.table = paf.table, force = TRUE)
#'
flipPaf <- function(paf.table, majority.strand = NULL, force = FALSE, flip.seqnames = NULL) {
    ## Helper function
    flipPafAln <- function(paf.table = NULL, coordinates = "query") {
        paf.flip <- paf.table
        ## Flip alignment orientation
        paf.flip$strand[paf.table$strand == "+"] <- "-"
        paf.flip$strand[paf.table$strand == "-"] <- "+"
        ## Keep info on original query coordinates [To test]
        # paf.flip$q.end <- paf.table$q.start
        # paf.flip$q.start <- paf.table$q.end
        if (coordinates == "query") {
            ## Flip query coordinates
            paf.flip$q.end <- paf.table$q.len - paf.table$q.start
            paf.flip$q.start <- paf.table$q.len - paf.table$q.end
            ## Add info if the PAF alignment was flipped
            paf.flip$query.flip <- TRUE
        } else if (coordinates == "target") {
            ## Flip target coordinates
            paf.flip$t.end <- paf.table$q.len - paf.table$t.start
            paf.flip$t.start <- paf.table$q.len - paf.table$t.end
            ## Add info if the PAF alignment was flipped
            paf.flip$target.flip <- TRUE
        } else {
            warning("Parameter 'coordinates' can only take values 'query' or 'target' !!!")
            paf.flip <- paf.table
        }
        ## Make sure q.start position is always smaller than q.end
        # q.start <- pmin(paf.flip$q.start, paf.flip$q.end)
        # q.end <- pmax(paf.flip$q.start, paf.flip$q.end)
        # paf.flip$q.start <- q.start
        # paf.flip$q.end <- q.end

        return(paf.flip)
    }

    ## Check user input ##
    ## Make sure submitted paf.table has at least 12 mandatory fields
    if (ncol(paf.table) >= 12) {
        paf <- paf.table
        paf$query.flip <- FALSE
    } else {
        stop("Submitted PAF alignments do not contain a minimum of 12 mandatory fields, see PAF file format definition !!!")
    }

    ptm <- startTimedMessage("[flipPaf] Flipping orientation of PAF alignments")

    ## Get unique alignment ID
    paf$seq.pair <- paste0(paf$q.name, "__", paf$t.name)
    ## Add flip info
    paf$query.flip <- FALSE
    paf$target.flip <- FALSE

    ## Force alignment flip
    if (force) {
        paf.l <- split(paf, paf$seq.pair)
        for (i in seq_along(paf.l)) {
            paf.ctg <- paf.l[[i]]
            paf.l[[i]] <- flipPafAln(paf.table = paf.ctg)
        }
        paf <- do.call(rbind, paf.l)
    }

    ## Flip specific sequence names either query or target
    if (!is.null(flip.seqnames)) {
        ## Flip query specific sequences
        if (any(flip.seqnames %in% paf$q.name)) {
            sub.paf <- paf[paf$q.name %in% flip.seqnames, ]
            sub.paf.l <- split(sub.paf, sub.paf$seq.pair)
            for (i in seq_along(sub.paf.l)) {
                paf.ctg <- sub.paf.l[[i]]
                sub.paf.l[[i]] <- flipPafAln(paf.table = paf.ctg, coordinates = "query")
            }
            sub.paf <- do.call(rbind, sub.paf.l)
            paf[paf$q.name %in% flip.seqnames, ] <- sub.paf
        }
        ## Flip target specific sequences
        if (any(flip.seqnames %in% paf$t.name)) {
            sub.paf <- paf[paf$t.name %in% flip.seqnames, ]
            sub.paf.l <- split(sub.paf, sub.paf$seq.pair)
            for (i in seq_along(sub.paf.l)) {
                paf.ctg <- sub.paf.l[[i]]
                sub.paf.l[[i]] <- flipPafAln(paf.table = paf.ctg, coordinates = "target")
            }
            sub.paf <- do.call(rbind, sub.paf.l)
            paf[paf$t.name %in% flip.seqnames, ] <- sub.paf
        }
    }

    ## Sync by majority strand directionality
    if (!is.null(majority.strand) & force == FALSE) {
        ## Define majority and minority strand
        if (majority.strand == "+") {
            minority.strand <- "-"
        } else if (majority.strand == "-") {
            minority.strand <- "+"
        } else {
            stop("Parameter 'majority.strand' can only be defined as '+' or '-' !!!")
        }

        ## Make sure alignment is defined before deciding on majority directionality
        if (!is.na(sum(paf$aln.len))) {
            if (sum(paf$aln.len) > 0) {
                paf.l <- split(paf, paf$seq.pair)
                for (i in seq_along(paf.l)) {
                    paf.ctg <- paf.l[[i]]
                    ## Flip directionality based to make sure majority strand covers the most bases
                    majority.bases <- sum(paf.ctg$aln.len[paf.ctg$strand == majority.strand])
                    minority.bases <- sum(paf.ctg$aln.len[paf.ctg$strand == minority.strand])
                    if (majority.bases < minority.bases) {
                        paf.l[[i]] <- flipPafAln(paf.table = paf.ctg)
                    } else {
                        paf.l[[i]] <- paf.ctg
                    }
                }
                paf <- do.call(rbind, paf.l)
            }
        }
    }

    stopTimedMessage(ptm)
    ## Export data object
    if (tibble::is_tibble(paf)) {
        return(paf)
    } else {
        return(tibble::as_tibble(paf))
    }
}
