#' Prepare PAF alignments for plotting.
#'
#' This function takes loaded PAF alignments using \code{\link{readPaf}} function. Such alignment could be post-processed
#' using \code{\link{filterPaf}}, \code{\link{breakPaf}} and \code{\link{flipPaf}} functions. Subsequently such alignments are
#' expanded in a set of x and y coordinates ready to be plotted by \code{\link{geom_miropeats}} function.
#'
#' @param offset.alignments Set to \code{TRUE} if subsequent target alignments should be offsetted below and above the midline.
#' @param add.col A user defined column name present in `paf.table` to be added in returned coordinates table.
#' @param sync.x.coordinates If set to \code{TRUE} query coordinates will be adjusted to the limits of target coordinates. (Default : `TRUE`)
#' @inheritParams breakPaf
#' @inheritParams q2t
#' @return A \code{tibble} of PAF alignments reported as x and y coordinate values.
#' @importFrom tibble add_column
#' @importFrom dplyr bind_cols
#' @importFrom data.table ':='
#' @author David Porubsky
#' @export
#' @examples
#' ## Get PAF to process ##
#' paf.file <- system.file("extdata", "test1.paf", package = "SVbyEye")
#' ## Read in PAF
#' paf.table <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = "cg")
#' ## Convert PAF alignments to x and y coordinates needed for plotting
#' paf2coords(paf.table = paf.table)
#' ## Include optional an column to be exported as well (e.g. alignment specific GC content)
#' paf.table$GC.content <- round(runif(nrow(paf.table), min = 30, max = 60), digits = 2)
#' paf2coords(paf.table = paf.table, add.col = "GC.content")
#'
paf2coords <- function(paf.table, offset.alignments = FALSE, sync.x.coordinates = TRUE, q.range = NULL, t.range = NULL, add.col = NULL) {
    ## Check user input ##
    ## Make sure PAF has at least 12 mandatory fields
    if (ncol(paf.table) >= 12) {
        paf <- paf.table
    } else {
        stop("Submitted PAF alignments do not contain a minimum of 12 mandatory fields, see PAF file format definition !!!")
    }
    ## Set required fields with default values if not defined
    if (!"aln.id" %in% colnames(paf)) {
        paf$aln.id <- seq_len(nrow(paf))
    }
    if (!"ID" %in% colnames(paf)) {
        paf$ID <- "M"
    }

    ## Flip start-end if strand == '-'
    paf[paf$strand == "-", c("t.start", "t.end")] <- rev(paf[paf$strand == "-", c("t.start", "t.end")])
    # paf[paf$strand == '-', c('q.start','q.end')] <- rev(paf[paf$strand == '-', c('q.start','q.end')])

    ## Add continuous scale to PAF to make sure multiple query and target sequences are position next to each other
    paf <- paf2continuousScale(paf.table = paf)
    paf$q.start.shift <- paf$q.start + paf$q.shift
    paf$q.end.shift <- paf$q.end + paf$q.shift
    paf$t.start.shift <- paf$t.start + paf$t.shift
    paf$t.end.shift <- paf$t.end + paf$t.shift

    ## Get unique alignment ID
    # if (!"seq.pair" %in% colnames(paf)) {
    #     paf$seq.pair <- paste0(paf$q.name, "__", paf$t.name)
    # }

    ## Sync scales between alignments [per region id]
    if (sync.x.coordinates) {
      # paf.l <- split(paf, paf$seq.pair)
      # for (i in seq_along(paf.l)) {
      #     paf.sub <- paf.l[[i]]
      #     q.range <- range(c(paf.sub$q.start, paf.sub$q.end))
      #     t.range <- range(c(paf.sub$t.start, paf.sub$t.end))
      #     ## Adjust target ranges given the size difference with respect to query ranges
      #     range.offset <- diff(q.range) - diff(t.range)
      #     t.range[2] <- t.range[2] + range.offset ## Make a start position as offset and change only end position
      #     ## Covert query to target coordinates
      #     paf.sub$q.start.trans <- q2t(x = paf.sub$q.start, q.range = q.range, t.range = t.range)
      #     paf.sub$q.end.trans <- q2t(x = paf.sub$q.end, q.range = q.range, t.range = t.range)
      #     # q.range <- range(c(paf$q.start, paf$q.end))
      #     # t.range <- range(c(paf$t.start, paf$t.end))
      #     # paf$q.start.trans <- q2t(x = paf$q.start, q.range = q.range, t.range = t.range)
      #     # paf$q.end.trans <- q2t(x = paf$q.end, q.range = q.range, t.range = t.range)
      #     paf.l[[i]] <- paf.sub
      # }
      # paf <- do.call(rbind, paf.l)

      ## Use user defined q.range if defined
      if (is.null(q.range) | !is.numeric(q.range) | length(q.range) != 2) {
        q.range <- range(c(paf$q.start.shift, paf$q.end.shift))
      }
      ## Use user defined t.range if defined
      if (is.null(t.range) | !is.numeric(t.range) | length(t.range) != 2) {
        t.range <- range(c(paf$t.start.shift, paf$t.end.shift))
      }
      ## Adjust target ranges given the size difference with respect to query ranges
      range.offset <- diff(q.range) - diff(t.range)
      t.range[2] <- t.range[2] + range.offset ## Make a start position as offset and change only end position
      ## Covert query to target coordinates
      paf$q.start.trans <- q2t(x = paf$q.start.shift, q.range = q.range, t.range = t.range)
      paf$q.end.trans <- q2t(x = paf$q.end.shift, q.range = q.range, t.range = t.range)

    } else {
      paf$q.start.trans <- paf$q.start
      paf$q.end.trans <- paf$q.end
    }

    ## Vectorize data transformation ##
    ## Define x and y coordinates
    #x <- c(rbind(paf$q.start.trans, paf$t.start, paf$q.end.trans, paf$t.end))
    x <- c(rbind(paf$q.start.trans, paf$t.start.shift, paf$q.end.trans, paf$t.end.shift))
    y <- rep(c(1, 2, 1, 2), times = nrow(paf))

    ## If there are more than one query or target sequence use offset
    if (length(unique(paf$q.name)) > 1) {
      offset <- seq(from = 0, to = length(unique(paf$q.name))/10, by = 0.1)[seq_along(unique(paf$q.name))]
      names(offset) <- unique(paf$q.name)
      #q.offset <- rep(offset, table(paf$q.name) * 2)
      y[c(TRUE, FALSE, TRUE, FALSE)] <- y[c(TRUE, FALSE, TRUE, FALSE)] - rep(offset[paf$q.name], each = 2)
    }
    if (length(unique(paf$t.name)) > 1) {
      offset <- seq(from = 0, to = length(unique(paf$t.name))/10, by = 0.1)[seq_along(unique(paf$t.name))]
      names(offset) <- unique(paf$t.name)
      #t.offset <- rep(offset, table(paf$t.name) * 2)
      y[c(FALSE, TRUE, FALSE, TRUE)] <- y[c(FALSE, TRUE, FALSE, TRUE)] + rep(offset[paf$t.name], each = 2)
    }

    ## Offset overlapping TARGET alignments up&down based on start position [TODO rewrite to a helper function]
    if (offset.alignments) {
        if ("bin.id" %in% colnames(paf)) {
            ## Use an un-binned version of PAF alignments (group by bin ID if present)
            to.rep <- table(paf$bin.id)
            offset.l <- list()
            for (i in seq_along(to.rep)) {
                if (i %% 2 == 0) {
                    offset.l[[i]] <- rep(x = c(0, 0.05, 0, 0.05), times = to.rep[i])
                } else {
                    offset.l[[i]] <- rep(x = c(0, 0, 0, 0), times = to.rep[i])
                }
            }
            offset <- do.call(c, offset.l)
        } else {
            offset <- rep(c(0, 0, 0, 0, 0, 0.05, 0, 0.05), times = ceiling(nrow(paf) / 2))[seq_along(y)]
        }
        ## Add offset value to y-axis positions
        y <- y + offset
    }

    group <- rep(seq_len(nrow(paf)), each = 4)
    seq.name <- c(rbind(paf$q.name, paf$t.name, paf$q.name, paf$t.name))
    #seq.pos <- c(rbind(paf$q.start, paf$t.start, paf$q.end, paf$t.end))
    seq.pos <- c(rbind(paf$q.start.shift, paf$t.start.shift, paf$q.end.shift, paf$t.end.shift))
    pos.shift <- c(rbind(paf$q.shift, paf$t.shift, paf$q.shift, paf$t.shift))
    pos.genomic <- c(rbind(paf$q.start, paf$t.start, paf$q.end, paf$t.end))
    seq.id <- rep(c("query", "target", "query", "target"), times = nrow(paf))
    n.match <- rep(paf$n.match, each = 4)
    aln.len <- rep(paf$aln.len, each = 4)
    mapq <- rep(paf$mapq, each = 4)
    aln.id <- rep(paf$aln.id, each = 4)
    ID <- rep(paf$ID, each = 4)
    #seq.pair <- rep(paf$seq.pair, each = 4)
    direction <- rep(paf$strand, each = 4)

    ## Create coordinates table (for plotting)
    coords <- dplyr::bind_cols(
        x = x,
        y = y,
        group = group,
        seq.pos = seq.pos,
        pos.shift = pos.shift,
        pos.genomic = pos.genomic,
        direction = direction,
        seq.name = seq.name,
        seq.id = seq.id,
        n.match = n.match,
        aln.len = aln.len,
        mapq = mapq,
        aln.id = aln.id,
        ID = ID,
        #seq.pair = seq.pair
    )

    ## Add user defined column name if defined
    if (!is.null(add.col)) {
        if (nchar(add.col) > 0 & add.col %in% colnames(paf)) {
            new.col <- rep(paf[, eval(add.col), drop = TRUE], each = 4)
            coords <- coords %>% tibble::add_column(!!add.col := new.col)
        } else {
            warning("User defined column name in 'add.column' is not present in submitted 'paf.table', skipping !!!")
        }
    }
    return(coords)
}
