% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotMiro.R
\name{plotMiro}
\alias{plotMiro}
\title{Visualize PAF alignments.}
\usage{
plotMiro(
  paf.table,
  min.deletion.size = NULL,
  min.insertion.size = NULL,
  highlight.sv = NULL,
  binsize = NULL,
  color.by = "direction",
  perc.identity.breaks = c(90, 95, 99, 99.5, 99.6, 99.7, 99.8, 99.9),
  color.palette = NULL,
  outline.alignments = FALSE,
  offset.alignments = FALSE,
  add.alignment.arrows = TRUE,
  target.region = NULL,
  genomic.scale = "bp"
)
}
\arguments{
\item{paf.table}{A \code{data.frame} or \code{tibble} containing a single or multiple PAF record(s) with 12 mandatory columns
along with CIGAR string defined in 'cg' column.}

\item{min.deletion.size}{A minimum size (in base pairs) of a deletion to be retained.}

\item{min.insertion.size}{A minimum size (in base pairs) of an insertion to be retained.}

\item{highlight.sv}{Visualize alignment embedded structural variation either as an outlined ('outline') or filled ('fill') miropeats.}

\item{binsize}{A size of a bin in base pairs to split a PAF alignment into.}

\item{color.by}{Color alignments either by directionality ('direction'), fraction of matched base pairs ('identity'),
or a custom column name present in submitted `paf.table`.}

\item{perc.identity.breaks}{Set of percentage values to categorize alignment percent identity into a set of discrete colors.}

\item{color.palette}{A discrete color palette defined as named character vector (elements = colors, names = discrete levels)
to color alignment directionality, `[default: color.palette <- c('-' = 'cornflowerblue', '+' = 'forestgreen')]`.}

\item{outline.alignments}{Set to \code{TRUE} if boundaries of each alignment should be highlighted by gray outline.}

\item{offset.alignments}{Set to \code{TRUE} if subsequent target alignments should be offsetted below and above the midline.}

\item{add.alignment.arrows}{Set to \code{FALSE} if alignment arrows should not be added to the plot.}

\item{target.region}{A user defined target region either as character string ('chr:start-end') or as
a \code{\link{GRanges-class}} object containing a single genomic region to which PAF alignments will be
narrowed down.}

\item{genomic.scale}{Report genomic coordinates in base pairs ('bp') kilobase pairs ('kbp') or megabase pairs ('Mbp') `[default: 'bp']`.}
}
\value{
A \code{ggplot2} object
}
\description{
This function takes PAF output file from minimap2 alignments, and visualize the alignments
in a miropeat style.
}
\examples{
## Get PAF to plot ##
paf.file <- system.file("extdata", "test1.paf", package = "SVbyEye")
## Read in PAF
paf.table <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = "cg")
## Optional steps include PAF filtering and flipping query coordinates
## (see filterPaf and flipPaf function documentation)
## Make a plot ##
## Color by alignment directionality
plotMiro(paf.table = paf.table, color.by = "direction")
## Color by fraction of matched bases in each alignment
plotMiro(paf.table = paf.table, color.by = "identity")
## Use custom color palette to color alignment directionality
plotMiro(paf.table = paf.table, color.palette = c("+" = "azure3", "-" = "yellow3"))
## Change genomic scale to megabase pairs
plotMiro(paf.table = paf.table, genomic.scale = "Mbp")
## Outline PAF alignments
plotMiro(paf.table = paf.table, outline.alignments = TRUE)
## Offset target PAF alignments
plotMiro(paf.table = paf.table, offset.alignments = TRUE)
## Bin PAF alignments into user defined bin and color them by sequence identity (\% of matched bases)
plotMiro(paf.table = paf.table, binsize = 10000)
## Highlight structural variants
paf.file <- system.file("extdata", "test3.paf", package = "SVbyEye")
paf.table <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = "cg")
plotMiro(paf.table = paf.table, min.deletion.size = 50, highlight.sv = "outline")

}
\author{
David Porubsky
}
