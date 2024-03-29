% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/breakPafAlignment.R
\name{breakPaf}
\alias{breakPaf}
\title{A wrapper function for \code{\link{breakPafAlignment}} expanding multiple PAF alignments
into a set of matching bases between query and target sequence.}
\usage{
breakPaf(
  paf.table = NULL,
  min.deletion.size = 50,
  min.insertion.size = 50,
  collapse.mismatches = TRUE,
  report.sv = TRUE
)
}
\arguments{
\item{paf.table}{A \code{data.frame} or \code{tibble} containing a single or multiple PAF record(s) with 12 mandatory columns
along with CIGAR string defined in 'cg' column.}

\item{min.deletion.size}{A minimum size (in base pairs) of a deletion to be retained.}

\item{min.insertion.size}{A minimum size (in base pairs) of an insertion to be retained.}

\item{collapse.mismatches}{Set to \code{TRUE} if mismatches should be collapsed in order expand matched regions.}

\item{report.sv}{Set to \code{TRUE} if to report also ranges of deleted and inserted bases.}
}
\value{
A \code{list} of \code{tibble} objects storing matched ('M') alignments as well as structurally variable ('SV') bases if 'report.sv' is TRUE.
}
\description{
A wrapper function for \code{\link{breakPafAlignment}} expanding multiple PAF alignments
into a set of matching bases between query and target sequence.
}
\examples{
## Get PAF to break
paf.file <- system.file("extdata", "test1.paf", package = "SVbyEye")
## Read in PAF alignment
paf.table <- readPaf(paf.file = paf.file)
## Break PAF alignment at indels of 20 bp and longer
breakPaf(paf.table = paf.table, min.deletion.size = 20, min.insertion.size = 20)

}
\author{
David Porubsky
}
