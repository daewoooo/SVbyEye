% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processAssemblyAlignments.R
\name{collapseGaps}
\alias{collapseGaps}
\title{This function will take in a \code{\link{GRanges-class}} or \code{\link{GRangesList-class}} object of alignments of a
single contig to a reference and collapses alignment gaps based user specified maximum allowed gap.}
\usage{
collapseGaps(ranges, max.gap = 1e+05)
}
\arguments{
\item{ranges}{A \code{\link{GRanges-class}} or \code{\link{GRangesList-class}} object of regions of a single or multiple contigs aligned to a reference.}

\item{max.gap}{A maximum length of a gap within a single contig alignments to be collapsed.}
}
\value{
A \code{\link{GRanges-class}} object.
}
\description{
This function will take in a \code{\link{GRanges-class}} or \code{\link{GRangesList-class}} object of alignments of a
single contig to a reference and collapses alignment gaps based user specified maximum allowed gap.
}
\author{
David Porubsky
}
