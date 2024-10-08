% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pafToFasta.R
\name{paf2FASTA}
\alias{paf2FASTA}
\title{Export FASTA sequences from a set of alignments reported in PAF formatted file.}
\usage{
paf2FASTA(
  paf.table,
  alignment.space = "query",
  order.by = "query",
  bsgenome = NULL,
  asm.fasta = NULL,
  majority.strand = NULL,
  revcomp = NULL,
  report.longest.aln = FALSE,
  report.query.name = NULL,
  concatenate.aln = TRUE,
  fasta.save = NULL,
  return = "fasta"
)
}
\arguments{
\item{paf.table}{A \code{data.frame} or \code{tibble} containing a single or multiple PAF record(s) with 12 mandatory columns
along with CIGAR string defined in 'cg' column.}

\item{alignment.space}{What alignment coordinates should be exported as FASTA, either 'query' or 'target' (Default : `query`).}

\item{order.by}{Order alignment either by `query` or `target` coordinates.}

\item{bsgenome}{A \pkg{\link[BSgenome]{BSgenome-class}} object of reference genome to get the genomic sequence from.}

\item{asm.fasta}{An assembly FASTA file to extract DNA sequence from defined PAF alignments.}

\item{majority.strand}{A desired majority strand directionality to be reported.}

\item{revcomp}{If set to \code{TRUE} FASTA sequence will be reverse complemented regardless of value defined in `majority.strand`.}

\item{report.longest.aln}{If set to \code{TRUE} only the sequence with the most aligned bases will be reported in final FASTA file.}

\item{report.query.name}{A single query (contig) name/id to be reported as FASTA sequence.}

\item{concatenate.aln}{Set to \code{TRUE} if multiple aligned contigs should be concatenated by 100 N's in to a single FASTA sequence (Default : `TRUE`).}

\item{fasta.save}{A path to a filename where to store final FASTA file.}

\item{return}{Set to either 'fasta' or 'index' to return either FASTA in \code{\link{DNAStringSet-class}} object or region index in \code{\link{GRanges-class}} object is returned.}
}
\value{
A \code{\link{DNAStringSet-class}} object with exported sequence
}
\description{
Export FASTA sequences from a set of alignments reported in PAF formatted file.
}
\examples{
## Get PAF to process ##
paf.file <- system.file("extdata", "test4.paf", package = "SVbyEye")
## Read in PAF
paf.table <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = "cg")
## Get FASTA using query alignment coordinates ##
## Define assembly FASTA to get the sequence from
asm.fasta <- system.file("extdata", "test4_query.fasta", package = "SVbyEye")
paf2FASTA(paf.table = paf.table, alignment.space = "query", asm.fasta = asm.fasta)
\donttest{
## Get FASTA using target alignment coordinates ##
## Define BSgenome object to get the sequence from
paf2FASTA(
    paf.table = paf.table, alignment.space = "target",
    bsgenome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
)
}
}
\author{
David Porubsky
}
