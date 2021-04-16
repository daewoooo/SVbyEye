#' Plot genome-wide density of single nucleotide or larger variants across the genome or separate chromosomes. 
#'
#' @param infile A genome positions of variants to be visualized as binned density, thus far only BED file format is allowed.
#' @param layout An overall layout of the ideogram, either 'horizontal' or 'vertical'.
#' @param reverse.x Set to \code{TRUE} if x-axis should be horizontally reversed.
#' @param title A character string to be used as plot title.
#' @param blacklist A \code{\link{GRanges-class}} object of regions to be removed.
#' @param bsgenome A \code{BSgenome} object to provide chromosome lengths for plotting.
#' @param ref.fasta A user defined reference FASTA file from which chromosome lengths will be extracted.
#' @param gap.cent.bed A set of ranges in BED format to be highlighted as centromeres and/or gaps as white rectangles.
#' @inheritParams makeBins
#' @import ggplot2
#' @importFrom data.table fread
#' @importFrom BiocGenerics as.data.frame
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom utils read.table
#' @importFrom Rsamtools FaFile indexFa scanFaIndex
#' @importFrom IRanges subsetByOverlaps
#' @importFrom stats quantile
#' @importFrom wesanderson wes_palette
#' @return A \code{ggplot} object.
#' @author David Porubsky
#' @export
#' 
plotVariantDensity <- function(infile=NULL, chromosomes=NULL, binsize=20000, stepsize=200000, layout='horizontal', reverse.x=FALSE, title=NULL, blacklist=NULL, bsgenome=NULL, ref.fasta=NULL, gap.cent.bed=NULL) {

  ## Load the data
  if (!is.null(infile)) { 
    if (file.exists((infile))) {
      var.dt <- data.table::fread(infile, sep='\t', stringsAsFactors = FALSE)
      if (ncol(svs.dt) >= 3) {
        var.dt <- var.dt[,c(1:3)]
        colnames(var.dt) <- c('seqnames', 'start', 'end')
        var.gr <- GenomicRanges::makeGRangesFromDataFrame(var.dt, keep.extra.columns = FALSE)
      } else {
        stop("User defined input file, '", infile, "' doesn't contain required columns ('chrom.id', 'start', 'end') !!!")
      }
    } else {
      stop("User defined input file, '", infile, "' doesn't exists!!!")
    }
  } else {
    stop("No input file defined to plot!!!")
  }
  
  ## Keep only user defined chromosomes
  chroms.in.data <- GenomeInfoDb::seqlevels(var.gr)
  if (is.null(chromosomes)) {
    chromosomes <- chroms.in.data
  }
  chroms2use <- intersect(chromosomes, chroms.in.data)
  if (length(chroms2use) == 0) {
    stop("Defined 'chromosomes' not present in the submitted BED file!!!")
  }
  
  ## Check input format of blacklisted regions and load a BED file if needed
  if (!is.null(blacklist) & !class(blacklist) == 'GRanges') {
    if (file.exists(blacklist)) {
      blacklist.df <- utils::read.table(file = blacklist, header = FALSE, sep = '\t', stringsAsFactors = FALSE)
      if (ncol(blacklist.df) >= 3) {
        blacklist.df <- blacklist.df[,c(1:3)]
        colnames(blacklist.df) <-  c('seqnames', 'start', 'end')
        blacklist <- GenomicRanges::makeGRangesFromDataFrame(blacklist.df)
      } else {
        blacklist <- NULL
        warning("The BED file, '", blacklist, "' doesn't contain required fields ('chr.name', 'start', 'end').")
      }
    } else {
      blacklist <- NULL
      warning("The BED file, '", blacklist, "' doesn't exists.")
    }
  }
  
  ## Check input format of regions to be masked on the ideogram
  if (file.exists(gap.cent.bed)) {
    gap.cent.bed.df <- utils::read.table(file = gap.cent.bed, header = FALSE, sep = '\t', stringsAsFactors = FALSE)
    if (ncol(gap.cent.bed.df) >= 3) {
      gap.cent.bed.df <- gap.cent.bed.df[,c(1:3)]
    } else {
      gap.cent.bed <- NULL
      warning("The BED file, '", gap.cent.bed, "' doesn't contain required fields ('chr.name', 'start', 'end').")
    }
  } else {
    gap.cent.bed <- NULL
    warning("The BED file, '", gap.cent.bed, "' doesn't exists.")
  }
  
  ## Check user input
  ## Load BSgenome
  if (class(bsgenome) != 'BSgenome') {
    if (is.character(bsgenome)) {
      bsgenome <- tryCatch({
        suppressPackageStartupMessages(library(bsgenome, character.only=TRUE))
        bsgenome <- eval(parse(text=bsgenome)) ## replacing string by object
      }, error = function(e) {return(NULL)})
    } else {
      bsgenome <- NULL
    }
  }
  ## Load reference fasta
  if (!is.null(ref.fasta)) {
    if (class(Rsamtools::FaFile(ref.fasta)) != "FaFile") {
      ref.fasta <- NULL
      warning("User defined reference FASTA file, '", ref.fasta, "' is not in a proper FASTA format!!!")
    }
  }
  
  if (!is.null(bsgenome)) {
    ## Get chromosome length from bsgenome object
    seq.len <- GenomeInfoDb::seqlengths(bsgenome)[chroms2use]
  } else if (!is.null(ref.fasta)) {
    ## Check if submitted fasta file is indexed
    ref.fasta.idx <- paste0(ref.fasta, ".fai")
    if (!file.exists(ref.fasta.idx)) {
      message("Index file for '", ref.fasta, "' not available, indexing ...")
      fa.idx <- Rsamtools::indexFa(file = ref.fasta)
    }
    ## Get FASTA index
    fa.idx <- Rsamtools::scanFaIndex(ref.fasta)
    ## Get chromosome length from reference fasta
    seq.len <- GenomeInfoDb::seqlengths(fa.idx)[chroms2use]
  }

  ## Bin the data
  if (!is.null(bsgenome)) {
    bins <- makeBins(bsgenome = bsgenome, chromosomes = chroms2use, binsize = binsize, stepsize = stepsize)
  } else if (!is.null(ref.fasta)) {
    bins <- makeBins(fai = paste0(ref.fasta, ".fai"), chromosomes = chroms2use, binsize = binsize, stepsize = stepsize)
  }
  
  ## Filter variants from blacklisted regions
  if (!is.null(blacklist)) {
    var.gr <- suppressWarnings( IRanges::subsetByOverlaps(var.gr, blacklist, invert = TRUE) )
  }
    
  ## Count variants in bins
  bins$counts <- IRanges::countOverlaps(bins, var.gr)
  ## Filter bins with extremely high amount of reads
  outlier <- stats::quantile(bins$counts, 0.99)
  ## Set outlier bins to the limit
  bins$counts[bins$counts >= outlier] <- outlier
  
  ## Prepare the ideogram
  ideo.df <- data.frame(seqnames=names(seq.len), length=seq.len)
  ideo.df$seqnames <- factor(ideo.df$seqnames, levels=chroms2use)
  ideo.df$levels <- 1:length(seq.len)

  if (layout == 'horizontal') {
    ideo <- ggplot2::ggplot() + 
      ggplot2::geom_rect(data=ideo.df, ggplot2::aes(xmin=0, xmax=length, ymin=0, ymax=1), color='black', fill='white') 
  } else if (layout == 'vertical') {
    ideo <- ggplot2::ggplot() + 
      ggplot2::geom_rect(data=ideo.df, ggplot2::aes(ymin=0, ymax=length, xmin=0, xmax=1), color='black', fill='white') 
  }  
  
  ## Get color palette
  pal <- wesanderson::wes_palette("Zissou1", 100, type = "continuous")

  ## Plotting themes
  theme_horizontal <- theme(axis.title.y=element_blank(),
                            axis.text.y=element_blank(),
                            axis.ticks.y=element_blank(),
                            panel.grid.major = element_blank(), 
                            panel.grid.minor = element_blank(),
                            strip.text.y.left = element_text(angle = 0),
                            panel.background = element_blank(),
                            strip.background = element_rect(colour="white", fill="white"))
  
  theme_vertical <- theme(legend.position ="top",
                          axis.line = element_blank(),
                          axis.text.x=element_blank(), 
                          axis.ticks.x=element_blank(),   
                          strip.text.x = element_text(angle = 0),
                          panel.grid.major = element_blank(), 
                          panel.grid.minor = element_blank(),
                          panel.background = element_blank(),
                          strip.background = element_rect(colour="white", fill="white"))
  
  ## Prepare data for plotting
  bins.df <- BiocGenerics::as.data.frame(bins)
  
  ## Get chromosome breaks and labels
  max.len <- signif(max(ideo.df$length), digits = 2)
  breaks <- seq(from = 0, to = max.len, length.out = 6)
  labels <- breaks / 1000000
  labels <- paste0(labels, 'Mb')
  
  ## Make ideogram
  if (layout == 'horizontal') {
    plt <- ideo + 
      geom_rect(data=bins.df, aes(xmin=start, xmax=end, ymin=0, ymax=1, fill=counts)) +
      scale_fill_gradientn(colours = pal, guide="none") +
      facet_grid(seqnames ~ ., switch = 'y') +
      scale_x_continuous(expand = c(0,0), breaks = breaks, labels = labels, name="") +
      theme_horizontal
  } else if (layout == 'vertical') {
    plt <- ideo + 
      geom_rect(data=bins.df, aes(ymin=start, ymax=end, xmin=0, xmax=1, fill=counts)) +
      scale_fill_gradientn(colours = pal, guide="none") +
      facet_grid(. ~ seqnames, switch = 'x') +
      scale_y_continuous(expand = c(0,0), breaks = breaks, labels = labels, name="") +
      theme_vertical
  }  
  
  ## Overlay positions of centromeres and gaps as white rectangles
  if (!is.null(gap.cent.bed)) {
    if (layout == 'horizontal') {
      plt <- plt + geom_rect(data=cent.gap.mask.df, aes(xmin=start, xmax=end, ymin=0, ymax=1), fill='white')
    } else if (layout == 'vertical') {
      plt <- plt + geom_rect(data=cent.gap.mask.df, aes(ymin=start, ymax=end, xmin=0, xmax=1), fill='white')
    }  
  }
  ## Reverse x axis if required
  if (reverse.x & layout == 'horizontal') {
    plt <- suppressMessages( plt + scale_x_reverse() )
  }
  ## Add title
  if (!is.null(title) & nchar(title > 0)) {
    plt <- plt + ggplot2::ggtitle(title)
  }
  
  ## Highlight a user defined loci [TODO: allow region highlights]
  #hla.locus <- GRanges(seqnames = 'chr6', ranges=IRanges(start=28510120, end=33480577))
  #hla.locus.df <- as.data.frame(hla.locus)
  # + geom_rect(data=hla.locus.df, aes(xmin=start, xmax=end, ymin=-0.1, ymax=1.1), color="green", fill=NA) +
  
  ## Return final plot
  return(plt)  
}
