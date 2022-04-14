#' Read GAF from an input file
#' 
#' This function takes an GAF output file from minigraph and loads the file and 
#' parse all the alignments.
#'
#' @param gaf.file A path to a GAF file containing graph alignments.
#' @importFrom stringr str_split
#' @importFrom dplyr bind_cols bind_rows
#' @importFrom S4Vectors lapply
#' @author David Porubsky
#' @export
readGaf <- function(gaf.file = NULL) {
  ## Check user input ##
  if (is.null(gaf.file)) {
    stop("Path to a GAF file to load is not defined!!!")
  }
  
  ## Load GAF file ##
  con <- file(filepath, "r")
  paths <- list()
  gaf.aln <- NULL
  gaf <- NULL
  path.aligns <- list()
  while (TRUE) {
    line <- readLines(con, n = 1)
    if (length(line) != 0) {
      fields <- stringr::str_split(line, pattern = "\t")
      if (!grepl(line, pattern = '\\*')) {
        ## Store if path alignments completely loaded
        if (!is.null(gaf.aln)) {
          path.aligns <- S4Vectors::lapply(path.aligns, "[", 2:9)
          path.aligns.names <- c('s.name', 's.len', 'n.minimizers', 'divergence', 's.start', 's.end', 'p.start', 'p.end')
          for (i in seq_along(path.aligns)) {
            attr(path.aligns[[i]], "names") <- path.aligns.names
          }
          path.aligns.tab <- dplyr::bind_rows(path.aligns)
          ## Get strand
          dir.code <- gsub(path.aligns.tab$s.name, pattern = 's\\d+', replacement = '')
          strand <- dplyr::recode(dir.code, '>' = '+', '<' = '-')
          path.aligns.tab$s.strand <- strand
          ## Add path ids
          path.ids <- stringr::str_split(gaf.aln$path, pattern = "<|>")[[1]]
          path.ids <- path.ids[nchar(path.ids) > 0]
          ## Get continuous paths
          seq.ids <- as.numeric(gsub(path.aligns.tab$s.name, pattern = '>s|<s', replacement = ''))
          seq.ids.subseq <- cumsum(c(TRUE, abs(diff(seq.ids)) != 1))
          seq.ids.subseq.runs <- S4Vectors::runLength( S4Vectors::Rle(seq.ids.subseq) )
          path.aligns.tab$t.name <- rep(path.ids, seq.ids.subseq.runs)
          ## Merge tables
          gaf <- dplyr::bind_cols(gaf.aln, path.aligns.tab)
          gaf$record <- length(paths) + 1
          paths[[length(paths) + 1]] <- gaf
          path.aligns <- list()
        }
        ## Process alignment header
        paf.fields <- S4Vectors::lapply(fields, "[", 1:12)
        field.names <-  c('q.name', 'q.len', 'q.start', 'q.end', 'strand', 'path', 'path.len', 'path.start', 'path.end', 'n.match', 'aln.block.len', 'mapq')
        attr(paf.fields[[1]], "names") <- field.names
        gaf.aln <- dplyr::bind_rows(paf.fields)
      } else {
        path.aligns[[length(path.aligns) + 1]] <- fields[[1]]
      }
    } else {
      path.aligns <- S4Vectors::lapply(path.aligns, "[", 2:9)
      path.aligns.names <- c('s.name', 's.len', 'n.minimizers', 'divergence', 's.start', 's.end', 'p.start', 'p.end')
      for (i in seq_along(path.aligns)) {
        attr(path.aligns[[i]], "names") <- path.aligns.names
      }
      path.aligns.tab <- dplyr::bind_rows(path.aligns)
      ## Get strand
      dir.code <- gsub(path.aligns.tab$s.name, pattern = 's\\d+', replacement = '')
      strand <- dplyr::recode(dir.code, '>' = '+', '<' = '-')
      path.aligns.tab$s.strand <- strand
      ## Add path ids
      path.ids <- stringr::str_split(gaf.aln$path, pattern = "<|>")[[1]]
      path.ids <- path.ids[nchar(path.ids) > 0]
      ## Get continuous paths
      seq.ids <- as.numeric(gsub(path.aligns.tab$s.name, pattern = '>s|<s', replacement = ''))
      seq.ids.subseq <- cumsum(c(TRUE, abs(diff(seq.ids)) != 1))
      seq.ids.subseq.runs <- S4Vectors::runLength( S4Vectors::Rle(seq.ids.subseq) )
      path.aligns.tab$t.name <- rep(path.ids, seq.ids.subseq.runs)
      ## Merge tables
      gaf <- dplyr::bind_cols(gaf.aln, path.aligns.tab)
      gaf$record <- length(paths) + 1
      paths[[length(paths) + 1]] <- gaf
      break
    }
  }
  close(con)
  paths.tbl <- dplyr::bind_rows(paths)
  paths.tbl$s.name <- gsub(paths.tbl$s.name, pattern = '<|>', replacement = '')
  num.cols <- c(2, 3, 4, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19, 20)
  paths.tbl[,num.cols] <- dplyr::bind_cols(S4Vectors::lapply(paths.tbl[num.cols], as.numeric))
  return(paths.tbl)
}