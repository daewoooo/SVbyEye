

readPaf <- function(paf.file=NULL, include.paf.tags=TRUE, restrict.paf.tags=c('NM', 'cg')) {
  ## Check user input ##
  if (is.null(paf.file)) {
    stop("Path to a PAF file to load is not defined!!!")
  }
  
  ## Load PAF file ##
  if (file.exists(paf.file)) {
    message("Loading PAF file: ", paf.file)
    ## Read PAF lines
    paf.lines <- readLines(paf.file)
    if (include.paf.tags) {
      fields <- stringr::str_split(paf.lines, "\t")
    } else {
      fields <-  stringr::str_split(paf.lines, "\t", 13)
    }
    paf.fields <- S4Vectors::lapply(fields, "[", 1:12)
    field.names <-  c('q.name', 'q.len', 'q.start', 'q.end', 'strand', 't.name', 't.len', 't.start', 't.end', 'n.match', 'aln.len', 'mapq')
    
    for (i in seq_along(paf.fields)) {
      attr(paf.fields[[i]], "names") <- field.names
    }
    paf <- dplyr::bind_rows(paf.fields)
    cols.num <- c(2, 3, 4, 7:12)
    paf[cols.num] <- S4Vectors::sapply(paf[cols.num], as.numeric)
    
    if (include.paf.tags) {
      if (any(lengths(fields) > 12)) {
        paf.tags <- S4Vectors::sapply(fields, function(x) paste(x[13:length(x)]))
        paf <- dplyr::bind_cols(paf, processPafTags(paf.tags = paf.tags, restrict.paf.tags = restrict.paf.tags))
      }
    }
    return(paf)
  } else {
    stop(paste0("PAF file ", paf.file, " doesn't exists !!!"))
    return(NULL)
  }  
}  

processPafTags <- function(paf.tags, restrict.paf.tags=c('NM', 'cg')) {
  ## Make sure restrict.paf.tags takes only expected values
  allowed.tags <- c('tp', 'cm', 's1', 's2', 'NM', 'MD', 'AS', 'SA', 'ms', 'nn', 'ts', 'cg', 'cs', 'dv', 'de', 'rl')
  restrict.paf.tags <- restrict.paf.tags[restrict.paf.tags %in% allowed.tags]
  if (length(restrict.paf.tags) == 0) {
    message(paste0("Submitted 'restrict.paf.tags' are not present in the allowed set of PAF tags: ", paste(allowed.tags, collapse = '; ')))
    restrict.paf.tags <- NULL
  }
  
  tags <- character()
  to.numeric <- integer()
  res <- list()
  n <- length(paf.tags)
  t.idx <- 0
  for(ali.idx in seq_along(paf.tags)) {
    split.tags <- stringr::str_split(paf.tags[[ali.idx]], ":")
    for (tag in split.tags) {
      if (!is.null(restrict.paf.tags) & tag[1] %in% restrict.paf.tags) {
        if ( !(tag[1] %in% tags)) {
          t.idx <- t.idx + 1
          if ( tag[2] %in% c("f", "H", "i")) {
            to.numeric <- c(to.numeric, t.idx)
          }
          res[[ tag[1] ]] <- rep(NA, n)
          tags <- c(tags, tag[1])
        }
        res[[ tag[1] ]][ali.idx] <- tag[3]
      }  
    }
  }
  
  for(i in to.numeric){
    res[[i]] <- as.numeric(res[[i]])
  }
  dplyr::bind_rows(res)
}
