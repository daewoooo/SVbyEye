# =======================================================
# Some global variables that can be used in all functions
# =======================================================

#' @import utils
NULL

## Set global variables to avoid NSE notes in R CMD check
utils::globalVariables(c(
    "q.name", "q.len", "q.start", "q.end", "strand",
    "t.name", "t.len", "t.start", "t.end", "n.match",
    "aln.len", "mapq", "aln.id", "bin.id"
))
