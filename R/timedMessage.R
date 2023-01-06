messageU <- function(..., underline = "=", overline = "=") {
    x <- paste0(..., collapse = "")
    if (!is.null(overline)) {
        message(rep(overline, nchar(x)))
    }
    message(x)
    if (!is.null(underline)) {
        message(rep(underline, nchar(x)))
    }
}


startTimedMessage <- function(...) {
    x <- paste0(..., collapse = "")
    message(x, appendLF = FALSE)
    ptm <- proc.time()
    return(ptm)
}


stopTimedMessage <- function(ptm) {
    time <- proc.time() - ptm
    message(" ... ", round(time[3], 2), "s")
}
