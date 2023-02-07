#' Read nucmer coordinates file in to PAF formatted table
#'
#' This function takes alignments produced by nucmer (parameters: --mum --coords) and loads
#' reported coordinates file (suffix .coords) into a PAF formatted table (see PAF specification).
#'
#' @param nucmer.coords A path to a nucmer coordinates file containing alignments to be loaded.
#' @importFrom utils read.table
#' @importFrom tibble tibble
#' @author David Porubsky
#' @export
#' @examples
#' ## Get nucmer coordinates to read in ##
#' nucmer.coords <- system.file("extdata", "HG00438_1_nucmer.coords", package = "SVbyEye")
#' ## Read in nucmer coordinates
#' nucmer2PAF(nucmer.coords=nucmer.coords)
#' 
nucmer2PAF <- function(nucmer.coords=NULL) {
  ## Check is path to nucmer coordinates has been defined and if the file exists
  if (!is.null(nucmer.coords)) {
    if (file.exists(nucmer.coords)) {
      ## Read in coordinates from nucmer output
      coords <- utils::read.table(nucmer.coords, skip=5, stringsAsFactors = FALSE)
      coords.tb <- tibble::tibble(
        q.name = coords$V12,
        q.len = max(c(coords$V1, coords$V2)),
        q.start = coords$V1,
        q.end = coords$V2,
        strand = ifelse(coords$V4 < coords$V5, '+', '-'),
        t.name = coords$V13,
        t.len = max(c(coords$V4, coords$V5)),
        t.start = pmin(coords$V4, coords$V5),
        t.end = pmax(coords$V4, coords$V5),
        n.match = coords$V10,
        aln.len = coords$V8,
        mapq = NA
      )
    } else {
      coords.tb <- NULL
    }
    return(coords.tb)
  }
}
