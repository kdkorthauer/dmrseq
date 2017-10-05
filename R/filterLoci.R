#' Filter loci for coverage
#' 
#' Helper function to remove loci that have low coverage. 
#' 
#' @param bs a \code{BSseq} object
#' 
#' @param minCoverage the minimum coverage value allowed for any given loci.
#'   The default value (recommended) is 1. This means that any loci with zero
#'   reads mapping to it in at least one sample are removed.
#'   
#' @return a \code{BSseq} object with the offending loci removed.
#' 
#' @export
#' 
#' @examples 
#' 
#' # load example data 
#' data(BS.chr21)
#' 
#' BS.chr21 <- filterLoci(BS.chr21)
#' 
filterLoci <- function(bs, minCoverage = 1) {
    message(paste0("Filtering out loci with coverage less than ", minCoverage,
                   " read in at least one sample"))
    
    Cov <- as.matrix(getCoverage(bs, type = "Cov"))
    nLoci.original <- nrow(Cov)
    which.zero <- which(rowSums(Cov == 0) > 0)
    rm(Cov)
    nLoci.removed <- length(which.zero)
    if (length(which.zero) > 0) {
        bs <- bs[-which.zero]
    }
    message(paste0("Removed ", nLoci.removed, " out of ", 
                   nLoci.original, " loci"))
    return(bs)
}
