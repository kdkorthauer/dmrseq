#' Filter loci for coverage
#' 
#' Helper function to remove loci that have low coverage. 
#' 
#' @param bs a \code{BSseq} object
#' 
#' @param minCoverage the minimum coverage value allowed for any given loci.
#'   The default value (recommended) is 1. 
#'   
#' @param numSamples an option that specifies how many samples in each group 
#' must have \code{minCoverage} in order to keep a loci. The default option is
#' "all", which means that all samples in each group must have
#' \code{minCoverage} in order to keep a loci. The other option is "one" 
#' which keeps all loci that have have a coverage of at least \code{minCoverage}
#' in at least one sample per group. The \code{testCovariate} must be a two
#' group comparison if this argument is "one".
#' 
#' @param testCovariate Character value or vector indicating which variables
#' (column names) in \code{pData(bs)} to test
#'  for association of methylation levels. 
#'  Can alternatively specify an integer value or vector indicating
#'  which of columns of
#'  \code{pData(bs)} to use. This is used to construct the 
#'  design matrix for the test statistic calculation. Required if 
#'  \code{numSamples} is not "all", and required to be a two group comparison.
#'   
#' @return a \code{BSseq} object with the offending loci removed
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
filterLoci <- function(bs, minCoverage = 1, numSamples="all", 
                       testCovariate=NULL) {
    if (!(numSamples %in% c("all", "one"))){
      stop("numSamples must be either `all` or `one`")
    }
    
    if(numSamples=="one"){
      stopifnot(!is.null(testCovariate))
      
      coeff <- 2:(2 + length(testCovariate) - 1)
      testCov <- pData(bs)[, testCovariate]
      design <- model.matrix(~testCov)
        
      stopifnot(length(unique(testCov))==2)
      
      message("Filtering out loci with coverage less than ", minCoverage,
              " read in all samples within testCovariate group")
      Cov <- as.matrix(getCoverage(bs[,design[, coeff] == 0], type = "Cov"))
      nLoci.original <- nrow(Cov)
      which.zero <- which(rowSums(Cov == 0) == sum(design[, coeff] == 0))
      
      Cov <- as.matrix(getCoverage(bs[,design[, coeff] == 1], type = "Cov"))
      which.zero <- unique(c(which.zero, 
                         which(rowSums(Cov == 0) == sum(design[, coeff] == 1))))
      
      rm(Cov)
      nLoci.removed <- length(which.zero)
      if (length(which.zero) > 0) {
        bs <- bs[-which.zero]
      }
      message(paste0("Removed ", nLoci.removed, " out of ", 
                     nLoci.original, " loci"))
    }
    
    if (numSamples=="all"){ 
      message("Filtering out loci with coverage less than ", minCoverage,
              " read in at least one sample")
    
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
    }
    
    return(bs)
}
