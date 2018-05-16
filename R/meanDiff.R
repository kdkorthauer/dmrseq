#' Function to compute raw mean methylation differences
#' 
#' This function calculates raw mean methylation differences for the
#' covariate of interest over a set of DMRs (or regions of interest), 
#' assuming a simple two-group comparison.
#' 
#' @param bs a \code{BSseq} object
#' @param dmrs a data.frame with one row per DMR. This can be in the format
#' of \code{dmrseq} output, but at least should contain the indexStart and 
#' indexEnd values of the regions of interest.
#' @param testCovariate a character indicating the covariate of interest in
#' the \code{pData} slot of \code{bs}.
#' 
#' @return numeric vector of raw mean methylation differences.
#' 
#' @importFrom DelayedMatrixStats rowMeans2
#' 
#' @export
#' 
#' @examples
#' 
#' data(BS.chr21)
#' data(dmrs.ex)
#' rawDiff <- meanDiff(BS.chr21, dmrs=dmrs.ex, testCovariate="CellType")
#' 
meanDiff <- function(bs, dmrs, testCovariate) {
  # convert covariates to column numbers if characters
  if (is.character(testCovariate)) {
    testCovariate <- which(colnames(pData(bs)) == testCovariate)
    if (length(testCovariate) == 0) {
      stop("testCovariate not found in pData(). ", 
           "Please specify a valid testCovariate")
    }
  }
  
  coeff <- seq(2,(2 + length(testCovariate) - 1))
  testCov <- pData(bs)[, testCovariate]
  if (length(unique(testCov)) == 1) {
    message("Warning: only one unique value of the specified ", 
            "covariate of interest.  Assuming null comparison and ", 
            "splitting sample group into two equal groups")
    testCov <- rep(1, length(testCov))
    testCov[seq_len(round(length(testCov)/2))] <- 0
  }

  design <- model.matrix(~testCov)
  colnames(design)[coeff] <- colnames(pData(bs))[testCovariate]
  
  if (length(unique(design[, coeff])) != 2) {
    message("Not a two-group comparison. Can't compute simple mean ",
            "methylation differences. ", 
            "Returning beta estimates instead")
    return(dmrs$beta)
  } else {
    prop.mat <- getCoverage(bs, type = "M") / 
                getCoverage(bs, type = "Cov")
    levs <- unique(design[, coeff])
    
    indexRanges <- IRanges(start(dmrs$index), end(dmrs$index))
    prop.mat.dmr <- extractROWS(prop.mat, indexRanges)
    prop.mat1.means <- DelayedMatrixStats::rowMeans2(prop.mat.dmr[,
                                 design[, coeff] == levs[which.min(levs)]],
                                 na.rm=TRUE)
    prop.mat2.means <- DelayedMatrixStats::rowMeans2(prop.mat.dmr[,
                                 design[, coeff] == levs[which.max(levs)]],
                                 na.rm=TRUE)
    
    meanDiff <- IRanges::mean(IRanges::relist(prop.mat2.means - prop.mat1.means,
                                          indexRanges), na.rm=TRUE)
    
    return(meanDiff)
  }
}
