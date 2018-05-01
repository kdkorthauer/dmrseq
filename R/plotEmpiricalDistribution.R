#' Plot the empirical distribution of the methylation beta vals or coverage
#' 
#' Uses ggplot2 to plot smoothed density histograms of methylation
#' proportions (beta values), or coverage. Methylation proportion densities
#' are weighted by coverage.
#' The number of curves plotted 
#' will be equal to the number of different values of \code{testCovariate},
#' unless \code{bySample} is TRUE. This can take quite some time to 
#' execute for a large object, so it is recommended to first take a random
#' sample of loci (say one million) before plotting.
#' 
#' @param bs a BSseq object
#'
#' @param testCovariate character specifying the column name of the
#' \code{pData} slot of the BSseq object to include in the plot legend.
#' 
#' @param bySample logical whether to plot a separate line for each sample.
#'  Default value is FALSE (so samples with the same value of 
#'  \code{testCovariate} will be collapsed into the same line).
#' 
#' @param type a character indicating which type of density to plot - the 
#' methylation (beta) values ("M") or the coverage ("Cov"). Default is "M".
#' 
#' @param adj a numeric value for the \code{adjust} parameter to pass to the 
#' \code{geom_line} function. Specifies how smooth the make the function.
#' 
#' @return a ggplot object
#' 
#' @importFrom locfit locfit lp
#' 
#' @export
#' 
#' @examples 
#' 
#' data(BS.chr21)
#' 
#' # plot beta values by sample group
#' plotEmpiricalDistribution(BS.chr21, testCovariate="CellType")
#' 
plotEmpiricalDistribution <- function(bs, 
                                      testCovariate = NULL,
                                      bySample = FALSE,
                                      type="M",
                                      adj = 2.5) {
    #satisfy check
    M <- Cov <- group <- NULL
  
    if (!(type %in% c("M", "Cov"))){
      stop("type must be either M or Cov")
    }
  
    meth.mat <- getCoverage(bs, type = "M")
    unmeth.mat <- getCoverage(bs, type = "Cov") - meth.mat
    
    meth.levelsm <- data.frame(meth.mat/
                                 (meth.mat + unmeth.mat))
    cov.matm <- data.frame((meth.mat + unmeth.mat))
    
    if (!is.null(testCovariate)) {
        if (sum(grepl(testCovariate, colnames(pData(bs)))) == 0) {
            stop("Error: no column in pData() found ", 
                "that matches the testCovariate")
        } else if (length(grep(testCovariate, colnames(pData(bs)))) > 1) {
            stop("Error: testCovariate matches more ",
                "than one column in pData()")
        }
        mC <- grep(testCovariate, colnames(pData(bs)))
    }
    
    meth.levelsm <- utils::stack(meth.levelsm)
    colnames(meth.levelsm)[1] <- "M"
    meth.levelsm$Cov <- utils::stack(cov.matm)$values
    
    meth.levelsm$group <- 
      unlist(lapply(seq_len(ncol(bs)), function(x) 
        rep(pData(bs)[x,mC], nrow(bs))))
    meth.levelsm$sample <- sort(rep(seq_len(ncol(bs)), nrow(bs)))
    
    if (!bySample){
      if (type=="M"){
        # compute weights - sum over all samples in group ##
        covtots <- rep(NA, ncol(cov.matm))
        names(covtots) <- pData(bs)[,mC]
        for(l in unique(pData(bs)[,mC])){
          covtots[names(covtots) == l] <- sum(colSums(cov.matm)[pData(bs)[,mC] == l]) 
        }
        
        wt.matm <- data.frame(t(t(cov.matm) / covtots)) 
        meth.levelsm$wt <- utils::stack(wt.matm)$values 
        
        p1 <- ggplot(meth.levelsm, 
                 aes(M, colour = group, group = group, weight = wt)) + 
          geom_line(adjust = adj, alpha = 0.6, stat = "density", size = 1.3) + 
        xlab("Methylation Proportion") + 
        theme_bw()  # colour=condition
      }else{
        p1 <- ggplot(meth.levelsm, aes(Cov+0.1, 
                               colour = group, group = group)) + 
          geom_line(adjust = adj, alpha = 0.6, stat = "density", size = 1.3) +
          scale_x_continuous(trans="log2") +
        xlab("Coverage") + 
          theme_bw()  # colour=condition
      }
    }else{
      if (type=="M"){
        wt.matm <- data.frame(t(t(cov.matm) / colSums(cov.matm))) 
        meth.levelsm$wt <- utils::stack(wt.matm)$values 
        
        p1 <- ggplot(meth.levelsm, 
                     aes(M, colour = group, group = sample, weight = wt)) + 
          geom_line(adjust = adj, alpha = 0.6, stat = "density", size = 1.3) +
          xlab("Methylation Proportion") + 
          theme_bw()  # colour=condition
      }else{
        p1 <- ggplot(meth.levelsm, aes(Cov+0.1, 
                                       colour = group, group = sample)) + 
          geom_line(adjust = adj, alpha = 0.6, stat = "density", size = 1.3) +
          scale_x_continuous(trans="log2") +
          xlab("Coverage") + 
          theme_bw()  # colour=condition
      }
    }
    return(p1)
}

