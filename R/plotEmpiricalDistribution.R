#' Plot the empirical distribution of the methylation beta vals
#' 
#' Uses ggplot2 to plot smoothed density histograms of methylation
#' proportions (beta values)
#' 
#' @param bs a BSseq object
#' 
#' @param plotFile a character string specifying the filename (and full path,
#' or else file will be saved to the current working directory) for the 
#' pdf file that will contain the output plot.
#' 
#' @param tiss1 character indicating the label of population 1
#' 
#' @param tiss2 character indicating the lable of population 2
#' 
#' @param time character indicating a secondary label of population 1 (e.g.
#' time point if each population is measured at multiple time points)
#' 
#' @param time2 character indicating a secondary label of population 2 (e.g.
#' time point if each population is measured at multiple time points)
#' 
#' @param matchCovariate character specifying the column name of the
#' \code{pData} slot of the BSseq object to include in the plot legend.
#' 
#' @param genomeName character specifying the genome name if plots of
#' beta values by CpG category should be made. See package \code{annotatr}
#' for details of available genomes.
#' 
#' @return None. Produces a plot saved to a pdf file with filename 
#' \code{plotFile}
#' 
#' @importFrom reshape2 melt
#' @importFrom locfit locfit lp
#' 
#' @export
plotEmpiricalDistribution <- function(bs, plotFile, 
                                      tiss1, tiss2, 
                                      time=NULL, time2=NULL,
                                      matchCovariate=NULL,
                                      genomeName=NULL){
  value=NULL # satisfy R cmd check
  variable=NULL # satisfy R cmd check
  
  if (!is.null(genomeName)){								
    # annotate BS with island/shore/etc
    annot_CpG = paste0(c(genomeName, "_cpgs"), collapse="")
    
    # Build the annotations (a single GRanges object)
    
    # Download has a nonzero fail rate; allow up to 5 retries before
    # throwing an error
    for (attempt in 1:5){
      cpg = try(build_annotations(genome = genomeName, annotations = annot_CpG),
                silent=TRUE)  
      if(!is(cpg, 'try-error')){
        message("Download of CpG annotation successful!")
        break;
      }else{
        message(paste0("Download of CpG annotation failed. ",
                       5-attempt, " attempts left"))
      }
    }
    
    ols <- findOverlaps(bs, cpg)
    cpgCat <- mcols(cpg)$type[ols@to]
  }
  
  meth.mat = getCoverage(bs, type = "M")
  unmeth.mat = getCoverage(bs, type = "Cov") - meth.mat
  
  #plot empirical distribution of meth levels (betas)
  # random sample of 1M rows
  rows.plot <- sample(1:nrow((meth.mat)), min(nrow((meth.mat)),1e6), 
                      replace=FALSE)
  meth.levelsm <- data.frame(meth.mat[rows.plot,] / 
                               (meth.mat + unmeth.mat)[rows.plot,])
  cov.matm <- data.frame((meth.mat + unmeth.mat)[rows.plot,])
  
  if(!is.null(matchCovariate)){
    if(sum(grepl(matchCovariate, colnames(pData(bs))))==0){
      stop("Error: no column in pData() found that matches the matchCovariate")
    }else if(length(grep(matchCovariate, colnames(pData(bs)))) > 1){
      stop("Error: matchCovariate matches more than one column in pData()")
    }
    mC <- grep(matchCovariate, colnames(pData(bs)))
    
    if (length(grep(pData(bs)[1,mC], colnames(meth.mat)))==0){
      colnames(meth.levelsm) <- colnames(cov.matm) <- 
        paste0(colnames(meth.levelsm), "_", pData(bs)[,mC])
    }
  }
  
  if (!is.null(genomeName)){		
    cpgCat <- cpgCat[rows.plot]
    meth.levelsm$cpg <- cpgCat
    cov.matm$cpg <- cpgCat	 
  }
  
  meth.levelsm <- melt(meth.levelsm)
  cov.matm <- melt(cov.matm)
  
  #overlay version
  g1 <- grep(tiss1, gsub("\\.", " ", meth.levelsm$variable))
  g2 <- grep(tiss2, gsub("\\.", " ", meth.levelsm$variable))
  
  if(!exists("time")){
    time <- NULL
  }else if(is.null(time)){
    time <- time
  }
  if(!exists("time2")){
    time2 <- NULL
  }else if(is.null(time2)){
    time2 <- time
  }
  
  if (class(time)!="function"){
    if (length(grep("_", time)) < 0){
      t1 <- strsplit(time, "_")[[1]][[2]]
      t1 <- grep(t1, meth.levelsm$variable)
      c1 <- g1[g1 %in% t1]
    }
  }
  if (class(time)!="function"){
    if (length(grep("_", time2)) < 0){
      t2 <- strsplit(time2, "_")[[1]][[2]]
      t2 <- grep(t2, meth.levelsm$variable)
      c2 <- g2[g2 %in% t2]
    }
  }
  if (exists("t1") & exists("t2")){
    meth.levelsm$condition <- paste0(time, "_", tiss1)
    meth.levelsm$condition[t2] <- paste0(time2, "_", tiss2)
    
    meth.levelsm$time <- time
    meth.levelsm$time[t2] <- time2
    
    cov.matm$condition <- paste0(time, "_", tiss1)
    cov.matm$condition[t2] <- paste0(time2, "_", tiss2)
    
    cov.matm$time <- time
    cov.matm$time[t2] <- time2
  }
  
  meth.levelsm$tissue <- tiss1
  meth.levelsm$tissue[g2] <- tiss2
  meth.levelsm$condition <- meth.levelsm$tissue
  
  cov.matm$tissue <- tiss1
  cov.matm$tissue[g2] <- tiss2
  cov.matm$condition <- meth.levelsm$tissue
  
  pdf(file=plotFile, 
      height=3, width=6)	
  
  p0 <- ggplot(meth.levelsm, aes(value, colour=variable, group=variable)) +
    geom_line(adjust=2.5, alpha=0.6, stat="density", size=1.3)+
    xlab("Methylation Proportion") + 
    theme_bw() # colour=condition
  p1 <- ggplot(cov.matm, aes(log(value), colour=variable, group=variable)) +
    geom_line(adjust=2.5, alpha=0.6, stat="density", size=1.3)+
    xlab("Coverage (log scale)") + 
    theme_bw() # colour=condition
  print(p0)
  print(p1)
  
  if (!is.null(genomeName)){					
    p2 <- ggplot(meth.levelsm, aes(value, colour=cpg, group=cpg)) +
      geom_line(adjust=2.5, alpha=0.6, stat="density", size=1.3)+
      xlab("Methylation Proportion") + 
      theme_bw() # colour=condition
    p3 <- ggplot(cov.matm, aes(log(value), colour=cpg, group=cpg)) +
      geom_line(adjust=2.5, alpha=0.6, stat="density", size=1.3)+
      xlab("Coverage (log scale)") + 
      theme_bw()  # colour=condition	
    
    print(p2)
    print(p3)
  }	
  
  dev.off()
  
  rm(meth.levelsm) 
  rm(rows.plot)
}

