#' BS.chr21: Whole-genome bisulfite sequencing for chromosome 21 
#' from Lister et al.
#'
#' @description This dataset represents chromosome 21 
#' from the IMR90 and H1 cell lines sequenced in Lister et al.  
#' Only CpG methylation are included. The two samples from 
#' each cell line are two different extractions (ie. technical replicates), 
#' and are pooled in the analysis in the original paper.
#' @usage data(BS.chr21)
#' @format An object of class \code{BSseq}.
#' 
#' @details All coordinates are in hg18.
#' @source Obtained from 
#' \url{http://neomorph.salk.edu/human_methylome/data.html} specifically,
#'  the files \url{mc_h1_r1.tar.gz}, \url{mc_h1_r2.tar.gz}, 
#'  \url{mc_imr90_r1.tar.gz}, \url{mc_imr90_r2.tar.gz}
#'  A script which downloads these files and constructs the \code{BS.chr21} 
#'  object may be found in \file{inst/scripts/get_BS.chr21.R} - this was 
#'  based off of and modified from the get_BS.chr22.R script in the 
#'  \code{bsseq} package. The object constructed here contains a 
#'  different chromosome (22 instead of 21), and two additional samples
#'  (h1 and imr90 instead of just imr90) to enable identification of
#'  cell type-DMRs for examples.
#' @references R Lister et al. \emph{Human DNA methylomes at base 
#'  resolution show widespread epigenomic differences}. Nature (2009) 462,
#'   315-322.
#' @examples
#'   data(BS.chr21)
#'   BS.chr21
"BS.chr21"

#' annot.chr21: Annotation information for chromosome 21, hg38 genome
#'
#' @description This is the annotation information returned from 
#' \code{\link{getAnnot}}, subsetted for chromosome 21 for convenience
#' in running the examples. The annotation is obtained using the 
#' \code{annotatr} package.
#' 
#' @usage data(annot.chr21)
#' 
#' @format a \code{GRangesList} object with two elements returned
#' by \code{\link{getAnnot}}. The first
#' contains CpG category information in the first element (optional)
#' coding gene sequence information in the second element (optional).
#' At least one of these elements needs to be non-null in order for 
#' any annotation to be plotted, but it is not necessary to contain
#' both.
#' 
#' @source Obtained from running
#' \code{annoTrack} function and then subsetting the results to 
#' only include chromosome 21. A script which executes these steps 
#' and constructs the \code{annot.chr21} 
#'  object may be found in \file{inst/scripts/get_annot.chr21.R}
#'  
#' @examples
#' data(annot.chr21)
"annot.chr21"

#' dmrs.ex: Example results of DMRs 
#'
#'@description Example output from \code{dmrseq} function run on the 
#' example dataset \code{BS.chr21}.
#' @usage data(dmrs.ex)
#' @format a data.frame that contains the results of the inference. The
#'    data.frame contains one row for each candidate region, and 
#'    10 columns, in the following order: 1. chr = 
#'    chromosome, 2. start = 
#'    start basepair position of the region, 3. end = end basepair position
#'    of the region,
#'    4. indexStart = the index of the region's first CpG, 
#'    5. indexEnd = the index of the region's last CpG,
#'    6. L = the number of CpGs contained in the region,
#'    7. area = the sum of the smoothed beta values
#'    8. beta = the coefficient value for the condition difference,
#'    9. stat = the test statistic for the condition difference,
#'    10. pval = the permutation p-value for the significance of the test
#'    statistic, and 
#'    11. qval = the q-value for the test statistic (adjustment
#'    for multiple comparisons to control false discovery rate).
#' @source Obtained from running the examples in \code{\link{dmrseq}}. 
#' A script which executes these steps 
#' and constructs the \code{dmrs.ex} 
#'  object may be found in \file{inst/scripts/get_dmrs.ex.R}
#' @examples
#' data(dmrs.ex)
"dmrs.ex"

getEstimatePooled <- function(meth.mat = meth.mat, cov.mat = cov.mat, pos = pos,
    chr = chr, design, coeff) {
    # check whether the covariate of interest is a two group comparison if not
    # (covariate is a multi-level factor or a continuous variable) then use
    # single-loci estimate of methylation effect instead of pooled mean diff
    
    if (length(unique(design[, coeff])) == 2) {
        lev1 <- 1
        lev2 <- 0
        
        # pooled
        est <- rowSums(meth.mat[, design[, coeff] == lev1])/
          rowSums(cov.mat[, design[,coeff] == lev1]) - 
          rowSums(meth.mat[, design[, coeff] == lev2])/
          rowSums(cov.mat[, design[, coeff] == lev2])
        
        sd <- cbind(rowMads(meth.mat[, design[, coeff] == lev1]/
                  cov.mat[, design[, coeff] == lev1], na.rm=TRUE)^2, 
                  rowMads(meth.mat[, design[, coeff] == lev2]/
                  cov.mat[, design[, coeff] == lev2], na.rm=TRUE)^2)
        # when sd is zero because one of the groups had only a single sample
        # with nonzero coverage, impute with the other group's sd est
        sd[rowSums(cov.mat[, design[, coeff] == lev1] > 0) == 1 &
           rowSums(cov.mat[, design[, coeff] == lev2] > 0) == 1,] <- 1
        sd[rowSums(cov.mat[, design[, coeff] == lev1] > 0) == 1,1] <- 
          sd[rowSums(cov.mat[, design[, coeff] == lev1] > 0) == 1,2]
        sd[rowSums(cov.mat[, design[, coeff] == lev2] > 0) == 1,2] <- 
          sd[rowSums(cov.mat[, design[, coeff] == lev2] > 0) == 1,1]
        sd <- 1.4826 * sqrt(rowSums(sd))
        
        return(list(rawBeta = est, sd = sd))
    } else {
        # continuous or multi-level factor case
        stop("Error: don't use pooled estimate when there ",
                    "are more than 2 groups")
    }
}

bumphunt <- function(bs, 
    design, coeff = 2, coeff.adj = 3,
    minInSpan = 30, minNumRegion = 5, 
    cutoff = NULL, maxGap = 1000, maxGapSmooth = 2500, smooth = FALSE, 
    bpSpan = 1000, verbose = TRUE, parallel = FALSE, ...) {
    
    # calculate smoothing span from minInSpan
    bpSpan2 <- NULL
    chrs <- unique(as.character(seqnames(bs)))
    for (ch in chrs) {
      pos <- chrSelectBSseq(bs, ch)
      bpSpan2 <- c(bpSpan2, minInSpan * (max(start(pos)) - 
                      min(start(pos)) + 1)/sum(seqnames(pos) == ch))
    }
    bpSpan2 <- mean(bpSpan2, na.rm = TRUE)
    
    # get 75th percentile of mean coverage genomewide before subsetting
    # by chromosome
    covQ75 <- quantile(rowMeans(as.matrix(getCoverage(bs, type="Cov"))), 0.75)
  
    tab <- NULL
    for (chromosome in chrs) {
      meth.mat <- as.matrix(getCoverage(chrSelectBSseq(bs, chromosome), 
                                        type = "M"))
      cov.mat <- as.matrix(getCoverage(chrSelectBSseq(bs, chromosome), 
                                       type = "Cov"))
      pos <- start(chrSelectBSseq(bs, chromosome))
      chr <- as.character(seqnames(chrSelectBSseq(bs, chromosome)))
    
      if (verbose) 
        message("...Chromosome ", chromosome, ": ", appendLF = FALSE)
    
      cov.means <- rowMeans(cov.mat)
    
      if (length(unique(design[, coeff])) == 2 &&
          length(coeff) == 1 &&
          all.equal(sort(unique(as.vector(design[,coeff])))[seq_len(2)],
                    c(0,1))) {
          tmp <- getEstimatePooled(meth.mat = meth.mat, 
              cov.mat = cov.mat, pos = pos, 
              chr = chr, design, coeff)
          rawBeta <- tmp$rawBeta
          sd.raw <- tmp$sd
      } else {
          tmp <- estim(meth.mat = meth.mat, cov.mat = cov.mat, 
              pos = pos, chr = chr, 
              design = design, coeff = coeff)
          rawBeta <- tmp$meth.diff
          sd.raw <- tmp$sd.meth.diff
      }
    
      sd.raw[sd.raw < 1e-05] <- 1e-05
    
      # truncate coverage at 75th percentile
      # minimum sd proportional to median coverage
      
      weights <- pmin(cov.means, covQ75) / 
        pmax(sd.raw, 1/pmax(cov.means, 5))
    
      if (smooth) {
        
        beta <- vector("list", 2)
        beta[[1]] <- beta[[2]] <- rep(NA, length(pos))
        
        beta.tmp <- smoother(y = rawBeta, 
                                 x = pos, 
                                 chr = chromosome, maxGapSmooth = maxGapSmooth,
                                 weights = weights, 
                                 minNumRegion = minNumRegion, 
                                 minInSpan = minInSpan, 
                bpSpan = bpSpan, bpSpan2 = bpSpan2, verbose = verbose, 
                parallel = parallel)
        beta[[1]] <- beta.tmp[[1]]
        beta[[2]] <- beta.tmp[[2]]

        names(beta) <- names(beta.tmp)
        Index <- which(beta$smoothed)
        beta <- beta$fitted
        
    } else {
        beta <- rawBeta
        Index <- seq(along = beta)
    }
    
    if (length(unique(as.vector(design[, coeff]))) == 2 &&
        length(coeff) == 1 &&
        all.equal(sort(unique(as.vector(design[,coeff])))[seq_len(2)],c(0,1))) {
      ngroups <- length(unique(design[,coeff]))
      rawBeta <- rawBeta/(sd.raw * ngroups / sqrt(nrow(design)))
    }else{
      rawBeta <- rawBeta/sd.raw
    }
    
    beta[-Index] <- rawBeta[-Index]
    
    tab <- rbind(tab, 
        regionScanner(meth.mat = meth.mat, cov.mat = cov.mat, pos = pos, 
        chr = chr, x = beta, y = rawBeta, maxGap = maxGap, cutoff = cutoff, 
        minNumRegion = minNumRegion, design = design, coeff = coeff, 
        coeff.adj = coeff.adj,
        verbose = verbose, parallel = parallel,
        pDat=pData(bs)))
    }
    
    if (length(tab) == 0) {
        if (verbose) 
            message("No regions found.")
        return(NULL)
    }
    return(tab)
}


refineEdges <- function(y, candidates = NULL, 
    cutoff = qt(0.975, nrow(design) - 2), verbose = FALSE, 
    minNumRegion, design) {
    
  stopifnot(length(cutoff) <= 2)
  stopifnot(is.list(candidates) && length(candidates) == 2)
    
  if (verbose) 
      message("refineEdges: refining")
  direction <- as.integer(bumphunter.greaterOrEqual(y, cutoff))
  direction[y <= -cutoff] <- -1L
    
  refineOne <- function(x, sig) {
    idx <- which(direction[x] == sig)
    if (length(idx) > 0) {
      if (length(seq(min(idx),max(idx))) >= minNumRegion) {
        x[seq(min(idx),max(idx))]
      } else {
        x
      }
    } else {
      x
    }
  }
    
  refineLong <- function(candidates, direction, minNumRegion, sig){
    if (length(candidates) > 0) {
      which.long <- which(lengths(candidates) > minNumRegion)
      if (length(which.long) > 1) {
        candidates[which.long] <- lapply(candidates[which.long], 
                                        FUN=refineOne, sig=sig)
        candidates[vapply(candidates, is.null, FUN.VALUE=logical(1))] <- NULL
      } else if (length(which.long) == 1) {
        candidates[[which.long]] <- lapply(candidates[which.long], 
                                           FUN=refineOne, sig=sig)
      }
    }
    return(candidates)
  }
    
  trimmed <- vector("list", 2)
  trimmed[[1]] <- refineLong(candidates[[1]], direction, 
                             minNumRegion, sig = 1)
  trimmed[[2]] <- refineLong(candidates[[2]], direction, 
                             minNumRegion, sig = -1)
    
  return(trimmed)
}

trimEdges <- function(x, candidates = NULL, verbose = FALSE, minNumRegion) {
    
  stopifnot(is.list(candidates) && length(candidates) == 2)
    
  if (verbose) 
      message("trimEdges: trimming")
  
  trimOne <- function(w, x, sig) {
      mid <- which.max(x[w])
      new.start <- 1
      new.end <- length(w)
      
      if (x[w[mid]]/min(x[w]) > 4/3) {
        if (w[mid] - w[1] + 1 > 4) {
          fit1 <- lm(x[w[seq_len(mid)]] ~ w[seq_len(mid)])
          if (length(summary(fit1)) > 0) {
            if (nrow(summary(fit1)$coef) == 2) {
              if (summary(fit1)$coef[2, 1] > 0 && 
                  summary(fit1)$coef[2, 4] < 0.01) {
                new.cut <- (0.5 * (x[w[mid]] - min(x[w])) + 
                              min(x[w]) +  0.75 * mean(x[w]))/2
                new.start <- min(max(1, round(mid - 0.125 * length(w))), 
                                 max(1, mid - 2), 
                                 (seq_len(mid))[min(which(x[w[seq_len(mid)]] >=
                                                            new.cut))])
              }
            }
          }
        }
        
        if (w[length(w)] - w[mid] + 1 > 4) {
          fit2 <- lm(x[w[seq(mid,length(w))]] ~ w[seq(mid,length(w))])
          if (length(fit2) > 0) {
            if (nrow(summary(fit2)$coef) == 2) {
              if (summary(fit2)$coef[2, 1] < 0 && 
                  summary(fit2)$coef[2, 4] < 0.01) {
                new.cut <- (0.5 * (x[w[mid]] - min(x[w])) + 
                              min(x[w]) + 0.75 * mean(x[w]))/2
                new.end <- max(min(round(mid + 0.125 * length(w)),
                                   length(w)), min(mid + 2, length(w)), 
                               (seq(mid,length(w)))[max(which(x[w[seq(mid,
                                                  length(w))]] >= new.cut))])
              }
            }
          }
        }
      }
      
      if (length(seq(new.start,new.end)) >= minNumRegion) {
        return(w[seq(new.start,new.end)])
      } else {
        return(w)
      }
  }

  trimLong <- function(x, candidates, minNumRegion, sig){
    if (length(candidates) > 0) {
      x <- x * sig
      which.long <- which(lengths(candidates) > minNumRegion)
      if (length(which.long) > 1) {
        candidates[which.long] <- lapply(candidates[which.long], 
                                         FUN=trimOne, x=x, sig=sig)
        candidates[vapply(candidates, is.null, FUN.VALUE=logical(1))] <- NULL
      } else if (length(which.long) == 1) {
        candidates[[which.long]] <- lapply(candidates[which.long], 
                                           FUN=trimOne, x=x, sig=sig)
      }
    }
    return(candidates)
  }
  
  trimmed <- vector("list", 2)
  trimmed[[1]] <- trimLong(x, candidates[[1]], minNumRegion, sig = 1)
  trimmed[[2]] <- trimLong(x, candidates[[2]], minNumRegion, sig = -1)
  
  return(trimmed)
}

regionScanner <- function(meth.mat = meth.mat, cov.mat = cov.mat, pos = pos, 
    chr = chr, x, y = x, ind = seq(along = x), order = TRUE, minNumRegion = 5, 
    maxGap = 300, cutoff = quantile(abs(x), 0.99), assumeSorted = FALSE, 
    verbose = verbose, design = design, coeff = coeff, coeff.adj = coeff.adj,
    parallel = parallel, pDat) {
    if (any(is.na(x[ind]))) {
        warning("NAs found and removed. ind changed.")
        ind <- intersect(which(!is.na(x)), ind)
    }
    
    cluster <- bumphunter::clusterMaker(chr, pos, maxGap = maxGap, 
                                        assumeSorted = assumeSorted)
    Indexes <- bumphunter::getSegments(x = x[ind], f = cluster[ind], 
                                       cutoff = cutoff, 
                                       assumeSorted = assumeSorted, 
                                       verbose = FALSE)
    
    # only keep up and down indices
    Indexes <- Indexes[seq_len(2)]
    
    if (sum(lengths(Indexes)) == 0) {
        message("No candidates found. ")
        return(NULL)
    }
    
    # refine edges -> start = first (stop = last) position with a raw difference
    # that meets the threshold
    Indexes <- refineEdges(y = y[ind], candidates = Indexes, cutoff = cutoff, 
        verbose = FALSE, minNumRegion = minNumRegion, design = design)
    
    if (sum(lengths(Indexes)) == 0) {
        message("No candidates found. ")
        return(NULL)
    }
    
    # refine edges II -> for larger regions, when effect size changes over the
    # region, remove portions at the beginning and end where effect size is
    # drastically different than overall effect size
    Indexes <- trimEdges(x = x[ind], candidates = Indexes, verbose = FALSE, 
                         minNumRegion = minNumRegion)
    
    if (sum(lengths(Indexes)) == 0) {
        message("No candidates found. ")
        return(NULL)
    }
    
    for (i in seq_len(2)) {
        # get number of loci in region
        lns <- lengths(Indexes[[i]])
        Indexes[[i]] <- Indexes[[i]][lns >= minNumRegion]
    }
  
    asin.gls.cov <- function(ix, design, coeff, 
        correlation = corAR1(form = ~1 |sample), 
        correlationSmall = corCAR1(form = ~L | sample), 
        weights = varPower(form = ~1/MedCov, fixed = 0.5)) {
        dat <- data.frame(g.fac = rep(pDat[,colnames(design)[coeff[1]]], 
                                      each = length(ix)),
                          sample = factor(rep(seq_len(nrow(design)), 
                                              each=length(ix))),
                          meth = as.vector(meth.mat[ix, ]),
                          cov = as.vector(cov.mat[ix, ]), 
                          L = as.vector(rep(pos[ix], nrow(design))))
        
        if(length(coeff.adj) > 0){
          dat$a.fac <- rep(pDat[,colnames(design)[coeff.adj]], 
                                  each=length(ix))
        }
        
        if(length(unique(dat$g.fac)) == 2){
          dat$g.fac <- as.factor(dat$g.fac)
        }
        
        # condition to remove regions with constant methylation / unmeth values
        if (!((length(unique(dat$meth)) == 1 && dat$meth[1] == 0) || 
              (length(unique(dat$cov - dat$meth)) == 1 && 
               (dat$cov - dat$meth)[1] == 0))) {
            
            dat$pos <- as.numeric(factor(dat$L))
            if (length(coeff.adj)==0){
              X <- model.matrix( ~ dat$g.fac + dat$L)
              mm <- formula(Z ~ g.fac + factor(L))
            }else{
              X <- model.matrix( ~ dat$g.fac + dat$L + dat$a.fac)
              mm <- formula(Z ~ g.fac + factor(L) + a.fac)
            }
            
            Y <- as.matrix(dat$meth)
            N <- as.matrix(dat$cov)
            
            dat$MedCov <- rep(as.numeric(by(dat$cov, dat$pos, median)), 
                              nrow(design))
            # remove rows with zero coverage
            whichZ <- which(dat$cov==0)
            if (length(whichZ)>0){
              dat <- dat[-whichZ,]
            }
            
            ## small constants to bound p and phi pick this such that min value
            ##  of z[m>0] is greater than all values of z[m==0]
            c0 <- 0.05
            c1 <- 0.001
            
            ## check to make sure data is complete
            ixn <- N > 0
            if (mean(ixn) < 1) {
                ## has missing entries
                X <- X[ixn, , drop = FALSE]
                Y <- Y[ixn]
                N <- N[ixn]
                ## check design not enough df for regression
                if (nrow(X) < ncol(X) + 1) {
                  message("Not enough degree of freedom to fit the ", 
                          "model. Drop some terms in formula")
                  return(data.frame(beta = NA, stat = NA, constant = FALSE))
                }
                ## design is not of full rank because of missing. Skip
                if (any(abs(svd(X)$d) < 1e-08)) 
                  return(data.frame(beta = NA, stat = NA, constant = FALSE))
            }
            
            ## Transform the methylation levels.  Add a small constant to bound 
            ## away from 0/1.
            dat$Z <- asin(2 * (Y + c0)/(N + 2 * c0) - 1)
            
            # Add a tiny amt of jitter to avoid numerically constant Z vals 
            # across a sample over the entire region
            if (max(table(dat$Z, dat$sample)) >= length(ix) - 1) {
                dat$Z <- asin(2 * (Y + c0)/(N + 2 * c0) - 1) + 
                  runif(length(Y), -c1, c1)
            }
            
            
            if (length(ix) >= 40) {
                fit <- tryCatch({
                  summary(gls(mm, weights = weights, 
                              data = dat, correlation = correlation))
                }, error = function(e) {
                  return(NA)
                })
            } else {
                # check for presence of 1-2 coverage outliers that could end up
                # driving the difference between the groups
                if (length(unique(dat$MedCov[seq_len(length(ix))])) > 1 && 
                    length(ix) <= 10) {
                  grubbs.one <- suppressWarnings(grubbs.test(
                    dat$MedCov[seq_len(length(ix))])$p.value)
                  grubbs.two <- suppressWarnings(
                    grubbs.test(dat$MedCov[seq_len(length(ix))], 
                    type = 20)$p.value)
                } else {
                  grubbs.one <- grubbs.two <- 1
                }
                
                if (grubbs.one < 0.01 || grubbs.two < 0.01) {
                  weights <- varIdent(form = ~1)
                }
                
                fit <- tryCatch({
                  summary(gls(mm, weights = weights,
                              data = dat, correlation = correlationSmall))
                }, error = function(e) {
                  return(NA)
                })
                
                # error handling in case of false convergence (don't include 
                # first variance weighting, and then corr str)
                if (sum(is.na(fit)) == length(fit)) {
                  fit <- tryCatch({
                    summary(gls(mm, data = dat, 
                                correlation = correlationSmall))
                  }, error = function(e) {
                    return(NA)
                  })
                  if (sum(is.na(fit)) == length(fit)) {
                    fit <- tryCatch({
                      summary(gls(mm, data = dat))
                    }, error = function(e) {
                      return(NA)
                    })
                  }
                }
            }
            
            if (!(sum(is.na(fit)) == length(fit))) {
              if (length(coeff)==1){
                stat <- fit$tTable[2, 3]
                beta <- fit$tTable[2, 1]
              }else{
                stat <- sqrt(anova(fit)$`F-value`[2])
                beta <- fit$tTable[2:(2+length(coeff)-1), 1]
              }
            } else {
                stat <- beta <- NA
            }
            
            df <- data.frame(stat = stat, constant=FALSE)
            nms <- gsub("g.fac", "", 
                        rownames(fit$tTable)[2:(2+length(coeff)-1)])

            if (length(beta) > 1){
              for(b in seq_len(length(beta))){
                df[[paste0("beta_", nms[b])]] <- beta[b]
              }
            }else{
              df$beta <- beta
            }
            
            df$stat = stat
            
            return(df)
        } else {
            return(data.frame(beta = NA, stat = NA, constant = TRUE))
        }
    }
    
    numCandidates <- sum(lengths(Indexes))
    
    if (numCandidates == 0) {
      message("No candidates found. ")
      return(NULL)
    }
  
    Indexes <- c(Indexes[[1]], Indexes[[2]])
    
    maxLength <- max(lengths(Indexes))
    if (maxLength > 1000) {
        message("Note: candidate regions with more than 1000 ",
                  "CpGs detected.", 
            " It is recommended to decrease the value of maxGap to ",
            "increase computational efficiency.")
    }
    
    t1 <- proc.time()
    
    if (parallel) {
        ret <- do.call("rbind", bplapply(Indexes, 
            function(Index) asin.gls.cov(ix = ind[Index], 
            design = design, coeff = coeff)))
    } else {
        ret <- do.call("rbind", lapply(Indexes, 
            function(Index) asin.gls.cov(ix = ind[Index], 
            design = design, coeff = coeff)))
    }
    
    df <- S4Vectors::DataFrame(ind, x = x[ind], chr = chr[ind], pos = pos[ind])
    res <- as.data.frame(S4Vectors::aggregate(df, S4Vectors::List(Indexes), 
                                   chr = unlist(IRanges::heads(chr, 1L)),
                                   START = min(pos), END = max(pos),
                                   indexStart = min(ind), indexEnd = max(ind),
                                   L = lengths(chr), area = abs(sum(x))))[,-1]
    colnames(res)[colnames(res)=="START"] <- "start"
    colnames(res)[colnames(res)=="END"] <- "end"
    
    res <- cbind(res, ret[,grepl("beta", colnames(ret)), drop = FALSE])
    res$stat <- ret$stat
    
    # remove regions that had constant meth values
    res <- res[!ret$constant, ]
    
    t2 <- proc.time()
    if (verbose) {
        message(numCandidates, " regions scored (", round((t2 - t1)[3]/60, 2), 
                " min). ")
    }
    
    if (order && nrow(res) > 0) 
        res <- res[order(-abs(res$stat)), ]
    
    return(res)
}

smoother <- function(y, x = NULL, weights = NULL, chr = chr, 
    maxGapSmooth = maxGapSmooth, minNumRegion = 5, verbose = TRUE, 
    minInSpan = minInSpan, bpSpan = bpSpan, bpSpan2 = bpSpan2, 
    parallel = parallel) {
    
    locfitByCluster2 <- function(ix) {
        
        ## if y is vector change to matrix
        yi <- matrix(y[ix], ncol = 1)
        xi <- x[ix]
        if (!is.null(weights)) {
            weightsi <- matrix(weights[ix], ncol = 1)
        } else {
            weightsi <- NULL
        }
        clusteri <- clusterC[ix]
        
        if (is.null((yi))) 
            stop("y (rawBeta) is missing")
        if (is.null(xi)) 
            stop("x (pos) is missing")
        if (is.null(clusteri)) 
            stop("cluster is missing")
        if (is.null(weightsi)) 
            weightsi <- matrix(1, nrow = nrow(yi), ncol = ncol(yi))
        
        Indexes <- split(seq(along = clusteri), clusteri)
        clusterL <- lengths(Indexes)
        smoothed <- rep(TRUE, nrow(yi))
        
        for (i in seq(along = Indexes)) {
            Index <- Indexes[[i]]
            if (clusterL[i] >= minNumRegion && 
                sum(rowSums(is.na(yi[Index, , drop = FALSE])) == 
                0) >= minNumRegion) {
                nn <- min(bpSpan2/(max(xi[Index]) - min(xi[Index]) + 1), 
                          minInSpan/length(Index), 0.75)
                
                for (j in seq_len(ncol(yi))) {
                  sdata <- data.frame(posi = xi[Index], yi = yi[Index, j], 
                                      weightsi = weightsi[Index, j])
                  fit <- locfit(yi ~ lp(posi, nn = nn, h = bpSpan), 
                    data = sdata, weights = weightsi, family = "gaussian", 
                    maxk = 10000)
                  pp <- preplot(fit, where = "data", band = "local", 
                    newdata = data.frame(posi = xi[Index]))
                  yi[Index, j] <- pp$trans(pp$fit)
                }
            } else {
                yi[Index, ] <- NA
                smoothed[Index] <- FALSE
            }
        }
        return(data.frame(fitted = as.vector(yi), smoothed = smoothed))
    }
    
    
    if (is.null(dim(y))) 
        y <- matrix(y, ncol = 1)  ##need to change this to be more robust
    if (!is.null(weights) && is.null(dim(weights))) 
        weights <- matrix(weights, ncol = 1)
    
    ret.all <- NULL
    
    t1 <- proc.time()
    clusterC <- bumphunter::clusterMaker(rep(chr, length(x)), x, 
                                         maxGap = maxGapSmooth)
    
    Indexes <- split(seq(along = clusterC), clusterC)
    
    idx <- NULL  ## for R CMD check
    if (parallel) {
        ret <- do.call("rbind", bplapply(Indexes, 
                                         function(idx) locfitByCluster2(idx)))
    } else {
        ret <- do.call("rbind", lapply(Indexes, 
                                       function(idx) locfitByCluster2(idx)))
    }
    
    if (verbose) {
        t2 <- proc.time()
        message("Smoothed (",
                round((t2 - t1)[3]/60, 2), " min). ",
                appendLF = FALSE)
    }
    
    return(ret)
}



## Function to compute the coefficient estimates for regression of the 
## methylation levels on the group indicator variable, at each CpG site 
## - only used if not a standard two-group comparison
getEstimate <- function(mat, design, coeff) {
    vv <- design[, coeff, drop = FALSE]
    QR <- qr(design)
    Q <- qr.Q(QR)
    R <- qr.R(QR)
    df.residual <- ncol(mat) - QR$rank
    bhat <- t(tcrossprod(backsolve(R, t(Q)), mat))
    
    if (!is.matrix(bhat)) 
      b <- matrix(b, ncol = 1)
    if (!is.matrix(mat)) 
      mat <- matrix(mat, nrow = 1)
    
    meth.diff <- rowSums(exp(bhat))/(1 + rowSums(exp(bhat))) - 
                 exp(bhat[,1])/(1 + exp(bhat[,1]))
    
    X <- rep(design, each = nrow(mat))
    X1 <- design

    res <- mat - t(tcrossprod(design, bhat))
    se2 <- rowSums(res ^ 2) / df.residual
    vb <- chol2inv(R)[coeff[1],coeff[1]] * se2
    
    out <- list(meth.diff = meth.diff, 
                sd.meth.diff = sqrt(vb), 
                df.residual = df.residual)
    
    return(out)
}

# *Experimental* - only used if not a standard two-group comparison
estim <- function(meth.mat = meth.mat, cov.mat = cov.mat, pos = pos, 
                  chr = chr, design, coeff) {
    
    ## For a linear regression of logit(meth.level) on the biological group 
    ## indicator variable at each CpG site
    mat <- meth.mat/cov.mat
    
    eps <- min(mat[mat != 0])
    mat[mat == 0] <- eps
    mat[mat == 1] <- 1 - eps
    logit.mat <- log(mat/(1 - mat))
    
    lm.fit <- getEstimate(mat = logit.mat, design = design, coeff = coeff)
    meth.diff <- lm.fit$meth.diff
    sd.meth.diff <- lm.fit$sd.meth.diff
    
    ## Returning the final list of estimated mean methylation differences and 
    ## their SDs at all CpG sites
    out <- list(meth.diff = meth.diff, sd.meth.diff = sd.meth.diff)
    return(out)
}


# pasting bumphunter's greaterOrEqual function since not exported
bumphunter.greaterOrEqual <- function(x, y) {
    precision <- sqrt(.Machine$double.eps)
    (x >= y) | (abs(x - y) <= precision)
}

