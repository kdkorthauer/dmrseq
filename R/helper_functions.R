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
        stop(paste0("Error: don't use pooled estimate when there ",
                    "are more than 2 groups"))
    }
}

bumphunt <- function(meth.mat = meth.mat, cov.mat = cov.mat, pos = pos, 
    chr = chr, design, sampleSize, coeff = 2, minInSpan = 30, minNumRegion = 5, 
    cutoff = NULL, maxGap = 1000, maxGapSmooth = 2500, smooth = FALSE, 
    bpSpan = 1000, verbose = TRUE, parallel = FALSE, ...) {
    
    # calculate smoothing span from minInSpan
    bpSpan2 <- NULL
    for (ch in unique(chr)) {
        bpSpan2 <- c(bpSpan2, minInSpan * (max(pos[chr == ch]) - 
                              min(pos[chr == ch]) + 1)/sum(chr == ch))
    }
    bpSpan2 <- mean(bpSpan2, na.rm = TRUE)
    
    if (verbose) 
        message(".....Computing coefficients.")
    
    cov.means <- rowMeans(cov.mat)
    
    if (length(unique(design[, coeff])) == 2) {
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
 
    weights <- pmin(cov.means, quantile(cov.means, 0.75)) / 
      pmax(sd.raw, 1/pmax(cov.means, 5))
    
    if (smooth) {
        if (verbose) 
            message(".....Smoothing coefficients.")
        
        beta <- vector("list", 2)
        beta[[1]] <- beta[[2]] <- rep(NA, length(pos))
        for (chromosome in unique(chr)) {
            beta.tmp <- smoother(y = rawBeta[chr == chromosome], 
                                 x = pos[chr == chromosome], 
                                 chr = chromosome, maxGapSmooth = maxGapSmooth,
                                 weights = weights[chr == chromosome], 
                                 minNumRegion = minNumRegion, 
                                 minInSpan = minInSpan, 
                bpSpan = bpSpan, bpSpan2 = bpSpan2, verbose = verbose, 
                parallel = parallel)
            beta[[1]][chr == chromosome] <- beta.tmp[[1]]
            beta[[2]][chr == chromosome] <- beta.tmp[[2]]
        }
        
        names(beta) <- names(beta.tmp)
        Index <- which(beta$smoothed)
        beta <- beta$fitted
        
    } else {
        beta <- rawBeta
        Index <- seq(along = beta)
    }
    
    rawBeta <- rawBeta/(sd.raw * sqrt(2/sampleSize))
    
    beta[-Index] <- rawBeta[-Index]
    
    tab <- regionScanner(meth.mat = meth.mat, cov.mat = cov.mat, pos = pos, 
        chr = chr, x = beta, y = rawBeta, maxGap = maxGap, cutoff = cutoff, 
        minNumRegion = minNumRegion, design = design, coeff = coeff, 
        verbose = verbose, sampleSize = sampleSize, parallel = parallel)

    if (nrow(tab) == 0) {
        if (verbose) 
            message("No regions found.")
        return(NA)
    }
    return(table = tab)
}


refineEdges <- function(y, candidates = NULL, 
    cutoff = qt(0.975, 2 * sampleSize - 2), verbose = FALSE, 
    minNumRegion, sampleSize) {
    stopifnot(length(cutoff) <= 2)
    stopifnot(is.list(candidates) & length(candidates) == 2)
    
    if (verbose) 
        message("refineEdges: refining")
    direction <- as.integer(bumphunter.greaterOrEqual(y, cutoff))
    direction[y <= -cutoff] <- -1L
    
    trimmed <- candidates
    for (ii in 1:2) {
        if (ii == 1) {
            sig <- 1
        } else {
            sig <- -1
        }
        if (length(candidates[[ii]]) > 0) {
            which.long <- which(sapply(candidates[[ii]], length) > minNumRegion)
            if (length(which.long) > 1) {
                trimmed[[ii]][which.long] <-sapply(candidates[[ii]][which.long],
                  function(x) {
                    idx <- which(direction[x] == sig)
                    if (length(idx) > 0) {
                      if (length(min(idx):max(idx)) >= minNumRegion) {
                        x[min(idx):max(idx)]
                      } else {
                        x
                      }
                    } else {
                      x
                    }
                  })
                trimmed[[ii]][sapply(trimmed[[ii]], is.null)] <- NULL
            } else if (length(which.long) == 1) {
                trimmed[[ii]][[which.long]] <- sapply(
                  candidates[[ii]][which.long], 
                  function(x) {
                    idx <- which(direction[x] == sig)
                    if (length(idx) > 0) {
                      if (length(min(idx):max(idx)) >= minNumRegion) {
                        x[min(idx):max(idx)]
                      } else {
                        x
                      }
                    } else {
                      x
                    }
                  })
            }
        }
    }
    
    return(trimmed)
}

trimEdges <- function(x, candidates = NULL, verbose = FALSE, minNumRegion) {
    
    stopifnot(is.list(candidates) & length(candidates) == 2)
    
    if (verbose) 
        message("trimEdges: trimming")
    
    trimmed <- candidates
    for (ii in 1:2) {
        if (length(candidates[[ii]]) > 0) {
            if (ii == 1) {
                sig <- 1
            } else {
                sig <- -1
            }
            x <- x * sig
            which.long2 <- which(sapply(candidates[[ii]], length) > 
                                   minNumRegion)
            if (length(which.long2) > 0) {
                repl <- sapply(candidates[[ii]][which.long2], function(w) {
                  mid <- which.max(x[w])
                  new.start <- 1
                  new.end <- length(w)
                  
                  if (x[w[mid]]/min(x[w]) > 4/3) {
                    if (w[mid] - w[1] + 1 > 4) {
                      fit1 <- lm(x[w[1:mid]] ~ w[1:mid])
                      if (length(summary(fit1)) > 0) {
                        if (nrow(summary(fit1)$coef) == 2) {
                          if (summary(fit1)$coef[2, 1] > 0 & 
                              summary(fit1)$coef[2, 
                            4] < 0.01) {
                            new.cut <- (0.5 * (x[w[mid]] - min(x[w])) + 
                                          min(x[w]) +  0.75 * mean(x[w]))/2
                            new.start <- min(max(1, 
                                          round(mid - 0.125 * length(w))), 
                              max(1, mid - 2), 
                              (1:mid)[min(which(x[w[1:mid]] >= new.cut))])
                            
                          }
                        }
                      }
                    }
                    
                    if (w[length(w)] - w[mid] + 1 > 4) {
                      fit2 <- lm(x[w[mid:length(w)]] ~ w[mid:length(w)])
                      if (length(fit2) > 0) {
                        if (nrow(summary(fit2)$coef) == 2) {
                          if (summary(fit2)$coef[2, 1] < 0 & 
                              summary(fit2)$coef[2, 
                            4] < 0.01) {
                            new.cut <- (0.5 * (x[w[mid]] - min(x[w])) + 
                                          min(x[w]) + 0.75 * mean(x[w]))/2
                            new.end <- max(min(round(mid + 0.125 * length(w)),
                                      length(w)), min(mid + 2, length(w)), 
                                (mid:length(w))[max(which(x[w[mid:length(w)]] >=
                                new.cut))])
                          }
                        }
                      }
                    }
                  }
                  if (length(new.start:new.end) >= minNumRegion) {
                    return(w[new.start:new.end])
                  } else {
                    return(w)
                  }
                })
                if (length(which.long2) > 1) {
                  trimmed[[ii]][which.long2] <- repl
                } else {
                  trimmed[[ii]][[which.long2]] <- repl
                }
            }
        }
    }
    return(trimmed)
}

# function to compute raw mean methylation differences
meanDiff <- function(bs, dmrs, testCovariate, adjustCovariate) {
    # convert covariates to column numbers if characters
    if (is.character(testCovariate)) {
        testCovariate <- which(colnames(pData(bs)) == testCovariate)
        if (length(testCovariate) == 0) {
            stop(paste0("testCovariate not found in pData(). ", 
                        "Please specify a valid testCovariate"))
        }
    }
    
    if (ncol(pData(bs)) < max(testCovariate, adjustCovariate)) {
        stop(paste0("Error: pData(bs) has too few columns.  ", 
                    "Please specify valid ", 
            "covariates to use in the analysis"))
    }
    
    coeff <- 2:(2 + length(testCovariate) - 1)
    testCov <- pData(bs)[, testCovariate]
    if (length(unique(testCov)) == 1) {
        message(paste0("Warning: only one unique value of the specified ", 
                       "covariate of interest.  Assuming null comparison and ", 
            "splitting sample group into two equal groups"))
        testCov <- rep(1, length(testCov))
        testCov[1:round(length(testCov)/2)] <- 0
    }
    if (!is.null(adjustCovariate)) {
        adjustCov <- pData(bs)[, adjustCovariate]
        design <- model.matrix(~testCov + adjustCov)
        colnames(design)[coeff] <- colnames(pData(bs))[testCovariate]
        colnames(design)[, (max(coeff)+1):ncol(design)] <- colnames(pData(bs))[
          adjustCovariate]
    } else {
        design <- model.matrix(~testCov)
        colnames(design)[coeff] <- colnames(pData(bs))[testCovariate]
    }
    
    if (length(unique(design[, coeff])) != 2) {
        message(paste0("Not a two-group comparison. Can't compute simple mean ",
            "methylation differences. ", "Returning beta estimates instead"))
        return(dmrs$beta)
    } else {
        prop.mat <- as.matrix(getCoverage(bs, type = "M")/getCoverage(bs, type = "Cov"))
        levs <- unique(design[, coeff])
        
        meanDiff <- sapply(1:nrow(dmrs), function(x) {
            return(mean(rowMeans(prop.mat[(dmrs$indexStart[x]:dmrs$indexEnd[x]),
                which(design[, coeff] == levs[2])]) - rowMeans(
                  prop.mat[(dmrs$indexStart[x]:dmrs$indexEnd[x]), 
                which(design[, coeff] == levs[1])]), na.rm = TRUE))
        })
        return(meanDiff)
    }
}

regionScanner <- function(meth.mat = meth.mat, cov.mat = cov.mat, pos = pos, 
    chr = chr, x, y = x, ind = seq(along = x), order = TRUE, minNumRegion = 5, 
    maxGap = 300, cutoff = quantile(abs(x), 0.99), assumeSorted = FALSE, 
    verbose = verbose, design = design, coeff = coeff, sampleSize = sampleSize, 
    parallel = parallel) {
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
    Indexes <- Indexes[1:2]
    
    if (sum(sapply(Indexes, length)) == 0) {
        return(NULL)
    }
    
    # refine edges -> start = first (stop = last) position with a raw difference
    # that meets the threshold
    Indexes <- refineEdges(y = y[ind], candidates = Indexes, cutoff = cutoff, 
        verbose = FALSE, minNumRegion = minNumRegion, sampleSize = sampleSize)
    
    if (sum(sapply(Indexes, length)) == 0) {
        return(NULL)
    }
    
    # refine edges II -> for larger regions, when effect size changes over the
    # region, remove portions at the beginning and end where effect size is
    # drastically different than overall effect size
    Indexes <- trimEdges(x = x[ind], candidates = Indexes, verbose = FALSE, 
                         minNumRegion = minNumRegion)
    
    if (sum(sapply(Indexes, length)) == 0) {
        return(NULL)
    }
    
    for (i in 1:2) {
        # get number of loci in region
        lns <- sapply(Indexes[[i]], length)
        Indexes[[i]] <- Indexes[[i]][lns >= minNumRegion]
    }
    
    asin.gls.cov <- function(ix, design, coeff, correlation = corAR1(form = ~1 |
        sample), correlationSmall = corCAR1(form = ~L | sample), 
        weights = varPower(form = ~1/MedCov, 
        fixed = 0.5)) {
        sampleSize <- nrow(design)/2
        dat <- data.frame(g.fac = factor(as.vector(sapply(design[, coeff], 
                                                          function(x) rep(x, 
            length(ix))))), sample = factor(as.vector(sapply(1:(sampleSize * 2),
            function(x) rep(x, length(ix))))), 
            meth = melt(meth.mat[ix, ])$value,
            cov = melt(cov.mat[ix, ])$value, L = as.vector(rep(pos[ix], 
                                                               nrow(design))))
        
        # condition to remove regions with constant methylation / unmeth values
        if (!((length(unique(dat$meth)) == 1 & dat$meth[1] == 0) | 
              (length(unique(dat$cov - 
            dat$meth)) == 1 & (dat$cov - dat$meth)[1] == 0))) {
            
            dat$pos <- as.numeric(factor(dat$L))
            X <- model.matrix(~dat$g.fac)
            colnames(X)[2] <- "grp"
            
            if (nrow(X) <= ncol(X)){
                message("Not enough degree of freedom to fit the linear ", 
                            "model. Drop some terms in formula")
                return(data.frame(beta = NA, stat = NA, constant = FALSE))
            }
            
            Y <- as.matrix(dat$meth)
            N <- as.matrix(dat$cov)
            
            dat$MedCov <- rep(as.numeric(by(dat$cov, dat$pos, median)), 
                              sampleSize * 2)
            # remove rows with zero coverage
            whichZ <- which(dat$cov==0)
            if (length(whichZ)>0){
              dat <- dat[-whichZ,]
            }
            
            # tol is the convergence criteria (for difference Between two 
            # iterations of the dispersion parameter estimate
            
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
                if (nrow(X) < ncol(X) + 1) 
                  return(data.frame(beta = NA, stat = NA, constant = FALSE))
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
                  summary(gls(Z ~ g.fac + factor(L), weights = weights, 
                              data = dat, correlation = correlation))
                }, error = function(e) {
                  return(NA)
                })
            } else {
                # check for presence of 1-2 coverage outliers that could end up
                # driving the difference between the groups
                if (length(unique(dat$MedCov[1:length(ix)])) > 1 & length(ix) <=
                  10) {
                  grubbs.one <- suppressWarnings(grubbs.test(
                    dat$MedCov[1:length(ix)])$p.value)
                  grubbs.two <- suppressWarnings(
                    grubbs.test(dat$MedCov[1:length(ix)], 
                    type = 20)$p.value)
                } else {
                  grubbs.one <- grubbs.two <- 1
                }
                
                if (grubbs.one < 0.01 | grubbs.two < 0.01) {
                  weights <- varIdent(form = ~1)
                }
                
                fit <- tryCatch({
                  summary(gls(Z ~ g.fac + factor(L), weights = weights,
                              data = dat, correlation = correlationSmall))
                }, error = function(e) {
                  return(NA)
                })
                
                # error handling in case of false convergence (don't include 
                # first variance weighting, and then corr str)
                if (sum(is.na(fit)) == length(fit)) {
                  fit <- tryCatch({
                    summary(gls(Z ~ g.fac + factor(L), data = dat, 
                                correlation = correlationSmall))
                  }, error = function(e) {
                    return(NA)
                  })
                  if (sum(is.na(fit)) == length(fit)) {
                    fit <- tryCatch({
                      summary(gls(Z ~ g.fac + factor(L), data = dat))
                    }, error = function(e) {
                      return(NA)
                    })
                  }
                }
            }
            
            if (!(sum(is.na(fit)) == length(fit))) {
                stat <- fit$tTable[2, 3]
                beta <- fit$tTable[2, 1]
            } else {
                stat <- beta <- NA
            }
            
            return(data.frame(beta = beta, stat = stat, constant = FALSE))
        } else {
            return(data.frame(beta = NA, stat = NA, constant = TRUE))
        }
    }
    
    numCandidates <- sum(sapply(Indexes, length))
    
    if (verbose) {
        message(paste0(".....Evaluating ", numCandidates, 
                       " candidate regions."))
    }
    
    if (numCandidates == 0) {
        return(NULL)
    }
    
    Indexes <- c(Indexes[[1]], Indexes[[2]])
    
    maxLength <- max(sapply(Indexes, length))
    if (maxLength > 1000) {
        message(paste0("Note: candidate regions with more than 1000 ",
                       "CpGs detected.", 
            " It is recommended to decrease the value of maxGap to ",
            "increase computational efficiency."))
    }
    
    t1 <- proc.time()
    res <- data.frame(chr = sapply(Indexes, function(Index) chr[ind[Index[1]]]),
        start = sapply(Indexes, function(Index) min(pos[ind[Index]])),
        end = sapply(Indexes, 
            function(Index) max(pos[ind[Index]])), 
        indexStart = sapply(Indexes, function(Index) min(ind[Index])), 
        indexEnd = sapply(Indexes, function(Index) max(ind[Index])),
        L = sapply(Indexes, 
            length),
        area = sapply(Indexes, function(Index) abs(sum(x[ind[Index]]))), 
        stringsAsFactors = FALSE)
    
    
    if (parallel) {
        ret <- do.call("rbind", bplapply(Indexes, 
            function(Index) asin.gls.cov(ix = ind[Index], 
            design = design, coeff = coeff)))
    } else {
        ret <- do.call("rbind", lapply(Indexes, 
            function(Index) asin.gls.cov(ix = ind[Index], 
            design = design, coeff = coeff)))
    }
    res$beta <- ret$beta
    res$stat <- ret$stat
    
    # remove regions that had constant meth values
    res <- res[!ret$constant, ]
    
    t2 <- proc.time()
    if (verbose) {
        message(paste0(".....Took ", round((t2 - t1)[3]/60, 2), 
                       " min to score candidate regions."))
    }
    
    if (order & nrow(res) > 0) 
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
        clusterL <- sapply(Indexes, length)
        smoothed <- rep(TRUE, nrow(yi))
        
        for (i in seq(along = Indexes)) {
            Index <- Indexes[[i]]
            if (clusterL[i] >= minNumRegion & 
                sum(rowSums(is.na(yi[Index, , drop = FALSE])) == 
                0) >= minNumRegion) {
                nn <- min(bpSpan2/(max(xi[Index]) - min(xi[Index]) + 1), 
                          minInSpan/length(Index), 0.75)
                
                for (j in 1:ncol(yi)) {
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
        ret <- do.call("rbind", mclapply(Indexes, function(idx) {
            sm <- locfitByCluster2(idx)
            return(sm)
        }, mc.cores = bpparam()$workers))
    } else {
        ret <- do.call("rbind", lapply(Indexes, function(idx) {
            sm <- locfitByCluster2(idx)
            return(sm)
        }))
    }
    
    if (verbose) {
        t2 <- proc.time()
        message(paste0("..........Done Smoothing Chromosome ", chr, ". Took ",
                       round((t2 - 
            t1)[3]/60, 2), " minutes"))
    }
    
    return(ret)
}



## Function to compute the coefficient estimates for regression of the 
## methylation levels on the group indicator variable, at each CpG site 
## *Experimental* - only used if not a standard two-group comparison
getEstimate <- function(mat, design, coeff) {
    v <- design[, coeff]
    A <- design[, -coeff, drop = FALSE]
    qa <- qr(A)
    S <- diag(nrow(A)) - tcrossprod(qr.Q(qa))
    vv <- matrix(v, ncol = 1)
    
    sv <- S %*% vv
    vsv <- diag(crossprod(vv, sv))
    b <- (mat %*% crossprod(S, vv)/vsv)
    a <- mat %*% as.matrix(rep(1/nrow(design), nrow(design))) - b * mean(v)
    
    x.cov <- design[, 2]
    deno <- sum((x.cov - mean(x.cov))^2)
    vcov.mat <- matrix(c(mean(x.cov^2)/deno, rep(-mean(x.cov)/deno, 2), 1/deno),
        2, 2)
    A1 <- exp(a + b)/(1 + exp(a + b))^2 - exp(a)/(1 + exp(a))^2
    A2 <- exp(a + b)/(1 + exp(a + b))^2
    var.coeff <- A1^2 * vcov.mat[1, 1] + 2 * A1 * A2 * vcov.mat[1, 2] + A2^2 * 
      vcov.mat[2, 2]
    meth.diff <- exp(a + b)/(1 + exp(a + b)) - exp(a)/(1 + exp(a))
    
    if (!is.matrix(b)) 
        b <- matrix(b, ncol = 1)
    if (!is.matrix(mat)) 
        mat <- matrix(mat, nrow = 1)
    
    sy <- mat %*% S
    df.residual <- ncol(mat) - qa$rank - 1
    sigma <- matrix(sqrt(rowSums((sy - tcrossprod(b, sv))^2)/df.residual), 
                    ncol = 1)
    sd.meth.diff <- sqrt(var.coeff) * sigma
    out <- list(meth.diff = meth.diff, sd.meth.diff = sd.meth.diff, 
                stdev.unscaled = sqrt(1/vsv), 
        df.residual = df.residual)
    out$stdev <- as.numeric(out$stdev)
    
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

# pasting bumphunter's reduceIt function since not exported
bumphunter.reduceIt <- function(x, elem, bind = rbind) {
    if (missing(elem)) 
        elem <- names(x[[1]])
    ret <- lapply(elem, function(el) {
        xx <- lapply(x, "[[", el)
        if (is.matrix(xx[[1]])){
            return(do.call(bind, xx)) 
        }else if (is.vector(xx[[1]])){
            return(do.call("c", xx))
        }else{
          stop("reduce can only handle matrices or vectors")
        }
    })
    names(ret) <- elem
    ret
}

# pasting bumphunter's greaterOrEqual function since not exported
bumphunter.greaterOrEqual <- function(x, y) {
    precision <- sqrt(.Machine$double.eps)
    (x >= y) | (abs(x - y) <= precision)
}

