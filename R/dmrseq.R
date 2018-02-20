#' Main function for detecting and evaluating significance of DMRs.
#' 
#' Performs a two-step approach that (1) detects candidate regions, and
#' (2) scores candidate regions with an exchangeable (across the genome)
#' statistic and evaluates statistical significance using a 
#' permuation test on the pooled null distribution of scores.
#' 
#' @param bs bsseq object containing the methylation values as well as the 
#'   phenotype matrix that contains sample level covariates
#' @param testCovariate Character value or vector indicating which variables
#' (column names) in \code{pData(bs)} to test
#'  for association of methylation levels. 
#'  Can alternatively specify an integer value or vector indicating
#'  which of columns of
#'  \code{pData(bs)} to use. This is used to construct the 
#'  design matrix for the test statistic calculation.
#' @param adjustCovariate an (optional) character value or vector 
#' indicating which variables (column names) in \code{pData(bs)} 
#' will be adjusted for when 
#'  testing for the association of methylation value with the 
#'  \code{testCovariate}. 
#' Can alternatively specify an
#' integer value or vector indicating
#'  which of the columns of \code{pData(bs)} to adjust for.
#'  If not NULL (default), then this is also used to 
#'  construct the design matrix for the test statistic calculation.
#' @param matchCovariate an (optional) character value or vector 
#' indicating which variables (column names) of \code{pData(bs)} 
#' will be blocked for when 
#'  constructing the permutations in order to
#'  test for the association of methylation value with the 
#'  \code{testCovariate}. 
#'  Alternatively, you can specify an integer value or vector indicating
#'  which columns of \code{pData(bs)} to block for.
#'  Blocking means that only permutations with balanced
#'  composition of \code{testCovariate} values will be used (for example if
#'  you have samples from different gender and this is not your covariate of
#'  interest, 
#'  it is recommended to use gender as a matching covariate to avoid one 
#'  of the permutations testing entirely males versus females; this violates
#'  the null hypothesis and will decrease power).
#'  If not NULL (default), then no blocking is performed.
#' @param minInSpan positive integer that represents the minimum number of
#'    CpGs in a smoothing span window if \code{smooth} is TRUE.  
#'    Default value is 30.
#' @param minNumRegion positive integer that represents the minimum number of
#'    nucleotides to consider for a candidate region. Default value is 5.
#' @param cutoff scalar value that represents the absolute value (or a vector 
#'    of two numbers representing a lower and upper bound) for the cutoff of 
#'    the single CpG coefficient that is used to discover 
#'    candidate regions. Default value is 0.10.
#' @param smooth logical value that indicates whether or not to smooth the 
#'    CpG level signal when discovering candidate regions.
#'    Defaults to TRUE.
#' @param bpSpan a positive integer that represents the length in basepairs
#'    of the smoothing span window if \code{smooth} is TRUE.  Default value is 
#'    1000.
#' @param verbose logical value that indicates whether progress messages
#'    should be printed to stdout. Defaults value is TRUE.
#' @param BPPARAM a \code{BiocParallelParam} object to specify the parallel 
#'    backend. The default 
#'    option is \code{BiocParallel::bpparam()} which will automatically creates
#'    a cluster appropriate for the operating system. 
#' @param maxPerms a positive integer that represents the maximum number 
#'    of permutations that will be used to generate the global null 
#'    distribution of test statistics.  Default value is 20.
#' @param maxGap integer value representing maximum number of basepairs in 
#' between neighboring CpGs to be included in the same DMR. 
#' @param maxGapSmooth integer value representing maximum number of basepairs  
#' in between neighboring CpGs to be included in the same 
#' cluster when performing smoothing (should generally be larger than
#' \code{maxGap})
#' @param stat a character vector indicating the name of the column of the 
#'   output to use as the region-level test statistic. Default value is 'stat'
#'   which is the region level-statistic designed to be comparable across the
#'   genome.
#'   It is not recommended to change this argument, but it can be done for
#'   experimental purposes. Possible values are: 'L' - the number of loci
#'   in the region, 'area' - the sum of the smoothed loci statistics,
#'   'beta' - the effect size of the region, 'stat' - the test statistic for
#'   the region, or 'avg' - the average smoothed loci statistic.
#' @return a data.frame that contains the results of the inference. The
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
#' @keywords inference
#' @importFrom outliers grubbs.test
#' @importFrom bumphunter clusterMaker getSegments
#' @importFrom matrixStats colMedians
#' @importFrom matrixStats rowMads
#' 
#' @importClassesFrom bsseq BSseq 
#' @importMethodsFrom bsseq pData seqnames sampleNames start width 
#' 
#' @importFrom grDevices col2rgb colorRampPalette dev.off pdf rgb
#' @importFrom graphics axis layout legend lines mtext par 
#' plot points rect rug text
#' @importFrom methods is
#' @importFrom stats approxfun lm loess median model.matrix p.adjust
#' predict preplot qt quantile rbeta rbinom runif
#' @importFrom utils combn
#' @importFrom BiocParallel bplapply register MulticoreParam bpparam
#' @importFrom parallel mclapply
#'
#' @import bsseq 
#' @import GenomicRanges
#' @import nlme
#' @import annotatr
#' @import ggplot2
#' 
#' @export
#' 
#' @examples
#' 
#' # load example data 
#' data(BS.chr21)
#' 
#' # the covariate of interest is the 'CellType' column of pData(BS.chr21)
#' testCovariate <- 'CellType'
#' 
#' # run dmrseq on a subset of the chromosome (20K CpGs)
#' regions <- dmrseq(bs=BS.chr21[240001:260000,],
#'                  cutoff = 0.05,
#'                  testCovariate=testCovariate)
#' 
dmrseq <- function(bs, testCovariate, adjustCovariate = NULL, cutoff = 0.1, 
                   minNumRegion = 5, smooth = TRUE, bpSpan = 1000, 
                   minInSpan = 30, maxGapSmooth = 2500, maxGap = 1000, 
                   verbose = TRUE,  
                   maxPerms = 10, matchCovariate = NULL, 
                   BPPARAM = bpparam(), stat = "stat") {
    
    stopifnot(class(bs) == "BSseq")
    
    if (!(is.null(cutoff) || length(cutoff) %in% seq_len(2))) 
        stop("'cutoff' has to be either NULL or a vector of length 1 or 2")
    if (length(cutoff) == 2) 
        cutoff <- sort(cutoff)
    if (is.null(cutoff) | abs(cutoff) > 1 | abs(cutoff) == 0) 
        stop("Must specify a value for cutoff between 0 and 1")
    subverbose <- max(as.integer(verbose) - 1L, 0)
    
    # check statistic name
    if (!(stat %in% c("L", "area", "beta", "stat", "avg"))) {
        stop("Specified '", stat, 
            "' as the test statistic which is not ", 
            "in the results. Please specify a valid name from one of ",
            "L, area, beta, stat, or avg")
    }
    
    # convert covariates to column numbers if characters
    if (is.character(testCovariate)) {
        testCovariate <- which(colnames(pData(bs)) == testCovariate)
        if (length(testCovariate) == 0) {
            stop("testCovariate not found in pData(). ",
                 "Please specify a valid testCovariate")
        }
    }
    
    if (is.character(adjustCovariate)) {
        adjustCovariate <- which(colnames(pData(bs)) == adjustCovariate)
        if (length(adjustCovariate) == 0) {
            stop("adjustCovariate not found in pData(). ",
                "Please specify a valid adjustCovariate")
        }
    }
    
    
    # construct the design matrix using the pData of bs
    if (ncol(pData(bs)) < max(testCovariate, adjustCovariate)) {
        stop("Error: pData(bs) has too few columns.  ","
              Please specify valid ", 
            "covariates to use in the analysis")
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
    
    if (length(unique(testCov)) > 2) {
        message("Warning! testCovariate has more than two groups. ", 
                "Functionality is *experimental*!")
    }
    
    # check sampleSize is even in both conditions
    sampleSize <- ncol(bs)/2
    if (length(unique(table(testCov))) > 1 |
        !(sum(table(testCov) == rep(sampleSize, 
        length(unique(testCov)))) == length(unique(testCov)))) {
        stop("Error: testCov is not balanced. Need to specify an equal",
             "number of samples at each level")
    }
    
    if (!is.null(adjustCovariate)) {
        adjustCov <- pData(bs)[, adjustCovariate]
        design <- model.matrix(~testCov + adjustCov)
        colnames(design)[coeff] <- colnames(pData(bs))[testCovariate]
        colnames(design)[,seq((max(coeff) + 1),ncol(design))] <- colnames(pData(bs))[
          adjustCovariate]
    } else {
        design <- model.matrix(~testCov)
        colnames(design)[coeff] <- colnames(pData(bs))[testCovariate]
    }
    
    if (length(unique(testCov)) == 2 && 
        (is.character(pData(bs)[, testCovariate]) | 
        is.factor(pData(bs)[, testCovariate]))) {
        message("Condition ",
            unique(pData(bs)[, testCovariate][which(design[, coeff] == 1)]), 
            " vs ", 
            unique(pData(bs)[, testCovariate][which(design[, coeff] == 0)]))
    }
    if (!is.null(matchCovariate)) {
        if (is.character(matchCovariate)) {
            if (sum(grepl(matchCovariate, colnames(pData(bs)))) == 0) {
                stop("Error: no column in pData() found that matches ",
                      "the matchCovariate")
            } else if (length(grep(matchCovariate, colnames(pData(bs)))) > 1) {
                stop("Error: matchCovariate matches more than one ",
                      "column in pData()")
            }
            mC <- grep(matchCovariate, colnames(pData(bs)))
        } else {
            stopifnot(matchCovariate <= ncol(pData(bs)))
        }
    }
    
    # check for loci with missing data
    if (length(unique(testCov)) == 2){
      which.zero <- which(rowSums(as.matrix(getCoverage(
                    bs[,which(design[, coeff] == 0)], type = "Cov")) == 0) == 
                      sum(design[, coeff] == 0))
      which.zero <- unique(which.zero, which(rowSums(as.matrix(getCoverage(
                    bs[,which(design[, coeff] == 1)], type = "Cov")) == 0) ==
                      sum(design[, coeff] == 1)))
    
      if (length(which.zero) > 0) {
        stop(length(which.zero), " loci have zero coverage in all samples ",
             "of at least one condition. Please remove these loci ", 
             "before running dmrseq")
      }
    }
    
    # register the parallel backend
    BiocParallel::register(BPPARAM)
    backend <- "BiocParallel"
    
    if (bpparam()$workers == 1) {
      if (verbose) {
        mes <- "Using a single core (backend: %s).
                  "
        message(sprintf(mes, backend))
      }
      parallel <- FALSE
    } else {
      if (verbose) {
        mes <- paste0("Parallelizing using %s workers/cores ", 
                      "(backend: %s).
                        ")
        message(sprintf(mes, bpparam()$workers, backend))
      }
      parallel <- TRUE
    }
    
    meth.mat <- as.matrix(getCoverage(bs, type = "M"))
    cov.mat <- as.matrix(getCoverage(bs, type = "Cov"))
    pos <- start(bs)
    chr <- as.character(seqnames(bs))
    meta <- pData(bs)
    
    message("Detecting candidate regions with coefficient larger than ",
                   unique(abs(cutoff)), 
        " in magnitude.")
    OBS <- bumphunt(meth.mat = meth.mat, cov.mat = cov.mat, pos = pos, 
                    chr = chr, design = design, sampleSize = sampleSize, 
                    coeff = coeff, minInSpan = minInSpan, 
                    minNumRegion = minNumRegion, cutoff = cutoff, 
                    maxGap = maxGap, maxGapSmooth = maxGapSmooth, 
                    smooth = smooth, bpSpan = bpSpan, verbose = verbose, 
                    parallel = parallel)
    # check that at least one candidate region was found; if there were none 
    # there is no need to go on to compute permutation tests...
    
    if (nrow(OBS) > 0) {
        FLIP <- NULL
        # configure the permutation matrix first consider balanced, 
        # two group comparisons
        if (nrow(design)%%2 == 0 && length(unique(design[, coeff])) == 2) {
            if (verbose) {
                message("Performing balanced permutations of ",
                        "condition across samples ", 
                  "to generate a null distribution of region test statistics")
            }
            perms <- combn(seq(1, nrow(design)), sampleSize)
            perms <- perms[, seq(2,(ncol(perms)/2))]
            
            if (maxPerms < ncol(perms)) {
                # subset on 'balanced perms'
                if (sampleSize > 3 && sampleSize < 6) {
                  sg <- apply(perms, 2, function(x) sum(x > sampleSize))
                  perms <- perms[, sg < (sampleSize - 1) & sg >= 2]
                  maxPerms <- min(maxPerms, ncol(perms))
                } else if (sampleSize >= 6) {
                  sg <- apply(perms, 2, function(x) sum(x > sampleSize))
                  perms <- perms[, sg >= floor(sampleSize/2) & 
                                   sg <= ceiling(sampleSize/2)]
                }
                perms <- perms[, sort(sample(seq_len(ncol(perms)), maxPerms, 
                                             replace = FALSE))]
            }
        } else if (length(unique(design[, coeff])) > 2) {
            # Next consider a multilevel, or continuous covariate where the 
            # covariate will be permuted in an unrestricted manner
            if (verbose) {
                message("Performing unrestricted permutation of", 
                  " covariate of interest across samples ", 
                  "to generate a null distribution of region test statistics")
            }
            perms <- as.matrix(sample(seq(seq_len(nrow(design))), nrow(design)))
            
            for (p in seq_len(maxPerms - 1)) {
                tries <- 0
                candidate <- sample(seq(seq_len(nrow(design))), nrow(design))
                # check that the permutation is not a duplicate
                while (sum(apply(perms, 2, function(x) all.equal(x, 
                                                                 candidate)) == 
                  TRUE) > 0 && tries <= 20) {
                  candidate <- sample(seq(seq_len(nrow(design))), nrow(design))
                  tries <- tries + 1
                }
                # save the permutation to the permutation matrix
                perms <- cbind(perms, candidate)
            }
        } else {
            stop("Error: Currently only balanced designs ", 
                        "supported for 2-group comparisons")
        }
        
        # Now rerun on flipped designs and concatenate results
        for (j in seq_len(ncol(perms))) {
            if (verbose) {
                message("Beginning permutation ", j)
            }
            reorder <- perms[, j]
            designr <- design
            
            if (length(unique(design[, coeff])) == 2) {
                designr[, 2] <- 0
                designr[reorder, 2] <- 1
            } else {
                designr[, coeff] <- designr[reorder, coeff]
            }
            
            # if matchCovariate is not null, restrict permutations such that 
            # null comparisons are balanced for the values of 
            # pData$matchCovariate this avoids comparison of,
            # say two different individuals in the null, that the comparison of 
            # interest is tissue type. Not matching would mean the null is 
            # really not null
            if (!is.null(matchCovariate)) {
                permLabel <- paste0(paste0(meta[designr[, coeff] == 1, mC], 
                                           collapse = "_"), 
                  "vs", paste0(meta[(1 - designr[, coeff]) == 1, mC], 
                               collapse = "_"))
                
                c1 <- unlist(strsplit(permLabel, "vs"))[1]
                c2 <- unlist(strsplit(permLabel, "vs"))[2]
                
                c1 <- unlist(strsplit(c1, "_"))
                c2 <- unlist(strsplit(c2, "_"))
                
                keepPerm <- 1 * (sum(c1 %in% c2) == length(c1))
                
                if (keepPerm == 0) {
                  if (verbose) {
                    message(paste0("Skipping permutation ", permLabel))
                  }
                  next
                }
            } else {
                permLabel <- j
            }
            
            res.flip.p <- bumphunt(meth.mat = meth.mat, cov.mat = cov.mat, 
                                   pos = pos, chr = chr, design = designr, 
                                   sampleSize = sampleSize, coeff = coeff, 
                                   minInSpan = minInSpan, 
                                   minNumRegion = minNumRegion, cutoff = cutoff,
                                   maxGap = maxGap, maxGapSmooth = maxGapSmooth,
                                   smooth = smooth, bpSpan = bpSpan, 
                                   verbose = verbose, parallel = parallel)
            
            if (!is.null(res.flip.p)) {
                res.flip.p$permNum <- permLabel
                FLIP <- rbind(FLIP, res.flip.p)
            }
            
            if (verbose) {
                message("* ", j, " out of ", ncol(perms), 
                     " permutations completed
                     ")
            }
        }
        
        # if there are more than 1 million candidate null regions, 
        # take a random sample
        # of 1 million of them
        if (nrow(FLIP) > 1e+06) {
            rs <- sample(seq_len(nrow(FLIP)), 1e+06, replace = FALSE)
            FLIP <- FLIP[rs, ]
        }
        
        # which column of results to use as test statistic ?  
        # check statistic name
        if (!(stat %in% c(colnames(OBS), "avg"))) {
            stop("Specified '", stat, 
                        "' as the test statistic which is not ", 
                "in the results. Please specify a valid name from one of ",
                "L, area, beta, or stat")
        } else if (stat == "avg") {
            OBS$avg <- OBS$area/OBS$L
            FLIP$avg <- FLIP$area/FLIP$L
        }
        
        whichStatO <- which(colnames(OBS) == stat)
        whichStatF <- which(colnames(FLIP) == stat)
        
        # Faster way to compute the p-values that doesn't use multiple cores 
        # Step 1: sort the permuted statistics vector
        perm.ordered <- c(sort(abs(FLIP[, whichStatF]), method = "quick"), Inf)
        
        # Step 2: find the first instance in the sorted vector where the 
        # permuted value is greater than the observed and use this to determine 
        # the number of permuted values that are greater than or equal to the 
        # observed
        pval <- rep(NA, nrow(OBS))
        pval[!is.na(OBS[, whichStatO])] <- (1 + 
                    vapply(abs(OBS[!is.na(OBS[, whichStatO]), whichStatO]),
          function(x) length(perm.ordered) - min(which(x <= perm.ordered)),
          numeric(1))) /
          (1 + sum(!is.na(FLIP[, whichStatF])))
        
        # missing test statistics cause Inf for the p-value calculation instead,
        # propagate the missing values
        pval[abs(pval) == Inf] <- NA
        
        pval <- data.frame(x = pval, y = p.adjust(pval, method = "BH"))
        
        OBS$pval <- pval$x
        OBS$qval <- pval$y
        
        return(OBS)
    } else {
        message("No candidate regions pass the cutoff of ", unique(abs(cutoff)))
    }
}
