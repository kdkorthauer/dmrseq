#' Main function for detecting and evaluating significance of DMRs.
#' 
#' Performs a two-step approach that (1) detects candidate regions, and
#' (2) scores candidate regions with an exchangeable (across the genome)
#' statistic and evaluates statistical significance using a 
#' permuation test on the pooled null distribution of scores.
#' 
#' @param bs bsseq object containing the methylation values as well as the 
#'   phenotype matrix that contains sample level covariates
#' @param testCovariate Character value indicating which variable
#'  (column name) in \code{pData(bs)} to test
#'  for association of methylation levels. 
#'  Can alternatively specify an integer value indicating
#'  which of column of
#'  \code{pData(bs)} to use. This is used to construct the 
#'  design matrix for the test statistic calculation. To run using a 
#'  continuous or categorial covariate with more than two groups, simply pass in
#'  the name of a column in `pData` that contains this covariate. A continuous
#'  covariate is assmued if the data type in the `testCovariate` slot is 
#'  continuous, with the exception of if there are only two unique values 
#'  (then a two group comparison is carried out).
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
#' @param matchCovariate An (optional) character value 
#'  indicating which variable (column name) of \code{pData(bs)} 
#'  will be blocked for when 
#'  constructing the permutations in order to
#'  test for the association of methylation value with the 
#'  \code{testCovariate}, only to be used when \code{testCovariate}
#'  is a two-group factor and the number of permutations possible is less
#'  than 500000.
#'  Alternatively, you can specify an integer value indicating
#'  which column of \code{pData(bs)} to block for.
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
#'    CpGs to consider for a candidate region. Default value is 5. 
#'    Minimum value is 3.
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
#'    distribution of test statistics.  Default value is 10.
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
#' @param block logical indicating whether to search for large-scale (low
#'  resolution) blocks of differential methylation (default is FALSE, which
#'  means that local DMRs are desired). If TRUE, the parameters for 
#'  \code{bpSpan}, \code{minInSpan}, and \code{maxGapSmooth} should be adjusted
#'  (increased) accordingly. This setting will also merge
#'  candidate regions that (1) are in the same direction and (2) are less than 
#'  1kb apart with no covered CpGs separating them. The region-level model used 
#'  is also slightly modified - instead of a loci-specific intercept for each 
#'  CpG in theregion, the intercept term is modeled as a natural spline with 
#'  one interior knot per each 10kb of length (up to 10 interior knots).
#' @param blockSize numeric value indicating the minimum number of basepairs 
#'  to be considered a block (only used if \code{block}=TRUE). Default is 
#'  5000 basepairs.
#' @param chrsPerChunk a positive integer value indicating the number of 
#'  chromosomes per chunk. The default is 1, meaning that the data will be 
#'  looped through one chromosome at a time. When pairing up multiple 
#'  chromosomes per chunk, sizes (in terms of numbers of CpGs) will be taken
#'  into consideration to balance the sizes of each chunk.
#' @return a \code{GRanges} object that contains the results of the inference. 
#'    The object contains one row for each candidate region, sorted by q-value
#'    and then chromosome. The standard 
#'    \code{GRanges} chr, start, and end are included, along with at least
#'    7 metadata
#'    columns, in the following order: 
#'    1. L = the number of CpGs contained in the region,
#'    2. area = the sum of the smoothed beta values
#'    3. beta = the coefficient value for the condition difference (there 
#'       will be more than one column here if a multi-group comparison
#'       was performed),
#'    4. stat = the test statistic for the condition difference,
#'    5. pval = the permutation p-value for the significance of the test
#'    statistic, and 
#'    6. qval = the q-value for the test statistic (adjustment
#'    for multiple comparisons to control false discovery rate).
#'    7. index = an \code{IRanges} containing the indices of the region's 
#'       first CpG to last CpG.
#'
#' @keywords inference
#' @importFrom outliers grubbs.test
#' @importFrom bumphunter clusterMaker getSegments
#' @importFrom DelayedMatrixStats colMedians rowMads rowSums2 rowMeans2 rowDiffs
#' @importFrom matrixStats rowRanges
#' @importFrom stats formula anova as.formula
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
#' @importFrom splines ns
#'
#' @import bsseq 
#' @import GenomicRanges
#' @import nlme
#' @import annotatr
#' @import ggplot2
#' @import S4Vectors
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
#' # run dmrseq on a subset of the chromosome (10K CpGs)
#' regions <- dmrseq(bs=BS.chr21[240001:250000,],
#'                  cutoff = 0.05,
#'                  testCovariate=testCovariate)
#' 
dmrseq <- function(bs, testCovariate, adjustCovariate = NULL, cutoff = 0.1, 
                   minNumRegion = 5, smooth = TRUE, bpSpan = 1000, 
                   minInSpan = 30, maxGapSmooth = 2500, maxGap = 1000, 
                   verbose = TRUE,  
                   maxPerms = 10, matchCovariate = NULL, 
                   BPPARAM = bpparam(), stat = "stat", 
                   block = FALSE, blockSize = 5000,
                   chrsPerChunk = 1) {
    
    stopifnot(is(bs, "BSseq"))
    
    if (!(is.null(cutoff) || length(cutoff) %in% seq_len(2))) 
        stop("'cutoff' has to be either NULL or a vector of length 1 or 2")
    if (length(cutoff) == 2) 
        cutoff <- sort(cutoff)
    if (is.null(cutoff) | abs(cutoff) > 1 | abs(cutoff) == 0) 
        stop("Must specify a value for cutoff between 0 and 1")
    subverbose <- max(as.integer(verbose) - 1L, 0)
    
    if(minNumRegion < 3){
      stop("minNumRegion must be at least 3")
    }
    
    # check statistic name
    if (!(stat %in% c("L", "area", "beta", "stat", "avg"))) {
        stop("Specified '", stat, 
            "' as the test statistic which is not ", 
            "in the results. Please specify a valid name from one of ",
            "L, area, beta, stat, or avg")
    }
    
    # informative message about blocks if block=TRUE; check for increased
    # smoothing window
    if (block){
      message("Searching for large scale blocks with at least ",
              blockSize, " basepairs.")
      
      if(minInSpan < 100 && bpSpan < 2000 && maxGapSmooth < 1e5){
        warning("When block=TRUE, it is recommended to increase the values ",
                "of minInSpan, bpSpan, and maxGapSmooth in order to widen ",
                "the smoothing window")
      }
    }
    
    # convert covariates to column numbers if characters
    if (is.character(testCovariate)) {
        if(length(testCovariate) > 1)
          stop("Only one testCovariate can be specified")
        if(is.character(adjustCovariate)){
          if(sum(testCovariate %in% adjustCovariate) > 0)
            stop("adjustCovariate can't contain testCovariate")
        }
        if(is.character(matchCovariate)){
          if(sum(testCovariate %in% matchCovariate))
            stop("matchCovariate can't contain testCovariate")
        }
        testCovariate <- which(colnames(pData(bs)) == testCovariate)
        if (length(testCovariate) == 0) {
            stop("testCovariate not found in pData(). ",
                 "Please specify a valid testCovariate")
        }
    }
    
    if (is.character(adjustCovariate)) {
        if(is.character(matchCovariate)){
          if(matchCovariate == adjustCovariate)
            stop("matchCovariate can't be identical to adjustCovariate")
        }
        adjustCovariate <- which(colnames(pData(bs)) %in% adjustCovariate)
        if (length(adjustCovariate) == 0) {
            stop("adjustCovariate not found in pData(). ",
                "Please specify a valid adjustCovariate")
        }
    }
    
    # check that chrsPerChunk value makes sense
    if (chrsPerChunk != 1){
      if (chrsPerChunk%%1 != 0){
        stop("chrsPerChunk must be an integer")
      }else if(chrsPerChunk < 1){
        stop("chrsPerChunk must be strictly positive")
      }else if(chrsPerChunk > length(unique(seqnames(bs)))){
        stop("chrsPerChunk can't be larger than the total",
             " number of chromosomes")
      }
    }
    
    # check that bs object is sorted since `bsseq::BSseq()` no longer 
    # automatically sorts to ensure loci from same chr are indexed consecutively
    if (is.unsorted(bs)) {
      stop("'bs' must be sorted. Use 'sort(bs)'.")
    }
    
    # construct the design matrix using the pData of bs
    if (ncol(pData(bs)) < max(testCovariate, adjustCovariate)) {
        stop("Error: pData(bs) has too few columns.  ","
              Please specify valid ", 
            "covariates to use in the analysis")
    }
    
    coeff <- seq(2,(2 + length(testCovariate) - 1))
    testCov <- pData(bs)[, testCovariate]
    if (is.factor(testCov)) # drop unused levels of test
      testCov <- droplevels(testCov)
    fact <- TRUE
    sampleSize <- table(testCov)[names(table(testCov)) %in% pData(bs)[,testCovariate]]
    if (length(unique(testCov)) == 1) {
        message("Warning: only one unique value of the specified ", 
                "covariate of interest.  Assuming null comparison and ",
                "splitting sample group into two equal groups")
        testCov <- rep(1, length(testCov))
        testCov[seq_len(round(length(testCov)/2))] <- 0
    }else if (length(unique(testCov)) > 2 && !is.numeric(testCov)) {
        message("Performing a global test of H0: no difference among ",
                length(unique(testCov)), " groups (assuming the test ",
                "covariate ", colnames(pData(bs))[testCovariate],
                " is a factor).")
        coeff <- seq(coeff, coeff + length(unique(testCov)) - 2)
    }else if (length(unique(testCov)) > 2 && is.numeric(testCov)) {
        message("Assuming the test ",
              "covariate ", colnames(pData(bs))[testCovariate],
              " is continuous.")
        fact <- FALSE
    }else{
        message("Assuming the test ",
              "covariate ", colnames(pData(bs))[testCovariate],
              " is a factor.")
        if(min(sampleSize) < 2)
          stop("At least one group has only one sample! ",
               "Replicates are required to run dmrseq.")
        testCov <- as.factor(testCov)
    }
    
    if (!is.null(adjustCovariate)) {
        mmdat <- data.frame(testCov = testCov)
        adjustCov <- pData(bs)[, adjustCovariate, drop = FALSE]
        
        # check for number of unique values per adjust cov
        nunq <- apply(adjustCov, 2, function(x) length(unique(x)))
        if (any(nunq < 2))
          stop("At least one adjust covariate is constant across samples.",
               " Please remove this covariate from the model and try again.")
        
        # remove any empty factor levels
        for (f in 1:ncol(adjustCov)){
          if (is.factor(adjustCov[,f]))
            adjustCov[,f] <- droplevels(adjustCov[,f])
        }
        
        mmdat <- cbind(mmdat, adjustCov)
        frm <- paste0("~", paste0(colnames(mmdat), collapse = " + "))
        design <- model.matrix(as.formula(frm), data=mmdat)
        colnames(design)[coeff] <- colnames(pData(bs))[testCovariate]
        coeff.adj <- (max(coeff) + 1):(ncol(design))
    } else {
        design <- model.matrix(~testCov)
        colnames(design)[coeff] <- colnames(pData(bs))[testCovariate]
        coeff.adj <- NULL
    }
    
    # check model matrix is full rank
    e <- eigen(crossprod(as.matrix(design)), symmetric = TRUE, only.values = TRUE)$values
    if (! (e[1] > 0 && abs(e[length(e)]/e[1]) > 1e-13)){
      stop("Design matrix is not full rank")
    }
    
    # check for empty factor levels in design matrix
    if (sum(colSums(design) == 0) > 0){
      which.empty <- which(colSums(design) == 0)
      design <- design[,-which.empty]
    }
    
    # check that p <= n
    if (ncol(bs) < ncol(design) + 1)
      stop("Not enough degrees of freedom to estimate ", ncol(design)-1,
           " covariates using ", ncol(bs), " samples. Please use a larger ",
           "number of samples, or specify fewer adjust covariates.")
  
    # check for incompatible args
    if (fact && !is.null(matchCovariate) && length(unique(testCov)) > 2)
      stop("matchCovariate can't be used when testCovariate is not a 2-group ",
           "factor. Perhaps you'd like to add an adjustCovariate instead?")
    
    if(!is.null(matchCovariate) && choose(nrow(design), min(sampleSize)) >= 5e5)
      stop("matchCovariate can't be used when the sample size is large enough ",
           "to yield more than 500000 possible permutations. ",
           "Perhaps you'd like to add an adjustCovariate instead?")
    
    # check for interaction terms (not yet supported)
    if (length(coeff) > 1 && any(rowSums(design[,coeff]) > 1))
      stop("Interaction terms in testCovariate are not yet supported.")
    
    if (length(unique(testCov)) == 2) {
        message("Condition: ",
            unique(pData(bs)[, testCovariate][which(design[, coeff] == 1)]), 
            " vs ", 
            unique(pData(bs)[, testCovariate][which(design[, coeff] == 0)]))
    }
    if (!is.null(adjustCovariate)) {
      message("Adjusting for covariate (s): ", 
              paste(colnames(pData(bs))[adjustCovariate], collapse = ", "))
    }
    if (!is.null(matchCovariate)) {
        if (length(matchCovariate) > 1)
          stop("Covariate matching can only be carried out for one",
              " covariate")
        if (length(unique(testCov)) > 2)
          stop("Covariate matching can only be carried out for 2-group",
               " comparisons")
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
      message("Matching permutations on covariate: ", 
              colnames(pData(bs))[mC])
    }
    
    # check for loci with missing data
    if (fact){
      lev <- unique(pData(bs)[[testCovariate]])
      filter <- NULL
      for (l in seq_along(lev)){
        filter <- rbind(filter,
              1*(DelayedMatrixStats::rowSums2(getCoverage(bs)[,pData(bs)[[testCovariate]] == 
                                           lev[l], drop = FALSE]) == 0))
      }
      filter <- which( apply(filter, 2, max) > 0 )
  
      if (length(filter) > 0) {
        stop(length(filter), " loci have zero coverage in all samples ",
             "of at least one condition. Please remove these loci ", 
             "before running dmrseq")
      }
      
    }else{
      filter <- DelayedMatrixStats::rowSums2(getCoverage(bs)==0) >= ncol(bs) - 1
      if(sum(filter) > 0)
        stop(sum(filter), " loci have zero coverage in at least ",
             ncol(bs) - 1, " samples. Please remove these loci ", 
             "before running dmrseq")
    }
    
    # register the parallel backend
    BiocParallel::register(BPPARAM)
    backend <- paste0("BiocParallel:", class(bpparam())[1])
    
    if (bpparam()$workers == 1) {
      if (verbose) {
        mes <- "Using a single core (backend: %s)."
        message(sprintf(mes, backend))
      }
      parallel <- FALSE
    } else {
      if (verbose) {
        mes <- paste0("Parallelizing using %s workers/cores ", 
                      "(backend: %s).")
        message(sprintf(mes, bpparam()$workers, backend))
      }
      parallel <- TRUE
    }
    message("Computing on ", chrsPerChunk, 
            " chromosome(s) at a time.\n")
    
    message("Detecting candidate regions with coefficient larger than ",
                   unique(abs(cutoff)), 
           " in magnitude.")
    OBS <- bumphunt(bs=bs, design = design, 
                    coeff = coeff, coeff.adj = coeff.adj, minInSpan = minInSpan,
                    minNumRegion = minNumRegion, cutoff = cutoff, 
                    maxGap = maxGap, maxGapSmooth = maxGapSmooth, 
                    smooth = smooth, bpSpan = bpSpan, verbose = verbose, 
                    parallel = parallel, block = block, blockSize = blockSize,
                    chrsPerChunk = chrsPerChunk, fact = fact,
                    adjustCovariate = adjustCovariate)
   
    # check that at least one candidate region was found; if there were none 
    # there is no need to go on to compute permutation tests...
    
    if (length(OBS) > 0) {
        message("* ", nrow(OBS), " candidates detected")
        FLIP <- NULL
        # configure the permutation matrix for two group comparisons
        if (length(unique(design[, coeff[1]])) == 2 && 
            length(coeff) == 1 &&
            choose(nrow(design), min(sampleSize)) < 5e5 ) {
            if (verbose) {
                message("Performing balanced permutations of ",
                        "condition across samples ", 
                  "to generate a null distribution of region test statistics")
            }
            perms <- combn(seq(1, nrow(design)), min(sampleSize))
            
            # Remove redundant permutations (if balanced)
            if (length(unique(table(design[,coeff]))) == 1){
              perms <- perms[, seq_len(ncol(perms)/2)]
            }
            
            # restrict to unique permutations that don't include any 
            # groups consisting of all identical conditions
            rmv <- NULL
            for (p in seq_len(ncol(perms))){
              if (length(unique(design[perms[,p],coeff])) == 1){
                rmv <- c(rmv, p)
              }
            }
            if (length(rmv) > 0 )
              perms <- perms[,-rmv]
            
            # subsample permutations based on similarity to original partition
            # gives preference to those with the least similarity
            if (maxPerms < ncol(perms)) {
               similarity <- apply(perms, 2, function(x) {
                 max(table(design[x,coeff]))
               })
               perms.all <- perms
               perms <- NULL
               levs <- sort(unique(similarity))
               l <- 1
               num <- 0
               while(!(num == maxPerms) && l <= length(levs)) {
                 keep <- sample(which(similarity == levs[l]), 
                                min(maxPerms-num, sum(similarity == levs[l])) )
                 perms <- cbind(perms, perms.all[,keep])
                 l <- l + 1
                 num <- ncol(perms) 
               }
            }
        } else {
            # Next consider a multilevel, or continuous covariate where the 
            # covariate will be permuted in an unrestricted manner
            if (verbose) {
                message("Performing unrestricted permutation of", 
                  " covariate of interest across samples ", 
                  "to generate a null distribution of region test statistics")
            }
            perms <- as.matrix(seq_len(nrow(design)))

            for (p in seq_len(maxPerms)) {
                tries <- 0
                candidate <- sample(seq_len(nrow(design)), nrow(design))
                # check that the permutation is not a duplicate, and not 
                # equal to the original
                while ((sum(apply(perms, 2, function(x) 
                                all.equal(x, candidate)) == TRUE) > 0 || 
                       sum(apply(perms, 2, function(x) 
                         all.equal(x, rev(candidate))) == TRUE) > 0) &&
                       tries <= 20) {
                  candidate <- sample(seq(seq_len(nrow(design))), nrow(design))
                  tries <- tries + 1
                }
                # save the permutation to the permutation matrix
                if (tries <= 20){
                  perms <- cbind(perms, candidate)
                }
            }
            perms <- perms[,-1] # remove original
        }
        
        pData.orig <- pData(bs)
        levs <- unique(pData.orig[[testCovariate]])
        # Now rerun on permuted designs and concatenate results
        for (j in seq_len(ncol(perms))) {
            if (verbose) {
                message("\nBeginning permutation ", j)
            }
            reorder <- perms[, j]
            designr <- design
            
            if (length(unique(design[, coeff[1]])) == 2 && 
                length(coeff) == 1 && 
                !nrow(perms) == nrow(designr)) {
                designr[, coeff] <- 0
                designr[reorder, coeff] <- 1
                pData(bs)[[testCovariate]] <- levs[1]
                pData(bs)[[testCovariate]][reorder] <- levs[2]
                
                if (!all(sort(pData.orig[[testCovariate]]) ==
                              sort(pData(bs)[[testCovariate]]))){
                  designr[, coeff] <- 1
                  designr[reorder, coeff] <- 0
                  pData(bs)[[testCovariate]] <- levs[2]
                  pData(bs)[[testCovariate]][reorder] <- levs[1]
                }
                
                xr <- NULL
                for (rd in seq_len(nrow(pData.orig))) {
                  match <- which(pData.orig[[testCovariate]] %in%
                              pData(bs)[rd,][[testCovariate]])
                  taken <- which(match %in% xr)
                  if (length(taken) > 0)
                    match <- match[-taken]
                  if (length(match) > 0)
                    xr <- c(xr, match[1])
                }
                if(length(coeff.adj) > 0){
                  pData(bs)[,adjustCovariate] <- 
                    pData.orig[xr,adjustCovariate]
                }
            } else {
                designr[, coeff] <- designr[reorder, coeff]
                pData(bs) <- pData.orig[reorder, , drop = FALSE]
            }
            
            # if matchCovariate is not null, restrict permutations such that 
            # null comparisons are balanced for the values of 
            # pData$matchCovariate this avoids comparison of,
            # say two different individuals in the null, that the comparison of 
            # interest is tissue type. Not matching would mean the null is 
            # really not null
            if (!is.null(matchCovariate)) {
                permLabel <- paste0(paste0(pData(bs)[designr[, coeff[1]] == 1, 
                                                     mC], collapse = "_"), 
                  "vs", paste0(pData(bs)[(1 - designr[, coeff[1]]) == 1, 
                                           mC], collapse = "_"))
                
                c1 <- unlist(strsplit(permLabel, "vs"))[1]
                c2 <- unlist(strsplit(permLabel, "vs"))[2]
                
                c1 <- unlist(strsplit(c1, "_"))
                c2 <- unlist(strsplit(c2, "_"))
                
                keepPerm <- 1 * (sum(c1 %in% c2) > 0 &&
                                 sum(c2 %in% c1) > 0)
                
                if (keepPerm == 0) {
                  if (verbose) {
                    message(paste0("Skipping permutation ", 
                                   gsub("vs", " vs ", permLabel)))
                  }
                  next
                }
            } else {
                permLabel <- j
            }
            
            res.flip.p <- bumphunt(bs=bs, design = designr, 
                                   coeff = coeff, 
                                   coeff.adj = coeff.adj,
                                   minInSpan = minInSpan, 
                                   minNumRegion = minNumRegion, cutoff = cutoff,
                                   maxGap = maxGap, maxGapSmooth = maxGapSmooth,
                                   smooth = smooth, bpSpan = bpSpan, 
                                   verbose = verbose, parallel = parallel,
                                   block = block, blockSize = blockSize,
                                   chrsPerChunk = chrsPerChunk, fact = fact,
                                   adjustCovariate = adjustCovariate)
            
            if (verbose) {
              message("* ", j, " out of ", ncol(perms), 
                      " permutations completed (",
                      nrow(res.flip.p), " null candidates)")
            }
            
            if (!is.null(res.flip.p)) {
                res.flip.p$permNum <- permLabel
                FLIP <- rbind(FLIP, res.flip.p)
            }
        }
        
        # restore original pData
        pData(bs) <- pData.orig
        
        # if no candidates were found in permutation
        # provide informative error message
        if (is.null(FLIP)){
          warning("No candidate regions found in permutation, so inference ",
               "can't be carried out. ",
               "Try decreasing the cutoff, or running on a larger ",
               "dataset if you are currently using a subset.")
          OBS$pval <- NA
          OBS$qval <- NA
        }else if (nrow(FLIP) < 0.05*nrow(OBS)){
          message("Note: Very few null candidate regions were found.",
               "For more accurate and sensitive inference, ",
               "try decreasing the cutoff, or running on a larger ",
               "dataset if you are currently using a subset.")
        }
        
        if (!is.null(FLIP)){
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
          perm.ordered <- c(sort(abs(FLIP[, whichStatF]), 
                                 method = "quick"), Inf)
          
          # Step 2: find the first instance in the sorted vector where the 
          # permuted value is greater than the observed and use this to 
          # determine the number of permuted values that are greater than or  
          # equal to theobserved
          pval <- rep(NA, nrow(OBS))
          pval[!is.na(OBS[, whichStatO])] <- (1 + 
                        vapply(abs(OBS[!is.na(OBS[, whichStatO]), whichStatO]),
            function(x) length(perm.ordered) - min(which(x <= perm.ordered)),
                        numeric(1))) / (1 + sum(!is.na(FLIP[, whichStatF])))
          
          # missing test statistics cause Inf for the p-value calculation 
          # instead propagate the missing values
          pval[abs(pval) == Inf] <- NA
          
          pval <- data.frame(x = pval, y = p.adjust(pval, method = "BH"))
          
          OBS$pval <- pval$x
          OBS$qval <- pval$y
        }
        
        # convert output into GRanges, with indexStart/indexEnd as IRanges
        indexIR <- IRanges(OBS$indexStart, OBS$indexEnd)
        OBS.gr <- makeGRangesFromDataFrame(OBS[,-c(4:5)], 
                                           keep.extra.columns = TRUE)
        OBS.gr$index <- indexIR
        names(OBS.gr) <- NULL
        
        # sort on pval overall (currently sorted within chromsome)
        OBS.gr <- OBS.gr[order(OBS.gr$pval, -abs(OBS.gr$stat)),]
        
        return(OBS.gr)
    } else {
        message("No candidate regions pass the cutoff of ", unique(abs(cutoff)))
        return(NULL)
    }
}
