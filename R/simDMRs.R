
#' Simulate Differentially Methylated Regions
#' 
#' Add simulated DMRs to observed control data. Control data will be split
#' into two (artificial) populations.
#' 
#' @param bs a BSseq object containing only control samples (from the same
#' population) for which simulated DMRs will be added after dividing the 
#' population into two artificial groups.
#' 
#' @param num.dmrs an integer specifying how many DMRs to add.
#' 
#' @param delta.max0 a proportion value indicating the mode value for the
#' difference in proportion of methylated CpGs in the simulated DMRs (the
#' actual value will be drawn from a scaled Beta distribution centered at 
#' this value). Default value is 0.3.
#' 
#' @return A named list object with 5 elements: (1) 
#' \code{gr.dmrs} is a \code{GenomicRanges} object with \code{num.dmrs} 
#' ranges that represent the random DMRs added. (2) \code{dmr.mncov} is a 
#' numeric vector that contains the mean coverage in each simulated DMR. (3)
#' \code{dmr.L} is a numeric vector that contains the number of CpGs in each 
#' simulated DMR. (4) \code{bs} is the BSseq object that contains the 
#' simulated DMRs. (5) \code{deltas} is a numeric vector that contains the 
#' effect size used for each DMR. 
#' 
#' @importFrom IRanges IRanges
#' 
#' @export
#' 
#' @examples 
#' 
#' # Add simulated DMRs to a BSseq dataset
#' # This is just for illustrative purposes - ideally you would
#' # add DMRs to a set of samples from the same condition (in our
#' # example data, we have data from two different cell types)
#' # In this case, we shuffle the samples by cell type to create
#' # a null comparison.
#' 
#' data(BS.chr21)
#' 
#' BS.chr21.sim <- simDMRs(bs=BS.chr21[1:10000,c(1,3,2,4)], 
#'                         num.dmrs=50)
#' 
#' # show the simulated DMRs GRanges object
#' show(BS.chr21.sim$gr.dmrs)
#' 
#' # show the updated BSseq object that includes the simulated DMRs
#' show(BS.chr21.sim$bs)
#' 
#' # examine effect sizes of the DMRs
#' head(BS.chr21.sim$delta)
#' 
simDMRs <- function(bs, num.dmrs = 3000, delta.max0 = 0.3) {
    # check that all loci have coverage in every sample
    zero.cov <- which(rowSums(as.matrix(getCoverage(bs, 
                                             type="Cov")) == 0) > 0)
    if (length(zero.cov) > 0){
      stop("Zero coverage found for at least one sample in ", length(zero.cov),
           " loci. Please filter for loci with coverage at least one in ",
           "all samples before passing to 'simDMRs'")
    }
  
    sampleSize <- floor(nrow(pData(bs))/2)
    
    # code to simulate DMRs if some number of simulated dmrs was specified
    message("Simulating DMRs for ", sampleSize, " vs ", ncol(bs) - sampleSize, 
            " comparison")
    triwt <- function(x, amp = 1, base = 0, width = 1, center = 0, 
                      deg = 3, dir = 1) {
        y <- dir * (((width/2)^deg - abs(x - center)^deg)/
                      (width/2)^deg)^deg * amp +  base
        y[abs(x - center) > ceiling(width/2)] <- base[abs(x - center) >
                                                        ceiling(width/2)]
        return(y)
    }
    
    meth.mat <- as.matrix(getCoverage(bs, type = "M"))
    unmeth.mat <- as.matrix(getCoverage(bs, type = "Cov")) - meth.mat
    chr <- as.character(seqnames(bs))
    pos <- start(bs)
    
    cluster <- bumphunter::clusterMaker(chr, pos, maxGap = 500)
    Indexes <- split(seq(along = cluster), cluster)
    lns <- lengths(Indexes)
    Indexes <- Indexes[lns >= 5 & lns <= 500]
    
    # sample regions with intermediate methylation values preferentially
    prop.mat <- rowMeans(meth.mat/(meth.mat + unmeth.mat))
    prop.mat <- unlist(lapply(Indexes, function(x) median(prop.mat[x])))
  
    dmrs.ind <- sample(seq_len(length(Indexes)), num.dmrs, replace = FALSE, 
                       prob = pmax(1 - sqrt(2) * 
                                     abs(0.5 - prop.mat)^0.5, 0))
    dmrs.ind <- Indexes[dmrs.ind]
    fnc <- function(index) {
        gr.dmr <- GRanges(seqnames = unique(as.character(seqnames(bs)[index])),
                          IRanges(start = min(start(bs)[index]), 
                                  end = max(start(bs)[index])))
        return(gr.dmr)
    }
    ## GenomicRanges Object for the Simulated DMRs
    gr.dmrs <- suppressWarnings(Reduce("c", lapply(dmrs.ind, fnc)))
    
    ## Generating the Methylated and Unmethylated Read Counts 
    ## for the CpG sites in the
    ## DMRs and outside
    
    # set up null signal (smooth function of position) mcols(bs)$diff <- 0
    Diff <- Diff2 <- rep(0, length(bs))
    
    dmr.mncov <- dmr.L <- deltas <- rep(NA, num.dmrs)
    
    for (u in seq_len(num.dmrs)) {
        # coin flip for up or down
        up <- 1 - 2 * (rbinom(1, 1, 0.5) == 1)
        
        # let effect size change randomly
        delta.max <- delta.max0 + (rbeta(1, 2, 2) - 0.5)/3
        deltas[u] <- delta.max
        
        # grab loci in the dmr
        dmr.L[u] <- length(dmrs.ind[[u]])
        prop.mat <- meth.mat[dmrs.ind[[u]], ]/(meth.mat[dmrs.ind[[u]], ] + 
                                                 unmeth.mat[dmrs.ind[[u]], 
            ])
        
        # change direction if baseline mean is near boundary
        if (up == 1) {
            if (mean(prop.mat) > 1 - delta.max) {
                up <- -1
            }
        } else if (up == -1) {
            if (mean(prop.mat) < delta.max) {
                up <- 1
            }
        }
        
        # simulated mean as a smooth parabola added or subtracted 
        # from baseline mean
        last <- max(pos[dmrs.ind[[u]]])
        first <- min(pos[dmrs.ind[[u]]])
        width <- last - first
        
        # widen out so that first and last CpGs don't have a 
        # difference of zero
        last <- last + 0.2 * width
        first <- first - 0.2 * width
        width <- last - first
        
        mid <- round((last - first)/2 + first)
        
        # Diff.hit is the methylation percentage difference for each position in
        # spiked in DMR; restricted to between -1 and 1 (negative indicates that
        # difference is in the negative direction.
        
        Diff.hit <- round(triwt(pos[dmrs.ind[[u]]], amp = delta.max, 
                                base = Diff[dmrs.ind[[u]]], 
            width = width, center = mid, deg = 4, dir = up), 4)
        
        # calculate the mean coverage in the DMR over all the samples and save 
        # result in a vector dmr.mncov for exploring the characteristics of the
        # detected / missed DMRs in simulation results.
        mn.cov <- by(t(meth.mat[dmrs.ind[[u]], ] + unmeth.mat[dmrs.ind[[u]], ]),
            factor(paste0("Condition", c(rep(1, sampleSize), 
                                         rep(2, ncol(bs) - sampleSize)))), 
            colMeans)
        mn.cov <- rowMeans(cbind(mn.cov[[1]], mn.cov[[2]]))
        dmr.mncov[u] <- mean(mn.cov)
        
        
        # Conditional on the coverage for each site and sample combination,
        # sample the
        # number of methylated reads from a binomial distribution where the 
        # probability
        # parameter is taken as the observed estimate plus the Diff/Diff2 
        # (depending on
        # whether the sample is from condition 1 or condition 2 This will 
        # induce sampling
        # error into the number of reads observed in the simulation and 
        # thus will be more
        # realistic and a more practical evaluation of the methods
        cov <- meth.mat[dmrs.ind[[u]], ] + unmeth.mat[dmrs.ind[[u]], ]
        prop <- meth.mat[dmrs.ind[[u]], ]/cov
        grp <- runif(1) < 0.5
        ss <- ifelse(grp, sampleSize, ncol(bs) - sampleSize)
        for (samp in seq_len(ss)) {
            # randomly choose which condition is the one with the difference
            if (grp) {
                # first generate M counts for sample samp of condition 1
                meth.mat[dmrs.ind[[u]], samp] <- rbinom(n=length(dmrs.ind[[u]]),
                  size = cov[, samp], 
                  prob = pmax(pmin(prop[, samp] + Diff.hit, 1), 
                    0))
                
                # next assign the other Cov - M counts to unmethylated matrix
                unmeth.mat[dmrs.ind[[u]], samp] <- cov[, samp] -
                  meth.mat[dmrs.ind[[u]], samp]
                
            } else {
                # next generate M counts for sample samp of condition 2
                meth.mat[dmrs.ind[[u]], (sampleSize + samp)] <- 
                        rbinom(n = length(dmrs.ind[[u]]),
                  size = cov[, (sampleSize + samp)], 
                  prob = pmax(pmin(prop[, (sampleSize + samp)] + 
                                     Diff.hit, 1), 0))
                
                # next assign the other Cov - M counts to unmethylated matrix
                unmeth.mat[dmrs.ind[[u]], (sampleSize + samp)] <- cov[,
                  (sampleSize + samp)] - meth.mat[dmrs.ind[[u]],
                                                  (sampleSize + samp)]
            }
        }
    }

    # get everything in order to run bumphunter functions
    sampnames <- paste0("Condition", c(rep(1, sampleSize),
                                       rep(2, ncol(bs) - sampleSize)), 
                        "_Rep", c(seq_len(sampleSize), 
                                  seq_len(ncol(bs) - sampleSize)))
    colnames(meth.mat) <- colnames(unmeth.mat) <- sampnames 
    bsNew <- BSseq(pos = pos, chr = chr, M = meth.mat, 
                   Cov = (meth.mat + unmeth.mat), 
                   sampleNames = sampnames)
    
    sim.dat.red <- list(gr.dmrs = gr.dmrs, dmr.mncov = dmr.mncov, dmr.L = dmr.L,
        bs = bsNew, delta = deltas)
    return(sim.dat.red)
}

