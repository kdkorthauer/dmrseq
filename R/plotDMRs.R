#' Plot Differentially Methylated Regions
#'
#' Generates trace plots of methylation proportions by genomic position.
#' 
#' Creates aesthetially pleasing DMR plots. By default will plot individual
#' points with size proportional to coverage, along with a smoothed line 
#' for each sample. Elements will be colored by biological condition 
#' (\code{label}). Also has functionality to add annotations below the main
#' plot (CpG category, genes) if \code{annoTrack} is specified.
#' 
#' @param BSseq An object of class BSseq.
#' 
#' @param regions A data.frame containing the DMRs (output from the main
#' \code{dmrseq}) function.
#' 
#' @param extend Describes how much the plotting region should be extended in
#'  either direction. The total width of the plot is equal to the width of the 
#'  region plus twice extend.
#' 
#' @param main The plot title. The default is to construct a title with 
#' information about which genomic region is being plotted.
#' 
#' @param addRegions A set of additional regions to be highlighted on the 
#' plots. Same format as the \code{regions} argument.
#' 
#' @param annoTrack a \code{GRangesList} object with two elements returned
#' by \code{\link{getAnnot}}. The first
#' contains CpG category information in the first element (optional)
#' coding gene sequence information in the second element (optional).
#' At least one of these elements needs to be non-null in order for 
#' any annotation to be plotted, but it is not necessary to contain
#' both.
#' 
#' @param col The color of the methylation estimates. It is recommended to 
#' leave this value as default (NULL), and specify a value of 
#' \code{testCovariate} to indicate which column of \code{pData(bs)}
#' to use as a factor for coloring the points and lines of the plot.
#' Alternatively, you can specify particular colors by 
#' passing this information through the \code{pData} slot of the
#' object \code{BSseq} (a \code{data.frame} that houses metadata). To do 
#' so, place the color value for each sample in a column titled \code{col}, 
#' and leave this argument as its default value of NULL. Alternatively,
#' you may specify a vector of color values (one for each sample), but 
#' you *must* make sure that this vector is in the same order as the samples
#' are in the BSseq object. If NULL and no \code{col} column is found in 
#' \code{pData}, then estimates are plotted in black for all samples.
#' 
#' @param lty The line type of the methylation estimates. It is recommended to
#' pass this information through the \code{pData} slot of the
#' object \code{BSseq} (a \code{data.frame} that houses metadata). To do 
#' so, place the line type value for each sample in a column titled \code{lty}, 
#' and leave this argument as its default value of NULL. Alternatively,
#' you may specify a vector of line type values (one for each sample), but 
#' you *must* make sure that this vector is in the same order as the samples
#' are in the BSseq object. If NULL and no \code{lty} column is found in 
#' \code{pData}, then estimates are plotted with \code{lty=1} for all samples.
#' 
#' @param lwd The line width of the methylation estimates. It is recommended to
#' pass this information through the \code{pData} slot of the
#' object \code{BSseq} (a \code{data.frame} that houses metadata). To do 
#' so, place the line width value for each sample in a column titled \code{lwd}, 
#' and leave this argument as its default value of NULL. Alternatively,
#' you may specify a vector of line width values (one for each sample), but 
#' you *must* make sure that this vector is in the same order as the samples
#' are in the BSseq object. If NULL and no \code{lwd} column is found in 
#' \code{pData}, then estimates are plotted with \code{lwd=1} for all samples.
#' 
#' @param label The condition/population labels for the plot legend. These
#' should have a 1 to 1 correspondence with \code{col} (i.e. all samples of
#' one color should correspond to the same condition/population label). It
#' is recommended to leave this value as default (NULL), and s
#' It is recommended to specify a value of 
#' \code{testCovariate} to indicate which column of \code{pData(bs)}
#' to use as the condition labels. Alternatively, you can pass in 
#' arbitrary labels by 
#' adding this information through the \code{pData} slot of the
#' object \code{BSseq} (a \code{data.frame} that houses metadata). To do 
#' so, place the labels for each sample in a column titled \code{label}, 
#' and leave this argument as its default value of NULL. Alternatively,
#' you may specify a vector of labels (one for each sample), but 
#' you *must* make sure that this vector is in the same order as the samples
#' are in the BSseq object. If NULL and no \code{label} column is found in 
#' \code{pData}, then no legend is created.
#' 
#' @param mainWithWidth logical value indicating whether the default title
#' should include information about width of the plot region.
#' 
#' @param regionCol The color used for highlighting the region.
#' 
#' @param addTicks logical value indicating whether tick marks showing the
#'  location of methylation loci should be added. Default is TRUE.
#' 
#' @param addPoints logical value indicating whether the individual 
#' methylation estimates be plotted as points.
#' 
#' @param pointsMinCov The minimum coverage a methylation loci need in 
#' order for the raw methylation estimates to be plotted. Useful for filtering
#' out low coverage loci. Only used if addPoints = TRUE. Default value is 1 
#' (no filtering).
#' 
#' @param highlightMain logical value indicating whether the plot region 
#' should be highlighted? Should the function be verbose?
#' 
#' @param stat logical value indicating whether the region statistic 
#' should be displayed in the plot title. The value is extracted from the
#' \code{regions} argument.
#' 
#' @param qval logical value indicating whether the region FDR estimate  
#' (q-value) should be displayed in the plot title. The value is extracted 
#' from the \code{regions} argument.
#' 
#' @param verbose logical value indicating whether progress messages 
#' should be printed to the screen.
#' 
#' @param testCovariate integer value or vector indicating which of columns of
#'  \code{pData(bs)} contains the covariate of interest that will be tested 
#'  for association of methylation levels. This is used to construct the 
#'  design matrix for the test statistic calculation.
#' 
#' @export
#' 
#' @return None (generates a plot)
#' 
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices hcl
#' 
#' @examples
#' 
#' # load the example data
#' data(BS.chr21)
#' 
#' # load example results (computed with dmrseq function)
#' data(dmrs.ex)
#' 
#' # get annotation information (using getAnnot function)
#' # here we'll load the example annotation from chr21
#' data(annot.chr21)
#' 
#' # plot the 1st DMR
#' plotDMRs(BS.chr21, regions=dmrs.ex[1,], testCovariate=1,
#'    annoTrack=annot.chr21)
#' 
plotDMRs <- function(BSseq, regions = NULL, testCovariate=NULL,
                            extend = (regions$end-regions$start+1)/2, 
                            main = "", 
                            addRegions = regions, annoTrack = NULL, 
                            col = NULL, lty = NULL, lwd = NULL, label=NULL,
                            mainWithWidth = TRUE, 
                            regionCol = .alpha("#C77CFF", 0.2), 
                            addTicks = TRUE, addPoints = TRUE, 
                            pointsMinCov = 1, 
                            highlightMain = FALSE, 
                            qval=TRUE, stat=TRUE,
                            verbose = TRUE) {
  # adapted from plotManyRegions from bsseq plot to 
  # take in a vector of qval values 
  # (1 per region in regions 
  # argument) to be displayed in the plot title.
  # set addPoints = TRUE to plot individual points sized by coverage 
  # and one smooth (loess)
  # line per sample instead of a uniform-sized verbatim line going through 
  # each observation
  if (verbose) message(paste0("[plotDMRs] Plotting ", nrow(regions), 
                              " DMRs"))
    if(!is.null(regions)) {
        if(is(regions, "data.frame"))
            gr <- data.frame2GRanges(regions, keepColumns = FALSE)
        else
            gr <- regions
        if(!is(gr, "GRanges"))
          stop(paste0("'regions' needs to be either a 'data.frame' ", 
                        "(with a single row) or a 'GRanges' ", 
                        "(with a single element)")) 
    } else {
        gr <- granges(BSseq)
    }
    gr <- resize(gr, width = 2*extend + width(gr), fix = "center")
    BSseq <- subsetByOverlaps(BSseq, gr)
    
    if(length(start(BSseq)) == 0)
        stop("No overlap between BSseq data and regions")
    if(!is.null(main) && length(main) != length(gr))
        main <- rep(main, length = length(gr))
    
    opar <- par(mar = c(0,4.1,0,0), oma = c(0,0,2.5,1))
    on.exit(par(opar))
   	
   	if (length(extend)==1){
   		extend <- rep(extend, length(gr))
   	}
    
    if(!is.null(testCovariate)){
      coeff <- 2:(2+length(testCovariate)-1)
      design <- model.matrix( ~ pData(BSseq)[,testCovariate])

      if(is.null(col)){
        cov.unique <- unique(design[,coeff])
        colors <- gg_color_hue(length(cov.unique))
        if (length(cov.unique) == 2){
          colors <- c("mediumblue", "deeppink1")
        }
        colors <- cbind(cov.unique, colors[rank(as.numeric(cov.unique))])
        z <- colors[,2][match(design[,coeff], colors[,1])]
        pData(BSseq)$col <- as.character(z)
      }
      
      if(is.null(label)){
        pData(BSseq)$label <- paste0(pData(BSseq)[,testCovariate])
      }
    }
    
    for(ii in seq(along = gr)) {
        if(verbose & ii%%100==0){
          cat(sprintf("..... Plotting region %d (out of %d)\n",
                      ii, nrow(regions)))
        } 

        .plotSingleDMR(BSseq = BSseq, region = regions[ii,], 
                      extend = extend[ii], main = main[ii], 
                      col = col, lty = lty, lwd = lwd, label=label,
                      addRegions = addRegions, regionCol = regionCol, 
                      mainWithWidth = mainWithWidth,
                      annoTrack = annoTrack, addTicks = addTicks, 
                      addPoints = addPoints,
                      pointsMinCov = pointsMinCov, 
                      highlightMain = highlightMain,
                      qval=qval, stat=stat)
    }
}


