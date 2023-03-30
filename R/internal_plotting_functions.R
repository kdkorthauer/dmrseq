#' @title Add annotations to DMR plots
#' 
#' @description Function to add visual representation of CpG categories 
#' and/or coding 
#' sequences to DMR plot
#' 
#' @details An internal function that takes an annotation 
#' \code{SimpleGRangesList}
#' object that
#' contains CpG category information in the first element (optional) and / or
#' coding gene sequence information in the second element (optional). If neither
#' of these are present, then nothing will be plotted.
#' 
#' @param gr a \code{GRanges} object that contains the DMRs to be
#'  plotted
#' 
#' @param annoTrack a \code{SimpleGRangesList} object with two elements.
#' The first contains CpG category information in the first element (optional)
#' coding gene sequence information in the second element (optional).
#' At least one of these elements needs to be non-null in order for 
#' any annotation to be plotted, but it is not necessary to contain
#' both.
#' 
#' @return None
#' 
dmrPlotAnnotations <- function(gr, annoTrack) {
    # Code adapted from bsseq package
    
    ## check may need to be modified
    if (!is(annoTrack, "SimpleGRangesList"))
        stop("'annoTrack' needs to be a 'SimpleGRangesList'")
    plot(start(gr), 1, type = "n", xaxt = "n", yaxt = "n", bty = "n", 
        ylim = c(0, length(annoTrack) + 0.5), xlim = c(start(gr), end(gr)), 
        xlab = "", ylab = "")
    
    # add legend
    vars <- list(Island = "Island   ", Shore = "Shore   ", Shelf = "Shelf   ", 
                 OpenSea1 = "Open ", OpenSea2 = "Sea")
    cols <- c("forestgreen", "goldenrod2", "dodgerblue", "blue3", "blue3")
    for (i in seq_len(length(vars))) {
        tmpvars <- vars
        tmpvars[-i] <- paste("phantom('", tmpvars[-i], "')", sep = "")
        expr <- paste(tmpvars, collapse = "*")
        text(start(gr), 1.8, parse(text = expr), col = cols[i], adj = c(0, 1),
             cex = 0.85)
    }
    
    lapply(seq(along = annoTrack), function(ii) {
        jj <- length(annoTrack) + 1 - ii
        ir <- subsetByOverlaps(annoTrack[[ii]], gr)
        start(ir) <- pmax(start(ir), start(gr))
        end(ir) <- pmin(end(ir), end(gr))
        
        if (ii == 2) {
            jj <- jj - 0.15
        }
        
        if (length(ir) > 0) {
            if (ii == 2) {
                colourCount <- length(unique(ir$symbol))
                getPalette <- colorRampPalette(.alpha(brewer.pal(max(
                  length(unique(ir$symbol)), 
                  3), "Dark2"), 0.4))
                color.pal <- getPalette(colourCount)
                names(color.pal) <- unique(ir$symbol)
                map <- match(ir$symbol, names(color.pal))
                color <- color.pal[map]
                
                bord <- "black"
            } else if (ii == 1) {
                color <- ir$type
                colvec <- c("blue3", "dodgerblue", "goldenrod2", "forestgreen")
                names(colvec) <- c("inter", "shelves", "shores", "islands")
                
                for (ucol in unique(color)){
                  ix <- agrep(ucol, names(colvec))
                  if (length(ix) > 1)
                    stop("Ambiguous CpG annotation labels")
                  color[color == ucol] <- colvec[ix]
                }
                bord <- color
                
                rect(start(ir), jj - 0.06, end(ir), jj + 0.17, col = color, 
                     border = bord)
            }
            
            if (ii == 2) {
                lastPos <- rep(NA, length(unique(ir$symbol)) - 1)
                used <- NULL
                for (k in seq_len(length(unique(ir$symbol)))) {
                  irk <- ir[ir$symbol == unique(ir$symbol)[k], ]
                  rect(min(start(irk)), jj - 0.065, max(end(irk)), jj + 0.065, 
                       col = .alpha("black", 
                    0.1), border = .alpha("black", 0.1))
                  if (!(unique(ir$symbol)[k] %in% used)) {
                    rwidth <- end(gr) - start(gr)
                    sg <- pmin(start(irk), end(irk))
                    eg <- pmax(start(irk), end(irk))
                    gwidth <- min(max(eg), end(gr)) - max(min(sg), start(gr))
                    textPos <- max(min(sg), start(gr)) + gwidth/2
                    jj.orig <- jj
                    if (sum(!is.na(lastPos)) > 0) {
                      separation <- (textPos - lastPos[k - 1])/rwidth
                      if (abs(separation) <= 0.2 && k < 3) {
                        jj <- jj - 0.29
                      } else {
                        separation <- min(abs((textPos - lastPos)/rwidth), 
                                          na.rm = TRUE)
                        if (abs(separation) <= 0.2) {
                          jj <- jj - 0.29
                        }
                      }
                    }
                    lastPos[k] <- textPos
                    text(textPos, jj - 0.375, labels = unique(irk$symbol), 
                         cex = 0.85, 
                      col = unique(color)[k])
                    jj <- jj.orig
                    used <- c(used, unique(ir$symbol)[k])
                  }
                  rect(sg, jj - 0.11, eg, jj + 0.12, col = unique(color)[k], 
                       border = bord)
                }
            }
            
        }
        mtext(names(annoTrack)[ii], side = 2, at = jj, las = 1, line = 1)
    })
    
}

.isColor <- function(x)
{
  res <- try(col2rgb(x),silent=TRUE)
  return(!"try-error"%in%class(res))
}

.dmrGetMeta <- function(object, col, lty, lwd, label) {
    ## Assumes that object has pData and sampleNames methods Code adapted from 
    ## bsseq package
    
    ## extract col
    if (is.null(col)) {
        if ("col" %in% names(pData(object))) 
            col <- pData(object)[["col"]] else col <- rep("black", 
                                                          nrow(pData(object)))
    }else if (length(col) == 1){
      if (col %in% names(pData(object)))
        col <- pData(object)[[col]] else col <- rep("black", ncol(object))
      if (!.isColor(col))
        col <- rainbow(length(unique(col)))[as.numeric(as.factor(col))]
    }
    if (length(col) != ncol(object)) 
        col <- rep(col, length.out = ncol(object))
    if (is.null(names(col))) 
        names(col) <- sampleNames(object)
    
    ## extract lty
    if (is.null(lty)) {
      if ("lty" %in% names(pData(object)))
        lty <- pData(object)[["lty"]] else lty <- rep(1, ncol(object))
    }else if (length(lty) == 1){
      if (lty %in% names(pData(object)))
        lty <- pData(object)[[lty]] else lty <- rep(1, ncol(object))
      if (!is.numeric(lty))
        lty <- as.numeric(as.factor(lty))
    }
    if (length(lty) != ncol(object)) 
        lty <- rep(lty, length.out = ncol(object))
    if (is.null(names(lty))) 
        names(lty) <- sampleNames(object)
    
    # extract lwd
    if (is.null(lwd)) {
        if ("lwd" %in% names(pData(object))) 
            lwd <- pData(object)[["lwd"]] else lwd <- rep(1.5, 
                                                          nrow(pData(object)))
    }else if (length(lwd) == 1){
      if (lwd %in% names(pData(object)))
        lwd <- pData(object)[[lwd]] else lwd <- rep(1, ncol(object))
      if (!is.numeric(lwd))
        lwd <- as.numeric(as.factor(lwd))
    }
    if (length(lwd) != ncol(object)) 
        lwd <- rep(lwd, length.out = ncol(object))
    if (is.null(names(lwd))) 
        names(lwd) <- sampleNames(object)
    
    ## extract label
    if (is.null(label)) {
        if ("label" %in% names(pData(object))) 
            label <- pData(object)[["label"]] else label <- rep(NA, 
                                                                ncol(object))
    }else if (length(label) == 1){
      if (label %in% names(pData(object)))
        label <- pData(object)[[label]] else label <- rep(NA, ncol(object))
        if (!is.character(label))
          label <- as.character(label)
    }
    if (length(label) != ncol(object)) 
        label <- rep(label, length.out = ncol(object))
    if (is.null(names(label))) 
        names(label) <- sampleNames(object)
    
    return(list(col = col, lty = lty, lwd = lwd, label = label))
}


.dmrPlotTitle <- function(gr, extend, main, mainWithWidth, 
                          qval = NULL, stat = NULL) {
    # this function creates the main title for DMR plots Code adapted from bsseq
    # package
    if (is.data.frame(gr)) 
        gr <- data.frame2GRanges(gr)
    if (length(gr) > 1) {
        warning("plotTitle: gr has more than one element")
        gr <- gr[1]
    }
    plotChr <- as.character(seqnames(gr))
    plotRange <- c(start(gr), end(gr))
    regionCoord <- sprintf("%s: %s - %s", plotChr, format(plotRange[1], 
                                                          big.mark = ",", 
        scientific = FALSE), format(plotRange[2], big.mark = ",", 
                                    scientific = FALSE))
    if (mainWithWidth) {
        regionWidth <- sprintf("width = %s", format(width(gr), 
                                                    big.mark = ",",
                                                    scientific = FALSE))
        # add optional labels to plot titles
        if (!is.null(qval) && !is.null(stat)) {
            regionStat <- sprintf("Stat: %s", format(stat, big.mark = ",",
                                                     scientific = FALSE))
            regionFDR <- sprintf("FDR: %s", format(qval, big.mark = ",",
                                                   scientific = FALSE))
            regionCoord <- sprintf(paste0("%s (%s)\n%s, %s"), regionCoord,
                                   regionWidth, regionStat, regionFDR)
        } else if (!is.null(stat)) {
            regionStat <- sprintf("Stat: %s", format(stat, big.mark = ",",
                                                     scientific = FALSE))
            regionCoord <- sprintf(paste0("%s (%s)\n%s"), regionCoord,
                                   regionWidth, regionStat)
        } else if (!is.null(qval)) {
            regionFDR <- sprintf("FDR: %s", format(qval, big.mark = ",", 
                                                   scientific = FALSE))
            regionCoord <- sprintf(paste0("%s (%s)\n%s"), regionCoord, 
                                   regionWidth, regionFDR)
        } else {
            regionCoord <- sprintf("%s (%s)", regionCoord, regionWidth)
        }
    } else {
        # add optional labels to plot titles
        if (!is.null(qval) && !is.null(stat)) {
            regionStat <- sprintf("Stat: %s", format(stat, big.mark = ",",
                                                     scientific = FALSE))
            regionFDR <- sprintf("FDR: %s", format(qval, big.mark = ",", 
                                                   scientific = FALSE))
            regionCoord <- sprintf(paste0("%s\n%s, %s"), regionCoord, 
                                   regionStat, regionFDR)
        } else if (!is.null(stat)) {
            regionStat <- sprintf("Stat: %s", format(stat, big.mark = ",",
                                                     scientific = FALSE))
            regionCoord <- sprintf(paste0("%s\n%s"), regionCoord, regionStat)
        } else if (!is.null(qval)) {
            regionFDR <- sprintf("FDR: %s", format(qval, big.mark = ",",
                                                   scientific = FALSE))
            regionCoord <- sprintf(paste0("%s\n%s"), regionCoord, regionFDR)
        } else {
            regionCoord <- sprintf("%s", regionCoord)
        }
    }
    if (main != "") {
        main <- sprintf("%s\n%s", main, regionCoord)
    } else {
        main <- regionCoord
    }
    main
}

.dmrPlotLegend <- function(plotRange, col, label, horizLegend) {
    # this function plots a legend to the right of the plot to indicate
    # whichcolor corresponds to which samples
    
    numUnique <- length(unique(paste0(col, label, sep = "")))
    if (numUnique < length(col)) {
        col <- unique(col)
        label <- unique(label)
    }
    
    if (!horizLegend){
      for (lg in seq_len(length(label))) {
         mtext(label[lg], side = 4, line = lg - 1, col = .darken(col[lg]), 
               cex = 0.9, las = 0)
      }
    }else{
      for (lg in seq_len(length(label))) {
        mtext(label[lg], side = 4, line = 0.5, col = .darken(col[lg]), 
              cex = 0.9, las = 1, at = 1.02 - 0.08*(lg-1))
      }
    }
}


# function to transform a given color specified by a character object to a
# transparent version of that color
.makeTransparent <- function(someColor, alpha = 130) {
    newColor <- col2rgb(someColor)
    apply(newColor, 2, function(curcoldata) {
        rgb(red = curcoldata[1], green = curcoldata[2], blue = curcoldata[3], 
            alpha = alpha, maxColorValue = 255)
    })
}

.alpha <- function(col, alpha = 1) {
    if (missing(col)) 
        stop("Please provide a vector of colours.")
    apply(vapply(col, grDevices::col2rgb, matrix(rep(0,3)))/255, 2, 
          function(x) grDevices::rgb(x[1], x[2], x[3], alpha = alpha))
}


.dmrPlotPoints <- function(x, y, z, col, pointsMinCov, maxCov, 
                           regionWidth) {
    # modified from .bsPlotPoints in bsseq added functionality for point size 
    # to vary with coverage
    
    lwd <- 1.5
    # make color of points semi-transparent so that overlapping points can 
    # still be seen
    col.points <- .makeTransparent(col)
    
    # if there are a lot of CpGs to plot (the case for a block-level analysis
    # decrease the size of the plotted points since these can get very crowded
    c1 <- pmax(-0.25 * atan(3 * (length(x) - 80)/80 * pi)/atan(3 * 
               (80)/80 * pi) +  0.75, 0.5)
    ptSize <- c1 * (sqrt(z)/sqrt(maxCov) + 0.25)
    
    points(x[z >= pointsMinCov], y[z >= pointsMinCov], col = col.points, pch = 16,
        cex = ptSize)
}

.darken <- function(color, factor = 1.4) {
    col <- col2rgb(color)
    col <- col/factor
    col <- rgb(t(col), maxColorValue = 255)
    col
}



.dmrPlotLines <- function(x, y, z, col, lwd, linesMinCov, maxCov, 
                          regionWidth, lty) {
  
  if (length(x) > 100 && !is.null(lwd)) {
    lwd <- lwd + 1
  } else if (length(x) > 100) {
    lwd <- 2
  }
  
  spn <- max(1 - (1/160)*sum(z >= linesMinCov), 0.75)

  y[y==1] <- 0.99
  y[y==0] <- 0.01
  logit <- function(p){ log(p/(1-p))}
  inv.logit <- function(l){ exp(l) / (1 + exp(l)) }
  
  # don't interpolate smooth lines if there are fewer than 10 cpgs
  if (length(x) >= 10) {
    loess_fit <- loess(logit(y[z >= linesMinCov]) ~ x[z >= linesMinCov],
                       weights = z[z >= linesMinCov], span = spn)
    
    xl <- seq(min(x[z >= linesMinCov], na.rm=TRUE), 
              max(x[z >= linesMinCov], na.rm=TRUE), 
              (max(x[z >= linesMinCov], na.rm=TRUE) - 
               min(x[z >= linesMinCov], na.rm=TRUE))/500)
    lines(xl, inv.logit(predict(loess_fit,xl)), 
          col = .makeTransparent(.darken(col), 175), lwd = lwd,
          lty = lty)
    
  }else{
    lines(x[z >= linesMinCov], y[z >= linesMinCov], 
          col = .makeTransparent(.darken(col), 175), lwd = lwd,
          lty = lty)
  }
}

.dmrPlotSmoothData <- function(BSseq, region, extend, addRegions, col, lty, lwd,
    label, regionCol, addTicks, addPoints, pointsMinCov, highlightMain, 
    includeYlab = TRUE, horizLegend, addLines=TRUE, linesMinCov) {
    # modified from .plotSmoothData in bsseq to allow non-smoothed regions
    
    gr <- bsseq.bsGetGr(BSseq, region, extend)
    BSseq <- subsetByOverlaps(BSseq, gr)
    BSseq2 <- subsetByOverlaps(BSseq, bsseq.bsGetGr(BSseq, region, extend = 0))
    ## Extract basic information
    sampleNames <- sampleNames(BSseq)
    names(sampleNames) <- sampleNames
    positions <- start(BSseq)
    positions2 <- start(BSseq2)
    rawPs <- as.matrix(bsseq::getMeth(BSseq, type = "raw"))
    coverage <- as.matrix(bsseq::getCoverage(BSseq))
    
    ## get col, lwd, lty these are extracted from the pData data.frame that is 
    ## part of the bsseq object colEtc is a list object that contains col, lty,
    ## lwd and label which are used as plotting parameters label is a condition
    ##label to use for adding a legend with condition names
    
    colEtc <- .dmrGetMeta(object = BSseq, col = col, lty = lty, lwd = lwd, 
                          label = label)
    
    if (includeYlab) {
        yl <- "Methylation"
    } else {
        yl <- ""
    }
    
    ## The actual plotting starts here
    plot(positions[1], 0.5, type = "n", xaxt = "n", yaxt = "n", ylim = c(0, 1),
         xlim = c(start(gr), 
        end(gr)), xlab = "", ylab = yl)
    axis(side = 2, at = c(0, 0.25, 0.5, 0.75, 1), las = 1)
    if (addTicks) 
        rug(positions)
    
    if (is.list(addRegions) && !is.data.frame(addRegions)) {
        if (length(addRegions) > 2) {
            stop("Only two sets of regions can be highlighted")
        }
        if (length(regionCol) == 1) {
            regionCol <- c(regionCol, .alpha("blue", 0.2))
        }
        bsseq.bsHighlightRegions(regions = addRegions[[1]], gr = gr, ylim = c(0,
            1), regionCol = regionCol[1], highlightMain = highlightMain)
        bsseq.bsHighlightRegions(regions = addRegions[[2]], gr = gr, ylim = c(0,
            1), regionCol = regionCol[2], highlightMain = highlightMain)
        
    } else {
        bsseq.bsHighlightRegions(regions = addRegions, gr = gr, ylim = c(0, 1), 
                                 regionCol = regionCol, 
            highlightMain = highlightMain)
    }
    
    # add points first to avoid lines getting hidden by plotting many cpg points
    if (addPoints) {
        for(sampIdx in seq_len(ncol(BSseq))){
          .dmrPlotPoints(positions, rawPs[, sampIdx], coverage[, sampIdx], 
                         col = colEtc$col[sampIdx], 
                         pointsMinCov = pointsMinCov, 
                         maxCov = quantile(coverage, 0.95), 
                         regionWidth = end(gr) - 
                           start(gr))
        }
    }   
    
    if (addLines){
       for(sampIdx in seq_len(ncol(BSseq))){
            if (sum(coverage[, sampIdx] >= linesMinCov) > 1){
              .dmrPlotLines(positions, rawPs[, sampIdx], coverage[, sampIdx], 
                        col = colEtc$col[sampIdx], 
                        lwd = colEtc$lwd[sampIdx],
                        linesMinCov = linesMinCov, 
                        maxCov = quantile(coverage, 0.95), 
                        regionWidth = end(gr) - 
                          start(gr),
                        lty = colEtc$lty[sampIdx])
            }
        }
    }
    
    # if colEtc$label contains characters that are not null or missing, then 
    # create a legend which houses the labels as well as the colors that 
    # correspond to them -> pass in both colEtc$label as well as colEtc$col
    if (sum(!is.na(colEtc$label)) == length(colEtc$label)) {
        .dmrPlotLegend(plotRange = c(start(gr), end(gr)), 
                       colEtc$col, colEtc$label, horizLegend)
    }
    
}

# function doesn't need to be exported; not a user-level function since a single
# DMR can be plotted just fine with plotDMRs.
.plotSingleDMR <- function(BSseq, region = NULL, extend = 0, main = "", 
    addRegions = NULL, annoTrack = NULL, col = NULL, lty = NULL, lwd = NULL, 
    label = NULL, mainWithWidth = TRUE, regionCol = .alpha("orchid1", 0.2), 
    addTicks = TRUE, addPoints = FALSE, pointsMinCov = 5, highlightMain = FALSE,
    qval = NULL, stat = NULL, includeYlab = TRUE, compareTrack = NULL, 
    labelCols = NULL, horizLegend = FALSE, addLines = TRUE, linesMinCov = 1) {
    
    if(!is.null(annoTrack) || !is.null(compareTrack)){
      layout(matrix(seq_len(2), ncol = 1), heights = c(2, 1.5))
    }else{
      layout(matrix(seq_len(2), ncol = 1), heights = c(2, 0.2))
    }
    .dmrPlotSmoothData(BSseq = BSseq, region = region, extend = extend, 
        addRegions = addRegions, 
        col = col, lty = lty, lwd = lwd, label = label, 
        regionCol = regionCol, addTicks = addTicks, 
        addPoints = addPoints, pointsMinCov = pointsMinCov, 
        highlightMain = highlightMain, 
        includeYlab = includeYlab, 
        horizLegend = horizLegend, 
        addLines = addLines,
        linesMinCov = linesMinCov)
    gr <- bsseq.bsGetGr(BSseq, region, extend)
    
    if (!is.null(main)) {
        if (qval && stat) {
            qval <- round(region$qval, 4)
            stat <- round(region$stat, 3)
            main <- .dmrPlotTitle(gr = region, extend = extend, main = main, 
                                  mainWithWidth = mainWithWidth, 
                                  qval = qval, stat = stat)
        } else if (stat) {
            stat <- round(region$stat, 3)
            main <- .dmrPlotTitle(gr = region, extend = extend, main = main, 
                                  mainWithWidth = mainWithWidth, stat = stat)
        } else if (qval) {
            qval <- round(region$qval, 4)
            main <- .dmrPlotTitle(gr = region, extend = extend, main = main,
                                  mainWithWidth = mainWithWidth, qval = qval)
        } else {
            main <- .dmrPlotTitle(gr = region, extend = extend, main = main,
                                  mainWithWidth = mainWithWidth)
        }
        mtext(side = 3, text = main, outer = FALSE, cex = 0.8, line = 0)
    }
    
    if (!is.null(annoTrack)) {
        dmrPlotAnnotations(gr, annoTrack)
    } else if (!is.null(compareTrack)) {
        dmrPlotComparisons(gr, compareTrack, labelCols = labelCols)
    }
}

# pasting bsseq's .bsGetGr function since not exported
bsseq.bsGetGr <- function(object, region, extend) {
    if (is.null(region)) {
        gr <- GRanges(seqnames = seqnames(object)[1], 
                      ranges = IRanges(start = min(start(object)), 
            end = max(start(object))))
    } else {
        if (is(region, "data.frame")){
              gr <- data.frame2GRanges(region, keepColumns = FALSE)
            }else{ 
              gr <- region
            }
        if (!is(gr, "GRanges") || length(gr) != 1) 
            stop("'region' needs to be either a 'data.frame' ",
                "(with a single row) or a 'GRanges' (with a single element)")
        gr <- resize(gr, width = 2 * extend + width(gr), fix = "center")
    }
    gr
}

# pasting bsseq's .bsHighlightRegions function since not exported
bsseq.bsHighlightRegions <- function(regions, gr, ylim, regionCol, 
                                     highlightMain) {
    if (is.data.frame(regions)) 
        regions <- data.frame2GRanges(regions)
    if (highlightMain) 
        regions <- c(regions, gr)
    if (is.null(regions)) 
        return(NULL)
    regions <- subsetByOverlaps(regions, gr)
    regions <- pintersect(regions, rep(gr, length(regions)))
    if (length(regions) == 0) 
        return(NULL)
    rect(xleft = start(regions), xright = end(regions), ybottom = ylim[1], 
         ytop = ylim[2], col = regionCol, border = NA)
}

gg_color_hue <- function(n) {
    hues <- seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[seq_len(n)]
}

# function to draw nonoverlapping comparison regions (up to 3) below the main
# region plot (instead of annotations)
dmrPlotComparisons <- function(gr, annoTrack, labelCols = NULL) {
    
    if (!is(annoTrack, "SimpleGRangesList"))
      stop("'annoTrack' needs to be a 'SimpleGRangesList'")
  
    if (length(annoTrack) > 4) 
        stop("Can't plot more than 4 tracks")
    
    plot(start(gr), 1, type = "n", xaxt = "n", yaxt = "n", bty = "n", 
         ylim = c(0, length(annoTrack) + 0.5), xlim = c(start(gr), end(gr)),
         xlab = "", ylab = "")
    
    colourCount <- length(annoTrack)
    getPalette <- colorRampPalette(.alpha(brewer.pal(max(length(annoTrack), 3),
                                                     "Dark2"), 
        0.4))
    color <- getPalette(colourCount)
    bord <- "black"
    
    lapply(seq(along = annoTrack), function(ii) {
        jj <- length(annoTrack) + 1 - ii
        ir <- subsetByOverlaps(annoTrack[[ii]], gr)
        start(ir) <- pmax(start(ir), start(gr))
        end(ir) <- pmin(end(ir), end(gr))
        
        top.pos <- 4.1 - (4 - length(annoTrack))
        bot.pos <- 0 + (4 - length(annoTrack)) * 0.1
        jj.between <- (top.pos - bot.pos)/length(annoTrack)
        
        if (length(ir) > 0) {
            jj <- top.pos - (ii - 1) * jj.between
            
            arrows(start(ir), jj, end(ir), jj, code = 3, length = 0.05, 
                   angle = 90, col = .makeTransparent(color[ii], alpha = 185))
            if (is.null(labelCols)) {
                text((end(gr) + start(gr))/2, jj - 0.4, names(annoTrack)[[ii]])
            }else if(sum(labelCols %in% colnames(mcols(annoTrack[[ii]]))) > 0) {
                whichLabelCols <- match(labelCols, 
                                        colnames(mcols(annoTrack[[ii]])))
                notmiss <- !is.na(whichLabelCols)
                whichLabelCols <- whichLabelCols[notmiss]
                if (class(unlist(mcols(annoTrack[[ii]])[1, 
                                            whichLabelCols, drop = TRUE])) == 
                  "numeric") {
                  comps <- round(unlist(mcols(ir)[1, 
                                            whichLabelCols, drop = TRUE]), 
                    3)
                } else {
                  comps <- unlist(mcols(ir)[1, whichLabelCols, drop = TRUE])
                }
                if (sum(grepl("\\.", labelCols)) > 0) {
                  labelCols <- gsub("\\.", " ", labelCols)
                }
                Label <- paste0(names(annoTrack)[ii], ": ", 
                                paste0(labelCols[notmiss], 
                  "=", comps, collapse = ", "))
                text((end(gr) + start(gr))/2, jj - 0.4, Label, 
                     cex = 0.85, col = color[ii])
            } else {
                text((end(gr) + start(gr))/2, jj - 0.4, names(annoTrack)[[ii]])
            }
        }
    })
}
