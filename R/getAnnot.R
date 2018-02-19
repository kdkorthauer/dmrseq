#' Retrieve annotation information 
#' 
#' Uses the \code{annotatr} package to retrieve annotation information (
#' CpG category and gene coding sequences) for the \code{annoTrack} argument
#' of \code{\link{plotDMRs}}. Allows for 5 
#' re-tries if download fails (to allow for a spotty internet connection).
#' 
#' @details Note that this package needs to attach the \code{annotatr} package, 
#' and will
#' return NULL if this cannot be done. You can still use the 
#' \code{\link{plotDMRs}} function without this optional annotation step, 
#' just by leaving the \code{annoTrack} argument as NULL.
#' 
#' @param genomeName a character object that indicates which organism is 
#' under study. Use the function \code{builtin_genomes()} to see
#' a character vector of available genome names to choose from (see 
#' \code{annotatr} documentation for more details).
#' 
#' @return a \code{GRangesList} object with two elements returned
#' by \code{\link{getAnnot}}. The first
#' contains CpG category information in the first element (optional)
#' coding gene sequence information in the second element (optional).
#' At least one of these elements needs to be non-null in order for 
#' any annotation to be plotted, but it is not necessary to contain
#' both.
#' 
#' @export
#' 
#' @import annotatr
#' @importFrom AnnotationHub AnnotationHub query
#' @importFrom rtracklayer liftOver
#' @importFrom GenomeInfoDb genome
#' 
#' @examples
#' 
#' # get annotation information for hg19
#' annoTrack <- getAnnot('hg19')
#' 
#' 
getAnnot <- function(genomeName) {
    requireNamespace("annotatr")
    liftTo <- NULL
    if(genomeName == 'hg18'){
        message("Genome ", genomeName, " will be built by lifting over ",
                "hg19 annotations from annotatr")
        liftTo <- 'hg18'
        genomeName <- 'hg19'
    }else if (!genomeName %in% annotatr::builtin_genomes()) {
        message(paste0("Genome ", genomeName, " is not supported by ", 
                       "annotatr at this time"))
        return(NULL)
    }
    
    if (is.null(genomeName)) {
        return(NULL)
    } else {
        annot_CpG <- paste0(c(genomeName, "_cpgs"), collapse = "")
        annot_genes <- paste0(c(genomeName, "_genes_cds"), collapse = "")
        
        # Build the annotations (a single GRanges object)
        
        # annotatr downloads files for genomes that aren't natively supported
        # (e.g. hg38)
        # Download has a nonzero fail rate; allow up to 5 retries before 
        # throwing an error
        
        for (attempt in 1:5) {
            cpg <- try(annotatr::build_annotations(genome = genomeName, 
                                          annotations = annot_CpG), 
                silent = TRUE)
            if (!is(cpg, "try-error")) {
                message("Download of CpG annotation successful!")
                fail1 <- 0
                break
            } else {
                message(paste0("Download of CpG annotation failed. ",
                  5 - attempt, 
                  " attempts left"))
                fail1 <- 1
            }
        }
        
        for (attempt in 1:5) {
            genes <- try(annotatr::build_annotations(genome = genomeName,
                                            annotations = annot_genes), 
                silent = TRUE)
            if (!is(genes, "try-error")) {
                message("Download of Gene annotation successful!")
                fail2 <- 0
                break
            } else {
                message(paste0("Download of Gene annotation failed. ",
                  5 - attempt, 
                  " attempts left"))
                fail2 <- 1
            }
        }
        

        if (fail1 == 0 & fail2 == 0) {
            # lift over hg19 coordinates to hg18 using annotationHub
            if (!is.null(liftTo)){
              ah = AnnotationHub()
              chainfiles <- query(ah , c(genomeName, liftTo, "chainfile"))
              cf <- which(grepl(paste0(genomeName, "To", liftTo), 
                                chainfiles$title, 
                                ignore.case=TRUE))
              if (length(cf) == 0){
                stop("LiftOver from ", genomeName, " to ", liftTo, 
                     " was unsucccessful")
              }else if(length(cf) > 1){
                # take the first matching chain if more than one
                cf <- cf[1]
              }
              
              chain <- chainfiles[[names(chainfiles)[cf]]]
              
              cpg.new <- unlist(liftOver(cpg, chain))
              genes.new <- unlist(liftOver(genes, chain))
              
              GenomeInfoDb::genome(cpg) <- liftTo
              GenomeInfoDb::genome(genes) <- liftTo
            }
           
            keep <- which(!is.na(genes$symbol))
            genes <- genes[keep, ]
            cpg$type <- substr(cpg$type, 10, nchar(cpg$type))
            
            annot <- list(CpGs = cpg, Exons = genes)
            return(annot)
        } else {
            stop("Annotation could not be retrieved.")
        }
    }
}
