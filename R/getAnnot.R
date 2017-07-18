#' Retrieve annotation information 
#' 
#' Uses the \code{annotatr} package to retrieve annotation information (
#' CpG category and gene coding sequences) for the \code{annoTrack} argument
#' of \code{\link{plotDMRs}}. Allows for 5 
#' re-tries if download fails (to allow for a spotty internet connection).
#' 
#' @param genomeName a character object that indicates which organism is 
#' under study. Use the function \code{annotatr:::builtin_genomes()} to see
#' a character vector of available genome names to choose from (see 
#' \code{annotatr} documentation for more details).
#' 
#' @export
getAnnot <- function(genomeName){
  if (is.null(genomeName)){
    return(NULL)
  }else{
    annot_CpG = paste0(c(genomeName, "_cpgs"), collapse="")
    annot_genes = paste0(c(genomeName, "_genes_cds"), collapse="")
    
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
    
    for (attempt in 1:5){
      genes = try(build_annotations(genome = genomeName, annotations = annot_genes), 
                  silent=TRUE)  
      if(!is(genes, 'try-error')){
        message("Download of Gene annotation successful!")
        break;
      }else{
        message(paste0("Download of Gene annotation failed. ",
                       5-attempt, " attempts left"))
      }
    }
    
    keep <- which(!is.na(genes$symbol))
    genes <- genes[keep,]
    cpg$type <- substr(cpg$type, 10, nchar(cpg$type))
    
    annot <- list("CpGs"=cpg,
                  "Exons"=genes)
    return(annot)
  }}
