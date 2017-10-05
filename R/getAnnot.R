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
#' under study. Use the function \code{builtin_genomes2()} (modified from
#' the annotatr
#' package) to see
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
#' @import TxDb.Hsapiens.UCSC.hg18.knownGene
#' @importFrom AnnotationDbi mappedkeys select
#' @importFrom GenomicFeatures transcripts cdsBy promoters exonsBy 
#' @importFrom GenomicFeatures fiveUTRsByTranscript threeUTRsByTranscript
#' @importFrom GenomicFeatures intronsByTranscript
#' @importFrom AnnotationHub AnnotationHub query
#' @importFrom GenomeInfoDb Seqinfo
#' @importFrom S4Vectors elementNROWS Rle queryHits
#' @importFrom readr read_tsv
#' 
#' @examples
#' 
#' # need to load annotatr package with "library(annotatr)" before 
#' # running this command to retrieve annotations
#' 
#' # get annotation information for hg18
#' annoTrack <- getAnnot("hg18")
#' 
#' 
getAnnot <- function(genomeName){
  if (!requireNamespace("annotatr", quietly = TRUE)) {
    message(paste0("annotatr could not be loaded. Please make sure it is ",
                "installed, or skip the annotation step and leave as NULL",
                "(default value) in plotDMRs."))
    return(NULL)
  }
  
  if (!genomeName %in% builtin_genomes2()){
    message(paste0("Genome ", genomeName, " is not supported by ",
                   "annotatr at this time"))
    return(NULL)
  }
  
  if (is.null(genomeName)){
    return(NULL)
  }else{
    annot_CpG = paste0(c(genomeName, "_cpgs"), collapse="")
    annot_genes = paste0(c(genomeName, "_genes_cds"), collapse="")
    
    # Build the annotations (a single GRanges object)
    
    # Download has a nonzero fail rate; allow up to 5 retries before
    # throwing an error
    
    for (attempt in 1:5){
      cpg = try( build_annotations2(genome = genomeName, 
                                            annotations = annot_CpG), 
                silent=TRUE)  
      if(!is(cpg, 'try-error')){
        message("Download of CpG annotation successful!")
        fail1 <- 0
        break;
      }else{
        message(paste0("Download of CpG annotation failed. ",
                       5-attempt, " attempts left"))
        fail1 <- 1
      }
    }
    
    for (attempt in 1:5){
      genes = try( build_annotations2(genome = genomeName, 
                                               annotations = annot_genes), 
                  silent=TRUE)  
      if(!is(genes, 'try-error')){
        message("Download of Gene annotation successful!")
        fail2 <- 0
        break;
      }else{
        message(paste0("Download of Gene annotation failed. ",
                       5-attempt, " attempts left"))
        fail2 <- 1
      }
    }
    
    if (fail1==0 & fail2 ==0){
      keep <- which(!is.na(genes$symbol))
      genes <- genes[keep,]
      cpg$type <- substr(cpg$type, 10, nchar(cpg$type))
      
      annot <- list("CpGs"=cpg,
                    "Exons"=genes)
      return(annot)
    }else{
      stop("Annotation could not be retrieved.")
    }
  }}

##############################################################################
# annotatr functions modified to include the hg18 genome
##############################################################################

#' A helper function to build CpG related annotations.
#'
#' Using the \code{AnnotationHub} package, extract CpG island track 
#' for the appropriate \code{genome} and construct the shores, shelves, and 
#' interCGI annotations as desired. ((This function is a modified version
#' of that found in the annotatr package, simply to include the hg18 version
#' of the genome for CpG categories))
#'
#' @param genome The genome assembly.
#' @param annotations A character vector with entries of the form 
#' \code{[genome]_cpg_{islands,shores,shelves,inter}}.
#'
#' @return A list of \code{GRanges} objects.
build_cpg_annots2 = function(genome = builtin_genomes2(), 
                            annotations = builtin_annotations2()) {
  # Ensure valid arguments
  genome = match.arg(genome)
  annotations = match.arg(annotations, several.ok = TRUE)
  
  annot_codes = data.frame(
    code = c(sprintf('%s_cpg_islands', genome),
             sprintf('%s_cpg_shores', genome),
             sprintf('%s_cpg_shelves', genome),
             sprintf('%s_cpg_inter', genome)),
    var = c('islands','shores','shelves','inter_cgi'),
    stringsAsFactors = FALSE)
  
  # Decide whether to use URL or AnnotationHub
  ah_genomes = c('hg19','mm9','rn5','rn4')
  if(genome == 'hg19' || genome == 'mm9' || genome == 'rn5' || genome == 'rn4') {
    use_ah = TRUE
  } else if (genome == 'hg38') {
    use_ah = FALSE
    con = 'http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cpgIslandExt.txt.gz'
  } else if (genome == 'hg18') {
    use_ah = FALSE
    con = 'http://hgdownload.cse.ucsc.edu/goldenpath/hg18/database/cpgIslandExt.txt.gz'
  } else if (genome == 'mm10') {
    use_ah = FALSE
    con = 'http://hgdownload.cse.ucsc.edu/goldenpath/mm10/database/cpgIslandExt.txt.gz'
  } else if (genome == 'rn6') {
    use_ah = FALSE
    con = 'http://hgdownload.cse.ucsc.edu/goldenpath/rn6/database/cpgIslandExt.txt.gz'
  } else {
    stop(sprintf('CpG features are not supported for genome %s', genome))
  }
  
  if(use_ah) {
    # Create AnnotationHub connection
    ah = AnnotationHub::AnnotationHub()
    
    # And do the query for available CpG Islands
    query = AnnotationHub::query(ah, c('CpG Islands'))
    
    # Determine the correct ID to extract data from AnnotationHub
    ID = row.names(GenomicRanges::mcols(query)[GenomicRanges::mcols(query)$genome == genome, ])
  }
  
  if(any(grepl('islands', annotations)) || any(grepl('shores', annotations)) ||
     any(grepl('shelves', annotations)) || any(grepl('inter', annotations))) {
    message('Building CpG islands...')
    ### Islands
    # Extract and sort the islands based on use_ah
    if(use_ah) {
      islands = ah[[ID]]
    } else {
      # Read from URL. There is surprisingly nothing in base that
      # does this as easily, so here we are with readr again.
      if (genome == "hg18"){
        islands_tbl = readr::read_tsv(con,
                                      col_names = c('chr','start','end'),
                                      col_types = 'cii-------')
      }else{
        islands_tbl = readr::read_tsv(con,
                                    col_names = c('chr','start','end'),
                                    col_types = '-cii-------')
      }
      # Convert to GRanges
      islands = GenomicRanges::GRanges(
        seqnames = islands_tbl$chr,
        ranges = IRanges::IRanges(start = islands_tbl$start, end = islands_tbl$end),
        strand = '*',
        seqinfo = GenomeInfoDb::Seqinfo(genome=genome)
      )
    }
    islands = GenomicRanges::sort(islands)
    
    
    # Rename the islands
    GenomicRanges::mcols(islands)$id = paste0('island:', seq_along(islands))
    
    # Add tx_name, gene_id, and symbol columns
    GenomicRanges::mcols(islands)$tx_id = NA
    GenomicRanges::mcols(islands)$gene_id = NA
    GenomicRanges::mcols(islands)$symbol = NA
    
    # Give it the correct type
    GenomicRanges::mcols(islands)$type = sprintf('%s_cpg_islands', genome)
    
    GenomicRanges::mcols(islands) = GenomicRanges::mcols(islands)[, c('id','tx_id','gene_id','symbol','type')]
    
    if(any(grepl('shores', annotations)) || 
       any(grepl('shelves', annotations)) || 
       any(grepl('inter', annotations))) {
      message('Building CpG shores...')
      ### Shores
      # Construct the shores based on:
      # upstream from the island start and downstream from the island end
      up_shores = GenomicRanges::flank(islands, width = 2000, 
                                       start = TRUE, both = FALSE)
      down_shores = GenomicRanges::flank(islands, width = 2000, 
                                         start = FALSE, both = FALSE)
      
      # Combine, sort, trim, and reduce combined up_shores and down_shores
      shores = c(up_shores, down_shores)
      shores = GenomicRanges::sort(shores)
      shores = GenomicRanges::trim(shores)
      shores = GenomicRanges::reduce(shores)
      
      # Remove islands from the shores
      shores = GenomicRanges::setdiff(shores, islands)
      
      # Rename the shores
      GenomicRanges::mcols(shores)$id = paste0('shore:', seq_along(shores))
      
      # Add tx_name, gene_id, and symbol columns
      GenomicRanges::mcols(shores)$tx_id = NA
      GenomicRanges::mcols(shores)$gene_id = NA
      GenomicRanges::mcols(shores)$symbol = NA
      
      # Give it the correct type
      GenomicRanges::mcols(shores)$type = sprintf('%s_cpg_shores', genome)
      
      GenomicRanges::mcols(shores) = GenomicRanges::mcols(shores)[, c('id','tx_id','gene_id','symbol','type')]
      
      if(any(grepl('shelves', annotations)) || any(grepl('inter', annotations))) {
        message('Building CpG shelves...')
        ### Shelves
        # Construct the shelves based on:
        # upstream from the up_shores start and downstream from the down_shores end
        up_shelves = GenomicRanges::flank(shores, width = 2000, 
                                          start = TRUE, both = FALSE)
        down_shelves = GenomicRanges::flank(shores, width = 2000, 
                                            start = FALSE, both = FALSE)
        
        # Combine, sort, trim, and reduce combined up_shelves and down_shelves
        shelves = c(up_shelves, down_shelves)
        shelves = GenomicRanges::sort(shelves)
        shelves = GenomicRanges::trim(shelves)
        shelves = GenomicRanges::reduce(shelves)
        
        # Remove islands and/or shores from the shelves
        shelves = GenomicRanges::setdiff(shelves, islands)
        shelves = GenomicRanges::setdiff(shelves, shores)
        
        # Rename the shelves
        GenomicRanges::mcols(shelves)$id = paste0('shelf:', seq_along(shelves))
        
        # Add gene_id, and symbol columns
        GenomicRanges::mcols(shelves)$tx_id = NA
        GenomicRanges::mcols(shelves)$gene_id = NA
        GenomicRanges::mcols(shelves)$symbol = NA
        
        # Give it the correct type
        GenomicRanges::mcols(shelves)$type = sprintf('%s_cpg_shelves', genome)
        
        GenomicRanges::mcols(shelves) = GenomicRanges::mcols(shelves)[, c('id','tx_id','gene_id','symbol','type')]
        
        if(any(grepl('inter', annotations))) {
          message('Building inter-CpG-islands...')
          ### interCGI
          # Take the union of all the objects so far
          extended_cgi = Reduce(
            function(x,y){
              GenomicRanges::union(x,y)
            },
            list(islands, shores, shelves)
          )
          
          # A quirk of GenomicRanges::gaps() on unstranded ranges gives the whole
          # chroms on the + and - strands in addition to the * gaps. Remove +/- gaps.
          inter_cgi = GenomicRanges::gaps(extended_cgi)
          inter_cgi = inter_cgi[GenomicRanges::strand(inter_cgi) == '*']
          
          # Rename the interCGI
          GenomicRanges::mcols(inter_cgi)$id = paste0('inter:', seq_along(inter_cgi))
          
          # Add gene_id, and symbol columns
          GenomicRanges::mcols(inter_cgi)$tx_id = NA
          GenomicRanges::mcols(inter_cgi)$gene_id = NA
          GenomicRanges::mcols(inter_cgi)$symbol = NA
          
          # Give it the correct type
          GenomicRanges::mcols(inter_cgi)$type = sprintf('%s_cpg_inter', genome)
          
          GenomicRanges::mcols(inter_cgi) = GenomicRanges::mcols(inter_cgi)[, c('id','tx_id','gene_id','symbol','type')]
        }
      }
    }
  }
  
  ### Put it all together
  mgets = annot_codes[annot_codes$code %in% annotations, 'var']
  cpgs = do.call('GRangesList', mget(mgets))
  names(cpgs) = annotations
  
  return(cpgs)
}

#' A function to build annotations from TxDb.* and AnnotationHub resources
#'
#' Create a \code{GRanges} object consisting of all the desired 
#' \code{annotations}. Supported annotation codes are listed by 
#' \code{builtin_annotations2()}. 
#' The basis for enhancer annotations are FANTOM5 data, 
#' the basis for CpG related annotations are CpG island 
#' tracks from \code{AnnotationHub}, and the basis for genic 
#' annotations are from the \code{TxDb.*} and \code{org.db} group of packages.
#' Modified from annotatr package to include hg18
#'
#' @param genome The genome assembly.
#' @param annotations A character vector of annotations to build. Valid annotation codes are listed with \code{builtin_annotations2()}. The "basicgenes" shortcut builds the following regions: 1-5Kb upstream of TSSs, promoters, 5UTRs, exons, introns, and 3UTRs. The "cpgs" shortcut builds the following regions: CpG islands, shores, shelves, and interCGI regions. NOTE: Shortcuts need to be appended by the genome, e.g. \code{hg19_basicgenes}.
#' Custom annotations whose names are of the form \code{[genome]_custom_[name]} should also be included. Custom annotations should be read in and converted to \code{GRanges} with \code{read_annotations()}. They can be for a \code{supported_genome()}, or for an unsupported genome.
#'
#' @return A \code{GRanges} object of all the \code{annotations} combined. The \code{mcols} are \code{id, tx_id, gene_id, symbol, type}. The \code{id} column is a unique name, the \code{tx_id} column is either a UCSC knownGene transcript ID (genic annotations) or a Ensembl transcript ID (lncRNA annotations), the \code{gene_id} is the Entrez ID, the \code{symbol} is the gene symbol from the \code{org.*.eg.db} mapping from the Entrez ID, and the \code{type} is of the form \code{[genome]_[type]_[name]}.
#'
#' @examples
#' # Example with hg19 gene promoters
#' annots = c('hg19_genes_promoters')
#' annots_gr = build_annotations2(genome = 'hg19', annotations = annots)
#'
#' # See vignette for an example with custom annotation
#'
#' @export
build_annotations2 = function(genome, annotations) {
  # Expand annotations first
  annotations = annotatr:::expand_annotations(annotations)
  
  # Collect built-in annotations
  hmm_annotations = grep('_chromatin_', annotations, value=TRUE)
  enh_annotations = grep('_enhancers_', annotations, value=TRUE)
  gene_annotations = grep('_genes_', annotations, value=TRUE)
  cpg_annotations = grep('_cpg_', annotations, value=TRUE)
  lncrna_annotations = grep('_lncrna_', annotations, value=TRUE)
  builtin_annotations = c(hmm_annotations, enh_annotations, 
                          gene_annotations, cpg_annotations, 
                          lncrna_annotations)
  
  # Check builtin_annotations
  if(length(builtin_annotations) > 0) {
    check_annotations2(builtin_annotations)
  }
  
  # Collect custom and AnnotationHub annotations
  other_annotations = setdiff(annotations, builtin_annotations)
  
  # Container for the annotations
  annots_grl = GenomicRanges::GRangesList()
  
  # Do the other_annotations first because we don't want to fail unexpectedly at the end
  if(length(other_annotations) > 0) {
    annots_grl = c(annots_grl, 
                   GenomicRanges::GRangesList(sapply(other_annotations, 
                            function(ca){annotatr::annotatr_cache$get(ca)})))
  }
  # Take the builtin_annotations piece by piece
  if(length(enh_annotations) != 0) {
    annots_grl = c(annots_grl, 
                   GenomicRanges::GRangesList(enhancers_fantom = 
          suppressWarnings(annotatr:::build_enhancer_annots(genome = genome))))
  }
  if(length(hmm_annotations) != 0) {
    annots_grl = c(annots_grl, 
                   GenomicRanges::GRangesList(chromatin = 
                 suppressWarnings(annotatr:::build_hmm_annots(genome = genome, 
                                          annotations = hmm_annotations))))
  }
  if(length(gene_annotations) != 0) {
    annots_grl = c(annots_grl, 
                suppressWarnings(build_gene_annots2(genome = genome,
                                      annotations = gene_annotations)))
  }
  if(length(cpg_annotations) != 0) {
    annots_grl = c(annots_grl, 
      suppressWarnings(build_cpg_annots2(genome = genome, 
                                    annotations = cpg_annotations)))
  }
  if(length(lncrna_annotations) != 0) {
    annots_grl = c(annots_grl, 
              GenomicRanges::GRangesList(lncrna_gencode = 
          suppressWarnings(annotatr:::build_lncrna_annots(genome = genome))))
  }
  
  return(unlist(annots_grl, use.names=FALSE))
}

#' Function returning supported TxDb.* genomes
#' modified from annotatr to include hg18
#'
#' @return A character vector of genomes for supported TxDb.* packages
#'
#' @examples
#' builtin_genomes2()
#'
#' @export
builtin_genomes2 = function() {
  return(c('dm3','dm6','hg19','hg38','mm9','mm10','rn4','rn5','rn6',
           'hg18'))
}

#' Function listing which annotations are available.
#' Modified from annotatr to include hg18
#'
#' This includes the shortcuts. The \code{expand_annotations()} function helps
#' handle the shortcuts.
#'
#' @return A character vector of available annotations.
#'
#' @examples
#' builtin_annotations2()
#'
#' @export
builtin_annotations2 = function() {
  # Create annotation code endings
  shortcut_ends = c('basicgenes','cpgs')
  
  # Gene codes
  gene_genomes = builtin_genomes2()
  gene_ends = c('1to5kb', 'promoters', 'cds', '5UTRs', 'exons', 'firstexons', 'introns', 'intronexonboundaries', 'exonintronboundaries', '3UTRs', 'intergenic')
  
  # CpG codes
  cpg_genomes = base::setdiff(builtin_genomes2(),c('dm3','dm6'))
  cpg_ends = c('islands', 'shores', 'shelves', 'inter')
  
  # Chromatin state codes
  # Remove numbers, and underscores, and take unique
  chromatin_recode = unique(annotatr:::reformat_hmm_codes(annotatr:::HMMCODES))
  
  chromatin_ends = apply(
    expand.grid(annotatr:::HMMCELLLINES, chromatin_recode, stringsAsFactors = FALSE),
    1, paste, collapse='-')
  
  chromatin_shortcut_ends = apply(
    expand.grid(annotatr:::HMMCELLLINES, 'chromatin', stringsAsFactors = FALSE),
    1, paste, collapse='-')
  
  # Create full annotation codes
  gene_codes = apply(
    expand.grid(gene_genomes, 'genes', gene_ends, stringsAsFactors = FALSE),
    1, paste, collapse='_')
  cpg_codes = apply(
    expand.grid(cpg_genomes, 'cpg', cpg_ends, stringsAsFactors= FALSE),
    1, paste, collapse='_')
  chromatin_codes = apply(
    expand.grid('hg19', 'chromatin', chromatin_ends, stringsAsFactors=FALSE),
    1, paste, collapse='_')
  
  enhancer_codes = c('hg19_enhancers_fantom','hg38_enhancers_fantom','mm9_enhancers_fantom','mm10_enhancers_fantom')
  lncrna_codes = c('hg19_lncrna_gencode','hg38_lncrna_gencode','mm10_lncrna_gencode')
  
  gene_shortcut_codes = apply(
    expand.grid(gene_genomes, 'basicgenes', stringsAsFactors = FALSE),
    1, paste, collapse='_')
  cpg_shortcut_codes = apply(
    expand.grid(cpg_genomes, 'cpgs', stringsAsFactors = FALSE),
    1, paste, collapse='_')
  chromatin_shortcut_codes = paste('hg19', chromatin_shortcut_ends, sep='_')
  
  # Create the big vector of supported annotations
  annots = c(gene_codes, cpg_codes, chromatin_codes, enhancer_codes, lncrna_codes,
             gene_shortcut_codes, cpg_shortcut_codes, chromatin_shortcut_codes)
  
  return(annots)
}


#' Function to check for valid annotations
#' Modified from annotatr to include hg18
#'
#' Gives errors if any annotations are not in builtin_annotations2() 
#' (and they are not in the required custom format), basicgenes are used,
#'  or the genome prefixes are not the same for all annotations.
#'
#' @param annotations A character vector of annotations 
#' possibly using the shortcuts
#' @return If all the checks on the annotations pass, 
#' returns NULL to allow code to move forward.
check_annotations2 = function(annotations) {
  # Pull out any custom annotations before checking
  custom_annotations = grep('custom', annotations, value = TRUE)
  annotations = base::setdiff(annotations, custom_annotations)
  
  # Check that the annotations are supported, tell the user which are unsupported
  if( !all(annotations %in% builtin_annotations2()) ) {
    unsupported = base::setdiff(annotations, builtin_annotations2())
    
    stop(sprintf('Error: "%s" is(are) not supported. See builtin_annotations2().',
                 paste(unsupported, collapse=', ')))
  }
  
  # Recombine annotations and custom_annotations or you get failure when
  # there are only custom annotations
  annotations = c(custom_annotations, annotations)
  
  genomes = sapply(annotations, function(a){
    unlist(strsplit(a, '_'))[1]
  }, USE.NAMES = FALSE)
  
  # Check for same genome on all annotations
  if( length(unique(genomes)) != 1 ){
    stop('Error: genome prefix on all annotations must be the same.')
  }
  
  return(NULL)
}

#' A helper function to build genic annotations.
#' modified from annotatr to include hg18
#'
#' Using the \code{TxDb.*} group of packages, construct genic annotations consisting of any combination of 1-5kb upstream of a TSS, promoters (< 1kb from TSS), 5UTRs, CDS, exons, first exons, introns, intron/exon and exon/intron boundaries, 3UTRs, and intergenic.
#'
#' @param genome The genome assembly.
#' @param annotations A character vector with entries of the form \code{[genome]_genes_{1to5kb,promoters,5UTRs,cds,exons,firstexons,introns,intronexonboundaries,exonintronboundaries,3UTRs,intergenic}}.
#'
#' @return A list of \code{GRanges} objects with unique \code{id} of the form \code{[type]:i}, \code{tx_id} being the UCSC knownGene transcript name, \code{gene_id} being the Entrez Gene ID, \code{symbol} being the gene symbol from the Entrez ID to symbol mapping in \code{org.db} for that species, and \code{type} being the annotation type.
build_gene_annots2 = function(genome = builtin_genomes2(), 
                             annotations = builtin_annotations2()) {
  # Ensure valid arguments
  genome = match.arg(genome)
  annotations = match.arg(annotations, several.ok = TRUE)
  
  annot_codes = data.frame(
    code = c(sprintf('%s_genes_promoters', genome),
             sprintf('%s_genes_1to5kb', genome),
             sprintf('%s_genes_cds', genome),
             sprintf('%s_genes_5UTRs', genome),
             sprintf('%s_genes_exons', genome),
             sprintf('%s_genes_firstexons', genome),
             sprintf('%s_genes_introns', genome),
             sprintf('%s_genes_intronexonboundaries', genome),
             sprintf('%s_genes_exonintronboundaries', genome),
             sprintf('%s_genes_3UTRs', genome),
             sprintf('%s_genes_intergenic', genome)),
    var = c('promoters_gr','onetofive_gr','cds_gr','fiveUTRs_gr','exons_gr',
            'firstexons_gr','introns_gr','intronexon_gr','exonintron_gr',
            'threeUTRs_gr','intergenic_gr'),
    stringsAsFactors = FALSE)
  
  # Load the appropriate TxDb.* library and get the txdb
  txdb_name = get_txdb_name2(genome)
  library(txdb_name, character.only = TRUE)
  txdb = get(txdb_name)
  
  # Get the org.XX.eg.db mapping from Entrez ID to gene symbol
  # First element returned is package name, second is eg2SYMBOL name
  orgdb_name = get_orgdb_name2(genome)
  library(sprintf('org.%s.eg.db', orgdb_name), character.only = TRUE)
  x = get(sprintf('org.%s.egSYMBOL', orgdb_name))
  mapped_genes = mappedkeys(x)
  eg2symbol = as.data.frame(x[mapped_genes])
  
  # Build the base transcripts
  tx_gr = transcripts(txdb, columns = c('TXID','GENEID','TXNAME'))
  # Create TSS GRanges for later use with intronexon boundaries
  tss_gr = GenomicRanges::GRanges(
    seqnames = seqnames(tx_gr),
    ranges = IRanges(start = start(tx_gr), end = start(tx_gr)),
    strand = '*'
  )
  GenomicRanges::mcols(tss_gr) = GenomicRanges::mcols(tx_gr)
  seqinfo(tss_gr) = seqinfo(tx_gr)
  
  # Create TES GRanges for later use with exonintron boundaries
  tes_gr = GenomicRanges::GRanges(
    seqnames = seqnames(tx_gr),
    ranges = IRanges(start = end(tx_gr), end = end(tx_gr)),
    strand = '*'
  )
  GenomicRanges::mcols(tes_gr) = GenomicRanges::mcols(tx_gr)
  seqinfo(tes_gr) = seqinfo(tx_gr)
  
  # Build tables to map TXID to TXNAME and GENEID
  id_maps = AnnotationDbi::select(txdb, keys = as.character(GenomicRanges::mcols(tx_gr)$TXID), columns = c('TXNAME','GENEID'), keytype = 'TXID')
  
  # Each annotation should be a GRanges object with the following mcols:
  # id, tx_id, gene_id, symbol, type
  
  if(any(grepl('promoters', annotations)) || 
     any(grepl('1to5kb', annotations)) || 
     any(grepl('intergenic', annotations))) {
    message('Building promoters...')
    ### promoters
    promoters_gr = GenomicFeatures::promoters(txdb, 
                                              upstream = 1000, downstream = 0)
    # Add Entrez ID, symbol, and type
    GenomicRanges::mcols(promoters_gr)$gene_id =
      id_maps[match(GenomicRanges::mcols(promoters_gr)$tx_id,
                    id_maps$TXID), 'GENEID']
    GenomicRanges::mcols(promoters_gr)$symbol =
      eg2symbol[match(GenomicRanges::mcols(promoters_gr)$gene_id, 
                      eg2symbol$gene_id), 'symbol']
    GenomicRanges::mcols(promoters_gr)$type =
      sprintf('%s_genes_promoters', genome)
    GenomicRanges::mcols(promoters_gr)$id =
      paste0('promoter:', seq_along(promoters_gr))
    
    GenomicRanges::mcols(promoters_gr) = 
      GenomicRanges::mcols(promoters_gr)[, c('id','tx_name','gene_id','symbol','type')]
    colnames(GenomicRanges::mcols(promoters_gr)) = c('id','tx_id','gene_id','symbol','type')
    
    if(any(grepl('1to5kb', annotations)) || any(grepl('intergenic', annotations))) {
      message('Building 1to5kb upstream of TSS...')
      ### 1-5kb
      onetofive_gr = GenomicRanges::flank(promoters_gr, width = 4000, start = TRUE, both = FALSE)
      onetofive_gr = GenomicRanges::trim(onetofive_gr)
      # Add Entrez ID, symbol, and type (all but type are inherited from promoters_gr)
      GenomicRanges::mcols(onetofive_gr)$id = paste0('1to5kb:', seq_along(onetofive_gr))
      GenomicRanges::mcols(onetofive_gr)$type = sprintf('%s_genes_1to5kb', genome)
      
      if(any(grepl('intergenic', annotations))) {
        message('Building intergenic...')
        ### intergenic
        genic_gr = c(GenomicRanges::granges(tx_gr), GenomicRanges::granges(promoters_gr), GenomicRanges::granges(onetofive_gr))
        GenomicRanges::strand(genic_gr) = '*'
        intergenic_gr = GenomicRanges::reduce(genic_gr)
        intergenic_gr = GenomicRanges::gaps(intergenic_gr)
        
        # A quirk in gaps gives the entire + and - strand of a chromosome, ignore those
        intergenic_gr = intergenic_gr[GenomicRanges::strand(intergenic_gr) == '*']
        
        GenomicRanges::mcols(intergenic_gr)$id = paste0('intergenic:', seq_along(intergenic_gr))
        GenomicRanges::mcols(intergenic_gr)$tx_id = NA
        GenomicRanges::mcols(intergenic_gr)$gene_id = NA
        GenomicRanges::mcols(intergenic_gr)$symbol = NA
        GenomicRanges::mcols(intergenic_gr)$type = sprintf('%s_genes_intergenic', genome)
      }
    }
  }
  
  if(any(grepl('cds', annotations))) {
    message('Building cds...')
    ### cds
    cds_grl = GenomicFeatures::cdsBy(txdb, by = 'tx', use.names = TRUE)
    # Create Rle of the tx_names
    cds_txname_rle = S4Vectors::Rle(names(cds_grl), S4Vectors::elementNROWS(cds_grl))
    cds_txname_vec = as.character(cds_txname_rle)
    # Unlist and add the tx_names
    cds_gr = unlist(cds_grl, use.names = FALSE)
    GenomicRanges::mcols(cds_gr)$tx_name = cds_txname_vec
    # Add Entrez ID, symbol, and type
    GenomicRanges::mcols(cds_gr)$gene_id = id_maps[match(GenomicRanges::mcols(cds_gr)$tx_name, id_maps$TXNAME), 'GENEID']
    GenomicRanges::mcols(cds_gr)$symbol = eg2symbol[match(GenomicRanges::mcols(cds_gr)$gene_id, eg2symbol$gene_id), 'symbol']
    GenomicRanges::mcols(cds_gr)$type = sprintf('%s_genes_cds', genome)
    GenomicRanges::mcols(cds_gr)$id = paste0('CDS:', seq_along(cds_gr))
    
    GenomicRanges::mcols(cds_gr) = GenomicRanges::mcols(cds_gr)[, c('id','tx_name','gene_id','symbol','type')]
    colnames(GenomicRanges::mcols(cds_gr)) = c('id','tx_id','gene_id','symbol','type')
  }
  
  if(any(grepl('5UTR', annotations))) {
    message('Building 5UTRs...')
    ### fiveUTRs
    fiveUTRs_grl = GenomicFeatures::fiveUTRsByTranscript(txdb, use.names = TRUE)
    # Create Rle of the tx_names
    fiveUTRs_txname_rle = S4Vectors::Rle(names(fiveUTRs_grl), S4Vectors::elementNROWS(fiveUTRs_grl))
    fiveUTRs_txname_vec = as.character(fiveUTRs_txname_rle)
    # Unlist and add the tx_names
    fiveUTRs_gr = unlist(fiveUTRs_grl, use.names = FALSE)
    GenomicRanges::mcols(fiveUTRs_gr)$tx_name = fiveUTRs_txname_vec
    # Add Entrez ID, symbol, and type
    # NOTE: here we match on the tx_name because the tx_id is not given
    GenomicRanges::mcols(fiveUTRs_gr)$gene_id = id_maps[match(GenomicRanges::mcols(fiveUTRs_gr)$tx_name, id_maps$TXNAME), 'GENEID']
    GenomicRanges::mcols(fiveUTRs_gr)$symbol = eg2symbol[match(GenomicRanges::mcols(fiveUTRs_gr)$gene_id, eg2symbol$gene_id), 'symbol']
    GenomicRanges::mcols(fiveUTRs_gr)$type = sprintf('%s_genes_5UTRs', genome)
    GenomicRanges::mcols(fiveUTRs_gr)$id = paste0('5UTR:', seq_along(fiveUTRs_gr))
    
    GenomicRanges::mcols(fiveUTRs_gr) = GenomicRanges::mcols(fiveUTRs_gr)[, c('id','tx_name','gene_id','symbol','type')]
    colnames(GenomicRanges::mcols(fiveUTRs_gr)) = c('id','tx_id','gene_id','symbol','type')
    
  }
  
  if(any(grepl('3UTR', annotations))) {
    message('Building 3UTRs...')
    ### threeUTRs
    threeUTRs_grl = GenomicFeatures::threeUTRsByTranscript(txdb, use.names = TRUE)
    # Create Rle of the tx_names
    threeUTRs_txname_rle = S4Vectors::Rle(names(threeUTRs_grl), S4Vectors::elementNROWS(threeUTRs_grl))
    threeUTRs_txname_vec = as.character(threeUTRs_txname_rle)
    # Unlist and add the tx_names
    threeUTRs_gr = unlist(threeUTRs_grl, use.names = FALSE)
    GenomicRanges::mcols(threeUTRs_gr)$tx_name = threeUTRs_txname_vec
    # NOTE: here we match on the tx_name because the tx_id is not given
    GenomicRanges::mcols(threeUTRs_gr)$gene_id = id_maps[match(GenomicRanges::mcols(threeUTRs_gr)$tx_name, id_maps$TXNAME), 'GENEID']
    GenomicRanges::mcols(threeUTRs_gr)$symbol = eg2symbol[match(GenomicRanges::mcols(threeUTRs_gr)$gene_id, eg2symbol$gene_id), 'symbol']
    GenomicRanges::mcols(threeUTRs_gr)$type = sprintf('%s_genes_3UTRs', genome)
    GenomicRanges::mcols(threeUTRs_gr)$id = paste0('3UTR:', seq_along(threeUTRs_gr))
    
    GenomicRanges::mcols(threeUTRs_gr) = GenomicRanges::mcols(threeUTRs_gr)[, c('id','tx_name','gene_id','symbol','type')]
    colnames(GenomicRanges::mcols(threeUTRs_gr)) = c('id','tx_id','gene_id','symbol','type')
  }
  
  if(any(grepl('exon', annotations)) || any(grepl('intron', annotations))) {
    
    message('Building exons...')
    ### exons
    exons_grl = GenomicFeatures::exonsBy(txdb, by = 'tx', use.names = TRUE)
    # Create Rle of the tx_names
    exons_txname_rle = S4Vectors::Rle(names(exons_grl), S4Vectors::elementNROWS(exons_grl))
    exons_txname_vec = as.character(exons_txname_rle)
    # Unlist and add the tx_names
    exons_gr = unlist(exons_grl, use.names = FALSE)
    GenomicRanges::mcols(exons_gr)$tx_name = exons_txname_vec
    # Add Entrez ID, symbol, and type
    GenomicRanges::mcols(exons_gr)$gene_id = id_maps[match(GenomicRanges::mcols(exons_gr)$tx_name, id_maps$TXNAME), 'GENEID']
    GenomicRanges::mcols(exons_gr)$symbol = eg2symbol[match(GenomicRanges::mcols(exons_gr)$gene_id, eg2symbol$gene_id), 'symbol']
    GenomicRanges::mcols(exons_gr)$type = sprintf('%s_genes_exons', genome)
    GenomicRanges::mcols(exons_gr)$id = paste0('exon:', seq_along(exons_gr))
    
    # This needs to be here before we remove the exon_rank mcol in exons_gr
    if(any(grepl('firstexons', annotations))) {
      message('Building first exons...')
      ### first exons
      firstexons_gr = exons_gr[sapply(exons_gr$exon_rank, function(er){1 %in% er})]
      # NOTE: The mcol() contents are CharacterLists and IntegerLists, which requires a different approach from previous
      GenomicRanges::mcols(firstexons_gr)$type = sprintf('%s_genes_firstexons', genome)
      GenomicRanges::mcols(firstexons_gr)$id = paste0('firstexon:', seq_along(firstexons_gr))
      
      GenomicRanges::mcols(firstexons_gr) = GenomicRanges::mcols(firstexons_gr)[, c('id','tx_name','gene_id','symbol','type')]
      colnames(GenomicRanges::mcols(firstexons_gr)) = c('id','tx_id','gene_id','symbol','type')
    }
    
    GenomicRanges::mcols(exons_gr) = GenomicRanges::mcols(exons_gr)[, c('id','tx_name','gene_id','symbol','type')]
    colnames(GenomicRanges::mcols(exons_gr)) = c('id','tx_id','gene_id','symbol','type')
    
    message('Building introns...')
    ### introns
    introns_grl = GenomicFeatures::intronsByTranscript(txdb, use.names = TRUE)
    # Create Rle of the tx_names
    introns_txname_rle = S4Vectors::Rle(names(introns_grl), S4Vectors::elementNROWS(introns_grl))
    introns_txname_vec = as.character(introns_txname_rle)
    # Unlist and add the tx_names
    introns_gr = unlist(introns_grl, use.names = FALSE)
    GenomicRanges::mcols(introns_gr)$tx_name = introns_txname_vec
    # NOTE: here we match on the tx_name because the tx_id is not given
    GenomicRanges::mcols(introns_gr)$gene_id = id_maps[match(GenomicRanges::mcols(introns_gr)$tx_name, id_maps$TXNAME), 'GENEID']
    GenomicRanges::mcols(introns_gr)$symbol = eg2symbol[match(GenomicRanges::mcols(introns_gr)$gene_id, eg2symbol$gene_id), 'symbol']
    GenomicRanges::mcols(introns_gr)$type = sprintf('%s_genes_introns', genome)
    GenomicRanges::mcols(introns_gr)$id = paste0('intron:', seq_along(introns_gr))
    
    GenomicRanges::mcols(introns_gr) = GenomicRanges::mcols(introns_gr)[, c('id','tx_name','gene_id','symbol','type')]
    colnames(GenomicRanges::mcols(introns_gr)) = c('id','tx_id','gene_id','symbol','type')
    
    if(any(grepl('intronexonboundaries', annotations))) {
      message('Building intron exon boundaries...')
      ### Intron/exon boundary
      # Intron to exon transition will be the starts of exons_gr for + strand
      # and ends of exons_gr for - strand
      split_exons = split(exons_gr, GenomicRanges::strand(exons_gr))
      
      intronexon_plus_gr = GenomicRanges::GRanges(
        seqnames = seqnames(split_exons[['+']]),
        ranges = IRanges(start = start(split_exons[['+']]), end = start(split_exons[['+']])),
        strand = GenomicRanges::strand(split_exons[['+']]))
      GenomicRanges::mcols(intronexon_plus_gr) = GenomicRanges::mcols(split_exons[['+']])
      
      intronexon_minus_gr = GenomicRanges::GRanges(
        seqnames = seqnames(split_exons[['-']]),
        ranges = IRanges(start = end(split_exons[['-']]), end = end(split_exons[['-']])),
        strand = GenomicRanges::strand(split_exons[['-']]))
      GenomicRanges::mcols(intronexon_minus_gr) = GenomicRanges::mcols(split_exons[['-']])
      
      intronexon_gr = sort(c(intronexon_plus_gr, intronexon_minus_gr))
      seqinfo(intronexon_gr) = seqinfo(exons_gr)
      
      # Need to remove those boundaries that are TSSs
      tss_idx = unique(S4Vectors::queryHits(GenomicRanges::findOverlaps(intronexon_gr, tss_gr)))
      intronexon_gr = intronexon_gr[-tss_idx]
      
      # Expand 200bp up and down
      intronexon_gr = GenomicRanges::flank(intronexon_gr, width = 200, both = TRUE)
      GenomicRanges::mcols(intronexon_gr)$type = sprintf('%s_genes_intronexonboundaries', genome)
      GenomicRanges::mcols(intronexon_gr)$id = paste0('intronexonboundary:', seq_along(intronexon_gr))
      
      GenomicRanges::mcols(intronexon_gr) = GenomicRanges::mcols(intronexon_gr)[, c('id','tx_id','gene_id','symbol','type')]
    }
    
    if(any(grepl('exonintronboundaries', annotations))) {
      message('Building exon intron boundaries...')
      ### Exon/intron boundary
      if(!exists('split_exons')) {
        split_exons = split(exons_gr, GenomicRanges::strand(exons_gr))
      }
      
      # Exon to intron transition will be the ends of exons_gr for + strand
      # and starts of exons_gr for - strand
      exonintron_plus_gr = GenomicRanges::GRanges(
        seqnames = seqnames(split_exons[['+']]),
        ranges = IRanges(start = end(split_exons[['+']]), end = end(split_exons[['+']])),
        strand = GenomicRanges::strand(split_exons[['+']]))
      GenomicRanges::mcols(exonintron_plus_gr) = GenomicRanges::mcols(split_exons[['+']])
      
      exonintron_minus_gr = GenomicRanges::GRanges(
        seqnames = seqnames(split_exons[['-']]),
        ranges = IRanges(start = start(split_exons[['-']]), end = start(split_exons[['-']])),
        strand = GenomicRanges::strand(split_exons[['-']]))
      GenomicRanges::mcols(exonintron_minus_gr) = GenomicRanges::mcols(split_exons[['-']])
      
      exonintron_gr = sort(c(exonintron_plus_gr, exonintron_minus_gr))
      seqinfo(exonintron_gr) = seqinfo(exons_gr)
      
      # Need to remove those boundaries that are TSSs
      tss_idx = unique(S4Vectors::queryHits(GenomicRanges::findOverlaps(exonintron_gr, tss_gr)))
      exonintron_gr = exonintron_gr[-tss_idx]
      
      # Expand 200bp up and down
      exonintron_gr = GenomicRanges::flank(exonintron_gr, width = 200, both = TRUE)
      GenomicRanges::mcols(exonintron_gr)$type = sprintf('%s_genes_exonintronboundaries', genome)
      GenomicRanges::mcols(exonintron_gr)$id = paste0('exonintronboundary:', seq_along(exonintron_gr))
      
      GenomicRanges::mcols(exonintron_gr) = GenomicRanges::mcols(exonintron_gr)[, c('id','tx_id','gene_id','symbol','type')]
    }
  }
  
  ### Put it all together
  mgets = annot_codes[annot_codes$code %in% annotations, 'var']
  genes = do.call('GRangesList', mget(mgets))
  names(genes) = annotations
  
  return(genes)
}

#' Function to get correct TxDb.* package name based on genome
#' modified from annotatr to include hg18
#'
#' @param genome A string giving the genome assembly.
#'
#' @return A string giving the name of the correct TxDb.* package name based on \code{genome}.
get_txdb_name2 = function(genome = builtin_genomes2()) {
  # Ensure valid arguments
  genome = match.arg(genome)
  
  db = grep(genome, c(annotatr:::TXDBS,
                      "TxDb.Hsapiens.UCSC.hg18.knownGene"), value = TRUE)
  
  return(db)
}


ORGDBS2 = data.frame(
  genome = c('dm3','dm6','hg19','hg38','mm9','mm10','rn4','rn5','rn6', 'hg18'),
  org = c('Dm','Dm','Hs','Hs','Mm','Mm','Rn','Rn','Rn', 'Hs'),
  stringsAsFactors = FALSE)

#' Function to get correct org.* package name based on genome
#' modified from annotatr to include hg18
#'
#' @param genome A string giving the genome assembly.
#'
#' @return A string giving the correct org for org.db packages. e.g. hg19 -> Hs.
get_orgdb_name2 = function(genome = builtin_genomes2()) {
  # Ensure valid arguments
  genome = match.arg(genome)
  # org.* family of packages
  
  org = ORGDBS2[ORGDBS2$genome == genome, 'org']
  
  return(org)
}
