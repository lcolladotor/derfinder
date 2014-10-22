#' Change naming style for a set of sequence names
#'
#' If available, use the information from GenomeInfoDb for your species of 
#' interest to map the sequence names from the style currently used to another 
#' valid style. For example, for Homo sapiens map '2' (NCBI style) to 'chr2' 
#' (UCSC style). If the information from GenomeInfoDb is not available, the 
#' original sequence names will be returned.
#'
#' @param seqnames A character vector with the sequence names.
#' @param style A single character vector specifying the naming style to use
#' for renaming the sequence names. 
#' @param species A single character vector with the species of interest: it has
#' to match the valid species names supported in GenomeInfoDb. See
#' \code{names(GenomeInfoDb::genomeStyles())}. If \code{species = NULL}, 
#' a guess will be made using the available information in GenomeInfoDb.
#' @param currentStyle A single character vector with the currently used
#' naming style. If \code{NULL}, a guess will be made from the naming styles
#' supported by \code{species}.
#' @param ... Arguments passed to other methods and/or advanced arguments.
#'
#' @return A vector of sequence names using the specified naming \code{style}.
#'
#' @author L. Collado-Torres
#' @export
#' @aliases extended_map_seqlevels
#'
#' @details This function is inspired from \link[GenomeInfoDb]{mapSeqlevels}
#' with the difference that it will return the original sequence names if
#' the species, current naming style or target naming style are not supported
#' in GenomeInfoDb.
#'
#' @examples
#'
#'
#' ## Without guessing any information
#' extendedMapSeqlevels('2', 'UCSC', 'Homo sapiens', 'NCBI')
#' 
#' ## Guessing the current naming style
#' extendedMapSeqlevels('2', 'UCSC', 'Homo sapiens')
#'
#' ## Guess species and current style
#' extendedMapSeqlevels('2', 'NCBI')
#'
#' ## Guess species while supplying the current style.
#' ## Probably an uncommon use-case
#' extendedMapSeqlevels('2', 'NCBI', currentStyle = 'TAIR10')
#'
#' ## Sequence names are unchanged when using an unsupported species
#' extendedMapSeqlevels('seq2', 'NCBI', 'toyOrganism')
#' 
#' \dontrun{
#' ## Set global species and style option
#' options('chrsStyle' = 'UCSC')
#' options('species' = 'homo_sapiens')
#' 
#' ## Run using global options
#' extendedMapSeqlevels('2')
#' }

extendedMapSeqlevels <- function(seqnames, style = getOption('chrsStyle',
    'UCSC'), species = getOption('species', 'homo_sapiens'), 
    currentStyle = NULL, ...) {
    
#' @param verbose If \code{TRUE} basic status updates will be printed along the 
#' way.
    verbose <- .advanced_argument('verbose', TRUE, ...)
    
    if(missing(style)) {
#' @param chrsStyle The naming style of the chromosomes. By default, UCSC. See 
#' \link[GenomeInfoDb]{seqlevelsStyle}.   
        style <- .advanced_argument('chrsStyle', 'UCSC', ...)
    }
    if(is.null(style)) {
        ## If style is NULL return seqnames un-changed
        return(seqnames)
    }
    
    ## Check inputs
    if(!is.character(seqnames))
        seqnames <- is.character(seqnames)
    stopifnot(is.character(style) & length(style) == 1)
    if(!is.null(species)) stopifnot(is.character(species) & length(species) == 1)
    if(!is.null(currentStyle)) stopifnot(is.character(currentStyle) & length(currentStyle) == 1)
    
    guessed <- is.null(species) | is.null(currentStyle)

    ## Guess species
    if(is.null(species)) {
        species <- GenomeInfoDb:::.guessSpeciesStyle(seqnames)[1]
        if(is.na(species)) {
            if(verbose)
                message("extendedMapSeqlevels: the 'seqnames' you supplied are currently not supported in GenomeInfoDb. Consider adding your genome by following the information at http://www.bioconductor.org/packages/release/bioc/vignettes/GenomeInfoDb/inst/doc/Accept-organism-for-GenomeInfoDb.pdf")
            return(seqnames)
        }
        
    }
    
    ## Format species name
    species <- tolower(gsub(" ", "_", species))
    
    ## Extract valid mappings
    supported <- GenomeInfoDb:::.supportedSeqnameMappings()
    
    ## Check species is supported
    i <- which(tolower(names(supported)) == species)
    if(length(i) != 1){
        if(verbose)
            message(paste("extendedMapSeqlevels: the 'species'", species, "is currently not supported in GenomeInfoDb. Check valid 'species' by running names(GenomeInfoDb::genomeStyles()). If it's not present, consider adding your genome by following the information at http://www.bioconductor.org/packages/release/bioc/vignettes/GenomeInfoDb/inst/doc/Accept-organism-for-GenomeInfoDb.pdf"))
        return(seqnames)
    }
    ## Valid mapping for the species
    mapping <- supported[[i]]
    
    ## Check style is supported
    j <- which(tolower(colnames(mapping)) == tolower(style))
    if(length(j) != 1) {
        if(verbose)
            message(paste("extendedMapSeqlevels: the 'style'", style, "is currently not supported for the 'species'", species, "in GenomeInfoDb. Check valid naming styles by running GenomeInfoDb::genomeStyles(species). If it's not present, consider adding your genome by following the information at http://www.bioconductor.org/packages/release/bioc/vignettes/GenomeInfoDb/inst/doc/Accept-organism-for-GenomeInfoDb.pdf"))
        return(seqnames)
    }
    
    ## Guess current style for the current species
    if(is.null(currentStyle)) {
        if(ncol(mapping) < 2) {
            stop("extendedMapSeqlevels: the mapping for your specified 'species' is incomplete: there should be at least 2 different styles.")
        }
        ## Proceed as if using mapSeqlevels(best.only = TRUE)
        possible <- sapply(mapping, function(styleNames) {
            sum(tolower(seqnames) %in% tolower(styleNames)) })
        if(max(possible) == 0) {
            if(verbose)
                message("extendedMapSeqlevels: the current naming style could not be guessed. Consider supplying 'currentStyle'.")
            return(seqnames)
        } else if (max(possible) < length(seqnames)) {
            warning("extendedMapSeqlevels: not all the 'seqnames' match the best current naming style guess.")
        }
        currentStyle <- colnames(mapping)[ which.max(possible) ]
    }
    ## Find current naming map
    k <- which(tolower(colnames(mapping)) == tolower(currentStyle))
    
    ## Check currentStyle if it was supplied
    if(!is.null(currentStyle)) {
        if(length(k) != 1) {
            if(verbose)
                message(paste("extendedMapSeqlevels: the 'currentStyle'", currentStyle, "is currently not supported for the 'species'", species, "in GenomeInfoDb. Check valid naming styles by running GenomeInfoDb::genomeStyles(species). If it's not present, consider adding your genome by following the information at http://www.bioconductor.org/packages/release/bioc/vignettes/GenomeInfoDb/inst/doc/Accept-organism-for-GenomeInfoDb.pdf"))
            return(seqnames)
        }
    } 
    
    ## Is the target style the same as the current style?
    if(tolower(style) == tolower(currentStyle)) {
        return(seqnames)
    }
    
    ## Make the map
    seqnames <- mapping[match(seqnames, mapping[, k]), j]
    if(verbose & guessed)
        message(paste('extendedMapSeqlevels: sequence names mapped from', currentStyle, 'to', style, 'for species', species))
    return(seqnames)
}

#' @export
extended_map_seqlevels <- extendedMapSeqlevels
