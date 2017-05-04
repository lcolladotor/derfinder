#' Change naming style for a set of sequence names
#'
#' If available, use the information from GenomeInfoDb for your species of 
#' interest to map the sequence names from the style currently used to another 
#' valid style. For example, for Homo sapiens map '2' (NCBI style) to 'chr2' 
#' (UCSC style). If the information from GenomeInfoDb is not available, the 
#' original sequence names will be returned. To disable this functionality
#' specify the same \code{style} and \param{currentStyle}.
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
#' Advanced arguments:
#' \describe{
#' \item{verbose }{ If \code{TRUE} basic status updates will be printed along 
#' the way.}
#' \item{chrsStyle }{ The naming style of the chromosomes. By default, 
#' \code{UCSC}. See \link[GenomeInfoDb]{seqlevelsStyle}.}
#' }
#'
#' @return A vector of sequence names using the specified naming \code{style}.
#'
#' @author L. Collado-Torres
#' @export
#'
#' @details This function is inspired from \link[GenomeInfoDb]{mapSeqlevels}
#' with the difference that it will return the original sequence names if
#' the species, current naming style or target naming style are not supported
#' in GenomeInfoDb.
#'
#' If you want to disable this function, specify the same \code{style} and
#' \code{currentStyle}. From other functions in derfinder that pass the
#' \code{...} argument to this function, use \code{currentStyle = 'UCSC'}
#' to match the \code{style}'s default. This can be useful when working
#' with organisms that are absent from GenomeInfoDb as documented in
#' \url{https://support.bioconductor.org/p/95521/}.
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
    
# @param verbose If \code{TRUE} basic status updates will be printed along the 
# way.
    verbose <- .advanced_argument('verbose', TRUE, ...)
    
    if(missing(style)) {
# @param chrsStyle The naming style of the chromosomes. By default, UCSC. See 
# \link[GenomeInfoDb]{seqlevelsStyle}.   
        style <- .advanced_argument('chrsStyle', 'UCSC', ...)
    }
    if(is.null(style)) {
        ## If style is NULL return seqnames un-changed
        return(seqnames)
    }
    
    ## Check inputs
    if(!is.character(seqnames))
        seqnames <- as.character(seqnames)
    stopifnot(is.character(style) & length(style) == 1)
    if(!is.null(species)) stopifnot(is.character(species) & length(species) == 1)
    if(!is.null(currentStyle)) stopifnot(is.character(currentStyle) & length(currentStyle) == 1)
    
    ## Was the species and/or the currentStyle guessed?
    guessedSpecies <- is.null(species)
    guessedCurrent <- is.null(currentStyle)
    guessed <- guessedSpecies | guessedCurrent
    if(guessed) {
        guesses <- GenomeInfoDb:::.guessSpeciesStyle(seqnames)
        if(any(is.na(guesses))) {
            if(verbose)
                message("extendedMapSeqlevels: the 'seqnames' you supplied are currently not supported in GenomeInfoDb. Consider adding your genome by following the information at http://www.bioconductor.org/packages/release/bioc/vignettes/GenomeInfoDb/inst/doc/Accept-organism-for-GenomeInfoDb.pdf")
            return(seqnames)
        }
    }
    
    ## Extract valid mappings
    supported <- GenomeInfoDb:::.supportedSeqnameMappings()
    names(supported) <- tolower(names(supported))

    ## Guess species and/or currentStyle
    if(guessedSpecies) {
        
        if(guessedCurrent) {
            ## Guess current style
            best <- .selectBestGuess(seqnames, guesses$species, guesses$style,
                supported)
            currentStyle <- best$style
        } else {
            ## Check current style
            idx <- grep(currentStyle, guesses$style)
            if(length(idx) == 0) {
                best$species <- NA
            } else {
                best < .selectBestGuess(seqnames, guesses$species[idx],
                    guesses$style[idx], supported)
            }           
            
        }
        species <- best$species
        
        ## Was able to guess the current species?
        if(is.na(species)) {
            if(verbose)
                message("extendedMapSeqlevels: the species could not be guessed. Consider supplying 'species'.")
            return(seqnames)
        }
            
    } else {
        ## Format species name
        species <- tolower(gsub(" ", "_", species))
        
        ## Check species is supported
        if(length(which(names(supported) == species)) != 1){
            if(verbose)
                message(paste("extendedMapSeqlevels: the 'species'", species, "is currently not supported in GenomeInfoDb. Check valid 'species' by running names(GenomeInfoDb::genomeStyles()). If it's not present, consider adding your genome by following the information at http://www.bioconductor.org/packages/release/bioc/vignettes/GenomeInfoDb/inst/doc/Accept-organism-for-GenomeInfoDb.pdf"))
            return(seqnames)
        }
        
        if (guessedCurrent) {
            ## Guess current style
            idx <- which(tolower(gsub(" ", "_", guesses$species)) == species)
            if(length(idx) == 0) {
                currentStyle <- NA
            } else {
                currentStyle <- .selectBestGuess(seqnames, guesses$species[idx],
                    guesses$style[idx], supported)$style
            }   
        }
    } 
    
    ## Was able to guesss the current style?
    if(is.na(currentStyle) & guessedCurrent) {
        if(verbose)
            message("extendedMapSeqlevels: the current naming style could not be guessed. Consider supplying 'currentStyle'.")
        return(seqnames)
    }

    ## Valid mapping for the species
    mapping <- supported[[which(names(supported) == species)]]
    
    ## Check style is supported
    j <- grep(tolower(style), tolower(colnames(mapping)))
    if(length(j) != 1) {
        if(verbose)
            message(paste("extendedMapSeqlevels: the 'style'", style, "is currently not supported for the 'species'", species, "in GenomeInfoDb. Check valid naming styles by running GenomeInfoDb::genomeStyles(species). If it's not present, consider adding your genome by following the information at http://www.bioconductor.org/packages/release/bioc/vignettes/GenomeInfoDb/inst/doc/Accept-organism-for-GenomeInfoDb.pdf"))
        return(seqnames)
    }
    
    ## Is the target style the same as the current style?
    if(tolower(style) == tolower(currentStyle)) {
        return(seqnames)
    }
    
    ## Find current naming map
    k <- grep(tolower(currentStyle), tolower(colnames(mapping)))
    
    ## Check currentStyle if it was supplied
    if(!guessedCurrent) {
        if(length(k) != 1) {
            if(verbose)
                message(paste("extendedMapSeqlevels: the 'currentStyle'", currentStyle, "is currently not supported for the 'species'", species, "in GenomeInfoDb. Check valid naming styles by running GenomeInfoDb::genomeStyles(species). If it's not present, consider adding your genome by following the information at http://www.bioconductor.org/packages/release/bioc/vignettes/GenomeInfoDb/inst/doc/Accept-organism-for-GenomeInfoDb.pdf"))
            return(seqnames)
        }
    }
    
    ## Make the map
    seqnames <- mapping[match(seqnames, mapping[, k]), j]
    if(verbose & guessed)
        message(paste('extendedMapSeqlevels: sequence names mapped from', currentStyle, 'to', style, 'for species', species))
    return(seqnames)
}

.selectBestGuess <- function(seqnames, organisms, styles, supported) {
    organisms <- tolower(gsub(" ", "_", organisms))
    seqnames <- tolower(seqnames)
    if(length(organisms) == 1 & length(styles) == 1) {
        res <- list('species' = organisms, 'style' = styles)
    } else {
        scores <- mapply(function(org, sty) {
            i <- which(names(supported) == org)
            if(length(i) != 1) return(0)
            k <- grep(sty, colnames(supported[[i]]))
            sum(seqnames %in% tolower(supported[[i]][, k[1]]))            
        }, organisms, styles)    
        if(max(scores) == 0) {
            res <- list('species' = NA, 'style' = NA)
        } else {
            res <- list('species' = organisms[which.max(scores)], 'style' = styles[which.max(scores)])
        }
    }
    return(res)
}
