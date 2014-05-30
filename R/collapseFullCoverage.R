#' Collapse full coverage information for efficient quantile computations
#'
#' For a given data set this function collapses the full coverage information 
#' for each sample from all the chromosomes. The resulting information per 
#' sample is the number of bases with coverage 0, 1, etc. It is similar to 
#' using table() on a regular vector. This information is then used by 
#' \link{sampleDepth} for calculating the sample depth adjustments. The data 
#' set can loaded to R using (see \link{fullCoverage} and optionally filtered 
#' using \link{filterData}. 
#' 
#' @param fullCov A list where each element is the result from 
#' \link[derfinder]{loadCoverage} used with \code{cutoff=NULL}. The elements of 
#' the list should be named according to the chromosome number. Can be 
#' generated using \link{fullCoverage}.
#' @param colsubset Which columns of \code{coverageInfo$coverage} to use.
#' @param save If \code{TRUE}, the result is saved as 'collapsedFull.Rdata'.
#' @param verbose If \code{TRUE} basic status updates will be printed along the 
#' way.
#'
#' @return 
#' A list with one element per sample. Then per sample, a list with two vector 
#' elements: \code{values} and \code{weights}. The first one is the coverage 
#' value and the second one is the number of bases with that value.
#'
#' @author Leonardo Collado-Torres
#' @seealso \link{fullCoverage}, \link{sampleDepth}
#' @export
#' @importMethodsFrom IRanges names '[' sort
#' @importFrom IRanges runValue runLength
#' @examples
#' ## Collapse the coverage information for the filtered data
#' collapsedFull <- collapseFullCoverage(list(genomeData$coverage), 
#'     verbose=TRUE)
#' collapsedFull
#' 
#' \dontrun{
#' ## You can also collapsed the raw data
#' collapsedFullRaw <- collapseFullCoverage(list(genomeDataRaw), verbose=TRUE)
#' }

collapseFullCoverage <- function(fullCov, colsubset = NULL, save = FALSE, 
    verbose = FALSE) {
    ## Remove un-used columns
    if (!is.null(colsubset)) {
        fullCov <- lapply(fullCov, function(x) {
            if ("coverage" %in% names(x)) {
                ## Extract coverage info if fullCov is the result from using
                ## lapply filterData() on the output of fullCoverage()
                res <- x$coverage[, colsubset]
            } else {
                res <- x[, colsubset]
            }
            return(res)
        })
    } else if ("coverage" %in% names(fullCov[[1]])) {
        ## Extract coverage info if fullCov is the result from using
        ## lapply filterData() on the output of fullCoverage()
        fullCov <- lapply(fullCov, function(x) x$coverage)
    }
    
    ## Sort to reduce run lengths
    if (verbose) 
        message(paste(Sys.time(), "sampleDepth: Sorting fullCov"))
    sortedFull <- lapply(fullCov, function(x) {
        sapply(x, sort)
    })
    
    ## Collapse chrs by sample
    if (verbose) 
        message(paste(Sys.time(),
            "sampleDepth: Collapsing chromosomes information by sample"))
    samples <- names(sortedFull[[1]])
    collapsedFull <- vector("list", length(samples))
    names(collapsedFull) <- samples
    for (sample in samples) {
        ## Extract the data
        chrdata <- lapply(sortedFull, function(x) x[[sample]])
        
        ## Collapse it
        values <- unlist(lapply(chrdata, runValue), use.names = FALSE)
        weights <- as.numeric(unlist(lapply(chrdata, runLength), 
            use.names = FALSE))
        
        ## Save it
        collapsedFull[[sample]] <- list(values = values, weights = weights)
    }
    
    ## Save for future use if you want to use another quantile
    if (save) {
        if (verbose) 
            message(paste(Sys.time(), "sampleDepth: Saving collapsedFull"))
        save(collapsedFull, file = "collapsedFull.Rdata")
    }
    
    return(collapsedFull)
} 
