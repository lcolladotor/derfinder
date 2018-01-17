#' Export coverage to BigWig files
#'
#' Using output from \link{fullCoverage}, export the coverage from all the
#' samples to BigWig files using \link{createBwSample}.
#' 
#' @param fullCov A list where each element is the result from 
#' \link{loadCoverage} used with \code{returnCoverage = TRUE}. Can be generated 
#' using \link{fullCoverage}.
#' @param path The path where the BigWig files will be created.
#' @param keepGR If \code{TRUE}, the \link[GenomicRanges:GRanges-class]{GRanges} objects 
#' created by \link{coerceGR} grouped into a \link[GenomicRanges:GRangesList-class]{GRangesList} 
#' are returned. Otherwise they are discarded.
#' @param ... Arguments passed to \link{createBwSample}.
#'
#' @return If \code{keepGR = TRUE}, then a \link[GenomicRanges:GRangesList-class]{GRangesList}
#' with the output for \link{coerceGR} for each of the samples.
#'
#' @details Use at most one core per chromosome.
#'
#' @author Leonardo Collado-Torres
#' @seealso \link[GenomicRanges:GRangesList-class]{GRangesList}, \link[rtracklayer]{export}, 
#' \link{createBwSample}, \link{coerceGR}
#' @export
#'
#' @importMethodsFrom GenomicRanges names
#' @importFrom GenomicRanges GRangesList
#' 
#' @examples
#' ## Create a small fullCov object with data only for chr21
#' fullCov <- list('chr21' = genomeDataRaw)
#'
#' ## Keep only 2 samples
#' fullCov$chr21$coverage <- fullCov$chr21$coverage[c(1, 31)]
#' 
#' ## Create the BigWig files for all samples in a test dir
#' dir.create('createBw-example')
#' bws <- createBw(fullCov, 'createBw-example')
#' 
#' ## Explore the output
#' bws
#'
#' ## First sample
#' bws[[1]]
#'
#' ## Note that if a sample has no bases with coverage > 0, the GRanges object
#' ## is empty and no BigWig file is created for that sample.
#' bws[[2]]
#'

## Exports fullCoverage() output to BigWig files
createBw <- function(fullCov, path = '.', keepGR = TRUE, ...) {
    
    ## Determine sample names
    samples <- names(fullCov[[1]])
    if('coverage' %in% samples) {
        samples <- names(fullCov[[1]]$coverage)
    }
    
    ## Coerce to GR
    gr.samples <- lapply(samples, createBwSample, path = path,
        fullCov = fullCov, keepGR = keepGR, ...)
    
    ## Done
    if(keepGR) {
        gr.samples <- GRangesList(gr.samples)
        names(gr.samples) <- samples
        return(invisible(gr.samples))
    } else {
        return(invisible(NULL))
    }
    
}
