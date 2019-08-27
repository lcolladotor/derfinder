#' Create a BigWig file with the coverage information for a given sample
#'
#' Given the output of [fullCoverage], this function coerces the coverage
#' to a [GRanges][GenomicRanges::GRanges-class] object using [coerceGR] and then 
#' exports the coverage to a BigWig file using [export][rtracklayer::export].
#' 
#' @param sample The name or integer index of the sample of interest to coerce
#' to a `GRanges` object.
#' @param path The path where the BigWig file will be created.
#' @param fullCov A list where each element is the result from 
#' [loadCoverage] used with `returnCoverage = TRUE`. Can be generated 
#' using [fullCoverage].
#' @param keepGR If `TRUE`, the [GRanges][GenomicRanges::GRanges-class] object 
#' created by [coerceGR] is returned. Otherwise it is discarded.
#' @param ... Arguments passed to other methods and/or advanced arguments.
#' Advanced arguments:
#' \describe{
#' \item{verbose }{ If `TRUE` basic status updates will be printed along 
#' the way.}
#' }
#' Passed to [coerceGR].
#'
#' @return Creates a BigWig file with the coverage information (regions with
#' coverage greater than zero) for a given sample. If `keepGR` it returns
#' the output from [coerceGR].
#'
#' @author Leonardo Collado-Torres
#' @seealso [GRanges][GenomicRanges::GRanges-class], [export][rtracklayer::export], 
#' link{coerceGR}
#' @export
#'
#' @importMethodsFrom rtracklayer export
#' 
#' @examples
#' ## Create a small fullCov object with data only for chr21
#' fullCov <- list('chr21' = genomeDataRaw)
#'
#' ## Create a BigWig for the first sample in a test directory
#' dir.create('createBwSample-example')
#' bw <- createBwSample('ERR009101', 'createBwSample-example', 
#'     fullCov = fullCov, seqlengths = c('chr21' = 48129895))
#'
#' ## Explore the output
#' bw

## Exports a single sample to a BigWig file
createBwSample <- function(sample, path = '.', fullCov, keepGR = TRUE, ...) {

    ## Advanced arguments
# @param verbose If \code{TRUE} basic status updates will be printed along the 
# way.
    verbose <- .advanced_argument('verbose', TRUE, ...)

    
    ## Coerce to GRanges
    gr.sample <- coerceGR(sample = sample, fullCov = fullCov, ...)

    ## Export bw file
    if(verbose) 
        message(paste(Sys.time(), 'createBwSample: exporting bw for sample',
            sample))
    
    ## Check that there is something to write
    if(length(gr.sample) > 0) {
        if(.Platform$OS.type == 'windows') {
            warning('As of rtracklayer 1.25.16, BigWig is not supported on Windows. Thus exporting to BigWig will most likely fail!')
        }
        export(gr.sample, file.path(path, paste0(sample, '.bw')))
    } else {
        warning(paste0('There are no bases with coverage > 0 for sample ',
            sample, '. Thus no bw file will be created for this sample.'))
    }    
    
    ## Keep the GRanges?
    if(keepGR) {
        return(invisible(gr.sample))
    } else {
        return(invisible(NULL))
    }
}
