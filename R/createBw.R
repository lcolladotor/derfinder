#' Export coverage to BigWig files
#'
#' Using output from \link{fullCoverage}, export the coverage from all the
#' samples to BigWig files using \link{createBwSample}.
#' 
#' @param fullCov A list where each element is the result from 
#' \link{loadCoverage} used with \code{returnCoverage = TRUE}. Can be generated 
#' using \link{fullCoverage}.
#' @param path The path where the BigWig files will be created.
#' @param keepGR If \code{TRUE}, the \link[GenomicRanges]{GRanges} objects 
#' created by \link{coerceGR} grouped into a \link[GenomicRanges]{GRangesList} 
#' are returned. Otherwise they are discarded.
#' @param mc.cores This argument is passed to \link[BiocParallel]{SnowParam} 
#' to define the number of \code{workers}. Use at most one core per chromosome.
#' @param mc.outfile This argument is passed to \link[BiocParallel]{SnowParam} 
#' to specify the \code{outfile} for any output from the workers.
#' @param verbose If \code{TRUE} basic status updates will be printed along the 
#' way.
#'
#' @return If \code{keepGR = TRUE}, then a \link[GenomicRanges]{GRangesList}
#' with the output for \link{coerceGR} for each of the samples.
#'
#' @author Leonardo Collado-Torres
#' @seealso \link[GenomicRanges]{GRangesList}, \link[rtracklayer]{export}, 
#' \link{createBwSample}, \link{coerceGR}
#' @export
#' @aliases create_bw
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
createBw <- function(fullCov, path = '.', keepGR = TRUE, 
    mc.cores = getOption('mc.cores', 1L), 
    mc.outfile = Sys.getenv('SGE_STDERR_PATH'), verbose = TRUE) {
    
    ## Determine seqlengths and sample names
    samples <- names(fullCov[[1]])
    if('coverage' %in% samples) {
        seqlengths <- sapply(fullCov, function(x) { nrow(x$coverage )})
        samples <- names(fullCov[[1]]$coverage)
    } else {
        seqlengths <- sapply(fullCov, nrow)    
    }
    
    ## Coerce to GR
    gr.samples <- lapply(samples, createBwSample, path = path,
        fullCov = fullCov, seqlengths = seqlengths, keepGR = keepGR,
        mc.cores = mc.cores, mc.outfile = mc.outfile, verbose = verbose)
    
    ## Done
    if(keepGR) {
        gr.samples <- GRangesList(gr.samples)
        names(gr.samples) <- samples
        return(invisible(gr.samples))
    } else {
        return(invisible(NULL))
    }
    
}

#' @export
create_bw <- createBw
