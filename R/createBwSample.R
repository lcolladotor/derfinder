#' Create a BigWig file with the coverage information for a given sample
#'
#' Given the output of \link{fullCoverage}, this function coerces the coverage
#' to a \link[GenomicRanges]{GRanges} object using \link{coerceGR} and then 
#' exports the coverage to a BigWig file using \link[rtracklayer]{export}.
#' 
#' @param sample The name or integer index of the sample of interest to coerce
#' to a \code{GRanges} object.
#' @param path The path where the BigWig file will be created.
#' @param fullCov A list where each element is the result from 
#' \link{loadCoverage} used with \code{returnCoverage = TRUE}. Can be generated 
#' using \link{fullCoverage}.
#' @param seqlengths A named vector with the sequence lengths of the 
#' chromosomes. This argument is passed to \link[GenomicRanges]{GRanges}.
#' @param keepGR If \code{TRUE}, the \link[GenomicRanges]{GRanges} object 
#' created by \link{coerceGR} is returned. Otherwise it is discarded.
#' @param mc.cores This argument is passed to \link[BiocParallel]{SnowParam} 
#' to define the number of \code{workers}. Use at most one core per chromosome.
#' @param mc.outfile This argument is passed to \link[BiocParallel]{SnowParam} 
#' to specify the \code{outfile} for any output from the workers.
#' @param verbose If \code{TRUE} basic status updates will be printed along the 
#' way.
#'
#' @return Creates a BigWig file with the coverage information (regions with
#' coverage greater than zero) for a given sample. If \code{keepGR} it returns
#' the output from \link{coerceGR}.
#'
#' @author Leonardo Collado-Torres
#' @seealso \link[GenomicRanges]{GRanges}, \link[rtracklayer]{export}, 
#' link{coerceGR}
#' @export
#' @aliases create_bw_sample
#'
#' @importMethodsFrom rtracklayer export
#' 
#' @examples
#' \dontrun{
#' }
#'

## Exports a single sample to a BigWig file
createBwSample <- function(sample, path = '.', fullCov, seqlengths,
    keepGR = TRUE, mc.cores = getOption("mc.cores", 1L), 
    mc.outfile = Sys.getenv('SGE_STDERR_PATH'), verbose = TRUE) {
    ## Coerce to GRanges
    gr.sample <- coerceGRanges(sample = sample, fullCov = fullCov,
        seqlengths = seqlengths, mc.cores = mc.cores, mc.outfile = mc.outfile, 
        verbose = verbose)
    
    ## Export bw file
    if(verbose) 
        message(paste(Sys.time(), "create_bw: exporting bw for sample", sample))
    
    ## Check that there is something to write
    if(length(gr.sample) > 0) {
        export(gr.sample, file.path(path, paste0(sample, ".bw")))
    } else {
        warning(paste0("There are no bases with coverage > 0 for sample ",
            sample, ". Thus no bw file will be created for this sample."))
    }    
    
    ## Keep the GRanges?
    if(keepGR) {
        return(invisible(gr.sample))
    } else {
        return(invisible(NULL))
    }
}

#' @export
create_bw_sample <- createBwSample
