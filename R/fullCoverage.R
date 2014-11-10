#' Load the unfiltered coverage information from a group of BAM files and a 
#' list of chromosomes
#'
#' For a group of samples this function reads the coverage information for 
#' several chromosomes directly from the BAM files. Per chromosome, it merges 
#' the unfiltered coverage by sample into a DataFrame. The end result is a list 
#' with one such DataFrame objects per chromosome.
#' 
#' @param files A character vector with the full path to the sample BAM files
#' (or BigWig files). 
#' The names are used for the column names of the DataFrame. Check 
#' \link{rawFiles} for constructing \code{files}. \code{files} can also be a 
#' \code{BamFileList} object created with \link[Rsamtools]{BamFileList} or a
#' \code{BigWigFileList} object created with \link[rtracklayer]{BigWigFileList}.
#' @param chrs The chromosome of the files to read. The format has to match the
#' one used in the input files.
#' @param bai The full path to the BAM index files. If \code{NULL} it is 
#' assumed that the BAM index files are in the same location as the BAM files 
#' and that they have the .bai extension. Ignored if \code{files} is a 
#' \code{BamFileList} object.
#' @param chrlens The chromosome lengths in base pairs. If it's \code{NULL}, 
#' the chromosome length is extracted from the BAM files. Otherwise, it should 
#' have the same length as \code{chrs}.
#' @param outputs This argument is passed to the \code{output} argument of 
#' \link{loadCoverage}. If \code{NULL} or \code{'auto'} it is then recycled.
#' @inheritParams loadCoverage
#' @param ... Arguments passed to other methods and/or advanced arguments.
#'
#' @return A list with one element per chromosome.
#' \describe{ Each element is a DataFrame with the coverage information 
#' produced by \link{loadCoverage}.
#' }
#'
#' @seealso \link{loadCoverage}, \link{filterData}
#'
#' @author Leonardo Collado-Torres
#' @export
#' @aliases full_coverage
#' @importFrom BiocParallel bplapply
#' @importFrom GenomeInfoDb mapSeqlevels
#'
#' @examples
#' datadir <- system.file('extdata', 'genomeData', package='derfinder')
#' files <- rawFiles(datadir=datadir, samplepatt='*accepted_hits.bam$', 
#'     fileterm=NULL)
#' ## Shorten the column names
#' names(files) <- gsub('_accepted_hits.bam', '', names(files))
#' 
#' ## Read and filter the data, only for 1 file
#' fullCov <- fullCoverage(files=files[1], chrs=c('21', '22'))
#' fullCov
#'
#' \dontrun{
#' ## You can then use filterData() to filter the data if you want to. 
#' ## Use bplapply() if you want to do so with multiple cores as shown below.
#' library('BiocParallel')
#' p <- SnowParam(2L, outfile = Sys.getenv('SGE_STDERR_PATH'))
#' bplapply(fullCov, function(x) {
#'     library('derfinder'); filterData(x, cutoff=0) }, BPPARAM = p)
#' }

fullCoverage <- function(files, chrs, bai = NULL, chrlens = NULL, 
    outputs = NULL, cutoff = NULL, ...) {
        
    stopifnot(length(chrlens) == length(chrs) | is.null(chrlens))
    stopifnot(is.character(chrs))
    if (!is.null(outputs)) {
        stopifnot(length(outputs) == length(chrs) | outputs == 'auto')
        if (outputs == 'auto') {
            outputs <- rep('auto', length(chrs))
        }
    }
    
    ## Advanged argumentsa
#' @param verbose If \code{TRUE} basic status updates will be printed along the 
#' way.
    verbose <- .advanced_argument('verbose', TRUE, ...)


#' @param mc.cores.load Controls the number of cores to be used per chr for 
#' loading the data in chunks. Only used when \code{tilewidth} is specified.
    mc.cores.load <- .advanced_argument('mc.cores.load',
        .advanced_argument('mc.cores', getOption('mc.cores', 1L), ...), ...)

    
    ## Define cluster
    BPPARAM <- .define_cluster(...)
    
    ## Subsetting function that runs loadCoverage
    loadChr <- function(idx, files, chrs, bai, chrlens, outputs, cutoff,
        mc.cores.load, ...) {
        
        if (verbose) 
            message(paste(Sys.time(), 'fullCoverage: processing chromosome', 
                chrs[idx]))
        if (is.null(cutoff)) {
            res <- loadCoverage(files = files, chr = chrs[idx], cutoff = NULL, 
                bai = bai, chrlen = chrlens[idx], output = outputs[idx], 
                mc.cores = mc.cores.load, ...)$coverage
        } else {
            res <- loadCoverage(files = files, chr = chrs[idx],
                cutoff = cutoff, bai = bai, chrlen = chrlens[idx],
                output = outputs[idx], mc.cores = mc.cores.load, ...)
        }
        return(res)        
    }
    result <- bplapply(seq_len(length(chrs)), loadChr, files = files, 
        chrs = chrs, bai = bai, chrlens = chrlens, outputs = outputs, 
        cutoff = cutoff, mc.cores.load = mc.cores.load,
        ..., BPPARAM = BPPARAM)
    names(result) <- extendedMapSeqlevels(chrs, ...)
    
    ## Done
    return(result)
} 

#' @export
full_coverage <- fullCoverage
