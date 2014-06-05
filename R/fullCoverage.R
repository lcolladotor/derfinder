#' Load the unfiltered coverage information from a group of BAM files and a 
#' list of chromosomes
#'
#' For a group of samples this function reads the coverage information for 
#' several chromosomes directly from the BAM files. Per chromosome, it merges 
#' the unfiltered coverage by sample into a DataFrame. The end result is a list 
#' with one such DataFrame objects per chromosome.
#' 
#' @param dirs A character vector with the full path to the sample BAM files
#' (or bigWig files). 
#' The names are used for the column names of the DataFrame. Check 
#' \link{makeBamList} for constructing \code{dirs}. \code{dirs} can also be a 
#' \code{BamFileList} object created with \link[Rsamtools]{BamFileList} or a
#' \code{BigWigFileList} object created with \link[rtracklayer]{BigWigFileList}.
#' @param chrs The chromosome of the files to read. The format has to match the
#' one used in the input files.
#' @param bai The full path to the BAM index files. If \code{NULL} it is 
#' assumed that the BAM index files are in the same location as the BAM files 
#' and that they have the .bai extension. Ignored if \code{dirs} is a 
#' \code{BamFileList} object.
#' @param chrlens The chromosome lengths in base pairs. If it's \code{NULL}, 
#' the chromosome length is extracted from the BAM files. Otherwise, it should 
#' have the same length as \code{chrs}.
#' @param outputs This argument is passed to the \code{output} argument of 
#' \link{loadCoverage}. If \code{NULL} or \code{'auto'} it is then recycled.
#' @param mc.cores This argument is passed to \link[BiocParallel]{SnowParam} 
#' to define the number of \code{workers}. You should use at most one core per 
#' chromosome.
#' @param mc.outfile This argument is passed to \link[BiocParallel]{SnowParam} 
#' to specify the \code{outfile} for any output from the workers.
#' @param inputType Has to be either \code{bam} or \code{bigWig}. It specifies
#' the format of the raw data files.
#' @param isMinusStrand Use \code{TRUE} for negative strand alignments only, 
#' \code{FALSE} for positive strands and \code{NA} for both. This argument is 
#' passed to \link[Rsamtools]{scanBamFlag} when \code{inputType='bam'}.
#' @param chrsStyle The naming style of the chromosomes. By default, UCSC. See 
#' \link[GenomeInfoDb]{seqlevelsStyle}.
#' @param verbose If \code{TRUE} basic status updates will be printed along the 
#' way.
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
#' @importFrom BiocParallel SnowParam SerialParam bplapply
#' @importFrom GenomeInfoDb mapSeqlevels
#'
#' @examples
#' datadir <- system.file('extdata', 'genomeData', package='derfinder')
#' dirs <- makeBamList(datadir=datadir, samplepatt='*accepted_hits.bam$', 
#'     bamterm=NULL)
#' ## Shorten the column names
#' names(dirs) <- gsub('_accepted_hits.bam', '', names(dirs))
#' 
#' ## Read and filter the data, only for 1 file
#' fullCov <- fullCoverage(dirs=dirs[1], chrs=c('21', '22'))
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

fullCoverage <- function(dirs, chrs, bai = NULL, chrlens = NULL, 
    outputs = NULL, mc.cores = getOption("mc.cores", 1L), 
    mc.outfile = Sys.getenv('SGE_STDERR_PATH'), inputType = "bam", 
    isMinusStrand = NA, chrsStyle = "UCSC", verbose = TRUE) {
        
    stopifnot(length(chrlens) == length(chrs) | is.null(chrlens))
    if (!is.null(outputs)) {
        stopifnot(length(outputs) == length(chrs) | outputs == "auto")
        if (outputs == "auto") {
            outputs <- rep("auto", length(chrs))
        }
    }
    
    ## Define cluster
    if(mc.cores > 1) {
        BPPARAM <- SnowParam(workers = mc.cores, outfile = mc.outfile)
    } else {
        BPPARAM <- SerialParam()
    }
    
    ## Subsetting function that runs loadCoverage
    loadChr <- function(idx, dirs, chrs, bai, chrlens, outputs, inputType,
        isMinusStrand, verbose) {
        
        if (verbose) 
            message(paste(Sys.time(), "fullCoverage: processing chromosome", 
                chrs[idx]))
        loadCoverage(dirs = dirs, chr = chrs[idx], cutoff = NULL, 
            bai = bai, chrlen = chrlens[idx], output = outputs[idx], 
            inputType = inputType, isMinusStrand = isMinusStrand,
            verbose = verbose)$coverage
    }
    result <- bplapply(seq_len(length(chrs)), loadChr, dirs = dirs, 
        chrs = chrs, bai = bai, chrlens = chrlens, outputs = outputs, 
        verbose = verbose, inputType = inputType, 
        isMinusStrand = isMinusStrand, BPPARAM = BPPARAM)
    names(result) <- mapSeqlevels(chrs, chrsStyle)
    
    ## Done
    return(result)
} 
