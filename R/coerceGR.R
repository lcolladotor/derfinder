#' Coerce the coverage to a GRanges object for a given sample
#' 
#' Given the output of \link{fullCoverage}, coerce the coverage to a 
#' \link[GenomicRanges]{GRanges} object.
#'
#' @param sample The name or integer index of the sample of interest to coerce
#' to a \code{GRanges} object.
#' @param fullCov A list where each element is the result from 
#' \link{loadCoverage} used with \code{returnCoverage = TRUE}. Can be generated 
#' using \link{fullCoverage}.
#' @param seqlengths A named vector with the sequence lengths of the 
#' chromosomes. This argument is passed to \link[GenomicRanges]{GRanges}.
#' @param mc.cores This argument is passed to \link[BiocParallel]{SnowParam} 
#' to define the number of \code{workers}. Use at most one core per chromosome.
#' @param mc.outfile This argument is passed to \link[BiocParallel]{SnowParam} 
#' to specify the \code{outfile} for any output from the workers.
#' @param verbose If \code{TRUE} basic status updates will be printed along the 
#' way.
#'
#' @return A \link[GenomicRanges]{GRanges} object with \code{score} metadata 
#' vector containing the coverage information for the specified sample. The 
#' ranges reported are only those for regions of the genome with coverage 
#' greater than zero.
#'
#' @author Leonardo Collado-Torres
#' @seealso \link[GenomicRanges]{GRanges}
#' @export
#' @aliases coerce_gr
#'
#' @importFrom BiocParallel SnowParam SerialParam bpmapply
#' @importFrom GenomicRanges GRangesList
#' 
#' @examples
#' ## Create a small fullCov object with data only for chr21
#' fullCov <- list('chr21' = genomeDataRaw)
#'
#' ## Coerce to a GRanges the first sample
#' gr <- createBwSample('ERR009101', fullCov = fullCov,
#'     seqlengths = c('chr21' = 48129895))
#'
#' ## Explore the output
#' gr
#'

## Coerces fullCoverage() output to GRanges for a given sample
coerceGR <- function(sample, fullCov, seqlengths = sapply(fullCov, nrow), 
    mc.cores = getOption("mc.cores", 1L), 
    mc.outfile = Sys.getenv('SGE_STDERR_PATH'), verbose = TRUE) {
    
    if(verbose) 
        message(paste(Sys.time(), "coerceGR: coercing sample", sample))
    
    ## Define cluster
    if(mc.cores > 1) {
        BPPARAM <- SnowParam(workers = mc.cores, outfile = mc.outfile)
    } else {
        BPPARAM <- SerialParam()
    }
            
    ## Coerce to a list of GRanges (1 element per chr)
    gr.sample <- bpmapply(function(chr, DF, sample, seqlengths) {
        ## Extract sample Rle info
        if("coverage" %in% names(DF)) {
            rle <- DF$coverage[[sample]]
        } else {
            rle <- DF[[sample]]
        }
        
        ## Rle values
        vals <- runValue(rle)
        idx <- which(vals > 0)
    
        if(length(idx) == 0) {
            ## Nothing found
            res <- GRanges(seqlengths = seqlengths)
            
        } else {
            ## Construct IRanges
            lens <- runLength(rle)
            IR <- IRanges(start = cumsum(c(1, lens))[idx], width = lens[idx])
            
            ## Finish
            res <- GRanges(seqnames = rep(chr, length(IR)), ranges = IR, strand = "*", seqlengths = seqlengths, score = vals[idx])
        }  
    
        return(res)
    }, names(fullCov), fullCov, MoreArgs = list(sample = sample, seqlengths = seqlengths), BPPARAM = BPPARAM)
    gr.sample <- unlist(GRangesList(gr.sample))
    
    ## Done
    return(gr.sample)
}

#' @export
coerce_gr <- coerceGR
