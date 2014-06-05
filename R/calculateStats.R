#' Calculate F-statistics at base pair resolution from a loaded BAM files
#'
#' After defining the models of interest (see \link{makeModels}) and 
#' pre-processing the data (see \link{preprocessCoverage}), use 
#' \link{calculateStats} to calculate the F-statistics at base-pair resolution.
#' 
#' @param coveragePrep A list with \code{$coverageProcessed}, 
#' \code{$mclapplyIndex}, and \code{$position} normally generated using 
#' \link{preprocessCoverage}.
#' @param models A list with \code{$mod} and \code{$mod0} normally generated 
#' using \link{makeModels}.
#' @param mc.cores This argument is passed to \link[BiocParallel]{SnowParam} 
#' to define the number of \code{workers} used for running \link{fstats.apply}.
#' @param mc.outfile This argument is passed to \link[BiocParallel]{SnowParam} 
#' to specify the \code{outfile} for any output from the workers.
#' @param adjustF A single value to adjust that is added in the denominator of 
#' the F-stat calculation. Useful when the Residual Sum of Squares of the 
#' alternative model is very small.
#' @param lowMemDir The directory where the processed chunks are saved when 
#' using \link{preprocessCoverage} with a specified \code{lowMemDir}.
#' @param verbose If \code{TRUE} basic status updates will be printed along the 
#' way.
#'
#' @return A numeric Rle with the F-statistics per base pair that passed the 
#' cutoff.
#'
#' @author Leonardo Collado-Torres
#' @export
#' @seealso \link{makeModels}, \link{preprocessCoverage}
#' @importFrom BiocParallel SnowParam SerialParam bplapply
#' @importMethodsFrom IRanges ncol '[[' length unlist
#' @importFrom IRanges RleList
#' @examples
#' ## Collapse the coverage information
#' collapsedFull <- collapseFullCoverage(list(genomeData$coverage), 
#'     verbose=TRUE)
#' 
#' ## Calculate library size adjustments
#' sampleDepths <- sampleDepth(collapsedFull, probs=c(0.5), nonzero=TRUE, 
#'     verbose=TRUE)
#' 
#' ## Build the models
#' group <- genomeInfo$pop
#' adjustvars <- data.frame(genomeInfo$gender)
#' models <- makeModels(sampleDepths, testvars=group, adjustvars=adjustvars)
#'
#' ## Preprocess the data
#' prep <- preprocessCoverage(genomeData, cutoff=0, scalefac=32, chunksize=1e3, 
#'     colsubset=NULL)
#' 
#' ## Run the function
#' fstats <- calculateStats(prep, models, mc.cores=1, verbose=TRUE)
#' fstats
#'
#' \dontrun{
#' ## Compare vs pre-packaged F-statistics
#' summary(fstats - genomeFstats)
#' }

calculateStats <- function(coveragePrep, models, 
    mc.cores = getOption("mc.cores", 2L),
    mc.outfile = Sys.getenv('SGE_STDERR_PATH'), adjustF = 0, lowMemDir = NULL,
    verbose = TRUE) {
    
    stopifnot(length(intersect(names(coveragePrep), c("coverageProcessed", 
        "mclapplyIndex", "position"))) == 3)
    stopifnot(length(intersect(names(models), c("mod", "mod0"))) == 
        2)
    
    coverageProcessed <- coveragePrep$coverageProcessed
    if (is.null(lowMemDir) & is.null(coverageProcessed)) 
        stop("preprocessCoverage() was used with a non-null 'lowMemDir', so please specify 'lowMemDir'.")
    mclapplyIndex <- coveragePrep$mclapplyIndex
    rm(coveragePrep)
    
    if (is.null(lowMemDir)) {
        ## Check that the columns match
        numcol <- ncol(coverageProcessed)
        if (numcol != dim(models$mod)[1]) {
            stop("The alternative model 'models$mod' is not compatible with the number of samples in 'coveragePrep$coverageProcessed'. Check the dimensions of the alternative model.")
        }
        
        chunks <- length(coverageProcessed)
        if (chunks < mc.cores) {
            warning("The number of chunks in coveragePrep$coverageProcessed is smaller than the number of cores selected. For using all the cores specified consider splitting the data into more chunks.")
        }
    }
    
    ## Define cluster
    if(mc.cores > 1) {
        BPPARAM <- SnowParam(workers = mc.cores, outfile = mc.outfile)
    } else {
        BPPARAM <- SerialParam()
    }
    
    ## Fit a model to each row (chunk) of database:
    if (verbose) 
        message(paste(Sys.time(), "calculateStats: calculating the F-statistics"))
    fstats.output <- bplapply(mclapplyIndex, fstats.apply, 
        data = coverageProcessed, mod = models$mod, mod0 = models$mod0,
        adjustF = adjustF, lowMemDir = lowMemDir, BPPARAM = BPPARAM)
    ## Using mclapply is as fast as using lapply if mc.cores=1, so
    ## there is no damage in setting the default mc.cores=1.
    ## Specially since parallel is included in R 3.0.x More at
    ## http://stackoverflow.com/questions/16825072/deprecation-of-multicore-mclapply-in-r-3-0
    result <- unlist(RleList(fstats.output), use.names = FALSE)
    
    ## Done =)
    return(result)    
} 
