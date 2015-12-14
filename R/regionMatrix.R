#' Identify regions data by a coverage filter and get a count matrix
#'
#' Given a set of un-filtered coverage data (see \link{fullCoverage}), create
#' candidate regions by applying a cutoff on the coverage values,
#' and obtain a count matrix where the number of rows corresponds to the number
#' of candidate regions and the number of columns corresponds to the number of 
#' samples. The values are the mean coverage for a given sample for a given 
#' region.
#' 
#' @param fullCov A list where each element is the result from 
#' \link{loadCoverage} used with \code{returnCoverage = TRUE}. Can be generated 
#' using \link{fullCoverage}. If \code{runFilter = FALSE}, then 
#' \code{returnMean = TRUE} must have been used.
#' @inheritParams filterData
#' @param L The width of the reads used. Either a vector of length 1 or length
#' equal to the number of samples.
#' @param runFilter This controls whether to run \link{filterData} or not. If 
#' set to \code{FALSE} then \code{returnMean = TRUE} must have been used to 
#' create each element of \code{fullCov}.
#' @param returnBP If \code{TRUE}, returns \code{$bpCoverage} explained below.
#' @param ... Arguments passed to other methods and/or advanced arguments.
#'
#' @return A list with one entry per chromosome. Then per chromosome, a list 
#' with three components.
#' \describe{
#' \item{regions }{ A set of regions based on the coverage filter cutoff as
#' returned by \link{findRegions}.}
#' \item{bpCoverage }{ A list with one element per region. Each element is a matrix with numbers of rows equal to the number of base pairs in the region and number of columns equal to the number of samples. It contains the base-level coverage information for the regions. Only returned when \code{returnBP = TRUE}.}
#' \item{coverageMatrix }{  A matrix with the mean coverage by sample for each
#' candidate region.}
#' }
#'
#' @details This function uses several other \link{derfinder-package} 
#' functions. Inspect the code if interested.
#'
#' You should use at most one core per chromosome.
#'
#' @author Leonardo Collado-Torres
#' @export
#'
#' @importMethodsFrom IRanges nrow '$<-'
#' @importFrom BiocParallel bpmapply
#' @importFrom GenomeInfoDb 'seqlengths<-'
#'
#' @examples
#' ## Create some toy data
#' library('IRanges')
#' x <- Rle(round(runif(1e4, max=10)))
#' y <- Rle(round(runif(1e4, max=10)))
#' z <- Rle(round(runif(1e4, max=10)))
#' fullCov <- list('chr21' = DataFrame(x, y, z))
#'
#' ## Calculate a proxy of library size
#' libSize <- sapply(fullCov$chr21, sum)
#'
#' ## Run region matrix normalizing the coverage
#' regionMat <- regionMatrix(fullCov = fullCov, maxRegionGap = 10L, 
#'     maxClusterGap = 300L, L = 36, totalMapped = libSize, targetSize = 4e4)
#'
#' \dontrun{
#' ## You can alternatively use filterData() on fullCov to reduce the required
#' ## memory before using regionMatrix(). This can be useful when mc.cores > 1
#' filteredCov <- lapply(fullCov, filterData, returnMean=TRUE, filter='mean', 
#'     cutoff=5, totalMapped = libSize, targetSize = 4e4)
#' regionMat2 <- regionMatrix(filteredCov, maxRegionGap = 10L, 
#'     maxClusterGap = 300L, L = 36, runFilter=FALSE)
#' identical(regionMat2, regionMat)
#' }

regionMatrix <- function(fullCov, cutoff = 5, filter = 'mean', L,
    runFilter = TRUE, returnBP = TRUE, ...) {
        
    ## Have to filter by something
    stopifnot(!is.null(cutoff))
    
    ## fullCov has to be named
    stopifnot(!is.null(names(fullCov)))
    
    ## If filtering has been run, check that the information is there
    if(!runFilter) {
        if(!all(sapply(fullCov, function(x) { all(c('coverage', 'position', 'meanCoverage') %in% names(x)) })))
            stop("When 'runFilter' = FALSE, all the elements of 'fullCov' are expected to be the output from filterData(..., returnMean=TRUE) with some non-NULL 'cutoff'.")
    }
        
    ## Define args
    moreArgs <- list(cutoff = cutoff, filter = filter, L = L,
        runFilter = runFilter, returnBP = returnBP, ...)
    
    ## Define cluster
    BPPARAM <- .define_cluster(...)
        
    ## Get regions per chr
    bpmapply(.regionMatrixByChr, fullCov, names(fullCov), MoreArgs = moreArgs, 
        SIMPLIFY = FALSE, BPPARAM = BPPARAM)
}

.regionMatrixByChr <- function(covInfo, chr, cutoff, filter, 
    maxClusterGap = 300L, L, runFilter, returnBP, ...) {

    
    verbose <- .advanced_argument('verbose', TRUE, ...)
    chrsStyle <- .advanced_argument('chrsStyle', 'UCSC', ...)


    if (verbose) 
        message(paste(Sys.time(), 'regionMatrix: processing', chr))
        
    ## Filter by 'one' or 'mean' and get mean coverage
    if(runFilter) {
        if(all(c('coverage', 'position') %in% names(covInfo))) {
            filt <- filterData(data = covInfo$coverage, cutoff = cutoff,
                filter = filter, returnMean = TRUE, returnCoverage = TRUE,
                index = covInfo$position, ...)
        } else {
            filt <- filterData(data = covInfo, cutoff = cutoff, filter = filter,
                returnMean = TRUE, returnCoverage = TRUE, ...)
        }
        
        ## Identify regions
        regs <- findRegions(position = filt$position, 
            fstats = filt$meanCoverage, chr = chr, cutoff = 0, ...)
            
        ## Prepare for getRegionCoverage
        fullCovTmp <- list(filt)
        seqlengths <- length(filt$position)
    } else {        
        ## Identify regions
        regs <- findRegions(position = covInfo$position, 
            fstats = covInfo$meanCoverage, chr = chr, cutoff = 0, ...)
        
        ## Prepare for getRegionCoverage
        fullCovTmp <- list(covInfo)
        seqlengths <- length(covInfo$position)
    }
    
    ## If there are no regions, return NULL
    if(is.null(regs)) {
        if(returnBP) {
            return(list(regions = NULL, coverageMatrix = NULL, bpCoverage = NULL))
        } else {
            return(list(regions = NULL, coverageMatrix = NULL))
        }
    }
    
    ## Format appropriately
    names(regs) <- seq_len(length(regs))
    
    ## Set the length
    names(seqlengths) <- chr
    seqlengths(regs) <- seqlengths

    ## Prepare for getRegionCoverage
    names(fullCovTmp) <- chr
        
    ## Get region coverage    
    regionCov <- getRegionCoverage(fullCov = fullCovTmp, regions = regs, 
        totalMapped = NULL, targetSize = 0, verbose = verbose,
        chrsStyle = chrsStyle, mc.cores = 1L)
        
    if(verbose) message(paste(Sys.time(), "regionMatrix: calculating coverageMatrix"))
    covMat <- lapply(regionCov, colSums)
    covMat <- do.call(rbind, covMat)
    
    if(verbose) message(paste(Sys.time(), "regionMatrix: adjusting coverageMatrix for 'L'"))
    if(length(L) == 1) {
        covMat <- covMat / L
    } else if (length(L) == ncol(covMat)) {
        for(i in length(L)) covMat[, i] <- covMat[, i] / L[i]
    } else {
        warning("Invalid 'L' value so it won't be used. It has to either be a integer/numeric vector of length 1 or length equal to the number of samples.")
    }

    ## Finish
    if(returnBP) {
        res <- list(regions = regs, coverageMatrix = covMat, 
            bpCoverage = regionCov)
    } else {
        res <- list(regions = regs, coverageMatrix = covMat)
    }
    
    return(res)
    
}
