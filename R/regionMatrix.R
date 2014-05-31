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
#' \link{loadCoverage} used with \code{cutoff=NULL}. Can be generated using 
#' \link{fullCoverage}.
#' @param cutoff Per base pair, at least one sample has to have coverage 
#' strictly greater than \code{cutoff} to be included in the result. This
#' argument is passed to \link{filterData}.
#' @param filter Has to be either \code{"one"} (default) or \code{"mean"}. In 
#' the first case, at least one sample has to have coverage above \code{cutoff}.
#' In the second case, the mean coverage has to be greater than \code{cutoff}.
#' This argument is passed to \link{filterData}.
#' @param maxRegionGap This determines the maximum number of gaps between two 
#' genomic positions to be considered part of the same candidate Differentially 
#' Expressed Region (candidate DER). This argument is passed to 
#' \link{findRegions}.
#' @param maxClusterGap This determines the maximum gap between candidate DERs. 
#' It should be greater than \code{maxRegionGap}. This argument is passed to
#' \link{findRegions}.
#' @param L The width of the reads used. This argument is passed to
#' \link{coverageToExon}.
#' @param totalMapped The total number of reads mapped for each sample. 
#' Providing this data adjusts the coverage to reads in \code{targetSize} 
#' library. By default, to reads per 80 million reads.
#' @param targetSize The target library size to adjust the coverage to. Used
#' only when \code{totalMapped} is specified.
#' @param mc.cores The number of cores to use. Maximum, 1 core per chr.
#' @param verbose If \code{TRUE} basic status updates will be printed along the 
#' way.
#'
#' @return A list with one entry per chromosome. Then per chromosome, a list 
#' with two components.
#' \describe{
#' \item{regions }{ A set of regions based on the coverage filter cutoff as
#' returned by \link{findRegions}.}
#' \item{coverageMatrix }{  A matrix with the mean coverage by sample for each
#' candidate region.}
#' }
#'
#' @details This function uses several other \link{derfinder-package} functions. Inspect
#' the code if interested.
#'
#' @author Leonardo Collado-Torres
#' @export
#' @importFrom GenomeInfoDb 'seqlengths<-'
#' @importMethodsFrom IRanges nrow '$<-'
#' @importFrom parallel mcmapply
#' @examples
#' library('IRanges')
#' x <- Rle(round(runif(1e4, max=10)))
#' y <- Rle(round(runif(1e4, max=10)))
#' z <- Rle(round(runif(1e4, max=10)))
#' fullCov <- list("chr21" = DataFrame(x, y, z))
#' regionMat <- regionMatrix(fullCov = fullCov, maxRegionGap = 10L, 
#'     maxClusterGap = 300L, L = 36)

regionMatrix <- function(fullCov, cutoff = 5, filter = "mean", 
    maxRegionGap = 0L, maxClusterGap = 300L, L, totalMapped = NULL,
    targetSize = 80e6, mc.cores = 1, verbose = TRUE) {
        
    ## library size adjustments
    if(!is.null(totalMapped)) {
        mappedPerXM <- totalMapped / targetSize
    } else {
        mappedPerXM <- NULL
    }
    
    ## Define args
    moreArgs <- list(cutoff = cutoff, filter = filter,
        maxRegionGap = maxRegionGap, maxClusterGap = maxClusterGap, L = L,
        mappedPerXM = mappedPerXM, verbose = verbose)
    
    ## Get regions per chr
    mcmapply(.regionMatrixByChr, fullCov, names(fullCov), MoreArgs = moreArgs, 
        SIMPLIFY = FALSE, mc.cores = mc.cores)
}

.regionMatrixByChr <- function(covInfo, chr, cutoff, filter, maxRegionGap = 0L, 
    maxClusterGap = 300L, L, mappedPerXM, verbose) {
    
    if (verbose) 
        message(paste(Sys.time(), "regionMatrix: processing", chr))
        
    ## Normalize to a given size
    if(!is.null(mappedPerXM)) {
        if (verbose) 
            message(paste(Sys.time(), "regionMatrix: normalizing coverage"))
        covInfo <- DataFrame(mapply(function(x, d) x / d, covInfo, mappedPerXM))
    }    
        
    ## Filter by 'one' or 'mean' and get mean coverage
    filt <- filterData(data = covInfo, cutoff = cutoff, filter = filter,
        returnMean = TRUE, returnCoverage = FALSE, verbose = verbose)
    
    ## Identify regions
    regs <- findRegions(position = filt$position, fstats = filt$meanCoverage,
        chr = chr, cutoff = 0, maxRegionGap = maxRegionGap,
        maxClusterGap = maxClusterGap, verbose = verbose)
    
    ## Format appropriately
    seqlengths(regs) <- nrow(covInfo)
    names(regs) <- seq_len(length(regs))
    
    ## Get coverage
    fullCovTmp <- list(covInfo)
    names(fullCovTmp) <- chr
    
    regionCov <- getRegionCoverage(fullCov = fullCovTmp, regions = regs, 
        totalMapped = NULL, mc.cores = 1, verbose = verbose)
        
    covMat <- lapply(regionCov, colSums)
    covMat <- do.call(rbind, covMat) / L

    ## Finish
    res <- list(regions = regs, coverageMatrix = covMat)
    return(res)
    
}