#' Extract coverage information for a set of regions
#'
#' This function extracts the raw coverage information calculated by 
#' \link{fullCoverage} at each base for a set of regions found with 
#' \link{calculatePvalues}. It can further calculate the mean coverage per 
#' sample for each region.
#' 
#' @param fullCov A list where each element is the result from 
#' \link{loadCoverage} used with \code{returnCoverage = TRUE}. Can be generated 
#' using \link{fullCoverage}. Alternatively, specify \code{files} to extract
#' the coverage information from the regions of interest. This can be 
#' helpful if you do not wish to store \code{fullCov} for memory reasons.
#' @param regions The \code{$regions} output from \link{calculatePvalues}. It 
#' is important that the seqlengths information is provided.
#' @param totalMapped The total number of reads mapped for each sample. 
#' Providing this data adjusts the coverage to reads in \code{targetSize} 
#' library. By default, to reads per 80 million reads.
#' @param targetSize The target library size to adjust the coverage to. Used
#' only when \code{totalMapped} is specified.
#' @inheritParams fullCoverage
#' @param ... Arguments passed to other methods and/or advanced arguments.
#'
#' @return a list of data.frame where each data.frame has the coverage 
#' information (nrow = width of region, ncol = number of samples) for a given 
#' region. The names of the list correspond to the region indexes in 
#' \code{regions}
#' 
#' @author Andrew Jaffe, Leonardo Collado-Torres
#' @seealso \link{fullCoverage}, \link{calculatePvalues}
#' @export
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomeInfoDb seqlevels renameSeqlevels
#' mapSeqlevels seqlevelsInUse
#' @importMethodsFrom GenomicRanges names 'names<-' length '[' coverage sort 
#' width c '$'
#' @importMethodsFrom IRanges subset as.data.frame
#' @importFrom IRanges IRanges
#' @import S4Vectors
#' @importFrom BiocParallel bpmapply
#'
#' @details When \code{fullCov} is the output of \link{loadCoverage} with
#' \code{cutoff} non-NULL, \link{getRegionCoverage} assumes that the regions
#' come from the same data. Meaning that \link{filterData} was not used again.
#' This ensures that the regions are a subset of the data available in 
#' \code{fullCov}.
#'
#' If \code{fullCov} is \code{NULL} and \code{files} is specified, this function
#' will attempt to read the coverage from the files. Note that if you used
#' 'totalMapped' and 'targetSize' before, you will have to specify them again
#' to get the same results. 
#'
#' See also \link{advancedArg} with \code{fun='loadCoverage'} for other details.
#'
#' You should use at most one core per chromosome.
#' 
#'
#' @examples
#' ## Obtain fullCov object
#' fullCov <- list('21'=genomeDataRaw$coverage)
#' 
#' ## Assign chr lengths using hg19 information, use only first two regions
#' library('GenomicRanges')
#' data(hg19Ideogram, package = 'biovizBase', envir = environment())
#' regions <- genomeRegions$regions[1:2]
#' seqlengths(regions) <- seqlengths(hg19Ideogram)[names(seqlengths(regions))]
#'
#' ## Finally, get the region coverage
#' regionCov <- getRegionCoverage(fullCov=fullCov, regions=regions)

getRegionCoverage <- function(fullCov = NULL, regions, totalMapped = NULL, 
    targetSize = 80e6, files = NULL, ...) {
        
    ## Advanged arguments
# @param verbose If \code{TRUE} basic status updates will be printed along the 
# way.
    verbose <- .advanced_argument('verbose', TRUE, ...)


    names(regions) <- seq_len(length(regions))  # add names
    
    ## Use UCSC names for homo_sapiens by default
    regions <- renameSeqlevels(regions, 
        extendedMapSeqlevels(seqlevels(regions), ...))
    
    ## TODO check seqlengths are properly given in 'regions'
    
    ## Load data if 'fullCov' is not specified
    if(is.null(fullCov)) {
        fullCov <- .load_fullCov(files = files, regs = regions,
            fun = 'getRegionCoverage', ...)        
    }
    ## Fix naming style
    names(fullCov) <- extendedMapSeqlevels(names(fullCov), ...)
        
    # split by chromosome
    regions.chrs <- as.factor(seqnames(regions))
    ## Make sure the order matches the one from the fullCov names
    regions.chrs <- factor(regions.chrs, levels = names(fullCov))
    grl <- split(regions, regions.chrs)  
    
    ## Define args
    moreArgs <- list(totalMapped = totalMapped, verbose = verbose)
    
    ## Define cluster
    BPPARAM <- .define_cluster(...)
    
    counts <- bpmapply(function(chr, covInfo, g, totalMapped, verbose) {
        
        ## Parallel by chr, so no point in using mc.cores beyond the number of chrs
        if (verbose) 
            message(paste(Sys.time(), 'getRegionCoverage: processing', chr))
        
        ## Check whether fullCov has been filtered, then subset
        if(all(c('coverage', 'position') %in% names(covInfo))) {
            if(!is.null(g$indexStart) & !is.null(g$indexEnd)) {
                ## Subset if appropriate
                yy <- covInfo$coverage[IRanges(start = g$indexStart, end = g$indexEnd), ]
            } else if (is.null(covInfo$position)){
                yy <- covInfo$coverage[ranges(g), ]
            } else {
                stop("It seems that you have filtered the coverage but your 'regions' object is missing the 'indexStart' and 'indexEnd' information produced by findRegions().")
            }
            
        } else {
            yy <- covInfo[ranges(g), ]  # better subset
        }
            
        # depth-adjust, like for plotting
        if (!is.null(totalMapped) & targetSize != 0) {
            yy <- DataFrame(mapply(function(x, d) x / d, 
                yy, totalMapped / targetSize), check.names = FALSE)
        }
        
        ind <- rep(names(g), width(g))  # to split along
        ind <- factor(ind, levels = unique(ind))  # make factor in order
        # split(yy,ind) # 'CompressedSplitDataFrameList', faster but
        # less clear how to unlist below, so leave out
        res <- split(as.data.frame(yy), ind)
        
        if (verbose) 
            message(paste(Sys.time(), 'getRegionCoverage: done processing', 
                chr))
                
        ## Done
        return(res)
    }, names(fullCov), fullCov, grl, MoreArgs = moreArgs, BPPARAM = BPPARAM,
    SIMPLIFY = FALSE)
    covList <- do.call('c', counts)  # collect list elements into one large list
    
    # put in original order
    names(covList) <- sapply(strsplit(names(covList), '\\.'), 
        '[', 2)
    theData <- covList[order(as.numeric(names(covList)))]
    
    if (sum(sapply(theData, nrow)) != sum(width(regions))) {
        stop('The total width of the regions did not match with the dimensions of the extracted coverage data.')
    }
    
    return(theData)
}


## Helper function that runs fullCoverage
.load_fullCov <- function(files, regs, fun, ...) {
    stopifnot(!is.null(files))
    if(.advanced_argument('verbose', TRUE, ...))
        message(paste(Sys.time(), fun, ": attempting to load coverage data from 'files'."))

    ## If no protection was specified for calculating the coverage, then
    ## specify it. Details in loadCoverage()
    protectWhich <- .advanced_argument('protectWhich', NULL, ...)
    chrs <- seqlevelsInUse(regs)
    
        
    if(is.null(protectWhich)) {
        fullCov <- fullCoverage(files = files, chrs = chrs, protectWhich = 3e4,
            which = regs, ...)
    } else {
        fullCov <- fullCoverage(files = files, chrs = chrs, which = regs, ...)
    }
    return(fullCov)
}
