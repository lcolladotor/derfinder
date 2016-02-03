#' Transform and split the data
#'
#' This function takes the coverage data from \link{loadCoverage}, scales the 
#' data, does the log2 transformation, and splits it into appropriate chunks 
#' for using \link{calculateStats}.
#' 
#' @param coverageInfo A list containing a DataFrame --\code{$coverage}-- with 
#' the coverage data and a logical Rle --\code{$position}-- with the positions 
#' that passed the cutoff. This object is generated using \link{loadCoverage}.
#' @param groupInfo A factor specifying the group membership of each sample. If 
#' \code{NULL} no group mean coverages are calculated. If the factor has more 
#' than one level, the first one will be used to calculate the log2 fold change 
#' in \link{calculatePvalues}.
#' @inheritParams filterData
#' @param colsubset Optional vector of column indices of 
#' \code{coverageInfo$coverage} that denote samples you wish to include in 
#' analysis. 
#' @param lowMemDir If specified, each chunk is saved into a separate Rdata 
#' file under \code{lowMemDir} and later loaded in 
#' \link[derfinderHelper]{fstats.apply} when 
#' running \link{calculateStats} and \link{calculatePvalues}. Using this option 
#' helps reduce the memory load as each fork in \link[BiocParallel]{bplapply} 
#' loads only the data needed for the chunk processing. The downside is a bit 
#' longer computation time due to input/output.
#' @param ... Arguments passed to other methods and/or advanced arguments.
#'
#' @details If \code{chunksize} is \code{NULL}, then \code{mc.cores} is used to 
#' determine the \code{chunksize}. This is useful if you want to split the data 
#' so each core gets the same amount of data (up to rounding).
#'
#' Computing the indexes and using those for \link[parallel]{mclapply} reduces 
#' memory copying as described by Ryan Thompson and illustrated in approach #4 
#' at \url{http://lcolladotor.github.io/2013/11/14/Reducing-memory-overhead-when-using-mclapply}
#'
#' If \code{lowMemDir} is specified then \code{$coverageProcessed} is NULL and 
#' \code{$mclapplyIndex} is a vector with the chunk identifiers.
#'
#' @return A list with five components.
#' \describe{
#' \item{coverageProcessed }{ contains the processed coverage information in a 
#' DataFrame object. Each column represents a sample and the coverage 
#' information is scaled and log2 transformed. Note that if \code{colsubset} is 
#' not \code{NULL} the number of columns will be less than those in 
#' \code{coverageInfo$coverage}. The total number of rows depends on the number 
#' of base pairs that passed the \code{cutoff} and the information stored is 
#' the coverage at that given base. Further note that \link{filterData} is 
#' re-applied if \code{colsubset} is not \code{NULL} and could thus lead to 
#' fewer rows compared to \code{coverageInfo$coverage}. }
#' \item{mclapplyIndex }{ is a list of logical Rle objects. They contain the 
#' partioning information according to \code{chunksize}.}
#' \item{position }{  is a logical Rle with the positions of the chromosome 
#' that passed the cutoff.}
#' \item{meanCoverage }{ is a numeric Rle with the mean coverage at each 
#' filtered base.}
#' \item{groupMeans }{ is a list of Rle objects containing the mean coverage at 
#' each filtered base calculated by group. This list has length 0 if 
#' \code{groupInfo=NULL}.}
#' }
#'
#' @author Leonardo Collado-Torres
#' @seealso \link{filterData}, \link{loadCoverage}, \link{calculateStats}
#' @export
#'
#' @importMethodsFrom IRanges ncol nrow '[' '[[' '[[<-' c
#' @import S4Vectors
#'
#' @examples
#' ## Split the data and transform appropriately before using calculateStats()
#' dataReady <- preprocessCoverage(genomeData, cutoff = 0, scalefac = 32, 
#'     chunksize = 1e3, colsubset = NULL, verbose = TRUE)
#' names(dataReady)
#' dataReady

preprocessCoverage <- function(coverageInfo, groupInfo = NULL, cutoff = 5, 
    colsubset = NULL, lowMemDir = NULL, ...) {
    ## Check that the input is from loadCoverage()
    stopifnot(length(intersect(names(coverageInfo), c('coverage', 
        'position'))) == 2)
    stopifnot(is.factor(groupInfo) | is.null(groupInfo))
    
    ## Advanged arguments
# @param verbose If \code{TRUE} basic status updates will be printed along the 
# way.
    verbose <- .advanced_argument('verbose', FALSE, ...)


# @param toMatrix Determines whether the data in the chunk should already be 
# saved as a Matrix object, which can be useful to reduce the computation time of the F-statistics. Only used when \code{lowMemDir} is not NULL.
    toMatrix <- .advanced_argument('toMatrix', !is.null(lowMemDir), ...) 


# @param mc.cores Number of cores you will use for calculating the statistics.
    mc.cores <- .advanced_argument('mc.cores', getOption('mc.cores', 1L), ...)


# @param scalefac A log transformation is used on the count tables, so zero 
# counts present a problem.  What number should we add to the entire matrix?
    scalefac <- .advanced_argument('scalefac', 32, ...)


# @param chunksize How many rows of \code{coverageInfo$coverage} should be 
# processed at a time?
    chunksize <- .advanced_argument('chunksize', 5e+06, ...)


    ## Use pre-calculated mean coverage if available
    means <- coverageInfo$meanCoverage
    
    coverage <- coverageInfo$coverage
    if (is.null(colsubset)) {
        stopifnot(length(groupInfo) == length(coverage) | is.null(groupInfo))
    } else {
        stopifnot(length(groupInfo) == length(colsubset) | is.null(groupInfo))
    }
    
    position <- coverageInfo$position
    
    ## Subset the DataFrame to use only the columns of interest
    if (!is.null(colsubset)) {
        ## Re-filter
        if (verbose) 
            message(paste(Sys.time(), 'preprocessCoverage: filtering the data'))
        coverageInfo <- filterData(data = coverageInfo$coverage[, 
            colsubset], cutoff = cutoff, index = coverageInfo$position, 
            returnMean=TRUE, ...)
        coverage <- coverageInfo$coverage
        position <- coverageInfo$position
        means <- coverageInfo$meanCoverage
    }
    rm(coverageInfo)
    
    ## Get the positions and shorter variables
    numrow <- nrow(coverage)
    
    ## Automatic chunksize depending on the number of cores
    if (is.null(chunksize)) {
        chunksize <- ceiling(numrow/mc.cores)
        if (verbose) 
            message(paste(Sys.time(), 'preprocessCoverage: using chunksize', 
                chunksize))
    }
    
    ## Determine number of loops
    lastloop <- trunc(numrow / chunksize)
    
    ## Fix the lastloop in case that the N is a factor of
    ## chunksize
    if (numrow %% chunksize == 0 & lastloop > 0) {
        lastloop <- lastloop - 1
    }
    
    ## Find the overall mean coverage
    if(is.null(means)) {
        means <- Reduce('+', coverage)/length(coverage)
    }
    
    
    ## Find the by group mean coverage
    if (is.null(groupInfo)) {
        groupMeans <- list()
    } else {
        groupMeans <- vector('list', length(levels(groupInfo)))
        names(groupMeans) <- levels(groupInfo)
        for (group in levels(groupInfo)) {
            groupMeans[[group]] <- Reduce('+', coverage[groupInfo == 
                group])/sum(groupInfo == group)
        }
    }
    
    ## Log2 transform and scale
    numcol <- ncol(coverage)
    for (i in seq_len(numcol)) {
        coverage[[i]] <- log2(coverage[[i]] + scalefac)
    }
    
    ## Split the data into appropriate chunks
    if (verbose) 
        message(paste(Sys.time(), 'preprocessCoverage: splitting the data'))
    if (lastloop == 0) {
        split.len <- numrow
    } else {
        split.len <- rep(chunksize, lastloop)
        split.len.sum <- numrow - sum(split.len)
        if (split.len.sum > 0) {
            split.len <- c(split.len, split.len.sum)
        }
    }
    split.idx <- Rle(seq_len(lastloop + 1), split.len)
    if (!is.null(lowMemDir)) {
        chunks <- split(coverage, split.idx)
        dir.create(lowMemDir, showWarnings = FALSE, recursive = TRUE)
        
        for (i in seq_len(length(chunks))) {
            if(toMatrix) {
                chunkProcessed <- derfinderHelper:::.transformSparseMatrix(
                    data = chunks[[i]], scalefac = scalefac)
            } else {
                chunkProcessed <- chunks[[i]]
            }
            save(chunkProcessed, file = file.path(lowMemDir, 
                paste0('chunk', i, '.Rdata')))
        }
        coverage <- NULL
        coverage.split <- seq_len(lastloop + 1)
    } else {
        coverage.split <- lapply(seq_len(lastloop + 1), function(x) {
            split.idx == x
        })
    }
    
    ## Done =)
    result <- list(coverageProcessed = coverage, 
        mclapplyIndex = coverage.split, position = position,
        meanCoverage = means, groupMeans = groupMeans)
    return(result)
    
}
