#' Filter the positions of interest
#'
#' For a group of samples this function reads the coverage information for a 
#' specific chromosome directly from the BAM files. It then merges them into a 
#' DataFrame and removes the bases that do not pass the cutoff. This is a 
#' helper function for \link{loadCoverage} and \link{preprocessCoverage}.
#' 
#' @param data Either a list of Rle objects or a DataFrame with the coverage 
#' information.
#' @param cutoff The base-pair level cutoff to use. It's behavior is controlled 
#' by \code{filter}.
#' @param index A logical Rle with the positions of the chromosome that passed 
#' the cutoff. If \code{NULL} it is assumed that this is the first time using 
#' \link{filterData} and thus no previous index exists.
#' @param filter Has to be either \code{'one'} (default) or \code{'mean'}. In 
#' the first case, at least one sample has to have coverage above \code{cutoff}.
#' In the second case, the mean coverage has to be greater than \code{cutoff}.
#' @param totalMapped A vector with the total number of reads mapped for each 
#' sample. The vector should be in the same order as the samples in \code{data}.
#' Providing this data adjusts the coverage to reads in \code{targetSize} 
#' library prior to filtering. See \link{getTotalMapped} for
#' calculating this vector.
#' @param targetSize The target library size to adjust the coverage to. Used
#' only when \code{totalMapped} is specified. By default, it adjusts to 
#' libraries with 80 million reads.
#' @param ... Arguments passed to other methods and/or advanced arguments.
#' Advanced arguments:
#' \describe{
#' \item{verbose }{ If \code{TRUE} basic status updates will be printed along 
#' the way.}
#' \item{returnMean }{ If \code{TRUE} the mean coverage is included in the 
#' result. \code{FALSE} by default.}
#' \item{returnCoverage }{ If \code{TRUE}, the coverage DataFrame is returned.
#' \code{TRUE} by default.}
#' }
#'
#' @return A list with up to three components.
#'
#' \describe{
#' \item{coverage }{ is a DataFrame object where each column represents a 
#' sample. The number of rows depends on the number of base pairs that passed 
#' the cutoff and the information stored is the coverage at that given base. 
#' Included only when \code{returnCoverage = TRUE}.}
#' \item{position }{  is a logical Rle with the positions of the chromosome 
#' that passed the cutoff.}
#' \item{meanCoverage }{ is a numeric Rle with the mean coverage at each base. 
#' Included only when \code{returnMean = TRUE}.}
#' \item{colnames }{ Specifies the column names to be used for the results 
#' DataFrame. If \code{NULL}, names from \code{data} are used.}
#' \item{smoothMean }{ Whether to smooth the mean. Used only when
#' \code{filter = 'mean'}. This option is used internally by 
#' \link{regionMatrix}.}
#' }
#' Passed to the internal function \code{.smootherFstats}, see 
#' \link{findRegions}.
#'
#' @details If \code{cutoff} is \code{NULL} then the data is grouped into 
#' DataFrame without applying any cutoffs. This can be useful if you want to 
#' use \link{loadCoverage} to build the coverage DataFrame without applying any 
#' cutoffs for other downstream purposes like plotting the coverage values of a 
#' given region. You can always specify the \code{colsubset} argument in 
#' \link{preprocessCoverage} to filter the data before calculating the F 
#' statistics.
#'
#' @author Leonardo Collado-Torres
#' @export
#'
#' @importMethodsFrom IRanges '[' '[<-' '[[' colnames 'colnames<-' lapply
#' @import S4Vectors
#' @importFrom methods is
#'
#' @seealso \link{loadCoverage}, \link{preprocessCoverage},
#' \link{getTotalMapped}
#' @examples
#' ## Construct some toy data
#' library('IRanges')
#' x <- Rle(round(runif(1e4, max=10)))
#' y <- Rle(round(runif(1e4, max=10)))
#' z <- Rle(round(runif(1e4, max=10)))
#' DF <- DataFrame(x, y, z)
#'
#' ## Filter the data
#' filt1 <- filterData(DF, 5)
#' filt1
#'
#' ## Filter again but only using the first two samples
#' filt2 <- filterData(filt1$coverage[, 1:2], 5, index=filt1$position)
#' filt2
#'

filterData <- function(data, cutoff = NULL, index = NULL, filter = 'one',
    totalMapped = NULL, targetSize = 80e6, ...) {
        
    ## Check filter
    stopifnot(filter %in% c('one', 'mean'))
    
    ## Advanged arguments
# @param verbose If \code{TRUE} basic status updates will be printed along the 
# way.
    verbose <- .advanced_argument('verbose', TRUE, ...)



# @param returnMean If \code{TRUE} the mean coverage is included in the result.
    returnMean <- .advanced_argument('returnMean', FALSE, ...)



# @param returnCoverage If \code{TRUE}, the coverage DataFrame is returned.
    returnCoverage <- .advanced_argument('returnCoverage', TRUE, ...)


# @param colnames Specifies the column names to be used for the results 
# DataFrame. If \code{NULL}, names from \code{data} are used.
    colnames <- .advanced_argument('colnames', NULL, ...)

    ## Initialize meanCov
    meanCov <- NULL
    
    ## library size adjustments
    if(!is.null(totalMapped) & targetSize != 0) {
        mappedPerXM <- totalMapped / targetSize
        
        ## Normalize to a given library size
        if (verbose) 
            message(paste(Sys.time(), 'filterData: normalizing coverage'))
        data <- mapply(function(x, d) x / d, data, mappedPerXM)
        if (verbose) 
            message(paste(Sys.time(), 'filterData: done normalizing coverage'))
    }
    
    ## If there is no cutoff to apply, just build the DataFrame
    if (is.null(cutoff)) {
        newindex <- NULL
        finalidx <- index
    } else {
        ## Construct the filtering index
        if(filter == 'one') {
            for (i in seq_len(length(data))) {
                if (i == 1) {
                    newindex <- data[[i]] > cutoff
                } else {
                    newindex <- newindex | data[[i]] > cutoff
                }
            }
        } else if (filter == 'mean') {
            meanCov <- Reduce('+', data) / length(data)
            
            smoothMean <- .advanced_argument('smoothMean', FALSE, ...)
            if(smoothMean) {
                if(verbose) message(paste(Sys.time(), 'filterData: smoothing mean coverage'))
                newindex <- Rle(TRUE, length(meanCov))
                meanCov <- derfinder:::.smootherFstats(fstats = meanCov, position = newindex, ...)
            }            
            newindex <- meanCov > cutoff
        }
        
        ## Build the final index
        if (!is.null(index)) {
            finalidx <- index
            finalidx[index] <- newindex
        } else {
            finalidx <- newindex
        }
    }
    
    ## Keep only bases that pass the cutoff
    if(returnMean) {
        if(is.null(meanCov)) {
            ## Calculate the mean if needed
            meanCov <- Reduce('+', data) / length(data)
        }
        if(!is.null(newindex)) {
            meanCovFiltered <- meanCov[newindex]
        } else {
            meanCovFiltered <- meanCov
        }
    }
    
    if(returnCoverage) {
        if (is(data, 'DataFrame')) {
            if (!is.null(newindex)) {
                DF <- data[newindex, ]
            } else {
                DF <- data
            }
        } else {
            ## Subset the data and group into DataFrame
            if (!is.null(newindex)) {
                DF <- DataFrame(lapply(data, function(x) {
                    x[newindex]
                }), check.names = FALSE)
            } else {
                DF <- DataFrame(data, check.names = FALSE)
            }
        }
    }
        
    ## Info for the user
    if (verbose) {
        if(returnCoverage) {
            message(paste(Sys.time(), 'filterData: originally there were', 
                length(data[[1]]), 'rows, now there are', nrow(DF), 
                'rows. Meaning that', 100 - round(nrow(DF)/length(data[[1]]) * 
                    100, 2), 'percent was filtered.'))
        } else if (returnMean) {
            message(paste(Sys.time(), 'filterData: originally there were', 
                length(data[[1]]), 'rows, now there are', 
                length(meanCovFiltered), 'rows. Meaning that', 100 - 
                round(length(meanCovFiltered)/length(data[[1]]) * 
                    100, 2), 'percent was filtered.'))
        }
    }
    
    ## Assign column names
    if(returnCoverage) {
        if (!is.null(colnames)) {
            colnames(DF) <- colnames
        }
    }
    
    ## Make the final resulting object.
    if(returnMean & returnCoverage) {
        res <- list(coverage = DF, position = finalidx, 
            meanCoverage = meanCovFiltered)   
    } else if (!returnMean & returnCoverage){
        res <- list(coverage = DF, position = finalidx)
    } else if (returnMean & !returnCoverage) {
        res <- list(position = finalidx, meanCoverage = meanCovFiltered)
    } else {
        res <- list(position = finalidx)
    }
    
    return(res)
}
