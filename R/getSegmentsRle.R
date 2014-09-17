#' Segment a Rle into positive, zero, and negative regions
#'
#' Given two cutoffs, L and U, this function slices a numerical Rle into up and 
#' down sections. It is a wrapper for \link[IRanges]{slice} with functionality 
#' inspired from \link[bumphunter]{getSegments}.
#'
#' 
#' @param x A numeric Rle.
#' @param cutoff A numeric vector of length either 1 or 2. If length is 1, U 
#' will be cutoff and L will be -cutoff. Otherwise it specifies L and U. The 
#' function will furthermore always use the minimum of cutoff for L and the 
#' maximum for U.
#' @param verbose If \code{TRUE} basic status updates will be printed along the 
#' way.
#'
#' @return A list of IRanges objects, one for the up segments and one for the 
#' down segments.
#'
#' @seealso \link[bumphunter]{getSegments}, \link[IRanges]{slice}, 
#' \link{clusterMakerRle}, \link{findRegions}
#'
#' @author Leonardo Collado-Torres
#' @export
#' @aliases get_segments_rle
#' @importMethodsFrom IRanges quantile
#' @importFrom IRanges slice
#' @importMethodsFrom S4Vectors as.numeric
#' @examples
#' library("IRanges")
#' set.seed(20130725)
#' pos <- Rle(sample(c(TRUE, FALSE), 1e5, TRUE, prob=c(0.05, 0.95)))
#' data <- Rle(rnorm(sum(pos)))
#' cutoff <- quantile(data, .99)
#'
#' ## It's quite fast
#' system.time(segs <- getSegmentsRle(data, cutoff, verbose=TRUE))
#' 
#' \dontrun{
#' ## The output is different in look than the one from getSegments() but it's 
#' ## use is similar.
#' ## Plus it can be transformed into the same format as the ouptut from 
#' ## getSegmentsRle().
#' library("bumphunter")
#' cluster <- clusterMakerRle(pos, 100L)
#' foo <- function() {
#'     segs2 <- getSegments(as.numeric(data), as.integer(cluster), cutoff, 
#'     assumeSorted=TRUE)[c("upIndex", "dnIndex")]
#'     segs.ir <- lapply(segs2, function(ind) {
#'         tmp <- lapply(ind, function(segment) {
#'             c("start"=min(segment), "end"=max(segment))
#'         })
#'         info <- do.call(rbind, tmp)
#'         IRanges(start=info[,"start"], end=info[,"end"])
#'     })
#'     return(segs.ir)
#' }
#' identical(foo(), segs) 
#'
#' }
#'

getSegmentsRle <- function(x, cutoff = quantile(x, 0.99), verbose = FALSE) {
    
    ## Select the cutoff
    if (verbose) message(paste(Sys.time(),
        "getSegmentsRle: segmenting with cutoff(s)",
        paste(cutoff, collapse=", ")))
    stopifnot(length(cutoff) <= 2)
    if (length(cutoff) == 1) {
        cutoff <- c(-cutoff, cutoff)
    }
    cutoff <- sort(cutoff)
    
    ## Find the segments
    result <- lapply(c("upIndex", "dnIndex"), function(ind) {
        if(ind == "upIndex") {
            fcut <- slice(x=x, lower=cutoff[2], rangesOnly=TRUE)
        } else {
            fcut <- slice(x=x, upper=cutoff[1], rangesOnly=TRUE)
        }
        return(fcut)
    })
    names(result) <- c("upIndex", "dnIndex")

    ## Done!
    return(result)
}

#' @export
get_segments_rle <- getSegmentsRle
