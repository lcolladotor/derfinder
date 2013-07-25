#' Segment a Rle into positive, zero, and negative regions
#'
#' Given two cutoffs, L and U, this function divides a numerical Rle into contiguous parts that are above U, between L and U, and below L. This is very similar to \link[bumphunter]{getSegments}.
#'
#' @details
#' \link[bumphunter]{getSegments} adapted to Rle world.
#'
#' It is a a helper function for \link{findRegions}.
#' 
#' @param x A numeric Rle.
#' @param f An integer Rle used to pre-divide x into pieces. Each piece is then segmented based on the cutoff. Setting this to NULL says that there is no pre-division. Often, \link{clusterMakerRle} is used to define this integer Rle.
#' @param cutoff A numeric vector of length either 1 or 2. If length is 1, U will be cutoff and L will be -cutoff. Otherwise it specifies L and U. The function will furthermore always use the minimum of cutoff for L and the maximum for U.
#' @param verbose If \code{TRUE} basic status updates will be printed along the way.
#' @param zero If \code{TRUE} the zero index is computed. If \code{FALSE} only the up and down indexes are computed.
#'
#' @return A list of IRanges, one for the up segments, one for the down segments, and if \code{zero} is set to \code{TRUE} then one for the zero segments.
#'
#' @seealso \link[bumphunter]{getSegments}, \link{clusterMakerRle}, \link{findRegions}
#' @author Leonardo Collado-Torres
#' @export
#' @importFrom IRanges Rle runValue "runValue<-" nrun IRanges start end Views ranges
#' @importMethodsFrom IRanges quantile length "<=" ">=" "[" "[<-" ">" "<" "==" sapply
#' @examples
#' library("IRanges")
#' set.seed(20130725)
#' pos <- Rle(sample(c(TRUE, FALSE), 1e5, TRUE, prob=c(0.05, 0.95)))
#' cluster <- clusterMakerRle(pos, 100L)
#' data <- Rle(rnorm(sum(pos)))
#' segs <- getSegmentsRle(data, cluster, verbose=TRUE)
#' segs
#' 
#' \dontrun{
#' ## The output is different in look than the one from getSegments() but it's use is similar.
#' ## Plus it can be transformed into the same format as the ouptut from getSegmentsRle().
#' library("bumphunter")
#' foo <- function() {
#' 	segs2 <- getSegments(as.numeric(data), as.integer(cluster))
#' 	segs.ir <- lapply(segs2, function(ind) {
#' 		tmp <- lapply(ind, function(segment) {
#' 			c("start"=min(segment), "end"=max(segment))
#' 		})
#' 		info <- do.call(rbind, tmp)
#' 		IRanges(start=info[,"start"], end=info[,"end"])
#' 	})
#' 	return(segs.ir)
#' }
#' identical(foo(), segs) 
#'
#' ## getSegmentsRle() is slower, yet hopefully less memory intense.
#' library("microbenchmark")
#' micro <- microbenchmark(getSegmentsRle(data, cluster), foo())
#' micro
#' }

## This version is slower (on small data) than the original getSegments() but because everything is performed in Rle world, it should require less memory
getSegmentsRle <- function(x, f, cutoff = quantile(abs(x), 0.99), verbose = FALSE, zero = TRUE) {
	## Setup steps
	if (is.null(f))  {
		f <- Rle(1L, length(x))
	}
    stopifnot(length(x) == length(f))
    stopifnot(length(cutoff) <= 2)
    
    if (length(cutoff) == 1) 
        cutoff <- c(-cutoff, cutoff)
    cutoff <- sort(cutoff)
	
	## Construct segments
    if (verbose) 
        message("getSegmentsRle: segmenting")

	## Apply the cutoff
    direction <- x >= cutoff[2]
    direction[x <= cutoff[1]] <- -1L
	
	## Find the segments
	segments <- direction
	runValue(segments) <- seq_len(nrun(segments))
    segments <- segments + f
	
	## Original segments
	segs <- IRanges(start=start(segments), end=end(segments))
	
	## Construct the output
	if (verbose)
		message("getSegmentsRle: constructing the output")

	## Define the main states
	up <- direction > 0
	down <- direction < 0
		
	## Construct the individual indexes		
	segs.up <- Views(up, segs)
	upIndex <- ranges(segs.up)[sapply(segs.up, all)]
	
	segs.down <- Views(down, segs)
	dnIndex <- ranges(segs.down[sapply(segs.down, all)])
	
	## Construct the final output
	if(zero) {
		zero <- direction == 0
		segs.zero <- Views(zero, segs)
		zeroIndex <- ranges(segs.zero[sapply(segs.zero, all)])
		res <- list(upIndex = upIndex, dnIndex = dnIndex, zeroIndex = zeroIndex)
	} else {
		res <- list(upIndex = upIndex, dnIndex = dnIndex)
	}
	
	## Done!
    return(res)	
}
