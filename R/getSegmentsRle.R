#' Segment a Rle into positive, zero, and negative regions
#'
#' Given two cutoffs, L and U, this function divides a numerical Rle into contiguous parts that are above U, between L and U, and below L. This is wrapper for \link[bumphunter]{getSegments}.
#'
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
#' @references Rafael A. Irizarry, Martin Aryee, Hector Corrada Bravo, Kasper D. Hansen and Harris A. Jaffee. bumphunter: Bump Hunter. R package version 1.1.10.
#' @author Leonardo Collado-Torres
#' @export
#' @importMethodsFrom IRanges quantile as.vector as.numeric
#' @importFrom bumphunter getSegments
#' @examples
#' library("IRanges")
#' set.seed(20130725)
#' pos <- Rle(sample(c(TRUE, FALSE), 1e5, TRUE, prob=c(0.05, 0.95)))
#' cluster <- clusterMakerRle(pos, 100L)
#' data <- Rle(rnorm(sum(pos)))
#' cutoff <- quantile(data, .99)
#'
#' ## It's quite fast
#' system.time(segs <- getSegmentsRle(data, cluster, cutoff, verbose=TRUE))
#' 
#' \dontrun{
#' ## The output is different in look than the one from getSegments() but it's use is similar.
#' ## Plus it can be transformed into the same format as the ouptut from getSegmentsRle().
#' library("bumphunter")
#' foo <- function() {
#' 	segs2 <- getSegments(as.numeric(data), as.integer(cluster), cutoff, assumeSorted=TRUE)
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
#' ## Pretty much the same speed if not a tiny bit faster
#' library("microbenchmark")
#' micro <- microbenchmark(getSegmentsRle(data, cluster, cutoff), foo())
#' micro
#' }

getSegmentsRle <- function(x, f, cutoff = quantile(x, 0.99), verbose = FALSE, zero = TRUE) {
	## Construct segments
    if (verbose) message(paste(Sys.time(), "getSegmentsRle: segmenting with cutoff(s)", paste(cutoff, collapse=", ")))

	## Basically uses getSegments from bumphunter and then constructs the relevant ranges
	Indexes <- getSegments(as.numeric(x), as.vector(f), cutoff=cutoff, assumeSorted=TRUE, verbose=verbose)
	if(!zero) {
		Indexes <- Indexes[c("upIndex", "dnIndex")]
	}
	
	## Built the IRanges object
	if (verbose) message(paste(Sys.time(), "getSegmentsRle: constructing the output"))
	res <- lapply(Indexes, function(ind) {
		if(length(ind) == 0) {
			result <- IRanges()
		} else {
			tmp <- lapply(ind, function(segment) {
				c("start"=min(segment), "end"=max(segment))
			})
			info <- do.call(rbind, tmp)
			result <- IRanges(start=info[,"start"], end=info[,"end"])
		}
		return(result)
	})

	## Done!
    return(res)	
}
