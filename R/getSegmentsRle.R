#' Segment a Rle into positive, zero, and negative regions
#'
#' Given two cutoffs, L and U, this function divides a numerical Rle into contiguous parts that are above U, between L and U, and below L. This is very similar to \link[bumphunter]{getSegments}.
#'
#' 
#' @param x A numeric Rle.
#' @param f An integer Rle used to pre-divide x into pieces. Each piece is then segmented based on the cutoff. Setting this to NULL says that there is no pre-division. Often, \link{clusterMakerRle} is used to define this integer Rle.
#' @param cutoff A numeric vector of length either 1 or 2. If length is 1, U will be cutoff and L will be -cutoff. Otherwise it specifies L and U. The function will furthermore always use the minimum of cutoff for L and the maximum for U.
#' @param verbose If \code{TRUE} basic status updates will be printed along the way.
#' @param zero If \code{TRUE} the zero index is computed. If \code{FALSE} only the up and down indexes are computed.
#' @param method If \code{speed} then the original \link[bumphunter]{getSegments} is used and the final IRanges is built from the results. If \code{memory} then the function operates in Rle world and attempts to reduce the memory overhead. 
#'
#' @details Note that with an Rle of length 5 millon that resulted in a total of 7319 segments (upIndex and downIndex only), \link{getSegmentsRle} with \code{method="speed"} runs in around 12 seconds versus 24 seconds with \code{method="memory"}. The memory difference was less than 1Gb of RAM. Thus unless you have a much larger data set that has many more segments and are running into memory problems you should use \code{method="speed"}.
#'
#' @return A list of IRanges, one for the up segments, one for the down segments, and if \code{zero} is set to \code{TRUE} then one for the zero segments.
#'
#' @seealso \link[bumphunter]{getSegments}, \link{clusterMakerRle}, \link{findRegions}
#' @references Rafael A. Irizarry, Martin Aryee, Hector Corrada Bravo, Kasper D. Hansen and Harris A. Jaffee. bumphunter: Bump Hunter. R package version 1.1.10.
#' @author Leonardo Collado-Torres
#' @export
#' @importFrom IRanges Rle runValue "runValue<-" nrun IRanges start end Views ranges
#' @importMethodsFrom IRanges quantile length "<=" ">=" "[" "[<-" ">" "<" "==" sapply as.vector as.numeric
#' @importFrom bumphunter getSegments
#' @examples
#' library("IRanges")
#' set.seed(20130725)
#' pos <- Rle(sample(c(TRUE, FALSE), 1e5, TRUE, prob=c(0.05, 0.95)))
#' cluster <- clusterMakerRle(pos, 100L)
#' data <- Rle(rnorm(sum(pos)))
#' cutoff <- quantile(data, .99)
#'
#' ## Time differences between methods
#' system.time(segs <- getSegmentsRle(data, cluster, cutoff, verbose=TRUE, method="speed"))
#' system.time(segs <- getSegmentsRle(data, cluster, cutoff, verbose=TRUE, method="memory"))
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
#' ## getSegmentsRle(method="memory") is slower, yet hopefully less memory heavy.
#' library("microbenchmark")
#' micro <- microbenchmark(getSegmentsRle(data, cluster, cutoff, method="memory"), getSegmentsRle(data, cluster, cutoff, method="speed"), foo())
#' micro
#' }

## This version is slower (on small data) than the original getSegments() but because everything is performed in Rle world, it should require less memory
getSegmentsRle <- function(x, f, cutoff = quantile(x, 0.99), verbose = FALSE, zero = TRUE, method="speed") {
	stopifnot(method %in% c("speed", "memory"))
	
	## Construct segments
    if (verbose) message(paste(Sys.time(), "getSegmentsRle: segmenting with cutoff(s)", paste(cutoff, collapse=", ")))

	if (method == "memory") {
		## Setup steps
		if (is.null(f))  {
			f <- Rle(1L, length(x))
		}
	    stopifnot(length(x) == length(f))
	    stopifnot(length(cutoff) <= 2)		
    
	    if (length(cutoff) == 1) 
	        cutoff <- c(-cutoff, cutoff)
	    cutoff <- sort(cutoff)
		
		## Apply the cutoff
	    direction <- x >= cutoff[2]
	    direction[x <= cutoff[1]] <- -1L
		rm(x)
	
	
		## Find the segments
		segments <- direction
		runValue(segments) <- seq_len(nrun(segments))
	    segments <- segments + f
		rm(f)
	
	
		## Original segments
		segs <- IRanges(start=start(segments), end=end(segments))
		rm(segments)
	
	
		## Construct the output
		if (verbose) message(paste(Sys.time(), "getSegmentsRle: constructing the output"))

		## Define the main states
		up <- direction > 0
		down <- direction < 0
		
		## Construct the individual indexes		
		segs.up <- Views(up, segs)
		upIndex <- ranges(segs.up)[sapply(segs.up, all)]
		rm(up, segs.up)
	
	
		segs.down <- Views(down, segs)
		dnIndex <- ranges(segs.down[sapply(segs.down, all)])
		rm(down, segs.down)
	
	
		## Construct the final output
		if(zero) {
			zero <- direction == 0
			segs.zero <- Views(zero, segs)
			zeroIndex <- ranges(segs.zero[sapply(segs.zero, all)])
			rm(zero, segs.zero)
		
			res <- list(upIndex = upIndex, dnIndex = dnIndex, zeroIndex = zeroIndex)
		} else {
			res <- list(upIndex = upIndex, dnIndex = dnIndex)
		}
	} else {
		## This method basically uses getSegments from bumphunter and then constructs the relevant ranges
		segs2 <- getSegments(as.numeric(x), as.vector(f), cutoff=cutoff, assumeSorted=TRUE, verbose=verbose)
		if(!zero) {
			segs2 <- segs2[c("upIndex", "dnIndex")]
		}
		
		## Built the IRanges object
		if (verbose) message(paste(Sys.time(), "getSegmentsRle: constructing the output"))
		res <- lapply(segs2, function(ind) {
			tmp <- lapply(ind, function(segment) {
				c("start"=min(segment), "end"=max(segment))
			})
			info <- do.call(rbind, tmp)
			IRanges(start=info[,"start"], end=info[,"end"])
		})
	}
	
	
	## Done!
    return(res)	
}
