\name{getSegmentsRle}
\alias{getSegmentsRle}
\title{Segment a Rle into positive, zero, and negative regions}
\usage{
  getSegmentsRle(x, f, cutoff = quantile(abs(x), 0.99),
    verbose = FALSE, zero = TRUE)
}
\arguments{
  \item{x}{A numeric Rle.}

  \item{f}{An integer Rle used to pre-divide x into pieces.
  Each piece is then segmented based on the cutoff. Setting
  this to NULL says that there is no pre-division. Often,
  \link{clusterMakerRle} is used to define this integer
  Rle.}

  \item{cutoff}{A numeric vector of length either 1 or 2.
  If length is 1, U will be cutoff and L will be -cutoff.
  Otherwise it specifies L and U. The function will
  furthermore always use the minimum of cutoff for L and
  the maximum for U.}

  \item{verbose}{If \code{TRUE} basic status updates will
  be printed along the way.}

  \item{zero}{If \code{TRUE} the zero index is computed. If
  \code{FALSE} only the up and down indexes are computed.}
}
\value{
  A list of IRanges, one for the up segments, one for the
  down segments, and if \code{zero} is set to \code{TRUE}
  then one for the zero segments.
}
\description{
  Given two cutoffs, L and U, this function divides a
  numerical Rle into contiguous parts that are above U,
  between L and U, and below L. This is very similar to
  \link[bumphunter]{getSegments}.
}
\details{
  \link[bumphunter]{getSegments} adapted to Rle world.

  It is a a helper function for \link{findRegions}.
}
\examples{
library("IRanges")
set.seed(20130725)
pos <- Rle(sample(c(TRUE, FALSE), 1e5, TRUE, prob=c(0.05, 0.95)))
cluster <- clusterMakerRle(pos, 100L)
data <- Rle(rnorm(sum(pos)))
segs <- getSegmentsRle(data, cluster, verbose=TRUE)
segs

\dontrun{
## The output is different in look than the one from getSegments() but it's use is similar.
## Plus it can be transformed into the same format as the ouptut from getSegmentsRle().
library("bumphunter")
foo <- function() {
	segs2 <- getSegments(as.numeric(data), as.integer(cluster))
	segs.ir <- lapply(segs2, function(ind) {
		tmp <- lapply(ind, function(segment) {
			c("start"=min(segment), "end"=max(segment))
		})
		info <- do.call(rbind, tmp)
		IRanges(start=info[,"start"], end=info[,"end"])
	})
	return(segs.ir)
}
identical(foo(), segs)

## getSegmentsRle() is slower, yet hopefully less memory intense.
library("microbenchmark")
micro <- microbenchmark(getSegmentsRle(data, cluster), foo())
micro
}
}
\author{
  Leonardo Collado-Torres
}
\seealso{
  \link[bumphunter]{getSegments}, \link{clusterMakerRle},
  \link{findRegions}
}

