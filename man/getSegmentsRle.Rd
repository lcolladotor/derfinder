\name{getSegmentsRle}
\alias{getSegmentsRle}
\title{Segment a Rle into positive, zero, and negative regions}
\usage{
getSegmentsRle(x, cutoff = quantile(x, 0.99), verbose = FALSE)
}
\arguments{
  \item{x}{A numeric Rle.}

  \item{cutoff}{A numeric vector of length either 1 or 2.
  If length is 1, U will be cutoff and L will be -cutoff.
  Otherwise it specifies L and U. The function will
  furthermore always use the minimum of cutoff for L and
  the maximum for U.}

  \item{verbose}{If \code{TRUE} basic status updates will
  be printed along the way.}
}
\value{
A list of IRanges objects, one for the up segments and one
for the down segments.
}
\description{
Given two cutoffs, L and U, this function slices a
numerical Rle into up and down sections. It is a wrapper
for \link[IRanges]{slice} with functionality inspired from
\link[bumphunter]{getSegments}.
}
\examples{
library("IRanges")
set.seed(20130725)
pos <- Rle(sample(c(TRUE, FALSE), 1e5, TRUE, prob=c(0.05, 0.95)))
data <- Rle(rnorm(sum(pos)))
cutoff <- quantile(data, .99)

## It's quite fast
system.time(segs <- getSegmentsRle(data, cutoff, verbose=TRUE))

\dontrun{
## The output is different in look than the one from getSegments() but it's
## use is similar.
## Plus it can be transformed into the same format as the ouptut from
## getSegmentsRle().
library("bumphunter")
cluster <- clusterMakerRle(pos, 100L)
foo <- function() {
    segs2 <- getSegments(as.numeric(data), as.integer(cluster), cutoff,
    assumeSorted=TRUE)[c("upIndex", "dnIndex")]
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

## It's a bit faster
library("microbenchmark")
micro <- microbenchmark(getSegmentsRle(data, cutoff), foo())
micro
}
}
\author{
Leonardo Collado-Torres
}
\seealso{
\link[bumphunter]{getSegments}, \link[IRanges]{slice},
\link{clusterMakerRle}, \link{findRegions}
}

