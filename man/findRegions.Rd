\name{findRegions}
\alias{findRegions}
\title{Find non-zero regions in a Rle}
\usage{
  findRegions(statsInfo, chr, y = statsInfo$fstats,
    oneTable = TRUE, maxGap = 300L,
    cutoff = quantile(abs(statsInfo$fstats), 0.99),
    verbose = TRUE)
}
\arguments{
  \item{statsInfo}{A list with \code{$position} and
  \code{$fstats} components where the first one is a
  logical Rle of genomic positions and the second one is a
  numeric Rle.}

  \item{chr}{A single element character vector specifying
  the chromosome name.}

  \item{y}{A numeric Rle of the same length as
  \code{statsInfo$fstats} containing values to be averaged
  for the region.}

  \item{oneTable}{If \code{TRUE} only one results GRanges
  is returned. Otherwise, a GRangesList with two components
  is returned: one for the regions with positive values and
  one for the negative values.}

  \item{maxGap}{This argument is passed to
  \link{clusterMakerRle}.}

  \item{cutoff}{This argument is passed to
  \link{getSegmentsRle}.}

  \item{verbose}{If \code{TRUE} basic status updates will
  be printed along the way.}
}
\value{
  Either a GRanges or a GRangesList as determined by
  \code{oneTable}. Each of them has the following metadata
  variables. \describe{ \item{value }{ The mean of the
  values of \code{y} for the given region.} \item{area }{
  The absolute value of the sum of the values of \code{y}
  for the given region.} \item{indexStart }{ The start
  position of the region in terms of the index for
  \code{y}.} \item{indexEnd }{ The end position of the
  region in terms of the index for \code{y}.} \item{cluster
  }{ The cluser ID.} \item{clusterL }{ The total length of
  the cluster.} }
}
\description{
  Find genomic regions for which a numeric vector is above
  (or below) predefined thresholds. This is similar to
  \link[bumphunter]{regionFinder}.
}
\details{
  \link[bumphunter]{regionFinder} adapted to Rle world.
}
\examples{
## Get the statistics
group <- brainInfo$outcome
adjustvars <- brainInfo[, c("sex", "age", "left.hemisph", "pmi", "brainpH")]
stats <- calculateStats(brainData, group, adjustvars=adjustvars, mc.cores=1, verbose=TRUE)
regs <- findRegions(stats, "chr21", verbose=TRUE)
regs

\dontrun{
## Compare vs bumphunter
library("bumphunter")
regs2 <- regionFinder(as.numeric(stats$fstats), rep("chr21", length(stats$fstats)), which(stats$position), cluster=NULL, assumeSorted=TRUE, verbose=TRUE, order=FALSE)
regs2
## Note that regs$L can be calculated with width(regs)
identical(width(regs), as.integer(regs2$L))
## Time comparison
library("microbenchmark")
micro <- microbenchmark(findRegions(stats, "chr21", verbose=FALSE), regionFinder(as.numeric(stats$fstats), rep("chr21", length(stats$fstats)), which(stats$position), cluster=NULL, assumeSorted=TRUE, verbose=FALSE, order=FALSE))
levels(micro$expr) <- c("new", "original")
micro
## The bumphunter function regionFinder() is faster in small data sets.

## Once you have the regions you can proceed to annotate them
annotation <- annotateNearest(regs, "hg19")
annotation
}
}
\author{
  Leonardo Collado-Torres
}
\seealso{
  \link[bumphunter]{regionFinder}
}

