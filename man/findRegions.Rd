\name{findRegions}
\alias{findRegions}
\title{Find non-zero regions in a Rle}
\usage{
  findRegions(position = NULL, fstats, chr,
    oneTable = TRUE, maxRegionGap = 0L,
    maxClusterGap = 300L, cutoff = quantile(fstats, 0.99),
    segmentIR = NULL, basic = FALSE, verbose = TRUE)
}
\arguments{
  \item{position}{A logical Rle of genomic positions. This
  is generated in \link{loadCoverage}. Note that it gets
  updated in \link{preprocessCoverage} if \code{colsubset}
  is not \code{NULL}.}

  \item{fstats}{A numeric Rle with the F-statistics.
  Normally obtained using \link{calculateStats}.}

  \item{chr}{A single element character vector specifying
  the chromosome name.}

  \item{oneTable}{If \code{TRUE} only one results GRanges
  is returned. Otherwise, a GRangesList with two components
  is returned: one for the regions with positive values and
  one for the negative values.}

  \item{maxRegionGap}{This determines the maximum number of
  gaps between two genomic positions to be considered part
  of the same candidate Differentially Expressed Region
  (candidate DER).}

  \item{maxClusterGap}{This determines the maximum gap
  between candidate DERs. It should be greater than
  \code{maxRegionGap}.}

  \item{cutoff}{This argument is passed to
  \link{getSegmentsRle}.}

  \item{segmentIR}{An IRanges object with the genomic
  positions that are potentials DERs, normally given by
  \link{clusterMakerRle} with \code{maxGap=maxRegionGap}
  and \code{ranges=TRUE}. This is used in
  \link{calculatePvalues} to speed up permutation
  calculations.}

  \item{basic}{If \code{TRUE} a DataFrame is returned that
  has only basic information on the candidate DERs. This is
  used in \link{calculatePvalues} to speed up permutation
  calculations.}

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
  (or below) predefined thresholds. In other words, this
  function finds the candidate Differentially Expressed
  Regions (candidate DERs). This is similar to
  \link[bumphunter]{regionFinder} and is a helper function
  for \link{calculatePvalues}.
}
\details{
  \link[bumphunter]{regionFinder} adapted to Rle world.
}
\examples{
## Collapse the coverage information
collapsedFull <- collapseFullCoverage(list(genomeData$coverage), verbose=TRUE)

## Calculate library size adjustments
sampleDepths <- sampleDepth(collapsedFull, probs=c(0.5), nonzero=TRUE, verbose=TRUE)

## Build the models
group <- genomeInfo$pop
adjustvars <- data.frame(genomeInfo$gender)
models <- makeModels(sampleDepths, testvars=group, adjustvars=adjustvars)

## Preprocess the data
prep <- preprocessCoverage(genomeData, cutoff=0, scalefac=32, chunksize=1e3, colsubset=NULL)

## Get the F statistics
fstats <- calculateStats(prep, models, mc.cores=1, verbose=TRUE)

## Find the regions
regs <- findRegions(prep$position, fstats, "chr21", verbose=TRUE)
regs

\dontrun{
## Compare vs bumphunter
library("bumphunter")
regs2 <- regionFinder(as.numeric(fstats), rep("chr21", length(fstats)), which(prep$position), cluster=NULL, assumeSorted=TRUE, verbose=TRUE, order=FALSE, maxGap=1)
regs2
## Note that regs$L can be calculated with width(regs)
identical(width(regs), as.integer(regs2$L))
## Time comparison
library("microbenchmark")
micro <- microbenchmark(findRegions(prep$position, fstats, "chr21", verbose=FALSE, basic=TRUE), regionFinder(as.numeric(fstats), rep("chr21", length(fstats)), which(prep$position), cluster=NULL, assumeSorted=TRUE, verbose=FALSE, order=FALSE, maxGap=1))
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
\references{
  Rafael A. Irizarry, Martin Aryee, Hector Corrada Bravo,
  Kasper D. Hansen and Harris A. Jaffee. bumphunter: Bump
  Hunter. R package version 1.1.10.
}
\seealso{
  \link[bumphunter]{regionFinder}, \link{calculatePvalues}
}

