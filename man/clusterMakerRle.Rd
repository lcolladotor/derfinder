\name{clusterMakerRle}
\alias{clusterMakerRle}
\title{Make clusters of genomic locations based on distance in Rle() world}
\usage{
  clusterMakerRle(pos, maxGap = 300L)
}
\arguments{
  \item{pos}{A logical Rle indicating the chromosome
  positions.}

  \item{maxGap}{An integer. Genomic locations within
  \code{maxGap} from each other are placed into the same
  cluster.}
}
\value{
  An integer Rle with the cluster IDs.
}
\description{
  Genomic locations are grouped into clusters based on
  distance: locations that are close to each other are
  assigned to the same cluster. The operation is performed
  on each chromosome independently. This is very similar to
  \link[bumphunter]{clusterMaker}.
}
\details{
  \link[bumphunter]{clusterMaker} adapted to Rle world.
  Assumes that the data is sorted and that everything is in
  a single chromosome. It is also almost as fast as the
  original version with the advantage that everything is in
  Rle() world.

  It is a a helper function for \link{findRegions}.
}
\examples{
library("IRanges")
set.seed(20130725)
pos <- Rle(sample(c(TRUE, FALSE), 1e5, TRUE, prob=c(0.05, 0.95)))
cluster <- clusterMakerRle(pos, 100L)
cluster

\dontrun{
## clusterMakerRle() is comparable in speed if you start from the Rle world.
library("bumphunter")
library("microbenchmark")
micro <- microbenchmark(clusterMakerRle(pos, 100L), clusterMaker(chr=rep("chr21", sum(pos)), pos=which(pos)))
micro
}
}
\author{
  Leonardo Collado-Torres
}
\seealso{
  \link[bumphunter]{clusterMaker}, \link{findRegions}
}

