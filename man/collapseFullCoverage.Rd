\name{collapseFullCoverage}
\alias{collapseFullCoverage}
\title{Collapse full coverage information for efficient quantile computations}
\usage{
collapseFullCoverage(fullCov, colsubset = NULL, save = FALSE,
  verbose = FALSE)
}
\arguments{
  \item{fullCov}{A list where each element is the result
  from \link[derfinder]{loadCoverage} used with
  \code{cutoff=NULL}. The elements of the list should be
  named according to the chromosome number. Can be
  generated using \link{fullCoverage}.}

  \item{colsubset}{Which columns of
  \code{coverageInfo$coverage} to use.}

  \item{save}{If \code{TRUE}, the result is saved as
  "collapsedFull.Rdata".}

  \item{verbose}{If \code{TRUE} basic status updates will
  be printed along the way.}
}
\value{
A list with one element per sample. Then per sample, a list
with two vector elements: \code{values} and \code{weights}.
The first one is the coverage value and the second one is
the number of bases with that value.
}
\description{
For a given data set this function collapses the full
coverage information for each sample from all the
chromosomes. The resulting information per sample is the
number of bases with coverage 0, 1, etc. It is similar to
using table() on a regular vector. This information is then
used by \link{sampleDepth} for calculating the sample depth
adjustments. The data set can loaded to R using (see
\link{fullCoverage} and optionally filtered using
\link{filterData}.
}
\examples{
## Choose the adjusting variables and define all the parameters for makeModels()
coverageInfo <- genomeData

## Collapse the coverage information
collapsedFull <- collapseFullCoverage(list(coverageInfo$coverage), verbose=TRUE)
collapsedFull
}
\author{
Leonardo Collado-Torres
}
\seealso{
\link{fullCoverage}, \link{sampleDepth}
}

