\name{mergeResults}
\alias{mergeResults}
\title{Merge results from different chromosomes}
\usage{
  mergeResults(chrnums = c(1:22, "X", "Y"), prefix = ".",
    significantCut = c(0.05, 0.1), verbose = TRUE)
}
\arguments{
  \item{chrnums}{The chromosome numbers of the files to be
  merged.}

  \item{prefix}{The main data directory path, which can be
  useful if \link{analyzeChr} is used for several
  parameters and the results are saved in different
  directories.}

  \item{significantCut}{A vector of length two specifiying
  the cutoffs used to determine significance. The first
  element is used to determine significance for the
  p-values and the second element is used for the q-values
  just like in \link{calculatePvalues}.}

  \item{verbose}{If \code{TRUE} basic status updates will
  be printed along the way.}
}
\value{
  Four Rdata files. \describe{ \item{fullFstats.Rdata }{
  Full F-statistics from all chromosomes in a list of Rle
  objects.} \item{fullTime.Rdata }{ Timing information from
  all chromosomes.} \item{fullNullSummary.Rdata}{ A
  DataFrame with the null region information: statistic,
  width, chromosome and permutation identifier. It's
  ordered by the statistics} \item{fullRegions.Rdata}{
  GRanges object with regions found and with full
  annotation from \link[bumphunter]{annotateNearest}. Note
  that the column \code{strand} from
  \link[bumphunter]{annotateNearest} is renamed to
  \code{annoStrand} to comply with GRanges specifications.
  } }
}
\description{
  This function merges the results from running
  \link{analyzeChr} on several chromosomes. It
  re-calculates the p-values and q-values using the pooled
  areas from the null regions from all chromosomes.
}
\examples{
\dontrun{
mergeResults(prefix="run1")
}
}
\author{
  Leonardo Collado-Torres
}
\seealso{
  \link{analyzeChr}, \link{calculatePvalues}
}

