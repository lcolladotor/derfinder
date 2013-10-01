\name{sampleDepth}
\alias{sampleDepth}
\title{Calculate adjustments for library size}
\usage{
  sampleDepth(fullCov, prob = 0.9, nonzero = TRUE,
    center = FALSE, colsubset = NULL, verbose = FALSE)
}
\arguments{
  \item{fullCov}{A list where each element is the result
  from \link[derfinder]{loadCoverage} used with
  \code{cutoff=NULL}. The elements of the list should be
  named according to the chromosome number. Can be
  generated using \link{fullCoverage}.}

  \item{prob}{A number between 0 and 1 representing the
  quantile of interest. For example, 0.5 is the median.}

  \item{nonzero}{If \code{TRUE} only the nonzero counts are
  used to calculate the library size adjustment.}

  \item{center}{If \code{TRUE} the sample quantiles are
  centered by the mean of the sample quantiles across
  samples. This can be helpful for interpretation
  purposes.}

  \item{colsubset}{Which colums of
  \code{coverageInfo$coverage} to use.}

  \item{verbose}{If \code{TRUE} basic status updates will
  be printed along the way.}
}
\value{
  A vector with the library size depth adjustments per
  sample to be used in \link{makeModels}.
}
\description{
  For a given data set calculate the per-sample Q quantile
  of the base coverage. The data set can loaded to R using
  (see \link{fullCoverage} and optionally filtered using
  \link{filterData}. One the per-sample Q quantile has been
  calculated, you could be interested in centering them
  accross samples for interpretability concerns. This
  information is then used in \link{makeModels} for
  construcing the null and alternative models.
}
\examples{
## Choose the adjusting variables and define all the parameters for makeModels()
coverageInfo <- genomeData

## Calculate library size adjustments
sampleDepths <- sampleDepth(list(coverageInfo$coverage), prob=0.5, nonzero=TRUE, center=TRUE, verbose=TRUE)
sampleDepths
}
\author{
  Leonardo Collado-Torres
}
\seealso{
  \link{fullCoverage}, \link{makeModels}
}

