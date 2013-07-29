\name{calculateStats}
\alias{calculateStats}
\title{Calculate F-statistics at base pair resolution from a loaded BAM files}
\usage{
  calculateStats(coveragePrep, models,
    mc.cores = getOption("mc.cores", 2L), verbose = TRUE)
}
\arguments{
  \item{coveragePrep}{A list with \code{$coverageSplit} and
  \code{$position} normally generated using
  \link{preprocessCoverage}.}

  \item{models}{A list with \code{$mod} and \code{$mod0}
  normally generated using \link{makeModels}.}

  \item{mc.cores}{This argument is passed to
  \link[parallel]{mclapply} to run \link{fstats.apply}.}

  \item{verbose}{If \code{TRUE} basic status updates will
  be printed along the way.}
}
\value{
  A numeric Rle with the F-statistics per base pair that
  passed the cutoff.
}
\description{
  After defining the models of interest (see
  \link{makeModels}) and pre-processing the data (see
  \link{preprocessCoverage}), use \link{calculateStats} to
  calculate the F-statistics at base-pair resolution.
}
\examples{
## Construct the models
group <- genomeInfo$pop
adjustvars <- data.frame(genomeInfo$gender)
models <- makeModels(coverageInfo=genomeData, group=group, adjustvars=adjustvars, nonzero=TRUE)

## Preprocess the data
prep <- preprocessCoverage(genomeData, cutoff=0, scalefac=32, chunksize=1e3, colsubset=NULL)

## Run the function
fstats <- calculateStats(prep, models, mc.cores=1, verbose=TRUE)
fstats
}
\author{
  Leonardo Collado-Torres
}
\seealso{
  \link{makeModels}, \link{preprocessCoverage}
}

