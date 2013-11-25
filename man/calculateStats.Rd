\name{calculateStats}
\alias{calculateStats}
\title{Calculate F-statistics at base pair resolution from a loaded BAM files}
\usage{
  calculateStats(coveragePrep, models,
    mc.cores = getOption("mc.cores", 2L), adjustF = 0,
    verbose = TRUE)
}
\arguments{
  \item{coveragePrep}{A list with
  \code{$coverageProcessed}, \code{$mclapplyIndex}, and
  \code{$position} normally generated using
  \link{preprocessCoverage}.}

  \item{models}{A list with \code{$mod} and \code{$mod0}
  normally generated using \link{makeModels}.}

  \item{mc.cores}{This argument is passed to
  \link[parallel]{mclapply} to run \link{fstats.apply}.}

  \item{adjustF}{A single value to adjust that is added in
  the denominator of the F-stat calculation. Useful when
  the Residual Sum of Squares of the alternative model is
  very small.}

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

