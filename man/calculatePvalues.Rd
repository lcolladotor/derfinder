\name{calculatePvalues}
\alias{calculatePvalues}
\title{Calculate p-values and identify regions}
\usage{
  calculatePvalues(statsInfo, nPermute = 1L,
    seeds = as.integer(gsub("-", "", Sys.Date())) + seq_len(nPermute),
    chunksize = 5e+06, chr, maxGap = 300L,
    cutoff = quantile(abs(statsInfo$fstats), 0.99),
    mc.cores = getOption("mc.cores", 2L), verbose = TRUE)
}
\arguments{
  \item{statsInfo}{A list with \code{$coverage},
  \code{$position}, \code{$fstats}, \code{$mod}, and
  \code{mod0} components as generated using
  \link{calculateStats}.}

  \item{nPermute}{The number of permutations. Note that for
  a full chromosome, a small amount (10) of permutations is
  sufficient.}

  \item{seeds}{An integer vector of length \code{nPermute}
  specifying the seeds to be used for each permutation. If
  \code{NULL} no seeds are used.}

  \item{chunksize}{How many rows of
  \code{statsInfo$coverage} should be processed at a time?}

  \item{chr}{A single element character vector specifying
  the chromosome name. This argument is passed to
  \link{findRegions}.}

  \item{maxGap}{This argument is passed to
  \link{clusterMakerRle}.}

  \item{cutoff}{This argument is passed to
  \link{getSegmentsRle}.}

  \item{mc.cores}{This argument is passed to
  \link[parallel]{mclapply}) to run \link{fstats.apply}.}

  \item{verbose}{If \code{TRUE} basic status updates will
  be printed along the way.}
}
\value{
  A GRanges with metadata columns given by
  \link{findRegions} and \describe{ \item{pvalues }{
  p-value of the region calculated via permutations of the
  samples.} }
}
\description{
  First, this function clusters the genomic positions and
  finds the regions of interest according to specified
  cutoffs. Then it permutes the samples and re-calculates
  the F-statistics. The F-statistics are segmented using
  the original clusters and cutoffs. The mean of the
  statistics from these segments are then used to calculate
  p-values for the original regions.
}
\examples{
## Get the statistics
group <- brainInfo$outcome
adjustvars <- brainInfo[, c("sex", "age", "left.hemisph", "pmi", "brainpH")]
statsInfo <- calculateStats(brainData, group, adjustvars=adjustvars, mc.cores=1, verbose=TRUE)
## Calculate the p-values and define the regions of interest.
regsWithP <- calculatePvalues(statsInfo, nPermute=10, seeds=NULL, chr="chr21", cutoff=c(2, 5), mc.cores=1)
regsWithP
hist(regsWithP$pvalues)

\dontrun{
## Annotate the results
library("bumphunter")
annotation <- annotateNearest(regsWithP, "hg19")
head(annotation)

## Compare speed between 1 and 4 cores (must have them!)
## The chunksize is artifically reduced just to actually need to run mclapply
library("microbenchmark")
micro <- microbenchmark(
calculatePvalues(statsInfo, nPermute=10, seeds=NULL, chr="chr21", cutoff=c(2, 5), chunksize=1e3, mc.cores=1, verbose=FALSE),
calculatePvalues(statsInfo, nPermute=10, seeds=NULL, chr="chr21", cutoff=c(2, 5), chunksize=1e3, mc.cores=4, verbose=FALSE),
times=10)
levels(micro$expr) <- c("one", "four")
micro
## Doesn't seem to help with this toy data.
}
}
\author{
  Leonardo Collado-Torres
}
\seealso{
  \link{findRegions}, \link{clusterMakerRle},
  \link{getSegmentsRle}, \link{fstats.apply}
}

