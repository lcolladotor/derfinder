\name{calculatePvalues}
\alias{calculatePvalues}
\title{Calculate p-values and identify regions}
\usage{
  calculatePvalues(statsInfo, nPermute = 1L,
    seeds = as.integer(gsub("-", "", Sys.Date())) + seq_len(nPermute),
    chr, maxGap = 300L,
    cutoff = quantile(abs(statsInfo$fstats), 0.99),
    mc.cores = getOption("mc.cores", 2L), verbose = TRUE)
}
\arguments{
  \item{statsInfo}{A list with \code{$coverageSplit},
  \code{$position}, \code{$fstats}, \code{$mod}, and
  \code{mod0} components as generated using
  \link{calculateStats}.}

  \item{nPermute}{The number of permutations. Note that for
  a full chromosome, a small amount (10) of permutations is
  sufficient.}

  \item{seeds}{An integer vector of length \code{nPermute}
  specifying the seeds to be used for each permutation. If
  \code{NULL} no seeds are used.}

  \item{chr}{A single element character vector specifying
  the chromosome name. This argument is passed to
  \link{findRegions}.}

  \item{maxGap}{This argument is passed to
  \link{clusterMakerRle}.}

  \item{cutoff}{This argument is passed to
  \link{getSegmentsRle}.}

  \item{mc.cores}{This argument is passed to
  \link[parallel]{mclapply} to run \link{fstats.apply}.}

  \item{verbose}{If \code{TRUE} basic status updates will
  be printed along the way.}
}
\value{
  A list with two components. The first one
  --\code{$regions}-- is a GRanges with metadata columns
  given by \link{findRegions} and \describe{ \item{pvalues
  }{ p-value of the region calculated via permutations of
  the samples.} } The second one --\code{$nullstats}-- is a
  numeric Rle with the mean of the null statistics by
  segment.
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
\details{
  Partially based on \link[derfinder]{get.pvals.DF}.
}
\examples{
## Get the statistics
group <- brainInfo$outcome
adjustvars <- brainInfo[, c("sex", "age", "left.hemisph", "pmi", "brainpH")]
statsInfo <- calculateStats(brainData, group, adjustvars=adjustvars, mc.cores=1, verbose=TRUE)

## Determine a cutoff from the F-distribution.
## This step is very important and you should consider using quantiles from the observed F statistics
n <- dim(statsInfo$coverageSplit[[1]])[2]
df1 <- dim(statsInfo$mod)[2]
df0 <- dim(statsInfo$mod0)[2]
cutoff <- qf(0.95, df1-df0, n-df1)

## Calculate the p-values and define the regions of interest.
regsWithP <- calculatePvalues(statsInfo, nPermute=10, seeds=NULL, chr="chr21", cutoff=cutoff, mc.cores=1)
regsWithP

## Histogram of the original F statistics
f.ori <- as.numeric(statsInfo$fstats)
hist(f.ori, main="Distribution original F-stats", freq=FALSE)

## Histogram of the null F statistics
f.null <- as.numeric(regsWithP$nullstats)
hist(f.null, main="Distribution null F-stats by region", freq=FALSE)

## Histogram of the original p-values
hist(pf(f.ori, df1-df0, n-df1), main="Distribution original p-values", freq=FALSE)

## Histogram of the null p-values
hist(pf(f.null, df1-df0, n-df1), main="Distribution null p-values by region", freq=FALSE)

## By region

## Histogram of the original F statistics by region
hist(regsWithP$regions$value, main="Distribution original F-stats by region", freq=FALSE)

## Histogram of the original p-values by region
hist(pf(regsWithP$regions$value, df1-df0, n-df1), main="Distribution original p-values by region", freq=FALSE)

## Histogram of the p-values by region
hist(regsWithP$regions$pvalues, main="Distribution permutted p-values by region", freq=FALSE)

\dontrun{
## Annotate the results
library("bumphunter")
annotation <- annotateNearest(regsWithP$regions, "hg19")
head(annotation)

## Compare speed between 1 and 4 cores (must have them!)
## The chunksize is artifically reduced just to actually need to run mclapply
library("microbenchmark")
micro <- microbenchmark(
calculatePvalues(statsInfo, nPermute=10, seeds=NULL, chr="chr21", cutoff=c(2, 5), mc.cores=1, verbose=FALSE),
calculatePvalues(statsInfo, nPermute=10, seeds=NULL, chr="chr21", cutoff=c(2, 5), mc.cores=4, verbose=FALSE),
times=10)
levels(micro$expr) <- c("one", "four")
micro
## Doesn't seem to help much with this toy data.
}
}
\author{
  Leonardo Collado-Torres
}
\references{
  Frazee et al. Biostatistics in review.
}
\seealso{
  \link{findRegions}, \link{clusterMakerRle},
  \link{getSegmentsRle}, \link{fstats.apply}
}

