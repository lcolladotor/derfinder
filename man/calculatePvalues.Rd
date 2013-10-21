\name{calculatePvalues}
\alias{calculatePvalues}
\title{Calculate p-values and identify regions}
\usage{
  calculatePvalues(coveragePrep, models, fstats,
    nPermute = 1L,
    seeds = as.integer(gsub("-", "", Sys.Date())) + seq_len(nPermute),
    chr, maxRegionGap = 0L, maxClusterGap = 300L,
    cutoff = quantile(fstats, 0.99),
    mc.cores = getOption("mc.cores", 2L), verbose = TRUE,
    significantCut = c(0.05, 0.1), adjustF = 0)
}
\arguments{
  \item{coveragePrep}{A list with \code{$coverageSplit} and
  \code{$position} normally generated using
  \link{preprocessCoverage}.}

  \item{models}{A list with \code{$mod} and \code{$mod0}
  normally generated using \link{makeModels}.}

  \item{fstats}{A numerical Rle with the F-statistics
  normally generated using \link{calculateStats}.}

  \item{nPermute}{The number of permutations. Note that for
  a full chromosome, a small amount (10) of permutations is
  sufficient. If set to 0, no permutations are performed
  and thus no null regions are used, however, the
  \code{$regions} component is created.}

  \item{seeds}{An integer vector of length \code{nPermute}
  specifying the seeds to be used for each permutation. If
  \code{NULL} no seeds are used.}

  \item{chr}{A single element character vector specifying
  the chromosome name. This argument is passed to
  \link{findRegions}.}

  \item{maxRegionGap}{This argument is passed to
  \link{findRegions}.}

  \item{maxClusterGap}{This argument is passed to
  \link{findRegions}.}

  \item{cutoff}{This argument is passed to
  \link{getSegmentsRle}.}

  \item{mc.cores}{This argument is passed to
  \link[parallel]{mclapply} to run \link{fstats.apply}.}

  \item{verbose}{If \code{TRUE} basic status updates will
  be printed along the way.}

  \item{significantCut}{A vector of length two specifiying
  the cutoffs used to determine significance. The first
  element is used to determine significance for the
  p-values and the second element is used for the
  q-values.}

  \item{adjustF}{A single value to adjust that is added in
  the denominator of the F-stat calculation. Useful when
  the Residual Sum of Squares of the alternative model is
  very small.}
}
\value{
  A list with four components: \describe{ \item{regions }{
  is a GRanges with metadata columns given by
  \link{findRegions} with the additional metadata column
  \code{pvalues}: p-value of the region calculated via
  permutations of the samples; \code{padj}: the qvalues
  calculated using \link[qvalue]{qvalue};
  \code{significant}: whether the p-value is less than 0.05
  (by default); \code{significantPadj}: whether the q-value
  is less than 0.10 (by default). It also includes the mean
  coverage of the region (mean from the mean coverage at
  each base calculated in \link{preprocessCoverage}).
  Furthermore, if \code{groupInfo} was not \code{NULL} in
  \link{preprocessCoverage}, then the group mean coverage
  is calculated as well as the log 2 fold change (using
  group 1 as the reference). } \item{nullStats}{ is a
  numeric Rle with the mean of the null statistics by
  segment.} \item{nullWidths}{ is a numeric Rle with the
  length of each of the segments in the null distribution.
  The area can be obtained by multiplying the absolute
  \code{nullstats} by the corresponding lengths.}
  \item{nullPermutation}{ is a Rle with the permutation
  number from which the null region originated from.} }
}
\description{
  First, this function finds the regions of interest
  according to specified cutoffs. Then it permutes the
  samples and re-calculates the F-statistics. The area of
  the statistics from these segments are then used to
  calculate p-values for the original regions.
}
\examples{
## Calculate library size adjustments
sampleDepths <- sampleDepth(list(genomeData$coverage), prob=0.5, nonzero=TRUE, center=TRUE, verbose=TRUE)
sampleDepths

## Build the models
group <- genomeInfo$pop
adjustvars <- data.frame(genomeInfo$gender)
models <- makeModels(sampleDepths, testvars=group, adjustvars=adjustvars)

## Preprocess the data
## Automatic chunksize used to then compare 1 vs 4 cores in the 'do not run' section
prep <- preprocessCoverage(genomeData, groupInfo=group, cutoff=0, scalefac=32, chunksize=NULL, colsubset=NULL, mc.cores=4)

## Get the F statistics
fstats <- calculateStats(prep, models, mc.cores=1, verbose=TRUE)

## Determine a cutoff from the F-distribution.
## This step is very important and you should consider using quantiles from the observed F statistics
n <- dim(prep$coverageSplit[[1]])[2]
df1 <- dim(models$mod)[2]
df0 <- dim(models$mod0)[2]
cutoff <- qf(0.95, df1-df0, n-df1)

## Low cutoff used for illustrative purposes
cutoff <- 1

## Calculate the p-values and define the regions of interest.
regsWithP <- calculatePvalues(prep, models, fstats, nPermute=10, seeds=NULL, chr="chr21", cutoff=cutoff, mc.cores=1)
regsWithP

## Histogram of the theoretical p-values by region
hist(pf(regsWithP$regions$value, df1-df0, n-df1), main="Distribution original p-values by region", freq=FALSE)

## Histogram of the permutted p-values by region
hist(regsWithP$regions$pvalues, main="Distribution permutted p-values by region", freq=FALSE)

## MA style plot
library("ggplot2")
ma <- data.frame(mean=regsWithP$regions$meanCoverage, foldChange=regsWithP$regions$foldChangeYRIvsCEU)
ggplot(ma, aes(x=log2(mean), y=foldChange)) + geom_point() + ylab("Fold Change (log2)") + xlab("Mean coverage (log2)") + labs(title="MA style plot")

\dontrun{
## Annotate the results
library("bumphunter")
annotation <- annotateNearest(regsWithP$regions, "hg19")
head(annotation)

## Compare speed between 1 and 4 cores (must have them!)
library("microbenchmark")
micro <- microbenchmark(
calculatePvalues(prep, models, fstats, nPermute=10, seeds=NULL, chr="chr21", cutoff=c(2, 5), mc.cores=1, verbose=FALSE),
calculatePvalues(prep, models, fstats, nPermute=10, seeds=NULL, chr="chr21", cutoff=c(2, 5), mc.cores=4, verbose=FALSE),
times=10)
levels(micro$expr) <- c("one", "four")
micro
## Using 4 cores doesn't help with this toy data, but it will (at the expense of more RAM) if you have a larger data set.
}
}
\author{
  Leonardo Collado-Torres
}
\seealso{
  \link{findRegions}, \link{clusterMakerRle},
  \link{getSegmentsRle}, \link{fstats.apply},
  \link[qvalue]{qvalue}
}

