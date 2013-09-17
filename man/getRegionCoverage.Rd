\name{getRegionCoverage}
\alias{getRegionCoverage}
\title{Extract coverage information for a set of regions}
\usage{
  getRegionCoverage(fullCov, regions,
    calculateMeans = TRUE, verbose = TRUE)
}
\arguments{
  \item{fullCov}{A list where each element is the result
  from \link{loadCoverage} used with \code{cutoff=NULL}.
  The elements of the list should be named according to the
  chromosome number. Can be generated using
  \link{fullCoverage}.}

  \item{regions}{The \code{$regions} output from
  \link{calculatePvalues}. It is important that the
  seqlengths information is provided.}

  \item{calculateMeans}{If \code{TRUE} the mean coverage
  per sample for each region is calculated.}

  \item{verbose}{If \code{TRUE} basic status updates will
  be printed along the way.}
}
\value{
  A list with elements \code{coverageData} and
  \code{coverageMeans} (only if
  \code{calculateMeans=TRUE}). \describe{
  \item{coverageData }{This is a list of data.frame where
  each data.frame has the coverage information (nrow =
  width of region, ncol = number of samples) for a given
  region. The names of the list correspond to the region
  indexes in \code{regions}.} \item{coverageMeans }{This is
  a matrix (nrow = number of regions, ncol = number of
  samples) with the mean coverage per sample for all the
  regions.} }
}
\description{
  This function extracts the raw coverage information
  calculated by \link{fullCoverage} at each base for a set
  of regions found with \link{calculatePvalues}. It can
  further calculate the mean coverage per sample for each
  region.
}
\examples{
\dontrun{
## Construct the models
group <- genomeInfo$pop
adjustvars <- data.frame(genomeInfo$gender)
models <- makeModels(coverageInfo=genomeData, testvars=group, adjustvars=adjustvars, nonzero=TRUE)

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

## Obtain fullCov object
datadir <- system.file("extdata", "genomeData", package="derfinder2")
dirs <- makeBamList(datadir=datadir, samplepatt="*accepted_hits.bam$", bamterm=NULL)
## Shorten the column names
names(dirs) <- gsub("_accepted_hits.bam", "", names(dirs))

## Reading the data and filtering it is quite fast.
fullCov <- fullCoverage(dirs=dirs, chrnums="21", mc.cores=1)

## Assign chr lengths using hg19 information
library("GenomicRanges")
data(hg19Ideogram, package = "biovizBase", envir = environment())
regions <- regsWithP$regions
seqlengths(regions) <- seqlengths(hg19Ideogram)[names(seqlengths(regions))]

## Finally, get the region coverage
regionCov <- getRegionCoverage(fullCov=fullCov, regions=regions)
}
}
\author{
  Andrew Jaffe, Leonardo Collado-Torres
}
\seealso{
  \link{fullCoverage}, \link{calculatePvalues}
}

