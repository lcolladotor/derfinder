\name{fullCoverage}
\alias{fullCoverage}
\title{Load the unfiltered coverage information from a group of BAM files and a list of chromosomes}
\usage{
  fullCoverage(dirs, chrnums, chrlens = NULL,
    outputs = NULL, mc.cores = getOption("mc.cores", 2L),
    verbose = TRUE)
}
\arguments{
  \item{dirs}{A character vector with the full path to the
  sample BAM files. The names are used for the column names
  of the DataFrame. Check \link{makeBamList} for
  constructing \code{dirs}.}

  \item{chrnums}{The chromosome numbers of the files to
  read.}

  \item{chrlens}{The chromosome lengths in base pairs. If
  it's \code{NULL}, the chromosome length is extracted from
  the BAM files. Otherwise, it should have the same length
  as \code{chrnums}.}

  \item{outputs}{This argument is passed to the
  \code{output} argument of \link{loadCoverage}. If
  \code{NULL} or \code{"auto"} it is then recycled.}

  \item{mc.cores}{This argument is passed to
  \link[parallel]{mclapply}. You should use at most one
  core per chromosome.}

  \item{verbose}{If \code{TRUE} basic status updates will
  be printed along the way.}
}
\value{
  A list with one element per chromosome. \describe{ Each
  element is a DataFrame with the coverage information
  produced by \link{loadCoverage}. }
}
\description{
  For a group of samples this function reads the coverage
  information for several chromosomes directly from the BAM
  files. Per chromosome, it merges the unfiltered coverage
  by sample into a DataFrame. The end result is a list with
  one such DataFrame objects per chromosome.
}
\examples{
datadir <- system.file("extdata", "genomeData", package="derfinder")
dirs <- makeBamList(datadir=datadir, samplepatt="*accepted_hits.bam$", bamterm=NULL)
## Shorten the column names
names(dirs) <- gsub("_accepted_hits.bam", "", names(dirs))

## Reading the data and filtering it is quite fast.
system.time(fullCov <- fullCoverage(dirs=dirs, chrnums=c("21", "22"), mc.cores=2))
fullCov

## You can then use filterData to filter the data if you want to.
## Use mclapply if you want to do so with multiple cores.
library("parallel")
mclapply(fullCov, filterData, cutoff=0, mc.cores=2)
}
\author{
  Leonardo Collado-Torres
}
\seealso{
  \link{loadCoverage}, \link{filterData},
  \link{generateReport}
}

