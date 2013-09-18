\name{preprocessCoverage}
\alias{preprocessCoverage}
\title{Transform and split the data}
\usage{
  preprocessCoverage(coverageInfo, groupInfo = NULL,
    cutoff = 5, scalefac = 32, chunksize = 5e+06,
    colsubset = NULL, mc.cores = getOption("mc.cores", 1L),
    verbose = FALSE)
}
\arguments{
  \item{coverageInfo}{A list containing a DataFrame
  --\code{$coverage}-- with the coverage data and a logical
  Rle --\code{$position}-- with the positions that passed
  the cutoff. This object is generated using
  \link{loadCoverage}.}

  \item{groupInfo}{A factor specifying the group membership
  of each sample. If \code{NULL} no group mean coverages
  are calculated. If the factor has more than one level,
  the first one will be used to calculate the log2 fold
  change in \link{calculatePvalues}.}

  \item{cutoff}{This argument is passed to
  \link{filterData}.}

  \item{colsubset}{Optional vector of column indices of
  \code{coverageInfo$coverage} that denote samples you wish
  to include in analysis.}

  \item{scalefac}{A log transformation is used on the count
  tables, so zero counts present a problem.  What number
  should we add to the entire matrix?}

  \item{chunksize}{How many rows of
  \code{coverageInfo$coverage} should be processed at a
  time?}

  \item{verbose}{If \code{TRUE} basic status updates will
  be printed along the way.}

  \item{mc.cores}{This argument is passed to
  \link[parallel]{mclapply} to run \link{fstats.apply}.}
}
\value{
  A list with four components. \describe{
  \item{coverageSplit }{ is a list of DataFrames. Each
  DataFrame column represents a sample and the data is
  partioned according to \code{chunksize}. The coverage
  information is scaled and log2 transformed. Note that if
  \code{colsubset} is not \code{NULL} the number of columns
  will be less than those in \code{coverageInfo$coverage}.
  The total number of rows depends on the number of base
  pairs that passed the \code{cutoff} and the information
  stored is the coverage at that given base. Further note
  that \link{filterData} is re-applied if \code{colsubset}
  is not \code{NULL} and could thus lead to fewer rows
  compared to \code{coverageInfo$coverage}. }
  \item{position }{ is a logical Rle with the positions of
  the chromosome that passed the cutoff.}
  \item{meanCoverage }{ is a numeric Rle with the mean
  coverage at each filtered base.} \item{groupMeans }{ is a
  list of Rle objects containing the mean coverage at each
  filtered base calculated by group. This list has length 0
  if \code{groupInfo=NULL}.} }
}
\description{
  This function takes the coverage data from
  \link{loadCoverage}, scales the data, does the log2
  transformation, and splits it into appropriate chunks for
  using \link{calculateStats}.
}
\details{
  If \code{chunksize} is \code{NULL}, then \code{mc.cores}
  is used to determine the \code{chunksize}. This is useful
  if you want to split the data so each core gets the same
  amount of data (up to rounding).
}
\examples{
## Split the data and transform appropriately before using calculateStats()
dataReady <- preprocessCoverage(genomeData, cutoff=0, scalefac=32, chunksize=1e3, colsubset=NULL, verbose=TRUE)
names(dataReady)
dataReady
}
\author{
  Leonardo Collado-Torres
}
\seealso{
  \link{filterData}, \link{loadCoverage},
  \link{calculateStats}
}

