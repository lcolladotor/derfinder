\name{calculateStats}
\alias{calculateStats}
\title{Calculate F-statistics at base pair resolution from a loaded BAM files}
\usage{
  calculateStats(coverageInfo, group,
    comparison = "group differences", colsubset = NULL,
    adjustvars = NULL, cutoff = 5, scalefac = 32,
    nonzero = FALSE, chunksize = 5e+06,
    mc.cores = getOption("mc.cores", 2L), verbose = TRUE)
}
\arguments{
  \item{coverageInfo}{A list containing a DataFrame
  --\code{$coverage}-- with the coverage data and a logical
  Rle --\code{$position}-- with the positions that passed
  the cutoff. This object is generated using
  \link{loadCoverage}.}

  \item{group}{A factor vector specifying the sample
  groups. It's length should match the number of columns
  used from \code{coverageInfo$coverage}.}

  \item{comparison}{Whether you are comparing if there is
  \code{expression} (intercept vs no intercept) or
  \code{group differences} (model with group vs no group).}

  \item{colsubset}{Optional vector of column indices of
  \code{coverageInfo$coverage} that denote samples you wish
  to include in analysis.}

  \item{adjustvars}{Optional matrix of adjustment variables
  (e.g. measured confounders, output from SVA, etc.) to use
  in fitting linear models to each nucleotide. These
  variables have to be specified by sample and the number
  of rows must match the number of columns used.}

  \item{cutoff}{Per base pair, at least one sample has to
  have coverage greater than \code{cutoff} to be included
  in the result.}

  \item{scalefac}{A log transformation is used on the count
  tables, so zero counts present a problem.  What number
  should we add to the entire matrix before running the
  models?}

  \item{nonzero}{If \code{TRUE}, use the median of only the
  nonzero counts as the library size adjustment.}

  \item{chunksize}{How many rows of
  \code{coverageInfo$coverage} should be processed at a
  time?}

  \item{mc.cores}{This argument is passed to
  \link[parallel]{mclapply}) to run \link{fstats.apply}.}

  \item{verbose}{If \code{TRUE} basic status updates will
  be printed along the way.}
}
\value{
  A list with five components. \describe{ \item{coverage }{
  is a DataFrame object where each column represents a
  sample. The coverage information is scaled and log2
  transformed. Note that if \code{colsubset} is not
  \code{NULL} the number of columns will be less than those
  in \code{coverageInfo$coverage}. The number of rows
  depends on the number of base pairs that passed the
  cutoff and the information stored is the coverage at that
  given base. Further note that \link{filterData} is
  re-applied if \code{colsubset} is not \code{NULL} and
  could thus lead to fewer rows compared to
  \code{coverageInfo$coverage}. } \item{position }{ is a
  logical Rle with the positions of the chromosome that
  passed the cutoff.} \item{fstats }{ is a numeric Rle with
  the F-statistics per base pair that passed the cutoff.}
  \item{mod }{ The alternative model matrix.} \item{mod0 }{
  The null model matrix.} }
}
\description{
  After defining the models of interest and adjusting for
  confounders, this function extracts the data from loaded
  BAM files (\link{loadCoverage}) and calculates the
  F-statistics.
}
\examples{
## Choose the adjusting variables and define all the parameters for calculateStats()
coverageInfo <- brainData
group <- brainInfo$outcome
colsubset <- NULL
adjustvars <- brainInfo[, c("sex", "age", "left.hemisph", "pmi", "brainpH")]
cutoff <- 5
scalefac <- 32
nonzero <- TRUE
chunksize <- 1e+03
mc.cores <- 1
verbose <- TRUE
## Run the function
stats <- calculateStats(coverageInfo, group, comparison="group differences", colsubset, adjustvars, cutoff, scalefac, nonzero, chunksize, mc.cores, verbose)
names(stats)
stats
}
\author{
  Leonardo Collado-Torres
}

