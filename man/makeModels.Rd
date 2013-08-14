\name{makeModels}
\alias{makeModels}
\title{Build model matrices for differential expression}
\usage{
  makeModels(coverageInfo, group, adjustvars = NULL,
    nonzero = FALSE, verbose = FALSE)
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

  \item{adjustvars}{Optional matrix of adjustment variables
  (e.g. measured confounders, output from SVA, etc.) to use
  in fitting linear models to each nucleotide. These
  variables have to be specified by sample and the number
  of rows must match the number of columns used. It will
  also work if it is a vector of the correct length.}

  \item{nonzero}{If \code{TRUE}, use the median of only the
  nonzero counts as the library size adjustment.}

  \item{verbose}{If \code{TRUE} basic status updates will
  be printed along the way.}
}
\value{
  A list with two components. \describe{ \item{mod }{ The
  alternative model matrix.} \item{mod0 }{ The null model
  matrix.} }
}
\description{
  Builds the model matrices for testing for differential
  expression by comparing a model with a grouping factor
  versus one without it. It adjusts for the confounders
  specified and the median coverage of each sample. The
  resulting models can be used in \link{calculateStats}.
}
\examples{
## Choose the adjusting variables and define all the parameters for makeModels()
coverageInfo <- genomeData
group <- genomeInfo$pop
colsubset <- NULL
adjustvars <- data.frame(genomeInfo$gender)
nonzero <- TRUE
verbose <- TRUE

## Run the function
models <- makeModels(coverageInfo, group, adjustvars, nonzero, verbose)
names(models)
models
}
\author{
  Leonardo Collado-Torres
}
\seealso{
  \link{calculateStats}
}

