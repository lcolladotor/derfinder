\name{fstats.apply}
\alias{fstats.apply}
\title{Calculate F-statistics per base by extracting chunks from a DataFrame.}
\usage{
  fstats.apply(i, data, comparison, chunksize, lastloop,
    numrow, scalefac, mod, mod0)
}
\arguments{
  \item{i}{The chunk number identifier.}

  \item{data}{The DataFrame containing the coverage
  information. Normally stored in
  \code{coverageInfo$coverage} from \link{makeCoverage}.}

  \item{comparison}{Whether you are comparing if there is
  \code{expression} (intercept vs no intercept) or
  \code{group differences} (model with group vs no group).}

  \item{chunksize}{How many rows of \code{data} should be
  processed at a time?}

  \item{lastloop}{The last chunk number.}

  \item{numrow}{Total number of rows in \code{data}.}

  \item{scalefac}{A log transformation is used on the count
  tables, so zero counts present a problem.  What number
  should we add to the entire matrix before running the
  models?}

  \item{mod}{The design matrix for the alternative model.
  Should be m by p where p is the number of covariates
  (normally also including the intercept).}

  \item{mod0}{The deisgn matrix for the null model. Should
  be m by p_0.}
}
\value{
  A Rle with the F-statistics per base for the chunk in
  question.
}
\description{
  Extract chunks from a DataFrame, apply the scaling
  factor, log2 transform and then get the F-statistics.
  This is a helper function for \link{calculateStats}.
}
\examples{
## Create the model matrices
mod <- model.matrix(~ brainInfo$outcome)
mod0 <- model.matrix(~ 0 + rep(1, nrow(brainInfo)))
## Run the function
fstats.output <- fstats.apply(1, brainData$coverage, 1000, 5, nrow(brainData$coverage), 32, mod, mod0)
fstats.output
}
\author{
  Leonardo Collado-Torres
}
\seealso{
  \link{calculateStats}, \link{fstats}.
}

