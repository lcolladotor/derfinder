\name{fstats.apply}
\alias{fstats.apply}
\title{Calculate F-statistics per base by extracting chunks from a DataFrame}
\usage{
  fstats.apply(index = Rle(TRUE, nrow(data)), data, mod,
    mod0, adjustF = 0, lowMemDir = NULL)
}
\arguments{
  \item{index}{An index (logical Rle is the best for saving
  memory) indicating which rows of the DataFrame to use.}

  \item{data}{The DataFrame containing the coverage
  information. Normally stored in
  \code{coveragePrep$coverageProcessed} from
  \link{preprocessCoverage}. Could also be the full data
  from \link{loadCoverage}.}

  \item{mod}{The design matrix for the alternative model.
  Should be m by p where p is the number of covariates
  (normally also including the intercept).}

  \item{mod0}{The design matrix for the null model. Should
  be m by p_0.}

  \item{adjustF}{A single value to adjust that is added in
  the denominator of the F-stat calculation. Useful when
  the Residual Sum of Squares of the alternative model is
  very small.}

  \item{lowMemDir}{The directory where the processed chunks
  are saved when using \link{preprocessCoverage} with a
  specified \code{lowMemDir}.}
}
\value{
  A numeric Rle with the F-statistics per base for the
  chunk in question.
}
\description{
  Extract chunks from a DataFrame and get the F-statistics
  on the rows of \code{data}, comparing the models
  \code{mod} (alternative) and \code{mod0} (null). This is
  a helper function for \link{calculateStats} and
  \link{calculatePvalues}.
}
\examples{
## Create the model matrices
mod <- model.matrix(~ genomeInfo$pop)
mod0 <- model.matrix(~ 0 + rep(1, nrow(genomeInfo)))
## Run the function
fstats.output <- fstats.apply(data=genomeData$coverage, mod=mod, mod0=mod0)
fstats.output
}
\author{
  Jeff Leek, Leonardo Collado-Torres
}
\seealso{
  \link{calculateStats}, \link{calculatePvalues}
}

