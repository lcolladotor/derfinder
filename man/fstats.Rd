\name{fstats}
\alias{fstats}
\title{Calculate F-statistics per base}
\usage{
  fstats(dat, mod, mod0)
}
\arguments{
  \item{dat}{A numeric matrix of size n by m. There should
  be n bases of the genome and m samples.}

  \item{mod}{The design matrix for the alternative model.
  Should be m by p where p is the number of covariates
  (normally also including the intercept).}

  \item{mod0}{The deisgn matrix for the null model. Should
  be m by p_0.}
}
\value{
  A vector with the F-statistics per base (length n).
}
\description{
  A function for calculating F-statistics on the rows of
  dat, comparing the models mod (alternative) and mod0
  (null). This is a helper function for
  \link{calculateStats}.
}
\examples{
## Generate random data
set.seed(20130722)
dat <- matrix(rnorm(1000, 5), ncol=10)

## Define the model matrices
group <- sample(letters[1:2], 10, TRUE)
mod <- model.matrix(~ group)
mod0 <- model.matrix(~ 0 + rep(1, 10))

## Get the F-statistics
result <- fstats(dat, mod, mod0)
hist(result)
}
\author{
  Jeff Leek, Leonardo Collado-Torres
}
\seealso{
  \link{calculateStats}
}

