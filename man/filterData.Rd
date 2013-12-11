\name{filterData}
\alias{filterData}
\title{Filter the positions of interest}
\usage{
filterData(data, cutoff = NULL, index = NULL, colnames = NULL,
  verbose = TRUE)
}
\arguments{
  \item{data}{Either a list of Rle objects or a DataFrame
  with the coverage information.}

  \item{cutoff}{Per base pair, at least one sample has to
  have coverage strictly greater than \code{cutoff} to be
  included in the result.}

  \item{index}{A logical Rle with the positions of the
  chromosome that passed the cutoff. If \code{NULL} it is
  assumed that this is the first time using
  \link{filterData} and thus no previous index exists.}

  \item{colnames}{Specifies the column names to be used for
  the results DataFrame. If \code{NULL}, no names are
  assigned.}

  \item{verbose}{If \code{TRUE} it will report how many
  rows are remaining out of the original ones.}
}
\value{
A list with two components. \describe{ \item{coverage }{ is
a DataFrame object where each column represents a sample.
The number of rows depends on the number of base pairs that
passed the cutoff and the information stored is the
coverage at that given base.} \item{position }{ is a
logical Rle with the positions of the chromosome that
passed the cutoff.} }
}
\description{
For a group of samples this function reads the coverage
information for a specific chromosome directly from the BAM
files. It then merges them into a DataFrame and removes the
bases that do not pass the cutoff. This is a helper
function for \link{loadCoverage} and
\link{preprocessCoverage}.
}
\details{
If \code{cutoff} is \code{NULL} then the data is grouped
into DataFrame without applying any cutoffs. This can be
useful if you want to use \link{loadCoverage} to build the
coverage DataFrame without applying any cutoffs for other
downstream purposes like plotting the coverage values of a
given region. You can always specify the \code{colsubset}
argument in \link{preprocessCoverage} to filter the data
before calculating the F statistics.
}
\examples{
library("IRanges")
x <- Rle(round(runif(1e4, max=10)))
y <- Rle(round(runif(1e4, max=10)))
z <- Rle(round(runif(1e4, max=10)))
DF <- DataFrame(x, y, z)
filt1 <- filterData(DF, 5)
filt1
filt2 <- filterData(filt1$coverage[, 1:2], 5, index=filt1$position)
filt2
## The number of TRUE values in 'position' is the same as the number of rows as in 'coverage'.
identical(sum(filt2$pos), nrow(filt2$cov))
}
\author{
Leonardo Collado-Torres
}
\seealso{
\link{loadCoverage}, \link{preprocessCoverage}
}

