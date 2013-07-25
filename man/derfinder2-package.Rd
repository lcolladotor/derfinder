\docType{package}
\name{derfinder2-package}
\alias{derfinder2-package}
\title{Fast differential expression analysis of RNA-seq data at base-pair resolution}
\description{
  Fast differential expression analysis of RNA-seq data at
  base-pair resolution from multiple samples. The analysis
  pipeline involves loading the sample BAM files using
  \link{loadCoverage}, calculating the F-statistics
  according (while adjusting for some confounders) using
  \link{calculateStats}, finding the regions of interest
  using \link{findRegions}, and finally annotating them
  using \link[bumphunter]{annotateNearest} from the
  bumphunter package.
}
\author{
  Leonardo Collado-Torres <lcollado@jhsph.edu>
}
\references{
  Paper coming soon.
}
\keyword{package}

