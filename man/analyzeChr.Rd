\name{analyzeChr}
\alias{analyzeChr}
\title{Run the derfinder analysis on a chromosome}
\usage{
analyzeChr(chrnum, coverageInfo, models, cutoffPre = 5, colsubset = NULL,
  scalefac = 32, chunksize = NULL, adjustF = 0, cutoffFstat = 1e-08,
  cutoffType = "theoretical", nPermute = 1, seeds = as.integer(gsub("-",
  "", Sys.Date())) + seq_len(nPermute), maxRegionGap = 0L,
  maxClusterGap = 300L, groupInfo, subject = "hg19",
  mc.cores = getOption("mc.cores", 2L), writeOutput = TRUE,
  returnOutput = FALSE, runAnnotation = TRUE, lowMemDir = NULL,
  verbose = TRUE)
}
\arguments{
  \item{chrnum}{Used for naming the output files when
  \code{writeOutput=TRUE} and for
  \link[bumphunter]{annotateNearest}. Use '21' instead of
  'chr21'.}

  \item{coverageInfo}{The output from \link{loadCoverage}.}

  \item{models}{The output from \link{makeModels}.}

  \item{cutoffPre}{This argument is passed to
  \link{preprocessCoverage} (\code{cutoff}).}

  \item{colsubset}{This argument is passed to
  \link{preprocessCoverage}.}

  \item{scalefac}{This argument is passed to
  \link{preprocessCoverage}.}

  \item{chunksize}{This argument is passed to
  \link{preprocessCoverage}.}

  \item{adjustF}{A single value to adjust that is added in
  the denominator of the F-stat calculation. Useful when
  the Residual Sum of Squares of the alternative model is
  very small.}

  \item{cutoffFstat}{This is used to determine the cutoff
  argument of \link{calculatePvalues} and it's behaviour is
  determined by \code{cutoffType}.}

  \item{cutoffType}{If set to \code{empirical}, the
  \code{cutoffFstat} (example: 0.99) quantile is used via
  \link{quantile}. If set to \code{theoretical}, the
  theoretical \code{cutoffFstats} (example: 1e-08) is
  calculated via \link{qf}. If set to \code{manual},
  \code{cutoffFstats} is passed to \link{calculatePvalues}
  without any other calculation.}

  \item{nPermute}{This argument is passed to
  \link{calculatePvalues}.}

  \item{seeds}{This argument is passed to
  \link{calculatePvalues}.}

  \item{maxRegionGap}{This argument is passed to
  \link{calculatePvalues}.}

  \item{maxClusterGap}{This argument is passed to
  \link{calculatePvalues}.}

  \item{groupInfo}{A factor specifying the group membership
  of each sample that can later be used with
  \code{plotRegion}.}

  \item{subject}{This argument is passed to
  \link[bumphunter]{annotateNearest}. Note that only
  \code{hg19} works right now.}

  \item{mc.cores}{This argument is passed to
  \link{preprocessCoverage} (useful if
  \code{chunksize=NULL}), \link{calculateStats} and
  \link{calculatePvalues}.}

  \item{writeOutput}{If \code{TRUE}, output Rdata files are
  created at each step inside a directory with the
  chromosome name (example: 'chr21' if \code{chrnum="21"}).
  One Rdata files is created for each component described
  in the return section.}

  \item{returnOutput}{If \code{TRUE}, it returns a list
  with the results from each step. Otherwise, it returns
  \code{NULL}.}

  \item{runAnnotation}{If \code{TRUE}
  \link[bumphunter]{annotateNearest} is run. Otherwise this
  step is skipped.}

  \item{lowMemDir}{The directory where the processed chunks
  are saved when using \link{preprocessCoverage} with a
  specified \code{lowMemDir}.}

  \item{verbose}{If \code{TRUE} basic status updates will
  be printed along the way.}
}
\value{
If \code{returnOutput=TRUE}, a list with six components:
\describe{ \item{timeinfo }{ The wallclock timing
information for each step.} \item{optionsStats }{ The main
options used when running this function.}
\item{coveragePrep }{ The output from
\link{preprocessCoverage}.} \item{fstats}{ The output from
\link{calculateStats}.} \item{regions}{ The output from
\link{calculatePvalues}.} \item{annotation}{ The output
from \link[bumphunter]{annotateNearest}.} } These are the
same components that are written to Rdata files if
\code{writeOutput=TRUE}.
}
\description{
This is a major wrapper for running several key functions
from this package. It is meant to be used after
\link{loadCoverage} has been used for a specific
chromosome. The steps run include \link{makeModels},
\link{preprocessCoverage}, \link{calculateStats},
\link{calculatePvalues} and
\link[bumphunter]{annotateNearest}.
}
\examples{
## Collapse the coverage information
collapsedFull <- collapseFullCoverage(list(genomeData$coverage), verbose=TRUE)

## Calculate library size adjustments
sampleDepths <- sampleDepth(collapsedFull, probs=c(0.5), nonzero=TRUE, verbose=TRUE)

## Build the models
group <- genomeInfo$pop
adjustvars <- data.frame(genomeInfo$gender)
models <- makeModels(sampleDepths, testvars=group, adjustvars=adjustvars)

## Analyze the chromosome
results <- analyzeChr(chrnum="21", coverageInfo=genomeData, models=models, cutoffFstat=1, cutoffType="manual", groupInfo=group, mc.cores=1, writeOutput=FALSE, returnOutput=TRUE)
names(results)
}
\author{
Leonardo Collado-Torres
}
\seealso{
\link{makeModels}, \link{preprocessCoverage},
\link{calculateStats}, \link{calculatePvalues},
\link[bumphunter]{annotateNearest}
}

