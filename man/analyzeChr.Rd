\name{analyzeChr}
\alias{analyzeChr}
\title{Run the derfinder analysis on a chromosome}
\usage{
  analyzeChr(chrnum, coverageInfo, testvars, adjustvars,
    nonzero = TRUE, center = TRUE, testIntercept = FALSE,
    cutoffPre = 5, colsubset = NULL, scalefac = 32,
    chunksize = NULL, cutoffFstat = 1e-08,
    cutoffType = "theoretical", nPermute = 1,
    seeds = as.integer(gsub("-", "", Sys.Date())) + seq_len(nPermute),
    maxGap = 300L, groupInfo = testvars, subject = "hg19",
    mc.cores = getOption("mc.cores", 2L),
    writeOutput = TRUE, returnOutput = FALSE,
    verbose = TRUE)
}
\arguments{
  \item{chrnum}{Used for naming the output files when
  \code{writeOutput=TRUE} and for
  \link[bumphunter]{annotateNearest}. Use '21' instead of
  'chr21'.}

  \item{coverageInfo}{The output from \link{loadCoverage}.}

  \item{testvars}{This argument is passed to
  \link{makeModels}.}

  \item{adjustvars}{This argument is passed to
  \link{makeModels}.}

  \item{nonzero}{This argument is passed to
  \link{makeModels}.}

  \item{center}{This argument is passed to
  \link{makeModels}.}

  \item{testIntercept}{This argument is passed to
  \link{makeModels}.}

  \item{cutoffPre}{This argument is passed to
  \link{preprocessCoverage} (\code{cutoff}).}

  \item{colsubset}{This argument is passed to
  \link{preprocessCoverage}.}

  \item{scalefac}{This argument is passed to
  \link{preprocessCoverage}.}

  \item{chunksize}{This argument is passed to
  \link{preprocessCoverage}.}

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

  \item{maxGap}{This argument is passed to
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

  \item{verbose}{If \code{TRUE} basic status updates will
  be printed along the way.}
}
\value{
  If \code{returnOutput=TRUE}, a list with seven
  components: \describe{ \item{timeinfo }{ The wallclock
  timing information for each step.} \item{optionsStats }{
  The main options used when running this function.}
  \item{models }{ The models used in the analysis.}
  \item{coveragePrep }{ The output from
  \link{preprocessCoverage}.} \item{fstats}{ The output
  from \link{calculateStats}.} \item{regions}{ The output
  from \link{calculatePvalues}.} \item{annotation}{ The
  output from \link[bumphunter]{annotateNearest}.} } These
  are the same components that are written to Rdata files
  if \code{writeOutput=TRUE}.
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
group <- genomeInfo$pop
adjustvars <- data.frame(genomeInfo$gender)
results <- analyzeChr(chrnum="21", coverageInfo=genomeData, testvars=group, adjustvars=adjustvars, cutoffFstat=1, cutoffType="manual", mc.cores=1, writeOutput=FALSE, returnOutput=TRUE)
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

