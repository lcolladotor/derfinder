\name{generateReport}
\alias{generateReport}
\title{Generate a HTML report exploring the basic results}
\usage{
  generateReport(prefix, outdir = "basicExploration",
    output = "basicExploration.html", project = prefix,
    browse = interactive(), makeBestClusters = TRUE,
    nBestClusters = 2, fullCov = NULL, hg19 = TRUE,
    p.ideos = NULL, txdb = NULL, installMissing = TRUE,
    device = "CairoPNG", fullRegions = NULL,
    fullNullSummary = NULL, optionsStats = NULL)
}
\arguments{
  \item{prefix}{The main data directory path where
  \link{mergeResults} was run. It should be the same as
  \code{mergeResults(prefix)}.}

  \item{outdir}{The name of output directory relative to
  \code{prefix}.}

  \item{output}{The name output HTML file.}

  \item{project}{The title of the project.}

  \item{browse}{If \code{TRUE} the HTML report is opened in
  your browser once it's completed.}

  \item{makeBestClusters}{If \code{TRUE},
  \link{plotCluster} is used on the \code{nBestClusters}
  regions by area. These plots take some time to make.}

  \item{nBestClusters}{The number of best by area plots to
  make.}

  \item{fullCov}{A list where each element is the result
  from \link{loadCoverage} used with \code{cutoff=NULL}.
  The elements of the list should be named according to the
  chromosome number. Can be generated using
  \link{fullCoverage}.}

  \item{hg19}{If \code{TRUE} then the reference is assumed
  to be hg19 and chromosome lengths as well as the default
  transcription database
  (TxDb.Hsapiens.UCSC.hg19.knownGene) will be used.}

  \item{p.ideos}{A list where each element is the result of
  \link[ggbio]{plotIdeogram}. If it's \code{NULL} and
  \code{hg19=TRUE} then they are created for the hg19 human
  reference.}

  \item{txdb}{Specify the transcription database to use for
  making the plots for the top regions by area. If
  \code{NULL} and \code{hg19=TRUE} then
  TxDb.Hsapiens.UCSC.hg19.knownGene is used.}

  \item{installMissing}{If \code{TRUE} all missing required
  packages are installed. Note that some are development
  versions hosted in GitHub.}

  \item{device}{The graphical device used when knitting.
  See more at http://yihui.name/knitr/options (dev
  argument).}

  \item{fullRegions}{Part of the output of
  \link{mergeResults}. Specify it only if you have already
  loaded it in memory.}

  \item{fullNullSummary}{Part of the output of
  \link{mergeResults}. Specify it only if you have already
  loaded it in memory.}

  \item{optionsStats}{Part of the output of
  \link{analyzeChr}. Specify it only if you have already
  loaded it in memory.}
}
\value{
  An HTML report with a basic exploration of the results.
}
\description{
  This function generates a HTML report exploring the basic
  results.
}
\examples{
\dontrun{
generateReport(prefix="run1", makeBestClusters=FALSE)
}
}
\author{
  Leonardo Collado-Torres
}
\seealso{
  \link{mergeResults}, \link{analyzeChr},
  \link{fullCoverage}
}

