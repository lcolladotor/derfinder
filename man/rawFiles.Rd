% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rawFiles.R
\name{rawFiles}
\alias{rawFiles}
\title{Construct full paths to a group of raw input files}
\usage{
rawFiles(
  datadir = NULL,
  sampledirs = NULL,
  samplepatt = NULL,
  fileterm = "accepted_hits.bam"
)
}
\arguments{
\item{datadir}{The main directory where each of the \code{sampledirs} is a
sub-directory of \code{datadir}.}

\item{sampledirs}{A character vector with the names of the sample
directories. If \code{datadir} is \code{NULL} it is then assumed that
\code{sampledirs} specifies the full path to each sample.}

\item{samplepatt}{If specified and \code{sampledirs} is set to \code{NULL},
then the directories matching this pattern in \code{datadir} (set to
\code{.} if it's set to \code{NULL}) are used as the sample directories.}

\item{fileterm}{Name of the BAM or BigWig file used in each sample. By
default it is set to \code{accepted_hits.bam} since that is the automatic
name generated when aligning with TopHat. If \code{NULL} it is then ignored
when reading the rawfiles. This can be useful if all the raw files are
stored in a single directory.}
}
\value{
A vector with the full paths to the raw files and sample names
stored as the vector names.
}
\description{
For a group of samples this function creates the list of paths to the raw
input files which can then be used in \link{loadCoverage}. The raw input
files are either BAM files or BigWig files.
}
\details{
This function can also be used to identify a set of BigWig files.
}
\examples{
## Get list of BAM files included in derfinder
datadir <- system.file("extdata", "genomeData", package = "derfinder")
files <- rawFiles(
    datadir = datadir, samplepatt = "*accepted_hits.bam$",
    fileterm = NULL
)
files
}
\seealso{
\link{loadCoverage}
}
\author{
Leonardo Collado-Torres
}
