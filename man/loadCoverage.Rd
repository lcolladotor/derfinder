\name{loadCoverage}
\alias{loadCoverage}
\title{Load the coverage information from a group of BAM files}
\usage{
  loadCoverage(chr, datadir = NULL, sampledirs = NULL,
    samplepatt = NULL, cutoff = 5, chrlen = NULL,
    bamterm = "accepted_hits.bam", output = NULL,
    verbose = TRUE)
}
\arguments{
  \item{chr}{Chromosome to read. Should be in simple
  format. For example, use X and not chrX.}

  \item{datadir}{The main directory where each of the
  \code{sampledirs} is a sub-directory of \code{datadir}.}

  \item{sampledirs}{A character vector with the names of
  the sample directories. If \code{datadir} is \code{NULL}
  it is then assumed that \code{sampledirs} specifies the
  full path to each sample.}

  \item{samplepatt}{If specified and \code{sampledirs} is
  set to \code{NULL}, then the directories matching this
  pattern in \code{datadir} (set to \code{.} if it's set to
  \code{NULL}) are used as the sample directories.}

  \item{cutoff}{Per base pair, at least one sample has to
  have coverage greater than \code{cutoff} to be included
  in the result.}

  \item{chrlen}{The chromosome length in base pairs. If
  it's \code{NULL}, the chromosome length is extracted from
  the BAM files.}

  \item{bamterm}{Name of the BAM file used in each sample.
  By default it is set to \code{accepted_hits.bam} since
  that is the automatic name generated when aligning with
  TopHat. If \code{NULL} it is then ignored when reading
  the BAM files. This can be useful if all the BAM files
  are stored in a single directory.}

  \item{output}{If \code{NULL} then no output is saved in
  disk. If \code{auto} then an automatic name is
  constructed (chrXDF.Rdata for example). If another
  character is specified, then that name is used for the
  output file.}

  \item{verbose}{If \code{TRUE} basic status updates will
  be printed along the way.}
}
\value{
  A list with two components. \describe{ \item{coverage }{
  is a DataFrame object where each column represents a
  sample. The number of rows depends on the number of base
  pairs that passed the cutoff and the information stored
  is the coverage at that given base.} \item{position }{ is
  a logical Rle with the positions of the chromosome that
  passed the cutoff.} }
}
\description{
  For a group of samples this function reads the coverage
  information for a specific chromosome directly from the
  BAM files. It then merges them into a DataFrame and
  removes the bases that do not pass the cutoff.
}
\examples{
datadir <- system.file("extdata", "genomeData", package="derfinder2")
## Reading the data and filtering it is quite fast.
system.time(data <- loadCoverage(chr="21", datadir=datadir, samplepatt="*accepted_hits.bam$", cutoff=0, bamterm=NULL))

## Shorten the column names
colnames(data$coverage) <- gsub("_accepted_hits.bam", "", colnames(data$coverage))
data

## For other examples this data is included in the package
identical(data, genomeData)

## The data is compact enough to be loaded in memory
print(object.size(data), units="Kb")
}
\author{
  Leonardo Collado-Torres, Andrew Jaffe
}

