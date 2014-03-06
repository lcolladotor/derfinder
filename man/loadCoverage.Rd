\name{loadCoverage}
\alias{loadCoverage}
\title{Load the coverage information from a group of BAM files}
\usage{
loadCoverage(dirs, chr, cutoff = NULL, bai = NULL, chrlen = NULL,
  output = NULL, verbose = TRUE)
}
\arguments{
  \item{dirs}{A character vector with the full path to the
  sample BAM files. The names are used for the column names
  of the DataFrame. Check \link{makeBamList} for
  constructing \code{dirs}. \code{dirs} can also be a
  \code{BamFileList} object created with
  \link[Rsamtools]{BamFileList}.}

  \item{chr}{Chromosome to read. Should be in simple
  format. For example, use X and not chrX.}

  \item{cutoff}{This argument is passed to
  \link{filterData}.}

  \item{bai}{The full path to the BAM index files. If
  \code{NULL} it is assumed that the BAM index files are in
  the same location as the BAM files and that they have the
  .bai extension. Ignored if \code{dirs} is a
  \code{BamFileList} object.}

  \item{chrlen}{The chromosome length in base pairs. If
  it's \code{NULL}, the chromosome length is extracted from
  the BAM files.}

  \item{output}{If \code{NULL} then no output is saved in
  disk. If \code{auto} then an automatic name is
  constructed (chrXCovInfo.Rdata for example). If another
  character is specified, then that name is used for the
  output file.}

  \item{verbose}{If \code{TRUE} basic status updates will
  be printed along the way.}
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
bases that do not pass the cutoff.
}
\examples{
datadir <- system.file("extdata", "genomeData", package="derfinder")
dirs <- makeBamList(datadir=datadir, samplepatt="*accepted_hits.bam$", bamterm=NULL)
## Shorten the column names
names(dirs) <- gsub("_accepted_hits.bam", "", names(dirs))

## Reading the data and filtering it is quite fast.
system.time(data <- loadCoverage(dirs=dirs, chr="21", cutoff=0))
data

## For other examples this data is included in the package
identical(data, genomeData)

## The data is compact enough to be loaded in memory
print(object.size(data), units="Kb")

\dontrun{
## Read the coverage without applying any cutoff.
## This can be useful for downstream analysis including coverage plots.
system.time(data2 <- loadCoverage(dirs=dirs, chr="21", cutoff=NULL))

## Compare vs raw data provided in the package
identical(data2, genomeDataRaw)

## Note that the object size is pretty much the same due to the Rle compression
print(object.size(data2), units="Kb")
}
}
\author{
Leonardo Collado-Torres, Andrew Jaffe
}

