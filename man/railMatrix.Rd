% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/railMatrix.R
\name{railMatrix}
\alias{railMatrix}
\title{Identify regions data by a coverage filter and get a count matrix from
BigWig files}
\usage{
railMatrix(
  chrs,
  summaryFiles,
  sampleFiles,
  L,
  cutoff = NULL,
  maxClusterGap = 300L,
  totalMapped = NULL,
  targetSize = 4e+07,
  file.cores = 1L,
  chrlens = NULL,
  ...
)
}
\arguments{
\item{chrs}{A character vector with the names of the chromosomes.}

\item{summaryFiles}{A character vector (or BigWigFileList) with the paths to
the summary BigWig files created by Rail. Either mean or median files. These
are library size adjusted by default in Rail. The order of the files in this
vector must match the order in \code{chrs}.}

\item{sampleFiles}{A character vector with the paths to the BigWig files
by sample. These files are unadjusted for library size.}

\item{L}{The width of the reads used. Either a vector of length 1 or length
equal to the number of samples.}

\item{cutoff}{The base-pair level cutoff to use. It's behavior is controlled
by \code{filter}.}

\item{maxClusterGap}{This determines the maximum gap between candidate ERs.}

\item{totalMapped}{A vector with the total number of reads mapped for each
sample. The vector should be in the same order as the samples in \code{data}.
Providing this data adjusts the coverage to reads in \code{targetSize}
library prior to filtering. See \link{getTotalMapped} for
calculating this vector.}

\item{targetSize}{The target library size to adjust the coverage to. Used
only when \code{totalMapped} is specified. By default, it adjusts to
libraries with 40 million reads, which matches the default used in Rail.}

\item{file.cores}{Number of cores used for loading the BigWig files per chr.}

\item{chrlens}{The chromosome lengths in base pairs. If it's \code{NULL},
the chromosome length is extracted from the BAM files. Otherwise, it should
have the same length as \code{chrs}.}

\item{...}{Arguments passed to other methods and/or advanced arguments.
Advanced arguments:
\describe{
\item{verbose }{ If \code{TRUE} basic status updates will be printed along
the way. Default: \code{TRUE}.}
\item{verbose.load }{ If \code{TRUE} basic status updates will be printed
along the way when loading data. Default: \code{TRUE}.}
\item{BPPARAM.railChr }{ A BPPARAM object to use for the chr step. Set to
\link[BiocParallel:SerialParam-class]{SerialParam} when \code{file.cores = 1} and
\link[BiocParallel:SnowParam-class]{SnowParam} otherwise.}
\item{chunksize }{ Chunksize to use. Default: 1000.}
}
Passed to \link{filterData}, \link{findRegions} and \link{define_cluster}.}
}
\value{
A list with one entry per chromosome. Then per chromosome, a list
with two components.
\describe{
\item{regions }{ A set of regions based on the coverage filter cutoff as
returned by \link{findRegions}.}
\item{coverageMatrix }{  A matrix with the mean coverage by sample for each
candidate region.}
}
}
\description{
Rail (available at http://rail.bio) generates coverage BigWig files. These
files are faster to load in R and to process. Rail creates an un-adjusted
coverage BigWig file per sample and an adjusted summary coverage BigWig file
by chromosome (median or mean). \link{railMatrix} reads in the mean (or
median) coverage BigWig file and applies a threshold cutoff to identify
expressed regions (ERs). Then it goes back to the sample coverage BigWig
files and extracts the base level coverage for each sample. Finally it
summarizes this information in a matrix with one row per ERs and one column
per sample. This function is similar to \link{regionMatrix} but is faster
thanks to the advantages presented by BigWig files.
}
\details{
Given a set of un-filtered coverage data (see \link{fullCoverage}), create
candidate regions by applying a cutoff on the coverage values,
and obtain a count matrix where the number of rows corresponds to the number
of candidate regions and the number of columns corresponds to the number of
samples. The values are the mean coverage for a given sample for a given
region.

This function uses several other \link{derfinder-package}
functions. Inspect the code if interested.

You should use at most one core per chromosome.
}
\examples{

## BigWig files are not supported in Windows
if (.Platform$OS.type != "windows") {
    ## Get data
    library("derfinderData")

    ## Identify sample files
    sampleFiles <- rawFiles(system.file("extdata", "AMY",
        package =
            "derfinderData"
    ), samplepatt = "bw", fileterm = NULL)
    names(sampleFiles) <- gsub(".bw", "", names(sampleFiles))

    ## Create the mean bigwig file. This file is normally created by Rail
    ## but in this example we'll create it manually.
    library("GenomicRanges")
    fullCov <- fullCoverage(files = sampleFiles, chrs = "chr21")
    meanCov <- Reduce("+", fullCov$chr21) / ncol(fullCov$chr21)
    createBw(list("chr21" = DataFrame("meanChr21" = meanCov)),
        keepGR =
            FALSE
    )

    summaryFile <- "meanChr21.bw"

    ## Get the regions
    regionMat <- railMatrix(
        chrs = "chr21", summaryFiles = summaryFile,
        sampleFiles = sampleFiles, L = 76, cutoff = 5.1,
        maxClusterGap = 3000L
    )

    ## Explore results
    names(regionMat$chr21)
    regionMat$chr21$regions
    dim(regionMat$chr21$coverageMatrix)
}
}
\author{
Leonardo Collado-Torres
}
