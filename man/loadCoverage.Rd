% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loadCoverage.R
\name{loadCoverage}
\alias{loadCoverage}
\title{Load the coverage information from a group of BAM files}
\usage{
loadCoverage(
  files,
  chr,
  cutoff = NULL,
  filter = "one",
  chrlen = NULL,
  output = NULL,
  bai = NULL,
  ...
)
}
\arguments{
\item{files}{A character vector with the full path to the sample BAM files
(or BigWig files).
The names are used for the column names of the DataFrame. Check
\link{rawFiles} for constructing \code{files}. \code{files} can also be a
\link[Rsamtools:BamFile-class]{BamFileList}, \link[Rsamtools:BamFile-class]{BamFile},
\link[rtracklayer:BigWigFile]{BigWigFileList}, or \link[rtracklayer:BigWigFile]{BigWigFile} object.}

\item{chr}{Chromosome to read. Should be in the format matching the one used
in the raw data.}

\item{cutoff}{This argument is passed to \link{filterData}.}

\item{filter}{Has to be either \code{'one'} (default) or \code{'mean'}. In
the first case, at least one sample has to have coverage above \code{cutoff}.
In the second case, the mean coverage has to be greater than \code{cutoff}.}

\item{chrlen}{The chromosome length in base pairs. If it's \code{NULL}, the
chromosome length is extracted from the BAM files.}

\item{output}{If \code{NULL} then no output is saved in disk. If \code{auto}
then an automatic name is constructed using UCSC names (chrXCovInfo.Rdata
for example). If another character is specified, then that name is used for
the output file.}

\item{bai}{The full path to the BAM index files. If \code{NULL} it is
assumed that the BAM index files are in the same location as the BAM files
and that they have the .bai extension. Ignored if \code{files} is a
\code{BamFileList} object or if \code{inputType=='BigWig'}.}

\item{...}{Arguments passed to other methods and/or advanced arguments.
Advanced arguments:
\describe{
\item{verbose }{ If \code{TRUE} basic status updates will be printed along
the way.}
\item{inputType }{ Has to be either \code{bam} or \code{BigWig}. It specifies
the format of the raw data files. By default it's set to \code{bam} before
trying to guess it from the file names.}
\item{tilewidth }{ When specified, \link[GenomicRanges:tileGenome]{tileGenome} is used
to break up the chromosome into chunks. We don't recommend this for BAM
files as the coverage in the borders of the chunks might be slightly off.}
\item{which }{ \code{NULL} by default. When a \code{GRanges} is specified,
this specific region of the genome is loaded instead of the full chromosome.}
\item{fileStyle }{ The naming style of the chromosomes in the input files.
If the global option \code{chrsStyle} is not set, the naming style won't be
changed. This option is useful when you want to use a specific naming style
but the raw files use another style.}
\item{protectWhich }{ When not \code{NULL} and \code{which} is specified,
this argument specifies by how much the ranges in \code{which} will be
expanded. This helps get the same base level coverage you would get from
reading the coverage for a full chromosome from BAM files. Otherwise some
reads might be excluded and thus the base level coverage will be lower.
\code{NULL} by default.}
\item{drop.D }{ Whether to drop the bases with 'D' in the CIGAR strings
or to include them. Only used with BAM files. \code{FALSE} by default.}
\item{sampleNames }{ Column names to be used the samples. By default it's
\code{names(files)}.}
}
Passed to \link{extendedMapSeqlevels}, \link{define_cluster},
\link[Rsamtools:ScanBamParam-class]{scanBamFlag} and \link{filterData}.
Note that \link{filterData} is used internally
by \link{loadCoverage} and has the important arguments \code{totalMapped}
and \code{targetSize} which can be used to normalize the coverage by
library size. See \link{getTotalMapped} for calculating \code{totalMapped}.}
}
\value{
A list with two components.
\describe{
\item{coverage }{ is a DataFrame object where each column represents a
sample. The number of rows depends on the number of base pairs that passed
the cutoff and the information stored is the coverage at that given base.}
\item{position }{  is a logical Rle with the positions of the chromosome
that passed the cutoff.}
}
}
\description{
For a group of samples this function reads the coverage information for a
specific chromosome directly from the BAM files. It then merges them into a
DataFrame and removes the bases that do not pass the cutoff.
}
\details{
The \code{...} argument can be used to control which alignments to consider
when reading from BAM files. See \link[Rsamtools:ScanBamParam-class]{scanBamFlag}.

Parallelization for loading the data in chunks is used only used when
\code{tilewidth} is specified. You may use up to one core per tile.

If you set the advanced argument \code{drop.D = TRUE}, bases with CIGAR
string "D" (deletion from reference) will be excluded from the base-level
coverage calculation.

If you are working with data from an organism different from 'Homo sapiens'
specify so by setting the global 'species' and 'chrsStyle' options. For
example:
\code{options(species = 'arabidopsis_thaliana')}
\code{options(chrsStyle = 'NCBI')}
}
\examples{
datadir <- system.file("extdata", "genomeData", package = "derfinder")
files <- rawFiles(
    datadir = datadir, samplepatt = "*accepted_hits.bam$",
    fileterm = NULL
)
## Shorten the column names
names(files) <- gsub("_accepted_hits.bam", "", names(files))

## Read and filter the data, only for 2 files
dataSmall <- loadCoverage(files = files[1:2], chr = "21", cutoff = 0)
\dontrun{
## Export to BigWig files
createBw(list("chr21" = dataSmall))

## Load data from BigWig files
dataSmall.bw <- loadCoverage(c(
    ERR009101 = "ERR009101.bw", ERR009102 =
        "ERR009102.bw"
), chr = "chr21")

## Compare them
mapply(function(x, y) {
    x - y
}, dataSmall$coverage, dataSmall.bw$coverage)

## Note that the only difference is the type of Rle (integer vs numeric) used
## to store the data.
}

}
\seealso{
\link{fullCoverage}, \link{getTotalMapped}
}
\author{
Leonardo Collado-Torres, Andrew Jaffe
}
