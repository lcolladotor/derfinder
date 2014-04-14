\name{mergeResults}
\alias{mergeResults}
\title{Merge results from different chromosomes}
\usage{
mergeResults(chrnums = c(1:22, "X", "Y"), prefix = ".",
  significantCut = c(0.05, 0.1), genomicState, minoverlap = 20,
  fullOrCoding = "full", mergePrep = FALSE, verbose = TRUE)
}
\arguments{
  \item{chrnums}{The chromosome numbers of the files to be
  merged.}

  \item{prefix}{The main data directory path, which can be
  useful if \link{analyzeChr} is used for several
  parameters and the results are saved in different
  directories.}

  \item{significantCut}{A vector of length two specifiying
  the cutoffs used to determine significance. The first
  element is used to determine significance for the
  p-values and the second element is used for the q-values
  just like in \link{calculatePvalues}.}

  \item{genomicState}{This argument is passed to
  \link{annotateRegions}.}

  \item{minoverlap}{This argument is passed to
  \link{annotateRegions}.}

  \item{fullOrCoding}{This argument is passed to
  \link{annotateRegions}.}

  \item{mergePrep}{If \code{TRUE} the output from
  \link{preprocessCoverage} is merged.}

  \item{verbose}{If \code{TRUE} basic status updates will
  be printed along the way.}
}
\value{
Seven Rdata files. \describe{ \item{fullFstats.Rdata }{
Full F-statistics from all chromosomes in a list of Rle
objects.} \item{fullTime.Rdata }{ Timing information from
all chromosomes.} \item{fullNullSummary.Rdata}{ A DataFrame
with the null region information: statistic, width,
chromosome and permutation identifier. It's ordered by the
statistics} \item{fullRegions.Rdata}{ GRanges object with
regions found and with full annotation from
\link[bumphunter]{annotateNearest}. Note that the column
\code{strand} from \link[bumphunter]{annotateNearest} is
renamed to \code{annoStrand} to comply with GRanges
specifications. } \item{fullCoveragePrep.Rdata}{ A list
with the pre-processed coverage data from all chromosomes.}
\item{fullAnnotatedRegions.Rdata}{ A list as constructed in
\link{annotateRegions} with the assigned genomic states.}
\item{optionsMerge.Rdata}{ A list with the options used
when merging the results. Used in
\link[derfinderReport]{generateReport}.} }
}
\description{
This function merges the results from running
\link{analyzeChr} on several chromosomes and assigns
genomic states using \link{annotateRegions}. It
re-calculates the p-values and q-values using the pooled
areas from the null regions from all chromosomes. Once the
results have been merged,
\link[derfinderReport]{generateReport} can be used to
generate an HTML report of the results. The
\code{derfinderReport} package is available at
https://github.com/lcolladotor/derfinderReport.
}
\examples{
\dontrun{
## Merge the results
mergeResults(prefix="run1", genomicState=genomicState)

## You can then explore the wallclock time spent on each step
prefix <- "." ## Set to the same prefix you used in mergeResults()
load(file.path(prefix, "fullRegions.Rdata"))

## Process the time info
time <- lapply(fullTime, function(x) data.frame(diff(x)))
time <- do.call(rbind, time)
colnames(time) <- "sec"
time$sec <- as.integer(round(time$sec))
time$min <- time$sec / 60
time$chr <- paste0("chr", gsub("\\\\..*", "", rownames(time)))
time$step <- gsub(".*\\\\.", "", rownames(time))
rownames(time) <- seq_len(nrow(time))

## Make plot
library("ggplot2")
ggplot(time, aes(x=step, y=min, colour=chr)) + geom_point() + labs(title="Wallclock time by step") + scale_colour_discrete(limits=chrs) + scale_x_discrete(limits=names(fullTime[[1]])[-1]) + ylab("Time (min)") + xlab("Step")
}
}
\author{
Leonardo Collado-Torres
}
\seealso{
\link{analyzeChr}, \link{calculatePvalues},
\link{annotateRegions}
}

