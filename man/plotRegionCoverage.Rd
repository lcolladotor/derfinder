\name{plotRegionCoverage}
\alias{plotRegionCoverage}
\title{Makes plots for every region while summarizing the annotation}
\usage{
plotRegionCoverage(regions, regionCoverage, groupInfo, nearestAnnotation,
  annotatedRegions, whichRegions = seq_len(min(100, length(regions))),
  colors = NULL, scalefac = 32, ask = interactive(), ylab = "Coverage",
  verbose = TRUE)
}
\arguments{
  \item{regions}{The \code{$regions} output from
  \link{calculatePvalues}.}

  \item{regionCoverage}{The output from
  \link{getRegionCoverage} used on \code{regions}.}

  \item{groupInfo}{A factor specifying the group membership
  of each sample. It will be used to color the samples by
  group.}

  \item{nearestAnnotation}{The output from
  \link[bumphunter]{annotateNearest} used on
  \code{regions}.}

  \item{annotatedRegions}{The output from
  \link{annotateRegions} used on \code{regions}.}

  \item{whichRegions}{An integer vector with the index of
  the regions to plot.}

  \item{colors}{If \code{NULL} then
  \link[RColorBrewer]{brewer.pal} with the \code{'Dark2'}
  color scheme is used.}

  \item{scalefac}{The parameter used in
  \link{preprocessCoverage}.}

  \item{ask}{If \code{TRUE} then the user is prompted
  before each plot is made.}

  \item{ylab}{The name of the of the Y axis.}

  \item{verbose}{If \code{TRUE} basic status updates will
  be printed along the way.}
}
\value{
A plot for every region showing the coverage of each sample
at each base of the region as well as the summarized
annotation information.
}
\description{
This function takes the regions found in
\link{calculatePvalues} and assigns them genomic states
contructed with \link{makeGenomicState}. The main workhorse
functions are \link[IRanges]{countOverlaps} and
\link[IRanges]{findOverlaps}. For an alternative plot check
\link{plotCluster} which is much slower and we recommend
it's use only after quickly checking the results with this
function.
}
\examples{
## Annotate regions, first two regions only
regions <- genomeRegions$regions[1:2]
annotatedRegions <- annotateRegions(regions=regions,
    genomicState=genomicState, minoverlap=1)

## Find nearest annotation
library('bumphunter')
nearestAnnotation <- annotateNearest(regions, 'hg19')

## Obtain fullCov object
fullCov <- list('21'=genomeDataRaw$coverage)

## Assign chr lengths using hg19 information
library('GenomicRanges')
data(hg19Ideogram, package = 'biovizBase', envir = environment())
seqlengths(regions) <- seqlengths(hg19Ideogram)[names(seqlengths(regions))]

## Get the region coverage
regionCov <- getRegionCoverage(fullCov=fullCov, regions=regions)

## Make plots for the regions
plotRegionCoverage(regions=regions, regionCoverage=regionCov,
    groupInfo=genomeInfo$pop, nearestAnnotation=nearestAnnotation,
    annotatedRegions=annotatedRegions, whichRegions=1:2)

\dontrun{
## If you prefer, you can save the plots to a pdf file
pdf('ders.pdf', h = 6, w = 9)
plotRegionCoverage(regions=regions, regionCoverage=regionCov,
    groupInfo=genomeInfo$pop, nearestAnnotation=nearestAnnotation,
    annotatedRegions=annotatedRegions, whichRegions=1:2, ask=FALSE)
dev.off()
}
}
\author{
Andrew Jaffe, Leonardo Collado-Torres
}
\seealso{
\link{calculatePvalues}, \link{getRegionCoverage},
\link[bumphunter]{annotateNearest}, \link{annotateRegions},
\link{plotCluster}
}

