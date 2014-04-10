\name{plotRegionCoverage}
\alias{plotRegionCoverage}
\title{Makes plots for every region while summarizing the annotation}
\usage{
plotRegionCoverage(regions, regionCoverage, groupInfo, nearestAnnotation,
  annotatedRegions, whichRegions = seq_len(min(100, length(regions))),
  colors = NULL, scalefac = 32, ask = interactive(), verbose = TRUE,
  ylab = "Coverage")
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
  \link[RColorBrewer]{brewer.pal} with the \code{"Dark2"}
  color scheme is used.}

  \item{scalefac}{The parameter used in
  \link{preprocessCoverage}.}

  \item{ask}{If \code{TRUE} then the user is prompted
  before each plot is made.}

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
\dontrun{
## Create GenomicState object:
## Hsapiens.UCSC.hg19.knownGene GenomicState
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

## Creating this GenomicState object takes around 8 min
GenomicState.Hsapiens.UCSC.hg19.knownGene <- makeGenomicState(txdb=txdb)

## Annotate regions
annotatedRegions <- annotateRegions(regions=genomeRegions$regions, genomicState=GenomicState.Hsapiens.UCSC.hg19.knownGene, minoverlap=1)

## Find nearest annotation
library("bumphunter")
nearestAnnotation <- annotateNearest(genomeRegions$regions, "hg19")

## Assign chr lengths using hg19 information
library("GenomicRanges")
data(hg19Ideogram, package = "biovizBase", envir = environment())
regions <- genomeRegions$regions
seqlengths(regions) <- seqlengths(hg19Ideogram)[names(seqlengths(regions))]

## Obtain fullCov object
datadir <- system.file("extdata", "genomeData", package="derfinder")
dirs <- makeBamList(datadir=datadir, samplepatt="*accepted_hits.bam$", bamterm=NULL)
## Shorten the column names
names(dirs) <- gsub("_accepted_hits.bam", "", names(dirs))

## Reading the data and filtering it is quite fast.
fullCov <- fullCoverage(dirs=dirs, chrnums="21", mc.cores=1)

## Get the region coverage
regionCov <- getRegionCoverage(fullCov=fullCov, regions=regions)

## Make plots for the regions
plotRegionCoverage(regions=regions, regionCoverage=regionCov, groupInfo=genomeInfo$pop, nearestAnnotation=nearestAnnotation, annotatedRegions=annotatedRegions, N=2, ask=TRUE)

## If you prefer, you can save the plots to a pdf file
pdf("ders.pdf", h = 6, w = 9)
plotRegionCoverage(regions=regions, regionCoverage=regionCov, groupInfo=genomeInfo$pop, nearestAnnotation=nearestAnnotation, annotatedRegions=annotatedRegions)
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

