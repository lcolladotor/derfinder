\name{coverageToExon}
\alias{coverageToExon}
\title{Extract coverage information for exons}
\usage{
  coverageToExon(fullCov, genomicState,
    fullOrCoding = "full", L = 100, returnType = "raw",
    verbose = TRUE)
}
\arguments{
  \item{fullCov}{A list where each element is the result
  from \link{loadCoverage} used with \code{cutoff=NULL}.
  The elements of the list should be named according to the
  chromosome number. Can be generated using
  \link{fullCoverage}.}

  \item{genomicState}{The output from
  \link{makeGenomicState}.}

  \item{fullOrCoding}{If \code{full} then the
  \code{genomicState$fullGenome} genomic state information
  is used. If \code{coding}, then the
  \code{genomicState$codingGenome} genomic state
  information is used.}

  \item{L}{The width of the reads used.}

  \item{returnType}{If \code{raw}, then the raw coverage
  information per exon is returned. If \code{rpkm}, RPKM
  values are calculated for each exon.}

  \item{verbose}{If \code{TRUE} basic status updates will
  be printed along the way.}
}
\value{
  A matrix (nrow = number of exons in \code{genomicState}
  corresponding to the chromosomes in \code{fullCov}, ncol
  = number of samples) with the number of reads (or RPKM)
  per exon. The row names correspond to the row indexes of
  \code{genomicState$fullGenome} (if
  \code{fullOrCoding="full"}) or
  \code{genomicState$codingGenome} (if
  \code{fullOrCoding="coding"}).
}
\description{
  This function extracts the coverage information
  calculated by \link{fullCoverage} for a set of exons
  determined by \link{makeGenomicState}. The underlying
  code is similar to \link{getRegionCoverage} with
  additional tweaks for calculating RPKM values.
}
\examples{
\dontrun{
## Collapse the coverage information
collapsedFull <- collapseFullCoverage(list(genomeData$coverage), verbose=TRUE)

## Calculate library size adjustments
sampleDepths <- sampleDepth(collapsedFull, probs=c(0.5), nonzero=TRUE, verbose=TRUE)

## Build the models
group <- genomeInfo$pop
adjustvars <- data.frame(genomeInfo$gender)
models <- makeModels(sampleDepths, testvars=group, adjustvars=adjustvars)

## Preprocess the data
## Automatic chunksize used to then compare 1 vs 4 cores in the 'do not run' section
prep <- preprocessCoverage(genomeData, groupInfo=group, cutoff=0, scalefac=32, chunksize=NULL, colsubset=NULL, mc.cores=4)

## Get the F statistics
fstats <- calculateStats(prep, models, mc.cores=1, verbose=TRUE)

## Determine a cutoff from the F-distribution.
## This step is very important and you should consider using quantiles from the observed F statistics
n <- dim(prep$coverageSplit[[1]])[2]
df1 <- dim(models$mod)[2]
df0 <- dim(models$mod0)[2]
cutoff <- qf(0.95, df1-df0, n-df1)

## Low cutoff used for illustrative purposes
cutoff <- 1

## Calculate the p-values and define the regions of interest.
regsWithP <- calculatePvalues(prep, models, fstats, nPermute=10, seeds=NULL, chr="chr21", cutoff=cutoff, mc.cores=1)

## Obtain fullCov object
datadir <- system.file("extdata", "genomeData", package="derfinder")
dirs <- makeBamList(datadir=datadir, samplepatt="*accepted_hits.bam$", bamterm=NULL)
## Shorten the column names
names(dirs) <- gsub("_accepted_hits.bam", "", names(dirs))

## Reading the data and filtering it is quite fast.
fullCov <- fullCoverage(dirs=dirs, chrnums="21", mc.cores=1)

## Create GenomicState object:
## Hsapiens.UCSC.hg19.knownGene GenomicState
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

## Creating this GenomicState object takes around 8 min
GenomicState.Hsapiens.UCSC.hg19.knownGene <- makeGenomicState(txdb=txdb)

## Annotate regions
annotatedRegions <- annotateRegions(regions=regsWithP$regions, genomicState=GenomicState.Hsapiens.UCSC.hg19.knownGene, minoverlap=1)

## Finally, get the coverage information for each exon
exonCov <- coverageToExon(fullCov=fullCov, genomicState=GenomicState.Hsapiens.UCSC.hg19.knownGene)
}
}
\author{
  Andrew Jaffe, Leonardo Collado-Torres
}
\seealso{
  \link{fullCoverage}, \link{getRegionCoverage}
}

