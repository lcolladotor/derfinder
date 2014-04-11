\name{coverageToExon}
\alias{coverageToExon}
\title{Extract coverage information for exons}
\usage{
coverageToExon(fullCov, genomicState, fullOrCoding = "full", L,
  returnType = "raw", mc.cores = getOption("mc.cores", 2L),
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

  \item{mc.cores}{This argument is passed to
  \link[parallel]{mclapply} twice. First, it's is used by
  strand. Secondly, for processing the exons by chromosome.
  So there is no gain in using \code{mc.cores} greater than
  the maximum of the number of strands and number of
  chromosomes.}

  \item{verbose}{If \code{TRUE} basic status updates will
  be printed along the way.}
}
\value{
A matrix (nrow = number of exons in \code{genomicState}
corresponding to the chromosomes in \code{fullCov}, ncol =
number of samples) with the number of reads (or RPKM) per
exon. The row names correspond to the row indexes of
\code{genomicState$fullGenome} (if
\code{fullOrCoding="full"}) or
\code{genomicState$codingGenome} (if
\code{fullOrCoding="coding"}).
}
\description{
This function extracts the coverage information calculated
by \link{fullCoverage} for a set of exons determined by
\link{makeGenomicState}. The underlying code is similar to
\link{getRegionCoverage} with additional tweaks for
calculating RPKM values.
}
\examples{
\dontrun{
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

## Creating this GenomicState object takes around 8 min for all chrs and around 30 secs for chr21
GenomicState.Hsapiens.UCSC.hg19.knownGene.chr21 <- makeGenomicState(txdb=txdb, chrs="chr21")

## Finally, get the coverage information for each exon
exonCov <- coverageToExon(fullCov=fullCov, genomicState=GenomicState.Hsapiens.UCSC.hg19.knownGene.chr21)
}
}
\author{
Andrew Jaffe, Leonardo Collado-Torres
}
\seealso{
\link{fullCoverage}, \link{getRegionCoverage}
}

