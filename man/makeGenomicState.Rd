\name{makeGenomicState}
\alias{makeGenomicState}
\title{Obtain the genomic state per region from annotation}
\usage{
makeGenomicState(txdb, chrs = paste0("chr", c(1:22, "X", "Y")),
  addChrPrefix = FALSE)
}
\arguments{
  \item{txdb}{A \link[GenomicFeatures]{TranscriptDb}
  object.}

  \item{chrs}{The names of the chromosomes to use as
  denoted in the \code{txdb} object. Check
  \link[GenomicFeatures]{isActiveSeq}.}

  \item{addChrPrefix}{If \code{TRUE}, 'chr' is added as a
  prefix to the chromosome names (seqlevels).}
}
\value{
A \code{GRangesList} object with two elements:
\code{fullGenome} and \code{codingGenome}. Both have
metadata information for the type of region (theRegion),
transcript IDs (tx_id), transcript name (tx_name), and gene
ID (gene_id). \code{fullGenome} classifies each region as
either being exon, intron or intragenic.
\code{codingGenome} classfies the regions as being
promoter, exon, intro, 5UTR, 3UTR or intragenic.
}
\description{
This function summarizes the annotation contained in a
\link[GenomicFeatures]{TranscriptDb} at each given base of
the genome based on annotated transcripts. It groups
contiguous base pairs classified as the same type into
regions.
}
\examples{
\dontrun{
## Hsapiens.UCSC.hg19.knownGene GenomicState
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

## Creating this GenomicState object takes around 8 min for all chrs and around 30 secs for chr21
GenomicState.Hsapiens.UCSC.hg19.knownGene.chr21 <- makeGenomicState(txdb=txdb, chrs="chr21")

## Hsapiens ENSEMBL GRCh37
library("GenomicFeatures")
## Can take several minutes and speed will depend on your internet speed
xx <- makeTxDbPackageFromBiomart(version = "0.99", maintainer = "Your Name", author="Your Name")
txdb <- loadDb(file.path("TxDb.Hsapiens.BioMart.ensembl.GRCh37.p11", "inst", "extdata", "TxDb.Hsapiens.BioMart.ensembl.GRCh37.p11.sqlite"))

## Creating this GenomicState object takes around 13 min
GenomicState.Hsapiens.ensembl.GRCh37.p11 <- makeGenomicState(txdb=txdb, chrs=c(1:22, "X", "Y"), addChrPrefix=TRUE)

## Save for later use
save(GenomicState.Hsapiens.ensembl.GRCh37.p11, file="GenomicState.Hsapiens.ensembl.GRCh37.p11.Rdata")
}
}
\author{
Andrew Jaffe, Leonardo Collado-Torres
}
\seealso{
\link[GenomicFeatures]{TranscriptDb}
}

