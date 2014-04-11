\name{annotateRegions}
\alias{annotateRegions}
\title{Assign genomic states to regions}
\usage{
annotateRegions(regions, genomicState, minoverlap = 20,
  fullOrCoding = "full", annotate = TRUE, verbose = TRUE)
}
\arguments{
  \item{regions}{The \code{$regions} output from
  \link{calculatePvalues}.}

  \item{genomicState}{The output from
  \link{makeGenomicState}.}

  \item{minoverlap}{This parameter is passed to
  \link[IRanges]{countOverlaps} and determines the minimum
  overlap of a region to be assined a genomic state. Set to
  1 if you want all the genomic states that overlap the
  regions.}

  \item{fullOrCoding}{If \code{full} then the
  \code{genomicState$fullGenome} genomic state information
  is used. If \code{coding}, then the
  \code{genomicState$codingGenome} genomic state
  information is used.}

  \item{annotate}{If \code{TRUE} then the regions are
  annotated by the genomic state. Othewise, only the
  overlaps between the regions and the genomic states are
  computed.}

  \item{verbose}{If \code{TRUE} basic status updates will
  be printed along the way.}
}
\value{
A list with elements \code{countTable} and
\code{annotationList} (only if \code{annotate=TRUE}).
\describe{ \item{countTable }{This is a data.frame with the
number of overlaps from the regions vs the genomic states
with one type per column. For example, if
\code{fullOrCoding="full"} then the columns are
\code{exon}, \code{intragenic} and \code{intron}.}
\item{annotationList }{This is a \code{GRangesList} with
the genomic states that overlapped with the regions (if
any, depends on \code{minoverlap}). The names of this
\code{GRangesList} correspond to the region index in
\code{regions}.} }
}
\description{
This function takes the regions found in
\link{calculatePvalues} and assigns them genomic states
contructed with \link{makeGenomicState}. The main workhorse
functions are \link[IRanges]{countOverlaps} and
\link[IRanges]{findOverlaps}.
}
\examples{
\dontrun{
## Create GenomicState object:
## Hsapiens.UCSC.hg19.knownGene GenomicState
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

## Creating this GenomicState object takes around 8 min for all chrs and around 30 secs for chr21
GenomicState.Hsapiens.UCSC.hg19.knownGene.chr21 <- makeGenomicState(txdb=txdb, chrs="chr21")

## Annotate regions
annotatedRegions <- annotateRegions(regions=genomeRegions$regions, genomicState=GenomicState.Hsapiens.UCSC.hg19.knownGene.chr21, minoverlap=1)
}
}
\author{
Andrew Jaffe, Leonardo Collado-Torres
}
\seealso{
\link{makeGenomicState}, \link{calculatePvalues}
}

