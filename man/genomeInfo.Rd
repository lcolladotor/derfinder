\docType{data}
\name{genomeInfo}
\alias{genomeInfo}
\title{Genome samples information}
\format{A data.frame with 5 columns:
\describe{
\item{run }{ The short name used to identify the sample BAM file.}
\item{library.layout }{ Whether it was a single-end library or a paired-end library.}
\item{hapmap.id }{ The HapMap identifier of the person sequenced. Note that some were sequenced more than once.}
\item{gender }{ Whether the person sequence is a female or a male.}
\item{pop }{ The population the person belongs to.}
}}
\description{
Information for the 31 samples downloaded from the Short
Read Archive from studies comparing CEU and YRI
populations. This data is used to specify the adjustment
variables in \link{calculateStats}. The data is sorted
according to the BAM files identifiers. Reads were 36bp
long.
}
\details{
The samples are from: \describe{ \item{10 }{ unrelated
females from the YRI population.} \item{5 }{ unrelated
females from the CEU population.} \item{5 }{ unrelated
males (unrelated to the females too) from the CEU
population.} }
}
\references{
1. Pickrell JK, Marioni JC, Pai AA, Degner JF, Engelhardt
BE, Nkadori E, Veyrieras J-B, Stephens M, Gilad Y,
Pritchard JK. Understanding mechanisms underlying human
gene expression variation with RNA sequencing. Nature 2010
Apr.

2. Montgomery SB, Sammeth M, Gutierrez-Arcelus M, Lach RP,
Ingle C, Nisbett J, Guigo R, Dermitzakis ET. Transcriptome
genetics using second generation sequencing in a Caucasian
population. Nature 2010 Mar.
}
\seealso{
\link{genomeData}, \link{calculateStats}
}
\keyword{datasets}

