\name{plotCluster}
\alias{plotCluster}
\title{Plot the coverage information surrounding a region cluster}
\usage{
  plotCluster(idx, regions, annotation, coverageInfo,
    groupInfo, titleUse = "qval", txdb = NULL,
    p.ideogram = NULL, minExtend = 200, scalefac = 32)
}
\arguments{
  \item{idx}{A integer specifying the index number of the
  region of interest.}

  \item{regions}{The \code{$regions} output from
  \link{calculatePvalues}.}

  \item{annotation}{The output from running
  \link[bumphunter]{annotateNearest} on the output from
  \link{calculatePvalues}.}

  \item{coverageInfo}{A DataFrame resulting from
  \link{loadCoverage} using \code{cutoff=NULL}.}

  \item{groupInfo}{A factor specifying the group membership
  of each sample. It will be used to color the samples by
  group.}

  \item{titleUse}{Whether to show the p-value (\code{pval})
  or the q-value (\code{qval}) in the title. If
  \code{titleUse=none} then no p-value or q-value
  information is used; useful if no permutations were
  performed and thus p-value and q-value information is
  absent.}

  \item{txdb}{A transcript data base such as
  TxDb.Hsapiens.UCSC.hg19.knownGene If \code{NULL} then no
  annotation information is used.}

  \item{p.ideogram}{If \code{NULL}, the ideogram for hg19
  is built for the corresponding chromosome. Otherwise an
  ideogram resuling from \link[ggbio]{plotIdeogram}.}

  \item{minExtend}{The minimum number of base-pairs to
  extend the view before and after the region of interest.}

  \item{scalefac}{A log transformation is used on the count
  tables, so zero counts present a problem. Ideally, the
  same you provided to \link{preprocessCoverage}.}
}
\value{
  A ggplot2 plot that is ready to be printed out.
  Tecnically it is a ggbio object.
}
\description{
  For a given region found in \link{calculatePvalues}, plot
  the coverage for the cluster this region belongs to as
  well as some padding. The mean by group is shown to
  facilitate comparisons between groups. If annotation
  exists, you can plot the trancripts and exons (if any)
  overlapping in the vicinity of the region of interest.
}
\details{
  See the parameter \code{significantCut} in
  \link{calculatePvalues} for how the significance cutoffs
  are determined.
}
\examples{
## Get coverage info without any cutoff
datadir <- system.file("extdata", "genomeData", package="derfinder2")
dirs <- makeBamList(datadir=datadir, samplepatt="*accepted_hits.bam$", bamterm=NULL)
names(dirs) <- gsub("_accepted_hits.bam", "", names(dirs))
covInfo <- loadCoverage(dirs=dirs, chr="21", cutoff=NULL, verbose=FALSE)

## Construct the models
group <- genomeInfo$pop
adjustvars <- data.frame(genomeInfo$gender)
models <- makeModels(coverageInfo=genomeData, testvars=group, adjustvars=adjustvars, nonzero=TRUE)

## Preprocess the data
prep <- preprocessCoverage(genomeData, cutoff=0, scalefac=32, chunksize=NULL, colsubset=NULL, mc.cores=4)

## Get the F statistics
fstats <- calculateStats(prep, models, mc.cores=1, verbose=FALSE)

## Determine a cutoff from the F-distribution.
## This step is very important and you should consider using quantiles from the observed F statistics
\dontrun{
n <- dim(prep$coverageSplit[[1]])[2]
df1 <- dim(models$mod)[2]
df0 <- dim(models$mod0)[2]
cutoff <- qf(0.95, df1-df0, n-df1)
}
## Low cutoff used for illustrative purposes
cutoff <- 1

## Calculate the p-values and define the regions of interest.
regsWithP <- calculatePvalues(prep, models, fstats, nPermute=10, seeds=NULL, chr="chr21", cutoff=cutoff, mc.cores=1, verbose=FALSE)

## Annotate the results
suppressMessages(library("bumphunter"))
annotation <- annotateNearest(regsWithP$regions, "hg19")

## Make the plot
suppressMessages(library("TxDb.Hsapiens.UCSC.hg19.knownGene"))
plotCluster(idx=1, regions=regsWithP$regions, annotation=annotation, coverageInfo=covInfo$coverage, groupInfo=group, txdb=TxDb.Hsapiens.UCSC.hg19.knownGene)
## Resize the plot window and the labels will look good.

\dontrun{
## For a custom plot, check the ggbio and ggplot2 packages.
## Also feel free to look at the code for this function:
plotCluster
}
}
\author{
  Leonardo Collado-Torres
}
\seealso{
  \link{loadCoverage}, \link{calculatePvalues},
  \link[bumphunter]{annotateNearest},
  \link[ggbio]{plotIdeogram}
}

