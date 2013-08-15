\name{plotOverview}
\alias{plotOverview}
\title{Plot a karyotype overview of the genome with the identified regions}
\usage{
  plotOverview(regions, annotation = NULL, type = "pval",
    base_size = 12, areaRel = 4,
    significantCut = c(0.05, 0.1))
}
\arguments{
  \item{regions}{The \code{$regions} output from
  \link{calculatePvalues}.}

  \item{annotation}{The output from running
  \link[bumphunter]{annotateNearest} on the output from
  \link{calculatePvalues}. It is only required if
  \code{type="annotation"}.}

  \item{type}{Must be either \code{pval}, \code{qval} or
  \code{annotation}. It determines whether the plot
  coloring should be done according to significant p-values
  (<0.05), significant q-values (<0.10) or annotation
  regions.}

  \item{base_size}{Base point size of the plot. This
  argument is passed to \link[ggplot2]{element_text}
  (\code{size} argument).}

  \item{areaRel}{The relative size for the area label when
  \code{type="pval"} or \code{type="qval"}. Can be useful
  when making high resolution versions of these plots in
  devices like CairoPNG.}

  \item{significantCut}{A vector of length two specifiying
  the cutoffs used to determine significance. The first
  element is used to determine significance for the
  p-values and the second element is used for the
  q-values.}
}
\value{
  A ggplot2 plot that is ready to be printed out.
  Tecnically it is a ggbio object.
}
\description{
  Plots an overview of the genomic locations of the
  identified regions (see \link{calculatePvalues}) in a
  karyotype view. The coloring can be done either by
  significant regions according to their p-values,
  significant by adjusted p-values, or by annotated region
  if using \link[bumphunter]{annotateNearest}.
}
\examples{
## Construct the models
group <- genomeInfo$pop
adjustvars <- data.frame(genomeInfo$gender)
models <- makeModels(coverageInfo=genomeData, group=group, adjustvars=adjustvars, nonzero=TRUE)

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

## Overview with type pval
plotOverview(regions=regsWithP$regions, type="pval")

## Overview with type qval
plotOverview(regions=regsWithP$regions, type="qval")

\dontrun{
## Annotate the results
suppressMessages(library("bumphunter"))
annotation <- annotateNearest(regsWithP$regions, "hg19")

## Annotation overview
plotOverview(regions=regsWithP$regions, annotation=annotation, type="annotation")

## This function is meant to be an example of the plots you can make with the output from calculatePvalues()
## For more details or custom plots check the ggbio and ggplot2 packages
## as well as the code from this function:
plotOverview
}
}
\author{
  Leonardo Collado-Torres
}
\seealso{
  \link{calculatePvalues},
  \link[bumphunter]{annotateNearest}
}

