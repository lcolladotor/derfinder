\name{plotOverview}
\alias{plotOverview}
\title{Plot a karyotype overview of the genome with the identified regions}
\usage{
plotOverview(regions, annotation = NULL, type = "pval", base_size = 12,
  areaRel = 4, legend.position = c(0.85, 0.12), significantCut = c(0.05,
  0.1))
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

  \item{legend.position}{This argument is passed to
  \link[ggplot2]{theme}. From ggplot2: the position of
  legends. ("left", "right", "bottom", "top", or
  two-element numeric vector).}

  \item{significantCut}{A vector of length two specifiying
  the cutoffs used to determine significance. The first
  element is used to determine significance for the
  p-values and the second element is used for the
  q-values.}
}
\value{
A ggplot2 plot that is ready to be printed out. Tecnically
it is a ggbio object.
}
\description{
Plots an overview of the genomic locations of the
identified regions (see \link{calculatePvalues}) in a
karyotype view. The coloring can be done either by
significant regions according to their p-values,
significant by adjusted p-values, or by annotated region if
using \link[bumphunter]{annotateNearest}.
}
\examples{
## Construct toy data
chrs <- paste0("chr", c(1:22, "X", "Y"))
chrs <- factor(chrs, levels=chrs)
suppressMessages(library("GenomicRanges"))
regs <- GRanges(rep(chrs, 10), ranges=IRanges(runif(240, 1, 4e7), width=1e3), significant=sample(c(TRUE, FALSE), 240, TRUE, p=c(0.05, 0.95)), significantQval=sample(c(TRUE, FALSE), 240, TRUE, p=c(0.05, 0.95)), area=rnorm(240))
annotation <- data.frame(region=sample(c("upstream", "promoter", "overlaps 5'", "inside", "overlaps 3'", "close to 3'", "downstream"), 240, TRUE))

## Type pval
plotOverview(regs)

## Type qval
plotOverview(regs, type="qval")

## Annotation
plotOverview(regs, annotation, type="annotation")

## Resize the plots if needed.

\dontrun{
## You might prefer to leave the legend at ggplot2's default option: right
plotOverview(regs, legend.position="right")

## Although the legend looks better on the bottom
plotOverview(regs, legend.position="bottom")

## Example knitr chunk for higher res plot using the CairoPNG device
```{r overview, message=FALSE, fig.width=7, fig.height=9, dev="CairoPNG", dpi=300}
plotOverview(regs, base_size=30, areaRel=10, legend.position=c(0.95, 0.12))
```

## For more custom plots, take a look at the ggplot2 and ggbio packages
## and feel free to look at the code of this function:
plotOverview
}
}
\author{
Leonardo Collado-Torres
}
\seealso{
\link{calculatePvalues}, \link[bumphunter]{annotateNearest}
}

