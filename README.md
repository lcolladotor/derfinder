derfinder2
==========

Fast differential expression analysis of RNA-seq data at base-pair resolution.

This is a development version for a faster version of the derfinder core steps. The original implementation is available via GitHub at the [derfinder](https://github.com/alyssafrazee/derfinder) repository.

# Installation instructions

Get R 3.0.1 or newer from [CRAN](http://cran.r-project.org/).

```S
## If needed
install.packages("devtools")

## Pre-requisites from CRAN
install.packages(c("knitr", "Rcpp", "RcppArmadillo", "ggplot2", "reshape2", "plyr", "microbenchmark"))

## Pre-requisites from Bioconductor
source("http://bioconductor.org/biocLite.R")
biocLite(c("IRanges", "GenomicRanges", "Rsamtools", "bumphunter", "biovizBase", "ggbio", "qvalue", "TxDb.Hsapiens.UCSC.hg19.knownGene"))

## derfinder2 itself
library(devtools)
install_github("derfinder2", "lcolladotor")
```

Note that the current Bioconductor release version of __bumphunter__ for R 3.0.1 is a few versions before the one required by __derfinder2__. The version needed can be installed manually from http://bioconductor.org/packages/2.13/bioc/html/bumphunter.html You can download the source or other binaries.
