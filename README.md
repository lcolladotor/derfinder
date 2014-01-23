derfinder
==========

Fast differential expression analysis of RNA-seq data at base-pair resolution. You can generate HTML reports from the results using __derfinderReport__ available at https://github.com/lcolladotor/derfinderReport.

For a full example on how to use __derfinder__ check https://github.com/lcolladotor/derfinderExample.

# Installation instructions

Get R 3.0.2 or newer from [CRAN](http://cran.r-project.org/).

```S
## If needed
install.packages("devtools")

## Pre-requisites from CRAN
install.packages(c("Rcpp", "RcppArmadillo", "ggplot2", "reshape2", "plyr", "microbenchmark",
	 "RColorBrewer", "scales", "Hmisc"))

## Pre-requisites from Bioconductor
source("http://bioconductor.org/biocLite.R")
biocLite(c("IRanges", "GenomicRanges", "Rsamtools", "bumphunter", "biovizBase", "ggbio", "qvalue",
	 "TxDb.Hsapiens.UCSC.hg19.knownGene", "AnnotationDbi", "GenomicFeatures"))

## derfinder itself
library(devtools)
install_github("derfinder", "lcolladotor")
```

Note that the current Bioconductor release version of __bumphunter__ for R 3.0.2 is a few versions before the one required by __derfinder__. The version needed can be installed manually from http://bioconductor.org/packages/devel/bioc/html/bumphunter.html You can download the source or other binaries.

# 'Watch' for updates

This software is in development, so we highly recommend 'watching' the repository: Click on the top right under `Watch`. You will then receive notifications for issues, comments, and pull requests as described [here](https://help.github.com/articles/notifications).

You will need a GitHub account to be able to `Watch` the repository.

# Citation

Below is the citation output from using `citation("derfinder")` in R. Please run this yourself to check for any updates on how to cite __derfinder__.

---

To cite package __derfinder__ in publications use:

Leonardo Collado-Torres, Alyssa Frazee, Andrew Jaffe and Jeffrey Leek (2013). derfinder: Fast differential expression analysis of RNA-seq data at base-pair resolution. R package version 0.0.45. https://github.com/lcolladotor/derfinder

A BibTeX entry for LaTeX users is

@Manual{, title = {derfinder: Fast differential expression analysis of RNA-seq data at base-pair resolution}, author = {Leonardo Collado-Torres and Alyssa Frazee and Andrew Jaffe and Jeffrey Leek}, year = {2013}, note = {R package version 0.0.45}, url = {https://github.com/lcolladotor/derfinder}, }


# Branches

* __master__ This branch corresponds to the one that works with the latest Bioc-release.
* __BioC-2.13__ Version working with BioC 2.13
* __BioC-2.14__ Version working with BioC 2.14
* __devel__ Developer's corner for testing ideas. Any useful ideas will be added to the version corresponding to the latest BioC-devel.


# Origins

This is a development version for a faster version of the __derfinder__ core steps. The original implementation is available via GitHub at the [derfinder](https://github.com/alyssafrazee/derfinder) repository.
