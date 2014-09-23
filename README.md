derfinder [![Build Status](https://travis-ci.org/lcolladotor/derfinder.svg?branch=master)](https://travis-ci.org/lcolladotor/derfinder)
=========

Fast differential expression analysis of RNA-seq data at base-pair resolution. 
You can generate HTML reports from the results using __derfinderReport__ 
available [here](https://github.com/lcolladotor/derfinderReport).

For a full example on how to use __derfinder__ check 
[derfinderExample](https://github.com/lcolladotor/derfinderExample).

# Installation instructions

Get R 3.1.1 or newer from [CRAN](http://cran.r-project.org/).

```R
## If needed
install.packages('devtools')

## Pre-requisites from CRAN
install.packages(c('ggplot2', 'Hmisc', 'testthat', 'knitr', 'rmarkdown',
    'knitrBoostrap', 'knitcitations'))

## Pre-requisites from Bioconductor
source('http://bioconductor.org/biocLite.R')
biocLite(c('S4Vectors', 'IRanges', 'GenomicRanges', 'Rsamtools', 'bumphunter', 'biovizBase',
    'TxDb.Hsapiens.UCSC.hg19.knownGene', 'AnnotationDbi', 'GenomicFeatures', 'GenomeInfoDb',
    'rtracklayer', 'BiocParallel', 'qvalue', 'GenomicFiles'))
    
## derfinderHelper
library('devtools')
install_github('lcolladotor/derfinderHelper')

## derfinder itself

## If you are using BioC-devel use:
install_github('lcolladotor/derfinder@master')

## If you are using BioC-release use:
install_github('lcolladotor/derfinder@release')

## Suggested:
biocLite('ggbio')
install_github('lcolladotor/derfinderPlot')
```

# 'Watch' for updates

This software is in development, so we highly recommend 'watching' the 
repository: Click on the top right under `Watch`. You will then receive 
notifications for issues, comments, and pull requests as described 
[here](https://help.github.com/articles/notifications).

You will need a GitHub account to be able to `Watch` the repository.

# Citation

Below is the citation output from using `citation('derfinder')` in R. Please 
run this yourself to check for any updates on how to cite __derfinder__.

---

To cite package __derfinder__ in publications use:

Leonardo Collado-Torres, Alyssa Frazee, Andrew Jaffe and Jeffrey Leek (2014). 
derfinder: Fast differential expression analysis of RNA-seq data at base-pair 
resolution. R package version 0.0.74. https://github.com/lcolladotor/derfinder

A BibTeX entry for LaTeX users is

@Manual{,
    title = {derfinder: Fast differential expression analysis of RNA-seq data 
        at base-pair resolution},
    author = {Leonardo Collado-Torres and Alyssa Frazee and Andrew Jaffe 
        and Jeffrey Leek},
    year = {2014},
    note = {R package version 0.0.74},
    url = {https://github.com/lcolladotor/derfinder},
}


# Branches

* [__master__](https://github.com/lcolladotor/derfinder/tree/master) This 
branch corresponds to the one that works with the latest 
[Bioc-devel](http://master.bioconductor.org/packages/devel).
* [__release__](https://github.com/lcolladotor/derfinder/tree/release) This 
branch corresponds to the one that works with the latest [Bioc-release](http://master.bioconductor.org/packages/release). It is also
available at [derfinder-release](https://github.com/lcolladotor/derfinder-release) in order for the git-svn bridge (release version) to work as suggested by [Ilari Scheinin here](https://www.mail-archive.com/bioc-devel@r-project.org/msg01967.html).

## Old

* [__BioC-2.13__](https://github.com/lcolladotor/derfinder/tree/BioC-2.13) 
Version working with [BioC 2.13](http://master.bioconductor.org/packages/2.13)

## Travis CI

This package is automatically tested thanks to [Travis CI](travis-ci.org) and [r-travis](https://github.com/craigcitro/r-travis). If you want to add this to your own package use:

```R
## Use devtools to create the .travis.yml file
library('devtools')
use_travis('yourPackage')

## Read https://github.com/craigcitro/r-travis/wiki to configure .travis.yml appropriately

## Add a status image by following the info at http://docs.travis-ci.com/user/status-images/
```

# Origins

This is a development version for a faster version of the __derfinder__ core 
steps. The original implementation is available via GitHub at the 
[derfinder](https://github.com/alyssafrazee/derfinder) repository.
