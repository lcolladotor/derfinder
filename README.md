derfinder [![Build Status](https://travis-ci.org/lcolladotor/derfinder.svg?branch=master)](https://travis-ci.org/lcolladotor/derfinder)
=========

Annotation-agnostic fast differential expression analysis of RNA-seq data at base-pair resolution. For more information about `derfinder` check the vignettes [here](http://lcolladotor.github.io/derfinder/).


# Further documentation

You can generate HTML reports from the results using __derfinderReport__ 
available [here](https://github.com/lcolladotor/derfinderReport).

For a full example on how to use __derfinder__ check 
[derfinderExample](https://github.com/lcolladotor/derfinderExample). TODO: update this.

# Installation instructions

Get R 3.1.1 or newer from [CRAN](http://cran.r-project.org/).

```R
## From Bioconductor
source('http://bioconductor.org/biocLite.R')
biocLite('derfinder')

## Suggested:
biocLite(c('derfinderPlot', 'regionReport'))
```

# Vignette

The vignette for this package can be viewed [here](http://lcolladotor.github.io/derfinder/) or via [Bioconductor's website](http://www.bioconductor.org/packages/devel/bioc/html/derfinder.html).

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

Leonardo Collado-Torres, Alyssa C. Frazee, Andrew E. Jaffe and Jeffrey T. Leek (2014). derfinder: Annotation-agnostic differential expression analysis of RNA-seq data at base-pair resolution. R package version 1.1.1. https://github.com/lcolladotor/derfinder

A BibTeX entry for LaTeX users is

@Manual{,
    title = {derfinder: Annotation-agnostic differential expression analysis of RNA-seq
    data at base-pair resolution},
    author = {Leonardo Collado-Torres and Alyssa C. Frazee and Andrew E. Jaffe and Jeffrey T. Leek},
    year = {2014},
    note = {R package version 1.1.1},
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
