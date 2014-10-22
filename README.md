derfinder [![Build Status](https://travis-ci.org/lcolladotor/derfinder.svg?branch=master)](https://travis-ci.org/lcolladotor/derfinder)
=========

Annotation-agnostic fast differential expression analysis of RNA-seq data at base-pair resolution. For more information about `derfinder` check the vignettes [here](http://lcolladotor.github.io/derfinder/).


# Further documentation

You can generate HTML reports from the results using __regionReport__ 
available [here](https://github.com/lcolladotor/regionReport).

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

# Vignettes

The vignettes for this package can be viewed [here](http://lcolladotor.github.io/derfinder/) or via [Bioconductor's website](http://www.bioconductor.org/packages/devel/bioc/html/derfinder.html).

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

Leonardo Collado-Torres, Alyssa C. Frazee, Andrew E. Jaffe and Jeffrey T. Leek (2014). derfinder: Annotation-agnostic differential expression analysis of RNA-seq data at base-pair resolution. R package version 1.1.5. https://github.com/lcolladotor/derfinder

A BibTeX entry for LaTeX users is

@Manual{,
    title = {derfinder: Annotation-agnostic differential expression analysis of RNA-seq
    data at base-pair resolution},
    author = {Leonardo Collado-Torres and Alyssa C. Frazee and Andrew E. Jaffe and Jeffrey T. Leek},
    year = {2014},
    note = {R package version 1.0.3},
    url = {https://github.com/lcolladotor/derfinder},
}


# Travis CI

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
