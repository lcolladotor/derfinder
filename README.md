derfinder [![Build Status](https://travis-ci.org/lcolladotor/derfinder.svg?branch=master)](https://travis-ci.org/lcolladotor/derfinder)
=========

Annotation-agnostic fast differential expression analysis of RNA-seq data at base-pair resolution. For more information about `derfinder` check the vignettes [here](http://lcolladotor.github.io/derfinder/).


# Further documentation

You can generate HTML reports from the results using __regionReport__ 
available [here](https://github.com/lcolladotor/regionReport).

For a full example on how to use __derfinder__ check 
[derfinderExample](https://github.com/lcolladotor/derfinderExample). TODO: update this.

# Installation instructions

Get R 3.2.0 from [CRAN](http://cran.r-project.org/).

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

To cite package __derfinder__ in publications use:

Collado-Torres L, Frazee AC, Love MI, Irizarry RA, Jaffe AE and Leek JT (2015). "derfinder: Software for annotation-agnostic RNA-seq differential expression analysis". _bioRxiv_. <URL: http://dx.doi.org/10.1101/015370>, <URL:
http://www.biorxiv.org/content/early/2015/02/19/015370.abstract>.

Frazee AC, Sabunciyan S, Hansen KD, Irizarry RA and Leek JT (2014). “Differential expression analysis of RNA-seq data at
single-base resolution.” _Biostatistics_, *15 (3)*, pp. 413-426. <URL: http://dx.doi.org/10.1093/biostatistics/kxt053>, <URL:
http://biostatistics.oxfordjournals.org/content/15/3/413.long>.

A BibTeX entry for LaTeX users is

@Manual{,
    title = {derfinder: Software for annotation-agnostic RNA-seq differential expression analysis},
    author = {Leonardo Collado-Torres and Alyssa C. Frazee and Michael I. Love and Rafael A. Irizarry and Andrew E. Jaffe and Jeffrey T. Leek},
    year = {2015},
    journal = {bioRxiv},
    doi = {10.1101/015370},
    url = {http://www.biorxiv.org/content/early/2015/02/19/015370.abstract}
}

@Article{,
    title = {Differential expression analysis of RNA-seq data at single-base resolution},
    author = {Alyssa C. Frazee and Sarven Sabunciyan and Kasper D. Hansen and Rafael A. Irizarry and Jeffrey T. Leek},
    year = {2014},
    journal = {Biostatistics},
    volume = {15 (3)},
    pages = {413-426},
    doi = {10.1093/biostatistics/kxt053},
    url = {http://biostatistics.oxfordjournals.org/content/15/3/413.long},
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

Testing on R-devel for Bioc-devel is feasible thanks to [r-builder](https://github.com/metacran/r-builder).

# Origins

This is a development version for a faster version of the __derfinder__ core 
steps. The original implementation is available via GitHub at the 
[derfinder](https://github.com/alyssafrazee/derfinder) repository.
