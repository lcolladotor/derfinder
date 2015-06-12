<a href="http://www.bioconductor.org/packages/release/bioc/html/derfinder.html#archives"><img border="0" src="http://www.bioconductor.org/shields/availability/release/derfinder.svg" title="Whether the package is available on all platforms; click for details."></a> <a href="http://www.bioconductor.org/packages/release/bioc/html/derfinder.html#since"><img border="0" src="http://www.bioconductor.org/shields/years-in-bioc/derfinder.svg" title="How long since the package was first in a released Bioconductor version (or is it in devel only)."></a> <a href="http://bioconductor.org/packages/stats/bioc/derfinder.html"><img border="0" src="http://www.bioconductor.org/shields/downloads/derfinder.svg" title="Percentile (top 5/20/50% or 'available') of downloads over last 6 full months. Comparison is done across all package categories (software, annotation, experiment)."></a> <a href="https://support.bioconductor.org/t/derfinder/"><img border="0" src="http://www.bioconductor.org/shields/posts/derfinder.svg" title="Support site activity, last 6 months: tagged questions/avg. answers per question/avg. comments per question/accepted answers, or 0 if no tagged posts."></a> <a href="http://www.bioconductor.org/packages/release/bioc/html/derfinder.html#svn_source"><img border="0" src="http://www.bioconductor.org/shields/commits/bioc/derfinder.svg" title="average Subversion commits (to the devel branch) per month for the last 6 months"></a>

Build status: Travis [![Build Status](https://travis-ci.org/lcolladotor/derfinder.svg?branch=master)](https://travis-ci.org/lcolladotor/derfinder),
Bioc-release <a href="http://bioconductor.org/checkResults/release/bioc-LATEST/derfinder/"><img border="0" src="http://www.bioconductor.org/shields/build/release/bioc/derfinder.svg" title="build results; click for full report"></a>,
Bioc-devel <a href="http://bioconductor.org/checkResults/devel/bioc-LATEST/derfinder/"><img border="0" src="http://www.bioconductor.org/shields/build/devel/bioc/derfinder.svg" title="build results; click for full report"></a>

derfinder
=========

Annotation-agnostic fast differential expression analysis of RNA-seq data at base-pair resolution. For more information about `derfinder` check the vignettes [here](http://www.bioconductor.org/packages/derfinder).


# Further documentation

You can generate HTML reports from the results using __regionReport__ 
available [here](https://github.com/lcolladotor/regionReport).

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

The vignettes for this package can be viewed [here](http://lcolladotor.github.io/derfinder/) or via [Bioconductor's website](http://www.bioconductor.org/packages/derfinder).

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


# Testing

Testing on Bioc-devel is feasible thanks to [r-builder](https://github.com/metacran/r-builder) as well as Bioconductor's nightly build.

# Origins

This is a development version for a faster version of the __derfinder__ core 
steps. The original implementation is available via GitHub at the 
[derfinder](https://github.com/alyssafrazee/derfinder) repository.
