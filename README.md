
<!-- README.md is generated from README.Rmd. Please edit that file -->

# derfinder

<!-- badges: start -->

[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![Bioc release
status](http://www.bioconductor.org/shields/build/release/bioc/derfinder.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/derfinder)
[![Bioc devel
status](http://www.bioconductor.org/shields/build/devel/bioc/derfinder.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/derfinder)
[![Bioc downloads
rank](https://bioconductor.org/shields/downloads/release/derfinder.svg)](http://bioconductor.org/packages/stats/bioc/derfinder/)
[![Bioc
support](https://bioconductor.org/shields/posts/derfinder.svg)](https://support.bioconductor.org/tag/derfinder)
[![Bioc
history](https://bioconductor.org/shields/years-in-bioc/derfinder.svg)](https://bioconductor.org/packages/release/bioc/html/derfinder.html#since)
[![Bioc last
commit](https://bioconductor.org/shields/lastcommit/devel/bioc/derfinder.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/derfinder/)
[![Bioc
dependencies](https://bioconductor.org/shields/dependencies/release/derfinder.svg)](https://bioconductor.org/packages/release/bioc/html/derfinder.html#since)
[![Codecov test
coverage](https://codecov.io/gh/lcolladotor/derfinder/branch/master/graph/badge.svg)](https://codecov.io/gh/lcolladotor/derfinder?branch=master)
[![R build
status](https://github.com/lcolladotor/derfinder/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/lcolladotor/derfinder/actions)
[![GitHub
issues](https://img.shields.io/github/issues/lcolladotor/derfinder)](https://github.com/lcolladotor/derfinder/issues)
[![GitHub
pulls](https://img.shields.io/github/issues-pr/lcolladotor/derfinder)](https://github.com/lcolladotor/derfinder/pulls)
<!-- badges: end -->

Annotation-agnostic differential expression analysis of RNA-seq data at
base-pair resolution via the DER Finder approach. This package contains
two different implementations of this approach. The first one is the
single base-level F-statistics implementation and the second one is via
identifying expressed regions. For more information about `derfinder`
check the vignettes
[here](http://www.bioconductor.org/packages/derfinder).

## Documentation

For more information about `derfinder` check the vignettes [through
Bioconductor](http://bioconductor.org/packages/derfinder) or at the
[documentation website](http://lcolladotor.github.io/derfinder).

## Further documentation

You can generate HTML reports from the results using **regionReport**
available [here](https://github.com/lcolladotor/regionReport).

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `derfinder` from
[Bioconductor](http://bioconductor.org/) using the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("derfinder")
```

## Citation

Below is the citation output from using `citation('derfinder')` in R.
Please run this yourself to check for any updates on how to cite
**derfinder**.

``` r
print(citation("derfinder"), bibtex = TRUE)
#> To cite package 'derfinder' in publications use:
#> 
#>   Collado-Torres L, Nellore A, Frazee AC, Wilks C, Love MI, Langmead B,
#>   Irizarry RA, Leek JT, Jaffe AE (2017). "Flexible expressed region
#>   analysis for RNA-seq with derfinder." _Nucl. Acids Res._.
#>   doi:10.1093/nar/gkw852 <https://doi.org/10.1093/nar/gkw852>,
#>   <http://nar.oxfordjournals.org/content/early/2016/09/29/nar.gkw852>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {Flexible expressed region analysis for RNA-seq with derfinder},
#>     author = {Leonardo Collado-Torres and Abhinav Nellore and Alyssa C. Frazee and Christopher Wilks and Michael I. Love and Ben Langmead and Rafael A. Irizarry and Jeffrey T. Leek and Andrew E. Jaffe},
#>     year = {2017},
#>     journal = {Nucl. Acids Res.},
#>     doi = {10.1093/nar/gkw852},
#>     url = {http://nar.oxfordjournals.org/content/early/2016/09/29/nar.gkw852},
#>   }
#> 
#>   Frazee AC, Sabunciyan S, Hansen KD, Irizarry RA, Leek JT (2014).
#>   "Differential expression analysis of RNA-seq data at single-base
#>   resolution." _Biostatistics_, *15 (3)*, 413-426.
#>   doi:10.1093/biostatistics/kxt053
#>   <https://doi.org/10.1093/biostatistics/kxt053>,
#>   <http://biostatistics.oxfordjournals.org/content/15/3/413.long>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {Differential expression analysis of RNA-seq data at single-base resolution},
#>     author = {Alyssa C. Frazee and Sarven Sabunciyan and Kasper D. Hansen and Rafael A. Irizarry and Jeffrey T. Leek},
#>     year = {2014},
#>     journal = {Biostatistics},
#>     volume = {15 (3)},
#>     pages = {413-426},
#>     doi = {10.1093/biostatistics/kxt053},
#>     url = {http://biostatistics.oxfordjournals.org/content/15/3/413.long},
#>   }
#> 
#>   Collado-Torres L, Jaffe AE, Leek JT (2017). _derfinder:
#>   Annotation-agnostic differential expression analysis of RNA-seq data
#>   at base-pair resolution via the DER Finder approach_.
#>   doi:10.18129/B9.bioc.derfinder
#>   <https://doi.org/10.18129/B9.bioc.derfinder>,
#>   https://github.com/lcolladotor/derfinder - R package version 1.35.0,
#>   <http://www.bioconductor.org/packages/derfinder>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {derfinder: Annotation-agnostic differential expression analysis of RNA-seq data at base-pair resolution via the DER Finder approach},
#>     author = {Leonardo Collado-Torres and Andrew E. Jaffe and Jeffrey T. Leek},
#>     year = {2017},
#>     url = {http://www.bioconductor.org/packages/derfinder},
#>     note = {https://github.com/lcolladotor/derfinder - R package version 1.35.0},
#>     doi = {10.18129/B9.bioc.derfinder},
#>   }
```

Please note that the `derfinder` was only made possible thanks to many
other R and bioinformatics software authors, which are cited either in
the vignettes and/or the paper(s) describing this package.

## DER Finder versions

- The original implementation of the DER Finder approach as published in
  Frazee et al, Biostatistics 2014 is available via GitHub at
  [derfinder](https://github.com/leekgroup/derfinder).
- The version implementing the single base-level approach via
  calculating F-stastics as described in the pre-print Collado-Torres et
  al, Nucleic Acids Research 2017 is available via Bioconductor at
  [derfinder](http://bioconductor.org/packages/derfinder). The same
  package has the functions required for the expressed regions-level
  approach.

## Code of Conduct

Please note that the derfinder project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

## Development tools

- Continuous code testing is possible thanks to [GitHub
  actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)
  through *[usethis](https://CRAN.R-project.org/package=usethis)*,
  *[remotes](https://CRAN.R-project.org/package=remotes)*,
  *[sysreqs](https://github.com/r-hub/sysreqs)* and
  *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)* customized
  to use [Bioconductor’s docker
  containers](https://www.bioconductor.org/help/docker/) and
  *[BiocCheck](https://bioconductor.org/packages/3.17/BiocCheck)*.
- Code coverage assessment is possible thanks to
  [codecov](https://codecov.io/gh) and
  *[covr](https://CRAN.R-project.org/package=covr)*.
- The [documentation website](http://lcolladotor.github.io/derfinder) is
  automatically updated thanks to
  *[pkgdown](https://CRAN.R-project.org/package=pkgdown)*.
- The code is styled automatically thanks to
  *[styler](https://CRAN.R-project.org/package=styler)*.
- The documentation is formatted thanks to
  *[devtools](https://CRAN.R-project.org/package=devtools)* and
  *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

For more details, check the `dev` directory.

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.17/biocthis)*.
