# derfinder 1.19.8

* Added a `NEWS.md` file to track changes to the package.


# derfinder 1.19.4

BUG FIXES

* `railMatrix()` and `loadCoverage()` helper functions had an issue when
the input set of regions was duplicated. This could be reproduced with

```
sampleFile <- c('SRR387777' = 'http://duffel.rail.bio/recount/SRP009615/bw/SRR387777.bw')
regs <- GenomicRanges::GRanges('chrY', IRanges(start = c(1, 1), width = 10), strand = '-')
names(regs) <- c(1:2)
result <- rtracklayer::import(sampleFile, selection = regs, as = 'RleList')
```

This error affected recount and other reverse dependencies that use
derfinder for processing BigWig files.

# derfinder 1.19.2


BUG FIXES

* `railMatrix()` and `loadCoverage()` helper functions now attempt to import a
BigWig file 3 times before giving up. Based on
http://bioconductor.org/developers/how-to/web-query/ and
https://github.com/leekgroup/recount/commit/8da982b309e2d19638166f263057d9f85bb64e3f
which will make these functions more robust to occasional web access
issues.

# derfinder 1.17.3


NEW FEATURES

* Add ORCID's following changes at
http://bioconductor.org/developers/package-guidelines/#description

# derfinder 1.17.2


BUG FIXES

* Use R's random seeds from version 3.5.0 for the test
thanks to https://twitter.com/StrictlyStat/status/1103303028751372289

# derfinder 1.15.4


SIGNIFICANT USER-VISIBLE CHANGES

* Switched from `outfile` to `log` when invoking
`BiocParallel::SnowParam()`. Thus `define_cluster()`
now has a `mc.log` argument instead of `mc.outfile`.


# derfinder 1.15.2


BUG FIXES

* Fix a message regarding the deprecated `IRanges` `subset` method.


# derfinder 1.15.1


SIGNIFICANT USER-VISIBLE CHANGES

* Use `BiocManager`


# derfinder 1.13.8


BUG FIXES

* Fixed a unit test that was breaking version 1.13.7.

# derfinder 1.13.1


NEW FEATURES

* Added an extra example to `regionMatrix()` in response to
https://support.bioconductor.org/p/103591

# derfinder 1.11.8


BUG FIXES

* Improved the documentation regarding an error when coverageInfo$position
is NULL when running `analyzeChr()` and indirectly running
`preprocessCoverage()`. See https://support.bioconductor.org/p/99400/
for details.


# derfinder 1.11.7


SIGNIFICANT USER-VISIBLE CHANGES

* Vignette now uses the new `BiocStyle::html_document` that was recently
released.

# derfinder 1.11.4


BUG FIXES

* `regionMatrix()` will now pass the hidden arguments `species` and
`currentStyle` to `getRegionCoverage()` so they can be used by
`extendedMapSeqlevels()`. Related to
https://support.bioconductor.org/p/95721/.

# derfinder 1.11.2


BUG FIXES

* Improved documentation of `extendedMapSeqlevels()`. Related to
https://support.bioconductor.org/p/95521/.
* Improved `filterData()` based on
https://github.com/lcolladotor/derfinder/issues/38


# derfinder 1.9.6


BUG FIXES

* Fixed `define_cluster()` to match recent changes in BiocParallel and
fixed an if clause in `regionMatrix()` that could lead to warnings in
some situations.

# derfinder 1.9.3


SIGNIFICANT USER-VISIBLE CHANGES

* `regionMatrix()` now has explicit arguments `totalMapped` and `targetSize`
so that users will almost always normalize by library size when
using this function (if they see the help page) or in the steps
prior to using `regionMatrix()`.

# derfinder 1.9.2


BUG FIXES

* Clarified the documentation of `mc.cores` and `mc.cores.load` in
`fullCoverage()` thanks to feedback from Emily E Burke
https://github.com/emilyburke.


# derfinder 1.7.16


SIGNIFICANT USER-VISIBLE CHANGES

* Help pages now document advanced arguments.
* Deprecated `advancedArg()`.


# derfinder 1.7.14


NEW FEATURES

* Added the function `getTotalMapped()` for calculating the total number of
mapped reads for a BAM file or the area under the curve (AUC) for a
BigWig file. This information can then be used with fullCoverage(),
`filterData()` and other functions. Note that if you `totalMapped`
in `fullCoverage()` you should not use `totalMapped` again in
`filterData()`.


# derfinder 1.7.12


BUG FIXES

* Updated links to BrainSpan. Issue reported by Steve Semick
https://github.com/SteveSemick.


# derfinder 1.7.2


BUG FIXES

* Now derfinder uses `DataFrame(check.names = FALSE)` to avoid naming issues.


# derfinder 1.7.1


SIGNIFICANT USER-VISIBLE CHANGES

* Dropped defunct functions.


# derfinder 1.5.39


BUG FIXES

* `annotateRegions()` now ignores strand by default.

# derfinder 1.5.37


SIGNIFICANT USER-VISIBLE CHANGES

* The users guide vignette now has a short section explaining how to use
derfinder for differential binding analysis with ChIP-seq data.


# derfinder 1.5.27


SIGNIFICANT USER-VISIBLE CHANGES

* Added smoothing arguments for the single base-level approach based on
functions from the bumphunter package. These arguments are useful for
identifying differentially bounded ChIP-seq peaks.


# derfinder 1.5.13


BUG FIXES

* Fixed `railMatrix()`'s flexibility for defining the cluster used for loading
the BigWig files. You can now use `BPPARAM.railChr` which will take
priority over `file.cores`. Also, if `file.cores = 1L`, then the default
will be to use `SerialParam()`, which was the implementation available
prior to 1.5.11.


# derfinder 1.5.11


SIGNIFICANT USER-VISIBLE CHANGES

* Now `coverageToExon()`, `regionMatrix()` and `railMatrix()` can take an `L`
argument of length equal to the number of samples in case not all
samples have the same read length.
* `railMatrix()` has a new argument called `file.cores` for controlling how
many cores are used for loading the BigWig files. In theory this allows
using `railMatrix()` with `BPPARAM.custom` equal to a
`BiocParallel::BatchJobsParam()` to submit 1 job per chromosome, then
`file.cores` determines the number of cores for reading the files.
This is a highly experimental feature.

# derfinder 1.5.9


SIGNIFICANT USER-VISIBLE CHANGES

* Dropped the old introductory and advanced vignettes. We think that the
new vignettes are clearer. In particular, they do a better job at
highlighting the differences between the expressed regions-level and
single base-level F-statistics implementations of the DER Finder
approach to RNA-seq data.

# derfinder 1.5.8


SIGNIFICANT USER-VISIBLE CHANGES

* Added a users guide vignette which explains nearly every detail you would
want to know as a user.


# derfinder 1.5.7


SIGNIFICANT USER-VISIBLE CHANGES

* Added a quick start vignette.


# derfinder 1.5.6


NEW FEATURES

* Introduced `railMatrix()` which generates similar output to `regionMatrix()`
but is much faster and less memory intensive. It achieves this by
extracting the required information from BigWig files.

# derfinder 1.3.3


SIGNIFICANT USER-VISIBLE CHANGES

* Brought back the `mc.outfile` argument for specifying the `outfile`
argument in `SnowParam()`. See more details at
https://stat.ethz.ch/pipermail/bioc-devel/2015-May/007531.html

# derfinder 1.3.2


SIGNIFICANT USER-VISIBLE CHANGES

* Deprecated functions with underscores in their names in favor of
camelCase functions. This was done to simplify the package.

# derfinder 1.3.1


SIGNIFICANT USER-VISIBLE CHANGES

* Greatly increased the speed of the p-values calculation step. See
https://github.com/lcolladotor/derfinder/issues/29 for details.


# derfinder 1.1.18


BUG FIXES

* Updated to work with qvalue 1.99.0.


# derfinder 1.1.17


SIGNIFICANT USER-VISIBLE CHANGES

* Changed citation information to reference the bioRxiv pre-print.

BUG FIXES

* Polished the interaction with bumphunter >= 1.7.3.
* Updated the default cluster option now that `BiocParallel::SnowParam()`
no longer has an `outfile` argument.

# derfinder 1.1.16


SIGNIFICANT USER-VISIBLE CHANGES

* `analyzeChr()` now uses `annotateTranscripts()` and `matchGenes()` from
bumphunter version 1.7.3 (or greater). As announced at
https://support.bioconductor.org/p/63568/ these changes in bumphunter
allow straight forward use of non-human annotation. In `analyzeChr()`
using a different organism can be used by changing the `txdb` argument:
finer control can be achieved through `...`. For example, by specifying
the `annotationPackage` argument used in `annotateTranscripts()`.

# derfinder 1.1.15


BUG FIXES

* `makeGenomicState()` incorrectly labeled regions as intragenic. The correct
name is intergenic.


# derfinder 1.1.14


BUG FIXES

* Fixed an important bug on `calculatePvalues()`! Basically, internally
`maxRegionGap` was set to 300 instead of 0 in one step by default.
Thus the process of mapping regions to genomic coordinates was messed
up. If you have results prior to this fix you can try using
https://gist.github.com/bf85e2c7d5d1f8197707 to fix the results as
much as possible. Basically, regions will be correct but the
p-values will be approximated with the available information from the
null regions. Truly fixing the p-values can only be done by re-running
derfinder.


# derfinder 1.1.5


NEW FEATURES

* Introduced function `extendedMapSeqlevels()` for using `GenomeInfoDb` when
there is information regarding the species and naming style of interest.
otherwise sequence names are left unchanged. If used with
`verbose = TRUE`, a message is printed whenever `GenomeInfoDb` could not
be used or if some information had to be guessed.

BUG FIXES

* Fixes https://support.bioconductor.org/p/62136.

# derfinder 1.1.3


NEW FEATURES

* `loadCoverage()` and `fullCoverage()` now support `BamFile` and `BigWigFile`
objects.

BUG FIXES

* Fixed a bug in `loadCoverage()` when the input was a `BamFileList`.
Implemented tests based on the bug. Bug reported at
https://support.bioconductor.org/p/62073


# derfinder 0.99.0


NEW FEATURES

* Preparing to submit to Bioconductor.

# derfinder 0.0.81


NEW FEATURES

* Added an advanced vignette.

# derfinder 0.0.80


NEW FEATURES

* Introductory vignette completed.

# derfinder 0.0.79


SIGNIFICANT USER-VISIBLE CHANGES

* `mergeResults()` can now calculate FWER adjusted p-values when provided with
`optionsStats`. Updated `analyzeChr()` to supply the required information.


# derfinder 0.0.78


NEW FEATURES

* Added an introductory vignette.

# derfinder 0.0.77


NEW FEATURES

* Added `advancedArg()` and its alias `advanced_arg()` which links to the docs
for the advanced arguments by opening a browser window with the relevant
information from GitHub.
* `getRegionCoverage()` and `coverageToExon()` now have the `files` argument
which is used only when `fullCov` is `NULL`. Both functions will attempt
to extract the coverage data from the raw files for the regions of
interest in that case.

Special care has to be taken in order to guarantee that the coverage is
the same as some reads might be discarded if the region is too narrow.
See the advanced argument `protectWhich` in `loadCoverage()` for more
information. Also, if `totalMapped` and `targetSize` were used prior to
filtering, they should be used again.
* `loadCoverage()` has new advanced arguments that help when reading a
specific region (or regions) of the genome.

# derfinder 0.0.76


SIGNIFICANT USER-VISIBLE CHANGES

* `loadCoverage()` and `fullCoverage()` argument `dirs` has been renamed to
`files` for greater consistency with what it represents.

# derfinder 0.0.75


SIGNIFICANT USER-VISIBLE CHANGES

* `regionMatrix()` now returns the output of `getRegionCoverage()` so you don't
have to run it twice if you are interested in using
`derfinderPlot::plotRegionCoverage()`.
* `regionMatrix()$regions` now guesses the seqlengths

BUG FIXES

* Fixed `.advanced_argument()` to work in nested functions
* Fixed a case in `getRegionCoverage()` where `fullCov$position` was provided
but it was `NULL`.
* Fixed `regionMatrix(totalMapped, targetSize)` case which would previously
lead to an error in the `getRegionCoverage()` step.


# derfinder 0.0.74


SIGNIFICANT USER-VISIBLE CHANGES

* Most functions had their arguments changed to increase usability. Some
have advanced arguments inside the code.


# derfinder 0.0.71


SIGNIFICANT USER-VISIBLE CHANGES

* Introduced `coerceGR()`, `createBwSample()` and `createBw()` for exporting
output from `fullCoverage()` into BigWig files.


# derfinder 0.0.70


SIGNIFICANT USER-VISIBLE CHANGES

* All exported functions now have aliases with underscore names for those
that prefer them over camel case names. For example, `analyze_chr()` is
the new alias for `analyzeChr()`.
* `makeBamList()` has been renamed to `rawFiles()` since it can be used to
identify a list of BigWig files instead of BAM files.


# derfinder 0.0.69


SIGNIFICANT USER-VISIBLE CHANGES

* `loadCoverage()` and `fullCoverage()` now have a `tilewidth` argument. When
specified `GenomicFiles` is used to read the coverage in chunks. In
theory, this can lead to lower memory usage at the expense of time.

# derfinder 0.0.67


SIGNIFICANT USER-VISIBLE CHANGES

* `preprocessCoverage()` now has a `toMatrix` argument which is only used when
`lowMemDir` is not `NULL`. It controls whether to save the chunks as
`DataFrame` objects or `dgCMatrix` objects and the idea is that it can time
by just transforming the data once instead of doing so at each
permutation.

# derfinder 0.0.66


SIGNIFICANT USER-VISIBLE CHANGES

* `fstats.apply()` has been moved to it's own package: `derfinderHelper`. This
will speed up the run time when using `BiocParallel::SnowParam()` as
`derfinderHelper` takes much less time to load than derfinder.
* `plotCluster()`, `plotOverview()` and `plotRegionCoverage()` were all moved to
their new own package: `derfinderPlot`. This will make maintenance easier
as the dependency `ggbio` is still under active development.



# derfinder 0.0.65


* Re-organized code for `fstats.apply()`. Note that improving
`.transformSparseMatrix()` would speed up the `Matrix` method.
* Note that all parallel functions have some overhead from loading derfinder
on each worker. Check `system.time(library(derfinder))` to see how long
the overhead is. It only pays off to use more cores if the calculations
are taking longer than the overhead.


# derfinder 0.0.64


SIGNIFICANT USER-VISIBLE CHANGES

* `fullCoverage()` has several new arguments and now is a full parallel
implementation of `loadCoverage()`. These changes were introduced since
`fullCoverage()` no longer blows up in memory since version 0.0.62 and
thus the new recommended use case is to call `fullCoverage()` instead of
running one job with `loadCoverage()` per chromosome or using a `lapply()`
loop.

# derfinder 0.0.63


NEW FEATURES

* `fstats.apply()` now has `method` and `scalefac` arguments. The `method`
argument controls which of the 3 implementations to use. The old method
is called `regular` now. The new method `Rle` calculates the
F-statistics without de-compressing the data, which is good for memory
but gets considerably slower as the number of samples increases. The
default method is `Matrix` which uses the Matrix package and is both
faster (given that the coercion doesn't take long) and less memory
intensive than the `regular` method.

SIGNIFICANT USER-VISIBLE CHANGES

* Functions `analyzeChr()`, `calculatePvalues()` and `calculateStats()` now have
arguments `method` and `scalefac` to match the changes in `fstats.apply()`.

# derfinder 0.0.62


NEW FEATURES

* derfinder now uses `BiocParallel::blapply()` instead of `parallel::mclapply()`
When `mc.cores` is greater than 1, `BiocParallel::SnowParam()` is used to
construct the cluster. Otherwise, `BiocParallel::SerialParam()` is used.
This change reduces memory load when using the functions that have the
`mc.cores` argument greater than 1.

* Functions `analyzeChr()`, `calculatePvalues()`, `calculateStats()`,
`coverageToExon()`, `fullCoverage()`, `getRegionCoverage()`, `regionMatrix()`
all have a new argument `mc.output`. This is passed to
`BiocParallel::SnowParam(outfile)`.

SIGNIFICANT USER-VISIBLE CHANGES

* You may now use `fullCoverage()` without problems and should no longer
encounter errors due to longer vectors not being implemented.
* Functions like `fullCoverage()` now use much less memory and do not blow up
as you increase `mc.cores`. Note however that the memory does increase,
but now it`s close to linear.
* Examples might take longer to run with `mc.cores` greater than 1, but that
is due to the small setup overhead of `BiocParallel::SnowParam()` which
is minimal compared to the overall speed gains with real data sets.

# derfinder 0.0.61


SIGNIFICANT USER-VISIBLE CHANGES

* `filterData()` and `loadCoverage()` now have arguments `totalMapped` and
`targetSize`
* `getRegionCoverage()` and `regionMatrix()` can now work with list output from
`loadCoverage()` with a non-NULL cutoff
* `regionMatrix()` now has an argument `runFilter` so it can be used with
previous output from `loadCoverage()`/`filterData()` with `returnMean=TRUE`

# derfinder 0.0.60


SIGNIFICANT USER-VISIBLE CHANGES

* `analyzeChr()`, `annotateRegions()`, `calculatePvalues()`, `coverageToExon()`,
`findRegions()`, `fullCoverage()`, `getRegionCoverage()`,
`makeGenomicState()`, `mergeResults()`, `plotCluster()`, `plotOverview()`
now all have `chrsStyle` as an argument to specify the chromosome
naming convention used. Defaults to UCSC.
* `makeGenomicState()` no longer has the `addChrPrefix` argument. It has been
replaced by `chrsStyle` to use `GenomeInfoDb` to set the naming style.
* `chrnums` has been renamed to `chrs` in `fullCoverage()` and `mergeResults()`
* `chrnum` has been renamed to `chr` in `analyzeChr()`

# derfinder 0.0.59


NEW FEATURES

* Updated to work with BioC version 3.0

# derfinder 0.0.58


NEW FEATURES

* `loadCoverage()` and `fullCoverage()` can now import data from BigWig files.

# derfinder 0.0.57


SIGNIFICANT USER-VISIBLE CHANGES

* `regionMatrix()` now relies on `getRegionCoverage()` instead of
`coverageToExon()` making it faster and less memory intensive.

# derfinder 0.0.56


NEW FEATURES

* Added `regionMatrix()` for filtering coverage data and using the resulting
regions to construct a count matrix. Uses several derfinder functions.

SIGNIFICANT USER-VISIBLE CHANGES

* Made `coverageToExon()` more robust for different names in `fullCov`.
* `filterData()` and `loadCoverage()` have new arguments `filter`, `returnMean`,
and `returnCoverage` which allow speeding up `regionMatrix()`.
`preprocessCoverage()` was changed accordingly.
* `getRegionCoverage()` now internally uses USCS names.

BUG FIXES

* Fixed warnings in `coverageToExon()`.


# derfinder 0.0.55


SIGNIFICANT USER-VISIBLE CHANGES

* Added a `NEWS` file with curated information from the git commit history.

# derfinder 0.0.54


SIGNIFICANT USER-VISIBLE CHANGES

* Added example data for `mergeResults()`. Now all functions have examples.

# derfinder 0.0.53


SIGNIFICANT USER-VISIBLE CHANGES

* Now derfinder includes the genomic state output for hg19 chr21. This
allowed to implement examples for `annotateRegions()`,
`getRegionCoverage()`, and `coverageToExon()`.


# derfinder 0.0.52


SIGNIFICANT USER-VISIBLE CHANGES

* Several examples run much faster now

# derfinder 0.0.51


SIGNIFICANT USER-VISIBLE CHANGES

* `makeGenomicState()` now has an example


# derfinder 0.0.50


BUG FIXES

* Fixed warnings in `makeGenomicState()` related to changes in `AnnotationDbi`


# derfinder 0.0.49


BUG FIXES

* Updated the help for `plotRegionCoverage()` and fixed an issue.

# derfinder 0.0.48


SIGNIFICANT USER-VISIBLE CHANGES

* `loadCoverage()` now allows specifying which strand you want to load. More
at https://github.com/lcolladotor/derfinder/issues/16

# derfinder 0.0.47


SIGNIFICANT USER-VISIBLE CHANGES

* `getRegionCoverage()` now has a depth-adjustment argument

BUG FIXES

* Fixed several bugs as suggested by Andrew Jaffe


# derfinder 0.0.46


SIGNIFICANT USER-VISIBLE CHANGES

* Now requires R version 3.1
* Changed Github organization to match the Git-SVN Bioconductor bridge
guidelines.
* Dropped `Rcpp` and `RcppArmadillo` from F-stats calculation. More at
https://github.com/lcolladotor/derfinder/pull/17
* Stored the results from some of the example data to speed up other
examples. Check `?genomeDataRaw`, `?genomeFstats`, `?genomeRegions`

BUG FIXES

* Fixed `verbose` for `getRegionCoverage()`
* `plotRegionCoverage()` now matches latest `getRegionCoverage()` output

# derfinder 0.0.44


SIGNIFICANT USER-VISIBLE CHANGES

* Updated `getRegionCoverage()` with a new method for sub setting the coverage
matrices, allowing for coverage estimates from overlapping regions. Now
also uses `mclapply()`.

BUG FIXES

* Fixed `NAMESPACE` to match current bioc-devel (2.14) as suggested by Tim
Triche.

# derfinder 0.0.42


SIGNIFICANT USER-VISIBLE CHANGES

* A series of performance enhancements were made to reduce the memory load
 (albeit a very minor time increase).

BUG FIXES

* Updated `analyzeChr()` to handle correctly the new `lowMemDir` argument.

# derfinder 0.0.41


SIGNIFICANT USER-VISIBLE CHANGES

* Introduced the `lowMemDir` argument to `preprocessCoverage()`,
`calculateStats()`, `calculatePvalues()`, `fstats.apply()`, and `analyzeChr()`.
Reduces peak memory usage at the expense of some input-output.

# derfinder 0.0.40


SIGNIFICANT USER-VISIBLE CHANGES

* `mergeResults()` will not merge pre-processed data by default
* `coverageToExon()` now uses `mclapply()` when possible

# derfinder 0.0.39


SIGNIFICANT USER-VISIBLE CHANGES

* Changed how the data is pre-processed and pre-splitted.

# derfinder 0.0.38


SIGNIFICANT USER-VISIBLE CHANGES

* `preprocessCoverage()` now uses `Reduce()` instead of `.rowMeans()`

# derfinder 0.0.37


SIGNIFICANT USER-VISIBLE CHANGES

* Updated the example in `plotCluster()`

# derfinder 0.0.36


NEW FEATURES

* Added `collapseFullCoverage()`

SIGNIFICANT USER-VISIBLE CHANGES

* `sampleDepth()` has been greatly changed. It is now based on Hector
Corrada`s ideas implemented in metagenomeSeq.
* Updated several man pages.

# derfinder 0.0.35


BUG FIXES

* Fixed an issue with the bumphunter dependency.

# derfinder 0.0.34


BUG FIXES

* Merged changes suggested by Michael Love


# derfinder 0.0.33


BUG FIXES

* Implemented fixes suggested by Michael Love

# derfinder 0.0.32


BUG FIXES

* Fixed an important bug in the F-stats calculation
* Implemented fixes suggested by Michael Love

# derfinder 0.0.31


BUG FIXES

* `loadCoverage()` now uses `readGAlignmentsFromBam()`

# derfinder 0.0.30


SIGNIFICANT USER-VISIBLE CHANGES

* Added a `bai` argument to `fullCoverage()`
* `loadCoverage()` can now work with a pre-defined BamFileList object.

# derfinder 0.0.29


SIGNIFICANT USER-VISIBLE CHANGES

* Now requires R version 3.0.2
* Changed the default of `center` in `sampleDepth()` to `FALSE`.
* Added a `runAnnotation` argument to `analyzeChr()`.
* Added a `bai` argument to `loadCoverage()`.
* Added an `adjustF` argument to all stats functions. Useful for cases when
the RSS of the alternative model is very small.

BUG FIXES

* Fixed `plotRegionCoverage()` and `plotCluster()` for unexpected cases.


# derfinder 0.0.28


NEW FEATURES

* Added `sampleDepth()`

SIGNIFICANT USER-VISIBLE CHANGES

* `generateReport()` has been moved to it's own new package called
`derfinderReport`. It is available at
https://github.com/lcolladotor/derfinderReport
* Examples and `analyzeChr()` have been updated now that `sampleDepth()` was
added


# derfinder 0.0.26


SIGNIFICANT USER-VISIBLE CHANGES

* Renamed package from `derfinder2` to `derfinder` to comply with Bioconductor
guidelines.


# derfinder 0.0.25


SIGNIFICANT USER-VISIBLE CHANGES

* Changed how `makeModels()` deals with cases when mod and mod0 are not full
rank.
* `plotCluster() `now no longer depends on an active Internet connection for
`hg19 = TRUE`.

BUG FIXES

* Fixed minor graphical issues in `plotRegionCoverage()`, `plotCluster()` and
in `generateReport()`

# derfinder 0.0.24


BUG FIXES

* Fixed `plotRegionCoverage()` for cases when in `annotateRegions(minoverlap=x)`
lead to no overlaps being found between a region and annotation.

# derfinder 0.0.23


BUG FIXES

* Fixed some bugs in `calculatePvalues()` when no null regions or only some
were found.
* Fixed a bug using `colsubset` on `analyzeChr()`.
* Fixed an issue when the F-stat cutoff used is too high and no regions are
found.
* Fixed an issue when `testvars` in `makeModels()` had unused levels.
* Fixed an issue when `qvalue::qvalue()` fails due to incorrect estimation of
`pi0`
* Fixed a bug on how the regions were being clustered.

# derfinder 0.0.22


SIGNIFICANT USER-VISIBLE CHANGES

* MA-style plots in `generateReport()` now weight the mean by the number of
samples in each group. Also removed the mean coverage vs area section.
`generateReport()` now also has a `nBestClusters` argument.
* `plotCluster() `now uses scales and has a `forceLarge` argument

BUG FIXES

* Fixed an important bug in finding candidate DERs. The example
for `plotCluster()` now includes code that was used for visualizing the
bug.
* Fixed some NAMESPACE issues

# derfinder 0.0.21


NEW FEATURES

* Added `getRegionCoverage()`, `coverageToExon()`, `plotRegionCoverage()`

SIGNIFICANT USER-VISIBLE CHANGES

* Renamed `plotRegion()` to `plotCluster()` plus it no longer shows the exons
track as it is redundant information
* `mergeResults()` now also runs `annotateRegions()`
* `generateReport()` now uses `plotRegionCoverage()` and includes MA-style plots

# derfinder 0.0.20


NEW FEATURES

* Added a`nnotateRegions()`

# derfinder 0.0.19


NEW FEATURES

* Added `makeGenomicState()`

# derfinder 0.0.18


NEW FEATURES

* Added `fullCoverage()`

# derfinder 0.0.17


SIGNIFICANT USER-VISIBLE CHANGES

* `preprocessCoverage()` now uses the `groupInfo` argument
* `calculatePvalues()` now calculates log2 fold changes (without scaling or
adjusting for library size)
* Greatly improved `generateReport()`

# derfinder 0.0.16


NEW FEATURES

* Added `generateReport()`

BUG FIXES

* Fixed bugs in `mergeResults()`


# derfinder 0.0.15


SIGNIFICANT USER-VISIBLE CHANGES

* `getSegmentsRle()` was greatly simplified

BUG FIXES

* Fixed `analyzeChr()`, completed `mergeResults()`


# derfinder 0.0.14


NEW FEATURES

* Added `analyzeChr()` and `mergeResults()`

SIGNIFICANT USER-VISIBLE CHANGES

* Updated README.md
* `makeModels()` now uses `testvars` instead of `group` and has a new
arguments `groupInfo`, `center` and `testIntercept`
* `calculatePvalues()` now uses area of regions instead of mean to calculate
the p-values.
* `preprocessCoverage()` now calculates the mean coverage at each base

# derfinder 0.0.13


SIGNIFICANT USER-VISIBLE CHANGES

* Users can specify significance cutoffs for `plotOverview()` and `plotRegion()`

# derfinder 0.0.12


SIGNIFICANT USER-VISIBLE CHANGES

* `calculatePvalues()` now uses `qvalue::qvalue()` instead of `p.adjust()`

# derfinder 0.0.11


SIGNIFICANT USER-VISIBLE CHANGES

* Added `plotOverview()`

# derfinder 0.0.10


SIGNIFICANT USER-VISIBLE CHANGES

* Added `plotRegion()`

BUG FIXES

* `makeModels()` can now handle a vector for the `adjustvars` argument

# derfinder 0.0.9


SIGNIFICANT USER-VISIBLE CHANGES

* `calculatePvalues()` will adjust the `p.values` now using `p.adjust()`
* `makeModels()` now can handle a matrix for the `group` argument

BUG FIXES

* `getSegmentsRle()` will now work properly in the case that no segments are
found

# derfinder 0.0.8


SIGNIFICANT USER-VISIBLE CHANGES

* Attempted to reduce memory load in `calculateStats()` and `calculatePvalues()`

# derfinder 0.0.7


SIGNIFICANT USER-VISIBLE CHANGES

* `fstats.apply()` now uses `Rcpp` and `RcppArmadillo`

# derfinder 0.0.6


NEW FEATURES

* `preprocessCoverage()` can now automatically select the `chunksize`

SIGNIFICANT USER-VISIBLE CHANGES

* Dropped `fstats()` and is now part of `fstats.apply()`

# derfinder 0.0.5


NEW FEATURES

* Introduced `method` argument for getSegmetnsRle()

BUG FIXES

* No longer permuting the data matrix in `calculatePvalues()`

# derfinder 0.0.4


NEW FEATURES

* Added `makeBamList()`, `makeModels()`, and `preprocessCoverage()`
* Added example data from the Montogemery and Pickrell studies. Check
`?genomeData` and `?genomeInfo`

# derfinder 0.0.3


NEW FEATURES

* Added `calculatePvalues()`

# derfinder 0.0.2


NEW FEATURES

* Added `clusterMakerRle()`, `findRegions()`, and `getSegmentsRle()`

# derfinder 0.0.1


NEW FEATURES

* Added `calculateStats()`, `filterData()`, `fstats()`, and `fstats.apply()`

SIGNIFICANT USER-VISIBLE CHANGES

* Renamed `makeCoverage()` to `loadCoverage()`
* Improved `NAMESPACE`

# derfinder 0.0.0


NEW FEATURES

* Initialized the package (named `derfinder2`) from `derfinder`
https://github.com/alyssafrazee/derfinder version 1.0.2
This version is available at
https://github.com/alyssafrazee/derfinder/tree/d49f7b28c26f075da36a50ab67c9d192ab2fd63d

