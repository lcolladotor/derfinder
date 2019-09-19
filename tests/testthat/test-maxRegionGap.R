library('GenomicRanges')
library('bumphunter')

########### from ?findRegions

## Preprocess the data
prep <- preprocessCoverage(genomeData, cutoff=0, scalefac=32, chunksize=1e3,
    colsubset=NULL)

## Get the F statistics
fstats <- genomeFstats

## Find the regions
regs <- findRegions(prep$position, fstats, 'chr21', verbose=TRUE)


#### Scenario:
## Start from a regular GenomicRanges
## Use filterData() and findRegions() to collapse regions

## First, set as a GRanges
gr <- regs

## Filter the data to get positions above a cutoff
filt <- filterData(GenomicRanges::coverage(gr), cutoff = 0,
    returnCoverage = TRUE, returnMean = TRUE)

## Find regions with no maxRegionGap
nogap_derfinder <- findRegions(filt$position, filt$meanCoverage, chr = 'chr21',
    cutoff = 0)
nogap_bumphunter <- bumphunter::regionFinder(x = as.numeric(filt$meanCoverage),
    chr = rep('chr21', length(filt$meanCoverage)),
    pos = which(filt$position > 0), cutoff = 0, maxGap = 1)
## order the results
nogap_bumphunter <- nogap_bumphunter[order(nogap_bumphunter$start), ]

## Now with a gap
gap_derfinder <- findRegions(filt$position, filt$meanCoverage, chr = 'chr21',
    cutoff = 0, maxRegionGap = 50)
gap_bumphunter <- bumphunter::regionFinder(x = as.numeric(filt$meanCoverage),
    chr = rep('chr21', length(filt$meanCoverage)),
    pos = which(filt$position > 0), cutoff = 0, maxGap = 50)
gap_bumphunter <- gap_bumphunter[order(gap_bumphunter$start), ]

## Remove the excess stuff
remove_derfinder <- function(x) {
    mcols(x) <- mcols(x)[, c('value', 'area', 'indexStart', 'indexEnd')]
    return(x)
}
remove_bumphunter <- function(x) {
    x <- GRanges(x[, - which(colnames(x) %in% c('cluster', 'clusterL', 'L'))])
    names(x) <- rep('up', length(x))
    return(x)
}

## now check that the regions match (start and ends)
## as well as the rest except for the cluster stuff
test_that("maxRegionGap works", {
    expect_equal(
        remove_derfinder(nogap_derfinder),
        remove_bumphunter(nogap_bumphunter)
    )
    expect_equal(
        remove_derfinder(gap_derfinder),
        remove_bumphunter(gap_bumphunter)
    )
})
