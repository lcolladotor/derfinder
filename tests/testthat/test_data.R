context('Loading and exporting data')

# Setup
datadir <- system.file('extdata', 'genomeData', package='derfinder')
files <- rawFiles(datadir = datadir, samplepatt = '*accepted_hits.bam$', 
    fileterm = NULL)
## Shorten the column names
names(files) <- gsub('_accepted_hits.bam', '', names(files))

## Find files
bogus <- rawFiles(datadir = datadir, samplepatt = "bw")
test_that('Finding files', {
    expect_that(rawFiles(datadir = NULL, sampledirs = NULL),
        throws_error("Either 'samplepatt' or 'sampledirs' must be non-NULL."))
    expect_that(length(bogus), equals(0))    
})

## Load BAM data
dataFilt <- loadCoverage(files = files, chr = '21', cutoff = 0)
dataRaw <- loadCoverage(files = files, chr = '21', cutoff = NULL)
dataFilt.mean <- loadCoverage(files = files, chr = '21', cutoff = 0,
    returnMean = TRUE)
which <- GRanges('21', IRanges(47410303, 47418067))

test_that('Load BAM data', {
    expect_that(dataFilt, is_identical_to(genomeData))
    expect_that(dataRaw, is_identical_to(genomeDataRaw))
    expect_that(names(dataFilt.mean), is_identical_to(c('coverage', 'position',
        'meanCoverage')))
    expect_that(sum(dataFilt.mean$meanCoverage), equals(357.838709677))
    expect_that(loadCoverage(files = files[1:2], chr = '21', cutoff = 0), is_identical_to(loadCoverage(files = files[1:2], chr = '21', cutoff = 0, which = which)))
})


## Export BigWig
dir.create('bw')

test_that('Create BigWig files', {
    expect_that(bws <- createBw(list('chr21' = dataRaw), path = 'bw'), 
        gives_warning())
    expect_that(is(bws, "GRangesList"), is_identical_to(TRUE))
    expect_that(length(bws), equals(31))
    expect_that(length(bws[[1]]), equals(16))
    expect_that(names(mcols(bws[[1]])), is_identical_to('score'))
    expect_that(dir('bw'), is_identical_to(c('ERR009101.bw', 'ERR009102.bw', 'ERR009105.bw', 'ERR009107.bw', 'ERR009108.bw', 'ERR009112.bw', 'ERR009115.bw', 'ERR009116.bw', 'ERR009131.bw', 'ERR009138.bw', 'ERR009144.bw', 'ERR009145.bw', 'ERR009148.bw', 'ERR009151.bw', 'ERR009152.bw', 'ERR009153.bw', 'ERR009159.bw', 'ERR009161.bw', 'ERR009163.bw', 'ERR009164.bw', 'ERR009167.bw', 'SRR031812.bw', 'SRR031835.bw', 'SRR031904.bw')))
})


## Load BigWig
bigwigs <- rawFiles(datadir = "bw", samplepatt = "bw", fileterm = NULL)
names(bigwigs) <- gsub('.bw', '', names(bigwigs))
dataBW <- loadCoverage(bigwigs, chr = 'chr21', cutoff = NULL,
    inputType = "BigWig")

dataRaw.bw <- dataRaw
dataRaw.bw$coverage <- dataRaw.bw$coverage[, names(bigwigs)]

test_that('Load BigWig data', {
    expect_that(dataRaw.bw, equals(dataBW))
    expect_that(loadCoverage(bigwigs, chr = '21', cutoff = NULL,
        inputType = "BigWig", verbose = FALSE), throws_error())
    expect_that(loadCoverage(bigwigs, chr = 'chr21', cutoff = NULL, 
        verbose = FALSE), equals(dataRaw.bw))
})

## Loading with GenomicFiles
dataRaw.GF <- loadCoverage(files, chr = '21', tilewidth = 2e7, cutoff = NULL)
dataRaw.bw.GF <- loadCoverage(bigwigs, chr = 'chr21', tilewidth = 2e7,
    cutoff = NULL, inputType = 'BigWig')
test_that('Load using GenomicFiles', {
    expect_that(dataRaw, is_identical_to(dataRaw.GF))
    expect_that(dataBW, is_identical_to(dataRaw.bw.GF))
})

## Loading with fullCoverage
fullCov <- fullCoverage(files = files, chrs = '21')
fullCov.bw <- fullCoverage(files = bigwigs, chrs = 'chr21', inputType = 'BigWig')
test_that('Load with fullCoverage()', {
    expect_that(list('chr21' = dataRaw$coverage), is_identical_to(fullCov))
    expect_that(list('chr21' = dataBW$coverage), is_identical_to(fullCov.bw))
    expect_that(fullCoverage(files = files, chrs = '21', output = c('one',
        'two')), throws_error())
    expect_that(fullCoverage(files = files, chrs = '21', chrlens = c(4e7,
        5e7)), throws_error())
})

## Loading with BamFile, BamFileList

bam <- system.file('extdata', 'genomeData', 'ERR009101_accepted_hits.bam', package = 'derfinder')
names(bam) <- 'sample1'
library('Rsamtools')
bFile <- BamFile(bam, index = paste0(bam, '.bai'))

library('rtracklayer')
big1 <- BigWigFile(bigwigs[1])
big1.list <- BigWigFileList(bigwigs[1])

test_that('BamFile and BamFileList', {
    expect_that(fullCoverage(bam, chrs = '21', verbose = FALSE), is_identical_to(fullCoverage(bFile, chrs = '21', sampleNames = 'sample1', verbose = FALSE)))
    expect_that(fullCoverage(bam, chrs = '21', verbose = FALSE), is_identical_to(fullCoverage(BamFileList(bFile), chrs = '21', sampleNames = 'sample1', verbose = FALSE)))
    expect_that(fullCoverage(bam, chrs = '21', sampleNames = 'ERR009101', verbose = FALSE), equals(fullCoverage(big1.list, chrs = 'chr21', verbose = FALSE)))
    expect_that(fullCoverage(big1.list, chrs = 'chr21', verbose = FALSE), is_identical_to(fullCoverage(big1, chrs = 'chr21', verbose = FALSE)))
})

## Filtering
filt <- filterData(dataRaw$coverage, cutoff = 0)
test_that('Filtering', {
    expect_that(filterData(dataRaw$coverage, cutoff = -1,
        verbose = FALSE)$coverage, is_identical_to(dataRaw$coverage))
    expect_that(filterData(dataRaw), throws_error())
    expect_that(filterData(dataRaw$coverage, filter = 'two'), throws_error())
    expect_that(filterData(dataRaw$coverage, cutoff = 10,
        verbose = FALSE)$position, is_identical_to(Rle(FALSE, 48129895)))
    expect_that(dim(filterData(dataRaw$coverage, cutoff = 0.5, filter = 'mean', 
        verbose = FALSE)$coverage), equals(c(273, 31)))
    expect_that(nrow(filt$coverage), is_identical_to(sum(filt$position)))
    
})


## Test getRegionCoverage() and coverageToExon()

## Assign chr lengths using hg19 information, use only first two regions
library('GenomicRanges')
data(hg19Ideogram, package = 'biovizBase', envir = environment())
regions <- genomeRegions$regions[1:2]
seqlengths(regions) <- seqlengths(hg19Ideogram)[names(seqlengths(regions))]
## Get the region coverage
regionCov <- getRegionCoverage(fullCov=fullCov, regions=regions)

## Load using which
fullCov.regs <- fullCoverage(files, chrs = seqlevels(regions), fileStyle = 'NCBI', protectWhich = 3e4, which = regions, verbose = FALSE)

## Load without using the default 'protectWhich'
max.noP <- sapply(loadCoverage(files = files, chr = '21', verbose = FALSE, protectWhich = 0, cutoff = NULL, which = regions)$coverage, max)
max.wP <- sapply(fullCov$chr21, max)

## Use only the first two exons
smallGenomicState <- genomicState
smallGenomicState$fullGenome <- smallGenomicState$fullGenome[which(smallGenomicState$fullGenome$theRegion == 'exon')[1:2] ]

## Get the coverage information for each exon
exonCov <- coverageToExon(fullCov=fullCov,
    genomicState=smallGenomicState$fullGenome, L=36)

## Verify coverageToExon
exonCov.verify <- do.call(rbind, lapply(getRegionCoverage(regions = smallGenomicState$fullGenome, files = files, fileStyle = 'NCBI', verbose = FALSE), function(x) { sapply(x, sum)})) / 36

test_that('Obtaining region coverage', {
    expect_that(fullCov, is_identical_to(fullCov.regs))
    expect_that(regionCov, is_identical_to(getRegionCoverage(regions = regions, files = files, fileStyle = 'NCBI', verbose = FALSE)))
    expect_that(exonCov, is_identical_to(coverageToExon(genomicState = smallGenomicState$fullGenome, L=36, files = files, fileStyle = 'NCBI', verbose = FALSE)))
    expect_that(max.noP - max.wP, is_equivalent_to(rep(c(0,-1,0,-1,-2,-1,-3,0,-3,-1,-3,0,-1,0,-1,0), c(4,1,1,1,3,2,2,1,3,1,2,1,1,3,1,4))))
    expect_that(exonCov, equals(exonCov.verify))
}
)
