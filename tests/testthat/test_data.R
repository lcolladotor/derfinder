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
    expect_that(rawFiles(datadir = NULL, samplefiles = NULL),
        throws_error("Either 'samplepatt' or 'samplefiles' must be non-NULL."))
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
