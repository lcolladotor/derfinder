context('Finding regions')

library('GenomicRanges')

pos <- Rle(c(TRUE, FALSE, TRUE, FALSE, TRUE), c(1000, 1000, 1000, 200, 1000))
test_that('.clusterMakerRle', {
    expect_that(derfinder:::.clusterMakerRle(pos), equals(Rle(c(1L, 2L), c(1000, 2000))))
    expect_that(derfinder:::.clusterMakerRle(pos, ranges = TRUE), equals(IRanges(c(1, 1001), c(1000, 3000))))
    expect_that(derfinder:::.clusterMakerRle(pos), equals(derfinder:::.clusterMakerRle(pos, maxGap = 999)))
    expect_that(derfinder:::.clusterMakerRle(pos, maxGap = 1000), equals(Rle(1L, 3000)))
})


x <- Rle(c(-1, 0, 1, 2), rep(1000, 4))
test_that('.getSegmentsRle', {
    expect_that(derfinder:::.getSegmentsRle(x, 0.5), equals(list(upIndex = IRanges(2001, 4000), dnIndex = IRanges(1, 1000))))
    expect_that(derfinder:::.getSegmentsRle(x + 3, c(1, 1)), equals(list(upIndex = IRanges(1, 4000), dnIndex = IRanges())))
})


prep <- preprocessCoverage(genomeData, cutoff=0, scalefac=32, chunksize=1e3, 
    colsubset=NULL)
prep.small <- preprocessCoverage(genomeData, cutoff=1, scalefac=32, 
    chunksize=1e3, colsubset=1:2, groupInfo = factor(letters[1:2]))
data.small <- filterData(genomeData$coverage[1:2], 1, genomeData$position)
tmp <- sapply(data.small$coverage, identity)
names(tmp) <- c('a', 'b')

test_that('preprocessCoverage', {
    expect_that(genomeData$position, equals(prep$position))
    expect_that(DataFrame(sapply(genomeData$coverage, function(x) { log2(x + 32)})), equals(prep$coverageProcessed))
    expect_that(list(Rle(c(TRUE, FALSE), c(1000, 434)), Rle(c(FALSE, TRUE), c(1000, 434))), equals(prep$mclapplyIndex))
    expect_that(Reduce('+', genomeData$coverage) / ncol(genomeData$coverage), equals(prep$meanCoverage))
    expect_that(data.small$position, equals(prep.small$position))
    expect_that(DataFrame(sapply(data.small$coverage, function(x) { log2(x + 32)})), equals(prep.small$coverageProcessed))
    expect_that(tmp, equals(prep.small$groupMeans))
})


regs <- findRegions(prep$position, genomeFstats, 'chr21', verbose=TRUE)
regs.small <- findRegions(prep$position[Rle(TRUE, sum(47407536, 32, 1256, 36))], genomeFstats[seq_len(32 + 36)], 'chr21', verbose=TRUE, cutoff = 1)
regs.small.basic <- findRegions(prep$position[Rle(TRUE, sum(47407536, 32, 1256, 36))], genomeFstats[seq_len(32 + 36)], 'chr21', verbose=TRUE, cutoff = 1, basic = TRUE)
sums <- sapply(split(genomeFstats[seq_len(32 + 36)][genomeFstats[seq_len(32 + 36)] > 1], rep(letters[1:3], c(6, 15, 36))), sum)
basic.check <- DataFrame(area = Rle(sums), width = Rle(c(6, 15, 36)))
basic.check$stat <- basic.check$area / basic.check$width

library('bumphunter')
regs2 <- regionFinder(as.numeric(genomeFstats), rep('chr21', length(genomeFstats)), 
    which(prep$position), cluster=NULL, assumeSorted=TRUE, verbose=TRUE, 
    order=FALSE, maxGap=1)

test_that('findRegions', {
    expect_that(width(regs), equals(as.integer(regs2$L)))
    expect_that(width(regs.small[1:2]), equals(runLength(genomeFstats[seq_len(32)] > 1)[runValue(genomeFstats[seq_len(32)] > 1)]))
    expect_that(width(regs.small[3]), equals(runLength(genomeFstats[32 + seq_len(36)] > 1)[runValue(genomeFstats[32 + seq_len(36)] > 1)]))
    expect_that(regs.small.basic, equals(basic.check))
    expect_that(ranges(findRegions(prep$position, genomeFstats, 'chr21', verbose=TRUE, cutoff = 0, maxRegionGap = 1e5, maxClusterGap = 1e6)), equals(IRanges(47407537, 47408970, names = 'up')))
    expect_equal(findRegions(position = FALSE), NULL)
})


test_that('findRegions-smoother', {
    expect_equal(derfinder:::.smootherFstats(genomeFstats, prep$position), Rle(bumphunter::locfitByCluster(genomeFstats, which(prep$position), derfinder:::.clusterMakerRle(prep$position, maxGap = 300L))$fitted))
    expect_equal(derfinder:::.smootherFstats(genomeFstats, prep$position), Rle(smoother(genomeFstats, which(prep$position), derfinder:::.clusterMakerRle(prep$position, maxGap = 300L), smoothFunction = bumphunter::locfitByCluster)$fitted))
    expect_error(findRegions(fstats = genomeFstats, chr = 'chr21', smooth = TRUE))
    expect_warning(findRegions(prep$position, fstats = genomeFstats, chr = 'chr21', smooth = TRUE, basic = TRUE), "Ignoring 'smooth' = TRUE since 'basic' = TRUE")
})

## Find smoothed regions


# Works right now because bumphunter was loaded previously in the session,
# otherwise it would fail due to https://github.com/ririzarr/bumphunter/issues/7
regs_s1 <- findRegions(prep$position, genomeFstats, 'chr21', verbose=TRUE, smooth = TRUE)
regs_s2 <- findRegions(prep$position, genomeFstats, 'chr21', verbose=TRUE, smooth = TRUE, smoothFunction = bumphunter::runmedByCluster)
regs_s3 <- findRegions(prep$position, genomeFstats, 'chr21', verbose=TRUE, smooth = TRUE, minNum = 1435)

test_that('findRegions-smooth-regions', {
    expect_equal(length(regs_s1), 1)
    expect_equal(length(regs_s2), 2)
    expect_equal(regs, regs_s3)
})

## Tests with 2 cores, can't run on Travis
library('BiocParallel')
if('SnowParam' %in% names(registered())) {
    #if(bpworkers(registered()$SnowParam) >= 2) {
    if(FALSE) {
        ## Won't work now due to https://github.com/ririzarr/bumphunter/issues/7
        regs_s4 <- findRegions(prep$position, genomeFstats, 'chr21', verbose=TRUE, smooth = TRUE, mc.cores = 2, maxClusterGap = 1256)
        regs_s5 <- findRegions(prep$position, genomeFstats, 'chr21', verbose=TRUE, smooth = TRUE, mc.cores = 2)
        regs_s6 <- findRegions(prep$position, genomeFstats, 'chr21', verbose=TRUE, smooth = TRUE, smoothFunction = bumphunter::runmedByCluster, mc.cores = 2)
        regs_s7 <- findRegions(prep$position, genomeFstats, 'chr21', verbose=TRUE, smooth = TRUE, mc.cores = 2, maxClusterGap = 1256, minNum = 1435)
        ## The following code use to fail before version 1.5.34
        smoothF <- derfinder:::.smootherFstats(genomeFstats, prep$position, maxClusterGap = 1256, mc.cores = 2)

        test_that('findRegions-smooth-regions-2cores', {
            expect_equal(ranges(regs_s4), IRanges(start = 47409522, width = 15, names = 'up'))
            expect_equal(regs_s1, regs_s5)
            expect_equal(regs_s2, regs_s6)
            expect_equal(ranges(regs), ranges(regs_s7))
            expect_equal(max(smoothF - genomeFstats), 8.5196371283809764208)
            expect_equal(min(smoothF - genomeFstats), -15.547398327000998108)
        })
    }
}
