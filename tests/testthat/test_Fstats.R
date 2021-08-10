library("GenomicRanges")

## Check that the sample depths are correctly calculated
collapsedFull <- collapseFullCoverage(list(genomeData$coverage),
    verbose = TRUE
)

sampleDepths <- sampleDepth(collapsedFull, probs = c(0.5), verbose = TRUE)

sample.manual <- sapply(collapsedFull, function(x) {
    log2(sum(x$values * x$weights) + 32)
})
names(sample.manual) <- paste0(names(sample.manual), ".100%")


test_that("sampleDepths", {
    expect_that(lapply(genomeData$coverage[1:2], function(x) {
        vals <- unique(runValue(x))
        weights <- sapply(vals, function(y) {
            sum(x == y)
        })
        return(list(values = vals, weights = weights))
    }), equals(collapseFullCoverage(list(genomeData$coverage), colsubset = 1:2)))
    expect_that(sampleDepth(collapsedFull, probs = c(1), verbose = FALSE), equals(sample.manual))
})


## Check that the models are properly being built
group <- genomeInfo$pop
adjustvars <- data.frame(genomeInfo$gender)
models <- makeModels(sampleDepths, testvars = group, adjustvars = adjustvars)

testvars <- group
adjustVar1 <- adjustvars[[1]]
models.manual <- list(mod = model.matrix(~ testvars + sampleDepths + adjustVar1), mod0 = model.matrix(~ sampleDepths + adjustVar1))

test_that("makeModels", {
    expect_that(models, equals(models.manual))
    expect_that(makeModels(sampleDepths, testvars = group, adjustvars = adjustvars, testIntercept = TRUE), equals(list(mod = model.matrix(~ sampleDepths + adjustVar1), mod0 = model.matrix(~ 0 + sampleDepths + adjustVar1))))
})


prep <- preprocessCoverage(genomeData,
    cutoff = 0, scalefac = 32,
    chunksize = 1e3, colsubset = NULL
)

## Manually calculate the F-stats
fstats <- calculateStats(prep, models, verbose = TRUE, method = "regular")
prep$coverageProcessed

n <- ncol(prep$coverageProcessed)
Id <- diag(n)
df1 <- ncol(models$mod)
df0 <- ncol(models$mod0)
P1 <- (Id - models$mod %*% solve(t(models$mod) %*% models$mod) %*% t(models$mod))
P0 <- (Id - models$mod0 %*% solve(t(models$mod0) %*% models$mod0) %*% t(models$mod0))
resid0 <- as.matrix(as.data.frame(prep$coverageProcessed)) %*% P0
resid1 <- as.matrix(as.data.frame(prep$coverageProcessed)) %*% P1
rss0 <- (resid0 * resid0) %*% rep(1L, n)
rss1 <- (resid1 * resid1) %*% rep(1L, n)

fstats.manual <- Rle(drop(((rss0 - rss1) / (df1 - df0)) / (rss1 / (n - df1))))

test_that("calculateFstats", {
    expect_equivalent(fstats, fstats.manual)
    expect_equivalent(fstats, genomeFstats)
})

## Check that the F-stat cutoff is properly constructed
test_that("F-stat cutoff", {
    expect_that(derfinder:::.calcFstatCutoff("empirical", 0.5, Rle(1:2, rep(10, 2))), equals(c("50%" = 1.5)))
    expect_that(derfinder:::.calcFstatCutoff("theoretical", 0.05, models = models), equals(qf(0.05, 4 - 3, 31 - 4, lower.tail = FALSE)))
    expect_that(derfinder:::.calcFstatCutoff("manual", 1), equals(1))
})


## Check calculating the p-values. Specially check that the default arguments
## are the ones intended
suppressWarnings(RNGversion("3.5.0"))
regsWithP <- calculatePvalues(prep, models, fstats,
    nPermute = 10, seeds = 1:10,
    chr = "chr21", cutoff = 1, mc.cores = 1, method = "regular"
)

test_that("calculatePvalues", {
    expect_that(regsWithP, is_equivalent_to(genomeRegions))
    expect_that(calculatePvalues(prep, models, fstats, nPermute = 1, seeds = 11, chr = "chr21", cutoff = 1, mc.cores = 1, method = "regular"), equals(calculatePvalues(prep, models, fstats, nPermute = 1, seeds = 11, chr = "chr21", cutoff = 1, mc.cores = 1, method = "regular", maxRegionGap = 0L, maxClusterGap = 300L, verbose = TRUE, scalefac = 32, adjustF = 0)))
    expect_that(derfinder:::.calcPval(1:2, 1:4), equals(c(0.8, 0.6)))
    expect_that(derfinder:::.calculateFWER(10, 1:4, rep(c(1, 3), each = 2), 3, 1:2), equals(c(1, 0.75)))
})


## Check merging results
dir.create("generateReport-example", showWarnings = FALSE, recursive = TRUE)
file.copy(system.file(file.path("extdata", "chr21"),
    package = "derfinder",
    mustWork = TRUE
), "generateReport-example", recursive = TRUE)
mergeResults(
    chrs = "21", prefix = "generateReport-example",
    genomicState = genomicState$fullGenome
)
load(file.path("generateReport-example", "fullRegions.Rdata"))
fullRegions.included <- fullRegions

## Make a fresh copy
initialPath <- getwd()
dir.create("generateReport-example-rerun",
    showWarnings = FALSE,
    recursive = TRUE
)
setwd(file.path(initialPath, "generateReport-example-rerun"))
collapsedFull <- collapseFullCoverage(list(genomeData$coverage),
    verbose = TRUE
)
sampleDepths <- sampleDepth(collapsedFull,
    probs = c(0.5), nonzero = TRUE,
    verbose = TRUE
)
groupInfo <- genomeInfo$pop
adjustvars <- data.frame(genomeInfo$gender)
models <- makeModels(sampleDepths, testvars = groupInfo, adjustvars = adjustvars)
analyzeChr(
    chr = "21", coverageInfo = genomeData, models = models,
    cutoffFstat = 1, cutoffType = "manual", seeds = 20140330, groupInfo = groupInfo,
    mc.cores = 1, writeOutput = TRUE, returnOutput = FALSE
)

## Merge fresh results
setwd(initialPath)
mergeResults(
    chrs = "21", prefix = "generateReport-example-rerun",
    genomicState = genomicState$fullGenome
)

## Load results
load(file.path("generateReport-example-rerun", "fullRegions.Rdata"))
fullRegions.rerun <- fullRegions

test_that("mergeResults", {
    expect_that(ranges(fullRegions.included), is_equivalent_to(ranges(genomeRegions$regions)))
    expect_that(fullRegions.included, is_equivalent_to(fullRegions.rerun))
    expect_that(dir("generateReport-example"), equals(dir("generateReport-example-rerun")))
})


library("TxDb.Hsapiens.UCSC.hg19.knownGene")

results <- analyzeChr(
    chr = "21", coverageInfo = genomeData, models = models,
    cutoffFstat = 1, cutoffType = "manual", groupInfo = groupInfo, mc.cores = 1,
    writeOutput = FALSE, returnOutput = TRUE, method = "regular"
)
results.param <- analyzeChr(
    chr = "21", coverageInfo = genomeData, models = models,
    cutoffFstat = 1, cutoffType = "manual", groupInfo = groupInfo, mc.cores = 1,
    writeOutput = FALSE, returnOutput = TRUE, method = "regular", maxRegionGap = 0L,
    maxClusterGap = 300L, txdb = TxDb.Hsapiens.UCSC.hg19.knownGene, scalefac = 32
)

test_that("analyzeChr", {
    expect_that(results$regions, equals(results.param$regions))
    expect_that(results$annotation, equals(results.param$annotation))
    expect_that(results$coveragePrep, equals(results.param$coveragePrep))
    expect_that(c(results$optionsStats[-which(names(results$optionsStats) == "analyzeCall")], list(maxRegionGap = 0L, maxClusterGap = 300L, scalefac = 32)), is_equivalent_to(results.param$optionsStats[-which(names(results.param$optionsStats) == "analyzeCall")]))
})


## Clean up
unlink("generateReport-example", recursive = TRUE)
unlink("generateReport-example-rerun", recursive = TRUE)
unlink("chr21", recursive = TRUE)

