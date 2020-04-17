context("Annotation-related functions")


library("TxDb.Hsapiens.UCSC.hg19.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

GenomicState.Hsapiens.UCSC.hg19.knownGene.chr21 <- makeGenomicState(
    txdb = txdb,
    chrs = "chr21"
)

library("GenomicFeatures")
samplefile <- system.file("extdata", "hg19_knownGene_sample.sqlite",
    package = "GenomicFeatures"
)
txdb <- loadDb(samplefile)

## Generate genomic state object, only for chr6
sampleGenomicState <- makeGenomicState(txdb, chrs = "chr6")

test_that("makeGenomicState", {
    expect_that(length(sampleGenomicState$fullGenome), equals(49))
    expect_that(as.vector(table(sampleGenomicState$fullGenome$theRegion)), equals(c(24, 6, 19)))
    expect_that(GenomicState.Hsapiens.UCSC.hg19.knownGene.chr21, is_equivalent_to(genomicState))
})


copy <- GRanges(seqnames = "chr21", ranges = ranges(genomicState$fullGenome[1:2]), strand = "+")
annoRegs <- annotateRegions(regions = copy, genomicState = genomicState$fullGenome, minoverlap = 1)

test_that("annotateRegions", {
    expect_that(annoRegs$annotationList[[1]], is_equivalent_to(genomicState$fullGenome[1]))
    expect_that(annoRegs$countTable, is_equivalent_to(data.frame(exon = c(1, 1), intragenic = c(0, 0), intron = c(0, 0))))
    expect_that(annotateRegions(regions = copy, genomicState = genomicState$fullGenome, minoverlap = min(width(copy)), annotate = FALSE, verbose = FALSE), is_equivalent_to(list(count.table = data.frame(exon = c(1, 1), intragenic = c(0, 0), intron = c(0, 0)))))
    expect_that(annotateRegions(regions = copy, genomicState = genomicState$fullGenome, minoverlap = max(width(copy)) + 1, annotate = FALSE, verbose = FALSE), is_equivalent_to(list(count.table = data.frame(exon = c(0, 0), intragenic = c(0, 0), intron = c(0, 0)))))
})


## Strand info
regs <- genomeRegions$regions[1:2]
ann <- annotateRegions(regions = regs, genomicState = genomicState$fullGenome, minoverlap = 1)
library("GenomicRanges")
regs2 <- regs
strand(regs2) <- "-"
ann2 <- annotateRegions(regions = regs2, genomicState = genomicState$fullGenome, minoverlap = 1)
test_that("ignore strand", {
    expect_equal(ann, ann2)
})
