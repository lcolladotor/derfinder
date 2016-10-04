context('Miscellaneous tests')

test_that('Advanced args doc', {
    expect_that((advancedArg('loadCoverage', browse = FALSE)), is_identical_to('https://github.com/search?utf8=%E2%9C%93&q=advanced_argument+filename%3Aload+repo%3Alcolladotor%2Fderfinder+path%3A%2FR&type=Code'))
}
)

test_that('Mapping levels', {
    expect_that(extendedMapSeqlevels('seq2', 'UCSC', verbose = FALSE), equals('seq2'))
    expect_that(extendedMapSeqlevels('seq2', 'UCSC', NULL, verbose = TRUE), shows_message("extendedMapSeqlevels: the 'seqnames' you supplied are currently not supported in GenomeInfoDb. Consider adding your genome by following the information at http://www.bioconductor.org/packages/release/bioc/vignettes/GenomeInfoDb/inst/doc/Accept-organism-for-GenomeInfoDb.pdf"))
    expect_that(extendedMapSeqlevels('seq2', 'UCSC', 'toy', verbose = FALSE), equals('seq2'))
    expect_that(extendedMapSeqlevels('seq2', 'UCSC', 'toy', verbose = TRUE), shows_message())
    expect_that(extendedMapSeqlevels('2', 'toyStyle', 'Arabidopsis thaliana', verbose = FALSE), equals('2'))
    expect_that(extendedMapSeqlevels('2', 'toyStyle', 'Arabidopsis thaliana', verbose = TRUE), shows_message())
    expect_that(extendedMapSeqlevels('2', 'NCBI', 'Arabidopsis thaliana', 'toyStyle', verbose = FALSE), equals('2'))
    expect_that(extendedMapSeqlevels('2', 'NCBI', 'Arabidopsis thaliana', 'toyStyle', verbose = TRUE), shows_message(""))
    expect_that(extendedMapSeqlevels('2', 'UCSC', 'Homo sapiens', 'NCBI'), equals('chr2'))
    expect_that(extendedMapSeqlevels('2', 'UCSC', 'Homo sapiens', verbose = TRUE), shows_message("extendedMapSeqlevels: sequence names mapped from NCBI to UCSC for species homo_sapiens"))
    expect_that(extendedMapSeqlevels('2', 'UCSC', 'Homo sapiens', 'NCBI'), equals(extendedMapSeqlevels('2', 'UCSC', 'Homo sapiens', verbose = FALSE)))
    expect_that(extendedMapSeqlevels('2', 'UCSC', 'Homo sapiens', 'NCBI'), equals(extendedMapSeqlevels('2', 'UCSC', 'homo_saPiens', 'NcBI')))
    expect_that(extendedMapSeqlevels('28', NULL), equals('28'))
    expect_that(extendedMapSeqlevels('T', 'NCBI', 'homo_sapiens'), shows_message())
    expect_that(extendedMapSeqlevels('T', 'NCBI', 'homo_sapiens', verbose = FALSE), equals('T'))
    expect_that(extendedMapSeqlevels('T', 'NCBI', NULL, verbose = FALSE), equals('LGI'))
})


foo <- function(...) {
    derfinder::define_cluster(...)
}
fuu <- function(...) {
    cores <- derfinder:::.advanced_argument('cores', 1L, ...)
    derfinder::define_cluster(cores = 'test', test = cores)
}

library('BiocParallel')
test_that('Cluster setup', {
    expect_that(foo(), equals(SerialParam()))
    expect_that(foo(mc.cores = 2), equals(SnowParam(2, outfile = Sys.getenv('SGE_STDERR_PATH'))))
    expect_that(foo(), equals(fuu()))
    expect_that(fuu(cores = 4), equals(SnowParam(4, outfile = Sys.getenv('SGE_STDERR_PATH'))))
})
