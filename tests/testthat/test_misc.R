context('Miscellaneous tests')

test_that('Advanced args doc', {
    expect_that((advancedArg('loadCoverage', browse = FALSE)), is_identical_to('https://github.com/search?utf8=%E2%9C%93&q=advanced_argument+filename%3Aload+repo%3Alcolladotor%2Fderfinder+path%3A%2FR&type=Code'))
}
)
