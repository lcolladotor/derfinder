## Based on https://github.com/hadley/ggplot2/blob/master/R/zzz.r
.onAttach <- function(...) {
    if (!interactive() || stats::runif(1) > 0.99) return()
        
    tips <- c(
        "Found a bug? Report it at https://github.com/lcolladotor/derfinder/issues",
        "Want to contribute a new feature? Fork derfinder at\n https://github.com/lcolladotor/derfinder/fork\nThen submit a pull request =)",
        paste("Find out what's changed in derfinder with\n",
            "news(Version == \"", utils::packageVersion("derfinder"),
            "\", package = \"derfinder\")", sep = ""),
        "Create HTML reports from derfinder results using regionReport available at\nhttp://www.bioconductor.org/packages/devel/bioc/html/regionReport.html",
        "Use suppressPackageStartupMessages to eliminate package startup messages."
    )
    
    tip <- sample(tips, 1, prob = c(0.05, 0.05, 0.55, 0.2, 0.15))
    packageStartupMessage(tip)
}
