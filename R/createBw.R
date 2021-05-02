#' Export coverage to BigWig files
#'
#' Using output from [fullCoverage], export the coverage from all the
#' samples to BigWig files using [createBwSample].
#'
#' @param fullCov A list where each element is the result from
#' [loadCoverage] used with `returnCoverage = TRUE`. Can be generated
#' using [fullCoverage].
#' @param path The path where the BigWig files will be created.
#' @param keepGR If `TRUE`, the [GRanges][GenomicRanges::GRanges-class] objects
#' created by [coerceGR] grouped into a [GRangesList][GenomicRanges::GRangesList-class]
#' are returned. Otherwise they are discarded.
#' @param ... Arguments passed to [createBwSample].
#'
#' @return If `keepGR = TRUE`, then a [GRangesList][GenomicRanges::GRangesList-class]
#' with the output for [coerceGR] for each of the samples.
#'
#' @details Use at most one core per chromosome.
#'
#' @author Leonardo Collado-Torres
#' @seealso [GRangesList][GenomicRanges::GRangesList-class],
#' [export.bw][rtracklayer::BigWigFile-class],
#' [createBwSample], [coerceGR]
#' @export
#'
#' @importMethodsFrom GenomicRanges names
#' @importFrom GenomicRanges GRangesList
#'
#' @examples
#' ## Create a small fullCov object with data only for chr21
#' fullCov <- list("chr21" = genomeDataRaw)
#'
#' ## Keep only 2 samples
#' fullCov$chr21$coverage <- fullCov$chr21$coverage[c(1, 31)]
#'
#' ## Create the BigWig files for all samples in a test dir
#' dir.create("createBw-example")
#' bws <- createBw(fullCov, "createBw-example")
#'
#' ## Explore the output
#' bws
#'
#' ## First sample
#' bws[[1]]
#'
#' ## Note that if a sample has no bases with coverage > 0, the GRanges object
#' ## is empty and no BigWig file is created for that sample.
#' bws[[2]]
#'
#'
#' ## Exports fullCoverage() output to BigWig files
createBw <- function(fullCov, path = ".", keepGR = TRUE, ...) {

    ## Determine sample names
    samples <- names(fullCov[[1]])
    if ("coverage" %in% samples) {
        samples <- names(fullCov[[1]]$coverage)
    }

    ## Coerce to GR
    gr.samples <- lapply(samples, createBwSample,
        path = path,
        fullCov = fullCov, keepGR = keepGR, ...
    )

    ## Done
    if (keepGR) {
        gr.samples <- GRangesList(gr.samples)
        names(gr.samples) <- samples
        return(invisible(gr.samples))
    } else {
        return(invisible(NULL))
    }
}
