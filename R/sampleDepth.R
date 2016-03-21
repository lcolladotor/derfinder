#' Calculate adjustments for library size
#'
#' For a given data set calculate the per-sample coverage adjustments. Hector 
#' Corrada's group proposed calculating the sum of the coverage for genes below 
#' a given sample quantile. In this function, we calculate the sample quantiles 
#' of interest by sample, and then the sum of the coverage for bases below or 
#' equal to quantiles of interest. The resulting values are transformed {log2(x 
#' + scalefac)} to avoid very large numbers that could potentially affect the 
#' stability of the F-statistics calculation. The sample coverage adjustments 
#' are then used in \link{makeModels} for constructing the null and alternative 
#' models.
#' 
#' @param collapsedFull The full coverage data collapsed by sample as produced 
#' by \link{collapseFullCoverage}.
#' @param probs Number(s) between 0 and 1 representing the quantile(s) of 
#' interest. For example, 0.5 is the median.
#' @param scalefac Number added to the sample coverage adjustments before the 
#' log2 transformation.
#' @param ... Arguments passed to other methods and/or advanced arguments.
#'
#' @return 
#' A matrix (vector of \code{length(probs) == 1}) with the library size depth 
#' adjustments per sample to be used in \link{makeModels}. The number of rows 
#' corresponds to the number of quantiles used for the sample adjustments. 
#'
#' @references
#' Paulson, J. N., Stine, O. C., Bravo, H. C. & Pop, M. Differential abundance 
#' analysis for microbial marker-gene surveys. Nat. Methods (2013). 
#' doi:10.1038/nmeth.2658
#'
#' @author Leonardo Collado-Torres
#' @seealso \link{collapseFullCoverage}, \link{makeModels}
#' @export
#' @importFrom Hmisc wtd.quantile
#' @examples
#' ## Collapse the coverage information
#' collapsedFull <- collapseFullCoverage(list(genomeData$coverage), 
#'     verbose=TRUE)
#' 
#' ## Calculate library size adjustments
#' sampleDepths <- sampleDepth(collapsedFull, probs=c(0.5, 1), verbose=TRUE)
#' sampleDepths

sampleDepth <- function(collapsedFull, probs = c(0.5, 1), scalefac = 32, ...) {
    
    ## Check probs are valid
    stopifnot(all(probs >= 0) & all(probs <= 1))

    ## Advanged arguments
# @param verbose If \code{TRUE} basic status updates will be printed along the 
# way.
    verbose <- .advanced_argument('verbose', TRUE, ...)


# @param nonzero If \code{TRUE} only the nonzero counts are used to calculate 
# the library size adjustment.
    nonzero <- .advanced_argument('nonzero', TRUE, ...)


# @param center If \code{TRUE} the sample coverage adjustements are centered. 
# In some cases, this could be helpful for interpretation purposes.
    center <- .advanced_argument('center', FALSE, ...)


    if (verbose) 
        message(paste(Sys.time(), 'sampleDepth: Calculating sample quantiles'))
    sampleQuant <- lapply(collapsedFull, .calcQuantile, nonzero = nonzero, 
        probs = probs)
    
    if (verbose) 
        message(paste(Sys.time(),
            'sampleDepth: Calculating sample adjustments'))
    sampleDepths <- mapply(.sampleCorradaAdj, collapsedFull, 
        sampleQuant, MoreArgs = list(scalefac = scalefac, center = center))
    
    ## Done =)
    return(sampleDepths)
}

## Function for calculating the quantiles
.calcQuantile <- function(x, nonzero, probs) {
    values <- x$values
    weights <- x$weights
    if (nonzero) {
        idx <- which(values == 0)
        values <- values[-idx]
        weights <- weights[-idx]
    }
    wtd.quantile(values, weights = weights, probs = probs)
}

## Function for calculating the Hector Corrada-type sample
## depth adjustments
.sampleCorradaAdj <- function(collapsedF, quants, scalefac, center) {
    sapply(quants, function(z) {
        idx <- which(collapsedF$values <= z)
        adj <- sum(collapsedF$weights[idx] * collapsedF$values[idx])
        res <- log2(adj + scalefac)
        if (center) {
            res <- res - mean(res)
        }
        return(res)
    })
}
