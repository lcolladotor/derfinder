#' Build model matrices for differential expression
#'
#' Builds the model matrices for testing for differential expression by
#' comparing a model with a grouping factor versus one without it. It adjusts
#' for the confounders specified and the median coverage of each sample. The
#' resulting models can be used in [calculateStats].
#'
#' @param sampleDepths Per sample library size adjustments calculated with
#' [sampleDepth].
#' @param testvars A vector or matrix specifying the variables to test. For
#' example, a factor with the group memberships when testing for differences
#' across groups. It's length should match the number of columns used from
#' `coverageInfo$coverage`.
#' @param adjustvars Optional matrix of adjustment variables (e.g. measured
#' confounders, output from SVA, etc.) to use in fitting linear models to each
#' nucleotide. These variables have to be specified by sample and the number of
#' rows must match the number of columns used. It will also work if it is a
#' vector of the correct length.
#' @param testIntercept If `TRUE` then `testvars` is ignored and mod0
#' will contain the column medians and any adjusting variables specified, but
#' no intercept.
#'
#' @return A list with two components.
#' \describe{
#' \item{mod }{ The alternative model matrix.}
#' \item{mod0 }{ The null model matrix.}
#' }
#'
#' @author Leonardo Collado-Torres
#' @seealso [sampleDepth], [calculateStats]
#' @export
#'
#' @importMethodsFrom IRanges ncol median '['
#' @import S4Vectors
#' @examples
#' ## Collapse the coverage information
#' collapsedFull <- collapseFullCoverage(list(genomeData$coverage),
#'     verbose = TRUE
#' )
#'
#' ## Calculate library size adjustments
#' sampleDepths <- sampleDepth(collapsedFull,
#'     probs = c(0.5), nonzero = TRUE,
#'     verbose = TRUE
#' )
#'
#' ## Build the models
#' group <- genomeInfo$pop
#' adjustvars <- data.frame(genomeInfo$gender)
#' models <- makeModels(sampleDepths, testvars = group, adjustvars = adjustvars)
#' names(models)
#' models
makeModels <- function(sampleDepths, testvars, adjustvars = NULL,
    testIntercept = FALSE) {
    ## Drop unused levels in testvars if it is a factor
    if (is.factor(testvars)) {
        testvars <- droplevels(testvars)
    }

    ## Check that the columns match
    numcol <- length(sampleDepths)
    if (!testIntercept & numcol != NROW(testvars)) {
        stop("The length of 'testvars' and the number of sample library size adjustments do not match.")
    }
    if (!is.null(adjustvars) & NROW(adjustvars) != numcol) {
        stop("The dimensions of 'adjustvars' should match with the number of library size adjustments.")
    }

    ## To avoid a warning in R CMD check
    mod <- mod0 <- NULL

    ## Build the adjusted variables if needed
    string1 <- ""
    if (!is.null(adjustvars)) {
        if (NCOL(adjustvars) == 1) {
            ## Only if a vector was supplied
            adjustvars <- as.data.frame(adjustvars)
        }
        for (i in seq_len(NCOL(adjustvars))) {
            eval(parse(text = paste0(
                "adjustVar", i, " <- adjustvars[,",
                i, "]"
            )))
            string1 <- paste(string1, paste0("adjustVar", i),
                sep = "+"
            )
        }
    }
    if (!testIntercept) {
        eval(parse(text = paste0(
            "mod = model.matrix(~ testvars + sampleDepths",
            string1, ")"
        )))
        eval(parse(text = paste0(
            "mod0 = model.matrix(~ + sampleDepths",
            string1, ")"
        )))
    } else {
        eval(parse(text = paste0(
            "mod = model.matrix(~ sampleDepths",
            string1, ")"
        )))
        eval(parse(text = paste0(
            "mod0 = model.matrix(~ 0 + sampleDepths",
            string1, ")"
        )))
    }

    ## Check that the matrices are full rank
    if (qr(mod)$rank != ncol(mod)) {
        r <- qr(mod)$rank
        drop.mod.col <- colnames(mod)[(r + 1):ncol(mod)]
        warning(paste(
            "Dropping from the alternative model matrix (mod) the column(s)",
            paste(drop.mod.col, collapse = ", "),
            "as the matrix is not full rank."
        ))
        mod <- mod[, seq_len(r), drop = FALSE]
        stopifnot(ncol(mod) > 0)

        ## Drop the same columns from mod0 if present in mod0
        if (any(colnames(mod0) %in% drop.mod.col)) {
            drop.mod0.col <- colnames(mod0)[colnames(mod0) %in%
                drop.mod.col]
            mod0 <- mod0[, !colnames(mod0) %in% drop.mod0.col,
                drop = FALSE
            ]
            warning(paste(
                "Dropping from the null model matrix (mod0) the column(s)",
                paste(drop.mod0.col, collapse = ", "),
                "as they were dropped in the alternative model matrix (mod)."
            ))
        }
        stopifnot(ncol(mod0) > 0)
    }
    ## This should not happen, since mod0 is a subset of mod and
    ## mod has been truncated to be full rank.
    stopifnot(qr(mod0)$rank == ncol(mod0))

    ## Finish
    result <- list(mod = mod, mod0 = mod0)

    ## Done =)
    return(result)
}
