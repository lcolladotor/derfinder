#' Merge results from different chromosomes
#'
#' This function merges the results from running [analyzeChr] on several
#' chromosomes and assigns genomic states using [annotateRegions]. It
#' re-calculates the p-values and q-values using the pooled areas from the null
#' regions from all chromosomes. Once the results have been merged,
#' `derfinderReport::generateReport` can be used to generate an HTML
#' report of the results. The `derfinderReport` package is available at
#' https://github.com/lcolladotor/derfinderReport.
#'
#' @param chrs The chromosomes of the files to be merged.
#' @param prefix The main data directory path, which can be useful if
#' [analyzeChr] is used for several parameters and the results are saved
#' in different directories.
#' @param significantCut A vector of length two specifiying the cutoffs used to
#' determine significance. The first element is used to determine significance
#' for the P-values and FWER adjusted P-values, while the second element is
#' used for the Q-values (FDR adjusted P-values) similar to
#' [calculatePvalues].
#' @param minoverlap Determines the mininum overlap needed when annotating
#' regions with [annotateRegions].
#' @param mergePrep If `TRUE` the output from [preprocessCoverage] is
#' merged.
#' @param ... Arguments passed to other methods and/or advanced arguments.
#' Advanced arguments:
#' \describe{
#' \item{verbose }{ If `TRUE` basic status updates will be printed along
#' the way.}
#' \item{optionsStats }{ The options used in [analyzeChr]. By default
#' `NULL` and will be inferred from the output files.}
#' \item{cutoffFstatUsed }{ The actual F-statistic cutoff used. This can be
#' obtained from the logs or from the output of [analyzeChr]. If
#' `NULL` then this function will attempt to re-calculate it.}
#' }
#' Passed to [annotateRegions] and [extendedMapSeqlevels].
#' @inheritParams annotateRegions
#'
#'
#' @return Seven Rdata files.
#' \describe{
#' \item{fullFstats.Rdata }{ Full F-statistics from all chromosomes in a list
#' of Rle objects.}
#' \item{fullTime.Rdata }{ Timing information from all chromosomes.}
#' \item{fullNullSummary.Rdata}{ A DataFrame with the null region information:
#' statistic, width, chromosome and permutation identifier. It's ordered by the
#' statistics}
#' \item{fullRegions.Rdata}{ GRanges object with regions found and with full
#' annotation from [matchGenes][bumphunter::matchGenes]. Note that the column
#' `strand` from [matchGenes][bumphunter::matchGenes] is renamed to
#' `annoStrand` to comply with GRanges specifications. }
#' \item{fullCoveragePrep.Rdata}{ A list with the pre-processed coverage data
#' from all chromosomes.}
#' \item{fullAnnotatedRegions.Rdata}{ A list as constructed in
#' [annotateRegions] with the assigned genomic states.}
#' \item{optionsMerge.Rdata}{ A list with the options used when merging the
#' results. Used in `derfinderReport::generateReport`.}
#' }
#'
#' @author Leonardo Collado-Torres
#' @seealso [analyzeChr], [calculatePvalues], [annotateRegions]
#' @export
#'
#' @importFrom GenomicRanges GRangesList
#' @importMethodsFrom GenomicRanges '$' '$<-' '['
#' @importFrom IRanges RleList
#' @import S4Vectors
#' @importMethodsFrom IRanges cbind values 'values<-' '[' length
#' order unlist nrow
#' @importFrom qvalue qvalue
#'
#' @details If you want to calculate the FWER adjusted P-values, supply
#' `optionsStats` which is produced by [analyzeChr].
#'
#' @examples
#' ## The output will be saved in the 'generateReport-example' directory
#' dir.create("generateReport-example", showWarnings = FALSE, recursive = TRUE)
#'
#' ## For convenience, the derfinder output has been pre-computed
#' file.copy(system.file(file.path("extdata", "chr21"),
#'     package = "derfinder",
#'     mustWork = TRUE
#' ), "generateReport-example", recursive = TRUE)
#'
#' ## Merge the results from the different chromosomes. In this case, there's
#' ## only one: chr21
#' mergeResults(
#'     chrs = "21", prefix = "generateReport-example",
#'     genomicState = genomicState$fullGenome
#' )
#' \dontrun{
#' ## You can then explore the wallclock time spent on each step
#' load(file.path("generateReport-example", "fullRegions.Rdata"))
#'
#' ## Process the time info
#' time <- lapply(fullTime, function(x) data.frame(diff(x)))
#' time <- do.call(rbind, time)
#' colnames(time) <- "sec"
#' time$sec <- as.integer(round(time$sec))
#' time$min <- time$sec / 60
#' time$chr <- paste0("chr", gsub("\\..*", "", rownames(time)))
#' time$step <- gsub(".*\\.", "", rownames(time))
#' rownames(time) <- seq_len(nrow(time))
#'
#' ## Make plot
#' library("ggplot2")
#' ggplot(time, aes(x = step, y = min, colour = chr)) +
#'     geom_point() +
#'     labs(title = "Wallclock time by step") +
#'     scale_colour_discrete(limits = chrs) +
#'     scale_x_discrete(limits = names(fullTime[[1]])[-1]) +
#'     ylab("Time (min)") +
#'     xlab("Step")
#' }
#'
mergeResults <- function(chrs = c(seq_len(22), "X", "Y"), prefix = ".",
    significantCut = c(0.05, 0.1), genomicState, minoverlap = 20,
    mergePrep = FALSE, ...) {
    ## For R CMD check
    prep <- fstats <- regions <- annotation <- timeinfo <- NULL

    stopifnot(length(significantCut) == 2 & all(significantCut >=
        0 & significantCut <= 1))

    ## Advanged argumentsa
    # @param verbose If \code{TRUE} basic status updates will be printed along the
    # way.
    verbose <- .advanced_argument("verbose", TRUE, ...)

    # @param optionsStats The options used in \link{analyzeChr}.
    optionsStats <- .advanced_argument("optionsStats", NULL, ...)

    # @param cutoffFstatUsed The actual F-statistic cutoff used. This can be
    # obtained from the logs or from the output of \link{analyzeChr}. If
    # \code{NULL} then this function will attempt to re-calculate it.
    cutoffFstatUsed <- .advanced_argument("cutoffFstatUsed", optionsStats$cutoffFstatUsed, ...)



    ## Use UCSC names for homo_sapiens by default
    chrs <- extendedMapSeqlevels(chrs, ...)

    ## save merging options used
    optionsMerge <- list(
        chrs = chrs, significantCut = significantCut,
        minoverlap = minoverlap, mergeCall = match.call(),
        cutoffFstatUsed = cutoffFstatUsed, ...
    )
    if (verbose) {
          message(paste(Sys.time(), "mergeResults: Saving options used"))
      }
    save(optionsMerge, file = file.path(prefix, "optionsMerge.Rdata"))

    ## Initialize
    fullCoveragePrep <- fullTime <- fullNullPermutation <- fullNullWidths <-
        fullNullStats <- fullFstats <- fullAnno <- fullRegs <- vector(
            "list",
            length(chrs)
        )
    names(fullCoveragePrep) <- names(fullTime) <- names(fullNullPermutation) <-
        names(fullNullWidths) <- names(fullNullStats) <- names(fullFstats) <-
        names(fullAnno) <- names(fullRegs) <- chrs

    ## Actual processing
    for (chr in chrs) {
        if (verbose) {
              message(paste(Sys.time(), "Loading chromosome", chr))
          }

        ## Process the F-statistics
        load(file.path(prefix, chr, "fstats.Rdata"))
        fullFstats[[chr]] <- fstats

        ## Process the regions, nullstats and nullwidths
        load(file.path(prefix, chr, "regions.Rdata"))
        fullRegs[[chr]] <- regions$regions
        fullNullStats[[chr]] <- regions$nullStats
        fullNullWidths[[chr]] <- regions$nullWidths
        fullNullPermutation[[chr]] <- regions$nullPermutation

        ## Process the annotation results
        load(file.path(prefix, chr, "annotation.Rdata"))
        fullAnno[[chr]] <- annotation

        ## Process the timing information
        load(file.path(prefix, chr, "timeinfo.Rdata"))
        fullTime[[chr]] <- timeinfo

        ## Process the covPrep data
        if (mergePrep) {
            load(file.path(prefix, chr, "coveragePrep.Rdata"))
            fullCoveragePrep[[chr]] <- prep
        }
    }

    ## Merge regions
    fullRegions <- unlist(GRangesList(fullRegs), use.names = FALSE)

    ## Process the annotation
    fullAnnotation <- do.call(rbind, fullAnno)
    if (!is.null(fullAnnotation)) {
        colnames(fullAnnotation)[which(colnames(fullAnnotation) ==
            "strand")] <- "annoStrand"
        rownames(fullAnnotation) <- NULL

        ## For some reason, signature 'AsIs' does not work when
        ## assigning the values() <-
        fullAnnotation$name <- as.character(fullAnnotation$name)
        fullAnnotation$annotation <- as.character(fullAnnotation$annotation)

        ## Combine regions with annotation
        values(fullRegions) <- cbind(
            values(fullRegions),
            DataFrame(fullAnnotation, check.names = FALSE)
        )
    }


    ## Summarize the null regions
    nulls <- unlist(RleList(fullNullStats), use.names = FALSE)
    widths <- unlist(RleList(fullNullWidths), use.names = FALSE)
    permutations <- unlist(RleList(fullNullPermutation), use.names = FALSE)
    howMany <- unlist(lapply(fullNullStats, length))

    if (length(nulls) > 0) {
        ## Proceed only if there are null regions to work with
        fullNullSummary <- DataFrame(
            stat = nulls, width = widths,
            chr = Rle(names(fullNullStats), howMany),
            permutation = permutations, check.names = FALSE
        )
        rm(nulls, widths, howMany, permutations)

        fullNullSummary$area <- abs(fullNullSummary$stat) *
            fullNullSummary$width
        fullNullSummary <- fullNullSummary[order(fullNullSummary$area,
            decreasing = TRUE
        ), ]
    } else {
        fullNullSummary <- DataFrame(NULL)
    }

    ## Determine F-stat cutoff if not supplied
    if (is.null(cutoffFstatUsed)) {
        if (!is.null(optionsStats)) {
            stopifnot(all(c("cutoffType", "cutoffFstat", "models") %in% names(optionsStats)))
            if (verbose) {
                  message(paste(
                      Sys.time(),
                      "mergeResults: calculating F-stat cutoff"
                  ))
              }
            cutoffFstatUsed <- .calcFstatCutoff(
                cutoffType = optionsStats$cutoffType,
                cutoffFstat = optionsStats$cutoffFstat,
                fstats = unlist(fullFstats),
                models = optionsStats$models
            )
        } else {
            if (verbose) {
                  message("Neither 'cutoffFstatUsed' nor 'optionsStats' were supplied, so the FWER calculation step will be skipped.")
              }
        }
    }
    if (!is.null(cutoffFstatUsed)) {
        if (nrow(fullNullSummary) > 0) {
            if (verbose) {
                  message(paste(
                      Sys.time(),
                      "mergeResults: calculating FWER"
                  ))
              }
            stopifnot(!is.null(optionsStats))
            stopifnot("nPermute" %in% names(optionsStats))
            fullRegions$fwer <- .calculateFWER(
                cutoffFstatUsed = cutoffFstatUsed,
                areaNull = fullNullSummary$area,
                permutation = fullNullSummary$permutation,
                nPermute = optionsStats$nPermute,
                areaReg = fullRegions$area
            )
            fullRegions$significantFWER <- factor(fullRegions$fwer < significantCut[1],
                levels = c(TRUE, FALSE)
            )
        } else {
            if (verbose) {
                  warning("No null regions were found, so the FWER calculation step will be skipped.")
              }
        }
    }

    if (verbose) {
          message(paste(Sys.time(), "mergeResults: Saving fullNullSummary"))
      }
    save(fullNullSummary, file = file.path(prefix, "fullNullSummary.Rdata"))

    ## Re-calculate p-values and q-values
    if (verbose) {
          message(paste(Sys.time(), "mergeResults: Re-calculating the p-values"))
      }

    ## Sort by decreasing area
    fullRegions <- fullRegions[order(fullRegions$area, decreasing = TRUE)]

    if (nrow(fullNullSummary) > 0) {
        ## Actual calculation
        fullRegions$pvalues <- .calcPval(
            fullRegions$area,
            as.numeric(fullNullSummary$area)
        )

        ## Update info
        fullRegions$significant <- factor(fullRegions$pvalues <
            significantCut[1], levels = c(TRUE, FALSE))

        ## Sometimes qvalue() fails due to incorrect pi0 estimates
        qvalues <- tryCatch(qvalue(fullRegions$pvalues), error = function(e) NULL)
        if (!is.null(qvalues)) {
            qvalues <- qvalues$qvalues
            sigQval <- factor(qvalues < significantCut[2], levels = c(
                TRUE,
                FALSE
            ))
        } else {
            message(paste(Sys.time(), "mergeResults: skipping q-value calculation."))
            qvalues <- rep(NA, length(fullRegions$pvalues))
            sigQval <- rep(NA, length(fullRegions$pvalues))
        }
        fullRegions$qvalues <- qvalues
        fullRegions$significantQval <- sigQval
    }

    ## save GRanges version
    if (verbose) {
          message(paste(Sys.time(), "mergeResults: Saving fullRegions"))
      }
    save(fullRegions, file = file.path(prefix, "fullRegions.Rdata"))

    ## Assign genomic states
    if (verbose) {
          message(paste(Sys.time(), "mergeResults: assigning genomic states"))
      }
    fullAnnotatedRegions <- annotateRegions(
        regions = fullRegions,
        genomicState = genomicState, minoverlap = minoverlap,
        annotate = TRUE, ...
    )

    if (verbose) {
          message(paste(Sys.time(), "mergeResults: Saving fullAnnotatedRegions"))
      }
    save(fullAnnotatedRegions, file = file.path(
        prefix,
        "fullAnnotatedRegions.Rdata"
    ))

    ## Save Fstats, Nullstats, and time info
    if (verbose) {
          message(paste(Sys.time(), "mergeResults: Saving fullFstats"))
      }
    save(fullFstats, file = file.path(prefix, "fullFstats.Rdata"))

    if (verbose) {
          message(paste(Sys.time(), "mergeResults: Saving fullTime"))
      }
    save(fullTime, file = file.path(prefix, "fullTime.Rdata"))

    if (mergePrep) {
        if (verbose) {
              message(paste(Sys.time(), "mergeResults: Saving fullCoveragePrep"))
          }
        save(fullCoveragePrep, file = file.path(
            prefix,
            "fullCoveragePrep.Rdata"
        ))
    }

    ## Finish
    return(invisible(NULL))
}

## Calculate FWER adjustments as used in the brainDERs project
.calculateFWER <- function(cutoffFstatUsed, areaNull, permutation, nPermute, areaReg) {
    nullList <- split(as.numeric(areaNull), factor(as.numeric(permutation), levels = seq_len(nPermute)))
    nullList[sapply(nullList, length) == 0] <- cutoffFstatUsed
    maxArea <- sapply(nullList, max)
    fwer <- (sapply(areaReg, function(x) sum(maxArea > x)) + 1) / (length(maxArea) + 1)
    return(fwer)
}
