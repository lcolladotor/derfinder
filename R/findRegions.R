#' Find non-zero regions in a Rle
#'
#' Find genomic regions for which a numeric vector is above (or below)
#' predefined thresholds. In other words, this function finds the candidate
#' Differentially Expressed Regions (candidate DERs). This is similar to
#' [regionFinder][bumphunter::regionFinder] and is a helper function for
#' [calculatePvalues].
#'
#' @param position A logical Rle of genomic positions. This is generated in
#' [loadCoverage]. Note that it gets updated in [preprocessCoverage]
#' if `colsubset` is not `NULL`.
#' @param fstats A numeric Rle with the F-statistics. Usually obtained using
#' [calculateStats].
#' @param chr A single element character vector specifying the chromosome name.
#' @param oneTable If `TRUE` only one GRanges is returned.
#' Otherwise, a GRangesList with two components is returned: one for the
#' regions with positive values and one for the negative values.
#' @param maxClusterGap This determines the maximum gap between candidate DERs.
#' It should be greater than `maxRegionGap` (0 by default).
#' @param cutoff Threshold applied to the `fstats` used to determine the
#' regions.
#' @param segmentIR An IRanges object with the genomic positions that are
#' potentials DERs. This is used in [calculatePvalues] to speed up
#' permutation calculations.
#' @param smooth Whether to smooth the F-statistics (`fstats`) or not. This
#' is by default `FALSE`. For RNA-seq data we recommend using `FALSE`.
#' @param weights Weights used by the smoother as described in
#' [smoother][bumphunter::smoother].
#' @param smoothFunction A function to be used for smoothing the F-statistics.
#' Two functions are provided by the `bumphunter` package:
#' [loessByCluster][bumphunter::loessByCluster] and [runmedByCluster][bumphunter::runmedByCluster]. If
#' you are using your own custom function, it has to return a named list with
#' an element called `$fitted` that contains the smoothed F-statistics and
#' an element claled `$smoothed` that is a logical vector indicating
#' whether the F-statistics were smoothed or not. If they are not smoothed, the
#' original values will be used.
#' @param ... Arguments passed to other methods and/or advanced arguments.
#' Advanced arguments:
#' \describe{
#' \item{verbose }{ If `TRUE` basic status updates will be printed along
#' the way.}
#' \item{basic }{ If `TRUE` a DataFrame is returned that has only basic
#' information on the candidate DERs. This is used in [calculatePvalues]
#' to speed up permutation calculations. Default: `FALSE`.}
#' \item{maxRegionGap }{ This determines the maximum number of gaps between two
#' genomic positions to be considered part of the same candidate region. The
#' default is 0L.}
#' }
#' Passed to [extendedMapSeqlevels] and the internal function
#' `.getSegmentsRle` that has by default `verbose = FALSE`.
#'
#' When `smooth = TRUE`, `...` is passed to the internal function
#' `.smootherFstats`. This internal function has the advanced argument
#' `maxClusterGap` (same as above) and passes `...` to
#' [define_cluster] and the formal arguments of `smoothFun`.
#'
#' @return Either a GRanges or a GRangesList as determined by `oneTable`.
#' Each of them has the following metadata variables.
#' \describe{
#' \item{value }{ The mean of the values of `y` for the given region.}
#' \item{area }{  The absolute value of the sum of the values of `y` for
#' the given region.}
#' \item{indexStart }{ The start position of the region in terms of the index
#' for `y`.}
#' \item{indexEnd }{ The end position of the region in terms of the index for
#' `y`.}
#' \item{cluster }{ The cluser ID.}
#' \item{clusterL }{ The total length of the cluster.}
#' }
#'
#' @details [regionFinder][bumphunter::regionFinder] adapted to Rle world.
#'
#' @seealso [calculatePvalues]
#'
#' @references Rafael A. Irizarry, Martin Aryee, Hector Corrada Bravo, Kasper
#' D. Hansen and Harris A. Jaffee. bumphunter: Bump Hunter. R package version
#' 1.1.10.
#'
#' @author Leonardo Collado-Torres
#'
#' @export
#' @importFrom IRanges IRanges start end width Views ranges
#' @import S4Vectors
#' @importFrom GenomicRanges GRanges GRangesList
#' @importMethodsFrom IRanges quantile which length mean rbind
#' @importFrom BiocParallel bpworkers
#' @importFrom bumphunter locfitByCluster runmedByCluster
#' @examples
#' ## Preprocess the data
#' prep <- preprocessCoverage(genomeData,
#'     cutoff = 0, scalefac = 32, chunksize = 1e3,
#'     colsubset = NULL
#' )
#'
#' ## Get the F statistics
#' fstats <- genomeFstats
#'
#' ## Find the regions
#' regs <- findRegions(prep$position, fstats, "chr21", verbose = TRUE)
#' regs
#' \dontrun{
#' ## Once you have the regions you can proceed to annotate them
#' library("bumphunter")
#' genes <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene)
#' annotation <- matchGenes(regs, genes)
#' annotation
#' }
#'
#' # Find regions with smoothing the F-statistics by bumphunter::runmedByCluster
#' regs_smooth <- findRegions(prep$position, fstats, "chr21",
#'     verbose = TRUE,
#'     smoothFunction = bumphunter::runmedByCluster
#' )
#' ## Compare against the object regs obtained earlier
#' regs_smooth
findRegions <- function(position = NULL, fstats, chr, oneTable = TRUE,
    maxClusterGap = 300L, cutoff = quantile(fstats, 0.99, na.rm = TRUE),
    segmentIR = NULL, smooth = FALSE, weights = NULL,
    smoothFunction = bumphunter::locfitByCluster, ...) {

    ## Advanged arguments
    # @param basic If \code{TRUE} a DataFrame is returned that has only basic
    # information on the candidate DERs. This is used in \link{calculatePvalues} to speed up permutation calculations.
    basic <- .advanced_argument("basic", FALSE, ...)


    # @param maxRegionGap This determines the maximum number of gaps between two
    # genomic positions to be considered part of the same candidate Differentially Expressed Region (candidate DER).
    maxRegionGap <- .advanced_argument("maxRegionGap", 0L, ...)

    # @param verbose If \code{TRUE} basic status updates will be printed along the
    # way.
    verbose <- .advanced_argument("verbose", TRUE, ...)

    if (maxClusterGap < maxRegionGap) {
        warning("'maxClusterGap' is less than 'maxRegionGap' which nullifies it's intended use.")
    }

    if (!basic) {
        if (is.null(segmentIR) | smooth) {
            stopifnot(!is.null(position))
        }
        if (smooth) {
            if (verbose) {
                  message(paste(Sys.time(), "findRegions: smoothing"))
              }
            fstats <- .smootherFstats(fstats = fstats, position = position, weights = weights, smoothFunction = smoothFunction, ...)
        }
    } else {
        if (smooth) warning("Ignoring 'smooth' = TRUE since 'basic' = TRUE")
    }


    ## Identify the segments
    if (is.null(segmentIR)) {
        if (verbose) {
              message(paste(
                  Sys.time(),
                  "findRegions: identifying potential segments"
              ))
          }
        if (!any(position)) {
            warning("Found no regions")
            return(NULL)
        }
        segmentIR <- .clusterMakerRle(position,
            maxGap = maxRegionGap,
            ranges = TRUE
        )
    }

    ## Create the F-stats segments
    if (verbose) {
          message(paste(
              Sys.time(),
              "findRegions: segmenting information"
          ))
      }
    segments <- .getSegmentsRle(x = fstats, cutoff = cutoff, ...)

    ## Work only with those that have some information
    hasInfo <- sapply(segments, length) != 0

    ## Stop if there are no segments
    if (!any(hasInfo)) {
        if (verbose) {
              message(paste(
                  Sys.time(),
                  "findRegions: found no segments to work with!!"
              ))
          }
        return(NULL)
    }

    ## Proceed of there is some data to work with
    segments <- segments[hasInfo]

    ## Find the actual DERs
    if (verbose) {
          message(paste(
              Sys.time(),
              "findRegions: identifying candidate regions"
          ))
      }
    ders <- lapply(segments, function(fcut) {
        ## Merge with segment ranges
        all <- c(fcut, segmentIR)

        ## Find all the small pieces
        pieces <- disjoin(all)

        ## Find the actual DERs
        Views(fstats, pieces[queryHits(findOverlaps(pieces, fcut))])
    })

    ## Sadly, this is required to map the positions of the index
    ## to the chr positions.  It's 275 mb in RAM for a length of
    ## 72097604 instead of 4.7 Mb in Rle world.  The good thing is
    ## that it's temporary and the user will not need to save this
    if (!basic) {
        pos <- which(position)
    }

    ## Build the output shell
    res <- vector("list", sum(hasInfo))
    names(res) <- names(hasInfo)[hasInfo]

    ## Use UCSC names for homo_sapiens by default
    chr <- extendedMapSeqlevels(chr, ...)

    for (i in names(hasInfo)[hasInfo]) {
        if (!basic) {
            ## Define the chr ranges
            pos.ir <- IRanges(
                start = pos[start(ders[[i]])],
                end = pos[end(ders[[i]])]
            )

            ## Actually build the GRanges
            res[[i]] <- GRanges(
                seqnames = Rle(chr, length(ders[[i]])),
                ranges = pos.ir, value = mean(ders[[i]]),
                area = abs(sum(ders[[i]])), indexStart = start(ders[[i]]),
                indexEnd = end(ders[[i]])
            )

            ## Identify clusters
            if (verbose) {
                  message(paste(
                      Sys.time(),
                      "findRegions: identifying region clusters"
                  ))
              }
            regionPos <- coverage(res[[i]])[[chr]]
            runValue(regionPos) <- as.logical(runValue(regionPos))
            cluster <- .clusterMakerRle(regionPos, maxClusterGap)

            ## Extract DERs ranges and shift the IR to the cluster' scale
            derCWs <- cumsum(width(ranges(ders[[i]])))
            derIR <- IRanges(start = c(1, derCWs[-length(derCWs)] +
                1), end = derCWs)
            clus <- Views(cluster, derIR)

            ## Finally, identify the clusters
            clusterFinal <- as.integer(mean(clus))
            clusterWidth <- tapply(pos.ir, clusterFinal, function(x) {
                max(end(x)) - min(start(x)) + 1
            })

            res[[i]]$cluster <- Rle(clusterFinal)
            res[[i]]$clusterL <- Rle(clusterWidth[clusterFinal])
        } else {
            ## Actually build the GRanges
            res[[i]] <- DataFrame(
                area = Rle(abs(sum(ders[[i]]))),
                width = Rle(width(ders[[i]])), stat = Rle(mean(ders[[i]])),
                check.names = FALSE
            )
        }
    }

    if (!basic) {
        ## Fix names and format
        names(res) <- gsub("Index", "", names(res))
        res <- GRangesList(res)

        ## Finish up
        if (oneTable) {
            res <- unlist(res)
        }
    } else {
        res <- do.call(rbind, res)
    }

    return(res)
}







#' Segment a Rle into positive, zero, and negative regions
#'
#' Given two cutoffs, L and U, this function slices a numerical Rle into up and
#' down sections. It is a wrapper for [slice][IRanges::slice] with functionality
#' inspired from [getSegments][bumphunter::getSegments].
#'
#'
#' @param x A numeric Rle.
#' @param cutoff A numeric vector of length either 1 or 2. If length is 1, U
#' will be cutoff and L will be -cutoff. Otherwise it specifies L and U. The
#' function will furthermore always use the minimum of cutoff for L and the
#' maximum for U.
#' @param ... Arguments passed to other methods and/or advanced arguments.
#'
#' @return A list of IRanges objects, one for the up segments and one for the
#' down segments.
#'
#' @seealso [getSegments][bumphunter::getSegments], [slice][IRanges::slice],
#' [findRegions]
#'
#' @author Leonardo Collado-Torres
#'
#' @keywords internal
#' @importMethodsFrom IRanges quantile
#' @importFrom IRanges slice
#' @import S4Vectors
#' @examples
#' library("IRanges")
#' set.seed(20130725)
#' pos <- Rle(sample(c(TRUE, FALSE), 1e5, TRUE, prob = c(0.05, 0.95)))
#' data <- Rle(rnorm(sum(pos)))
#' cutoff <- quantile(data, .99, na.rm = TRUE)
#'
#' ## It's quite fast
#' system.time(segs <- derfinder:::.getSegmentsRle(data, cutoff, verbose = TRUE))
#' \dontrun{
#' ## The output is different in look than the one from getSegments() but it's
#' ## use is similar.
#' ## Plus it can be transformed into the same format as the ouptut from
#' ## .getSegmentsRle().
#' library("bumphunter")
#' cluster <- derfinder:::.clusterMakerRle(pos, 100L)
#' foo <- function() {
#'     segs2 <- getSegments(as.numeric(data), as.integer(cluster), cutoff,
#'         assumeSorted = TRUE
#'     )[c("upIndex", "dnIndex")]
#'     segs.ir <- lapply(segs2, function(ind) {
#'         tmp <- lapply(ind, function(segment) {
#'             c("start" = min(segment), "end" = max(segment))
#'         })
#'         info <- do.call(rbind, tmp)
#'         IRanges(start = info[, "start"], end = info[, "end"])
#'     })
#'     return(segs.ir)
#' }
#' identical(foo(), segs)
#' }
#'
#' @noRd
.getSegmentsRle <- function(x, cutoff = quantile(x, 0.99, na.rm = TRUE), ...) {

    ## Advanged arguments
    # @param verbose If \code{TRUE} basic status updates will be printed along the
    # way.
    verbose <- .advanced_argument("verbose", FALSE, ...)

    ## Select the cutoff
    if (verbose) {
        message(paste(
            Sys.time(),
            ".getSegmentsRle: segmenting with cutoff(s)",
            paste(cutoff, collapse = ", ")
        ))
    }
    stopifnot(length(cutoff) <= 2)
    if (length(cutoff) == 1) {
        cutoff <- c(-cutoff, cutoff)
    }
    cutoff <- sort(cutoff)

    ## Find the segments
    result <- lapply(c("upIndex", "dnIndex"), function(ind) {
        if (ind == "upIndex") {
            fcut <- slice(x = x, lower = cutoff[2], rangesOnly = TRUE)
        } else {
            fcut <- slice(x = x, upper = cutoff[1], rangesOnly = TRUE)
        }
        return(fcut)
    })
    names(result) <- c("upIndex", "dnIndex")

    ## Done!
    return(result)
}




#' Make clusters of genomic locations based on distance in Rle() world
#'
#' Genomic locations are grouped into clusters based on distance: locations
#' that are close to each other are assigned to the same cluster. The operation
#' is performed on each chromosome independently. This is very similar to
#' [clusterMaker][bumphunter::clusterMaker].
#'
#' @details
#' [clusterMaker][bumphunter::clusterMaker] adapted to Rle world. Assumes that the data
#' is sorted and that everything is in a single chromosome.
#' It is also almost as fast as the original version with the advantage that
#' everything is in Rle() world.
#'
#' It is a a helper function for [findRegions].
#'
#' @param position A logical Rle indicating the chromosome positions.
#' @param maxGap An integer. Genomic locations within `maxGap` from each
#' other are labeled as part of the same cluster.
#' @param ranges If `TRUE` then an IRanges object is returned instead of
#' the usual integer Rle.
#' @param ... Arguments passed to other methods and/or advanced arguments.
#'
#' @return An integer Rle with the cluster IDs. If `ranges=TRUE` then it
#' is an IRanges object with one range per cluster.
#'
#' @keywords internal
#' @seealso [clusterMaker][bumphunter::clusterMaker], [findRegions]
#' @references Rafael A. Irizarry, Martin Aryee, Hector Corrada Bravo, Kasper
#' D. Hansen and Harris A. Jaffee. bumphunter: Bump Hunter. R package version
#' 1.1.10.
#' @author Leonardo Collado-Torres
#'
#' @importFrom IRanges IRanges start end reduce Views runLength
#' @importMethodsFrom IRanges length sum
#' @import S4Vectors
#'
#' @examples
#' library("IRanges")
#' set.seed(20130725)
#' pos <- Rle(sample(c(TRUE, FALSE), 1e5, TRUE, prob = c(0.05, 0.95)))
#' cluster <- .clusterMakerRle(pos, 100L)
#' cluster
#' @noRd

.clusterMakerRle <- function(position, maxGap = 300L, ranges = FALSE, ...) {
    ## Instead of using which(), identify the regions of the chr
    ## with data
    ir <- IRanges(
        start = start(position)[runValue(position)],
        end = end(position)[runValue(position)]
    )

    ## Apply the gap reduction
    ir.red <- reduce(ir, min.gapwidth = maxGap + 1)

    ## Identify the clusters
    clusterIDs <- Rle(seq_len(length(ir.red)), sum(Views(
        position,
        ir.red
    )))
    ## Note that sum(Views(pos, ir.red)) is faster than
    ## sapply(ir.red, function(x) sum(pos[x]))

    ## Group the information into an IRanges object
    if (ranges) {
        csum <- cumsum(runLength(clusterIDs))
        result <- IRanges(start = c(1, csum[-length(csum)] +
            1), end = csum)
    } else {
        result <- clusterIDs
    }

    ## Done
    return(result)
}


.smootherFstats <- function(fstats, position, weights = NULL,
    smoothFunction = bumphunter::locfitByCluster, ...) {
    ## Based on bumphunter::smoother

    ## Advanced arguments
    # @param maxClusterGap This determines the maximum gap between candidate DERs.
    # It should be greater than \code{maxRegionGap} (0 by default).
    maxClusterGap <- .advanced_argument("maxClusterGap", 300L, ...)

    ## Identify clusters
    cluster <- .clusterMakerRle(position, maxGap = maxClusterGap)

    ## Define computing cluster
    BPPARAM <- define_cluster(...)
    cores <- bpworkers(BPPARAM)

    if (cores > 1) {
        IndexesChunks <- split(runValue(cluster), cut(runValue(cluster),
            breaks = cores
        ))
        iChunks <- rep(seq_len(cores), sapply(IndexesChunks, function(i) {
              sum(cluster %in% i)
          }))
    } else {
        iChunks <- rep(1, length(cluster))
    }

    ## Define the chunks of data to process in parallel
    fstatsChunks <- split(fstats, iChunks)
    posChunks <- split(which(position), iChunks)
    clusterChunks <- split(cluster, iChunks)

    if (is.null(weights)) {
        weightChunks <- vector("list", length = length(unique(iChunks)))
    } else {
        weightChunks <- split(weights, iChunks)
    }

    ## Run in parallel
    res <- bpmapply(.smoothFstatsFun, fstatsChunks, posChunks, clusterChunks,
        weightChunks,
        MoreArgs = list(smoothFun = smoothFunction, ...),
        BPPARAM = BPPARAM
    )

    ## Get back a Rle
    res <- unlist(RleList(res), use.names = FALSE)

    return(res)
}


## Helper function to smoothFstats() for running in parallel and coercing
## the result back to a Rle object

.smoothFstatsFun <- function(y, x, cluster, weights, smoothFun, ...) {
    hostPackage <- environmentName(environment(smoothFun))
    requireNamespace(hostPackage)
    smoothed <- .runFunFormal(smoothFun, y = y, x = x, cluster = cluster, weights = weights, ...)

    ## Use original values if they were not smoothed
    if (any(!smoothed$smoothed)) {
        if (is(y, "Rle")) {
            smoothed$fitted[!smoothed$smoothed] <- as.vector(y[!smoothed$smoothed])
        } else {
            smoothed$fitted[!smoothed$smoothed] <- y[!smoothed$smoothed]
        }
    }

    ## Extract only the smoothed data
    res <- Rle(smoothed$fitted)
    return(res)
}
