#' Find non-zero regions in a Rle
#'
#' Find genomic regions for which a numeric vector is above (or below) 
#' predefined thresholds. In other words, this function finds the candidate 
#' Differentially Expressed Regions (candidate DERs). This is similar to 
#' \link[bumphunter]{regionFinder} and is a helper function for 
#' \link{calculatePvalues}.
#' 
#' @param position A logical Rle of genomic positions. This is generated in 
#' \link{loadCoverage}. Note that it gets updated in \link{preprocessCoverage} 
#' if \code{colsubset} is not \code{NULL}.
#' @param fstats A numeric Rle with the F-statistics. Usually obtained using 
#' \link{calculateStats}.
#' @param chr A single element character vector specifying the chromosome name.
#' @param oneTable If \code{TRUE} only one GRanges is returned. 
#' Otherwise, a GRangesList with two components is returned: one for the 
#' regions with positive values and one for the negative values.
#' @param maxClusterGap This determines the maximum gap between candidate DERs. 
#' It should be greater than \code{maxRegionGap} (0 by default).
#' @param cutoff Threshold applied to the \code{fstats} used to determine the #' regions.
#' @param segmentIR An IRanges object with the genomic positions that are
#' potentials DERs. This is used in \link{calculatePvalues} to speed up
#' permutation calculations.
#' @param ... Arguments passed to other methods and/or advanced arguments.
#'
#' @return Either a GRanges or a GRangesList as determined by \code{oneTable}. 
#' Each of them has the following metadata variables.
#' \describe{
#' \item{value }{ The mean of the values of \code{y} for the given region.}
#' \item{area }{  The absolute value of the sum of the values of \code{y} for 
#' the given region.}
#' \item{indexStart }{ The start position of the region in terms of the index 
#' for \code{y}.}
#' \item{indexEnd }{ The end position of the region in terms of the index for 
#' \code{y}.}
#' \item{cluster }{ The cluser ID.}
#' \item{clusterL }{ The total length of the cluster.}
#' }
#'
#' @details \link[bumphunter]{regionFinder} adapted to Rle world.
#'
#' @seealso \link{calculatePvalues}
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
#' @examples
#' ## Preprocess the data
#' prep <- preprocessCoverage(genomeData, cutoff=0, scalefac=32, chunksize=1e3, 
#'     colsubset=NULL)
#' 
#' ## Get the F statistics
#' fstats <- genomeFstats
#'
#' ## Find the regions
#' regs <- findRegions(prep$position, fstats, 'chr21', verbose=TRUE)
#' regs
#'
#' \dontrun{
#' ## Once you have the regions you can proceed to annotate them
#' library('bumphunter')
#' library('TxDb.Hsapiens.UCSC.hg19.knownGene')
#' genes <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' annotation <- matchGenes(regs, genes)
#' annotation
#' }

findRegions <- function(position = NULL, fstats, chr, oneTable = TRUE, 
    maxClusterGap = 300L, cutoff = quantile(fstats, 0.99), segmentIR = NULL,
    ...){
    
    ## Advanged arguments
#' @param basic If \code{TRUE} a DataFrame is returned that has only basic
#' information on the candidate DERs. This is used in \link{calculatePvalues} to speed up permutation calculations.
    basic <- .advanced_argument('basic', FALSE, ...)


#' @param maxRegionGap This determines the maximum number of gaps between two 
#' genomic positions to be considered part of the same candidate Differentially Expressed Region (candidate DER). 
    maxRegionGap <- .advanced_argument('maxRegionGap', 0L, ...)

#' @param verbose If \code{TRUE} basic status updates will be printed along the 
#' way.
    verbose <- .advanced_argument('verbose', TRUE, ...)

    
    if (!basic) {
        if (is.null(segmentIR)) {
            stopifnot(!is.null(position))
        }
    }
    if (maxClusterGap < maxRegionGap) {
        warning("'maxClusterGap' is less than 'maxRegionGap' which nullifies it's intended use.")
    }
    
    ## Identify the segments
    if (is.null(segmentIR)) {
        if (verbose) 
            message(paste(Sys.time(),
                'findRegions: identifying potential segments'))
        segmentIR <- .clusterMakerRle(position, maxGap = maxRegionGap,
            ranges = TRUE)
    }
    
    ## Create the F-stats segments
    if (verbose) 
        message(paste(Sys.time(),
            'findRegions: segmenting F-stats information'))
    segments <- .getSegmentsRle(x = fstats, cutoff = cutoff, ...)
    
    ## Work only with those that have some information
    hasInfo <- sapply(segments, length) != 0
    
    ## Stop if there are no segments
    if (!any(hasInfo)) {
        if (verbose) 
            message(paste(Sys.time(),
                'findRegions: found no segments to work with!!'))
        return(NULL)
    }
    
    ## Proceed of there is some data to work with
    segments <- segments[hasInfo]
    
    ## Find the actual DERs
    if (verbose) 
        message(paste(Sys.time(),
            'findRegions: identifying candidate regions'))
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
    res <- vector('list', sum(hasInfo))
    names(res) <- names(hasInfo)[hasInfo]
    
    ## Use UCSC names for homo_sapiens by default
    chr <- extendedMapSeqlevels(chr, ...)
    
    for (i in names(hasInfo)[hasInfo]) {
        if (!basic) {
            ## Define the chr ranges
            pos.ir <- IRanges(start = pos[start(ders[[i]])], 
                width = width(ders[[i]]))
            
            ## Actually build the GRanges
            res[[i]] <- GRanges(seqnames = Rle(chr, length(ders[[i]])), 
                ranges = pos.ir, value = mean(ders[[i]]), 
                area = abs(sum(ders[[i]])), indexStart = start(ders[[i]]),
                indexEnd = end(ders[[i]]))
            
            ## Identify clusters
            if (verbose) 
                message(paste(Sys.time(),
                    'findRegions: identifying region clusters'))
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
            res[[i]] <- DataFrame(area = Rle(abs(sum(ders[[i]]))), 
                width = Rle(width(ders[[i]])), stat = Rle(mean(ders[[i]])))
        }
        
    }
    
    if (!basic) {
        ## Fix names and format
        names(res) <- gsub('Index', '', names(res))
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
#' down sections. It is a wrapper for \link[IRanges]{slice} with functionality 
#' inspired from \link[bumphunter]{getSegments}.
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
#' @seealso \link[bumphunter]{getSegments}, \link[IRanges]{slice}, 
#' \link{findRegions}
#'
#' @author Leonardo Collado-Torres
#'
#' @keywords internal 
#' @importMethodsFrom IRanges quantile
#' @importFrom IRanges slice
#' @import S4Vectors
#' @examples
#' library('IRanges')
#' set.seed(20130725)
#' pos <- Rle(sample(c(TRUE, FALSE), 1e5, TRUE, prob=c(0.05, 0.95)))
#' data <- Rle(rnorm(sum(pos)))
#' cutoff <- quantile(data, .99)
#'
#' ## It's quite fast
#' system.time(segs <- derfinder:::.getSegmentsRle(data, cutoff, verbose=TRUE))
#' 
#' \dontrun{
#' ## The output is different in look than the one from getSegments() but it's 
#' ## use is similar.
#' ## Plus it can be transformed into the same format as the ouptut from 
#' ## .getSegmentsRle().
#' library('bumphunter')
#' cluster <- derfinder:::.clusterMakerRle(pos, 100L)
#' foo <- function() {
#'     segs2 <- getSegments(as.numeric(data), as.integer(cluster), cutoff, 
#'     assumeSorted=TRUE)[c('upIndex', 'dnIndex')]
#'     segs.ir <- lapply(segs2, function(ind) {
#'         tmp <- lapply(ind, function(segment) {
#'             c('start'=min(segment), 'end'=max(segment))
#'         })
#'         info <- do.call(rbind, tmp)
#'         IRanges(start=info[,'start'], end=info[,'end'])
#'     })
#'     return(segs.ir)
#' }
#' identical(foo(), segs) 
#'
#' }
#'

.getSegmentsRle <- function(x, cutoff = quantile(x, 0.99), ...) {
    
    ## Advanged arguments
#' @param verbose If \code{TRUE} basic status updates will be printed along the 
#' way.
    verbose <- .advanced_argument('verbose', FALSE, ...)
    
    ## Select the cutoff
    if (verbose) message(paste(Sys.time(),
        '.getSegmentsRle: segmenting with cutoff(s)',
        paste(cutoff, collapse=', ')))
    stopifnot(length(cutoff) <= 2)
    if (length(cutoff) == 1) {
        cutoff <- c(-cutoff, cutoff)
    }
    cutoff <- sort(cutoff)
    
    ## Find the segments
    result <- lapply(c('upIndex', 'dnIndex'), function(ind) {
        if(ind == 'upIndex') {
            fcut <- slice(x=x, lower=cutoff[2], rangesOnly=TRUE)
        } else {
            fcut <- slice(x=x, upper=cutoff[1], rangesOnly=TRUE)
        }
        return(fcut)
    })
    names(result) <- c('upIndex', 'dnIndex')

    ## Done!
    return(result)
}




#' Make clusters of genomic locations based on distance in Rle() world
#'
#' Genomic locations are grouped into clusters based on distance: locations 
#' that are close to each other are assigned to the same cluster. The operation 
#' is performed on each chromosome independently. This is very similar to 
#' \link[bumphunter]{clusterMaker}.
#'
#' @details
#' \link[bumphunter]{clusterMaker} adapted to Rle world. Assumes that the data 
#' is sorted and that everything is in a single chromosome.
#' It is also almost as fast as the original version with the advantage that 
#' everything is in Rle() world.
#' 
#' It is a a helper function for \link{findRegions}.
#' 
#' @param position A logical Rle indicating the chromosome positions.
#' @param maxGap An integer. Genomic locations within \code{maxGap} from each 
#' other are labeled as part of the same cluster.
#' @param ranges If \code{TRUE} then an IRanges object is returned instead of 
#' the usual integer Rle.
#' @param ... Arguments passed to other methods and/or advanced arguments.
#'
#' @return An integer Rle with the cluster IDs. If \code{ranges=TRUE} then it 
#' is an IRanges object with one range per cluster.
#'
#' @keywords internal 
#' @seealso \link[bumphunter]{clusterMaker}, \link{findRegions}
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
#' library('IRanges')
#' set.seed(20130725)
#' pos <- Rle(sample(c(TRUE, FALSE), 1e5, TRUE, prob=c(0.05, 0.95)))
#' cluster <- .clusterMakerRle(pos, 100L)
#' cluster
#'

.clusterMakerRle <- function(position, maxGap = 300L, ranges = FALSE, ...) {
    ## Instead of using which(), identify the regions of the chr
    ## with data
    ir <- IRanges(start = start(position)[runValue(position)], 
        end = end(position)[runValue(position)])
    
    ## Apply the gap reduction
    ir.red <- reduce(ir, min.gapwidth = maxGap + 1)
    
    ## Identify the clusters
    clusterIDs <- Rle(seq_len(length(ir.red)), sum(Views(position, 
        ir.red)))
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

