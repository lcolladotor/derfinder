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
#' @param fstats A numeric Rle with the F-statistics. Normally obtained using 
#' \link{calculateStats}.
#' @param chr A single element character vector specifying the chromosome name.
#' @param oneTable If \code{TRUE} only one results GRanges is returned. 
#' Otherwise, a GRangesList with two components is returned: one for the 
#' regions with positive values and one for the negative values.
#' @param maxRegionGap This determines the maximum number of gaps between two 
#' genomic positions to be considered part of the same candidate Differentially 
#' Expressed Region (candidate DER).
#' @param maxClusterGap This determines the maximum gap between candidate DERs. 
#' It should be greater than \code{maxRegionGap}.
#' @param cutoff This argument is passed to \link{getSegmentsRle}.
#' @param segmentIR An IRanges object with the genomic positions that are 
#' potentials DERs, normally given by \link{clusterMakerRle} with 
#' \code{maxGap=maxRegionGap} and \code{ranges=TRUE}. This is used in 
#' \link{calculatePvalues} to speed up permutation calculations.
#' @param basic If \code{TRUE} a DataFrame is returned that has only basic 
#' information on the candidate DERs. This is used in \link{calculatePvalues} 
#' to speed up permutation calculations.
#' @param chrsStyle The naming style of the chromosomes. By default, UCSC. See 
#' \link[GenomeInfoDb]{seqlevelsStyle}.
#' @param verbose If \code{TRUE} basic status updates will be printed along the 
#' way.
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
#' @seealso \link[bumphunter]{regionFinder}, \link{calculatePvalues}
#' @references Rafael A. Irizarry, Martin Aryee, Hector Corrada Bravo, Kasper 
#' D. Hansen and Harris A. Jaffee. bumphunter: Bump Hunter. R package version 
#' 1.1.10.
#'
#' @author Leonardo Collado-Torres
#' @export
#' @aliases find_regions
#' @importFrom IRanges IRanges start end width Views ranges
#' @importFrom S4Vectors Rle runLength  runValue 'runValue<-' DataFrame
#' @importFrom GenomicRanges GRanges GRangesList
#' @importFrom GenomeInfoDb mapSeqlevels
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
#' ## Compare vs bumphunter
#' library('bumphunter')
#' regs2 <- regionFinder(as.numeric(fstats), rep('chr21', length(fstats)), 
#'     which(prep$position), cluster=NULL, assumeSorted=TRUE, verbose=TRUE, 
#'     order=FALSE, maxGap=1)
#' regs2
#' ## Note that regs$L can be calculated with width(regs)
#' identical(width(regs), as.integer(regs2$L))
#' ## Time comparison
#' library('microbenchmark')
#' micro <- microbenchmark(findRegions(prep$position, fstats, 'chr21', 
#'     verbose=FALSE, basic=TRUE), regionFinder(as.numeric(fstats), 
#'     rep('chr21', length(fstats)), which(prep$position), cluster=NULL, 
#'     assumeSorted=TRUE, verbose=FALSE, order=FALSE, maxGap=1))
#' levels(micro$expr) <- c('new', 'original')
#' micro
#' ## The bumphunter function regionFinder() is faster in small data sets.
#'
#' ## Once you have the regions you can proceed to annotate them
#' annotation <- annotateNearest(regs, 'hg19')
#' annotation
#' }

findRegions <- function(position = NULL, fstats, chr, oneTable = TRUE, 
    maxRegionGap = 0L, maxClusterGap = 300L, cutoff = quantile(fstats, 
        0.99), segmentIR = NULL, basic = FALSE, chrsStyle = "UCSC", 
        verbose = TRUE) {
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
                "findRegions: identifying potential segments"))
        segmentIR <- clusterMakerRle(position, maxRegionGap, 
            ranges = TRUE)
    }
    
    ## Create the F-stats segments
    if (verbose) 
        message(paste(Sys.time(),
            "findRegions: segmenting F-stats information"))
    segments <- getSegmentsRle(x = fstats, cutoff = cutoff, verbose = verbose)
    
    ## Work only with those that have some information
    hasInfo <- sapply(segments, length) != 0
    
    ## Stop if there are no segments
    if (!any(hasInfo)) {
        if (verbose) 
            message(paste(Sys.time(),
                "findRegions: found no segments to work with!!"))
        return(NULL)
    }
    
    ## Proceed of there is some data to work with
    segments <- segments[hasInfo]
    
    ## Find the actual DERs
    if (verbose) 
        message(paste(Sys.time(),
            "findRegions: identifying candidate DERs"))
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
    
    ## Use UCSC names by default
    chr <- mapSeqlevels(chr, chrsStyle)
    
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
                    "findRegions: identifying DER clusters"))
            regionPos <- coverage(res[[i]])[[chr]]
            runValue(regionPos) <- as.logical(runValue(regionPos))
            cluster <- clusterMakerRle(regionPos, maxClusterGap)
            
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

#' @export
find_regions <- findRegions
