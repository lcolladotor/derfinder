#' Extract coverage information for a set of regions
#'
#' This function extracts the raw coverage information calculated by 
#' \link{fullCoverage} at each base for a set of regions found with 
#' \link{calculatePvalues}. It can further calculate the mean coverage per 
#' sample for each region.
#' 
#' @param fullCov A list where each element is the result from 
#' \link{loadCoverage} used with \code{cutoff=NULL}. Can be generated using 
#' \link{fullCoverage}.
#' @param regions The \code{$regions} output from \link{calculatePvalues}. It 
#' is important that the seqlengths information is provided.
#' @param totalMapped The total number of reads mapped for each sample. 
#' Providing this data returns coverage adjusted totalMapped per million at 
#' each base.
#' @param mc.cores The number of cores to use for computing coverage. Default=1
#' @param chrsStyle The naming style of the chromosomes. By default, UCSC. See 
#' \link[GenomeInfoDb]{seqlevelsStyle}.
#' @param verbose If \code{TRUE} basic status updates will be printed along the 
#' way.
#'
#' @return a list of data.frame where each data.frame has the coverage 
#' information (nrow = width of region, ncol = number of samples) for a given 
#' region. The names of the list correspond to the region indexes in 
#' \code{regions}
#' 
#' @author Andrew Jaffe, Leonardo Collado-Torres
#' @seealso \link{fullCoverage}, \link{calculatePvalues}
#' @export
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomeInfoDb seqlevels seqlevelsStyle 'seqlevelsStyle<-'
#' mapSeqlevels
#' @importMethodsFrom GenomicRanges names 'names<-' length '[' coverage sort 
#' width c
#' @importMethodsFrom IRanges subset as.data.frame as.factor
#'
#' @examples
#' ## Obtain fullCov object
#' fullCov <- list('21'=genomeDataRaw$coverage)
#' 
#' ## Assign chr lengths using hg19 information, use only first two regions
#' library('GenomicRanges')
#' data(hg19Ideogram, package = 'biovizBase', envir = environment())
#' regions <- genomeRegions$regions[1:2]
#' seqlengths(regions) <- seqlengths(hg19Ideogram)[names(seqlengths(regions))]
#'
#' ## Finally, get the region coverage
#' regionCov <- getRegionCoverage(fullCov=fullCov, regions=regions)

getRegionCoverage <- function(fullCov, regions, totalMapped = NULL, 
    mc.cores = 1, chrsStyle = "UCSC", verbose = TRUE) {
    
    names(regions) <- seq_len(length(regions))  # add names
    
    ## Use UCSC style names by default
    names(fullCov) <- mapSeqlevels(names(fullCov), chrsStyle)
    seqlevelsStyle(regions) <- chrsStyle
    
    ## Warning when seqlengths are not specified
    if (any(is.na(seqlengths(regions)))) 
        warning("'regions' does not have seqlengths assigned! In some cases, this can lead to erroneous results. getRegionCoverage() will proceed, but please check for other warnings or errors.")
    
    # split by chromosome
    grl <- split(regions, as.factor(seqnames(regions)))  
    
    counts <- mclapply(grl, function(g) {
        # now can be parallel
        if (verbose) 
            cat(".")
        thechr <- as.character(unique(seqnames(g)))
        yy <- fullCov[[thechr]][ranges(g), ]  # better subset
        # depth-adjust, like for plotting
        if (!is.null(totalMapped)) {
            yy <- DataFrame(mapply(function(x, d) x/(d/1e+06), 
                yy, totalMapped))
        }
        ind <- rep(names(g), width(g))  # to split along
        ind <- factor(ind, levels = unique(ind))  # make factor in order
        # split(yy,ind) # 'CompressedSplitDataFrameList', faster but
        # less clear how to unlist below, so leave out
        split(as.data.frame(yy), ind)
    }, mc.cores = mc.cores)
    covList <- do.call("c", counts)  # collect list elements into one large list
    
    # put in original order
    names(covList) <- sapply(strsplit(names(covList), "\\."), 
        "[", 2)
    theData <- covList[order(as.numeric(names(covList)))]
    
    if (sum(sapply(theData, nrow)) != sum(width(regions))) {
        stop("The total width of the regions did not match with the dimensions of the extracted coverage data. Check that 'regions' has seqlengths specified correctly.")
    }
    
    return(theData)
} 
