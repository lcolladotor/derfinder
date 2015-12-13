#' Identify regions data by a coverage filter and get a count matrix from
#' BigWig files
#'
#' Rail (available at http://rail.bio) generates coverage BigWig files. These
#' files are faster to load in R and to process. Rail creates an un-adjusted
#' coverage BigWig file per sample and an adjusted summary coverage BigWig file 
#' by chromosome (median or mean). \link{railMatrix} reads in the mean (or 
#' median) coverage BigWig file and applies a threshold cutoff to identify 
#' expressed-regions (ERs). Then it goes back to the sample coverage BigWig 
#' files and extracts the base level coverage for each sample. Finally it
#' summarizes this information in a matrix with one row per ERs and one column
#' per sample. This function is similar to \link{regionMatrix} but is faster
#' thanks to the advantages presented by BigWig files.

#' Given a set of un-filtered coverage data (see \link{fullCoverage}), create
#' candidate regions by applying a cutoff on the coverage values,
#' and obtain a count matrix where the number of rows corresponds to the number
#' of candidate regions and the number of columns corresponds to the number of 
#' samples. The values are the mean coverage for a given sample for a given 
#' region.
#' 

#' @param chrs A character vector with the names of the chromosomes.
#' @param summaryFiles A character vector (or BigWigFileList) with the paths to 
#' the summary BigWig files created by Rail. Either mean or median files. These 
#' are library size adjusted by default in Rail. The order of the files in this 
#' vector must match the order in \code{chrs}.
#' @param sampleFiles A character vector with the paths to the BigWig files
#' by sample. These files are unadjusted for library size.
#' @param targetSize The target library size to adjust the coverage to. Used
#' only when \code{totalMapped} is specified. By default, it adjusts to 
#' libraries with 40 million reads, which matches the default used in Rail.
#' @inheritParams filterData
#' @param L The width of the reads used. Either a vector of length 1 or length
#' equal to the number of samples.
#' @param maxClusterGap This determines the maximum gap between candidate ERs.
#' @param returnBP If \code{TRUE}, returns \code{$bpCoverage} explained below.
#' @param file.cores Number of cores used for loading the BigWig files per chr.
#' @param ... Arguments passed to other methods and/or advanced arguments.
#'
#' @return A list with one entry per chromosome. Then per chromosome, a list 
#' with three components.
#' \describe{
#' \item{regions }{ A set of regions based on the coverage filter cutoff as
#' returned by \link{findRegions}.}
#' \item{bpCoverage }{ A list with one element per region. Each element is a matrix with numbers of rows equal to the number of base pairs in the region and number of columns equal to the number of samples. It contains the base-level coverage information for the regions. Only returned when \code{returnBP = TRUE}. Set to \code{FALSE} when the number of samples is large.}
#' \item{coverageMatrix }{  A matrix with the mean coverage by sample for each
#' candidate region.}
#' }
#'
#' @details This function uses several other \link{derfinder-package} 
#' functions. Inspect the code if interested.
#'
#' You should use at most one core per chromosome.
#'
#' @author Leonardo Collado-Torres
#' @export
#'
#' @importFrom BiocParallel bpmapply
#' @importFrom GenomeInfoDb 'seqlengths<-'
#'
#' @examples
#'
#' ## BigWig files are not supported in Windows
#' if(.Platform$OS.type != 'windows') {
#'     ## Get data
#'     library('derfinderData')
#'     
#'     ## Identify sample files
#'     sampleFiles <- rawFiles(system.file('extdata', 'AMY', package = 
#'         'derfinderData'), samplepatt = 'bw', fileterm = NULL)
#'     names(sampleFiles) <- gsub('.bw', '', names(sampleFiles))
#'     
#'     ## Create the mean bigwig file. This file is normally created by Rail
#'     ## but in this example we'll create it manually.
#'     library('GenomicRanges')
#'     fullCov <- fullCoverage(files = sampleFiles, chrs = 'chr21')
#'     meanCov <- Reduce('+', fullCov$chr21) / ncol(fullCov$chr21)
#'     createBw(list('chr21' = DataFrame('meanChr21' = meanCov)), keepGR = 
#'         FALSE)
#'     
#'     summaryFile <- 'meanChr21.bw'
#'     
#'     ## Get the regions
#'     regionMat <- railMatrix(chrs = 'chr21', summaryFiles = summaryFile, 
#'         sampleFiles = sampleFiles, L = 76, cutoff = 5, maxClusterGap = 3000L)
#'     
#'     ## Explore results
#'     names(regionMat$chr21)
#'     regionMat$chr21$regions
#'     dim(regionMat$chr21$coverageMatrix)
#' }


railMatrix <- function(chrs, summaryFiles, sampleFiles, L = NULL, cutoff = NULL, maxClusterGap = 300L, totalMapped = NULL, targetSize = 40e6, returnBP = TRUE, file.cores = 1L, ...) {
    stopifnot(length(chrs) == length(summaryFiles))
    ## In the future, summaryFiles might be a vector or a list, but the length
    ## check must be maintained. It might be a list if the mean bigwig info
    ## is split into different files.
    
    ## Have to filter by something
    stopifnot(!is.null(cutoff))
    
    ## Verbose value
    verbose <- .advanced_argument('verbose', TRUE, ...)
    
    ## Define cluster used for per chromosome
    BPPARAM <- .define_cluster(...)
    
    ## Define cluster used for loading BigWig files    
    if(file.cores == 1L) {
        BPPARAM.railChr <- .advanced_argument('BPPARAM.railChr', SerialParam(),
            ...)
    } else {
        mc.outfile <- .advanced_argument('mc.outfile', 
            Sys.getenv('SGE_STDERR_PATH'), ...)
        BPPARAM.railChr <- .advanced_argument('BPPARAM.railChr', 
            SnowParam(workers = file.cores, outfile = mc.outfile), ...)
    }
    
    if(!is.null(totalMapped) & targetSize != 0) {
        mappedPerXM <- totalMapped / targetSize
    } else {
        mappedPerXM <- NULL
    }
    
    ## Chunksize to use
    chunksize <- .advanced_argument('chunksize', 1000, ...)
    
    
    regionMat <- bpmapply(.railMatrixChr, chrs, summaryFiles, SIMPLIFY = FALSE, MoreArgs = list(sampleFiles = sampleFiles, L = L, maxClusterGap = maxClusterGap, cutoff = cutoff, mappedPerXM = mappedPerXM, returnBP = returnBP, chunksize = chunksize, BPPARAM.railChr = BPPARAM.railChr, verbose = verbose), BPPARAM = BPPARAM)
    return(regionMat)
}


.railMatrixChr <- function(chr, summaryFile, sampleFiles, L = NULL, cutoff = NULL,  maxClusterGap = 300L, mappedPerXM = mappedPerXM, returnBP = TRUE, chunksize = 1000, BPPARAM.railChr = BPPARAM.railChr, verbose = TRUE) {
    meanCov <- loadCoverage(files = summaryFile, chr = chr)    
    regs <- findRegions(position = Rle(TRUE, length(meanCov$coverage[[1]])), fstats = meanCov$coverage[[1]], chr = chr, maxClusterGap = maxClusterGap, cutoff = cutoff)
    
    ## Format appropriately
    names(regs) <- seq_len(length(regs))
    
    ## Set the length
    seqlengths(regs) <- length(meanCov$coverage[[1]])
    
    ## Get coverage matrix by chunks of regions
    nChunks <- length(regs) %/% chunksize
    if(length(regs) %% chunksize > 0) nChunks <- nChunks + 1
    
    ## Split regions into chunks
    if(nChunks == 1) {
        regs_split <- list(regs)
    } else {
        regs_split <- split(regs, cut(seq_len(length(regs)), breaks = nChunks, labels = FALSE))
    }
    
    ## Actually calculate the coverage matrix
    resChunks <- bplapply(regs_split, .railMatChrRegion, sampleFiles = sampleFiles, chr = chr, mappedPerXM = mappedPerXM, L = L, returnBP = returnBP, verbose = verbose, BPPARAM = BPPARAM.railChr)
    
    
    ## Finish
    if(returnBP) {
        regionCov <- lapply(resChunks, '[[', 'bpCoverage')
        names(regionCov) <- NULL
        regionCov <- do.call(c, regionCov)
        
        result <- list(
            regions = regs,
            coverageMatrix = do.call(rbind, lapply(resChunks, '[[', 'coverageMatrix')),
            bpCoverage = regionCov
        )
    } else {
        result <- list(
            regions = regs,
            coverageMatrix = do.call(rbind, lapply(resChunks, '[[', 'coverageMatrix'))
        )
    }
    return(result)
}

.railMatChrRegion <- function(regions, sampleFiles, chr, mappedPerXM, L, returnBP, verbose = TRUE) {
        
    ## Get the region coverage matrix
    fullCov <- fullCoverage(files = sampleFiles, chrs = chr,
        protectWhich = 0, which = regions)
    
    ## Normalize coverage
    if(!is.null(mappedPerXM)) {
        if(verbose) message(paste(Sys.time(), 'railMatrix: normalizing coverage'))
        for(i in ncol(fullCov[[1]])) fullCov[[1]][[i]] <- fullCov[[1]][[i]] / mappedPerXM[i]
        if(verbose) message(paste(Sys.time(), 'railMatrix: done normalizing coverage'))
    }

    regionCov <- getRegionCoverage(fullCov = fullCov, regions = regions,
        mc.cores = 1L)

    if(verbose) message(paste(Sys.time(), "railMatrix: calculating coverageMatrix"))
    covMat <- lapply(regionCov, colSums)
    covMat <- do.call(rbind, covMat)

    if(verbose) message(paste(Sys.time(), "railMatrix: adjusting coverageMatrix for 'L'"))
    if(length(L) == 1) {
        covMat <- covMat / L
    } else if (length(L) == ncol(covMat)) {
        for(i in length(L)) covMat[, i] <- covMat[, i] / L[i]
    } else {
        warning("Invalid 'L' value so it won't be used. It has to either be a integer/numeric vector of length 1 or length equal to the number of samples.")
    }

    ## Finish
    if(returnBP) {
        res <- list(coverageMatrix = covMat, bpCoverage = regionCov)
    } else {
        res <- list(coverageMatrix = covMat)
    }
    return(res)
}
