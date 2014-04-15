#' Extract coverage information for exons
#'
#' This function extracts the coverage information calculated by 
#' \link{fullCoverage} for a set of exons determined by 
#' \link{makeGenomicState}. The underlying code is similar to 
#' \link{getRegionCoverage} with additional tweaks for calculating RPKM values.
#' 
#' @param fullCov A list where each element is the result from 
#' \link{loadCoverage} used with \code{cutoff=NULL}. The elements of the list 
#' should be named according to the chromosome number. Can be generated using 
#' \link{fullCoverage}.
#' @param genomicState The output from \link{makeGenomicState}.
#' @param fullOrCoding If \code{full} then the \code{genomicState$fullGenome} 
#' genomic state information is used. If \code{coding}, then the 
#' \code{genomicState$codingGenome} genomic state information is used.
#' @param L The width of the reads used.
#' @param returnType If \code{raw}, then the raw coverage information per exon 
#' is returned. If \code{rpkm}, RPKM values are calculated for each exon.
#' @param mc.cores This argument is passed to \link[parallel]{mclapply} twice. 
#' First, it's is used by strand. Secondly, for processing the exons by 
#' chromosome. So there is no gain in using \code{mc.cores} greater than the 
#' maximum of the number of strands and number of chromosomes.
#' @param verbose If \code{TRUE} basic status updates will be printed along the 
#' way.
#'
#' @return A matrix (nrow = number of exons in \code{genomicState} 
#' corresponding to the chromosomes in \code{fullCov}, ncol = number of 
#' samples) with the number of reads (or RPKM) per exon. The row names 
#' correspond to the row indexes of \code{genomicState$fullGenome}  (if 
#' \code{fullOrCoding='full'}) or \code{genomicState$codingGenome} (if 
#' \code{fullOrCoding='coding'}).
#'
#' @author Andrew Jaffe, Leonardo Collado-Torres
#' @seealso \link{fullCoverage}, \link{getRegionCoverage}
#' @export
#' @importFrom GenomicRanges seqlevels seqnames
#' @importMethodsFrom GenomicRanges names 'names<-' length '[' coverage sort 
#' width c strand subset as.data.frame
#' @importMethodsFrom IRanges subset as.data.frame as.character runValue '%in%'
#' @importFrom parallel mclapply
#'
#' @examples
#' ## Obtain fullCov object
#' fullCov <- list('21'=genomeDataRaw$coverage)
#'
#' ## Use only the first two exons
#' smallGenomicState <- genomicState
#' smallGenomicState$fullGenome <- smallGenomicState$fullGenome[ 
#'     which(smallGenomicState$fullGenome$theRegion == 'exon')[1:2] ]
#' 
#' ## Finally, get the coverage information for each exon
#' exonCov <- coverageToExon(fullCov=fullCov, genomicState=smallGenomicState, 
#'     L=36)


coverageToExon <- function(fullCov, genomicState, fullOrCoding = "full", 
    L = NULL, returnType = "raw", mc.cores = getOption("mc.cores", 
        2L), verbose = TRUE) {
    stopifnot(length(intersect(names(genomicState), c("fullGenome", 
        "codingGenome"))) == 2)
    stopifnot(length(intersect(fullOrCoding, c("full", "coding"))) == 
        1)
    stopifnot(length(intersect(returnType, c("raw", "rpkm"))) == 
        1)
    if (is.null(L)) 
        stop("'L' has to be specified")
    
    if (fullOrCoding == "full") {
        gs <- genomicState$fullGenome
    } else if (fullOrCoding == "coding") {
        gs <- genomicState$codingGenome
    }
    
    # just the reduced exons
    etab <- gs[gs$theRegion == "exon"]
    
    ## Keep only the exons from the chromosomes in fullCov
    etab <- etab[seqnames(etab) %in% paste0("chr", names(fullCov))]
    
    # split by strand
    strandIndexes <- split(seq_len(length(etab)), as.character(strand(etab)))
    
    # count reads covering exon on each strand
    nCores <- min(mc.cores, length(unique(runValue(strand(etab)))))
    # Use at most n cores where n is the number of unique strands
    
    exonByStrand <- mclapply(strandIndexes, .coverageToExonStrandStep, 
        fullCov = fullCov, etab = etab, L = L, nCores = nCores, 
        verbose = verbose, mc.cores = nCores)
    
    # combine two strands
    exons <- do.call("rbind", exonByStrand)
    
    # put back in annotation order
    theExons <- exons[names(etab), ]
    
    if (returnType == "rpkm") {
        Mchr <- t(sapply(fullCov, function(z) sapply(z, function(xx)
            sum(as.numeric(runValue(xx))))))
        M <- colSums(Mchr)/L/1e+06
        theExons <- theExons/(width(etab)/1000)/M
    }
    return(theExons)
}

.coverageToExonStrandStep <- function(ii, fullCov, etab, L, nCores, 
    verbose) {
    e <- etab[ii]  # subset
    
    ## use logical rle to subset large coverage matrix
    cc <- coverage(e)  # first coverage
    for (i in seq(along = cc)) {
        # then convert to logical
        cc[[i]]@values <- ifelse(cc[[i]]@values > 0, TRUE, FALSE)
    }
    
    # now count exons
    exonList <- mclapply(seq(along = fullCov), .coverageToExonChrStep, 
        fullCov = fullCov, cc = cc, e = e, L = L, verbose = verbose, 
        mc.cores = nCores)
    
    out <- do.call("rbind", exonList)  # combine
    
    ## Clean up
    rm(e, cc, exonList)
    gc()
    return(out)
}


.coverageToExonChrStep <- function(i, fullCov, cc, e, L, verbose = verbose) {
    chrnum <- names(fullCov)[i]
    chr <- paste0("chr", chrnum)
    if (verbose) 
        message(paste(Sys.time(), "coverageToExon: processing chromosome", 
            chrnum))
    
    # subset using logical rle (fastest way)
    z.tmp <- subset(fullCov[[chrnum]], cc[[chr]])
    z <- as.data.frame(z.tmp)
    
    # only exons from this chr
    g <- e[seqnames(e) == chr]
    ind <- rep(names(g), width(g))  # to split
    tmpList <- split(z, ind)  # split
    res <- t(sapply(tmpList, colSums)/L)  # get # reads
    
    ## Clean up
    rm(z, g, ind, tmpList)
    gc()
    return(res)
} 
