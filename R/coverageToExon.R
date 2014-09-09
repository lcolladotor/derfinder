#' Extract coverage information for exons
#'
#' This function extracts the coverage information calculated by 
#' \link{fullCoverage} for a set of exons determined by 
#' \link{makeGenomicState}. The underlying code is similar to 
#' \link{getRegionCoverage} with additional tweaks for calculating RPKM values.
#' 
#' @param fullCov A list where each element is the result from 
#' \link{loadCoverage} used with \code{cutoff=NULL}. Can be generated using 
#' \link{fullCoverage}.
#' @param genomicState The output from \link{makeGenomicState}.
#' @param fullOrCoding If \code{full} then the \code{genomicState$fullGenome} 
#' genomic state information is used. If \code{coding}, then the 
#' \code{genomicState$codingGenome} genomic state information is used.
#' @param L The width of the reads used.
#' @param returnType If \code{raw}, then the raw coverage information per exon 
#' is returned. If \code{rpkm}, RPKM values are calculated for each exon.
#' @param mc.cores This argument is passed to \link[BiocParallel]{SnowParam} 
#' twice. First, it is used by strand. Second, for processing the exons by 
#' chromosome. So there is no gain in using \code{mc.cores} greater than the 
#' maximum of the number of strands and number of chromosomes.
#' @param mc.outfile This argument is passed to \link[BiocParallel]{SnowParam} 
#' to specify the \code{outfile} for any output from the workers.
#' @param chrsStyle The naming style of the chromosomes. By default, UCSC. See 
#' \link[GenomeInfoDb]{seqlevelsStyle}.
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
#' @aliases coverage_to_exon
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomeInfoDb seqlevels seqlevelsStyle 'seqlevelsStyle<-'
#' mapSeqlevels
#' @importMethodsFrom GenomicRanges names 'names<-' length '[' coverage sort 
#' width c strand subset as.data.frame
#' @importMethodsFrom IRanges subset as.data.frame
#' @importMethodsFrom S4Vectors as.character '%in%'
#' @importFrom S4Vectors runValue
#' @importFrom BiocParallel SnowParam SerialParam bplapply bpmapply
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
    L = NULL, returnType = "raw", mc.cores = getOption("mc.cores", 1L),
    mc.outfile = Sys.getenv('SGE_STDERR_PATH'), chrsStyle = "UCSC",
    verbose = TRUE) {
    
    ## Run some checks
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
    
    ## Use UCSC names by default
    seqlevelsStyle(gs) <- chrsStyle
    names(fullCov) <- mapSeqlevels(names(fullCov), chrsStyle)
    
    # just the reduced exons
    etab <- gs[gs$theRegion == "exon"]
    
    ## Check that the names are unique
    stopifnot(length(etab) == length(unique(names(etab))))
    
    ## Keep only the exons from the chromosomes in fullCov
    chrKeep <- names(fullCov)
    etab <- etab[seqnames(etab) %in% chrKeep]
    
    # split by strand
    strandIndexes <- split(seq_len(length(etab)), as.character(strand(etab)))
    
    # count reads covering exon on each strand
    nCores <- min(mc.cores, length(unique(runValue(strand(etab)))))

    ## Define cluster
    if(nCores > 1) {
        BPPARAM <- SnowParam(workers = nCores, outfile = mc.outfile)
    } else {
        BPPARAM <- SerialParam()
    }

    # Use at most n cores where n is the number of unique strands
    exonByStrand <- bplapply(strandIndexes, .coverageToExonStrandStep, 
        fullCov = fullCov, etab = etab, L = L, nCores = mc.cores,
        mcOut = mc.outfile, chrs = chrKeep,
        verbose = verbose, BPPARAM = BPPARAM)
    
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

.coverageToExonStrandStep <- function(ii, fullCov, etab, L, nCores, mcOut, chrs,
    verbose) {
        
    e <- etab[ii]  # subset
    
    ## use logical rle to subset large coverage matrix
    cc <- coverage(e)  # first coverage
    for (i in seq(along = cc)) {
        # then convert to logical
        cc[[i]]@values <- ifelse(cc[[i]]@values > 0, TRUE, FALSE)
    }
    
    ## Subset data
    # subset using logical rle (fastest way)
    subsets <- mapply(function(covInfo, chr) { subset(covInfo, cc[[chr]]) },
        fullCov, chrs, SIMPLIFY = FALSE)
    
    # now count exons
    moreArgs <- list(e = e, L = L, verbose = verbose)
    
    ## Define cluster
    nCores <- min(nCores, length(subsets))
    if(nCores > 1) {
        BPPARAM.chrStep <- SnowParam(workers = nCores, outfile = mcOut)
    } else {
        BPPARAM.chrStep <- SerialParam()
    }
    
    ## Define ChrStep function
    .coverageToExonChrStep <- function(z.DF, chr, e, L, verbose) {
        if (verbose) 
            message(paste(Sys.time(), "coverageToExon: processing chromosome", chr))
    
        ## Transform to regular data.frame   
        z <- as.data.frame(z.DF)
    
        # only exons from this chr
        g <- e[seqnames(e) == chr]
        ind <- rep(names(g), width(g))  # to split
        tmpList <- split(z, ind)  # split
        res <- t(sapply(tmpList, colSums)/L)  # get # reads
    
        # done
        return(res)
    }    
    
    ## Now run it
    exonList <- bpmapply(.coverageToExonChrStep, subsets, chrs, 
        MoreArgs = moreArgs, BPPARAM = BPPARAM.chrStep, SIMPLIFY = FALSE)
        
    # combine
    out <- do.call("rbind", exonList)
    
    # done
    return(out)
}

#' @export
coverage_to_exon <- coverageToExon
