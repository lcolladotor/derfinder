#' Extract coverage information for exons
#'
#' This function extracts the coverage information calculated by 
#' \link{fullCoverage} for a set of exons determined by 
#' \link{makeGenomicState}. The underlying code is similar to 
#' \link{getRegionCoverage} with additional tweaks for calculating RPKM values.
#' 
#' @inheritParams getRegionCoverage
#' @inheritParams annotateRegions
#' @param L The width of the reads used.
#' @param returnType If \code{raw}, then the raw coverage information per exon 
#' is returned. If \code{rpkm}, RPKM values are calculated for each exon.
#' @inheritParams fullCoverage
#' @param ... Arguments passed to other methods and/or advanced arguments.
#'
#' @return A matrix (nrow = number of exons in \code{genomicState} 
#' corresponding to the chromosomes in \code{fullCov}, ncol = number of 
#' samples) with the number of reads (or RPKM) per exon. The row names 
#' correspond to the row indexes of \code{genomicState$fullGenome}  (if 
#' \code{fullOrCoding='full'}) or \code{genomicState$codingGenome} (if 
#' \code{fullOrCoding='coding'}).
#'
#' @details
#' Parallelization is used twice.
#' First, it is used by strand. Second, for processing the exons by 
#' chromosome. So there is no gain in using \code{mc.cores} greater than the 
#' maximum of the number of strands and number of chromosomes.
#'
#' If \code{fullCov} is \code{NULL} and \code{files} is specified, this function
#' will attempt to read the coverage from the files. Note that if you used
#' 'totalMapped' and 'targetSize' before, you will have to specify them again
#' to get the same results. 
#'
#' See also \link{advanedArg} with \code{fun='loadCoverage'} for other details.
#'
#' @author Andrew Jaffe, Leonardo Collado-Torres
#' @seealso \link{fullCoverage}, \link{getRegionCoverage}
#' @export
#' @aliases coverage_to_exon
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomeInfoDb seqlevels seqlevelsStyle 'seqlevelsStyle<-'
#' mapSeqlevels seqlevelsInUse
#' @importMethodsFrom GenomicRanges names 'names<-' length '[' coverage sort 
#' width c strand subset as.data.frame
#' @importMethodsFrom IRanges subset as.data.frame
#' @importMethodsFrom S4Vectors as.character '%in%'
#' @importFrom S4Vectors runValue
#' @importFrom BiocParallel bplapply bpmapply
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
#' exonCov <- coverageToExon(fullCov=fullCov,
#'     genomicState=smallGenomicState$fullGenome, L=36)
#'


coverageToExon <- function(fullCov = NULL, genomicState, L = NULL, 
    returnType = "raw", files = NULL, ...) {
        
    ## Run some checks
    stopifnot(length(intersect(returnType, c("raw", "rpkm"))) == 
        1)
    stopifnot(is(genomicState, 'GRanges'))
    stopifnot(identical(names(mcols(genomicState)), c('theRegion', 'tx_id',
        'tx_name', 'gene')))

    if (is.null(L)) 
        stop("'L' has to be specified")
    
    ## Advanged argumentsa
#' @param chrsStyle The naming style of the chromosomes. By default, UCSC. See 
#' \link[GenomeInfoDb]{seqlevelsStyle}.    
    chrsStyle <- .advanced_argument('chrsStyle', 'UCSC', ...)


#' @param verbose If \code{TRUE} basic status updates will be printed along the 
#' way.
    verbose <- .advanced_argument('verbose', TRUE, ...)


    ## Use UCSC names by default
    seqlevelsStyle(genomicState) <- chrsStyle
    
    # just the reduced exons
    etab <- genomicState[genomicState$theRegion == "exon"]
    
    ## Load data if 'fullCov' is not specified
    if(is.null(fullCov)) {
        fullCov <- .load_fullCov(files = files, chrs = seqlevelsInUse(etab),
            fun = 'coverageToExon', verbose = verbose, ...)        
    }
    ## Fix naming style
    names(fullCov) <- mapSeqlevels(names(fullCov), chrsStyle)
    
    
    
    ## Check that the names are unique
    stopifnot(length(etab) == length(unique(names(etab))))
    
    ## Keep only the exons from the chromosomes in fullCov
    chrKeep <- names(fullCov)
    etab <- etab[seqnames(etab) %in% chrKeep]
    
    # split by strand
    strandIndexes <- split(seq_len(length(etab)), as.character(strand(etab)))
    
    # count reads covering exon on each strand
    strandCores <- min(.advanced_argument('mc.cores', getOption('mc.cores', 1L),
        ...), length(unique(runValue(strand(etab)))))

    ## Define cluster
    BPPARAM <- .advanced_argument('BPPARAM',
        .define_cluster(cores = 'strandCores', strandCores), ...)

    # Use at most n cores where n is the number of unique strands
    exonByStrand <- bplapply(strandIndexes, .coverageToExonStrandStep, 
        fullCov = fullCov, etab = etab, L = L,
        nCores = .advanced_argument('mc.cores', getOption('mc.cores', 1L), ...),
        chrs = chrKeep, verbose = verbose, BPPARAM = BPPARAM, ...)
    
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

.coverageToExonStrandStep <- function(ii, fullCov, etab, L, nCores, chrs,
    verbose, ...) {
        
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
    exonCores <- min(nCores, length(subsets))

    ## Define cluster
    BPPARAM.chrStep <- .advanced_argument('BPPARAM.chrStep', 
        .define_cluster(cores = 'exonCores', exonCores), ...)
    
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
