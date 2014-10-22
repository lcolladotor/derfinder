#' Assign genomic states to regions
#'
#' This function takes the regions found in \link{calculatePvalues} and assigns 
#' them genomic states contructed with \link{makeGenomicState}. The main 
#' workhorse functions are \link[IRanges]{countOverlaps} and 
#' \link[IRanges]{findOverlaps}.
#' 
#' @param regions The \code{$regions} output from \link{calculatePvalues}.
#' @param genomicState A GRanges object created with \link{makeGenomicState}. 
#' It can be either the \code{genomicState$fullGenome} or 
#' \code{genomicState$codingGenome} component.
#' @param annotate If \code{TRUE} then the regions are annotated by the genomic 
#' state. Otherwise, only the overlaps between the regions and the genomic 
#' states are computed.
#' @param ... Arguments passed to other methods and/or advanced arguments.
#'
#' @return A list with elements \code{countTable} and \code{annotationList} 
#' (only if \code{annotate=TRUE}). 
#' \describe{
#' \item{countTable }{This is a data.frame with the number of overlaps from the 
#' regions vs the genomic states with one type per column. For example, if 
#' \code{fullOrCoding='full'} then the columns are \code{exon}, 
#' \code{intragenic} and \code{intron}.}
#' \item{annotationList }{This is a \code{GRangesList} with the genomic states 
#' that overlapped with the regions (if any, depends on \code{minoverlap}). The 
#' names of this \code{GRangesList} correspond to the region index in 
#' \code{regions}.}
#' }
#'
#' @author Andrew Jaffe, Leonardo Collado-Torres
#' @seealso \link{makeGenomicState}, \link{calculatePvalues}
#' @export
#' @aliases annotate_regions
#' @importFrom IRanges queryHits subjectHits
#' @importFrom GenomeInfoDb renameSeqlevels seqlevels
#' @importMethodsFrom GenomicRanges names 'names<-' length '$' split mcols
#' countOverlaps findOverlaps '['
#' @importMethodsFrom S4Vectors sapply
#'
#' @details
#' You might want to specify arguments such as \code{minoverlap} to control
#' how the overlaps are determined. See \link[IRanges]{findOverlaps} for further
#' details.
#'
#' @examples
#' ## Annotate regions, first two only
#' annotatedRegions <- annotateRegions(regions=genomeRegions$regions[1:2], 
#'     genomicState=genomicState$fullGenome, minoverlap=1)
#' annotatedRegions

annotateRegions <- function(regions, genomicState, annotate = TRUE, ...) {
    stopifnot(is(genomicState, 'GRanges'))
    stopifnot(identical(names(mcols(genomicState)), c('theRegion', 'tx_id',
        'tx_name', 'gene')))

    ## Advanged arguments

#' @param verbose If \code{TRUE} basic status updates will be printed along the 
#' way.
    verbose <- .advanced_argument('verbose', TRUE, ...)


    ## Fix row names
    names(regions) <- seq_len(length(regions))

    ## Use UCSC names for homo_sapiens by default
    genomicState <- renameSeqlevels(genomicState,
        extendedMapSeqlevels(seqlevels(genomicState), ...))
    regions <- renameSeqlevels(regions,
        extendedMapSeqlevels(seqlevels(regions), ...))
    
    genomicState.list <- split(genomicState, genomicState$theRegion)
    
    if (verbose) 
        message(paste(Sys.time(), 'annotateRegions: counting'))
    
    countTable <- sapply(genomicState.list, function(x, ...) {
        .runFunFormal(countOverlaps, query = regions, subject = x, ...)
    }, ...)
    countTable <- data.frame(countTable)
    out <- list(countTable = countTable)
    
    if (annotate) {
        if (verbose) 
            message(paste(Sys.time(), 'annotateRegions: annotating'))
        
        oo <- findOverlaps(regions, genomicState, minoverlap = 1L)  # don't care about min overlap here
        glist <- split(genomicState[subjectHits(oo)], queryHits(oo))
        out$annotationList <- glist
    }
    return(out)
} 

#' @export
annotate_regions <- annotateRegions
