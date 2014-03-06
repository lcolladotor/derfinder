#' Assign genomic states to regions
#'
#' This function takes the regions found in \link{calculatePvalues} and assigns them genomic states contructed with \link{makeGenomicState}. The main workhorse functions are \link[IRanges]{countOverlaps} and \link[IRanges]{findOverlaps}.
#' 
#' @param regions The \code{$regions} output from \link{calculatePvalues}.
#' @param genomicState The output from \link{makeGenomicState}.
#' @param minoverlap This parameter is passed to \link[IRanges]{findOverlaps} and determines the minimum overlap of a region to be assined a genomic state. Set to 1 if you want all the genomic states that overlap the regions.
#' @param fullOrCoding If \code{full} then the \code{genomicState$fullGenome} genomic state information is used. If \code{coding}, then the \code{genomicState$codingGenome} genomic state information is used.
#' @param annotate If \code{TRUE} then the regions are annotated by the genomic state. Othewise, only the overlaps between the regions and the genomic states are computed.
#' @param verbose If \code{TRUE} basic status updates will be printed along the way.
#'
#' @return A list with elements \code{countTable} and \code{annotationList} (only if \code{annotate=TRUE}). 
#' \describe{
#' \item{countTable }{This is a data.frame with the number of overlaps from the regions vs the genomic states with one type per column. For example, if \code{fullOrCoding="full"} then the columns are \code{exon}, \code{intragenic} and \code{intron}.}
#' \item{annotationList }{This is a \code{GRangesList} with the genomic states that overlapped with the regions (if any, depends on \code{minoverlap}). The names of this \code{GRangesList} correspond to the region index in \code{regions}.}
#' }
#'
#' @author Andrew Jaffe, Leonardo Collado-Torres
#' @seealso \link{makeGenomicState}, \link{calculatePvalues}
#' @export
#' @importFrom IRanges queryHits subjectHits
#' @importMethodsFrom GenomicRanges names "names<-" length "$" split countOverlaps findOverlaps "["
#' @importMethodsFrom IRanges sapply
#'
#' @examples
#' \dontrun{
#' ## Create GenomicState object:
#' ## Hsapiens.UCSC.hg19.knownGene GenomicState
#' library("TxDb.Hsapiens.UCSC.hg19.knownGene")
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#'
#' ## Creating this GenomicState object takes around 8 min
#' GenomicState.Hsapiens.UCSC.hg19.knownGene <- makeGenomicState(txdb=txdb)
#' 
#' ## Annotate regions
#' annotatedRegions <- annotateRegions(regions=genomeRegions$regions, genomicState=GenomicState.Hsapiens.UCSC.hg19.knownGene, minoverlap=1)
#' }

annotateRegions <- function(regions, genomicState, minoverlap=20, fullOrCoding = "full", annotate=TRUE, verbose=TRUE) {
	stopifnot(length(intersect(names(genomicState), c("fullGenome", "codingGenome"))) == 2)
	stopifnot(length(intersect(fullOrCoding, c("full", "coding"))) == 1)
	
	## Fix row names
	names(regions) <- seq_len(length(regions))
	
	if(fullOrCoding == "full") {
		gs <- genomicState$fullGenome
	} else if (fullOrCoding == "coding") {
		gs <- genomicState$codingGenome
	}
	gsl <- split(gs, gs$theRegion)

	if(verbose) message(paste(Sys.time(), "annotateRegions: counting"))

	countTable <- sapply(gsl, function(x) countOverlaps(regions, x, minoverlap=minoverlap))
	countTable <- data.frame(countTable)
	out <- list(countTable=countTable)
	
	if(annotate) {
		if(verbose) message(paste(Sys.time(), "annotateRegions: annotating"))

		oo <- findOverlaps(regions, gs, minoverlap=minoverlap)
		glist <- split(gs[subjectHits(oo)], queryHits(oo))
		out$annotationList <- glist
	}
	return(out)
}