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
#' ## Collapse the coverage information
#' collapsedFull <- collapseFullCoverage(list(genomeData$coverage), verbose=TRUE)
#' 
#' ## Calculate library size adjustments
#' sampleDepths <- sampleDepth(collapsedFull, probs=c(0.5), nonzero=TRUE, verbose=TRUE)
#' 
#' ## Build the models
#' group <- genomeInfo$pop
#' adjustvars <- data.frame(genomeInfo$gender)
#' models <- makeModels(sampleDepths, testvars=group, adjustvars=adjustvars)
#'
#' ## Preprocess the data
#' ## Automatic chunksize used to then compare 1 vs 4 cores in the 'do not run' section
#' prep <- preprocessCoverage(genomeData, groupInfo=group, cutoff=0, scalefac=32, chunksize=NULL, colsubset=NULL, mc.cores=4)
#' 
#' ## Get the F statistics
#' fstats <- calculateStats(prep, models, mc.cores=1, verbose=TRUE)
#'
#' ## Determine a cutoff from the F-distribution.
#' ## This step is very important and you should consider using quantiles from the observed F statistics
#' n <- dim(prep$coverageProcessed)[2]
#' df1 <- dim(models$mod)[2]
#' df0 <- dim(models$mod0)[2]
#' cutoff <- qf(0.95, df1-df0, n-df1)
#' 
#' ## Low cutoff used for illustrative purposes
#' cutoff <- 1
#'
#' ## Calculate the p-values and define the regions of interest.
#' regsWithP <- calculatePvalues(prep, models, fstats, nPermute=10, seeds=NULL, chr="chr21", cutoff=cutoff, mc.cores=1)
#'
#' ## Create GenomicState object:
#' ## Hsapiens.UCSC.hg19.knownGene GenomicState
#' library("TxDb.Hsapiens.UCSC.hg19.knownGene")
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#'
#' ## Creating this GenomicState object takes around 8 min
#' GenomicState.Hsapiens.UCSC.hg19.knownGene <- makeGenomicState(txdb=txdb)
#' 
#' ## Annotate regions
#' annotatedRegions <- annotateRegions(regions=regsWithP$regions, genomicState=GenomicState.Hsapiens.UCSC.hg19.knownGene, minoverlap=1)
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