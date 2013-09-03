#' Merge results from different chromosomes
#'
#' This function merges the results from running \link{analyzeChr} on several chromosomes. It re-calculates the p-values and q-values using the pooled areas from the null regions from all chromosomes.
#' 
#' @param chrnums The chromosome numbers of the files to be merged.
#' @param prefix The main data directory path, which can be useful if \link{analyzeChr} is used for several parameters and the results are saved in different directories.
#' @param significantCut A vector of length two specifiying the cutoffs used to determine significance. The first element is used to determine significance for the p-values and the second element is used for the q-values just like in \link{calculatePvalues}.
#' @param verbose If \code{TRUE} basic status updates will be printed along the way.
#'
#' @return Six Rdata files.
#' \describe{
#' \item{fullFstats.Rdata }{ Full F-statistics from all chromosomes.}
#' \item{fullTime.Rdata }{ Timing information from all chromosomes.}
#' \item{fullNullstats.Rdata }{ Null F-statistics from all chromosomes.}
#' \item{fullNullwidths.Rdata }{ Null region widths from all chromosomes.}
#' \item{fullNullchrs.Rdata}{ Null region chromosome identifier.}
#' \item{fullRegions.Rdata}{ Regions found with full annotation from \link[bumphunter]{annotateNearest}. Note that the column \code{strand} from \link[bumphunter]{annotateNearest} is renamed to \code{annoStrand} to comply with GRanges specifications. }
#' }
#'
#' @author Leonardo Collado-Torres
#' @seealso \link{analyzeChr}, \link{calculatePvalues}
#' @export
#' @importFrom GenomicRanges GRangesList
#' @importMethodsFrom GenomicRanges unlist
#' @importFrom IRanges DataFrame
#' @importMethodsFrom IRanges cbind values "values<-" "[" "$" "$<-"
#' @importFrom qvalue qvalue
#'
#' @examples
#' \dontrun{
#' mergeResults(prefix="run1")
#' }

mergeResults <- function(chrnums=c(1:22, "X", "Y"), prefix=".", significantCut=c(0.05, 0.10), verbose=TRUE) {
	library("IRanges")
	library("GenomicRanges")
	
	## Initialize
	fullTime <- fullNullwidths <- fullNullstats <- fullFstats <- fullAnno <- fullRegs <- vector("list", length(chrnums))
	names(fullTime) <- names(fullNullwidths) <- names(fullNullstats) <- names(fullFstats) <- names(fullAnno) <- names(fullRegs) <- chrnums

	## Actual processing
	for(current in chrnums) {
		chr <- paste0("chr", current)
		if(verbose) message(paste(Sys.time(), "Loading chromosome", current))
	
		## Process the F-statistics
		load(file.path(prefix, chr, "fstats.Rdata"))
		fullFstats[[chr]] <- fstats
	
		## Process the regions, nullstats and nullwidths
		load(file.path(prefix, chr, "regions.Rdata"))
		fullRegs[[chr]] <- regs$regions
		fullNullstats[[chr]] <- regs$nullstats
		fullNullwidths[[chr]] <- regs$nullwidths
	
		## Process the annotation results
		load(file.path(prefix, chr, "annotation.Rdata"))
		fullAnno[[chr]] <- annotation
	
		## Process the timing information
		load(file.path(prefix, chr, "timeinfo.Rdata"))
		fullTime[[chr]] <- timeinfo
	}

	## Save Fstats, Nullstats, and time info
	if(verbose) message(paste(Sys.time(), "mergeResults: Saving fullFstats"))
	save(fullFstats, file=file.path(prefix, "fullFstats.Rdata"))
	
	if(verbose) message(paste(Sys.time(), "mergeResults: Saving fullTime"))
	save(fullTime, file=file.path(prefix, "fullTime.Rdata"))
	
	if(verbose) message(paste(Sys.time(), "mergeResults: Saving fullNullstats"))
	save(fullNullstats, file=file.path(prefix, "fullNullstats.Rdata"))
	
	if(verbose) message(paste(Sys.time(), "mergeResults: Saving fullNullwidths"))
	save(fullNullwidths, file=file.path(prefix, "fullNullwidths.Rdata"))
	
	if(verbose) message(paste(Sys.time(), "mergeResults: Saving fullNullchrs"))
	save(fullNullchrs, file=file.path(prefix, "fullNullchrs.Rdata"))

	## Process the annotation 
	fullAnnotation <- do.call(rbind, fullAnno)
	colnames(fullAnnotation)[which(colnames(fullAnnotation) == "strand")] <- "annoStrand"

	## Combine regions with annotation
	fullRegions <- unlist(GRangesList(fullRegs), use.names=FALSE)
    values(fullRegions) <- cbind( values(fullRegions), DataFrame(fullAnnotation))	
	
	## Re-calculate p-values and q-values
	if(verbose) message(paste(Sys.time(), "mergeResults: Re-calculating the p-values"))
		
	## Actual calculation
	nullareas <- as.numeric(fullNullstats * fullNullwidths)
	pvals <- sapply(fullRegions$area, function(x) { sum(nullareas > x) })
	
	## Update info
	fullRegions$pvalues <- (pvals + 1) / (length(nullareas) + 1)
	fullRegions$qvalues <- qvalue(fullRegions$pvalues)$qvalues
	fullRegions$significant <- factor(fullRegions$pvalues < significantCut[1], levels=c(TRUE, FALSE))
	fullRegions$significantQval <- factor(fullRegions$qvalues < significantCut[2], levels=c(TRUE, FALSE))
	
	## Sort by decreasing area
	fullRegions <- fullRegions[order(fullRegions$area, decreasing=TRUE), ]

	## save GRanges version
	if(verbose) message(paste(Sys.time(), "mergeResults: Saving fullRegions"))
	save(fullRegions, file=file.path(prefix, "fullRegions.Rdata"))
	
	## Finish
	return(invisible(NULL))
}


