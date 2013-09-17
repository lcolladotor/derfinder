#' Merge results from different chromosomes
#'
#' This function merges the results from running \link{analyzeChr} on several chromosomes and assigns genomic states using \link{annotateRegions}. It re-calculates the p-values and q-values using the pooled areas from the null regions from all chromosomes.
#' 
#' @param chrnums The chromosome numbers of the files to be merged.
#' @param prefix The main data directory path, which can be useful if \link{analyzeChr} is used for several parameters and the results are saved in different directories.
#' @param significantCut A vector of length two specifiying the cutoffs used to determine significance. The first element is used to determine significance for the p-values and the second element is used for the q-values just like in \link{calculatePvalues}.
#' @param genomicState This argument is passed to \link{annotateRegions}.
#' @param minoverlap This argument is passed to \link{annotateRegions}.
#' @param fullOrCoding This argument is passed to \link{annotateRegions}.
#' @param verbose If \code{TRUE} basic status updates will be printed along the way.
#'
#' @return Six Rdata files.
#' \describe{
#' \item{fullFstats.Rdata }{ Full F-statistics from all chromosomes in a list of Rle objects.}
#' \item{fullTime.Rdata }{ Timing information from all chromosomes.}
#' \item{fullNullSummary.Rdata}{ A DataFrame with the null region information: statistic, width, chromosome and permutation identifier. It's ordered by the statistics}
#' \item{fullRegions.Rdata}{ GRanges object with regions found and with full annotation from \link[bumphunter]{annotateNearest}. Note that the column \code{strand} from \link[bumphunter]{annotateNearest} is renamed to \code{annoStrand} to comply with GRanges specifications. }
#' \item{fullCoveragePrep.Rdata}{ A list with the pre-processed coverage data from all chromosomes.}
#' \item{fullAnnotatedRegions.Rdata}{ A list as constructed in \link{annotateRegions} with the assigned genomic states.}
#' }
#'
#' @author Leonardo Collado-Torres
#' @seealso \link{analyzeChr}, \link{calculatePvalues}, \link{annotateRegions}
#' @export
#' @importFrom GenomicRanges GRangesList
#' @importMethodsFrom GenomicRanges unlist
#' @importFrom IRanges DataFrame RleList
#' @importMethodsFrom IRanges cbind values "values<-" "[" "$" "$<-" length order unlist as.numeric
#' @importFrom qvalue qvalue
#'
#' @examples
#' \dontrun{
#' ## Hsapiens.UCSC.hg19.knownGene GenomicState
#' library("TxDb.Hsapiens.UCSC.hg19.knownGene")
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#'
#' ## Creating this GenomicState object takes around 8 min
#' GenomicState.Hsapiens.UCSC.hg19.knownGene <- makeGenomicState(txdb=txdb)
#'
#' ## Merge the results
#' mergeResults(prefix="run1", genomicState=GenomicState.Hsapiens.UCSC.hg19.knownGene)
#'
#' ## You can then explore the wallclock time spent on each step
#' prefix <- "." ## Set to the same prefix you used in mergeResults()
#' load(file.path(prefix, "fullRegions.Rdata"))
#' 
#' ## Process the time info
#' time <- lapply(fullTime, function(x) data.frame(diff(x)))
#' time <- do.call(rbind, time)
#' colnames(time) <- "sec"
#' time$sec <- as.integer(round(time$sec))
#' time$min <- time$sec / 60
#' time$chr <- paste0("chr", gsub("\\..*", "", rownames(time)))
#' time$step <- gsub(".*\\.", "", rownames(time))
#' rownames(time) <- seq_len(nrow(time))
#' 
#' ## Make plot
#' library("ggplot2")
#' ggplot(time, aes(x=step, y=min, colour=chr)) + geom_point() + labs(title="Wallclock time by step") + scale_colour_discrete(limits=chrs) + scale_x_discrete(limits=names(fullTime[[1]])[-1]) + ylab("Time (min)") + xlab("Step")
#' }

mergeResults <- function(chrnums=c(1:22, "X", "Y"), prefix=".", significantCut=c(0.05, 0.10), genomicState, minoverlap=20, fullOrCoding = "full", verbose=TRUE) {	
	## For R CMD check
	prep <- fstats <- regions <- annotation <- timeinfo <- NULL
	
	## Initialize
	fullCoveragePrep <- fullTime <- fullNullPermutation <- fullNullWidths <- fullNullStats <- fullFstats <- fullAnno <- fullRegs <- vector("list", length(chrnums))
	names(fullCoveragePrep) <- names(fullTime) <- names(fullNullPermutation) <- names(fullNullWidths) <- names(fullNullStats) <- names(fullFstats) <- names(fullAnno) <- names(fullRegs) <- paste0("chr", chrnums)

	## Actual processing
	for(current in chrnums) {
		chr <- paste0("chr", current)
		if(verbose) message(paste(Sys.time(), "Loading chromosome", current))
	
		## Process the F-statistics
		load(file.path(prefix, chr, "fstats.Rdata"))
		fullFstats[[chr]] <- fstats
	
		## Process the regions, nullstats and nullwidths
		load(file.path(prefix, chr, "regions.Rdata"))
		fullRegs[[chr]] <- regions$regions
		fullNullStats[[chr]] <- regions$nullStats
		fullNullWidths[[chr]] <- regions$nullWidths
		fullNullPermutation[[chr]] <- regions$nullPermutation
	
		## Process the annotation results
		load(file.path(prefix, chr, "annotation.Rdata"))
		fullAnno[[chr]] <- annotation
	
		## Process the timing information
		load(file.path(prefix, chr, "timeinfo.Rdata"))
		fullTime[[chr]] <- timeinfo
		
		## Process the covPrep data
		load(file.path(prefix, chr, "coveragePrep.Rdata"))
		fullCoveragePrep[[chr]] <- prep
	}

	## Save Fstats, Nullstats, and time info
	if(verbose) message(paste(Sys.time(), "mergeResults: Saving fullFstats"))
	save(fullFstats, file=file.path(prefix, "fullFstats.Rdata"))
	
	if(verbose) message(paste(Sys.time(), "mergeResults: Saving fullTime"))
	save(fullTime, file=file.path(prefix, "fullTime.Rdata"))
	
	if(verbose) message(paste(Sys.time(), "mergeResults: Saving fullCoveragePrep"))
	save(fullCoveragePrep, file=file.path(prefix, "fullCoveragePrep.Rdata"))
	
	## Process the annotation 
	fullAnnotation <- do.call(rbind, fullAnno)
	colnames(fullAnnotation)[which(colnames(fullAnnotation) == "strand")] <- "annoStrand"

	## Combine regions with annotation
	fullRegions <- unlist(GRangesList(fullRegs), use.names=FALSE)
    values(fullRegions) <- cbind( values(fullRegions), DataFrame(fullAnnotation))	
	
	## Re-calculate p-values and q-values
	if(verbose) message(paste(Sys.time(), "mergeResults: Re-calculating the p-values"))
		
	## Summarize the null regions
	nulls <- unlist(RleList(fullNullStats), use.names=FALSE)
	widths <- unlist(RleList(fullNullWidths), use.names=FALSE)
	permutations <- unlist(RleList(fullNullPermutation), use.names=FALSE)
	howMany <- unlist(lapply(fullNullStats, length))
		
	if(length(nulls) > 0) {
		## Proceed only if there are null regions to work with
		fullNullSummary <- DataFrame(stat=nulls, width=widths, chr=Rle(names(fullNullStats), howMany), permutation=permutations)
		rm(nulls, widths, howMany, permutations)
		fullNullSummary$area <- fullNullSummary$stat * fullNullSummary$width
		fullNullSummary <- fullNullSummary[order(fullNullSummary$area, decreasing=TRUE), ]
	} else {
		fullNullSummary <- DataFrame(NULL)	
	}
	
	if(verbose) message(paste(Sys.time(), "mergeResults: Saving fullNullSummary"))
	save(fullNullSummary, file=file.path(prefix, "fullNullSummary.Rdata"))
		
	if(nrow(fullNullSummary) > 0) {
		## Actual calculation
		nullareas <- as.numeric(abs(fullNullSummary$stat) * fullNullSummary$width)
		pvals <- sapply(fullRegions$area, function(x) { sum(nullareas > x) })
	
		## Update info
		fullRegions$pvalues <- (pvals + 1) / (length(nullareas) + 1)
		fullRegions$qvalues <- qvalue(fullRegions$pvalues)$qvalues
		fullRegions$significant <- factor(fullRegions$pvalues < significantCut[1], levels=c(TRUE, FALSE))
		fullRegions$significantQval <- factor(fullRegions$qvalues < significantCut[2], levels=c(TRUE, FALSE))
	}	
	## Sort by decreasing area
	fullRegions <- fullRegions[order(fullRegions$area, decreasing=TRUE), ]
	
	## save GRanges version
	if(verbose) message(paste(Sys.time(), "mergeResults: Saving fullRegions"))
	save(fullRegions, file=file.path(prefix, "fullRegions.Rdata"))
	
	## Assign genomic states
	if(verbose) message(paste(Sys.time(), "mergeResults: assigning genomic states"))
	fullAnnotatedRegions <- annotateRegions(regions=fullRegions, genomicState=genomicState, minoverlap=minoverlap, fullOrCoding = fullOrCoding, annotate=TRUE, verbose=verbose)
	
	if(verbose) message(paste(Sys.time(), "mergeResults: Saving fullAnnotatedRegions"))
	save(fullAnnotatedRegions, file=file.path(prefix, "fullAnnotatedRegions.Rdata"))		
	
	## Finish
	return(invisible(NULL))
}


