#' Extract coverage information for a set of regions
#'
#' This function extracts the raw coverage information calculated by \link{fullCoverage} at each base for a set of regions found with \link{calculatePvalues}. It can further calculate the mean coverage per sample for each region.
#' 
#' @param fullCov A list where each element is the result from \link{loadCoverage} used with \code{cutoff=NULL}. The elements of the list should be named according to the chromosome number. Can be generated using \link{fullCoverage}.
#' @param regions The \code{$regions} output from \link{calculatePvalues}. It is important that the seqlengths information is provided.
#' @param totalMapped The total number of reads mapped for each sample. Providing this data returns coverage adjusted totalMapped per million at each base.
#' @param mc.cores The number of cores to use for computing coverage. Default=1
#' @param verbose If \code{TRUE} basic status updates will be printed along the way.
#'
#' @return a list of data.frame where each data.frame has the coverage information (nrow = width of region, ncol = number of samples) for a given region. The names of the list correspond to the region indexes in \code{regions}
#' #'
#' @author Andrew Jaffe, Leonardo Collado-Torres
#' @seealso \link{fullCoverage}, \link{calculatePvalues}
#' @export
#' @importFrom GenomicRanges seqlevels seqnames seqlengths
#' @importMethodsFrom GenomicRanges names "names<-" length "[" coverage sort width c
#' @importMethodsFrom IRanges subset as.data.frame
#'
#' @examples
#' \dontrun{
#' ## Obtain fullCov object
#' datadir <- system.file("extdata", "genomeData", package="derfinder")
#' dirs <- makeBamList(datadir=datadir, samplepatt="*accepted_hits.bam$", bamterm=NULL)
#' ## Shorten the column names
#' names(dirs) <- gsub("_accepted_hits.bam", "", names(dirs))
#' 
#' ## Reading the data and filtering it is quite fast.
#' fullCov <- fullCoverage(dirs=dirs, chrnums="21", mc.cores=1)
#' 
#' ## Assign chr lengths using hg19 information
#' library("GenomicRanges")
#' data(hg19Ideogram, package = "biovizBase", envir = environment())
#' regions <- genomeRegions$regions
#' seqlengths(regions) <- seqlengths(hg19Ideogram)[names(seqlengths(regions))]
#'
#' ## Finally, get the region coverage
#' regionCov <- getRegionCoverage(fullCov=fullCov, regions=regions)
#' }


getRegionCoverage <- function(fullCov, regions, totalMapped = NULL,
	mc.cores=1, verbose=TRUE) {

	names(regions) <- seq_len(length(regions)) # add names
	
	if(sum(grepl("chr", names(fullCov))) == 0) {
		names(fullCov) = paste0("chr",names(fullCov))
	}
	
	## Warning when seqlengths are not specified
	if(any(is.na(seqlengths(regions)))) warning("'regions' does not have seqlengths assigned! In some cases, this can lead to erroneous results. getRegionCoverage() will proceed, but please check for other warnings or errors.")
	
	grl <- split(regions, as.character(seqnames(regions))) # split by chromosome
	counts <- mclapply(grl, function(g) { # now can be parallel
		if(verbose) cat(".")
		thechr <- as.character(unique(seqnames(g)))
		yy <- fullCov[[thechr]][ranges(g),] # better subset
		# depth-adjust, like for plotting
		if(!is.null(totalMapped)) {
			yy <- DataFrame(mapply(function(x,d) x/(d/1e6), yy, totalMapped))
		}
		ind <- rep(names(g), width(g)) # to split along
		ind <- factor(ind, levels = unique(ind)) # make factor in order
		# split(yy,ind) # "CompressedSplitDataFrameList", faster but less clear
						#   how to unlist below, so leave out
		split(as.data.frame(yy), ind) 
	}, mc.cores=mc.cores)
	covList <- do.call("c", counts) # collect list elements into one large list
	
	# put in original order
	names(covList) <- sapply(strsplit(names(covList), "\\."), "[", 2)
	theData <- covList[order(as.numeric(names(covList)))]	
	
	if(sum(sapply(theData, nrow)) != sum(width(regions))) {
		stop("The total width of the regions did not match with the dimensions of the extracted coverage data. Check that 'regions' has seqlengths specified correctly.")
	}
	
	return(theData)
}