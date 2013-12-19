#' Extract coverage information for a set of regions
#'
#' This function extracts the raw coverage information calculated by \link{fullCoverage} at each base for a set of regions found with \link{calculatePvalues}. It can further calculate the mean coverage per sample for each region.
#' 
#' @param fullCov A list where each element is the result from \link{loadCoverage} used with \code{cutoff=NULL}. The elements of the list should be named according to the chromosome number. Can be generated using \link{fullCoverage}.
#' @param regions The \code{$regions} output from \link{calculatePvalues}. It is important that the seqlengths information is provided.
#' @param calculateMeans If \code{TRUE} the mean coverage per sample for each region is calculated.
#' @param verbose If \code{TRUE} basic status updates will be printed along the way.
#'
#' @return A list with elements \code{coverageData} and \code{coverageMeans} (only if \code{calculateMeans=TRUE}). 
#' \describe{
#' \item{coverageData }{This is a list of data.frame where each data.frame has the coverage information (nrow = width of region, ncol = number of samples) for a given region. The names of the list correspond to the region indexes in \code{regions}.}
#' \item{coverageMeans }{This is a matrix (nrow = number of regions, ncol = number of samples) with the mean coverage per sample for all the regions.}
#' }
#'
#' @author Andrew Jaffe, Leonardo Collado-Torres
#' @seealso \link{fullCoverage}, \link{calculatePvalues}
#' @export
#' @importFrom GenomicRanges seqlevels seqnames seqlengths
#' @importMethodsFrom GenomicRanges names "names<-" length "[" coverage sort width c
#' @importMethodsFrom IRanges subset as.data.frame
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
#' regions <- regsWithP$regions
#' seqlengths(regions) <- seqlengths(hg19Ideogram)[names(seqlengths(regions))]
#'
#' ## Finally, get the region coverage
#' regionCov <- getRegionCoverage(fullCov=fullCov, regions=regions)
#' }


getRegionCoverage <- function(fullCov, regions, mc.cores=1, verbose=TRUE) {
	names(regions) <- seq_len(length(regions)) # add names
	fullCov <- fullCov[gsub("chr", "", seqlevels(regions))] ## added
	
	## Warning when seqlengths are not specified
	if(any(is.na(seqlengths(regions)))) warning("'regions' does not have seqlengths assigned! In some cases, this can lead to erroneous results. getRegionCoverage() will proceed, but please check for other warnings or errors.")
	
	grl = split(regions, seqnames(regions)) # split by chromosome
	counts = mclapply(grl, function(g) { # now can be parallel
		cat(".")
		thechr = as.character(unique(seqnames(g)))
		yy = y[[thechr]][ranges(g),] # better subset
		ind = rep(names(g), width(g)) # to split along
		ind = factor(ind, levels = unique(ind)) # make factor in order
		# split(yy,ind) # "CompressedSplitDataFrameList", faster but less clear
						#   how to unlist below, so leave out
		split(as.data.frame(yy),ind) 
	}, mc.cores=mc.cores)
	covList = do.call("c",counts) # collect list elements into one large list
	
	# put in original order
	names(covList) = sapply(strsplit(names(covList), "\\."), "[", 2)
	theData = covList[order(as.numeric(names(covList)))]	
	
	if(sum(sapply(theData, nrow)) != sum(width(regions))) {
		stop("The total width of the regions did not match with the dimensions of the extracted coverage data. Check that 'regions' has seqlengths specified correctly.")
	}
	
	return(theData)
}