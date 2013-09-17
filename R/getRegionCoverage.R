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
#' ## Construct the models
#' group <- genomeInfo$pop
#' adjustvars <- data.frame(genomeInfo$gender)
#' models <- makeModels(coverageInfo=genomeData, testvars=group, adjustvars=adjustvars, nonzero=TRUE)
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
#' n <- dim(prep$coverageSplit[[1]])[2]
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
#' datadir <- system.file("extdata", "genomeData", package="derfinder2")
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


getRegionCoverage <- function(fullCov, regions, calculateMeans=TRUE, verbose=TRUE) {
	names(regions) <- seq_len(length(regions)) # add names
	fullCov <- fullCov[gsub("chr", "", seqlevels(regions))] ## added
	
	## Warning when seqlengths are not specified
	if(all(is.na(seqlengths(regions)))) warning("'regions' does not have seqlengths assigned! In some cases, this can lead to erroneous results. getRegionCoverage() will proceed, but please check for other warnings or errors.")
	
	## use logical rle to subset large coverage matrix
	cc <- coverage(regions) ## The output of coverage() differs on whether seqlenths are provided or not. If absent, it assumes that the last base is the end of the chromosome.
	for(i in seq(along=cc)) {
		cc[[i]]@values <- ifelse(cc[[i]]@values > 0, TRUE, FALSE)
	}

	fullCovList <- list()
	for(i in seq(along=fullCov)) {
		if(verbose) message(paste(Sys.time(), "getRegionCoverage: processing chromosome", names(fullCov)[i]))
		z <- as.data.frame(subset(fullCov[[i]], cc[[i]]))
		g <- sort(regions[seqnames(regions) == seqlevels(regions)[i]])
		ind <- rep(names(g), width(g))
		tmpList <- split(z, ind)
		fullCovList[[i]] <- tmpList
	}
	tmp <- do.call(c, fullCovList)
	theData <- tmp[order(as.numeric(names(tmp)), decreasing=FALSE)]
	out <- list(coverageData = theData)
	
	if(sum(unlist(lapply(out$coverageData, nrow))) != sum(width(regions))) {
		stop("The total width of the regions did not match with the dimensions of the extracted coverage data. Check that 'regions' has seglengths specified correctly.")
	}
	
	if(calculateMeans) {
		theMeans <- t(sapply(theData, colMeans))
		out$coverageMeans <- theMeans
	}

	return(out)
}