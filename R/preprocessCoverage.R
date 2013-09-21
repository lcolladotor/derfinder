#' Transform and split the data
#'
#' This function takes the coverage data from \link{loadCoverage}, scales the data, does the log2 transformation, and splits it into appropriate chunks for using \link{calculateStats}.
#' 
#' @param coverageInfo A list containing a DataFrame --\code{$coverage}-- with the coverage data and a logical Rle --\code{$position}-- with the positions that passed the cutoff. This object is generated using \link{loadCoverage}.
#' @param groupInfo A factor specifying the group membership of each sample. If \code{NULL} no group mean coverages are calculated. If the factor has more than one level, the first one will be used to calculate the log2 fold change in \link{calculatePvalues}.
#' @param cutoff This argument is passed to \link{filterData}.
#' @param colsubset Optional vector of column indices of \code{coverageInfo$coverage} that denote samples you wish to include in analysis. 
#' @param scalefac A log transformation is used on the count tables, so zero counts present a problem.  What number should we add to the entire matrix?
#' @param chunksize How many rows of \code{coverageInfo$coverage} should be processed at a time?
#' @param verbose If \code{TRUE} basic status updates will be printed along the way.
#' @param mc.cores This argument is passed to \link[parallel]{mclapply} to run \link{fstats.apply}.
#'
#' @details If \code{chunksize} is \code{NULL}, then \code{mc.cores} is used to determine the \code{chunksize}. This is useful if you want to split the data so each core gets the same amount of data (up to rounding).
#'
#' @return A list with four components.
#' \describe{
#' \item{coverageSplit }{ is a list of DataFrames. Each DataFrame column represents a sample and the data is partioned according to \code{chunksize}. The coverage information is scaled and log2 transformed. Note that if \code{colsubset} is not \code{NULL} the number of columns will be less than those in \code{coverageInfo$coverage}. The total number of rows depends on the number of base pairs that passed the \code{cutoff} and the information stored is the coverage at that given base. Further note that \link{filterData} is re-applied if \code{colsubset} is not \code{NULL} and could thus lead to fewer rows compared to \code{coverageInfo$coverage}. }
#' \item{position }{  is a logical Rle with the positions of the chromosome that passed the cutoff.}
#' \item{meanCoverage }{ is a numeric Rle with the mean coverage at each filtered base.}
#' \item{groupMeans }{ is a list of Rle objects containing the mean coverage at each filtered base calculated by group. This list has length 0 if \code{groupInfo=NULL}.}
#' }
#'
#' @author Leonardo Collado-Torres
#' @seealso \link{filterData}, \link{loadCoverage}, \link{calculateStats}
#' @export
#' @importMethodsFrom IRanges ncol nrow sapply "[" "[[" "[[<-" c split
#' @importFrom IRanges Rle
#' @examples
#' ## Split the data and transform appropriately before using calculateStats()
#' dataReady <- preprocessCoverage(genomeData, cutoff=0, scalefac=32, chunksize=1e3, colsubset=NULL, verbose=TRUE)
#' names(dataReady)
#' dataReady

preprocessCoverage <- function(coverageInfo, groupInfo=NULL, cutoff = 5, scalefac = 32, chunksize = 5e6, colsubset = NULL, mc.cores=getOption("mc.cores", 1L), verbose=FALSE) {
	## Check that the input is from loadCoverage()
	stopifnot(length(intersect(names(coverageInfo), c("coverage", "position"))) == 2)
	stopifnot(is.factor(groupInfo) | is.null(groupInfo))
	
	coverage <- coverageInfo$coverage
	if(is.null(colsubset)) {
		stopifnot(length(groupInfo) == length(coverage) | is.null(groupInfo))
	} else {
		stopifnot(length(groupInfo) == length(colsubset) | is.null(groupInfo))
	}
		
	position <- coverageInfo$position
		
	## Subset the DataFrame to use only the columns of interest
	if(!is.null(colsubset)) {
		## Re-filter
		if(verbose) message(paste(Sys.time(), "preprocessCoverage: filtering the data"))
		coverageInfo <- filterData(data=coverageInfo$coverage[, colsubset], cutoff=cutoff, index=coverageInfo$position, verbose=verbose)
		coverage <- coverageInfo$coverage
		position <- coverageInfo$position
	}
	rm(coverageInfo)
	
	## Get the positions and shorter variables
	numrow <- nrow(coverage)
	
	## Automatic chunksize depending on the number of cores
	if(is.null(chunksize)) {
		chunksize <- ceiling(numrow / mc.cores)
		if(verbose) message(paste(Sys.time(), "preprocessCoverage: using chunksize", chunksize))
	}
	
	## Determine number of loops
	lastloop <- trunc(numrow / chunksize)
	
	## Fix the lastloop in case that the N is a factor of chunksize
	if(numrow %% chunksize == 0 & lastloop > 0)  {
		lastloop <- lastloop - 1
	}
	
	## Find the overall mean coverage
	means <- .rowMeansDF(coverage)
	
	## Find the by group mean coverage
	if(is.null(groupInfo)) {
		groupMeans <- list()
	} else {
		groupMeans <- vector("list", length(levels(groupInfo)))
		names(groupMeans) <- levels(groupInfo)
		for(group in levels(groupInfo)) {
			groupMeans[[group]] <- .rowMeansDF(coverage[groupInfo == group])
		}
	}
	
		
	## Log2 transform and scale
	numcol <- ncol(coverage)
	for(i in seq_len(numcol)) {
		coverage[[i]] <- log2(coverage[[i]] + scalefac)
	}
	
	## Split the data into appropriate chunks
	if(verbose) message(paste(Sys.time(), "preprocessCoverage: splitting the data"))
	if(lastloop == 0) {
		split.len <- numrow
	} else {
		split.len <- rep(chunksize, lastloop)
		split.len.sum <- numrow - sum(split.len)
		if(split.len.sum > 0) {
			split.len <- c(split.len, split.len.sum)
		}
	}
	split.idx <- Rle(seq_len(lastloop + 1), split.len)
	coverage.split <- as.list(split(coverage, split.idx))
	
	## Done =)
	result <- list("coverageSplit"=coverage.split, "position"=position, "meanCoverage"=means, "groupMeans"=groupMeans)
	return(result)	
	
}

## Auxiliary function that calculates the rowMeans of a DataFrame
.rowMeansDF <- function(coverage) {
	means <- Rle(0, nrow(coverage))
	for(i in seq_len(length(coverage))) {
		means <- means + coverage[[i]]
	}
	means <- means / length(coverage)
	return(means)
}
