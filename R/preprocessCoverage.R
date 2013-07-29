#' Transform and split the data
#'
#' This function takes the coverage data from \link{loadCoverage}, scales the data, does the log2 transformation, and splits it into appropriate chunks for using \link{calculateStats}.
#' 
#' @param coverageInfo A list containing a DataFrame --\code{$coverage}-- with the coverage data and a logical Rle --\code{$position}-- with the positions that passed the cutoff. This object is generated using \link{loadCoverage}.
#' @param cutoff Per base pair, at least one sample has to have coverage greater than \code{cutoff} to be included in the result.
#' @param colsubset Optional vector of column indices of \code{coverageInfo$coverage} that denote samples you wish to include in analysis. 
#' @param scalefac A log transformation is used on the count tables, so zero counts present a problem.  What number should we add to the entire matrix before running the models?
#' @param chunksize How many rows of \code{coverageInfo$coverage} should be processed at a time?
#' @param verbose If \code{TRUE} basic status updates will be printed along the way.
#'
#' @return A list with two components.
#' \describe{
#' \item{coverageSplit }{ is a SplitDataFrameList object where each column represents a sample and the data is partioned according to \code{chunksize}. The coverage information is scaled and log2 transformed. Note that if \code{colsubset} is not \code{NULL} the number of columns will be less than those in \code{coverageInfo$coverage}. The total number of rows depends on the number of base pairs that passed the \code{cutoff} and the information stored is the coverage at that given base. Further note that \link{filterData} is re-applied if \code{colsubset} is not \code{NULL} and could thus lead to fewer rows compared to \code{coverageInfo$coverage}. }
#' \item{position }{  is a logical Rle with the positions of the chromosome that passed the cutoff.}
#' }
#'
#' @author Leonardo Collado-Torres
#' @seealso \link{loadCoverage}, \link{calculateStats}
#' @export
#' @importMethodsFrom IRanges ncol nrow sapply "[" "[[" "[[<-" c split
#' @importFrom IRanges Rle
#' @examples
#' ## Split the data and transform appropriately before using calculateStats()
#' dataReady <- preprocessCoverage(genomeData, cutoff=0, scalefac=32, chunksize=1e3, colsubset=NULL, verbose=TRUE)
#' names(dataReady)
#' dataReady

preprocessCoverage <- function(coverageInfo, cutoff = 5, scalefac = 32, chunksize = 5e6, colsubset = NULL, verbose=FALSE) {
	## Check that the input is from loadCoverage()
	stopifnot(length(intersect(names(coverageInfo), c("coverage", "position"))) == 2)
		
	## Subset the DataFrame to use only the columns of interest
	if(!is.null(colsubset)) {
		## Re-filter
		if(verbose) message(paste(date(), "preprocessCoverage: filtering the data"))
		coverageInfo <- filterData(data=coverageInfo$coverage[, colsubset], cutoff=cutoff, index=coverageInfo$position, verbose=verbose)
	}
	
	## Get the positions and shorter variables
	data <- coverageInfo$coverage
	
	## Determine total and loop sizes
	numrow <- nrow(data)
	lastloop <- trunc(numrow / chunksize)
	
	## Fix the lastloop in case that the N is a factor of chunksize
	if(numrow %% chunksize == 0 & lastloop > 0)  {
		lastloop <- lastloop - 1
	}
		
	## Log2 transform and scale
	numcol <- ncol(coverageInfo$coverage)
	for(i in seq_len(numcol)) {
		data[[i]] <- log2(data[[i]] + scalefac)
	}
	
	## Split the data into appropriate chunks
	if(verbose) message(paste(date(), "preprocessCoverage: splitting the data"))
	if(lastloop == 0) {
		split.len <- numrow
	} else {
		split.len <- rep(chunksize, lastloop)
		split.len.sum <- numrow - sum(split.len)
		if(split.len.sum > 0) {
			split.len <- c(split.len, split.len.sum)
		}
	}
	split.idx <- Rle(0:lastloop, split.len)
	data.split <- split(data, split.idx)
	
	## Done =)
	result <- list("coverageSplit"=data.split, "position"=coverageInfo$position)
	return(result)	
	
}
