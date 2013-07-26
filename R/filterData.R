#' Filter the positions of interest
#'
#' For a group of samples this function reads the coverage information for a specific chromosome directly from the BAM files. It then merges them into a DataFrame and removes the bases that do not pass the cutoff. This is a helper function for \link{loadCoverage} and \link{calculateStats}.
#' 
#' @param data Either a list of Rle objects or a DataFrame with the coverage information.
#' @param cutoff Per base pair, at least one sample has to have coverage greater than \code{cutoff} to be included in the result.
#' @param index A logical Rle with the positions of the chromosome that passed the cutoff. If \code{NULL} it is assumed that this is the first time using \link{filterData} and thus no previous index exists.
#' @param colnames Specifies the column names to be used for the results DataFrame. If \code{NULL}, no names are assigned.
#' @param verbose If \code{TRUE} it will report how many rows are remaining out of the original ones.
#'
#' @return A list with two components. 
#' \describe{
#' \item{coverage }{ is a DataFrame object where each column represents a sample. The number of rows depends on the number of base pairs that passed the cutoff and the information stored is the coverage at that given base.}
#' \item{position }{  is a logical Rle with the positions of the chromosome that passed the cutoff.}
#' }
#'
#' @author Leonardo Collado-Torres
#' @export
#' @importFrom IRanges DataFrame
#' @importMethodsFrom IRanges "[" "[<-" colnames "colnames<-" lapply
#' @seealso \link{loadCoverage}
#' @examples
#' library("IRanges")
#' x <- Rle(round(runif(1e4, max=10)))
#' y <- Rle(round(runif(1e4, max=10)))
#' z <- Rle(round(runif(1e4, max=10)))
#' DF <- DataFrame(x, y, z)
#' filt1 <- filterData(DF, 5)
#' filt1
#' filt2 <- filterData(filt1$coverage[, 1:2], 5, index=filt1$position)
#' filt2
#' ## The number of TRUE values in 'position' is the same as the number of rows as in 'coverage'.
#' identical(sum(filt2$pos), nrow(filt2$cov))

filterData <- function(data, cutoff, index=NULL, colnames=NULL, verbose=TRUE) {
	## Construct the filtering index
	for(i in 1:length(data)) {
		if(i == 1) {
			newindex <- data[[i]] > cutoff
		} else {
			newindex <- newindex | data[[i]] > cutoff
		}
	}
	
	## Keep only bases that pass the cutoff
	if(is(data, "DataFrame")) {
		DF <- data[newindex, ]
	} else {
		## Subset the data
		newdata <- lapply(data, function(x) { x[newindex] })
		
		## Group into DataFrame, though the list might be more useful for downstream stuff
		DF <- DataFrame(newdata)	
	}	
	
	## Assign column names
	if(!is.null(colnames)) {
		colnames(DF) <- colnames
	}	
	
	## Build the final index	
	if(!is.null(index)) {
		finalidx <- index
		finalidx[index] <- newindex
	} else {
		finalidx <- newindex
	}
	
	## Info for the user
	if(verbose) {
		message(paste(date(), "filterData: originally there were", length(data[[1]]), "rows, now there are", nrow(DF), "rows. Meaning that", 100 - round(nrow(DF) / length(data[[1]]) * 100, 2), "percent was filtered."))
	}
	
	## Make the final resulting object.
	res <- list("coverage"=DF, "position"=finalidx)
	return(res)
	
}
