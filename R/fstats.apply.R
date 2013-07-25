#' Calculate F-statistics per base by extracting chunks from a DataFrame
#'
#' Extract chunks from a DataFrame, apply the scaling factor, log2 transform and then get the F-statistics. This is a helper function for \link{calculateStats}.
#'
#' @param i The chunk number identifier.
#' @param data The DataFrame containing the coverage information. Normally stored in \code{coverageInfo$coverage} from \link{loadCoverage}.
#' @param chunksize How many rows of \code{data} should be processed at a time?
#' @param lastloop The last chunk number.
#' @param numrow Total number of rows in \code{data}.
#' @param mod This argument is passed to \link{fstats}.
#' @param mod0 This argument is passed to \link{fstats}.
#'
#' @return A Rle with the F-statistics per base for the chunk in question.
#'
#' @author Leonardo Collado-Torres
#' @export
#' @importFrom IRanges Rle
#' @importMethodsFrom IRanges as.data.frame as.matrix "["
#' @seealso \link{calculateStats}, \link{fstats}.
#' @examples
#' ## Create the model matrices
#' mod <- model.matrix(~ brainInfo$outcome)
#' mod0 <- model.matrix(~ 0 + rep(1, nrow(brainInfo)))
#' ## Run the function
#' fstats.output <- fstats.apply(1, brainData$coverage, 1000, 5, nrow(brainData$coverage), mod, mod0)
#' fstats.output
#' 

fstats.apply <- function(i, data, chunksize, lastloop, numrow, mod, mod0) {
	## Define which bases to subset
	if(i!=lastloop) {
		index <- (chunksize * i + 1):(chunksize * (i + 1))
	} else  {
		index <- (chunksize * i + 1):numrow		
	}
	##  Subset the DataFrame to the current chunk and transform to a regular matrix
	mymat <- as.matrix(as.data.frame(data[index,]))
		
	## I don't think that we need the Amean for the F-stats.
	# Amean <- rowMeans(mymat) 
	
	## Get the Fstats
	stats <- Rle(fstats(mymat, mod, mod0))
	return(stats)
}
