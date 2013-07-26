#' Calculate F-statistics per base by extracting chunks from a DataFrame
#'
#' Extract chunks from a DataFrame and get the F-statistics. This is a helper function for \link{calculateStats} and \link{calculatePvalues}.
#'
#' @param data The DataFrame containing the coverage information. Normally stored in \code{coverageInfo$coverage} from \link{loadCoverage}.
#' @param mod This argument is passed to \link{fstats}.
#' @param mod0 This argument is passed to \link{fstats}.
#'
#' @return A Rle with the F-statistics per base for the chunk in question.
#'
#' @author Leonardo Collado-Torres
#' @export
#' @importFrom IRanges Rle
#' @importMethodsFrom IRanges as.data.frame as.matrix
#' @seealso \link{calculateStats}, \link{fstats}, \link{calculatePvalues}
#' @examples
#' ## Create the model matrices
#' mod <- model.matrix(~ brainInfo$outcome)
#' mod0 <- model.matrix(~ 0 + rep(1, nrow(brainInfo)))
#' ## Run the function
#' fstats.output <- fstats.apply(brainData$coverage, mod, mod0)
#' fstats.output
#' 

fstats.apply <- function(data, mod, mod0) {
	##  Subset the DataFrame to the current chunk and transform to a regular matrix
	mymat <- as.matrix(as.data.frame(data))
		
	## I don't think that we need the Amean for the F-stats.
	# Amean <- rowMeans(mymat) 
	
	## Get the Fstats
	stats <- Rle(fstats(mymat, mod, mod0))
	return(stats)
}
