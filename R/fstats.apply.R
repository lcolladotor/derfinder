#' Calculate F-statistics per base by extracting chunks from a DataFrame
#'
#' Extract chunks from a DataFrame and get the F-statistics on the rows of \code{data}, comparing the models \code{mod} (alternative) and \code{mod0} (null). This is a helper function for \link{calculateStats} and \link{calculatePvalues}.
#'
#' @param data The DataFrame containing the coverage information. Normally stored in \code{coveragePrep$coverageSplit} from \link{preprocessCoverage}. Could also be the full data from \link{loadCoverage}.
#' @param mod The design matrix for the alternative model. Should be m by p where p is the number of covariates (normally also including the intercept).
#' @param mod0 The design matrix for the null model. Should be m by p_0.
#'
#' @return A numeric Rle with the F-statistics per base for the chunk in question.
#'
#' @author Jeff Leek, Leonardo Collado-Torres
#' @export
#' @importMethodsFrom IRanges as.data.frame as.matrix
#' @useDynLib derfinder2
#' @seealso \link{calculateStats}, \link{calculatePvalues}
#' @examples
#' ## Create the model matrices
#' mod <- model.matrix(~ genomeInfo$pop)
#' mod0 <- model.matrix(~ 0 + rep(1, nrow(genomeInfo)))
#' ## Run the function
#' fstats.output <- fstats.apply(genomeData$coverage, mod, mod0)
#' fstats.output
#' 

fstats.apply <- function(data, mod, mod0) {
	##  Subset the DataFrame to the current chunk and transform to a regular matrix
	dat <- as.matrix(as.data.frame(data))
	rm(data)
	
	# A function for calculating F-statistics
	# on the rows of dat, comparing the models
	# mod (alternative) and mod0 (null). 	
	fstats <- Rle(drop(.Call("rcppFstats", t(dat), mod, mod0, package="derfinder2")))
	
	## Done
	return(fstats)
}
