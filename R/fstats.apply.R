#' Calculate F-statistics per base by extracting chunks from a DataFrame.
#'
#' Extract chunks from a DataFrame, apply the scaling factor, log2 transform and then get the F-statistics. This is a helper function for \link{calculateStats}.
#'
#' @param i The chunk number identifier.
#' @param data The DataFrame containing the coverage information. Normally stored in \code{coverageInfo$coverage} from \link{makeCoverage}.
#' @param comparison Whether you are comparing if there is \code{expression} (intercept vs no intercept) or \code{group differences} (model with group vs no group). 
#' @param chunksize How many rows of \code{data} should be processed at a time?
#' @param lastloop The last chunk number.
#' @param numrow Total number of rows in \code{data}.
#' @param scalefac A log transformation is used on the count tables, so zero counts present a problem.  What number should we add to the entire matrix before running the models?
#' @param mod The design matrix for the alternative model. Should be m by p where p is the number of covariates (normally also including the intercept).
#' @param mod0 The deisgn matrix for the null model. Should be m by p_0.
#'
#' @return A Rle with the F-statistics per base for the chunk in question.
#'
#' @author Leonardo Collado-Torres
#' @export
#' @import IRanges
#' @seealso \link{calculateStats}, \link{fstats}.
#' @examples
#' ## Create the model matrices
#' mod <- model.matrix(~ brainInfo$outcome)
#' mod0 <- model.matrix(~ 0 + rep(1, nrow(brainInfo)))
#' ## Run the function
#' fstats.output <- fstats.apply(1, brainData$coverage, 1000, 5, nrow(brainData$coverage), 32, mod, mod0)
#' fstats.output
#' 

fstats.apply <- function(i, data, comparison, chunksize, lastloop, numrow, scalefac, mod, mod0) {
	## Define which bases to subset
	if(i!=lastloop) {
		index <- (chunksize * i + 1):(chunksize * (i + 1))
	} else  {
		index <- (chunksize * i + 1):numrow		
	}
	##  Subset the DataFrame to the current chunk and transform to a regular matrix
	mymat <- log2(as.matrix(as.data.frame(data[index,])) + scalefac)
	if(comparison == "expression"){
		## Subtract log2(scalefac) since we want to test beta_0 = 0, not beta_0 = log2(scalefac + 0), which is what we see with 0 expression under the transformation.
		mymat <- mymat - log2(scalefac)
	}
	
	## I don't think that we need the Amean for the F-stats.
	# Amean <- rowMeans(mymat) 
	
	## Get the Fstats
	stats <- Rle(fstats(mymat, mod, mod0))
	return(stats)
}
