#' Calculate adjustments for library size
#'
#' For a given data set  calculate the per-sample Q quantile of the base coverage. The data set can loaded to R using (see \link{fullCoverage} and optionally filtered using \link{filterData}. One the per-sample Q quantile has been calculated, you could be interested in centering them accross samples for interpretability concerns. This information is then used in \link{makeModels} for construcing the null and alternative models.
#' 
#' @param fullCov A list where each element is the result from \link[derfinder]{loadCoverage} used with \code{cutoff=NULL}. The elements of the list should be named according to the chromosome number. Can be generated using \link{fullCoverage}.
#' @param prob A number between 0 and 1 representing the quantile of interest. For example, 0.5 is the median.
#' @param nonzero If \code{TRUE} only the nonzero counts are used to calculate the library size adjustment.
#' @param center If \code{TRUE} the sample quantiles are centered by the mean of the sample quantiles across sampls. This can be helpful for interpretation purposes.
#' @param colsubset Which colums of \code{coverageInfo$coverage} to use.
#' @param verbose If \code{TRUE} basic status updates will be printed along the way.
#'
#' @return 
#' A vector with the library size depth adjustments per sample to be used in \link{makeModels}.
#'
#' @author Leonardo Collado-Torres
#' @seealso \link{fullCoverage}, \link{makeModels}
#' @export
#' @importMethodsFrom IRanges sapply "[" rbind quantile
#' @examples
#' ## Choose the adjusting variables and define all the parameters for makeModels()
#' coverageInfo <- genomeData
#' 
#' ## Calculate library size adjustments
#' sampleDepths <- sampleDepth(list(coverageInfo$coverage), prob=0.5, nonzero=TRUE, center=TRUE, verbose=TRUE)
#' sampleDepths

sampleDepth <- function(fullCov, prob = 0.9, nonzero = TRUE, center=TRUE, colsubset=NULL, verbose=FALSE) {
	stopifnot(prob >= 0 & prob <= 1 & length(prob) == 1)
	## Remove un-used columns
	if(!is.null(colsubset)) {
		fullCov <- lapply(fullCov, function(x) {
			x[, colsubset]
		})
	}
	## Append the information from all chrs
	coverage <- do.call(rbind, fullCov)
	
	## Get the medians of the columns
	if(nonzero) {
		sampleDepths <- sapply(coverage, function(y) { 
			## Catch cases where the is no data points greater than 0
			tmp <- try(quantile(y[y > 0], prob), silent=TRUE)
			res <- ifelse(inherits(tmp, "try-error"), 0, tmp)
			return(res)
		})
	} else {
		sampleDepths <- sapply(coverage, function(y) {
			quantile(y, prob)
		})
	}
	
	if(any(is.na(sampleDepths))) {
		col.na <- is.na(sampleDepths)
		warning(paste0("Sample column(s) ", paste(which(col.na), collapse=", "), " have sample depth coverage (nonzero=", nonzero, ") as NA. Setting them to 0 (before centering if center=TRUE). Check for possible issues with this sample!"))
		sampleDepths[col.na] <- 0
	}
	
	if(verbose) message(paste0(Sys.time(), " sampleDepth: sample depth by sample (before centering): ", paste(sampleDepths, collapse=", "), "."))
	
	if(center) {
		sampleDepths <- sampleDepths - mean(sampleDepths)
		if(verbose) message(paste0(Sys.time(), " sampleDepth: sample depth by sample (after centering): ", paste(sampleDepths, collapse=", "), "."))
	}
	
	
		
	## Done =)
	return(sampleDepths)	
}
