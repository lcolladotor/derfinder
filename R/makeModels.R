#' Build model matrices for differential expression
#'
#' Builds the model matrices for testing for differential expression by comparing a model with a grouping factor versus one without it. It adjusts for the confounders specified and the median coverage of each sample. The resulting models can be used in \link{calculateStats}.
#' 
#' @param coverageInfo A list containing a DataFrame --\code{$coverage}-- with the coverage data and a logical Rle --\code{$position}-- with the positions that passed the cutoff. This object is generated using \link{loadCoverage}.
#' @param group A factor vector specifying the sample groups. It's length should match the number of columns used from \code{coverageInfo$coverage}.
#' @param adjustvars Optional matrix of adjustment variables (e.g. measured confounders, output from SVA, etc.) to use in fitting linear models to each nucleotide. These variables have to be specified by sample and the number of rows must match the number of columns used. It will also work if it is a vector of the correct length.
#' @param nonzero If \code{TRUE}, use the median of only the nonzero counts as the library size adjustment.
#' @param verbose If \code{TRUE} basic status updates will be printed along the way.
#'
#' @return A list with two components.
#' \describe{
#' \item{mod }{ The alternative model matrix.}
#' \item{mod0 }{ The null model matrix.}
#' }
#'
#' @author Leonardo Collado-Torres
#' @seealso \link{calculateStats}
#' @export
#' @importMethodsFrom IRanges ncol sapply median "["
#' @examples
#' ## Choose the adjusting variables and define all the parameters for makeModels()
#' coverageInfo <- genomeData
#' group <- genomeInfo$pop
#' colsubset <- NULL
#' adjustvars <- data.frame(genomeInfo$gender)
#' nonzero <- TRUE
#' verbose <- TRUE
#' 
#' ## Run the function
#' models <- makeModels(coverageInfo, group, adjustvars, nonzero, verbose)
#' names(models)
#' models

makeModels <- function(coverageInfo, group, adjustvars = NULL, nonzero = FALSE, verbose=FALSE) {
	## Check that the input is from loadCoverage()
	stopifnot(length(intersect(names(coverageInfo), c("coverage", "position"))) == 2)
	coverage <- coverageInfo$coverage
		
	## Check that the columns match
	numcol <- ncol(coverage)
	if(numcol != NROW(group)) {
		stop("The length of 'group' and the number of columns in 'coverageInfo$coverage' do not match.")
	} else if (!is.null(adjustvars) & NROW(adjustvars) != numcol) {
		stop("The dimensions of 'adjustvars' should match with the number of columns in 'coverageInfo$coverage'.")
	}
			
	## Get the medians of the columns
	if(nonzero) {
		colmedians <- sapply(coverage, function(y) { 
			## Catch cases where the is no data points greater than 0
			tmp <- try(median(y[y > 0]), silent=TRUE)
			res <- ifelse(inherits(tmp, "try-error"), 0, tmp)
			return(res)
		})
	} else {
		colmedians <- sapply(coverage, median)
	}
	rm(coverage)
	## Info for the user
	if(verbose) message(paste0(Sys.time(), " makeModels: these are the column medians used: ", paste(colmedians, collapse=", "), "."))
		
	## To avoid a warning in R CMD check
	mod <- mod0 <- NULL
	
	## Build the adjusted variables	if needed
	string1 <- ""
	if(!is.null(adjustvars)){
		if(NCOL(adjustvars) == 1) {
			## Only if a vector was supplied
			adjustvars <- as.data.frame(adjustvars)
		}		
		for(i in seq_len(NCOL(adjustvars))) {
			eval(parse(text=paste0("adjustVar", i, " <- adjustvars[,", i, "]")))
			string1 <- paste(string1, paste0("adjustVar", i), sep="+")
		}		
	} 
	eval(parse(text=paste0("mod = model.matrix(~ group + colmedians", string1, ")")))
	eval(parse(text=paste0("mod0 = model.matrix(~ + colmedians", string1, ")")))	
		
	## Finish
	result <- list(mod=mod, mod0=mod0)
	
	## Done =)
	return(result)	
	
}
