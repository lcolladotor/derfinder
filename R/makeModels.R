#' Build model matrices for differential expression
#'
#' Builds the model matrices for testing for differential expression by comparing a model with a grouping factor versus one without it. It adjusts for the confounders specified and the median coverage of each sample. The resulting models can be used in \link{calculateStats}.
#' 
#' @param coverageInfo A list containing a DataFrame --\code{$coverage}-- with the coverage data and a logical Rle --\code{$position}-- with the positions that passed the cutoff. This object is generated using \link{loadCoverage}.
#' @param testvars A vector or matrix specifying the variables to test. For example, a factor with the group memberships when testing for differences across groups. It's length should match the number of columns used from \code{coverageInfo$coverage}.
#' @param adjustvars Optional matrix of adjustment variables (e.g. measured confounders, output from SVA, etc.) to use in fitting linear models to each nucleotide. These variables have to be specified by sample and the number of rows must match the number of columns used. It will also work if it is a vector of the correct length.
#' @param nonzero If \code{TRUE}, use the median of only the nonzero counts as the library size adjustment.
#' @param verbose If \code{TRUE} basic status updates will be printed along the way.
#' @param center If \code{TRUE} the column medians are centered by the mean of the column medians. This is helpful for interpretation purposes.
#' @param testIntercept If \code{TRUE} then \code{testvars} is ignored and mod0 will contain the column medians and any adjusting variables specified, but no intercept.
#' @param colsubset Which colums of \code{coverageInfo$coverage} to use.
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
#' testvars <- genomeInfo$pop
#' colsubset <- NULL
#' adjustvars <- data.frame(genomeInfo$gender)
#' nonzero <- TRUE
#' verbose <- TRUE
#' 
#' ## Run the function
#' models <- makeModels(coverageInfo, testvars, adjustvars, nonzero, verbose)
#' names(models)
#' models

makeModels <- function(coverageInfo, testvars, adjustvars = NULL, nonzero = TRUE, verbose=FALSE, center=TRUE, testIntercept=FALSE, colsubset=NULL) {
	## Check that the input is from loadCoverage()
	stopifnot(length(intersect(names(coverageInfo), c("coverage", "position"))) == 2)
	coverage <- coverageInfo$coverage
	
	## Drop unused levels in testvars if it is a factor
	if(is.factor(testvars)) {
		testvars <- droplevels(testvars)
	}
	
	
	if(!is.null(colsubset)) {
		coverage <- coverage[, colsubset]
	}
		
	## Check that the columns match
	numcol <- ncol(coverage)
	if(!testIntercept & numcol != NROW(testvars)) {
		stop("The length of 'testvars' and the number of columns in 'coverageInfo$coverage' do not match.")
	} 
	if (!is.null(adjustvars) & NROW(adjustvars) != numcol) {
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
	
	if(any(is.na(colmedians))) {
		col.na <- is.na(colmedians)
		warning(paste0("Sample column(s) ", paste(which(col.na), collapse=", "), " median coverage (nonzero=", nonzero, ") are NA. Setting them to 0 (before centering if center=TRUE). Check for possible issues with this sample!"))
		colmedians[col.na] <- 0
	}
	
	if(center) {
		colmedians <- colmedians - mean(colmedians, na.rm=TRUE)
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
	if(!testIntercept) {
		eval(parse(text=paste0("mod = model.matrix(~ testvars + colmedians", string1, ")")))
		eval(parse(text=paste0("mod0 = model.matrix(~ + colmedians", string1, ")")))	
	} else {
		eval(parse(text=paste0("mod = model.matrix(~ colmedians", string1, ")")))
		eval(parse(text=paste0("mod0 = model.matrix(~ 0 + colmedians", string1, ")")))	
	}	
	
	## Check that the matrices are full rank
	if(qr(mod)$rank == ncol(mod)) {
		r <- qr(mod)$rank
		warning(paste("Dropping from the alternative model matrix (mod) the column(s)", paste(colnames(mod)[(r+1):ncol(mod)], collapse=", "), "as the matrix is not full rank."))
		mod <- mod[, seq_len(r), drop=FALSE]
		stopifnot(ncol(mod) > 0)
	}
	if(qr(mod0)$rank == ncol(mod0)) {
		r <- qr(mod0)$rank
		warning(paste("Dropping from the null model matrix (mod0) the column(s)", paste(colnames(mod0)[(r+1):ncol(mod0)], collapse=", "), "as the matrix is not full rank."))
		mod0 <- mod0[, seq_len(r), drop=FALSE]
		stopifnot(ncol(mod0) > 0)
	}
		
	## Finish
	result <- list(mod=mod, mod0=mod0)
	
	## Done =)
	return(result)	
	
}
