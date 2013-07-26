#' Calculate F-statistics at base pair resolution from a loaded BAM files
#'
#' After defining the models of interest and adjusting for confounders, this function extracts the data from loaded BAM files (\link{loadCoverage}) and calculates the F-statistics.
#' 
#' @param coverageInfo A list containing a DataFrame --\code{$coverage}-- with the coverage data and a logical Rle --\code{$position}-- with the positions that passed the cutoff. This object is generated using \link{loadCoverage}.
#' @param group A factor vector specifying the sample groups. It's length should match the number of columns used from \code{coverageInfo$coverage}.
#' @param comparison Whether you are comparing if there is \code{expression} (intercept vs no intercept) or \code{group differences} (model with group vs no group). 
#' @param colsubset Optional vector of column indices of \code{coverageInfo$coverage} that denote samples you wish to include in analysis. 
#' @param adjustvars Optional matrix of adjustment variables (e.g. measured confounders, output from SVA, etc.) to use in fitting linear models to each nucleotide. These variables have to be specified by sample and the number of rows must match the number of columns used.
#' @param cutoff Per base pair, at least one sample has to have coverage greater than \code{cutoff} to be included in the result.
#' @param scalefac A log transformation is used on the count tables, so zero counts present a problem.  What number should we add to the entire matrix before running the models?
#' @param nonzero If \code{TRUE}, use the median of only the nonzero counts as the library size adjustment.
#' @param chunksize How many rows of \code{coverageInfo$coverage} should be processed at a time?
#' @param mc.cores This argument is passed to \link[parallel]{mclapply} to run \link{fstats.apply}.
#' @param verbose If \code{TRUE} basic status updates will be printed along the way.
#'
#' @return A list with five components.
#' \describe{
#' \item{coverageSplit }{ is a SplitDataFrameList object where each column represents a sample and the data is partioned according to \code{chunksize}. The coverage information is scaled and log2 transformed. Note that if \code{colsubset} is not \code{NULL} the number of columns will be less than those in \code{coverageInfo$coverage}. The total number of rows depends on the number of base pairs that passed the \code{cutoff} and the information stored is the coverage at that given base. Further note that \link{filterData} is re-applied if \code{colsubset} is not \code{NULL} and could thus lead to fewer rows compared to \code{coverageInfo$coverage}. }
#' \item{position }{  is a logical Rle with the positions of the chromosome that passed the cutoff.}
#' \item{fstats }{ is a numeric Rle with the F-statistics per base pair that passed the cutoff.}
#' \item{mod }{ The alternative model matrix.}
#' \item{mod0 }{ The null model matrix.}
#' }
#' @details Partially based on \link[derfinder]{getLimmaInput.DF}.
#' @references Frazee et al. Biostatistics in review.
#'
#' @author Leonardo Collado-Torres
#' @export
#' @importFrom parallel mclapply
#' @importMethodsFrom IRanges ncol nrow sapply median "[[" "[[<-" c split unlist
#' @importFrom IRanges RleList
#' @examples
#' ## Choose the adjusting variables and define all the parameters for calculateStats()
#' coverageInfo <- brainData
#' group <- brainInfo$outcome
#' colsubset <- NULL
#' adjustvars <- brainInfo[, c("sex", "age", "left.hemisph", "pmi", "brainpH")]
#' cutoff <- 5
#' scalefac <- 32
#' nonzero <- TRUE
#' chunksize <- 1e+03
#' mc.cores <- 1	
#' verbose <- TRUE
#' 
#' ## Run the function
#' stats <- calculateStats(coverageInfo, group, comparison="group differences", colsubset, adjustvars, cutoff, scalefac, nonzero, chunksize, mc.cores, verbose)
#' names(stats)
#' stats

calculateStats <- function(coverageInfo, group, comparison = "group differences", colsubset = NULL, adjustvars = NULL, cutoff = 5, scalefac = 32, nonzero = FALSE, chunksize = 5e6, mc.cores=getOption("mc.cores", 2L), verbose=TRUE) {
	stopifnot(length(intersect(names(coverageInfo), c("coverage", "position"))) == 2)
	
	if(!comparison %in% c("expression", "group differences")) {
		stop("Invalid value for 'comparison'.")
	}
	
	## Subset the DataFrame to use only the columns of interest
	if(!is.null(colsubset)) {
		## Re-filter
		coverageInfo <- filterData(data=coverageInfo$coverage[, colsubset], cutoff=cutoff, index=coverageInfo$position, verbose=verbose)
	}
	
	## Get the positions and shorter variables
	data <- coverageInfo$coverage
	## I don't think that we need this anymore.
	#pos <- which(coverageInfo$position)
	
	## Check that the columns match
	numcol <- ncol(data)
	if(numcol != length(group)) {
		stop("The length of 'group' and the number of columns in 'coverageInfo$coverage[, colsubset]' do not match.")
	} else if (!is.null(adjustvars) & NROW(adjustvars) != numcol) {
		stop("The dimensions of 'adjustvars' should match with the number of columns in 'coverageInfo$coverage[, colsubset]'.")
	}
	
	## Determine total and loop sizes
	numrow <- nrow(data)
	lastloop <- trunc(numrow / chunksize)
	
	## Fix the lastloop in case that the N is a factor of chunksize
	if(numrow %% chunksize == 0 & lastloop > 0)  {
		lastloop <- lastloop - 1
	}
		
	## Get the medians of the columns
	if(nonzero) {
		colmeds <- sapply(data, function(y) { median(y[y > 0]) })
	} else {
		colmeds <- sapply(data, median)
	}
	## Info for the user
	if(verbose) message(paste0(date(), " calculateStats: these are the column medians used: ", paste(colmeds, collapse=", "), "."))

	## Create model matrices
	if(!is.null(adjustvars)){
		## Build the adjusted variables
		string1 <- ""
		for(i in 1:dim(adjustvars)[2]){
			eval(parse(text=paste0("av", i, " <- adjustvars[,", i, "]")))
			string1 <- paste(string1, paste0("av", i), sep="+")
		}
		if(comparison == "expression") {
			warning("Check that it makes sense! TODO")
			eval(parse(text=paste0("mod = model.matrix(~ as.factor(group) + colmeds", string1, ")")))
			eval(parse(text=paste0("mod0 = model.matrix(~ 0 + as.factor(group) + colmeds", string1, ")")))
		} else {
			eval(parse(text=paste0("mod = model.matrix(~ as.factor(group) + colmeds", string1, ")")))
			eval(parse(text=paste0("mod0 = model.matrix(~ + colmeds", string1, ")")))
		}		
	} else {
		if(comparison == "expression"){ 
			stop("Fstats for 1 column matrix vs 0? TODO")
			
			## For expression-test, add a column of 1s to the model matrix
			mod <- model.matrix(~ 0 + rep(1, numcol)) #intercept-only model, in model.matrix form.
			mod0 <- model.matrix(~ 0 + rep(0, numcol))
		} else	{
			## For differential expresion
			mod <- model.matrix(~ as.factor(group) + colmeds)
			mod0 <- model.matrix(~ colmeds)
		}
	}	
	
	## Log2 transform and scale
	for(i in seq_len(numcol)) {
		data[[i]] <- log2(data[[i]] + scalefac)
	}
	
	## Split the data into appropriate chunks
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
	
	## Fit a model to each row (chunk) of database:
	if(verbose) message(paste(date(), "calculateStats: calculating the F-statistics"))
	fstats.output <- mclapply(data.split, fstats.apply, mod=mod, mod0=mod0, mc.cores=mc.cores)
	## Using mclapply is as fast as using lapply if mc.cores=1, so there is no damage in setting the default mc.cores=1. Specially since parallel is included in R 3.0.x
	## More at http://stackoverflow.com/questions/16825072/deprecation-of-multicore-mclapply-in-r-3-0
	fstats.output <- unlist(RleList(fstats.output), use.names=FALSE)
	
	## Done =)
	result <- list("coverageSplit"=data.split, "position"=coverageInfo$position, "fstats"=fstats.output, "mod"=mod, "mod0"=mod0)
	return(result)	
	
}
