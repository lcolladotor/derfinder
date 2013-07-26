#' Calculate p-values and identify regions
#'
#' First, this function clusters the genomic positions and finds the regions of interest according to specified cutoffs. Then it permutes the samples and re-calculates the F-statistics. The F-statistics are segmented using the original clusters and cutoffs. The mean of the statistics from these segments are then used to calculate p-values for the original regions.
#' 
#' @param statsInfo A list with \code{$coverage}, \code{$position}, \code{$fstats}, \code{$mod}, and \code{mod0} components as generated using \link{calculateStats}.
#' @param nPermute The number of permutations. Note that for a full chromosome, a small amount (10) of permutations is sufficient.
#' @param seeds An integer vector of length \code{nPermute} specifying the seeds to be used for each permutation. If \code{NULL} no seeds are used.
#' @param chunksize How many rows of \code{statsInfo$coverage} should be processed at a time?
#' @param chr A single element character vector specifying the chromosome name. This argument is passed to \link{findRegions}.
#' @param maxGap This argument is passed to \link{clusterMakerRle}.
#' @param cutoff This argument is passed to \link{getSegmentsRle}.
#' @param mc.cores This argument is passed to \link[parallel]{mclapply}) to run \link{fstats.apply}.
#' @param verbose If \code{TRUE} basic status updates will be printed along the way.
#'
#' @return A GRanges with metadata columns given by \link{findRegions} and
#' \describe{
#' \item{pvalues }{ p-value of the region calculated via permutations of the samples.}
#' }
#'
#' @author Leonardo Collado-Torres
#' @seealso \link{findRegions}, \link{clusterMakerRle}, \link{getSegmentsRle}, \link{fstats.apply}
#' @export
#' @importMethodsFrom IRanges nrow ncol c mean
#' @importFrom IRanges Views
#' @importFrom parallel mclapply
#' @examples
#' ## Get the statistics
#' group <- brainInfo$outcome
#' adjustvars <- brainInfo[, c("sex", "age", "left.hemisph", "pmi", "brainpH")]
#' statsInfo <- calculateStats(brainData, group, adjustvars=adjustvars, mc.cores=1, verbose=TRUE)
#' ## Calculate the p-values and define the regions of interest.
#' regsWithP <- calculatePvalues(statsInfo, nPermute=10, seeds=NULL, chr="chr21", cutoff=c(2, 5), mc.cores=1)
#' regsWithP
#' hist(regsWithP$pvalues)
#'
#' \dontrun{
#' ## Annotate the results
#' library("bumphunter")
#' annotation <- annotateNearest(regsWithP, "hg19")
#' head(annotation)
#'
#' ## Compare speed between 1 and 4 cores (must have them!)
#' ## The chunksize is artifically reduced just to actually need to run mclapply
#' library("microbenchmark")
#' micro <- microbenchmark(
#' calculatePvalues(statsInfo, nPermute=10, seeds=NULL, chr="chr21", cutoff=c(2, 5), chunksize=1e3, mc.cores=1, verbose=FALSE),
#' calculatePvalues(statsInfo, nPermute=10, seeds=NULL, chr="chr21", cutoff=c(2, 5), chunksize=1e3, mc.cores=4, verbose=FALSE),
#' times=10)
#' levels(micro$expr) <- c("one", "four")
#' micro
#' ## Doesn't seem to help with this toy data.
#' }

calculatePvalues <- function(statsInfo, nPermute = 1L, seeds = as.integer(gsub("-", "", Sys.Date())) + seq_len(nPermute), chunksize = 5e6, chr, maxGap = 300L, cutoff = quantile(abs(statsInfo$fstats), 0.99), mc.cores=getOption("mc.cores", 2L), verbose=TRUE) {
	## Setup
	if(is.null(seeds)) {
		seeds <- rep(NA, nPermute)
	}
	stopifnot(nPermute == length(seeds))
	stopifnot(length(intersect(names(statsInfo), c("coverage", "position", "fstats", "mod", "mod0"))) == 5)
	
	## Identify the clusters
	if(verbose) message("calculatePvalues: identifying clusters")
	cluster <- clusterMakerRle(statsInfo$position, maxGap)
	
	## Find the regions
	regs <- findRegions(statsInfo, chr=chr, fstats=statsInfo$fstats, cluster=cluster, y=statsInfo$fstats, cutoff=cutoff, verbose=verbose) 
	
	## Determine total and loop sizes
	numrow <- nrow(statsInfo$coverage)
	lastloop <- trunc(numrow / chunksize)
	
	## Fix the lastloop in case that the N is a factor of chunksize
	if(numrow %% chunksize == 0 & lastloop > 0)  {
		lastloop <- lastloop - 1
	}
	
	## Pre-allocate memory
	nullstats <- vector("list", length(seeds) * 2)
	last <- 0
		
	for(i in seq_along(seeds)) {
		if(verbose) message(paste("calculatePvalues: calculating F-statistics for permutation", i))		
		
		if(!is.na(seeds[i])) {
			set.seed(seeds[i])
		}
		idx.permute <- sample(seq_len(ncol(statsInfo$coverage)))
	
		### !!! I don't think that the model matrices have to get permuted too. Otherwise the results should be the same (except for some tiny numerical differences).
		## New model matrices
		## mod.p <- statsInfo$mod[idx.permute, ]
		## mod0.p <- statsInfo$mod0[idx.permute, ]
		
		## Permuted data
		data.p <- statsInfo$coverage[, idx.permute]
		
		## Get the F-statistics
		fstats.output <- mclapply(0:lastloop, fstats.apply, data=data.p, chunksize=chunksize, lastloop=lastloop, numrow=numrow, mod=statsInfo$mod, mod0=statsInfo$mod0, mc.cores=mc.cores)
		fstats.output <- do.call(c, fstats.output)
			
		## Find the segments
		Indexes <- getSegmentsRle(x = fstats.output, f = cluster, cutoff = cutoff, verbose = verbose, zero=FALSE)
		
		## Calculate mean statistics
	    for (j in 1:2) {
			nullstats[[last + j]] <- mean(Views(fstats.output, Indexes[[j]]))
	    }		
		last <- last + 2
		
		## Finish loop
		rm(idx.permute, data.p, fstats.output, Indexes)
	}
	nullstats <- sort(abs(do.call(c, nullstats)))
	
	## Calculate pvalues
	if(verbose) message("calculatePvalues: calculating the p-values")
	pvals <- sapply(abs(regs$value), function(x) { sum(nullstats > x) })
	pvals <- (pvals + 1) / length(nullstats)
	
	## Finish up
	regs$pvalues <- pvals
	
	## Done =)
	return(regs)	
}
