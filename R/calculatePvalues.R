#' Calculate p-values and identify regions
#'
#' First, this function clusters the genomic positions and finds the regions of interest according to specified cutoffs. Then it permutes the samples and re-calculates the F-statistics. The F-statistics are segmented using the original clusters and cutoffs. The area of the statistics from these segments are then used to calculate p-values for the original regions.
#' 
#' @param coveragePrep A list with \code{$coverageSplit} and \code{$position} normally generated using \link{preprocessCoverage}.
#' @param models A list with \code{$mod} and \code{$mod0} normally generated using \link{makeModels}.
#' @param fstats A numerical Rle with the F-statistics normally generated using \link{calculateStats}.
#' @param nPermute The number of permutations. Note that for a full chromosome, a small amount (10) of permutations is sufficient.
#' @param seeds An integer vector of length \code{nPermute} specifying the seeds to be used for each permutation. If \code{NULL} no seeds are used.
#' @param chr A single element character vector specifying the chromosome name. This argument is passed to \link{findRegions}.
#' @param maxGap This argument is passed to \link{clusterMakerRle}.
#' @param cutoff This argument is passed to \link{getSegmentsRle}.
#' @param mc.cores This argument is passed to \link[parallel]{mclapply} to run \link{fstats.apply}.
#' @param verbose If \code{TRUE} basic status updates will be printed along the way.
#' @param significantCut A vector of length two specifiying the cutoffs used to determine significance. The first element is used to determine significance for the p-values and the second element is used for the q-values.
#'
#' @return A list with three components:
#' \describe{
#' \item{regions }{ is a GRanges with metadata columns given by \link{findRegions} with the additional metadata column \code{pvalues}: p-value of the region calculated via permutations of the samples; \code{padj}: the qvalues calculated using \link[qvalue]{qvalue}; \code{significant}: whether the p-value is less than 0.05 (by default); \code{significantPadj}: whether the q-value is less than 0.10 (by default). It also includes the mean coverage of the region (mean from the mean coverage at each base calculated in \link{preprocessCoverage}).}
#' \item{nullstats}{ is a numeric Rle with the mean of the null statistics by segment.}
#' \item{nullwidths}{ is a numeric Rle with the length of each of the segments in the null distribution. The area can be obtained by multiplying the absolute \code{nullstats} by the corresponding lengths.}
#' }
#'
#' @author Leonardo Collado-Torres
#' @seealso \link{findRegions}, \link{clusterMakerRle}, \link{getSegmentsRle}, \link{fstats.apply}, \link[qvalue]{qvalue}
#' @export
#' @importMethodsFrom IRanges quantile nrow ncol c mean lapply unlist as.numeric "$" "$<-"
#' @importFrom IRanges Views RleList Rle IRanges Views
#' @importFrom parallel mclapply
#' @importFrom qvalue qvalue
#'
#' @examples
#' ## Construct the models
#' group <- genomeInfo$pop
#' adjustvars <- data.frame(genomeInfo$gender)
#' models <- makeModels(coverageInfo=genomeData, testvars=group, adjustvars=adjustvars, nonzero=TRUE)
#'
#' ## Preprocess the data
#' ## Automatic chunksize used to then compare 1 vs 4 cores in the 'do not run' section
#' prep <- preprocessCoverage(genomeData, cutoff=0, scalefac=32, chunksize=NULL, colsubset=NULL, mc.cores=4)
#' 
#' ## Get the F statistics
#' fstats <- calculateStats(prep, models, mc.cores=1, verbose=TRUE)
#'
#' ## Determine a cutoff from the F-distribution.
#' ## This step is very important and you should consider using quantiles from the observed F statistics
#' n <- dim(prep$coverageSplit[[1]])[2]
#' df1 <- dim(models$mod)[2]
#' df0 <- dim(models$mod0)[2]
#' cutoff <- qf(0.95, df1-df0, n-df1)
#' 
#' ## Low cutoff used for illustrative purposes
#' cutoff <- 1
#'
#' ## Calculate the p-values and define the regions of interest.
#' regsWithP <- calculatePvalues(prep, models, fstats, nPermute=10, seeds=NULL, chr="chr21", cutoff=cutoff, mc.cores=1)
#' regsWithP
#'
#' ## Calculate areas of the null segments
#' abs(regsWithP$nullstats) * regsWithP$nullwidths
#'
#' ## Histogram of the original F statistics
#' f.ori <- as.numeric(fstats)
#' hist(f.ori, main="Distribution original F-stats", freq=FALSE)
#'
#' ## Histogram of the null F statistics
#' f.null <- as.numeric(regsWithP$nullstats)
#' hist(f.null, main="Distribution null F-stats by region", freq=FALSE)
#'
#' ## Histogram of the original p-values
#' hist(pf(f.ori, df1-df0, n-df1), main="Distribution original p-values", freq=FALSE)
#' 
#' ## Histogram of the null p-values
#' hist(pf(f.null, df1-df0, n-df1), main="Distribution null p-values by region", freq=FALSE)
#'
#' ## By region
#'
#' ## Histogram of the original F statistics by region
#' hist(regsWithP$regions$value, main="Distribution original F-stats by region", freq=FALSE)
#'
#' ## Histogram of the original p-values by region
#' hist(pf(regsWithP$regions$value, df1-df0, n-df1), main="Distribution original p-values by region", freq=FALSE)
#'
#' ## Histogram of the p-values by region
#' hist(regsWithP$regions$pvalues, main="Distribution permutted p-values by region", freq=FALSE)
#'
#' \dontrun{
#' ## Annotate the results
#' library("bumphunter")
#' annotation <- annotateNearest(regsWithP$regions, "hg19")
#' head(annotation)
#'
#' ## Compare speed between 1 and 4 cores (must have them!)
#' library("microbenchmark")
#' micro <- microbenchmark(
#' calculatePvalues(prep, models, fstats, nPermute=10, seeds=NULL, chr="chr21", cutoff=c(2, 5), mc.cores=1, verbose=FALSE),
#' calculatePvalues(prep, models, fstats, nPermute=10, seeds=NULL, chr="chr21", cutoff=c(2, 5), mc.cores=4, verbose=FALSE),
#' times=10)
#' levels(micro$expr) <- c("one", "four")
#' micro
#' ## Using 4 cores doesn't help with this toy data, but it will (at the expense of more RAM) if you have a larger data set.
#' }

calculatePvalues <- function(coveragePrep, models, fstats, nPermute = 1L, seeds = as.integer(gsub("-", "", Sys.Date())) + seq_len(nPermute), chr, maxGap = 300L, cutoff = quantile(fstats, 0.99), mc.cores=getOption("mc.cores", 2L), verbose=TRUE, significantCut=c(0.05, 0.10)) {
	## Setup
	if(is.null(seeds)) {
		seeds <- rep(NA, nPermute)
	}
	stopifnot(nPermute == length(seeds))
	stopifnot(length(intersect(names(coveragePrep), c("coverageSplit", "position", "meanCoverage"))) == 3)
	stopifnot(length(intersect(names(models), c("mod", "mod0"))) == 2)
	stopifnot(length(significantCut) == 2 & all(significantCut >=0 & significantCut <=1))
	
	## Identify the clusters
	if(verbose) message(paste(Sys.time(), "calculatePvalues: identifying clusters"))
	position <- coveragePrep$position
	means <- coveragePrep$meanCoverage
	cluster <- clusterMakerRle(position, maxGap)
	
	## Find the regions
	regs <- findRegions(position, chr=chr, fstats=fstats, cluster=cluster, y=fstats, cutoff=cutoff, verbose=verbose) 
	regs$meanCoverage <- mean(Views(means, IRanges(start=regs$indexStart, end=regs$indexEnd)))
	rm(fstats, position)
	
	
	## Pre-allocate memory
	nullwidths <- nullstats <- vector("list", length(seeds) * 2)
	last <- 0
	nSamples <- seq_len(nrow(models$mod))
	coverageSplit <- coveragePrep$coverageSplit
		
	for(i in seq_along(seeds)) {
		if(verbose) message(paste(Sys.time(), "calculatePvalues: calculating F-statistics for permutation", i))		
		
		if(!is.na(seeds[i])) {
			set.seed(seeds[i])
		}
		idx.permute <- sample(nSamples)
		
		## Permuted sample labels
		mod.p <- models$mod[idx.permute, , drop=FALSE]
		mod0.p <- models$mod0[idx.permute, , drop=FALSE]
		
		## Get the F-statistics
		fstats.output <- mclapply(coverageSplit, fstats.apply, mod=mod.p, mod0=mod0.p, mc.cores=mc.cores)
		fstats.output <- unlist(RleList(fstats.output), use.names=FALSE)	
			
		## Find the segments
		Indexes <- getSegmentsRle(x=fstats.output, f=cluster, cutoff=cutoff, verbose=verbose, zero=FALSE)
		
		## Calculate mean statistics
	    for (j in 1:2) {
			view <- Views(fstats.output, Indexes[[j]])
			nullstats[[last + j]] <- Rle(mean(view))
			nullwidths[[last + j]] <- Rle(width(view))
	    }		
		last <- last + 2
		
		## Finish loop
		rm(idx.permute, fstats.output, Indexes, view, mod.p, mod0.p)
		
	}
	nullstats <- do.call(c, nullstats)
	nullwidths <- do.call(c, nullwidths)
	
	if(length(nullstats) > 0) {
		## Proceed only if there is at least one null stats
		nullareas <- as.numeric(abs(nullstats) * nullwidths)
		rm(coveragePrep, coverageSplit)
	
	
		## Calculate pvalues
		if(verbose) message(paste(Sys.time(), "calculatePvalues: calculating the p-values"))
		pvals <- sapply(regs$area, function(x) { sum(nullareas > x) })
		regs$pvalues <- (pvals + 1) / (length(nullareas) + 1)
		regs$qvalues <- qvalue(regs$pvalues)$qvalues
		regs$significant <- factor(regs$pvalues < significantCut[1], levels=c(TRUE, FALSE))
		regs$significantQval <- factor(regs$qvalues < significantCut[2], levels=c(TRUE, FALSE))
		regs <- regs[order(regs$area, decreasing=TRUE), ]
	} else {
		if(verbose) message(paste(Sys.time(), "calculatePvalues: no null regions found. Skipping p-value calculation."))
	}
	## Save the nullstats too
	final <- list(regions=regs, nullstats=nullstats, nullwidths=nullwidths)
	
	
	## Done =)
	return(final)	
}
