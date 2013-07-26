#' Find non-zero regions in a Rle
#'
#' Find genomic regions for which a numeric vector is above (or below) predefined thresholds. This is similar to \link[bumphunter]{regionFinder} and is a helper function for \link{calculatePvalues}.
#' 
#' @param statsInfo A list with \code{$position} and \code{$fstats} components where the first one is a logical Rle of genomic positions and the second one is a numeric Rle. 
#' @param chr A single element character vector specifying the chromosome name.
#' @param fstats A numeric Rle with the F-statistics. Normally stored in \code{statsInfo$fstats}.
#' @param cluster The clusters of locations that are to be analyzed together, normally given by \link{clusterMakerRle}.
#' @param y A numeric Rle of the same length as \code{statsInfo$fstats} containing values to be averaged for the region. 
#' @param oneTable If \code{TRUE} only one results GRanges is returned. Otherwise, a GRangesList with two components is returned: one for the regions with positive values and one for the negative values.
#' @param maxGap This argument is passed to \link{clusterMakerRle}.
#' @param cutoff This argument is passed to \link{getSegmentsRle}.
#' @param verbose If \code{TRUE} basic status updates will be printed along the way.
#'
#' @return Either a GRanges or a GRangesList as determined by \code{oneTable}. Each of them has the following metadata variables.
#' \describe{
#' \item{value }{ The mean of the values of \code{y} for the given region.}
#' \item{area }{  The absolute value of the sum of the values of \code{y} for the given region.}
#' \item{indexStart }{ The start position of the region in terms of the index for \code{y}.}
#' \item{indexEnd }{ The end position of the region in terms of the index for \code{y}.}
#' \item{cluster }{ The cluser ID.}
#' \item{clusterL }{ The total length of the cluster.}
#' }
#'
#' @details \link[bumphunter]{regionFinder} adapted to Rle world.
#' @seealso \link[bumphunter]{regionFinder}, \link{calculatePvalues}
#' @references Rafael A. Irizarry, Martin Aryee, Hector Corrada Bravo, Kasper D. Hansen and Harris A. Jaffee. bumphunter: Bump Hunter. R package version 1.1.10.
#'
#' @author Leonardo Collado-Torres
#' @export
#' @importFrom IRanges IRanges start end width Views Rle runLength
#' @importFrom GenomicRanges GRanges GRangesList
#' @importMethodsFrom IRanges which length mean
#' @importMethodsFrom GenomicRanges unlist
#' @examples
#' ## Get the statistics
#' group <- brainInfo$outcome
#' adjustvars <- brainInfo[, c("sex", "age", "left.hemisph", "pmi", "brainpH")]
#' stats <- calculateStats(brainData, group, adjustvars=adjustvars, mc.cores=1, verbose=TRUE)
#'
#' ## Find the regions
#' regs <- findRegions(stats, "chr21", verbose=TRUE)
#' regs
#'
#' \dontrun{
#' ## Compare vs bumphunter
#' library("bumphunter")
#' regs2 <- regionFinder(as.numeric(stats$fstats), rep("chr21", length(stats$fstats)), which(stats$position), cluster=NULL, assumeSorted=TRUE, verbose=TRUE, order=FALSE)
#' regs2
#' ## Note that regs$L can be calculated with width(regs)
#' identical(width(regs), as.integer(regs2$L))
#' ## Time comparison
#' library("microbenchmark")
#' micro <- microbenchmark(findRegions(stats, "chr21", verbose=FALSE), regionFinder(as.numeric(stats$fstats), rep("chr21", length(stats$fstats)), which(stats$position), cluster=NULL, assumeSorted=TRUE, verbose=FALSE, order=FALSE))
#' levels(micro$expr) <- c("new", "original")
#' micro
#' ## The bumphunter function regionFinder() is faster in small data sets.
#'
#' ## Once you have the regions you can proceed to annotate them
#' annotation <- annotateNearest(regs, "hg19")
#' annotation
#' }

findRegions <- function(statsInfo, chr, fstats=statsInfo$fstats, cluster=NULL, y = fstats, oneTable = TRUE, maxGap = 300L, cutoff = quantile(abs(fstats), 0.99), verbose = TRUE) {
	stopifnot(length(intersect(names(statsInfo), c("position"))) == 1)
	
	## Identify the clusters
	if(is.null(cluster)) {
		if(verbose) message(paste(date(), "findRegions: identifying clusters"))
		cluster <- clusterMakerRle(statsInfo$position, maxGap)
	}	
	
	## Find the segments
	Indexes <- getSegmentsRle(x = statsInfo$fstats, f = cluster, cutoff = cutoff, verbose = verbose, zero=FALSE)
	
	## Sadly, this is required to map the positions of the index to the chr positions. It's 275 mb in RAM for a length of 72097604 instead of 4.7 Mb in Rle world.
	## The good thing is that it's temporary and the user will not need to save this
	pos <- which(statsInfo$position)
	
	## Build the output shell
	hasInfo <- sapply(Indexes, length) != 0
    res <- vector("list", sum(hasInfo))
	names(res) <- names(hasInfo)[hasInfo]
	
    for (i in names(hasInfo)[hasInfo]) {
		## Define the chr ranges
		pos.ir <- IRanges(start=pos[start(Indexes[[i]])], width = width(Indexes[[i]]) )
		
		## Extract info from y
		view <- Views(y, Indexes[[i]])
		clus <- mean(Views(cluster, Indexes[[i]]))
	
		## Actually build the GRanges
		res[[i]] <- GRanges(
			seqnames = Rle(chr, length(Indexes[[i]])),
			ranges = pos.ir,
			value = mean(view),
			area = abs(sum(view)),
            indexStart = start(Indexes[[i]]), 
            indexEnd = end(Indexes[[i]]),
			cluster = Rle(as.integer(clus)),
			clusterL = Rle(runLength(cluster)[clus])
		)
    }
	## Fix names and format
    names(res) <- gsub("Index", "", names(res))
	res <- GRangesList(res)
	
	## Finish up
    if (oneTable) {
        res <- unlist(res)
    }
    return(res)
	
}
