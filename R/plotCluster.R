#' Plot the coverage information surrounding a region cluster
#'
#' For a given region found in \link{calculatePvalues}, plot the coverage for the cluster this region belongs to as well as some padding. The mean by group is shown to facilitate comparisons between groups. If annotation exists, you can plot the trancripts and exons (if any) overlapping in the vicinity of the region of interest.
#'
#' 
#' @param idx A integer specifying the index number of the region of interest.
#' @param regions The \code{$regions} output from \link{calculatePvalues}.
#' @param annotation The output from running \link[bumphunter]{annotateNearest} on the output from \link{calculatePvalues}.
#' @param coverageInfo A DataFrame resulting from \link{loadCoverage} using \code{cutoff=NULL}.
#' @param groupInfo A factor specifying the group membership of each sample. It will be used to color the samples by group.
#' @param titleUse Whether to show the p-value (\code{pval}) or the q-value (\code{qval}) in the title. If \code{titleUse=none} then no p-value or q-value information is used; useful if no permutations were performed and thus p-value and q-value information is absent.
#' @param txdb A transcript data base such as TxDb.Hsapiens.UCSC.hg19.knownGene If \code{NULL} then no annotation information is used.
#' @param p.ideogram If \code{NULL}, the ideogram for hg19 is built for the corresponding chromosome. Otherwise an ideogram resuling from \link[ggbio]{plotIdeogram}.
#' @param maxExtend The maximum number of base-pairs to extend the view (on each side) before and after the region cluster of interest. For small region clusters, the one side extension is equal to the width of the region cluster.
#' @param colsubset Column subset in case that it was specified in \link{preprocessCoverage}.
#' @param forceLarge If \code{TRUE} then the data size limitations are ignored. The window size (region cluster width + 2 times \code{maxExtend}) has to be less than 100 kb. Note that a single plot at the 300kb range can take around 2 hours to complete.
#'
#' @return A ggplot2 plot that is ready to be printed out. Tecnically it is a ggbio object. The region with the red bar is the one whose information is shown in the title.
#'
#' @details See the parameter \code{significantCut} in \link{calculatePvalues} for how the significance cutoffs are determined.
#'
#' @seealso \link{loadCoverage}, \link{calculatePvalues}, \link[bumphunter]{annotateNearest}, \link[ggbio]{plotIdeogram}
#' @author Leonardo Collado-Torres
#' @export
#' @importFrom IRanges width resize
#' @importMethodsFrom IRanges "$" "[" as.matrix findOverlaps queryHits
#' @importFrom GenomicRanges seqnames
#' @importMethodsFrom GenomicRanges findOverlaps start end as.data.frame range
#' @importFrom ggbio plotIdeogram tracks theme_tracks_sunset
#' @importMethodsFrom ggbio autoplot
#' @importFrom reshape2 melt
#' @importFrom ggplot2 autoplot aes scale_fill_manual ggplot geom_line ylab guides scale_y_continuous geom_segment
#' @importFrom plyr ddply summarise
#' @importFrom scales log2_trans log_trans
#' @examples
#' ## Get coverage info without any cutoff
#' datadir <- system.file("extdata", "genomeData", package="derfinder")
#' dirs <- makeBamList(datadir=datadir, samplepatt="*accepted_hits.bam$", bamterm=NULL)
#' names(dirs) <- gsub("_accepted_hits.bam", "", names(dirs))
#' covInfo <- loadCoverage(dirs=dirs, chr="21", cutoff=NULL, verbose=FALSE)
#' 
#' ## Construct the models
#' group <- genomeInfo$pop
#' adjustvars <- data.frame(genomeInfo$gender)
#' models <- makeModels(coverageInfo=genomeData, testvars=group, adjustvars=adjustvars, nonzero=TRUE)
#'
#' ## Preprocess the data
#' prep <- preprocessCoverage(genomeData, cutoff=0, scalefac=32, chunksize=NULL, colsubset=NULL, mc.cores=4)
#' 
#' ## Get the F statistics
#' fstats <- calculateStats(prep, models, mc.cores=1, verbose=FALSE)
#'
#' ## Determine a cutoff from the F-distribution.
#' ## This step is very important and you should consider using quantiles from the observed F statistics
#' \dontrun{
#' n <- dim(prep$coverageSplit[[1]])[2]
#' df1 <- dim(models$mod)[2]
#' df0 <- dim(models$mod0)[2]
#' cutoff <- qf(0.95, df1-df0, n-df1)
#' }
#' ## Low cutoff used for illustrative purposes
#' cutoff <- 1
#'
#' ## Calculate the p-values and define the regions of interest.
#' regsWithP <- calculatePvalues(prep, models, fstats, nPermute=10, seeds=NULL, chr="chr21", cutoff=cutoff, mc.cores=1, verbose=FALSE)
#'
#' ## Annotate the results
#' suppressMessages(library("bumphunter"))
#' annotation <- annotateNearest(regsWithP$regions, "hg19")
#'
#' ## Make the plot
#' suppressMessages(library("TxDb.Hsapiens.UCSC.hg19.knownGene"))
#' plotCluster(idx=1, regions=regsWithP$regions, annotation=annotation, coverageInfo=covInfo$coverage, groupInfo=group, txdb=TxDb.Hsapiens.UCSC.hg19.knownGene)
#' ## Resize the plot window and the labels will look good.
#'
#' \dontrun{
#' ## For a custom plot, check the ggbio and ggplot2 packages.
#' ## Also feel free to look at the code for this function:
#' plotCluster
#' 
#' #### This is a detailed example for a specific cluster
#' ## The purpose is to illustrate how data filtering (and availability), 
#' ## F-stat cutoff, cluster cutoff interact into determing the actual DERs.
#'
#' ## Using as example cluster #7
#' ## Note how despite having data and using a very small F-stat cutoff, some regions with data are split into different DERs
#' plotCluster(idx=7, regions=regsWithP$regions, annotation=annotation, coverageInfo=covInfo$coverage, groupInfo=group, txdb=TxDb.Hsapiens.UCSC.hg19.knownGene)
#' 
#' ## Identify DERs clusters and regions of the genome where we have data
#' clusters <- clusterMakerRle(prep$position, ranges=TRUE)
#' dataRegions <- clusterMakerRle(prep$position, maxGap=0, ranges=TRUE)
#' 
#' ## Apply F-stat cutoff of 1
#' segs <- getSegmentsRle(fstats, 1)$upIndex
#'
#' ## Separate the sections that passed the F-stat cutoff by regions in the genome
#' library("IRanges")
#' pieces <- disjoin(c(segs, dataRegions))
#' 
#' ## The DERs are actually the following ones:
#' ders <- pieces[queryHits(findOverlaps(pieces, segs))]
#' ## You can very that this is the case:
#' identical(width(ders), width(sort(ranges(regsWithP$regions))))
#' 
#' ## Ranges plotting function (from IRanges documentation)
#' plotRanges <- function(x, xlim = x, main = deparse(substitute(x)), col = "black", 
#'    sep = 0.5, ...) {
#'    height <- 1
#'    if (is(xlim, "Ranges")) 
#'        xlim <- c(min(start(xlim)), max(end(xlim)))
#'    bins <- disjointBins(IRanges(start(x), end(x) + 1))
#'    plot.new()
#'    plot.window(xlim, c(0, max(bins) * (height + sep)))
#'    ybottom <- bins * (sep + height) - height
#'    rect(start(x) - 0.5, ybottom, end(x) + 0.5, ybottom + height, col = col, 
#'        ...)
#'    title(main)
#'    axis(1)
#' }
#'
#' ## Visualize the different DER clusters
#' plotRanges(clusters)
#' ## Note that region 7 is part of cluster #4.
#' regsWithP$regions$cluster[7]
#' 
#' ## Plot the different segmentation steps, the final DERs, and the fstats with the cutoff
#' par(mfrow=c(5, 1))
#' plotRanges(dataRegions, xlim=c(200, 617))
#' plotRanges(segs, xlim=c(200, 617))
#' plotRanges(pieces, xlim=c(200, 617))
#' plotRanges(ders, xlim=c(200, 617))
#' f <- as.numeric(fstats)
#' plot(f, type="l", xlim=c(200, 617))
#' abline(h=1, col="red")
#' 
#' ## We can see that the 10 different data regions match with how 
#' ## many sections of the genome have data in plotCluster(idx=7, ...)
#'
#' ## The F-stat cutoff applied to F-stats leads to 6 different segments passing the cutoff.
#' ## We can see this both in the F-stat panel as well as in the segs panel.
#'
#' ## Between the regions with data and the segments, we have lots
#' ## of different pieces to take into account. 
#'
#' ## From the pieces, only 12 of them correspond to unique regions in
#' ## the genome that passed the F-stat cutoff.
#' 
#' ## They are the 12 different DERs we see in cluster 4 as shown in plotCluster(idx=7, ...)
#'
#' ## So despite using a very low F-stat cutoff, with the intention of getting
#' ## anything that had data (for illustrative purposes, in reality you will want 
#' ## to use a higher cutoff), some regions were not considered to be candidate DERs.
#'
#' ## From plotCluster(idx=7, ...) we see that the 4th and 5th region with data 
#' ## had sections (or none for the 4th one) that were considered DERs but others were not.
#' ## This is due to the higher variability in the CEU samples, and of course, the low 
#' ## F-stat cutoff used.  
#' }

plotCluster <- function(idx, regions, annotation, coverageInfo, groupInfo, titleUse="qval", txdb=NULL, p.ideogram=NULL, maxExtend=300L, colsubset=NULL, forceLarge=FALSE) {
	stopifnot(titleUse %in% c("pval", "qval", "none"))
	stopifnot(is.factor(groupInfo))
	if(is.null(colsubset)) colsubset <- seq_len(length(groupInfo))
	current <- regions[idx]
	
	## Satisfying R CMD check
	x <- xend <- y <- meanCov <- significant <- significantQval <- position <- valueScaled <- variable <- group <- value <- meanScaled <- NULL
		
	## Select region and build title
	cluster <- regions[ seqnames(regions) == seqnames(current) & regions$cluster == current$cluster]
	if(length(cluster) > 1) {
		cluster <- range(cluster)
	}
	l <-  width(cluster) + 2 * min(maxExtend, width(cluster))
	
	if(l > 1e5 & !forceLarge) {
		message(paste("No plot will be made because the data is too large. The window size exceeds 100 kb."))
		return(invisible(NULL))
	}
	
	wh <- resize(cluster, l, fix="center")
	if(titleUse == "pval") {
		title <- paste0("Cluster for region with name ", annotation$name[idx], " and p-value ", signif(current$pvalues, 4))
	} else if (titleUse == "qval") {
		title <- paste0("Cluster for region with name ", annotation$name[idx], " and q-value ", signif(current$qvalues, 4))
	} else {
		title <- paste0("Cluster for region with name ", annotation$name[idx])
	}
	
	## Plot the ideogram if not supplied
	if(is.null(p.ideogram)) {
		chr <- as.character(seqnames(wh))
		## Now load the ideogram info
		hg19IdeogramCyto <- NULL
		load(system.file("data", "hg19IdeogramCyto.rda", package="biovizBase", mustWork=TRUE))
		p.ideogram <- plotIdeogram(hg19IdeogramCyto, chr)
	}
	
	## Regions found (from the view)
	neighbors <- regions[queryHits(findOverlaps(regions, wh))]
	neighbors$originalRegion <- neighbors == current
	ann_line <- data.frame(x=start(current), xend=end(current), y=1)
	if(titleUse == "pval") {
		p.region <- autoplot(neighbors, aes(fill=significant)) + 
		scale_fill_manual(values=c("chartreuse4", "wheat2"), limits=c("TRUE", "FALSE")) + 
		geom_segment(aes(x=x, xend=xend, y=y, yend=y, size=3), data=ann_line, colour="red") + guides(size=FALSE)
	} else if (titleUse == "qval" ){
		p.region <- autoplot(neighbors, aes(fill=significantQval)) +
		scale_fill_manual(values=c("chartreuse4", "wheat2"), limits=c("TRUE", "FALSE")) +
		geom_segment(aes(x=x, xend=xend, y=y, yend=y, size=3), data=ann_line, colour="red") + guides(size=FALSE)
	} else {
		p.region <- autoplot(neighbors) + 
		geom_segment(aes(x=x, xend=xend, y=y, yend=y, size=3), data=ann_line, colour="red") + guides(size=FALSE)
	}

	## Graphical parameters
	nGroups <- length(levels(groupInfo))
	
	## Construct the coverage plot
	pos <- start(wh):end(wh)
	rawData <- as.data.frame(coverageInfo[pos, colsubset])
	rawData$position <- pos
	covData <- melt(rawData, id.vars="position")
	covData$group <- rep(groupInfo, each=nrow(rawData))
	p.coverage <- ggplot(covData, aes(x=position, y=value, group=variable, colour=group)) + geom_line(alpha=1/nGroups) + scale_y_continuous(trans=log2_trans())
	
	## Construct mean by group coverage plot
	meanCoverage <- ddply(covData, c("position", "group"), summarise, meanCov=mean(value))
	p.meanCov <- ggplot(meanCoverage, aes(x=position, y=meanCov, colour=group)) + geom_line(alpha=1/max(1, 1/2 * nGroups)) + scale_y_continuous(trans=log2_trans())
	
	## Annotation info and final plot
	if(is.null(txdb)) {
		p.transcripts <- FALSE
	} else {
		## The tryCatch is needed because not all regions overlap a transcript
		p.transcripts <- tryCatch(autoplot(txdb, which = wh, names.expr = "tx_name(gene_id)"), warning=function(w) { FALSE })
	}	
	if(!is.logical(p.transcripts)) {
		result <- tracks(p.ideogram, "Coverage" = p.coverage, "Mean coverage" = p.meanCov, "Regions" = p.region, "tx_name\n(gene_id)" = p.transcripts, heights = c(2, 4, 4, 1.5, 3), xlim=wh, title=title) + ylab("") + theme_tracks_sunset()		
	} else {
		result <- tracks(p.ideogram, "Coverage" = p.coverage, "Mean coverage" = p.meanCov, "Regions" = p.region, heights = c(2, 5, 5, 2), xlim=wh, title=title) + ylab("") + theme_tracks_sunset()
	}
	return(result)	
}
