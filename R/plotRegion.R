#' Plot the coverage information surrounding a region
#'
#' For a given region found in \link{calculatePvalues}, plot the coverage in it's vicinity as well as the mean by group. If annotation exists, you can plot the trancripts and exons (if any) overlapping in the vicinity of the region of interest.
#'
#' 
#' @param idx A integer specifying the index number of the region of interest.
#' @param regions The \code{$regions} output from \link{calculatePvalues}.
#' @param annotation The output from running \link[bumphunter]{annotateNearest} on the output from \link{calculatePvalues}.
#' @param coverageInfo A DataFrame resulting from \link{loadCoverage} using \code{cutoff=NULL}.
#' @param groupInfo A factor specifying the group membership of each sample.
#' @param titleUse Whether to show the p-value (\code{pval}) or the q-value (\code{qval}) in the title.
#' @param txdb A transcript data base such as TxDb.Hsapiens.UCSC.hg19.knownGene If \code{NULL} then no annotation information is used.
#' @param p.ideogram If \code{NULL}, the ideogram for hg19 is built for the corresponding chromosome. Otherwise an ideogram resuling from \link[ggbio]{plotIdeogram}.
#' @param minExtend The minimum number of base-pairs to extend the view before and after the region of interest.
#' @param scalefac A log transformation is used on the count tables, so zero counts present a problem. Ideally, the same you provided to \link{preprocessCoverage}.
#'
#' @return A ggplot2 plot that is ready to be printed out. Tecnically it is a ggbio object.
#'
#' @details See the parameter \code{significantCut} in \link{calculatePvalues} for how the significance cutoffs are determined.
#'
#' @seealso \link{loadCoverage}, \link{calculatePvalues}, \link[bumphunter]{annotateNearest}, \link[ggbio]{plotIdeogram}
#' @author Leonardo Collado-Torres
#' @export
#' @importFrom IRanges width resize
#' @importMethodsFrom IRanges "$" "[" as.matrix findOverlaps
#' @importFrom GenomicRanges seqnames
#' @importMethodsFrom GenomicRanges findOverlaps start end as.data.frame
#' @importFrom ggbio plotIdeogram tracks theme_tracks_sunset
#' @importMethodsFrom ggbio autoplot
#' @importFrom reshape2 melt
#' @importFrom ggplot2 autoplot aes scale_fill_manual ggplot geom_line ylab
#' @importFrom plyr ddply summarise
#' @examples
#' ## Get coverage info without any cutoff
#' datadir <- system.file("extdata", "genomeData", package="derfinder2")
#' dirs <- makeBamList(datadir=datadir, samplepatt="*accepted_hits.bam$", bamterm=NULL)
#' names(dirs) <- gsub("_accepted_hits.bam", "", names(dirs))
#' covInfo <- loadCoverage(dirs=dirs, chr="21", cutoff=NULL, verbose=FALSE)
#' 
#' ## Construct the models
#' group <- genomeInfo$pop
#' adjustvars <- data.frame(genomeInfo$gender)
#' models <- makeModels(coverageInfo=genomeData, group=group, adjustvars=adjustvars, nonzero=TRUE)
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
#' plotRegion(idx=1, regions=regsWithP$regions, annotation=annotation, coverageInfo=covInfo$coverage, groupInfo=group, txdb=TxDb.Hsapiens.UCSC.hg19.knownGene)
#' ## Resize the plot window and the labels will look good.
#'
#' \dontrun{
#' ## For a custom plot, check the ggbio and ggplot2 packages.
#' ## Also feel free to look at the code for this function:
#' plotRegion
#' }

plotRegion <- function(idx, regions, annotation, coverageInfo, groupInfo, titleUse="pval", txdb=NULL, p.ideogram=NULL, minExtend=200, scalefac=32) {
	stopifnot(titleUse %in% c("pval", "qval"))
	
	## Satisfying R CMD check
	significant <- significantQval <- position <- valueScaled <- variable <- group <- value <- meanScaled <- NULL
		
	## Select region and build title
	l <-  width(regions[idx]) + 2 * max(minExtend, width(regions[idx]))
	wh <- resize(regions[idx], l, fix="center")
	if(titleUse == "pval") {
		title <- paste0("Annotated name ", annotation$name[idx], " with p-value ", regions$pvalues[idx],", sf=", scalefac)
	} else {
		title <- paste0("Annotated name ", annotation$name[idx], " with q-value ", regions$qvalues[idx],", sf=", scalefac)
	}
	
	## Plot the ideogram if not supplied
	if(is.null(p.ideogram)) {
		chr <- as.character(seqnames(wh))
		p.ideogram <- plotIdeogram(genome = "hg19", subchr=chr)
	}
	
	## Regions found (from the view)
	if(titleUse == "pval") {
		p.region <- autoplot(regions[as.matrix(findOverlaps(regions, wh))[, "queryHits"]], aes(fill=significant)) + scale_fill_manual(values=c("chartreuse4", "wheat2"), limits=c("TRUE", "FALSE"))
	} else {
		p.region <- autoplot(regions[as.matrix(findOverlaps(regions, wh))[, "queryHits"]], aes(fill=significantQval)) + scale_fill_manual(values=c("chartreuse4", "wheat2"), limits=c("TRUE", "FALSE"))
	}
	

	## Construct the coverage plot
	pos <- start(wh):end(wh)
	rawData <- as.data.frame(coverageInfo[pos, ])
	rawData$position <- pos
	covData <- melt(rawData, id.vars="position")
	covData$group <- rep(groupInfo, each=nrow(rawData))
	covData$valueScaled <- log2(covData$value + scalefac)
	p.coverage <- ggplot(covData, aes(x=position, y=valueScaled, group=variable, colour=group)) + geom_line(alpha=1/2, size=1.5)
	
	## Construct mean by group coverage plot
	meanCov <- ddply(covData, c("position", "group"), summarise, meanCov=mean(value))
	meanCov$meanScaled <- log2(meanCov$meanCov + scalefac)
	p.meanCov <- ggplot(meanCov, aes(x=position, y=meanScaled, colour=group)) + geom_line(alpha=1/2, size=1.5)
	
	## Annotation info and final plot
	if(is.null(txdb)) {
		p.transcripts <- FALSE
	} else {
		## The tryCatch is needed because not all regions overlap a transcript
		p.transcripts <- tryCatch(autoplot(txdb, which = wh, names.expr = "tx_name(gene_id)"), warning=function(w) { FALSE })
	}	
	if(!is.logical(p.transcripts)) {
		## tryCatch needed because if the region is smaller than an exon, the stat="reduce" fails
		p.exons <- tryCatch(autoplot(txdb, which = wh, stat = "reduce", color = "brown", fill = "brown"), error=function(w) { FALSE })
		if(!is.logical(p.exons)) {
			result <- tracks(p.ideogram, "Coverage\nlog2(x + sf)" = p.coverage, "Mean coverage\nlog2(xbar + sf)" = p.meanCov, "Regions" = p.region, "Known\ntx_name(gene_id)" = p.transcripts, "Exons" = p.exons, heights = c(2, 4, 4, 1.5, 3, 1), xlim=wh, title=title) + ylab("") + theme_tracks_sunset()
		} else {
			result <- tracks(p.ideogram, "Coverage\nlog2(x + sf)" = p.coverage, "Mean coverage\nlog2(xbar + sf)" = p.meanCov, "Regions" = p.region, "Known\ntx_name(gene_id)" = p.transcripts, heights = c(2, 4, 4, 1.5, 3), xlim=wh, title=title) + ylab("") + theme_tracks_sunset()
		}
		
		
	} else {
		result <- tracks(p.ideogram, "Coverage\nlog2(x + sf)" = p.coverage, "Mean coverage\nlog2(xbar + sf)" = p.meanCov, "Regions" = p.region, heights = c(1.5, 5, 5, 2), xlim=wh, title=title) + ylab("") + theme_tracks_sunset()
	}
	return(result)	
}
