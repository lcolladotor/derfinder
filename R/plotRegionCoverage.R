#' Makes plots for every region while summarizing the annotation
#'
#' This function takes the regions found in \link{calculatePvalues} and assigns them genomic states contructed with \link{makeGenomicState}. The main workhorse functions are \link[IRanges]{countOverlaps} and \link[IRanges]{findOverlaps}. For an alternative plot check \link{plotCluster} which is much slower and we recommend it's use only after quickly checking the results with this function.
#' 
#' @param regions The \code{$regions} output from \link{calculatePvalues}.
#' @param regionCoverage The output from \link{getRegionCoverage} used on \code{regions}.
#' @param groupInfo A factor specifying the group membership of each sample. It will be used to color the samples by group.
#' @param nearestAnnotation The output from \link[bumphunter]{annotateNearest} used on \code{regions}.
#' @param annotatedRegions The output from \link{annotateRegions} used on \code{regions}.
#' @param N The maximum number of regions to plot: will only make less than \code{N} plots only if \code{N} is greater than the number of regions.
#' @param colors If \code{NULL} then \link[RColorBrewer]{brewer.pal} with the \code{"Dark2"} color scheme is used.
#' @param scalefac The parameter used in \link{preprocessCoverage}.
#' @param ask If \code{TRUE} then the user is prompted before each plot is made.
#'
#' @return A plot for every region showing the coverage of each sample at each base of the region as well as the summarized annotation information.
#'
#' @author Andrew Jaffe, Leonardo Collado-Torres
#' @seealso \link{calculatePvalues}, \link{getRegionCoverage}, \link[bumphunter]{annotateNearest}, \link{annotateRegions}, \link{plotCluster}
#' @export
#' @importMethodsFrom GenomicRanges mcols names start end "$" "[[" as.data.frame
#'
#' @examples
#' \dontrun{
#' ## Construct the models
#' group <- genomeInfo$pop
#' adjustvars <- data.frame(genomeInfo$gender)
#' models <- makeModels(coverageInfo=genomeData, testvars=group, adjustvars=adjustvars, nonzero=TRUE)
#'
#' ## Preprocess the data
#' ## Automatic chunksize used to then compare 1 vs 4 cores in the 'do not run' section
#' prep <- preprocessCoverage(genomeData, groupInfo=group, cutoff=0, scalefac=32, chunksize=NULL, colsubset=NULL, mc.cores=4)
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
#'
#' ## Create GenomicState object:
#' ## Hsapiens.UCSC.hg19.knownGene GenomicState
#' library("TxDb.Hsapiens.UCSC.hg19.knownGene")
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#'
#' ## Creating this GenomicState object takes around 8 min
#' GenomicState.Hsapiens.UCSC.hg19.knownGene <- makeGenomicState(txdb=txdb)
#' 
#' ## Annotate regions
#' annotatedRegions <- annotateRegions(regions=regsWithP$regions, genomicState=GenomicState.Hsapiens.UCSC.hg19.knownGene, minoverlap=1)
#'
#' ## Find nearest annotation
#' library("bumphunter")
#' nearestAnnotation <- annotateNearest(regsWithP$regions, "hg19")
#'
#' ## Assign chr lengths using hg19 information
#' library("GenomicRanges")
#' data(hg19Ideogram, package = "biovizBase", envir = environment())
#' regions <- regsWithP$regions
#' seqlengths(regions) <- seqlengths(hg19Ideogram)[names(seqlengths(regions))]
#'
#' ## Obtain fullCov object
#' datadir <- system.file("extdata", "genomeData", package="derfinder2")
#' dirs <- makeBamList(datadir=datadir, samplepatt="*accepted_hits.bam$", bamterm=NULL)
#' ## Shorten the column names
#' names(dirs) <- gsub("_accepted_hits.bam", "", names(dirs))
#' 
#' ## Reading the data and filtering it is quite fast.
#' fullCov <- fullCoverage(dirs=dirs, chrnums="21", mc.cores=1)
#'
#' ## Get the region coverage
#' regionCov <- getRegionCoverage(fullCov=fullCov, regions=regions)
#'
#' ## Make plots for the regions
#' plotRegionCoverage(regions=regions, regionCoverage=regionCov, groupInfo=group, nearestAnnotation=nearestAnnotation, annotatedRegions=annotatedRegions, N=2, ask=TRUE)
#' 
#' ## If you prefer, you can save the plots to a pdf file
#' pdf("ders.pdf", h = 6, w = 9)
#' plotRegionCoverage(regions=regions, regionCoverage=regionCov, groupInfo=group, nearestAnnotation=nearestAnnotation, annotatedRegions=annotatedRegions)
#' dev.off()
#' }

plotRegionCoverage <- function(regions, regionCoverage, groupInfo, nearestAnnotation, annotatedRegions, N = 100, colors=NULL, scalefac = 32, ask = interactive()) {
	stopifnot(length(intersect(names(annotatedRegions), c("annotationList"))) == 1)
	stopifnot(length(intersect(names(regionCoverage), c("coverageData"))) == 1)
	stopifnot(is.data.frame(nearestAnnotation) | is(nearestAnnotation, "GRanges"))
	if(is.data.frame(nearestAnnotation)) {
		stopifnot(length(intersect(colnames(nearestAnnotation), c("name", "distance", "region"))) == 3)
	} else {
		stopifnot(length(intersect(names(mcols(nearestAnnotation)), c("name", "distance", "region"))) == 3)
	}
	stopifnot(is.factor(groupInfo))
	
	if(is.null(colors)) {
		library("RColorBrewer")
		palette(brewer.pal(length(levels(groupInfo)), "Dark2"))
	}
	
	anno <- annotatedRegions$annotationList

	layout(matrix(c(1, 1, 2), ncol = 1))
	N <- min(N, length(regions))
	for(i in seq_len(N)) {
		if(i %% 10 == 0) cat(".")
		y <- log2(regionCoverage$coverageData[[i]] + scalefac)
		x <- start(regions[i]):end(regions[i])
		
		if(ask) {
			devAskNewPage(TRUE)
		}
		
		par(mar=c(0, 4.5, 0.25, 1.1), oma=c(0,0,2,0))

		matplot(x, y, lty=1, col = as.numeric(groupInfo), type="l", yaxt="n", ylab="", xlab="", xaxt = "n", cex.lab = 1.7)
		m <- ceiling(max(y))
		axis(2, at = 5:m, labels = 2^(5:m) - scalefac, cex.axis = 1.5)
		
		legend("topleft", levels(groupInfo), pch = 15, col=seq(along=levels(groupInfo)), ncol = length(levels(groupInfo)), cex = 1.2, pt.cex = 1.5)
		mtext("Coverage", side = 2, line = 2.5, cex = 1.3)
		mtext(paste(nearestAnnotation$name[i], ",", nearestAnnotation$distance[i], "bp from tss:", nearestAnnotation$region[i]), outer = TRUE, cex = 1.3)
		## annotation
		par(mar=c(3.5, 4.5, 0.25, 1.1))
		plot(0, 0, type="n", xlim=range(x), ylim=c(-1.5, 1.5), yaxt="n", ylab="", xlab="", cex.axis = 1.5, cex.lab = 1.5)
		a <- as.data.frame(anno[[i]])
		Strand <- ifelse(a$strand == "+", 1, ifelse(a$strand == "-", -1, 0))
		Col <- ifelse(a$theRegion == "exon", "blue", ifelse(a$theRegion == "intron", "lightblue", "white"))
		Lwd <- ifelse(a$theRegion == "exon", 1, ifelse(a$theRegion == "intron", 0.5, 0))
		axis(2, c(-1, 1) , c("-", "+"), tick = FALSE, las = 1, cex.axis = 3)
		abline(h = 0, lty = 3)
		for(j in seq_len(nrow(a))) {
			polygon(c(a$start[j], a$end[j], a$end[j], a$start[j]), Strand[j]/2 + c(-0.3, -0.3, 0.3, 0.3) * Lwd[j], col = Col[j])
		}
		e <- a[a$theRegion=="exon",]
		s2 <- Strand[a$theRegion=="exon"]
		g  <- sapply(e$tx_name, paste, collapse="\n")
		if(length(g) > 0) text(x = e$start + e$width/2, y = s2 * -1, g, font = 2, pos = s2+2)
		mtext("Genes", side = 2, line = 2.5, cex = 1.3)
		mtext(unique(a$seqnames), side=1, line = 2.2, cex = 1.1)
	}
	return(invisible(NULL))
}
