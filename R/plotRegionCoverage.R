#' Makes plots for every region while summarizing the annotation
#'
#' This function takes the regions found in \link{calculatePvalues} and assigns them genomic states contructed with \link{makeGenomicState}. The main workhorse functions are \link[IRanges]{countOverlaps} and \link[IRanges]{findOverlaps}. For an alternative plot check \link{plotCluster} which is much slower and we recommend it's use only after quickly checking the results with this function.
#' 
#' @param regions The \code{$regions} output from \link{calculatePvalues}.
#' @param regionCoverage The output from \link{getRegionCoverage} used on \code{regions}.
#' @param groupInfo A factor specifying the group membership of each sample. It will be used to color the samples by group.
#' @param nearestAnnotation The output from \link[bumphunter]{annotateNearest} used on \code{regions}.
#' @param annotatedRegions The output from \link{annotateRegions} used on \code{regions}.
#' @param whichRegions An integer vector with the index of the regions to plot.
#' @param colors If \code{NULL} then \link[RColorBrewer]{brewer.pal} with the \code{"Dark2"} color scheme is used.
#' @param scalefac The parameter used in \link{preprocessCoverage}.
#' @param ask If \code{TRUE} then the user is prompted before each plot is made.
#' @param ylab The name of the of the Y axis.
#' @param verbose If \code{TRUE} basic status updates will be printed along the way.
#'
#' @return A plot for every region showing the coverage of each sample at each base of the region as well as the summarized annotation information.
#'
#' @author Andrew Jaffe, Leonardo Collado-Torres
#' @seealso \link{calculatePvalues}, \link{getRegionCoverage}, \link[bumphunter]{annotateNearest}, \link{annotateRegions}, \link{plotCluster}
#' @export
#' @importMethodsFrom GenomicRanges mcols names start end "$" "[[" as.data.frame
#'
#' @examples
#' ## Annotate regions, first two regions only
#' regions <- genomeRegions$regions[1:2]
#' annotatedRegions <- annotateRegions(regions=regions, genomicState=genomicState, minoverlap=1)
#'
#' ## Find nearest annotation
#' library("bumphunter")
#' nearestAnnotation <- annotateNearest(regions, "hg19")
#'
#' ## Obtain fullCov object
#' fullCov <- list("21"=genomeDataRaw$coverage)
#' 
#' ## Assign chr lengths using hg19 information
#' library("GenomicRanges")
#' data(hg19Ideogram, package = "biovizBase", envir = environment())
#' seqlengths(regions) <- seqlengths(hg19Ideogram)[names(seqlengths(regions))]
#'
#' ## Get the region coverage
#' regionCov <- getRegionCoverage(fullCov=fullCov, regions=regions)
#'
#' ## Make plots for the regions
#' plotRegionCoverage(regions=regions, regionCoverage=regionCov, groupInfo=genomeInfo$pop, nearestAnnotation=nearestAnnotation, annotatedRegions=annotatedRegions, whichRegions=1:2)
#' 
#' \dontrun{
#' ## If you prefer, you can save the plots to a pdf file
#' pdf("ders.pdf", h = 6, w = 9)
#' plotRegionCoverage(regions=regions, regionCoverage=regionCov, groupInfo=genomeInfo$pop, nearestAnnotation=nearestAnnotation, annotatedRegions=annotatedRegions, whichRegions=1:2, ask=FALSE)
#' dev.off()
#' }

plotRegionCoverage <- function(regions, regionCoverage, groupInfo, nearestAnnotation, annotatedRegions, whichRegions = seq_len(min(100, length(regions))), colors=NULL, scalefac = 32, ask = interactive(), ylab="Coverage", verbose=TRUE) {
	stopifnot(length(intersect(names(annotatedRegions), c("annotationList"))) == 1)
	stopifnot(is.data.frame(nearestAnnotation) | is(nearestAnnotation, "GRanges"))
	if(is.data.frame(nearestAnnotation)) {
		stopifnot(length(intersect(colnames(nearestAnnotation), c("name", "distance", "region"))) == 3)
	} else {
		stopifnot(length(intersect(names(mcols(nearestAnnotation)), c("name", "distance", "region"))) == 3)
	}
	stopifnot(is.factor(groupInfo))
	
	## Make sure that the user is not going beyond the number of regions
	N <- max(whichRegions)
	if(N > length(regions)) {
		warning("'whichRegions' exceeds the number of regions available. Dropping invalid indexes.")
		whichRegions <- whichRegions[whichRegions <= length(regions)]
	}
	
	## Color setup
	if(is.null(colors)) {
		library("RColorBrewer")
		palette(brewer.pal(max(3, length(levels(groupInfo))), "Dark2"))
	}
	
	## Annotation information
	anno <- annotatedRegions$annotationList

	layout(matrix(c(1, 1, 2), ncol = 1))
	for(i in whichRegions) {
		## Status update
		if((i - 1) %% 10 == 0 & verbose) {
			par(mar=c(5, 4, 4, 2) + 0.1, oma=c(0, 0, 0, 0))
			plot.new()
			text(0.5, 0.5, i, cex=20)
		}
		
		## For subsetting the named lists
		ichar <- as.character(i)
		
		## Obtain data
		y <- log2(regionCoverage[[ichar]] + scalefac)
		x <- start(regions[i]):end(regions[i])
		
		if(ask) {
			devAskNewPage(TRUE)
		}
		
		## Plot coverage
		par(mar=c(0, 4.5, 0.25, 1.1), oma=c(0, 0, 2, 0))
		if(length(x) > 1) {
			matplot(x, y, lty=1, col = as.numeric(groupInfo), type="l", yaxt="n", ylab="", xlab="", xaxt = "n", cex.lab = 1.7)
		} else {
			matplot(x, y, lty=1, col = as.numeric(groupInfo), type="p", yaxt="n", ylab="", xlab="", xaxt = "n", cex.lab = 1.7)
		}
		m <- ceiling(max(y))
		y.labs <- seq(from=0, to=log2(2^m - scalefac), by=1)
		axis(2, at = log2(scalefac + c(0, 2^y.labs)), labels = c(0, 2^y.labs), cex.axis = 1.5)
		
		## Coverage legend
		legend("topleft", levels(groupInfo), pch = 15, col=seq(along=levels(groupInfo)), ncol = length(levels(groupInfo)), cex = 1.5, pt.cex = 1.5)
		mtext(ylab, side = 2, line = 2.5, cex = 1.3)
		mtext(paste(nearestAnnotation$name[i], ",", nearestAnnotation$distance[i], "bp from tss:", nearestAnnotation$region[i]), outer = TRUE, cex = 1.3)
		
		## Plot annotation
		par(mar=c(3.5, 4.5, 0.25, 1.1))
		xrange <- range(x)
		if(length(x) == 1) {
			xrange <- xrange + c(-1, 1)
		}
		plot(0, 0, type="n", xlim=xrange, ylim=c(-1.5, 1.5), yaxt="n", ylab="", xlab="", cex.axis = 1.5, cex.lab = 1.5)
		gotAnno <- !is.null(anno[[ichar]])
		if(gotAnno) {
			a <- as.data.frame(anno[[ichar]])
			Strand <- ifelse(a$strand == "+", 1, ifelse(a$strand == "-", -1, 0))
			Col <- ifelse(a$theRegion == "exon", "blue", ifelse(a$theRegion == "intron", "lightblue", "white"))
			Lwd <- ifelse(a$theRegion == "exon", 1, ifelse(a$theRegion == "intron", 0.5, 0))
		}		
		axis(2, c(-1, 1) , c("-", "+"), tick = FALSE, las = 1, cex.axis = 3)
		abline(h = 0, lty = 3)
		if(gotAnno) {
			for(j in seq_len(nrow(a))) {
				polygon(c(a$start[j], a$end[j], a$end[j], a$start[j]), Strand[j]/2 + c(-0.3, -0.3, 0.3, 0.3) * Lwd[j], col = Col[j])
			}
			e <- a[a$theRegion=="exon",]
			s2 <- Strand[a$theRegion=="exon"]
			g <- unlist(e$symbol)
			if(!is.null(g)) {
				g[is.na(g)] <- ""
				if(length(g) > 0) text(x=e$start + e$width/2, y=s2*0.9, g,font=2,pos=s2+2,cex=2)
			}
			
		}		
		mtext("Genes", side = 2, line = 2.5, cex = 1.3)
		mtext(as.character(seqnames(regions)[i]), side=1, line = 2.2, cex = 1.1)

	}
	return(invisible(NULL))
}
