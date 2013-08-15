#' Plot a karyotype overview of the genome with the identified regions
#'
#' Plots an overview of the genomic locations of the identified regions (see \link{calculatePvalues}) in a karyotype view. The coloring can be done either by significant regions according to their p-values, significant by adjusted p-values, or by annotated region if using \link[bumphunter]{annotateNearest}.
#'
#' 
#' @param regions The \code{$regions} output from \link{calculatePvalues}.
#' @param annotation The output from running \link[bumphunter]{annotateNearest} on the output from \link{calculatePvalues}. It is only required if \code{type="annotation"}.
#' @param type Must be either \code{pval}, \code{qval} or \code{annotation}. It determines whether the plot coloring should be done according to significant p-values (<0.05), significant q-values (<0.10) or annotation regions.
#' @param base_size Base point size of the plot. This argument is passed to \link[ggplot2]{element_text} (\code{size} argument).
#' @param areaRel The relative size for the area label when \code{type="pval"} or \code{type="qval"}. Can be useful when making high resolution versions of these plots in devices like CairoPNG.
#'
#' @return A ggplot2 plot that is ready to be printed out. Tecnically it is a ggbio object.
#'
#' @seealso \link{calculatePvalues}, \link[bumphunter]{annotateNearest}
#' @author Leonardo Collado-Torres
#' @export
#' @importFrom GenomicRanges seqlengths "seqlengths<-" seqinfo
#' @importMethodsFrom ggbio autoplot layout_karyogram
#' @importFrom ggplot2 aes labs scale_colour_manual scale_fill_manual geom_text rel geom_segment xlab theme element_text
#' 
#' @examples
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
#' ## Overview with type pval
#' plotOverview(regions=regsWithP$regions, type="pval")
#'
#' ## Overview with type qval
#' plotOverview(regions=regsWithP$regions, type="qval")
#'
#' \dontrun{
#' ## Annotate the results
#' suppressMessages(library("bumphunter"))
#' annotation <- annotateNearest(regsWithP$regions, "hg19")
#'
#' ## Annotation overview
#' plotOverview(regions=regsWithP$regions, annotation=annotation, type="annotation")
#'
#' ## This function is meant to be an example of the plots you can make with the output from calculatePvalues()
#' ## For more details or custom plots check the ggbio and ggplot2 packages
#' ## as well as the code from this function:
#' plotOverview
#' }

plotOverview <- function(regions, annotation=NULL, type="pval", base_size=12, areaRel=4) {
	stopifnot(type %in% c("pval", "qval", "annotation"))
	
	## Keeping R CMD check happy
	hg19Ideogram <- significant <- midpoint <- area <- x <- y <- xend <- significantQval <- region <- NULL
	
	## Assign chr lengths if needed
	if(any(is.na(seqlengths(regions)))) {
		message(paste(Sys.time(), "plotOverview: assigning chromosome lengths from hg19!!!"))
		data(hg19Ideogram, package = "biovizBase", envir = environment())
		seqlengths(regions) <- seqlengths(hg19Ideogram)[names(seqlengths(regions))]
	}
		
	## Graphical setup
	ann_text <- data.frame(x=225e6, y=10, lab="Area", seqnames="chrX")
	ann_line <- data.frame(x=200e6, xend=215e6, y=10, seqnames="chrX")

	## Make the plot
	if(type == "pval") {
		## P-value plot
		result <- autoplot(seqinfo(regions)) +
			layout_karyogram(regions, aes(fill=significant, color=significant), geom="rect", base_size=30) +
			layout_karyogram(regions, aes(x=midpoint, y=area), geom="line", color="coral1", ylim=c(10, 20)) +
			labs(title="Overview of regions found in the genome; significant: p-value <0.05") +
			scale_colour_manual(values=c("chartreuse4", "wheat2"), limits=c("TRUE", "FALSE")) +
			scale_fill_manual(values=c("chartreuse4", "wheat2"), limits=c("TRUE", "FALSE")) +
			geom_text(aes(x=x, y=y), data=ann_text, label="Area", size=rel(areaRel)) +
			geom_segment(aes(x=x, xend=xend, y=y, yend=y), data=ann_line, colour="coral1") +
			xlab("Genomic coordinate") + 
			theme(text=element_text(size=base_size))
	} else if (type == "qval") {
		## Adjusted p-value plot
		result <- autoplot(seqinfo(regions)) +
			layout_karyogram(regions, aes(fill=significantQval, color=significantQval), geom="rect") +
			layout_karyogram(regions, aes(x=midpoint, y=area), geom="line", color="coral1", ylim=c(10, 20)) +
			labs(title="Overview of regions found in the genome; significant: q-value <0.10") +
			scale_colour_manual(values=c("chartreuse4", "wheat2"), limits=c("TRUE", "FALSE")) +
			scale_fill_manual(values=c("chartreuse4", "wheat2"), limits=c("TRUE", "FALSE")) +
			geom_text(aes(x=x, y=y), data=ann_text, label="Area", size=rel(areaRel)) +
			geom_segment(aes(x=x, xend=xend, y=y, yend=y), data=ann_line, colour="coral1") +
			xlab("Genomic coordinate") + 
			theme(text=element_text(size=base_size))
	} else {
		## Annotation region plot
		stopifnot(is.null(annotation) == FALSE)
		regions$region <- annotation$region
		result <- autoplot(seqinfo(regions)) +
			layout_karyogram(regions, aes(fill=region, color=region), geom="rect") +
			labs(title="Annotation region (if available)") +
			xlab("Genomic location") + 
			theme(text=element_text(size=base_size))
	}
	return(result)
}
