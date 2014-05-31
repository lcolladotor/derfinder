#' Plot a karyotype overview of the genome with the identified regions
#'
#' Plots an overview of the genomic locations of the identified regions (see 
#' \link{calculatePvalues}) in a karyotype view. The coloring can be done 
#' either by significant regions according to their p-values, significant by 
#' adjusted p-values, or by annotated region if using 
#' \link[bumphunter]{annotateNearest}.
#'
#' 
#' @param regions The \code{$regions} output from \link{calculatePvalues}.
#' @param annotation The output from running \link[bumphunter]{annotateNearest} 
#' on the output from \link{calculatePvalues}. It is only required if 
#' \code{type='annotation'}.
#' @param type Must be either \code{pval}, \code{qval} or \code{annotation}. It 
#' determines whether the plot coloring should be done according to significant 
#' p-values (<0.05), significant q-values (<0.10) or annotation regions.
#' @param base_size Base point size of the plot. This argument is passed to 
#' \link[ggplot2]{element_text} (\code{size} argument).
#' @param areaRel The relative size for the area label when \code{type='pval'} 
#' or \code{type='qval'}. Can be useful when making high resolution versions of 
#' these plots in devices like CairoPNG.
#' @param legend.position This argument is passed to \link[ggplot2]{theme}. 
#' From ggplot2: the position of legends. ('left', 'right', 'bottom', 'top', or 
#' two-element numeric vector).
#' @param significantCut A vector of length two specifiying the cutoffs used to 
#' determine significance. The first element is used to determine significance 
#' for the p-values and the second element is used for the q-values.
#' @param chrsStyle The naming style of the chromosomes. By default, UCSC. See 
#' \link[GenomeInfoDb]{seqlevelsStyle}.
#'
#' @return A ggplot2 plot that is ready to be printed out. Tecnically it is a 
#' ggbio object.
#'
#' @seealso \link{calculatePvalues}, \link[bumphunter]{annotateNearest}
#' @author Leonardo Collado-Torres
#' @export
#' @importFrom GenomicRanges seqinfo
#' @importFrom GenomeInfoDb seqlengths 'seqlengths<-' seqlevelsStyle 
#' 'seqlevelsStyle<-'
#' @importMethodsFrom ggbio autoplot layout_karyogram
#' @importFrom ggplot2 aes labs scale_colour_manual scale_fill_manual geom_text 
#' rel geom_segment xlab theme element_text element_blank
#' 
#' @examples
#' ## Construct toy data
#' chrs <- paste0('chr', c(1:22, 'X', 'Y'))
#' chrs <- factor(chrs, levels=chrs)
#' library('GenomicRanges')
#' regs <- GRanges(rep(chrs, 10), ranges=IRanges(runif(240, 1, 4e7), 
#'     width=1e3), significant=sample(c(TRUE, FALSE), 240, TRUE, p=c(0.05, 
#'     0.95)), significantQval=sample(c(TRUE, FALSE), 240, TRUE, p=c(0.1, 
#'     0.9)), area=rnorm(240))
#' annotation <- data.frame(region=sample(c('upstream', 'promoter', 
#'     "overlaps 5'", 'inside', "overlaps 3'", "close to 3'", 'downstream'), 
#'     240, TRUE))
#'
#' ## Type pval
#' plotOverview(regs)
#'
#' \dontrun{
#' ## Type qval
#' plotOverview(regs, type='qval')
#'
#' ## Annotation
#' plotOverview(regs, annotation, type='annotation')
#' 
#' ## Resize the plots if needed.
#'
#' ## You might prefer to leave the legend at ggplot2's default option: right
#' plotOverview(regs, legend.position='right')
#' 
#' ## Although the legend looks better on the bottom
#' plotOverview(regs, legend.position='bottom')
#'
#' ## Example knitr chunk for higher res plot using the CairoPNG device
#' ```{r overview, message=FALSE, fig.width=7, fig.height=9, dev='CairoPNG', dpi=300}
#' plotOverview(regs, base_size=30, areaRel=10, legend.position=c(0.95, 0.12))
#' ```
#' 
#' ## For more custom plots, take a look at the ggplot2 and ggbio packages
#' ## and feel free to look at the code of this function:
#' plotOverview
#' }

plotOverview <- function(regions, annotation = NULL, type = "pval", 
    base_size = 12, areaRel = 4, legend.position = c(0.85, 0.12), 
    significantCut = c(0.05, 0.1), chrsStyle = "UCSC") {
    stopifnot(type %in% c("pval", "qval", "annotation"))
    stopifnot(length(significantCut) == 2 & all(significantCut >= 
        0 & significantCut <= 1))
    
    ## Keeping R CMD check happy
    hg19Ideogram <- significant <- midpoint <- area <- x <- y <- xend <-
        significantQval <- region <- NULL
        
    ## Use UCSC names by default
    seqlevelsStyle(regions) <- chrsStyle
    
    ## Assign chr lengths if needed
    if (any(is.na(seqlengths(regions)))) {
        message(paste(Sys.time(), 
            "plotOverview: assigning chromosome lengths from hg19!!!"))
        data(hg19Ideogram, package = "biovizBase", envir = environment())
        seqlengths(regions) <- seqlengths(hg19Ideogram)[
            names(seqlengths(regions))]
    }
    
    ## Graphical setup
    ann_chr <- ifelse(any(seqnames(regs) == "chrX"), "chrX",
        levels(seqnames(regs))[length(levels(seqnames(regs)))])
    ann_text <- data.frame(x = 2.25e+08, y = 10, lab = "Area", 
        seqnames = ann_chr)
    ann_line <- data.frame(x = 2e+08, xend = 2.15e+08, y = 10, 
        seqnames = ann_chr)
    
    ## Make the plot
    if (type == "pval") {
        ## P-value plot
        result <- autoplot(seqinfo(regions)) + layout_karyogram(regions, 
            aes(fill = significant, color = significant), geom = "rect", 
            base_size = 30) + layout_karyogram(regions, aes(x = midpoint, 
            y = area), geom = "line", color = "coral1", ylim = c(10, 
            20)) + labs(title = paste0("Overview of regions found in the genome; significant: p-value <", 
            significantCut[1])) + scale_colour_manual(values = c("chartreuse4", 
            "wheat2"), limits = c("TRUE", "FALSE")) + scale_fill_manual(
                values = c("chartreuse4", "wheat2"), limits = c("TRUE", 
                "FALSE")) + geom_text(aes(x = x, y = y), data = ann_text,
                label = "Area", size = rel(areaRel)) + 
            geom_segment(aes(x = x, xend = xend, y = y, yend = y), 
                data = ann_line, colour = "coral1") + 
                xlab("Genomic coordinate") + theme(text = element_text(
                    size = base_size), legend.background = element_blank(), 
                legend.position = legend.position)
    } else if (type == "qval") {
        ## Adjusted p-value plot
        result <- autoplot(seqinfo(regions)) + layout_karyogram(regions, 
            aes(fill = significantQval, color = significantQval), 
            geom = "rect") + layout_karyogram(regions, aes(x = midpoint, 
            y = area), geom = "line", color = "coral1", ylim = c(10, 
            20)) + labs(title = paste0("Overview of regions found in the genome; significant: q-value <", 
            significantCut[2])) + scale_colour_manual(values = c("chartreuse4", 
            "wheat2"), limits = c("TRUE", "FALSE")) + scale_fill_manual(
                values = c("chartreuse4", "wheat2"), limits = c("TRUE",
                "FALSE")) + geom_text(aes(x = x, y = y), data = ann_text,
                label = "Area", size = rel(areaRel)) + 
            geom_segment(aes(x = x, xend = xend, y = y, yend = y), 
                data = ann_line, colour = "coral1") +
                xlab("Genomic coordinate") + theme(text = element_text(
                    size = base_size), legend.background = element_blank(), 
                legend.position = legend.position)
    } else {
        ## Annotation region plot
        stopifnot(is.null(annotation) == FALSE)
        regions$region <- annotation$region
        result <- autoplot(seqinfo(regions)) + layout_karyogram(regions, 
            aes(fill = region, color = region), geom = "rect") + 
            labs(title = "Annotation region (if available)") + 
            xlab("Genomic location") + theme(text = element_text(
                size = base_size), legend.background = element_blank(),
                legend.position = legend.position)
    }
    return(result)
} 
