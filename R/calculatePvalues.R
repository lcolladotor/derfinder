#' Calculate p-values and identify regions
#'
#' First, this function finds the regions of interest according to specified 
#' cutoffs. Then it permutes the samples and re-calculates the F-statistics. 
#' The area of the statistics from these segments are then used to calculate 
#' p-values for the original regions.
#' 
#' @param coveragePrep A list with \code{$coverageProcessed}, 
#' \code{mclapplyIndeex}, and \code{$position} normally generated using 
#' \link{preprocessCoverage}.
#' @param models A list with \code{$mod} and \code{$mod0} normally generated 
#' using \link{makeModels}.
#' @param fstats A numerical Rle with the F-statistics normally generated using 
#' \link{calculateStats}.
#' @param nPermute The number of permutations. Note that for a full chromosome, 
#' a small amount (10) of permutations is sufficient. If set to 0, no 
#' permutations are performed and thus no null regions are used, however, the 
#' \code{$regions} component is created.
#' @param seeds An integer vector of length \code{nPermute} specifying the 
#' seeds to be used for each permutation. If \code{NULL} no seeds are used.
#' @param chr A single element character vector specifying the chromosome name. 
#' This argument is passed to \link{findRegions}.
#' @param maxRegionGap This argument is passed to \link{findRegions}.
#' @param maxClusterGap This argument is passed to \link{findRegions}.
#' @param cutoff This argument is passed to \link{getSegmentsRle}.
#' @param mc.cores This argument is passed to \link[BiocParallel]{SnowParam} 
#' to define the number of \code{workers} used for running \link{fstats.apply}.
#' @param mc.outfile This argument is passed to \link[BiocParallel]{SnowParam} 
#' to specify the \code{outfile} for any output from the workers.
#' @param verbose If \code{TRUE} basic status updates will be printed along the 
#' way.
#' @param significantCut A vector of length two specifiying the cutoffs used to 
#' determine significance. The first element is used to determine significance 
#' for the p-values and the second element is used for the q-values.
#' @param adjustF A single value to adjust that is added in the denominator of 
#' the F-stat calculation. Useful when the Residual Sum of Squares of the 
#' alternative model is very small.
#' @param lowMemDir The directory where the processed chunks are saved when 
#' using \link{preprocessCoverage} with a specified \code{lowMemDir}.
#' @param method This argument is passed to \link{fstats.apply}. Check the 
#' details there for more information
#' @param scalefac This argument is passed to \link{fstats.apply} and should be 
#' the same as the one used in \link{preprocessCoverage}.
#' @param chrsStyle The naming style of the chromosomes. By default, UCSC. See 
#' \link[GenomeInfoDb]{seqlevelsStyle}.
#'
#' @return A list with four components:
#' \describe{
#' \item{regions }{ is a GRanges with metadata columns given by 
#' \link{findRegions} with the additional metadata column \code{pvalues}: 
#' p-value of the region calculated via permutations of the samples; 
#' \code{padj}: the qvalues calculated using \link[qvalue]{qvalue}; 
#' \code{significant}: whether the p-value is less than 0.05 (by default); 
#' \code{significantPadj}: whether the q-value is less than 0.10 (by default). 
#' It also includes the mean coverage of the region (mean from the mean 
#' coverage at each base calculated in \link{preprocessCoverage}). Furthermore,
#' if \code{groupInfo} was not \code{NULL} in \link{preprocessCoverage}, then 
#' the group mean coverage is calculated as well as the log 2 fold change 
#' (using group 1 as the reference). }
#' \item{nullStats}{ is a numeric Rle with the mean of the null statistics by 
#' segment.}
#' \item{nullWidths}{ is a numeric Rle with the length of each of the segments 
#' in the null distribution. The area can be obtained by multiplying the 
#' absolute \code{nullstats} by the corresponding lengths.}
#' \item{nullPermutation}{ is a Rle with the permutation number from which the 
#' null region originated from.}
#' }
#'
#' @author Leonardo Collado-Torres
#' @seealso \link{findRegions}, \link{clusterMakerRle}, \link{getSegmentsRle}, 
#' \link{fstats.apply}, \link[qvalue]{qvalue}
#' @export
#' @importMethodsFrom IRanges quantile nrow ncol c mean lapply unlist 
#' as.numeric '$' '$<-' cbind
#' @importFrom IRanges Views RleList Rle IRanges Views DataFrame values 
#' 'values<-' nrow
#' @importFrom BiocParallel SnowParam SerialParam bplapply
#' @importFrom qvalue qvalue
#'
#' @examples
#' ## Collapse the coverage information
#' collapsedFull <- collapseFullCoverage(list(genomeData$coverage), 
#'     verbose=TRUE)
#' 
#' ## Calculate library size adjustments
#' sampleDepths <- sampleDepth(collapsedFull, probs=c(0.5), nonzero=TRUE, 
#'     verbose=TRUE)
#' 
#' ## Build the models
#' group <- genomeInfo$pop
#' adjustvars <- data.frame(genomeInfo$gender)
#' models <- makeModels(sampleDepths, testvars=group, adjustvars=adjustvars)
#'
#' ## Preprocess the data
#' ## Automatic chunksize used to then compare 1 vs 4 cores in the 'do not run'
#' ## section
#' prep <- preprocessCoverage(genomeData, groupInfo=group, cutoff=0, 
#'     scalefac=32, chunksize=NULL, colsubset=NULL, mc.cores=4)
#' 
#' ## Get the F statistics
#' fstats <- genomeFstats
#'
#' ## Determine a cutoff from the F-distribution.
#' ## This step is very important and you should consider using quantiles from 
#' ## the observed F statistics
#' n <- dim(prep$coverageProcessed)[2]
#' df1 <- dim(models$mod)[2]
#' df0 <- dim(models$mod0)[2]
#' cutoff <- qf(0.95, df1-df0, n-df1)
#' 
#' ## Low cutoff used for illustrative purposes
#' cutoff <- 1
#'
#' ## Calculate the p-values and define the regions of interest.
#' regsWithP <- calculatePvalues(prep, models, fstats, nPermute=1, seeds=1, 
#'     chr='chr21', cutoff=cutoff, mc.cores=1, method='regular')
#' regsWithP
#'
#' \dontrun{
#' ## Calculate again, but with 10 permutations instead of just 1
#' regsWithP <- calculatePvalues(prep, models, fstats, nPermute=10, seeds=1:10, 
#'     chr='chr21', cutoff=cutoff, mc.cores=2, method='regular')
#' 
#' ## Check that they are the same as the previously calculated regions
#' identical(regsWithP, genomeRegions)
#'
#' ## Histogram of the theoretical p-values by region
#' hist(pf(regsWithP$regions$value, df1-df0, n-df1), main='Distribution 
#'     original p-values by region', freq=FALSE)
#'
#' ## Histogram of the permutted p-values by region
#' hist(regsWithP$regions$pvalues, main='Distribution permutted p-values by 
#'     region', freq=FALSE)
#'
#' ## MA style plot
#' library('ggplot2')
#' ma <- data.frame(mean=regsWithP$regions$meanCoverage, 
#'     log2FoldChange=regsWithP$regions$log2FoldChangeYRIvsCEU)
#' ggplot(ma, aes(x=log2(mean), y=log2FoldChange)) + geom_point() + 
#'     ylab('Fold Change (log2)') + xlab('Mean coverage (log2)') + 
#'     labs(title='MA style plot')
#'
#' ## Annotate the results
#' library('bumphunter')
#' annotation <- annotateNearest(regsWithP$regions, 'hg19')
#' head(annotation)
#'
#' ## Compare speed between 1 and 4 cores (must have them!)
#' library('microbenchmark')
#' micro <- microbenchmark(
#'     calculatePvalues(prep, models, fstats, nPermute=10, seeds=NULL, 
#'         chr='chr21', cutoff=c(2, 5), mc.cores=1, verbose=FALSE,
#'         method='regular'),
#'     calculatePvalues(prep, models, fstats, nPermute=10, seeds=NULL, 
#'         chr='chr21', cutoff=c(2, 5), mc.cores=4, verbose=FALSE,
#'         method='regular'),
#'     times=10)
#' levels(micro$expr) <- c('one', 'four')
#' micro
#' ## Using 4 cores doesn't help with this toy data, but it will (at the 
#' ## expense of more RAM) if you have a larger data set.
#' }

calculatePvalues <- function(coveragePrep, models, fstats, nPermute = 1L, 
    seeds = as.integer(gsub("-", "", Sys.Date())) + seq_len(nPermute), 
    chr, maxRegionGap = 0L, maxClusterGap = 300L, cutoff = quantile(fstats, 
        0.99), mc.cores = getOption("mc.cores", 2L),
    mc.outfile = Sys.getenv('SGE_STDERR_PATH'), verbose = TRUE, 
    significantCut = c(0.05, 0.1), adjustF = 0, lowMemDir = NULL,
    method = 'Matrix', scalefac = 32, chrsStyle = "UCSC") {
    ## Setup
    if (is.null(seeds)) {
        seeds <- rep(NA, nPermute)
    }
    
    ## Run some checks
    stopifnot(nPermute == length(seeds))
    stopifnot(length(intersect(names(coveragePrep), c("coverageProcessed", 
        "mclapplyIndex", "position", "meanCoverage", "groupMeans"))) == 
        5)
    stopifnot(length(intersect(names(models), c("mod", "mod0"))) == 
        2)
    stopifnot(length(significantCut) == 2 & all(significantCut >= 
        0 & significantCut <= 1))
    
    ## Identify the data segments
    if (verbose) 
        message(paste(Sys.time(),
            "calculatePvalues: identifying data segments"))
    
    ## Extract data
    position <- coveragePrep$position
    means <- coveragePrep$meanCoverage
    groupMeans <- coveragePrep$groupMeans
    mclapplyIndex <- coveragePrep$mclapplyIndex
    coverageProcessed <- coveragePrep$coverageProcessed
    if (is.null(lowMemDir) & is.null(coverageProcessed)) 
        stop("preprocessCoverage() was used with a non-null 'lowMemDir', so please specify 'lowMemDir'.")
    rm(coveragePrep)
    
    ## Avoid re-calculating possible candidate DERs for every
    ## permutation
    segmentIR <- clusterMakerRle(position, maxRegionGap, ranges = TRUE)
    
    ## Find the regions
    regs <- findRegions(position = position, chr = chr,
        maxRegionGap = maxRegionGap, maxClusterGap = maxClusterGap,
        fstats = fstats, cutoff = cutoff, segmentIR = segmentIR,
        chrsStyle = chrsStyle, verbose = verbose)
    if (is.null(regs)) {
        final <- list(regions = NULL, nullStats = NULL, nullWidths = NULL, 
            nullPermutation = NULL)
        return(final)
    }
    
    ## Assign mean coverage (overall)
    indexIR <- IRanges(start = regs$indexStart, end = regs$indexEnd)
    regs$meanCoverage <- mean(Views(means, indexIR))
    
    ## Calculate mean coverage by group and fold changes
    if (length(groupMeans) > 0) {
        regionGroupMean <- lapply(groupMeans, function(x) {
            mean(Views(x, indexIR))
        })
        
        ## Calculate fold coverage vs group 1
        if (length(regionGroupMean) > 1) {
            log2FoldChange <- vector("list", length(regionGroupMean) - 
                1)
            names(log2FoldChange) <- names(regionGroupMean)[-1]
            for (group in names(log2FoldChange)) {
                log2FoldChange[[group]] <- 
                    log2(regionGroupMean[[group]]/regionGroupMean[[1]])
            }
        }
        
        ## Finish up
        names(regionGroupMean) <- paste0("mean", names(regionGroupMean))
        values(regs) <- cbind(values(regs), DataFrame(regionGroupMean))
        if (length(regionGroupMean) > 1) {
            names(log2FoldChange) <- paste0("log2FoldChange", 
                names(log2FoldChange), "vs", names(groupMeans)[1])
            values(regs) <- cbind(values(regs), DataFrame(log2FoldChange))
            rm(log2FoldChange)
        }
        rm(regionGroupMean)
    }
    
    rm(fstats, position, means, groupMeans)
    
    ## Pre-allocate memory
    nullareas <- nullpermutation <- nullwidths <- nullstats <- vector("list", 
        length(seeds) * 2)
    last <- 0
    nSamples <- seq_len(nrow(models$mod))
    
    for (i in seq_along(seeds)) {
        if (verbose) 
            message(paste(Sys.time(),
                "calculatePvalues: calculating F-statistics for permutation", 
                i, "and seed", seeds[i]))
        
        if (!is.na(seeds[i])) {
            set.seed(seeds[i])
        }
        idx.permute <- sample(nSamples)
        
        ## Permuted sample labels
        mod.p <- models$mod[idx.permute, , drop = FALSE]
        mod0.p <- models$mod0[idx.permute, , drop = FALSE]
        
        ## Define cluster
        if(mc.cores > 1) {
            BPPARAM <- SnowParam(workers = mc.cores, outfile = mc.outfile)
        } else {
            BPPARAM <- SerialParam()
        }
        
        
        ## Get the F-statistics
        fstats.output <- bplapply(mclapplyIndex, fstats.apply, 
            data = coverageProcessed, mod = mod.p, mod0 = mod0.p, 
            adjustF = adjustF, lowMemDir = lowMemDir, method = method, 
            scalefac = scalefac, BPPARAM = BPPARAM)
        fstats.output <- unlist(RleList(fstats.output), use.names = FALSE)
        
        ## Find the segments
        regs.perm <- findRegions(chr = chr, maxRegionGap = maxRegionGap, 
            maxClusterGap = maxClusterGap, fstats = fstats.output, 
            cutoff = cutoff, segmentIR = segmentIR, basic = TRUE, 
            verbose = verbose)
        
        ## Calculate mean statistics
        if (!is.null(regs.perm)) {
            for (j in 1:2) {
                nullstats[[last + j]] <- regs.perm$stat
                nullwidths[[last + j]] <- regs.perm$width
                nullareas[[last + j]] <- regs.perm$area
                nullpermutation[[last + j]] <- Rle(i, nrow(regs.perm))
            }
        }
        last <- last + 2
        
        ## Finish loop
        rm(idx.permute, fstats.output, regs.perm, mod.p, mod0.p)        
    }
    nullstats <- do.call(c, nullstats[!sapply(nullstats, is.null)])
    nullwidths <- do.call(c, nullwidths[!sapply(nullwidths, is.null)])
    nullpermutation <- do.call(c, nullpermutation[!sapply(nullpermutation, 
        is.null)])
    nullareas <- do.call(c, nullareas[!sapply(nullareas, is.null)])
    
    if (length(nullstats) > 0) {
        ## Proceed only if there is at least one null stats
        
        ## Calculate pvalues
        if (verbose) 
            message(paste(Sys.time(),
                "calculatePvalues: calculating the p-values"))
        pvals <- sapply(regs$area, function(x) {
            sum(nullareas > x)
        })
        regs$pvalues <- (pvals + 1)/(length(nullareas) + 1)
        regs$significant <- factor(regs$pvalues < significantCut[1], 
            levels = c(TRUE, FALSE))
        
        ## Sometimes qvalue() fails due to incorrect pi0 estimates
        qvalues <- qvalue(regs$pvalues)
        if (is(qvalues, "qvalue")) {
            qvalues <- qvalues$qvalues
            sigQval <- factor(qvalues < significantCut[2], levels = c(TRUE, 
                FALSE))
        } else {
            qvalues <- rep(NA, length(regs$pvalues))
            sigQval <- rep(NA, length(regs$pvalues))
        }
        regs$qvalues <- qvalues
        regs$significantQval <- sigQval
        
        regs <- regs[order(regs$area, decreasing = TRUE), ]
    } else {
        if (verbose) 
            message(paste(Sys.time(),
                "calculatePvalues: no null regions found. Skipping p-value calculation."))
        regs$pvalues <- rep(NA, length(regs))
        regs$significant <- rep(NA, length(regs))
        regs$qvalues <- rep(NA, length(regs))
        regs$significantQval <- rep(NA, length(regs))
    }
    ## Save the nullstats too
    final <- list(regions = regs, nullStats = nullstats,
            nullWidths = nullwidths, nullPermutation = nullpermutation)
    
    
    ## Done =)
    return(final)
} 
