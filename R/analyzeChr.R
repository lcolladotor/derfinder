#' Run the derfinder analysis on a chromosome
#'
#' This is a major wrapper for running several key functions from this package. 
#' It is meant to be used after \link{loadCoverage} has been used for a 
#' specific chromosome. The steps run include \link{makeModels}, 
#' \link{preprocessCoverage}, \link{calculateStats}, \link{calculatePvalues} 
#' and \link[bumphunter]{annotateNearest}. 
#' 
#' @param chr Used for naming the output files when \code{writeOutput=TRUE} 
#' and for \link[bumphunter]{annotateNearest}.
#' @param models The output from \link{makeModels}.
#' @param cutoffPre This argument is passed to \link{preprocessCoverage} 
#' (\code{cutoff}).
#' @inheritParams preprocessCoverage
#' @param cutoffFstat This is used to determine the cutoff argument of 
#' \link{calculatePvalues} and it's behaviour is determined by 
#' \code{cutoffType}.
#' @param cutoffType If set to \code{empirical}, the \code{cutoffFstat} 
#' (example: 0.99) quantile is used via \link{quantile}. If set to 
#' \code{theoretical}, the theoretical \code{cutoffFstats} (example: 1e-08) is 
#' calculated via \link{qf}. If set to \code{manual}, \code{cutoffFstats} is 
#' passed to \link{calculatePvalues} without any other calculation.
#' @inheritParams calculatePvalues
#' @param groupInfo A factor specifying the group membership of each sample 
#' that can later be used with the plotting functions in the 
#' \code{derfinderPlot} package.
#' @param subject This argument is passed to 
#' \link[bumphunter]{annotateNearest}. Note that only \code{hg19} works right 
#' now.
#' @param writeOutput If \code{TRUE}, output Rdata files are created at each 
#' step inside a directory with the chromosome name (example: 'chr21' if 
#' \code{chrnum='21'}). One Rdata files is created for each component described 
#' in the return section.
#' @param runAnnotation If \code{TRUE} \link[bumphunter]{annotateNearest} is 
#' run. Otherwise this step is skipped.
#' @param ... Arguments passed to other methods and/or advanced arguments.
#'
#' @return If \code{returnOutput=TRUE}, a list with six components:
#' \describe{
#' \item{timeinfo }{ The wallclock timing information for each step.}
#' \item{optionsStats }{ The main options used when running this function.}
#' \item{coveragePrep }{ The output from \link{preprocessCoverage}.}
#' \item{fstats}{ The output from \link{calculateStats}.}
#' \item{regions}{ The output from \link{calculatePvalues}.}
#' \item{annotation}{ The output from \link[bumphunter]{annotateNearest}.}
#' }
#' These are the same components that are written to Rdata files if 
#' \code{writeOutput=TRUE}.
#'
#' @author Leonardo Collado-Torres
#' @seealso \link{makeModels}, \link{preprocessCoverage}, 
#' \link{calculateStats}, \link{calculatePvalues}, 
#' \link[bumphunter]{annotateNearest}
#' @export
#' @aliases analyze_chr
#' @importMethodsFrom S4Vectors as.numeric
#' @importFrom bumphunter annotateNearest
#' @importFrom GenomeInfoDb mapSeqlevels
#' 
#' @examples
#' ## Collapse the coverage information
#' collapsedFull <- collapseFullCoverage(list(genomeData$coverage), 
#'     verbose = TRUE)
#' 
#' ## Calculate library size adjustments
#' sampleDepths <- sampleDepth(collapsedFull, probs = c(0.5), nonzero=TRUE, 
#'     verbose=TRUE)
#' 
#' ## Build the models
#' group <- genomeInfo$pop
#' adjustvars <- data.frame(genomeInfo$gender)
#' models <- makeModels(sampleDepths, testvars=group, adjustvars=adjustvars)
#'
#' ## Analyze the chromosome
#' results <- analyzeChr(chr='21', coverageInfo=genomeData, models=models, 
#'     cutoffFstat=1, cutoffType='manual', groupInfo=group, mc.cores=1, 
#'     writeOutput=FALSE, returnOutput=TRUE, method='regular')
#' names(results)

analyzeChr <- function(chr, coverageInfo, models, cutoffPre = 5, 
    cutoffFstat = 1e-08, cutoffType = 'theoretical', nPermute = 1, 
    seeds = as.integer(gsub('-', '', Sys.Date())) + seq_len(nPermute), 
    groupInfo, subject = 'hg19', writeOutput = TRUE, runAnnotation = TRUE, ...){
        
    ## Run some checks
    stopifnot(length(intersect(cutoffType, c('empirical', 'theoretical', 
        'manual'))) == 1)
    stopifnot(is.factor(groupInfo))
    
    ## Advanged argumentsa
#' @param chrsStyle The naming style of the chromosomes. By default, UCSC. See 
#' \link[GenomeInfoDb]{seqlevelsStyle}.    
    chrsStyle <- .advanced_argument('chrsStyle', 'UCSC', ...)


    ## Use UCSC names by default
    chr <- mapSeqlevels(chr, chrsStyle)

#' @param verbose If \code{TRUE} basic status updates will be printed along the 
#' way.
    verbose <- .advanced_argument('verbose', TRUE, ...)



#' @param scalefac This argument is passed to \link{preprocessCoverage}.
    scalefac <- .advanced_argument('scalefac', 32, ...)



#' @param chunksize This argument is passed to \link{preprocessCoverage}.
    chunksize <- .advanced_argument('chunksize', NULL, ...)
    

#' @param lowMemDir The directory where the processed chunks are saved when
#' using \link{preprocessCoverage} with a specified \code{lowMemDir}.
    lowMemDir <- .advanced_argument('lowMemDir', file.path(chr, 'chunksDir'),
        ...)
        

#' @param returnOutput If \code{TRUE}, it returns a list with the results from 
#' each step. Otherwise, it returns \code{NULL}.
    returnOutput <- .advanced_argument('returnOutput', !writeOutput, ...)
        

    ## Begin timing
    timeinfo <- NULL
    ## Init
    timeinfo <- c(timeinfo, list(Sys.time()))
    
    ## Drop unused levels in groupInfo
    groupInfo <- droplevels(groupInfo)
    
    ## Setup
    timeinfo <- c(timeinfo, list(Sys.time()))
    
    ## pre-process the coverage data with automatic chunks
    ## depending on the number of cores
    if (verbose) 
        message(paste(Sys.time(),
            'analyzeChr: Pre-processing the coverage data'))
    prep <- preprocessCoverage(coverageInfo = coverageInfo,
        groupInfo = groupInfo, cutoff = cutoffPre, ...)
    rm(coverageInfo)
    
    ## prepData
    timeinfo <- c(timeinfo, list(Sys.time()))
    
    ## Save the prepared data
    if (writeOutput) {
        save(prep, file = file.path(chr, 'coveragePrep.Rdata'))
    }
    ## savePrep
    timeinfo <- c(timeinfo, list(Sys.time()))
    
    ## Run calculateStats
    if (verbose) 
        message(paste(Sys.time(), 'analyzeChr: Calculating statistics'))
    fstats <- calculateStats(coveragePrep = prep, models = models, ...)
    
    ## calculateStats
    timeinfo <- c(timeinfo, list(Sys.time()))
    
    ## Save the output from calculateStats
    if (writeOutput) {
        save(fstats, file = file.path(chr, 'fstats.Rdata'))
    }
    
    ## saveStats
    timeinfo <- c(timeinfo, list(Sys.time()))
    
    ## Choose the cutoff
    cutoff <- .calcFstatCutoff(cutoffType, cutoffFstat, fstats, models)
    
    ## Save parameters used for running calculateStats
    optionsStats <- list(models = models, cutoffPre = cutoffPre, 
        scalefac = scalefac, chunksize = chunksize, 
        cutoffFstat = cutoffFstat, cutoffType = cutoffType, 
        nPermute = nPermute, seeds = seeds, groupInfo = groupInfo,
        lowMemDir = lowMemDir, analyzeCall = match.call(), 
        cutoffFstatUsed = cutoff, ...)
        
    if (writeOutput) {
        dir.create(chr, showWarnings = FALSE, recursive = TRUE)
        save(optionsStats, file = file.path(chr, 'optionsStats.Rdata'))
    }
    
    ## saveStatsOpts
    timeinfo <- c(timeinfo, list(Sys.time()))
    
    ## Calculate p-values and find regions
    if (verbose) 
        message(paste(Sys.time(), 'analyzeChr: Calculating pvalues'))
    
    if (verbose) 
        message(paste(Sys.time(), 'analyzeChr: Using the following', 
            cutoffType, 'cutoff for the F-statistics', cutoff))
    
    regions <- calculatePvalues(coveragePrep = prep, models = models, 
        fstats = fstats, nPermute = nPermute, seeds = seeds, 
        chr = chr, cutoff = cutoff, ...)
    if (!returnOutput) {
        rm(prep)
    }
    
    ## calculatePValues
    timeinfo <- c(timeinfo, list(Sys.time()))
    
    ## Save the output from calculatePvalues
    if (writeOutput) {
        save(regions, file = file.path(chr, 'regions.Rdata'))
    }
    
    ## saveRegs
    timeinfo <- c(timeinfo, list(Sys.time()))
    
    ## Annotate
    if (verbose) 
        message(paste(Sys.time(), 'analyzeChr: Annotating regions'))
    
    if (!is.null(regions$regions) & runAnnotation) {
        annotation <- annotateNearest(regions$regions, subject)
    } else {
        annotation <- NULL
    }
    
    ## Annotate
    timeinfo <- c(timeinfo, list(Sys.time()))
    
    if (writeOutput) {
        save(annotation, file = file.path(chr, 'annotation.Rdata'))
    }
    
    ## saveAnnotation
    timeinfo <- c(timeinfo, list(Sys.time()))
    
    ## Save timing information
    timeinfo <- do.call(c, timeinfo)
    names(timeinfo) <- c('init', 'setup', 'prepData', 'savePrep',
        'calculateStats', 'saveStats', 'saveStatsOpts', 'calculatePvalues', 
        'saveRegs', 'annotate', 'saveAnno')
    if (writeOutput) {
        save(timeinfo, file = file.path(chr, 'timeinfo.Rdata'))
    }
    
    if (returnOutput) {
        result <- list(timeinfo = timeinfo, optionsStats = optionsStats, 
            coveragePrep = prep, fstats = fstats, regions = regions, 
            annotation = annotation)
    } else {
        result <- NULL
    }
    
    ## Done
    return(invisible(result))
} 

#' @export
analyze_chr <- analyzeChr

## Helper function for calculating the F-stat cutoff
.calcFstatCutoff <- function(cutoffType, cutoffFstat, fstats, models) {
    if (cutoffType == 'empirical') {
        if(cutoffFstat == 1e-08) {
            cutoffFstat <- 0.99
            warning("Switching 'cutoffFstat' to 0.99 as the user probably forgot to change its default value.")
        }
        cutoff <- quantile(as.numeric(fstats), cutoffFstat)
    } else if (cutoffType == 'theoretical') {
        n <- dim(models$mod)[1]
        df1 <- dim(models$mod)[2]
        df0 <- dim(models$mod0)[2]
        cutoff <- qf(cutoffFstat, df1 - df0, n - df1, lower.tail = FALSE)
    } else if (cutoffType == 'manual') {
        cutoff <- cutoffFstat
    }
    return(cutoff)
}
