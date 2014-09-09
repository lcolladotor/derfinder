#' Load the coverage information from a group of BAM files
#'
#' For a group of samples this function reads the coverage information for a 
#' specific chromosome directly from the BAM files. It then merges them into a 
#' DataFrame and removes the bases that do not pass the cutoff.
#' 
#' @param dirs A character vector with the full path to the sample BAM files
#' (or bigWig files). 
#' The names are used for the column names of the DataFrame. Check 
#' \link{rawFiles} for constructing \code{dirs}. \code{dirs} can also be a 
#' \code{BamFileList} object created with \link[Rsamtools]{BamFileList} or a
#' \code{BigWigFileList} object created with \link[rtracklayer]{BigWigFileList}.
#' @param chr Chromosome to read. Should be in the format matching the one used
#' in the raw data.
#' @param cutoff This argument is passed to \link{filterData}.
#' @param bai The full path to the BAM index files. If \code{NULL} it is 
#' assumed that the BAM index files are in the same location as the BAM files 
#' and that they have the .bai extension. Ignored if \code{dirs} is a 
#' \code{BamFileList} object or if \code{inputType=='bigWig'}.
#' @param chrlen The chromosome length in base pairs. If it's \code{NULL}, the 
#' chromosome length is extracted from the BAM files.
#' @param output If \code{NULL} then no output is saved in disk. If \code{auto} 
#' then an automatic name is constructed using UCSC names (chrXCovInfo.Rdata
#' for example). If another character is specified, then that name is used for #' the output file.
#' @param inputType Has to be either \code{bam} or \code{bigWig}. It specifies
#' the format of the raw data files.
#' @param isMinusStrand Use \code{TRUE} for negative strand alignments only, 
#' \code{FALSE} for positive strands and \code{NA} for both. This argument is 
#' passed to \link[Rsamtools]{scanBamFlag} when \code{inputType='bam'}.
#' @param filter This argument is passed to \link{filterData}.
#' @param returnMean This argument is passed to \link{filterData}.
#' @param returnCoverage This argument is passed to \link{filterData}.
#' @param totalMapped The total number of reads mapped for each sample. 
#' Providing this data adjusts the coverage to reads in \code{targetSize} 
#' library prior to filtering. By default, to reads per 80 million reads.
#' @param targetSize The target library size to adjust the coverage to. Used
#' only when \code{totalMapped} is specified.
#' @param tilewidth When specified, \link[GenomicRanges]{tileGenome} is used to
#' break up the chromosome into chunks.
#' @param mc.cores This argument is passed to \link[BiocParallel]{SnowParam} 
#' to define the number of \code{workers}. You may use up to one core per tile.
#' Only used when \code{tilewidth} is specified.
#' @param mc.outfile This argument is passed to \link[BiocParallel]{SnowParam} 
#' to specify the \code{outfile} for any output from the workers.
#' Only used when \code{tilewidth} is specified.
#' @param verbose If \code{TRUE} basic status updates will be printed along the 
#' way.
#'
#' @return A list with two components.
#' \describe{
#' \item{coverage }{ is a DataFrame object where each column represents a 
#' sample. The number of rows depends on the number of base pairs that passed 
#' the cutoff and the information stored is the coverage at that given base.}
#' \item{position }{  is a logical Rle with the positions of the chromosome 
#' that passed the cutoff.}
#' }
#'
#' @author Leonardo Collado-Torres, Andrew Jaffe
#' @export
#' @aliases load_coverage
#' @importFrom Rsamtools BamFileList scanBamHeader ScanBamParam path scanBamFlag
#' @importFrom GenomicAlignments readGAlignmentsFromBam
#' @importFrom IRanges IRanges RangesList
#' @importFrom rtracklayer BigWigFileList
#' @importFrom GenomeInfoDb mapSeqlevels
#' @importFrom GenomicRanges tileGenome
#' @importFrom GenomicFiles reduceByFile
#' @importFrom BiocParallel SnowParam SerialParam
#' @importMethodsFrom S4Vectors Reduce
#' @importMethodsFrom GenomicRanges coverage
#' @importMethodsFrom Rsamtools names
#' @importMethodsFrom rtracklayer import import.bw
#' @examples
#' datadir <- system.file('extdata', 'genomeData', package='derfinder')
#' dirs <- rawFiles(datadir=datadir, samplepatt='*accepted_hits.bam$', 
#'     fileterm=NULL)
#' ## Shorten the column names
#' names(dirs) <- gsub('_accepted_hits.bam', '', names(dirs))
#'  
#' ## Read and filter the data, only for 2 files
#' dataSmall <- loadCoverage(dirs=dirs[1:2], chr='21', cutoff=0)
#'
#' \dontrun{
#' ## Read all the data
#' dataAll <- loadCoverage(dirs=dirs, chr='21', cutoff=0)
#'
#' ## For other examples this data is included in the package
#' identical(dataAll, genomeData)
#'
#' ## Read the coverage without applying any cutoff.
#' ## This can be useful for downstream analysis including coverage plots.
#' dataRaw <- loadCoverage(dirs=dirs, chr='21', cutoff=NULL)
#'
#' ## Compare vs raw data provided in the package
#' identical(dataRaw, genomeDataRaw)
#' 
#' ## Note that the object size is pretty much the same due to the Rle
#' ## compression
#' print(object.size(dataRaw), units='Kb')
#' 
#'
#' #################################################################
#' 
#' ## The following code shows how to export the coverage to a BigWig file
#' library('rtracklayer')
#' sample1 <- RleList("chr21" = dataSmall$coverage[[1]])
#' export(sample1, "sample1.bw")
#' 
#' ## Re-load it in R
#' sample1.new <- import("sample1.bw", selection = GRanges("chr21", IRanges(1, 
#'    48129895)), as = "RleList")
#' 
#' ## Compare them
#' sample1 - sample1.new
#' ## Note that the original one is a RleList of integer-Rle while the new one 
#' ## is a RleList of numeric-Rle
#'
#' #################################################################
#' ## Below is an example of using loadCoverage(inputType='bigWig')
#' 
#' ## First we need to export the second sample too
#' sample2 <- RleList("chr21" = dataSmall$coverage[[2]])
#' export(sample2, "sample2.bw")
#'
#' ## Now we can use loadCoverage
#' dataSmall.bigWig <- loadCoverage(c(ERR009101 = 'sample1.bw', ERR009102 = 
#'     'sample2.bw'), chr="chr21", inputType="bigWig")
#' 
#' ## We can compare the results
#' mapply(function(x, y) { x - y}, dataSmall$coverage, 
#'     dataSmall.bigWig$coverage)
#' ## Note that the only difference is the type of Rle (integer vs numeric) used
#' ## to store the data.
#'
#' }

loadCoverage <- function(dirs, chr, cutoff = NULL, bai = NULL,
    chrlen = NULL, output = NULL, inputType = "bam", isMinusStrand = NA,
    filter = "one", returnMean = FALSE, returnCoverage = TRUE,
    totalMapped = NULL, targetSize = 80e6, tilewidth = NULL, mc.cores = 1,
    mc.outfile = Sys.getenv('SGE_STDERR_PATH'), verbose = TRUE) {
    stopifnot(inputType %in% c("bam", "bigWig"))
        
    ## Do the indexes exist?
    if (is(dirs, "BamFileList") & inputType == "bam") {
        bList <- dirs
    } else if (is(dirs, "BigWigFileList") & inputType == "bigWig") {
        bList <- dirs
    } else if (inputType == "bam"){
        if (is.null(bai)) {
            bai <- paste0(dirs, ".bai")
        }
        if (all(file.exists(bai))) {
            bList <- BamFileList(dirs, bai)
        } else {
            stop("Not all BAM files have a BAM index. If the BAM index files are in a separate directory from the BAM files or are not named as 'bamFile.bam.bai' then consider using the 'bai' argument.")
        }
    } else if (inputType == "bigWig") {
        bList <- BigWigFileList(dirs)
    }
    
    
    ## Determine the chromosome length
    if (is.null(chrlen)) {
        ## This assumes that all the BAM files are from the same
        ## organism.
        if (verbose) 
            message(paste(Sys.time(),
                "loadCoverage: finding chromosome lengths"))
        if (inputType == "bam") {
            clengths <- scanBamHeader(bList[[1]])$targets
        } else if (inputType == "bigWig") {
            clengths <- seqlengths(bList[[1]])
        }
       
        if (!chr %in% names(clengths)) {
            stop(paste("'chr' is not correctly specified. Valid options are:", 
                paste(names(clengths), collapse = ", ")))
        }
        chrlen <- clengths[chr]
    }
    
    ## Make tiles if using GenomicFiles
    if(!is.null(tilewidth)) {
        tiles <- tileGenome(chrlen, tilewidth = tilewidth)
        
        ## Define cluster
        if(mc.cores > 1) {
            BPPARAM <- SnowParam(workers = mc.cores, outfile = mc.outfile)
        } else {
            BPPARAM <- SerialParam()
        }
    }    
    
    ## Construct the objects so only the chr of interest is read
    ## from the BAM/bigWig file
    which <- GRanges(seqnames=chr, ranges=IRanges(1, chrlen))
    if(inputType == "bam") {
        param <- ScanBamParam(which = which, 
            flag = scanBamFlag(isMinusStrand = isMinusStrand))
        
        ## Read in the data for all the chrs
        if(is.null(tilewidth)) {
            data <- lapply(bList, .loadCoverageBAM, param = param, chr = chr,
                verbose=verbose)
        } else {
            data <- reduceByFile(tiles, bList, .bamMAPPER, .REDUCER, 
                chr = chr, verbose = verbose, isMinusStrand = isMinusStrand, 
                BPPARAM = BPPARAM)
        }
        
    } else if (inputType == "bigWig") {
        ## Read in the data for all the chrs
        if(is.null(tilewidth)) {
            data <- lapply(bList, .loadCoverageBigWig, range = which,
                chr = chr, verbose = verbose)
        } else {
            data <- reduceByFile(tiles, bList, .loadCoverageBigWig, .REDUCER, 
                chr = chr, verbose = verbose, BPPARAM = BPPARAM)
        }
        
    }
   
    ## Identify which bases pass the cutoff
    if (verbose) 
        message(paste(Sys.time(),
            "loadCoverage: applying the cutoff to the merged data"))
    
    ## Rename the object to a name that will make more sense later
    varname <- paste0(mapSeqlevels(chr, "UCSC"), "CovInfo")
    assign(varname, filterData(data = data, cutoff = cutoff, index = NULL, 
        colnames = names(dirs), filter = filter, returnMean = returnMean,
        returnCoverage = returnCoverage, totalMapped = totalMapped, 
        targetSize = targetSize, verbose = verbose))
    rm(data)    
    
    ## Save if output is specified
    if (!is.null(output)) {
        ## Automatic output name
        if (output == "auto") {
            output <- paste0(varname, ".Rdata")
        }
        
        ## Print output name
        if (verbose) 
            message(paste(Sys.time(), "loadCoverage: saving the output file", 
                output))
        
        ## Save the DataFrame
        save(list = varname, file = output, compress = "gzip")
    }
    
    ## Done
    return(get(varname))
} 

#' @export
load_coverage <- loadCoverage

## GenomicFiles functions for BAM/bigWig files
.bamMAPPER <- function(range, file, chr, verbose, isMinusStrand, ...) {
    param <- ScanBamParam(which = range, 
        flag = scanBamFlag(isMinusStrand = isMinusStrand))
    .loadCoverageBAM(file, param, chr, verbose)
}

.REDUCER <- function(mapped, ...) {
    Reduce('+', mapped)
}


.loadCoverageBAM <- function(file, param, chr, verbose) {
    if (verbose) 
        message(paste(Sys.time(), "loadCoverage: loading BAM file", 
            path(file)))
    
    ## Read the BAM file and get the coverage. Extract only the
    ## one for the chr in question.
    output <- coverage(readGAlignmentsFromBam(file, param = param))[[chr]]
        
    ## Done
    return(output)
}

.loadCoverageBigWig <- function(file, range, chr, verbose) {
    if (verbose) 
        message(paste(Sys.time(), "loadCoverage: loading BigWig file", 
            path(file)))
    
    ## Read the BAM file and get the coverage. Extract only the
    ## one for the chr in question.
    output <- import(file, selection = range, as = "RleList")[[chr]]
        
    ## Done
    return(output)
}
