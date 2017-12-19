#' Load the coverage information from a group of BAM files
#'
#' For a group of samples this function reads the coverage information for a 
#' specific chromosome directly from the BAM files. It then merges them into a 
#' DataFrame and removes the bases that do not pass the cutoff.
#' 
#' @param files A character vector with the full path to the sample BAM files
#' (or BigWig files). 
#' The names are used for the column names of the DataFrame. Check 
#' \link{rawFiles} for constructing \code{files}. \code{files} can also be a 
#' \link[Rsamtools]{BamFileList}, \link[Rsamtools]{BamFile}, 
#' \link[rtracklayer]{BigWigFileList}, or \link[rtracklayer]{BigWigFile} object.
#' @param chr Chromosome to read. Should be in the format matching the one used
#' in the raw data.
#' @param cutoff This argument is passed to \link{filterData}.
#' @inheritParams filterData
#' @param chrlen The chromosome length in base pairs. If it's \code{NULL}, the 
#' chromosome length is extracted from the BAM files.
#' @param output If \code{NULL} then no output is saved in disk. If \code{auto} 
#' then an automatic name is constructed using UCSC names (chrXCovInfo.Rdata
#' for example). If another character is specified, then that name is used for #' the output file.
#' @param bai The full path to the BAM index files. If \code{NULL} it is 
#' assumed that the BAM index files are in the same location as the BAM files 
#' and that they have the .bai extension. Ignored if \code{files} is a 
#' \code{BamFileList} object or if \code{inputType=='BigWig'}.
#' @param ... Arguments passed to other methods and/or advanced arguments.
#' Advanced arguments:
#' \describe{
#' \item{verbose }{ If \code{TRUE} basic status updates will be printed along 
#' the way.}
#' \item{inputType }{ Has to be either \code{bam} or \code{BigWig}. It specifies
#' the format of the raw data files. By default it's set to \code{bam} before
#' trying to guess it from the file names.}
#' \item{tilewidth }{ When specified, \link[GenomicRanges]{tileGenome} is used 
#' to break up the chromosome into chunks. We don't recommend this for BAM
#' files as the coverage in the borders of the chunks might be slightly off.}
#' \item{which }{ \code{NULL} by default. When a \code{GRanges} is specified, 
#' this specific region of the genome is loaded instead of the full chromosome.}
#' \item{fileStyle }{ The naming style of the chromosomes in the input files. 
#' If the global option \code{chrsStyle} is not set, the naming style won't be 
#' changed. This option is useful when you want to use a specific naming style
#' but the raw files use another style.}
#' \item{protectWhich }{ When not \code{NULL} and \code{which} is specified, 
#' this argument specifies by how much the ranges in \code{which} will be 
#' expanded. This helps get the same base level coverage you would get from 
#' reading the coverage for a full chromosome from BAM files. Otherwise some 
#' reads might be excluded and thus the base level coverage will be lower. 
#' \code{NULL} by default.}
#' \item{drop.D }{ Whether to drop the bases with 'D' in the CIGAR strings
#' or to include them. Only used with BAM files. \code{FALSE} by default.}
#' \item{sampleNames }{ Column names to be used the samples. By default it's
#' \code{names(files)}.}
#' }
#' Passed to \link{extendedMapSeqlevels}, \link{define_cluster}, 
#' \link[Rsamtools]{scanBamFlag} and \link{filterData}. 
#' Note that \link{filterData} is used internally 
#' by \link{loadCoverage} and has the important arguments \code{totalMapped} 
#' and \code{targetSize} which can be used to normalize the coverage by
#' library size. See \link{getTotalMapped} for calculating \code{totalMapped}.
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
#' @details
#' 
#' The \code{...} argument can be used to control which alignments to consider
#' when reading from BAM files. See \link[Rsamtools]{scanBamFlag}.
#'
#' Parallelization for loading the data in chunks is used only used when 
#' \code{tilewidth} is specified. You may use up to one core per tile.
#'
#' If you set the advanced argument \code{drop.D = TRUE}, bases with CIGAR 
#' string "D" (deletion from reference) will be excluded from the base-level
#' coverage calculation.
#'
#' If you are working with data from an organism different from 'Homo sapiens'
#' specify so by setting the global 'species' and 'chrsStyle' options. For 
#' example:
#' \code{options(species = 'arabidopsis_thaliana')}
#' \code{options(chrsStyle = 'NCBI')}
#' 
#' @seealso \link{fullCoverage}, \link{getTotalMapped}
#' 
#' @author Leonardo Collado-Torres, Andrew Jaffe
#' @export
#' @importFrom BiocGenerics path
#' @importFrom Rsamtools BamFileList scanBamHeader ScanBamParam
#' scanBamFlag BamFile
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom IRanges IRanges RangesList
#' @importFrom rtracklayer BigWigFileList BigWigFile
#' @importFrom GenomeInfoDb seqlevels
#' mapSeqlevels
#' @importFrom GenomicRanges tileGenome
#' @importFrom GenomicFiles reduceByFile
#' @import S4Vectors
#' @importMethodsFrom GenomicRanges coverage isDisjoint reduce width resize
#' @importMethodsFrom Rsamtools names
#' @importMethodsFrom rtracklayer import import.bw
#' @importFrom methods is
#' @examples
#' datadir <- system.file('extdata', 'genomeData', package='derfinder')
#' files <- rawFiles(datadir = datadir, samplepatt = '*accepted_hits.bam$', 
#'     fileterm = NULL)
#' ## Shorten the column names
#' names(files) <- gsub('_accepted_hits.bam', '', names(files))
#'  
#' ## Read and filter the data, only for 2 files
#' dataSmall <- loadCoverage(files = files[1:2], chr = '21', cutoff = 0)
#'
#' \dontrun{
#' ## Export to BigWig files
#' createBw(list('chr21' = dataSmall))
#'
#' ## Load data from BigWig files
#' dataSmall.bw <- loadCoverage(c(ERR009101 = 'ERR009101.bw', ERR009102 = 
#'     'ERR009102.bw'), chr = 'chr21')
#' 
#' ## Compare them
#' mapply(function(x, y) { x - y }, dataSmall$coverage, dataSmall.bw$coverage)
#'
#' ## Note that the only difference is the type of Rle (integer vs numeric) used
#' ## to store the data.
#'
#' }

loadCoverage <- function(files, chr, cutoff = NULL, filter = 'one', 
    chrlen = NULL, output = NULL, bai = NULL, ...) {
    
    stopifnot(is.character(chr))
    
    ## Advanged arguments
# @param verbose If \code{TRUE} basic status updates will be printed along the 
# way.
    verbose <- .advanced_argument('verbose', TRUE, ...)

    
# @param inputType Has to be either \code{bam} or \code{BigWig}. It specifies
# the format of the raw data files.
    inputType <- .advanced_argument('inputType', 'bam', ...)


    ## Guess the input type if it's not supplied
    if(is(files, 'BigWigFileList') | is(files, 'BigWigFile')) {
        inputType <- 'BigWig'
    } else if (is(files, 'BamFileList') | is(files, 'BamFile')) {
        inputType <- 'bam'
    } else if(all(grepl('bw$|bigwig$', tolower(files)))) {
        inputType <- 'BigWig'
    }    
    stopifnot(inputType %in% c('bam', 'BigWig'))

# @param tilewidth When specified, \link[GenomicRanges]{tileGenome} is used to
# break up the chromosome into chunks.
    tilewidth <- .advanced_argument('tilewidth', NULL, ...)

  
# @param which \code{NULL} by default. When a \code{GRanges} is specified, 
# this specific region of the genome is loaded instead of the full chromosome.
    which <- .advanced_argument('which', NULL, ...)


# @param fileStyle The naming style of the chromosomes in the input files. If 
# the global option 'chrsStyle' is not set, the naming style won't be changed.
    fileStyle <- .advanced_argument('fileStyle', getOption('chrsStyle',
        NULL), ...)


# @param protectWhich When not \code{NULL} and \code{which} is specified, this argument specifies by how much the ranges in \code{which} will be expanded.
# This helps get the same base level coverage you would get from reading the coverage for a full chromosome. Otherwise some reads might be excluded and thus the base level coverage will be lower.
    protectWhich <- .advanced_argument('protectWhich', NULL, ...)


    ## Assign naming style
    chr <- extendedMapSeqlevels(chr, style = fileStyle, ...)
    
    ## Fix 'which'
    if(!is.null(which)) {
        stopifnot(is(which, 'GRanges'))
        which <- renameSeqlevels(which, extendedMapSeqlevels(seqlevels(which),
            style = fileStyle, ...))
        
        if(!is.null(protectWhich)) {
            stopifnot(protectWhich >= 0)
            which <- resize(which, width(which) + protectWhich, fix = 'center')
        }     
        
        if(!isDisjoint(which)) {
            ## Prevent counting more than once the coverage for a given region
            which <- reduce(which)
        }
    }    

    ## Do the indexes exist?
    if (is(files, 'BamFileList')) {
        bList <- files
    } else if (is(files, 'BamFile')) {
        bList <- BamFileList(files)
    } else if (is(files, 'BigWigFileList')) {
        bList <- files
    } else if (is(files, 'BigWigFile')) {
        bList <- BigWigFileList(path(files))
    } else if (inputType == 'bam'){
        if (is.null(bai)) {
            bai <- paste0(files, '.bai')
        }
        if (all(file.exists(bai))) {
            bList <- BamFileList(files, bai)
        } else {
            stop("Not all BAM files have a BAM index. If the BAM index files are in a separate directory from the BAM files or are not named as 'bamFile.bam.bai' then consider using the 'bai' argument.")
        }
    } else if (inputType == 'BigWig') {
        bList <- BigWigFileList(files)
        if(.Platform$OS.type == 'windows') {
            warning('As of rtracklayer 1.25.16, BigWig is not supported on Windows. Thus loading data from BigWig files will most likely fail!')
        }
    }
    
    
    ## Determine the chromosome length
    if (is.null(chrlen)) {
        ## This assumes that all the BAM files are from the same
        ## organism.
        if (verbose) 
            message(paste(Sys.time(),
                'loadCoverage: finding chromosome lengths'))
        if (inputType == 'bam') {
            clengths <- scanBamHeader(bList[[1]])$targets       
        } else if (inputType == 'BigWig') {
            clengths <- seqlengths(bList[[1]])
        }
       
        if (!chr %in% names(clengths)) {
            stop(paste("'chr' is not correctly specified. Valid options are:", 
                paste(names(clengths), collapse = ', ')))
        }
        chrlen <- clengths[chr]
    }
    
    ## Make tiles if using GenomicFiles
    if(!is.null(tilewidth)) {
        tiles <- tileGenome(chrlen, tilewidth = tilewidth)
        
        ## Define cluster
        BPPARAM <- define_cluster(...)
    }    
    
    ## Construct the objects so only the chr of interest is read
    ## from the BAM/BigWig file
    if(is.null(which) | !is(which, 'GRanges')) {
        which <- GRanges(seqnames=chr, ranges=IRanges(1, chrlen))
    }    
    if(inputType == 'bam') {
        param <- ScanBamParam(which = which, 
            flag = .runFunFormal(scanBamFlag, ...))
            
        # @param drop.D Whether to drop the bases with 'D' in the CIGAR strings
        # or to include them.
        drop.D <- .advanced_argument('drop.D', FALSE, ...)
        
        ## Read in the data for all the chrs
        if(is.null(tilewidth)) {
            data <- lapply(bList, .loadCoverageBAM, param = param, chr = chr,
                verbose = verbose, drop.D.ranges = drop.D)
        } else {
            data <- reduceByFile(tiles, bList, .bamMAPPER, .REDUCER, 
                chr = chr, verbose = verbose, BPPARAM = BPPARAM, 
                drop.D.ranges = drop.D, ...)
        }
        
    } else if (inputType == 'BigWig') {
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
            'loadCoverage: applying the cutoff to the merged data'))
    
    ## Rename the object to a name that will make more sense later
    varname <- paste0(mapSeqlevels(chr, 'UCSC'), 'CovInfo')

# @param sampleNames Column names to be used the samples
    sampleNames <- .advanced_argument('sampleNames', names(files), ...)
    
    assign(varname, filterData(data = data, cutoff = cutoff, index = NULL, 
        colnames = sampleNames, filter = filter, ...))
    rm(data)    
    
    ## Save if output is specified
    if (!is.null(output)) {
        ## Automatic output name
        if (output == 'auto') {
            output <- paste0(varname, '.Rdata')
        }
        
        ## Print output name
        if (verbose) 
            message(paste(Sys.time(), 'loadCoverage: saving the output file', 
                output))
        
        ## Save the DataFrame
        save(list = varname, file = output, compress = 'gzip')
    }
    
    ## Done
    return(get(varname))
}


## GenomicFiles functions for BAM/BigWig files
.bamMAPPER <- function(range, file, chr, verbose, drop.D.ranges, ...) {
    param <- ScanBamParam(which = range, flag = .runFunFormal(scanBamFlag,
        ...))
    .loadCoverageBAM(file, param, chr, verbose, drop.D.ranges)
}

.REDUCER <- function(mapped, ...) {
    Reduce('+', mapped)
}


.loadCoverageBAM <- function(file, param, chr, verbose, drop.D.ranges) {
    if (verbose) 
        message(paste(Sys.time(), 'loadCoverage: loading BAM file', 
            path(file)))
    
    ## Read the BAM file and get the coverage. Extract only the
    ## one for the chr in question.
    output <- coverage(readGAlignments(file, param = param), 
        drop.D.ranges = drop.D.ranges)[[chr]]
        
    ## Done
    return(output)
}

.loadCoverageBigWig <- function(range, file, chr, verbose) {
    if (verbose) 
        message(paste(Sys.time(), 'loadCoverage: loading BigWig file', 
            path(file)))
    
    ## Read the BAM file and get the coverage. Extract only the
    ## one for the chr in question.
    output <- import(file, selection = range, as = 'RleList')[[chr]]
        
    ## Done
    return(output)
}
