#' Load the coverage information from a group of BAM files
#'
#' For a group of samples this function reads the coverage information for a specific chromosome directly from the BAM files. It then merges them into a DataFrame and removes the bases that do not pass the cutoff.
#' 
#' @param chr Chromosome to read. Should be in simple format. For example, use X and not chrX.
#' @param datadir The main directory where each of the \code{sampledirs} is a sub-directory of \code{datadir}.
#' @param sampledirs A character vector with the names of the sample directories. If \code{datadir} is \code{NULL} it is then assumed that \code{sampledirs} specifies the full path to each sample.
#' @param samplepatt If specified and \code{sampledirs} is set to \code{NULL}, then the directories matching this pattern in \code{datadir} (set to \code{.} if it's set to \code{NULL}) are used as the sample directories.
#' @param cutoff Per base pair, at least one sample has to have coverage greater than \code{cutoff} to be included in the result.
#' @param chrlen The chromosome length in base pairs.
#' @param bamterm Name of the BAM file used in each sample. By default it is set to \code{accepted_hits.bam} since that is the automatic name generated when aligning with TopHat. If \code{NULL} it is then ignored when reading the BAM files. This can be useful if all the BAM files are stored in a single directory.
#' @param output If \code{NULL} then no output is saved in disk. If \code{auto} then an automatic name is constructed (chrXDF.Rdata for example). If another character is specified, then that name is used for the output file.
#' @param verbose If \code{TRUE} basic status updates will be printed along the way.
#'
#' @return A list with two components.
#' \describe{
#' \item{coverage }{ is a DataFrame object where each column represents a sample. The number of rows depends on the number of base pairs that passed the cutoff and the information stored is the coverage at that given base.}
#' \item{position }{  is a logical Rle with the positions of the chromosome that passed the cutoff.}
#' }
#'
#' @author Leonardo Collado-Torres, Andrew Jaffe
#' @export
#' @importFrom Rsamtools BamFileList scanBamHeader ScanBamParam path readBamGappedAlignments
#' @importFrom IRanges IRanges RangesList
#' @importMethodsFrom GenomicRanges coverage
#' @importMethodsFrom Rsamtools names
#' @examples
#' datadir <- system.file("extdata", "brainData", package="derfinder2")
#' ## Reading the data and filtering it is quite fast.
#' system.time(data <- loadCoverage(chr="21", datadir=datadir, samplepatt="*accepted_hits.bam$", bamterm=NULL))
#' ## Shorten the column names
#' colnames(data$coverage) <- gsub("_accepted_hits.bam", "", colnames(data$coverage))
#' data
#' ## The data is compact enough to be loaded in memory
#' print(object.size(data), units="Kb")

loadCoverage <- function(chr, datadir=NULL, sampledirs=NULL, samplepatt=NULL, cutoff=5, chrlen=NULL, bamterm="accepted_hits.bam", output=NULL, verbose=TRUE) {
	## Determine the full paths to the sample directories
	if(!is.null(sampledirs)) {
		if(!is.null(datadir)) {
			## Using sampledirs with datadir
			dirs <- sapply(sampledirs, function(x) { file.path(datadir, x) })
			names(dirs) <- sampledirs
		} else {
			## Using only the sampledirs since datadir is NULL
			dirs <- sampledirs
			names(dirs) <- sampledirs
		}
	} else if (!is.null(samplepatt)) {
		if(is.null(datadir)) {
			## This case assumes that the datadir is the current directory
			datadir <- "."
		}
		## Identify the directories with this pattern
		dirs <- dir(path=datadir, pattern=samplepatt, full.names=TRUE)
		names(dirs) <- dir(path=datadir, pattern=samplepatt, full.names=FALSE)
	} else {
		stop("Either 'samplepatt' or 'sampledirs' must be non-NULL.")
	}
	
	## Tell R which are the BAM files
	if(!is.null(bamterm)) {
		dirs <- file.path(dirs, bamterm)
	}
	## Do the indexes exist?
	bai <- paste0(dirs, ".bai")
	if(all(file.exists(bai))) {
		bList <- BamFileList(dirs, bai)
	} else {
		bList <- BamFileList(dirs)
	}
		
	## Determine the chromosome length
	if(is.null(chrlen)) {
		## This assumes that all the BAM files are from the same organism.
		clengths <- scanBamHeader(bList[[1]])$targets
		if(!chr %in% names(clengths)) {
			stop(paste("'chr' is not correctly specified. Valid options are:", paste(names(clengths), collapse=", ")))
		}
		chrlen <- clengths[chr]
	}
	
	## Construct the objects so only the chr of interest is read from the BAM file
	which <- RangesList(IRanges(1, chrlen))
	names(which) <- chr
	param <- ScanBamParam(which=which)
	
	## Read in the data for all the chrs
	data <- lapply(bList, function(x) {
		if(verbose) message(paste("loadCoverage: loading BAM file", path(x)))
	
		## Read the BAM file and get the coverage. Extract only the one for the chr in question.
		output <- coverage(readBamGappedAlignments(x, param=param))[[chr]]
		
		## Done
		return(output)
	})
	
	## Identify which bases pass the cutoff
	if(verbose) message("loadCoverage: applying the cutoff to the merged data")
	
	res <- filterData(data=data, cutoff=cutoff, index=NULL, colnames=names(dirs), verbose=verbose)
	
	## Save if output is specified
	if(!is.null(output)) {
		## Rename the object to a name that will make more sense later
		varname <- paste0("chr", chr, "DF")
		assign(varname, res)
		
		## Automatic output name
		if(output=="auto") {
			output <- paste0(varname, ".Rdata")
		}
		
		## Print output name
		if(verbose) message(paste("loadCoverage: saving the output file", output))
		
		## Save the DataFrame
		save(list=varname, file=output, compress="gzip")
	}
	
	## Done
	return(res)	
}
