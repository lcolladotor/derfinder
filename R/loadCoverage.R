#' Load the coverage information from a group of BAM files
#'
#' For a group of samples this function reads the coverage information for a specific chromosome directly from the BAM files. It then merges them into a DataFrame and removes the bases that do not pass the cutoff.
#' 
#' @param dirs A character vector with the full path to the sample BAM files. The names are used for the column names of the DataFrame. Check \link{makeBamList} for constructing \code{dirs}. \code{dirs} can also be a \code{BamFileList} object created with \link[Rsamtools]{BamFileList}.
#' @param chr Chromosome to read. Should be in simple format. For example, use X and not chrX.
#' @param cutoff This argument is passed to \link{filterData}.
#' @param bai The full path to the BAM index files. If \code{NULL} it is assumed that the BAM index files are in the same location as the BAM files and that they have the .bai extension. Ignored if \code{dirs} is a \code{BamFileList} object.
#' @param chrlen The chromosome length in base pairs. If it's \code{NULL}, the chromosome length is extracted from the BAM files.
#' @param output If \code{NULL} then no output is saved in disk. If \code{auto} then an automatic name is constructed (chrXCovInfo.Rdata for example). If another character is specified, then that name is used for the output file.
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
#' @importFrom Rsamtools BamFileList scanBamHeader ScanBamParam path
#' @importFrom GenomicAlignments readGAlignmentsFromBam
#' @importFrom IRanges IRanges RangesList
#' @importMethodsFrom GenomicRanges coverage
#' @importMethodsFrom Rsamtools names
#' @examples
#' datadir <- system.file("extdata", "genomeData", package="derfinder")
#' dirs <- makeBamList(datadir=datadir, samplepatt="*accepted_hits.bam$", bamterm=NULL)
#' ## Shorten the column names
#' names(dirs) <- gsub("_accepted_hits.bam", "", names(dirs))
#' 
#' ## Reading the data and filtering it is quite fast.
#' system.time(data <- loadCoverage(dirs=dirs, chr="21", cutoff=0))
#' data
#'
#' ## For other examples this data is included in the package
#' identical(data, genomeData)
#'
#' ## The data is compact enough to be loaded in memory
#' print(object.size(data), units="Kb")
#'
#' \dontrun{
#' ## Read the coverage without applying any cutoff.
#' ## This can be useful for downstream analysis including coverage plots.
#' system.time(data2 <- loadCoverage(dirs=dirs, chr="21", cutoff=NULL))
#'
#' ## Compare vs raw data provided in the package
#' identical(data2, genomeDataRaw)
#' 
#' ## Note that the object size is pretty much the same due to the Rle compression
#' print(object.size(data2), units="Kb")
#' }

loadCoverage <- function(dirs, chr, cutoff=NULL, bai=NULL, chrlen=NULL, output=NULL, verbose=TRUE) {
	## Do the indexes exist?
	if(is(dirs, "BamFileList")) {
		bList <- dirs
	} else {
		if(is.null(bai)) {
			bai <- paste0(dirs, ".bai")
		}	
		if(all(file.exists(bai))) {
			bList <- BamFileList(dirs, bai)
		} else {
			stop("Not all BAM files have a BAM index. If the BAM index files are in a separate directory from the BAM files or are not named as 'bamFile.bam.bai' then consider using the 'bai' argument.")
		}
	}
	
		
	## Determine the chromosome length
	if(is.null(chrlen)) {
		## This assumes that all the BAM files are from the same organism.
		if(verbose) message(paste(Sys.time(), "loadCoverage: finding chromosome lengths"))
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
		if(verbose) message(paste(Sys.time(), "loadCoverage: loading BAM file", path(x)))
	
		## Read the BAM file and get the coverage. Extract only the one for the chr in question.
## In R-devel, readBamGappedAlignments is deprecated! Check ?readGAlignmentsFromBam
		#output <- coverage(readBamGappedAlignments(x, param=param))[[chr]]
		output <- coverage(readGAlignmentsFromBam(x, param=param))[[chr]]
		
		
		## Done
		return(output)
	})
	
	## Identify which bases pass the cutoff
	if(verbose) message(paste(Sys.time(), "loadCoverage: applying the cutoff to the merged data"))
	
	res <- filterData(data=data, cutoff=cutoff, index=NULL, colnames=names(dirs), verbose=verbose)
	rm(data)
	gc()
	
	
	## Save if output is specified
	if(!is.null(output)) {
		## Rename the object to a name that will make more sense later
		chrnum <- gsub("chr", "", chr)
		varname <- paste0("chr", chrnum, "CovInfo")
		assign(varname, res)
		
		## Automatic output name
		if(output=="auto") {
			output <- paste0(varname, ".Rdata")
		}
		
		## Print output name
		if(verbose) message(paste(Sys.time(), "loadCoverage: saving the output file", output))
		
		## Save the DataFrame
		save(list=varname, file=output, compress="gzip")
	}
	
	## Done
	return(res)	
}
