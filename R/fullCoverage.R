#' Load the unfiltered coverage information from a group of BAM files and a list of chromosomes
#'
#' For a group of samples this function reads the coverage information for several chromosomes directly from the BAM files. Per chromosome, it merges the unfiltered coverage by sample into a DataFrame. The end result is a list with one such DataFrame objects per chromosome.
#' 
#' @param dirs A character vector with the full path to the sample BAM files. The names are used for the column names of the DataFrame. Check \link{makeBamList} for constructing \code{dirs}.
#' @param chrnums The chromosome numbers of the files to read.
#' @param chrlens The chromosome lengths in base pairs. If it's \code{NULL}, the chromosome length is extracted from the BAM files. Otherwise, it should have the same length as \code{chrnums}.
#' @param outputs This argument is passed to the \code{output} argument of \link{loadCoverage}. If \code{NULL} or \code{"auto"} it is then recycled.
#' @param mc.cores This argument is passed to \link[parallel]{mclapply}. You should use at most one core per chromosome.
#' @param verbose If \code{TRUE} basic status updates will be printed along the way.
#'
#' @return A list with one element per chromosome.
#' \describe{ Each element is a DataFrame with the coverage information produced by \link{loadCoverage}.
#' }
#'
#' @seealso \link{loadCoverage}, \link{filterData}, \link{generateReport}
#'
#' @author Leonardo Collado-Torres
#' @export
#' @importFrom parallel mclapply
#'
#' @examples
#' datadir <- system.file("extdata", "genomeData", package="derfinder2")
#' dirs <- makeBamList(datadir=datadir, samplepatt="*accepted_hits.bam$", bamterm=NULL)
#' ## Shorten the column names
#' names(dirs) <- gsub("_accepted_hits.bam", "", names(dirs))
#' 
#' ## Reading the data and filtering it is quite fast.
#' system.time(fullCov <- fullCoverage(dirs=dirs, chrnums=c("21", "22"), mc.cores=2))
#' fullCov
#'
#' ## You can then use filterData to filter the data if you want to. 
#' ## Use mclapply if you want to do so with multiple cores.
#' library("parallel")
#' mclapply(fullCov, filterData, cutoff=0, mc.cores=2)


fullCoverage <- function(dirs, chrnums, chrlens=NULL, outputs=NULL, mc.cores=getOption("mc.cores", 2L), verbose=TRUE) {
	stopifnot(length(chrlens) == length(chrnums) | is.null(chrlens))
	if(!is.null(outputs)) {
		stopifnot(length(outputs) == length(chrnums) | outputs == "auto")
		if(outputs == "auto") {
			outputs <- rep("auto", length(chrnums))
		}
	}	
		
	## Subsetting function that runs loadCoverage
	loadChr <- function(idx, dirs, chrnums, chrlens, outputs, verbose) {
		if(verbose) message(paste(Sys.time(), "fullCoverage: processing chromosome", chrnums[idx]))
		loadCoverage(dirs=dirs, chr=chrnums[idx], cutoff=NULL, chrlen=chrlens[idx], output=outputs[idx], verbose=verbose)$coverage
	}
	result <- mclapply(seq_len(length(chrnums)), loadChr, dirs=dirs, chrnums=chrnums, chrlens=chrlens, outputs=outputs, verbose=verbose, mc.cores=mc.cores)
	names(result) <- chrnums
	
	## Done
	return(result)
}