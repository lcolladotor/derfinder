#' Construct full paths to a group of raw input files
#'
#' For a group of samples this function creates the list of paths to the raw 
#' input files which can then be used in \link{loadCoverage}. The raw input
#' files are either BAM files or BigWig files.
#' 
#' @param datadir The main directory where each of the \code{sampledirs} is a 
#' sub-directory of \code{datadir}.
#' @param sampledirs A character vector with the names of the sample 
#' directories. If \code{datadir} is \code{NULL} it is then assumed that 
#' \code{sampledirs} specifies the full path to each sample.
#' @param samplepatt If specified and \code{sampledirs} is set to \code{NULL}, 
#' then the directories matching this pattern in \code{datadir} (set to 
#' \code{.} if it's set to \code{NULL}) are used as the sample directories.
#' @param fileterm Name of the BAM or BigWig file used in each sample. By 
#' default it is set to \code{accepted_hits.bam} since that is the automatic 
#' name generated when aligning with TopHat. If \code{NULL} it is then ignored 
#' when reading the rawfiles. This can be useful if all the raw files are 
#' stored in a single directory.
#'
#' @return A vector with the full paths to the raw files and sample names 
#' stored as the vector names.
#'
#' @details This function can also be used to identify a set of BigWig files.
#'
#' @author Leonardo Collado-Torres
#' @export
#' @seealso \link{loadCoverage}
#' @examples
#' ## Get list of BAM files included in derfinder
#' datadir <- system.file('extdata', 'genomeData', package='derfinder')
#' files <- rawFiles(datadir=datadir, samplepatt='*accepted_hits.bam$', 
#'     fileterm=NULL)
#' files

rawFiles <- function(datadir = NULL, sampledirs = NULL, samplepatt = NULL, 
    fileterm = 'accepted_hits.bam') {
    ## Determine the full paths to the sample directories
    if (!is.null(sampledirs)) {
        if (!is.null(datadir)) {
            ## Using sampledirs with datadir
            files <- sapply(sampledirs, function(x) {
                file.path(datadir, x)
            })
            names(files) <- sampledirs
        } else {
            ## Using only the sampledirs since datadir is NULL
            files <- sampledirs
            names(files) <- sampledirs
        }
    } else if (!is.null(samplepatt)) {
        if (is.null(datadir)) {
            ## This case assumes that the datadir is the current directory
            datadir <- '.'
        }
        ## Identify the directories with this pattern
        files <- dir(path = datadir, pattern = samplepatt, full.names = TRUE)
        names(files) <- dir(path = datadir, pattern = samplepatt, 
            full.names = FALSE)
    } else {
        stop("Either 'samplepatt' or 'sampledirs' must be non-NULL.")
    }
    
    ## Tell R which are the BAM files
    if (!is.null(fileterm)) {
        tmp <- file.path(files, fileterm)
        names(tmp) <- names(files)
        files <- tmp
    }
    
    ## Done
    return(files)
}
