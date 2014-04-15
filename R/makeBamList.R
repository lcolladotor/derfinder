#' Construct full paths to a group of BAM files
#'
#' For a group of samples this function creates the list of paths to the BAM 
#' files which can then be used in \link{loadCoverage}.
#' 
#' @param datadir The main directory where each of the \code{sampledirs} is a 
#' sub-directory of \code{datadir}.
#' @param sampledirs A character vector with the names of the sample 
#' directories. If \code{datadir} is \code{NULL} it is then assumed that 
#' \code{sampledirs} specifies the full path to each sample.
#' @param samplepatt If specified and \code{sampledirs} is set to \code{NULL}, 
#' then the directories matching this pattern in \code{datadir} (set to 
#' \code{.} if it's set to \code{NULL}) are used as the sample directories.
#' @param bamterm Name of the BAM file used in each sample. By default it is 
#' set to \code{accepted_hits.bam} since that is the automatic name generated 
#' when aligning with TopHat. If \code{NULL} it is then ignored when reading 
#' the BAM files. This can be useful if all the BAM files are stored in a 
#' single directory.
#'
#' @return A vector with the full paths to the BAM files and sample names 
#' stored as the vector names.
#'
#' @author Leonardo Collado-Torres
#' @export
#' @seealso \link{loadCoverage}
#' @examples
#' ## Get list of BAM files included in derfinder
#' datadir <- system.file('extdata', 'genomeData', package='derfinder')
#' dirs <- makeBamList(datadir=datadir, samplepatt='*accepted_hits.bam$', 
#'     bamterm=NULL)
#' dirs

makeBamList <- function(datadir = NULL, sampledirs = NULL, samplepatt = NULL, 
    bamterm = "accepted_hits.bam") {
    ## Determine the full paths to the sample directories
    if (!is.null(sampledirs)) {
        if (!is.null(datadir)) {
            ## Using sampledirs with datadir
            dirs <- sapply(sampledirs, function(x) {
                file.path(datadir, x)
            })
            names(dirs) <- sampledirs
        } else {
            ## Using only the sampledirs since datadir is NULL
            dirs <- sampledirs
            names(dirs) <- sampledirs
        }
    } else if (!is.null(samplepatt)) {
        if (is.null(datadir)) {
            ## This case assumes that the datadir is the current directory
            datadir <- "."
        }
        ## Identify the directories with this pattern
        dirs <- dir(path = datadir, pattern = samplepatt, full.names = TRUE)
        names(dirs) <- dir(path = datadir, pattern = samplepatt, 
            full.names = FALSE)
    } else {
        stop("Either 'samplepatt' or 'sampledirs' must be non-NULL.")
    }
    
    ## Tell R which are the BAM files
    if (!is.null(bamterm)) {
        tmp <- file.path(dirs, bamterm)
        names(tmp) <- names(dirs)
        dirs <- tmp
    }
    
    ## Done
    return(dirs)
} 
