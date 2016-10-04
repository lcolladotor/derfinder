#' Define a cluster to use.
#'
#' @param cores The argument to use to define the number of cores. This is
#' useful for cases with nested parallelizations.
#' @param ... Advanced arguments are:
#' \describe{
#' \item{mc.cores }{ If 1 (default), then \link[BiocParallel]{SerialParam} will 
#' be used. If greater than 1, then it specifies the number of workers for
#' \link[BiocParallel]{SnowParam}.}
#' \item{mc.outfile }{ Passed to \code{outfile} when using
#' \link[BiocParallel]{SnowParam}.}
#' \item{BPPARAM.custom }{ If specified, that's the BPPARAM that will be used.}
#' }
#'
#' @importFrom BiocParallel SnowParam SerialParam
#' @export
#'
#' @author Leonardo Collado-Torres
#'
#' @details This function is used internally in many functions.
#'
#' @examples
#' ## Use SerialParam()
#' define_cluster(mc.cores = 1)
#'
#' ## Note that BPPARAM.custom takes precedence
#' define_cluster(mc.cores = 2, BPPARAM.custom = BiocParallel::SerialParam())
#'

define_cluster <- function(cores = 'mc.cores', ...) {
    args <- list(...)
    if('BPPARAM.custom' %in% names(args)) {
        return(args$BPPARAM.custom)
    } else {
        mc.cores <- .advanced_argument(cores, getOption('mc.cores', 1L), ...)
        if(mc.cores > 1) {
            mc.outfile <- .advanced_argument('mc.outfile',
                Sys.getenv('SGE_STDERR_PATH'), ...)
            BPPARAM <- SnowParam(workers = mc.cores, outfile = mc.outfile)
        } else {
            BPPARAM <- SerialParam()
        }
        return(BPPARAM)
    }
}
