#' Define a cluster to use.
#'
#' @param cores The argument to use to define the number of cores. This is
#' useful for cases with nested parallelizations.
#' @param ... Advanced arguments are:
#' \describe{
#' \item{mc.cores }{ If 1 (default), then [SerialParam][BiocParallel::SerialParam-class] will
#' be used. If greater than 1, then it specifies the number of workers for
#' [SnowParam][BiocParallel::SnowParam-class].}
#' \item{mc.log }{ Passed to `log` when using
#' [SnowParam][BiocParallel::SnowParam-class].}
#' \item{BPPARAM.custom }{ If specified, that's the BPPARAM that will be used.}
#' }
#'
#' @importFrom BiocParallel SnowParam SerialParam
#' @return A BiocParallel *Param object
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
define_cluster <- function(cores = "mc.cores", ...) {
    args <- list(...)
    if ("BPPARAM.custom" %in% names(args)) {
        return(args$BPPARAM.custom)
    } else {
        mc.cores <- .advanced_argument(cores, getOption("mc.cores", 1L), ...)
        if (mc.cores > 1) {
            mc.log <- .advanced_argument("mc.log", TRUE, ...)
            BPPARAM <- SnowParam(workers = mc.cores, log = mc.log)
        } else {
            BPPARAM <- SerialParam()
        }
        return(BPPARAM)
    }
}
