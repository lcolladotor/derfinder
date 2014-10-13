#' Set an advanced argument
#'
#' @param name Name of the advanced argument to look for in ...
#' @param value The default value of the advanged argument
#' @keywords internal 
.advanced_argument <- function(name, value, ...) {
    args <- list(...)
    if(!name %in% names(args)) {
        return(value)
    } else {
        return(args[[name]])
    }
}

#' Define a cluster to use.
#'
#' @param cores The argument to use to define the number of cores. This is
#' useful for cases with nested parallelizations.
#'
#' @importFrom BiocParallel SnowParam SerialParam
#' @keywords internal
#'
#' @details
#' Use If you define 'BPPARAM.custom' then this will be used instead of the
#' default SnowParam() [if the number of cores is > 1] or SerialParam()
.define_cluster <- function(cores = 'mc.cores', ...) {
    args <- list(...)
    if('BPPARAM.custom' %in% names(args)) {
        return(args$BPPARAM)
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
