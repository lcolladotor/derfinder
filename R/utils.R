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

#' Run a function only using it's formal arguments. That is, stop using ...
#' recursively.
#'
#' @param fun A function to use (the actual function, not the character)
#' @param ... A set of arguments. Only those matching the formal definition of
#' \code{fun} will be used when evaluating \code{fun}.
#' @param hiddenArgs A named list of additional arguments to use that are not
#' part of the formal definition of the function.
#'
#' @keywords internal
.runFunFormal <- function(fun, ..., hiddenArgs = NULL) {
    ## Identify the formal arguments and the supplied info
    formal <- formalArgs(fun)
    args <- list(...)
    
    ## Match any of the remaining formal arguments and drop any of the
    ## extra stuff in ... which doesn't match the formal arguments
    input <- args[names(args) %in% formal]
    if(!is.null(hiddenArgs)) input <- c(input, hiddenArgs)
    
    ## Evaluate the function
    result <- do.call(fun, input)
    return(result)
}
