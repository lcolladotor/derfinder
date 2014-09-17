#' Arrange the ... into a named list 
#'
#' @details Taken from https://github.com/jimhester/gmailr/blob/8b1c1042a3b69248edeb05d4ebd480008ee18914/R/utils.R#L45
#' @keywords internal 
.dots  <- function (...) { eval(substitute(alist(...))) }

#' Set an advanced argument
#'
#' @param name Name of the advanced argument to look for in ...
#' @param value The default value of the advanged argument
#' @keywords internal 
.advanced_argument <- function(name, value, ...) {
    args <- .dots(...)
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
.define_cluster <- function(cores = 'mc.cores', ...) {
    args <- .dots(...)
    if('BPPARAM' %in% names(args)) {
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
