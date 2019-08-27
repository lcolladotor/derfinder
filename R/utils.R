#' Set an advanced argument
#'
#' @param name Name of the advanced argument to look for in ...
#' @param value The default value of the advanged argument
#' @keywords internal
#' @noRd
.advanced_argument <- function(name, value, ...) {
    args <- list(...)
    if(!name %in% names(args)) {
        return(value)
    } else {
        return(args[[name]])
    }
}

#' Run a function only using it's formal arguments. That is, stop using ...
#' recursively.
#'
#' @param fun A function to use (the actual function, not the character)
#' @param ... A set of arguments. Only those matching the formal definition of
#' `fun` will be used when evaluating `fun`.
#' @param hiddenArgs A named list of additional arguments to use that are not
#' part of the formal definition of the function.
#'
#' @importFrom methods formalArgs
#' @keywords internal
#' @noRd
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
