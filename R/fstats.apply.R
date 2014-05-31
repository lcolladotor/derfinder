#' Calculate F-statistics per base by extracting chunks from a DataFrame
#'
#' Extract chunks from a DataFrame and get the F-statistics on the rows of 
#' \code{data}, comparing the models \code{mod} (alternative) and \code{mod0} 
#' (null). This is a helper function for \link{calculateStats} and 
#' \link{calculatePvalues}.
#'
#' @param index An index (logical Rle is the best for saving memory) indicating 
#' which rows of the DataFrame to use.
#' @param data The DataFrame containing the coverage information. Normally 
#' stored in \code{coveragePrep$coverageProcessed} from 
#' \link{preprocessCoverage}. Could also be the full data from 
#' \link{loadCoverage}.
#' @param mod The design matrix for the alternative model. Should be m by p 
#' where p is the number of covariates (normally also including the intercept).
#' @param mod0 The design matrix for the null model. Should be m by p_0.
#' @param adjustF A single value to adjust that is added in the denominator of 
#' the F-stat calculation. Useful when the Residual Sum of Squares of the 
#' alternative model is very small.
#' @param lowMemDir The directory where the processed chunks are saved when 
#' using \link{preprocessCoverage} with a specified \code{lowMemDir}.
#'
#' @details If \code{lowMemDir} is specified then \code{index} is expected to 
#' specify the chunk number.
#'
#' @return A numeric Rle with the F-statistics per base for the chunk in 
#' question.
#'
#' @author Jeff Leek, Leonardo Collado-Torres
#' @export
#' @importMethodsFrom IRanges as.data.frame as.matrix
#' @seealso \link{calculateStats}, \link{calculatePvalues}
#'
#' @examples
#' ## Create the model matrices
#' mod <- model.matrix(~ genomeInfo$pop)
#' mod0 <- model.matrix(~ 0 + rep(1, nrow(genomeInfo)))
#' ## Run the function
#' fstats.output <- fstats.apply(data=genomeData$coverage, mod=mod, mod0=mod0)
#' fstats.output
#' 

fstats.apply <- function(index=Rle(TRUE, nrow(data)), data, mod, mod0, 
    adjustF=0, lowMemDir=NULL) {
    ## Load the chunk file
    if(!is.null(lowMemDir)) {
        chunkProcessed <- NULL
        load(file.path(lowMemDir, paste0("chunk", index, ".Rdata")))
        data <- chunkProcessed
        
        ##  Transform to a regular matrix
        dat <- as.matrix(as.data.frame(data))
    } else{
        ##  Subset the DataFrame to the current chunk and transform to a regular matrix
        dat <- as.matrix(as.data.frame(data[index, ]))
    }
    rm(data)
    
    # A function for calculating F-statistics
    # on the rows of dat, comparing the models
    # mod (alternative) and mod0 (null).
    n <- dim(dat)[2]
    m <- dim(dat)[1]
    p <- rep(0L, m)
    Id <- diag(n)
    nVec <- rep(1L, n)

    ## Calculate rss1
    resid <- dat %*% (Id - mod %*% solve(t(mod) %*% mod) %*% t(mod))
    rss1 <- (resid*resid) %*% nVec
    rm(resid)


    ## Calculate rss2
    resid0 <- dat %*% (Id - mod0 %*% solve(t(mod0) %*% mod0) %*% t(mod0))
    rss0 <- (resid0*resid0) %*% nVec
    rm(resid0, nVec, Id)


    ## Get the F-stats
    df1 <- dim(mod)[2]
    df0 <- dim(mod0)[2]
    fstats <- Rle(drop(((rss0 - rss1) / (df1 - df0)) / (adjustF + rss1 /
        (n - df1))))
    
    ## Done
    return(fstats)
}
