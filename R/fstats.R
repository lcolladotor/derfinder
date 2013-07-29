#' Calculate F-statistics per base
#'
#' A function for calculating F-statistics on the rows of dat, comparing the models mod (alternative) and mod0 (null). This is a helper function for \link{calculateStats}.
#' 
#' @param dat A numeric matrix of size n by m. There should be n bases of the genome and m samples.
#' @param mod The design matrix for the alternative model. Should be m by p where p is the number of covariates (normally also including the intercept).
#' @param mod0 The deisgn matrix for the null model. Should be m by p_0.
#'
#' @return A vector with the F-statistics per base (length n).
#'
#' @author Jeff Leek, Leonardo Collado-Torres
#' @export
#' @seealso \link{calculateStats}
#' @examples
#' ## Generate random data
#' set.seed(20130722)
#' dat <- matrix(rnorm(1000, 5), ncol=10)
#' 
#' ## Define the model matrices
#' group <- sample(letters[1:2], 10, TRUE)
#' mod <- model.matrix(~ group)
#' mod0 <- model.matrix(~ 0 + rep(1, 10))
#' 
#' ## Get the F-statistics
#' result <- fstats(dat, mod, mod0)
#' hist(result, freq=FALSE, main="Distribution of F statistics")
#' 

fstats <- function(dat, mod, mod0){
	# A function for calculating F-statistics
	# on the rows of dat, comparing the models
	# mod (alternative) and mod0 (null). 	
	n <- dim(dat)[2]
	m <- dim(dat)[1]
	df1 <- dim(mod)[2]
	df0 <- dim(mod0)[2]
	p <- rep(0,m)
	Id <- diag(n)

	resid <- dat %*% (Id - mod %*% solve(t(mod) %*% mod) %*% t(mod))
	resid0 <- dat %*% (Id - mod0 %*% solve(t(mod0) %*% mod0) %*% t(mod0))

	rss1 <- (resid*resid) %*% rep(1,n)
	rss0 <- (resid0*resid0) %*% rep(1,n)

	fstats <- ((rss0 - rss1)/(df1-df0))/(rss1/(n-df1))
	return(drop(fstats))
}
