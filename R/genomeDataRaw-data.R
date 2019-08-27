#' Genome samples processed data
#'
#' 10kb region from chr21 processed for 31 RNA-seq samples described in 
#' [genomeInfo]. The TopHat BAM files are included in the package and this 
#' is the output of [loadCoverage] applied to it with `cutoff=NULL`. 
#' For more information check the example of [loadCoverage].
#'
#' @references
#' 1. Pickrell JK, Marioni JC, Pai AA, Degner JF, Engelhardt BE, Nkadori E, 
#' Veyrieras J-B, Stephens M, Gilad Y, Pritchard JK. Understanding mechanisms 
#' underlying human gene expression variation with RNA sequencing. Nature 2010 
#' Apr.
#'
#' 2. Montgomery SB, Sammeth M, Gutierrez-Arcelus M, Lach RP, Ingle C, Nisbett 
#' J, Guigo R, Dermitzakis ET. Transcriptome genetics using second generation 
#' sequencing in a Caucasian population. Nature 2010 Mar.
#'
#' @name genomeDataRaw
#' @docType data
#' @format  A list with two components.
#' \describe{
#' \item{coverage }{  is a DataFrame object where each column represents a 
#' sample.}
#' \item{position }{ is `NULL` because no bases were filtered.}
#' }
#' @keywords datasets
#' @seealso [loadCoverage], [genomeInfo]
NULL 
