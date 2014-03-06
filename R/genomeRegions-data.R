#' Candidate DERs for example data
#'
#' Candidate Differentially Expressed Regions (DERs) for the example data. For more information check \link{calculatePvalues}.
#'
#'
#' @name genomeRegions
#' @docType data
#' @format  A list with four components.
#' \describe{
#' \item{regions }{ a GRanges object with the candidate DERs.}
#' \item{nullStats }{ a numeric Rle with the mean F-statistics for the null DERs found from the permutations.}
#' \item{nullWidths }{ an integer Rle with the width of each null candidate DER.}
#' \item{nullPermutation }{ an integer Rle with the permutation number for each candidate DER. It identifies which permutation cycle created the null candidate DER.}
#' }
#'@keywords datasets
#'@seealso \link{calculatePvalues}
NULL
