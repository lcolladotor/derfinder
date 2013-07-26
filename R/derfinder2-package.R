#' Fast differential expression analysis of RNA-seq data at base-pair resolution
#'
#' Fast differential expression analysis of RNA-seq data at base-pair resolution from multiple samples. The analysis pipeline involves loading the sample BAM files using \link{loadCoverage}, calculating the F-statistics according (while adjusting for some confounders) using \link{calculateStats}, calculating the p-values and finding the regions of interest using \link{calculatePvalues}, and finally annotating them using \link[bumphunter]{annotateNearest} from the bumphunter package.
#'
#'@name derfinder2-package
#'@aliases derfinder2-package
#'@docType package
#'@author Leonardo Collado-Torres <lcollado@@jhsph.edu>
#'@references Paper coming soon.
#'@keywords package
NULL
