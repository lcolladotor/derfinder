#' Deprecated functions in package \sQuote{derfinder}
#'
#' These functions are provided for compatibility with older versions
#' of \sQuote{derfinder} only, and will be defunct at the next release.
#'
#' @details 
#'
#' The following functions are deprecated and will be made defunct; use
#'  the replacement indicated below:
#'  \itemize{
#'
#'    \item{advanced_arg: \code{\link{advancedArg}}}
#'    \item{analyze_chr: \code{\link{analyzeChr}}}
#'    \item{annotate_regions: \code{\link{annotateRegions}}}
#'    \item{calculate_pvalues: \code{\link{calculatePvalues}}}
#'    \item{calculate_stats: \code{\link{calculateStats}}}
#'    \item{coerce_gr: \code{\link{coerceGR}}}
#'    \item{collapse_full_coverage: \code{\link{collapseFullCoverage}}}
#'    \item{coverage_to_exon: \code{\link{coverageToExon}}}
#'    \item{create_bw: \code{\link{createBw}}}
#'    \item{create_bw_sample: \code{\link{createBwSample}}}
#'    \item{extended_map_seqlevels: \code{\link{extendedMapSeqlevels}}}
#'    \item{filter_data: \code{\link{filterData}}}
#'    \item{find_regions: \code{\link{findRegions}}}
#'    \item{full_coverage: \code{\link{fullCoverage}}}
#'    \item{get_region_coverage: \code{\link{getRegionCoverage}}}
#'    \item{load_coverage: \code{\link{loadCoverage}}}
#'    \item{make_genomic_state: \code{\link{makeGenomicState}}}
#'    \item{make_models: \code{\link{makeModels}}}
#'    \item{merge_results: \code{\link{mergeResults}}}
#'    \item{preprocess_coverage: \code{\link{preprocessCoverage}}}
#'    \item{raw_files: \code{\link{rawFiles}}}
#'    \item{region_matrix: \code{\link{regionMatrix}}}
#'    \item{sample_depth: \code{\link{sampleDepth}}}
#'
#'  }
#'
#' @name derfinder-deprecated
#' @aliases derfinder-deprecated

#' @export
advanced_arg <- function() { .Deprecated('advancedArg')}

#' @export
analyze_chr <- function() { .Deprecated('analyzeChr')}

#' @export
annotate_regions <- function() { .Deprecated('annotateRegions')}

#' @export
calculate_pvalues <- function() { .Deprecated('calculatePvalues')}

#' @export
calculate_stats <- function() { .Deprecated('calculateStats')}

#' @export
coerce_gr <- function() { .Deprecated('coerceGR')}

#' @export
collapse_full_coverage <- function() { .Deprecated('collapseFullCoverage')}

#' @export
coverage_to_exon <- function() { .Deprecated('coverageToExon')}

#' @export
create_bw <- function() { .Deprecated('createBw')}

#' @export
create_bw_sample <- function() { .Deprecated('createBwSample')}

#' @export
extended_map_seqlevels <- function() { .Deprecated('extendedMapSeqlevels')}

#' @export
filter_data <- function() { .Deprecated('filterData')}

#' @export
find_regions <- function() { .Deprecated('findRegions')}

#' @export
full_coverage <- function() { .Deprecated('fullCoverage')}

#' @export
get_region_coverage <- function() { .Deprecated('getRegionCoverage')}

#' @export
load_coverage <- function() { .Deprecated('loadCoverage')}

#' @export
make_genomic_state <- function() { .Deprecated('makeGenomicState')}

#' @export
make_models <- function() { .Deprecated('makeModels')}

#' @export
merge_results <- function() { .Deprecated('mergeResults')}

#' @export
preprocess_coverage <- function() { .Deprecated('preprocessCoverage')}

#' @export
raw_files <- function() { .Deprecated('rawFiles')}

#' @export
region_matrix <- function() { .Deprecated('regionMatrix')}

#' @export
sample_depth <- function() { .Deprecated('sampleDepth')}
