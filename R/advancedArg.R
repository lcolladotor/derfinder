#' List advanced arguments
#'
#' Find in GitHub the documentation for the advanced arguments of a given 
#' function.
#'
#' @param fun The name of a function(s) that has advanced arguments in
#' \code{package}.
#' @param package The name of the package where the function is stored. Only
#' 'derfinder', 'derfinderPlot', and 'regionReport' are accepted.
#' @param browse Whether to open the URLs in a browser.
#'
#' @author Leonardo Collado-Torres
#' @export
#'
#' @details If you are interested on the default options used for functions
#' that run on multiple cores, check 
#' https://github.com/lcolladotor/derfinder/blob/master/R/utils.R
#' Note that in general, \link[BiocParallel]{SnowParam} is more memory efficient
#' than link[BiocParallel]{MulticoreParam}. If you so desire, use your favorite
#' cluster type by specifying \code{BPPARAM}.
#'
#' 
#' @return A vector of URLs with the GitHub search queries.
#'
#' @examples
#' ## Open the advanced argument docs for loadCoverage()
#' if(interactive()) {
#'     advancedArg('loadCoverage')
#' } else {
#'     (advancedArg('loadCoverage', browse = FALSE))
#' }

advancedArg <- function(fun, package = 'derfinder', browse = interactive()) {
    stopifnot(package %in% c('derfinder', 'derfinderPlot', 'regionReport'))
    
    query_map <- data.frame(
        'fun' = c(
            'analyzeChr', 'analyze_chr',
            'annotateRegions', 'annotate_regions',
            'calculatePvalues', 'calculate_pvalues',
            'calculateStats', 'calculate_stats',
            'coerceGR', 'coerce_gr',
            'collapseFullCoverage', 'collapse_full_coverage',
            'coverageToExon', 'coverage_to_exon',
            'createBwSample', 'create_bw_sample',
            'filterData', 'filter_data',
            'findRegions', 'find_regions',
            'fullCoverage', 'full_coverage',
            'getRegionCoverage', 'get_region_coverage',
            'loadCoverage', 'load_coverage',
            'makeGenomicState', 'make_genomic_state',
            'mergeResults', 'merge_results',
            'preprocessCoverage', 'preprocess_coverage',
            'regionMatrix', 'region_matrix',
            'sampleDepth', 'sample_depth',
            'plotCluster', 'plot_cluster',
            'plotOverview', 'plot_overview',
            'derfinderReport', 'derfinder_report'
        ),
        'query' = rep(c(
            'analyze',
            'annotate',
            'pvalues',
            'stats',
            'coerce',
            'collapse',
            'exon',
            'createBw',
            'filter',
            'find',
            'full,load',
            'get',
            'load',
            'state',
            'merge',
            'preprocess',
            'matrix',
            'depth',
            'cluster',
            'overview',
            'derfinder'
        ), each = 2),
        'repo' = rep(c('derfinder', 'derfinderPlot', 'regionReport'), c(36, 4,
            2)),
        stringsAsFactors = FALSE
    )
    repo <- NULL
    query_map <- subset(query_map, repo == package)
    
    ## Find function
    i <- which(query_map$fun == fun)
    if(length(i) == 0) {
        stop(paste0("Invalid option. 'fun' should match one of '", paste(query_map$fun, collapse="', '"), "'."))
    }
    
    ## Build url
    url <- paste0(
    'https://github.com/search?utf8=%E2%9C%93&q=advanced_argument+filename%3A',
        query_map$query[i], '+repo%3Alcolladotor%2F', query_map$repo[i],
        '+path%3A%2FR&type=Code')
    
    if(browse) {
        for(u in seq_along(url)) {
            browseURL(url[u])
        }
    }    
    return(invisible(url))
}
