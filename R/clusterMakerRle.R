#' Make clusters of genomic locations based on distance in Rle() world
#'
#' Genomic locations are grouped into clusters based on distance: locations 
#' that are close to each other are assigned to the same cluster. The operation 
#' is performed on each chromosome independently. This is very similar to 
#' \link[bumphunter]{clusterMaker}.
#'
#' @details
#' \link[bumphunter]{clusterMaker} adapted to Rle world. Assumes that the data 
#' is sorted and that everything is in a single chromosome.
#' It is also almost as fast as the original version with the advantage that 
#' everything is in Rle() world.
#' 
#' It is a a helper function for \link{findRegions}.
#' 
#' @param position A logical Rle indicating the chromosome positions.
#' @param maxGap An integer. Genomic locations within \code{maxGap} from each 
#' other are placed into the same cluster.
#' @param ranges If \code{TRUE} then an IRanges object is returned instead of 
#' the usual integer Rle.
#'
#' @return An integer Rle with the cluster IDs. If \code{ranges=TRUE} then it 
#' is an IRanges object with one range per cluster.
#'
#' @seealso \link[bumphunter]{clusterMaker}, \link{findRegions}
#' @references Rafael A. Irizarry, Martin Aryee, Hector Corrada Bravo, Kasper 
#' D. Hansen and Harris A. Jaffee. bumphunter: Bump Hunter. R package version 
#' 1.1.10.
#' @author Leonardo Collado-Torres
#' @export
#' @importFrom IRanges IRanges start end reduce Views runLength
#' @importMethodsFrom IRanges length sum
#' @importFrom S4Vectors Rle runValue
#'
#' @examples
#' library('IRanges')
#' set.seed(20130725)
#' pos <- Rle(sample(c(TRUE, FALSE), 1e5, TRUE, prob=c(0.05, 0.95)))
#' cluster <- clusterMakerRle(pos, 100L)
#' cluster
#' 
#' \dontrun{
#' ## clusterMakerRle() is comparable in speed if you start from the Rle world.
#' library('bumphunter')
#' library('microbenchmark')
#' micro <- microbenchmark(clusterMakerRle(pos, 100L), 
#'     clusterMaker(chr=rep('chr21', sum(pos)), pos=which(pos)))
#' micro
#' }

clusterMakerRle <- function(position, maxGap = 300L, ranges = FALSE) {
    ## Instead of using which(), identify the regions of the chr
    ## with data
    ir <- IRanges(start = start(position)[runValue(position)], 
        end = end(position)[runValue(position)])
    
    ## Apply the gap reduction
    ir.red <- reduce(ir, min.gapwidth = maxGap + 1)
    rm(ir)
    
    ## Identify the clusters
    clusterIDs <- Rle(seq_len(length(ir.red)), sum(Views(position, 
        ir.red)))
    ## Note that sum(Views(pos, ir.red)) is faster than
    ## sapply(ir.red, function(x) sum(pos[x]))
    
    ## Group the information into an IRanges object
    if (ranges) {
        csum <- cumsum(runLength(clusterIDs))
        result <- IRanges(start = c(1, csum[-length(csum)] + 
            1), end = csum)
    } else {
        result <- clusterIDs
    }
    
    ## Done
    return(result)
} 
