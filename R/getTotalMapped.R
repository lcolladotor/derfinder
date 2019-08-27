#' Calculate the total number of mapped reads
#'
#' For a given BAM calculate the total number of mapped reads and for a BigWig
#' file calculate the area under the curve (AUC), which is related to the
#' number of mapped reads: the exact relationship depends on whether the
#' aligner soft clips reads and/or if the length of the reads is the same.
#' If you use the 'chrs' argument you can choose to only consider the
#' information for your chromosomes of interest. For example, you can exclude
#' the mitocondrial chromosome.
#'
#' @param rawFile Either a BAM file or a BigWig file.
#' @param chrs If `NULL`, all the chromosomes will be used. Otherwise,
#' only those in `chrs` will be used.
#'
#' @return The total number of mapped reads for a BAM file or the AUC for a
#' BigWig file in a single element vector.
#'
#' @author Leonardo Collado-Torres
#'
#' @export
#' @importFrom rtracklayer BigWigFileList BigWigFile
#' @importMethodsFrom rtracklayer import import.bw
#' @importFrom Rsamtools idxstatsBam BamFileList
#'
#' @examples
#'
#' ## Get the total number of mapped reads for an example BAM file:
#' bam <- system.file('extdata', 'genomeData', 'ERR009102_accepted_hits.bam',
#'     package='derfinder', mustWork=TRUE)
#' getTotalMapped(bam)
#'

getTotalMapped <- function(rawFile, chrs = NULL) {
    stopifnot(length(rawFile) == 1)
    
    ## Guess the input type
    if(is(rawFile, 'BigWigFileList') | is(rawFile, 'BigWigFile') | grepl('bw$|bigwig$', tolower(rawFile))) {
        inputType <- 'BigWig'
    } else if (is(rawFile, 'BamFileList') | is(rawFile, 'BamFile') | grepl('bam$', tolower(rawFile))) {
        inputType <- 'bam'
    }
    stopifnot(inputType %in% c('bam', 'BigWig'))
    
    if(inputType == 'bam') {
        info <- idxstatsBam(rawFile)
        if(!is.null(chrs)) {
            info <- subset(info, seqnames %in% chrs) 
        }
        res <- sum(info$mapped)
    } else if (inputType == 'BigWig') {
        if(.Platform$OS.type == 'windows') {
            warning('As of rtracklayer 1.25.16, BigWig is not supported on Windows. Thus loading data from BigWig files will most likely fail!')
        }
        gr <- import.bw(rawFile)
        
        if(!is.null(chrs)) {
            gr <- gr[seqnames(gr) %in% chrs]
        }
        
        res <- sum(width(gr) * gr$score) 
    }
    return(res)
}
