#' Obtain the genomic state per region from annotation
#'
#' This function summarizes the annotation contained in a 
#' \link[GenomicFeatures]{TxDb} at each given base of the genome based 
#' on annotated transcripts. It groups contiguous base pairs classified as the 
#' same type into regions.
#' 
#' @param txdb A \link[GenomicFeatures]{TxDb} object.
#' @param chrs The names of the chromosomes to use as denoted in the 
#' \code{txdb} object. Check \link[GenomicFeatures]{isActiveSeq}.
#' @param ... Arguments passed to other methods and/or advanced arguments.
#'
#' @return A \code{GRangesList} object with two elements: \code{fullGenome} and 
#' \code{codingGenome}. Both have metadata information for the type of region 
#' (theRegion), transcript IDs (tx_id), transcript name (tx_name), and gene ID 
#' (gene_id). \code{fullGenome} classifies each region as either being exon, 
#' intron or intergenic. \code{codingGenome} classfies the regions as being 
#' promoter, exon, intro, 5UTR, 3UTR or intergenic.
#'
#' @author Andrew Jaffe, Leonardo Collado-Torres
#' @seealso \link[GenomicFeatures]{TxDb}
#' @export
#'
#' @importFrom GenomicFeatures isActiveSeq 'isActiveSeq<-' intronsByTranscript 
#' fiveUTRsByTranscript threeUTRsByTranscript exonsBy
#' @importFrom IRanges CharacterList IntegerList
#' @import S4Vectors
#' @importFrom GenomicRanges GRangesList seqnames 
#' @importFrom GenomeInfoDb seqlengths seqlevels renameSeqlevels
#' @importMethodsFrom AnnotationDbi select
#' @importMethodsFrom GenomicRanges names 'names<-' reduce '$' 
#' '$<-' '[' '[<-' sort disjoin length findOverlaps
#' strand 'strand<-' gaps width
#' @importMethodsFrom IRanges names unlist relist length lapply '['
#' @import S4Vectors
#' @importMethodsFrom GenomicFeatures promoters
#'
#' @examples
#' ## Load the example data base from the GenomicFeatures vignette
#' library('GenomicFeatures')
#' samplefile <- system.file('extdata', 'hg19_knownGene_sample.sqlite', 
#'     package='GenomicFeatures')
#' txdb <- loadDb(samplefile)
#'
#' ## Generate genomic state object, only for chr6
#' sampleGenomicState <- makeGenomicState(txdb, chrs='chr6')
# 
#' \dontrun{
#' ## Create the GenomicState object for Hsapiens.UCSC.hg19.knownGene
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
#'
#' ## Creating this GenomicState object takes around 8 min for all chrs and 
#' ## around 30 secs for chr21
#' GenomicState.Hsapiens.UCSC.hg19.knownGene.chr21 <- 
#'     makeGenomicState(txdb=txdb, chrs='chr21')
#' 
#' ## For convinience, this object is already included in derfinder
#' library('testthat')
#' expect_that(GenomicState.Hsapiens.UCSC.hg19.knownGene.chr21, 
#'    is_equivalent_to(genomicState))
#'
#' ## Hsapiens ENSEMBL GRCh37
#' library('GenomicFeatures')
#' ## Can take several minutes and speed will depend on your internet speed
#' xx <- makeTxDbPackageFromBiomart(version = '0.99', maintainer = 'Your Name', 
#'     author='Your Name')
#' txdb <- loadDb(file.path('TxDb.Hsapiens.BioMart.ensembl.GRCh37.p11', 'inst', 
#'     'extdata', 'TxDb.Hsapiens.BioMart.ensembl.GRCh37.p11.sqlite'))
#' 
#' ## Creating this GenomicState object takes around 13 min
#' GenomicState.Hsapiens.ensembl.GRCh37.p11 <- makeGenomicState(txdb=txdb, 
#'     chrs=c(1:22, 'X', 'Y'))
#' 
#' ## Save for later use
#' save(GenomicState.Hsapiens.ensembl.GRCh37.p11, 
#'     file='GenomicState.Hsapiens.ensembl.GRCh37.p11.Rdata')
#' }

makeGenomicState <- function(txdb, chrs = c(1:22, 'X', 'Y'), ...) {
    stopifnot(is(txdb, "TxDb"))

    ## Use UCSC names for homo_sapiens by default
    chrs <- extendedMapSeqlevels(chrs, ...)

    ## Select chrs to use
    isActiveSeq(txdb) <- names(isActiveSeq(txdb)) %in% chrs
    
    ######### add introns
    introns <- intronsByTranscript(txdb)
    
    ## Get gene_id and tx_name associated with tx_id
    map <- select(txdb, keys = names(introns), keytype = 'TXID', 
        columns = c('TXID', 'TXNAME', 'GENEID'))  # for both R and R-devel
    
    ## Subset on transcripts that have a GENEID
    map <- map[!is.na(map$GENEID), ]
    
    ## Group introns by gene
    tx_names <- CharacterList(split(map$TXNAME, map$GENEID))
    tx_id <- CharacterList(split(map$TXID, map$GENEID))
    tx2gene <- data.frame(tx = unlist(unname(tx_id)), gene = rep(names(tx_id), 
        elementNROWS(tx_id)))
    txidRep <- rep(as.numeric(names(introns)), elementNROWS(introns))
    txnameRep <- rep(names(intronsByTranscript(txdb, use.names = TRUE)), 
        elementNROWS(introns))
    geneRep <- tx2gene$gene[match(txidRep, tx2gene$tx)]
    
    
    #### reduce introns
    intronsRed <- reduce(introns@unlistData, with.revmap = TRUE)
    introns2 <- introns@unlistData
    mcols(introns2) <- DataFrame(txidRep, txnameRep, geneRep,
        check.names = FALSE)
    
    mapping <- mcols(intronsRed)$revmap
    addTx <- relist(mcols(introns2)[unlist(mapping), 'txidRep'], 
        mapping)
    addTxName <- relist(mcols(introns2)[unlist(mapping), 'txnameRep'], 
        mapping)
    addGene <- relist(mcols(introns2)[unlist(mapping), 'geneRep'], 
        mapping)
    iRegion <- rep('intronic', length(addTx))
    values(intronsRed) <- DataFrame(iRegion, addTx, addTxName, 
        addGene, check.names = FALSE)
    names(mcols(intronsRed)) <- c('region', 'TxID', 'TxName', 
        'Gene')
    
    ## promoters next
    promoters <- promoters(txdb)
    mcols(promoters)$Gene <- tx2gene$gene[match(mcols(promoters)$tx_id, 
        tx2gene$tx)]
    
    promotersRed <- reduce(promoters, with.revmap = TRUE)
    mapping <- mcols(promotersRed)$revmap
    
    addTx <- relist(mcols(promoters)[unlist(mapping), 'tx_id'], 
        mapping)
    addTxName <- relist(mcols(promoters)[unlist(mapping), 'tx_name'], 
        mapping)
    addGene <- relist(mcols(promoters)[unlist(mapping), 'Gene'], 
        mapping)
    pRegion <- rep('promoter', length(addTx))
    
    values(promotersRed) <- DataFrame(pRegion, addTx, addTxName, 
        addGene, check.names = FALSE)
    names(mcols(promotersRed)) <- c('region', 'TxID', 'TxName', 
        'Gene')
    
    ##### UTRs 5'
    utr5 <- fiveUTRsByTranscript(txdb)
    txidRep5 <- rep(as.numeric(names(utr5)), elementNROWS(utr5))
    geneRep5 <- tx2gene$gene[match(txidRep5, tx2gene$tx)]
    txnameRep5 <- rep(names(fiveUTRsByTranscript(txdb, use.names = TRUE)), 
        elementNROWS(utr5))
    utr5a <- utr5@unlistData
    mcols(utr5a) <- DataFrame(txidRep5, txnameRep5, geneRep5,
        check.names = FALSE)
    
    utr5red <- reduce(utr5a, with.revmap = TRUE)
    
    mapping <- mcols(utr5red)$revmap
    addTx <- relist(mcols(utr5a)[unlist(mapping), 'txidRep5'], 
        mapping)
    addTxName <- relist(mcols(utr5a)[unlist(mapping), 'txnameRep5'], 
        mapping)
    addGene <- relist(mcols(utr5a)[unlist(mapping), 'geneRep5'], 
        mapping)
    uRegion <- rep('5UTR', length(addTx))
    values(utr5red) <- DataFrame(uRegion, addTx, addTxName, addGene,
        check.names = FALSE)
    names(mcols(utr5red)) <- c('region', 'TxID', 'TxName', 'Gene')
    
    # 3'
    utr3 <- threeUTRsByTranscript(txdb)
    txidRep3 <- rep(as.numeric(names(utr3)), elementNROWS(utr3))
    geneRep3 <- tx2gene$gene[match(txidRep3, tx2gene$tx)]
    txnameRep3 <- rep(names(threeUTRsByTranscript(txdb, use.names = TRUE)), 
        elementNROWS(utr3))
    utr3a <- utr3@unlistData
    mcols(utr3a) <- DataFrame(txidRep3, txnameRep3, geneRep3,
        check.names = FALSE)
    
    utr3red <- reduce(utr3a, with.revmap = TRUE)
    
    mapping <- mcols(utr3red)$revmap
    addTx <- relist(mcols(utr3a)[unlist(mapping), 'txidRep3'], 
        mapping)
    addTxName <- relist(mcols(utr3a)[unlist(mapping), 'txnameRep3'], 
        mapping)
    addGene <- relist(mcols(utr3a)[unlist(mapping), 'geneRep3'], 
        mapping)
    uRegion <- rep('3UTR', length(addTx))
    values(utr3red) <- DataFrame(uRegion, addTx, addTxName, addGene,
        check.names = FALSE)
    names(mcols(utr3red)) <- c('region', 'TxID', 'TxName', 'Gene')
    
    ### exons
    exons <- exonsBy(txdb)
    txidRep <- rep(as.numeric(names(exons)), elementNROWS(exons))
    geneRep <- tx2gene$gene[match(txidRep, tx2gene$tx)]
    txnameRep <- rep(names(exonsBy(txdb, use.names = TRUE)), 
        elementNROWS(exons))
    exonsa <- exons@unlistData
    mcols(exonsa) <- DataFrame(txidRep, txnameRep, geneRep, check.names = FALSE)
    
    exonsRed <- reduce(exonsa, with.revmap = TRUE)
    
    mapping <- mcols(exonsRed)$revmap
    addTx <- relist(mcols(exonsa)[unlist(mapping), 'txidRep'], 
        mapping)
    addTxName <- relist(mcols(exonsa)[unlist(mapping), 'txnameRep'], 
        mapping)
    addGene <- relist(mcols(exonsa)[unlist(mapping), 'geneRep'], 
        mapping)
    eRegion <- rep('exon', length(addTx))
    values(exonsRed) <- DataFrame(eRegion, addTx, addTxName, 
        addGene, check.names = FALSE)
    names(mcols(exonsRed)) <- c('region', 'TxID', 'TxName', 'Gene')
    
    ############ merge
    codingGR <- c(intronsRed, exonsRed, promotersRed, utr3red, 
        utr5red)
    codingGR <- sort(codingGR)
    
    ######## disjoin
    dGR <- disjoin(codingGR)
    
    theTx <- theGene <- IntegerList(vector('list', length(dGR)))
    theTxName <- CharacterList(vector('list', length(dGR)))
    theRegion <- rep(NA, length(dGR))
    
    ## introns
    ooIntrons <- findOverlaps(dGR, intronsRed, type = 'within')
    theRegion[unique(queryHits(ooIntrons))] <- 'intron'
    theTx[unique(queryHits(ooIntrons))] <-
        intronsRed$TxID[subjectHits(ooIntrons)]
    theTxName[unique(queryHits(ooIntrons))] <-
        intronsRed$TxName[subjectHits(ooIntrons)]
    theGene[unique(queryHits(ooIntrons))] <-
        intronsRed$Gene[subjectHits(ooIntrons)]
    
    ## promoters
    ooPromoters <- findOverlaps(dGR, promotersRed, type = 'within')
    theRegion[unique(queryHits(ooPromoters))] <- 'promoter'
    theTx[unique(queryHits(ooPromoters))] <-
        promotersRed$TxID[subjectHits(ooPromoters)]
    theTxName[unique(queryHits(ooPromoters))] <-
        promotersRed$TxName[subjectHits(ooPromoters)]
    theGene[unique(queryHits(ooPromoters))] <-
        promotersRed$Gene[subjectHits(ooPromoters)]
    
    ## Exons
    ooExons <- findOverlaps(dGR, exonsRed, type = 'within')
    theRegion[unique(queryHits(ooExons))] <- 'exon'
    theTx[unique(queryHits(ooExons))] <- exonsRed$TxID[subjectHits(ooExons)]
    theTxName[unique(queryHits(ooExons))] <-
        exonsRed$TxName[subjectHits(ooExons)]
    theGene[unique(queryHits(ooExons))] <- exonsRed$Gene[subjectHits(ooExons)]
    
    
    ## 5' UTRs
    oo5UTR <- findOverlaps(dGR, utr5red, type = 'within')
    theRegion[unique(queryHits(oo5UTR))] <- '5UTR'
    theTx[unique(queryHits(oo5UTR))] <- utr5red$TxID[subjectHits(oo5UTR)]
    theTxName[unique(queryHits(oo5UTR))] <- utr5red$TxName[subjectHits(oo5UTR)]
    theGene[unique(queryHits(oo5UTR))] <- utr5red$Gene[subjectHits(oo5UTR)]
    
    ## 3' UTRs
    oo3UTR <- findOverlaps(dGR, utr3red, type = 'within')
    theRegion[unique(queryHits(oo3UTR))] <- '3UTR'
    theTx[unique(queryHits(oo3UTR))] <- utr3red$TxID[subjectHits(oo3UTR)]
    theTxName[unique(queryHits(oo3UTR))] <- utr3red$TxName[subjectHits(oo3UTR)]
    theGene[unique(queryHits(oo3UTR))] <- utr3red$Gene[subjectHits(oo3UTR)]
    
    ## clean up vectors
    theTx <- IntegerList(lapply(theTx, function(x) sort(unique(x))))
    theTxName <- CharacterList(lapply(theTxName, function(x) sort(unique(x))))
    theGene <- IntegerList(lapply(theGene, function(x) sort(unique(x))))
    
    
    values(dGR) <- DataFrame(theRegion, theTx, theTxName, theGene,
        check.names = FALSE)
    names(values(dGR)) <- c('theRegion', 'tx_id', 'tx_name', 
        'gene')
    codingGR <- dGR
    
    ######## introns, exons, promoters
    fullGR <- c(intronsRed, exonsRed)
    ######## disjoin
    dGR <- disjoin(fullGR)
    
    theTx <- theGene <- IntegerList(vector('list', length(dGR)))
    theTxName <- CharacterList(vector('list', length(dGR)))
    theRegion <- rep(NA, length(dGR))
    
    ## introns
    ooIntrons <- findOverlaps(dGR, intronsRed, type = 'within')
    theRegion[unique(queryHits(ooIntrons))] <- 'intron'
    theTx[unique(queryHits(ooIntrons))] <-
        intronsRed$TxID[subjectHits(ooIntrons)]
    theTxName[unique(queryHits(ooIntrons))] <-
        intronsRed$TxName[subjectHits(ooIntrons)]
    theGene[unique(queryHits(ooIntrons))] <-
        intronsRed$Gene[subjectHits(ooIntrons)]
    
    ## Exons
    ooExons <- findOverlaps(dGR, exonsRed, type = 'within')
    theRegion[unique(queryHits(ooExons))] <- 'exon'
    theTx[unique(queryHits(ooExons))] <- exonsRed$TxID[subjectHits(ooExons)]
    theTxName[unique(queryHits(ooExons))] <-
        exonsRed$TxName[subjectHits(ooExons)]
    theGene[unique(queryHits(ooExons))] <- exonsRed$Gene[subjectHits(ooExons)]
    
    
    ## clean up vectors
    theTx <- IntegerList(lapply(theTx, function(x) sort(unique(x))))
    theTxName <- CharacterList(lapply(theTxName, function(x) sort(unique(x))))
    theGene <- IntegerList(lapply(theGene, function(x) sort(unique(x))))
    
    values(dGR) <- DataFrame(theRegion, theTx, theTxName, theGene,
        check.names = FALSE)
    names(values(dGR)) <- c('theRegion', 'tx_id', 'tx_name', 
        'gene')
    fullGR <- dGR
    
    ######## reduce segments of the same type
    rIndexes <- split(seq_len(length(fullGR)), fullGR$theRegion)
    grList <- lapply(rIndexes, function(i) {
        r <- fullGR[i]
        rr <- reduce(r, with.revmap = TRUE)
        
        mapping <- sapply(mcols(rr)$revmap, '[', 1)
        values(rr) <- mcols(r)[mapping, ]
        return(rr)
    })
    fullGR <- unlist(GRangesList(grList))
    fullGR <- sort(fullGR)
    
    rIndexes2 <- split(seq_len(length(codingGR)), codingGR$theRegion)
    grList2 <- lapply(rIndexes2, function(i) {
        r <- codingGR[i]
        rr <- reduce(r, with.revmap = TRUE)
        
        mapping <- sapply(mcols(rr)$revmap, '[', 1)
        values(rr) <- mcols(r)[mapping, ]
        return(rr)
    })
    codingGR <- unlist(GRangesList(grList2))
    codingGR <- sort(codingGR)
    
    #### gaps to get intergenic candidate regions
    
    ## full
    mcols(fullGR)$Gene_Strand <- strand(fullGR)
    strand(fullGR) <- Rle('*')
    
    fullInterGR <- gaps(fullGR)
    strand(fullGR) <- as.character(mcols(fullGR)$Gene_Strand)
    mcols(fullGR)$Gene_Strand <- NULL
    
    fList <- split(fullInterGR, seqnames(fullInterGR))
    for (i in seq(along = fList)) {
        x <- fList[[i]]
        fList[[i]] <- x[width(x) != seqlengths(fullInterGR)[i]]
    }
    fullInterGR <- unlist(fList)
    
    addTxName <- CharacterList(vector('list', length(fullInterGR)))
    addTx <- addGene <- IntegerList(vector('list', length(fullInterGR)))
    gRegion <- rep('intergenic', length(fullInterGR))
    values(fullInterGR) <- DataFrame(gRegion, addTx, addTxName, 
        addGene, check.names = FALSE)
    names(mcols(fullInterGR)) <- names(mcols(fullGR))
    
    
    fullGenome <- c(fullInterGR, fullGR)
    fullGenome <- sort(fullGenome)
    names(fullGenome) <- seq(along = fullGenome)
    
    ## coding
    mcols(codingGR)$Gene_Strand <- strand(codingGR)
    strand(codingGR) <- Rle('*')
    
    codingInterGR <- gaps(codingGR)
    strand(codingGR) <- as.character(mcols(codingGR)$Gene_Strand)
    mcols(codingGR)$Gene_Strand <- NULL
    
    fList <- split(codingInterGR, seqnames(codingInterGR))
    for (i in seq(along = fList)) {
        x <- fList[[i]]
        fList[[i]] <- x[width(x) != seqlengths(codingInterGR)[i]]
    }
    codingInterGR <- unlist(fList)
    
    addTxName <- CharacterList(vector('list', length(codingInterGR)))
    addTx <- addGene <- IntegerList(vector('list', length(codingInterGR)))
    gRegion <- rep('intergenic', length(codingInterGR))
    values(codingInterGR) <- DataFrame(gRegion, addTx, addTxName, 
        addGene, check.names = FALSE)
    names(mcols(codingInterGR)) <- names(mcols(codingGR))
    
    
    codingGenome <- c(codingInterGR, codingGR)
    codingGenome <- sort(codingGenome)
    names(codingGenome) <- seq(along = codingGenome)
    
    ## Use UCSC names for homo_sapiens by default
    codingGenome <- renameSeqlevels(codingGenome,
        extendedMapSeqlevels(seqlevels(codingGenome), ...))
    fullGenome <- renameSeqlevels(fullGenome,
        extendedMapSeqlevels(seqlevels(fullGenome), ...))
    
    GenomicState <- GRangesList(fullGenome = fullGenome,
        codingGenome = codingGenome)
    
    ## Done
    return(GenomicState)
}
