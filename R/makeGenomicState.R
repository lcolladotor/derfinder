#' Obtain the genomic state per region from annotation
#'
#' This function summarizes the annotation contained in a \link[GenomicFeatures]{TranscriptDb} at each given base of the genome based on annotated transcripts. It groups contiguous base pairs classified as the same type into regions.
#' 
#' @param txdb A \link[GenomicFeatures]{TranscriptDb} object.
#' @param chrs The names of the chromosomes to use as denoted in the \code{txdb} object. Check \link[GenomicFeatures]{isActiveSeq}.
#' @param addChrPrefix If \code{TRUE}, 'chr' is added as a prefix to the chromosome names (seqlevels).
#'
#' @return A \code{GRangesList} object with two elements: \code{fullGenome} and \code{codingGenome}. Both have metadata information for the type of region (theRegion), transcript IDs (tx_id), transcript name (tx_name), and gene ID (gene_id). \code{fullGenome} classifies each region as either being exon, intron or intragenic. \code{codingGenome} classfies the regions as being promoter, exon, intro, 5UTR, 3UTR or intragenic.
#'
#' @author Andrew Jaffe, Leonardo Collado-Torres
#' @seealso \link[GenomicFeatures]{TranscriptDb}
#' @export
#' @importFrom GenomicFeatures isActiveSeq "isActiveSeq<-" intronsByTranscript fiveUTRsByTranscript threeUTRsByTranscript exonsBy
#' @importFrom IRanges CharacterList elementLengths DataFrame IntegerList queryHits subjectHits Rle
#' @importFrom GenomicRanges GRangesList seqnames seqlengths seqlevels "seqlevels<-"
#' @importMethodsFrom AnnotationDbi select
#' @importMethodsFrom GenomicRanges names "names<-" reduce mcols "mcols<-" "$" "$<-" "[" "[<-" values "values<-" sort disjoin length findOverlaps split strand "strand<-" gaps width
#' @importMethodsFrom IRanges names unlist relist length lapply sapply as.character
#' @importMethodsFrom GenomicFeatures promoters
#'
#' @examples
#' \dontrun{
#' ## Hsapiens.UCSC.hg19.knownGene GenomicState
#' library("TxDb.Hsapiens.UCSC.hg19.knownGene")
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#'
#' ## Creating this GenomicState object takes around 8 min for all chrs and around 30 secs for chr21
#' GenomicState.Hsapiens.UCSC.hg19.knownGene.chr21 <- makeGenomicState(txdb=txdb, chrs="chr21")
#'
#' ## Hsapiens ENSEMBL GRCh37
#' library("GenomicFeatures")
#' ## Can take several minutes and speed will depend on your internet speed
#' xx <- makeTxDbPackageFromBiomart(version = "0.99", maintainer = "Your Name", author="Your Name")
#' txdb <- loadDb(file.path("TxDb.Hsapiens.BioMart.ensembl.GRCh37.p11", "inst", "extdata", "TxDb.Hsapiens.BioMart.ensembl.GRCh37.p11.sqlite"))
#' 
#' ## Creating this GenomicState object takes around 13 min
#' GenomicState.Hsapiens.ensembl.GRCh37.p11 <- makeGenomicState(txdb=txdb, chrs=c(1:22, "X", "Y"), addChrPrefix=TRUE)
#' 
#' ## Save for later use
#' save(GenomicState.Hsapiens.ensembl.GRCh37.p11, file="GenomicState.Hsapiens.ensembl.GRCh37.p11.Rdata")
#' }

makeGenomicState <- function(txdb, chrs=paste0("chr", c(1:22, "X", "Y")), addChrPrefix=FALSE) {
	## Select chrs to use
	isActiveSeq(txdb) <- names(isActiveSeq(txdb)) %in% chrs
	
	#########
	### add introns
	introns <- intronsByTranscript(txdb)

	## Get gene_id and tx_name associated with tx_id
	map <- select(txdb, keys=names(introns), keytype="TXID", columns=c("TXID", "TXNAME", "GENEID")) # for both R and R-devel

	## Subset on transcripts that have a GENEID
	map <- map[!is.na(map$GENEID),]

	## Group introns by gene
	tx_names <- CharacterList(split(map$TXNAME, map$GENEID))
	tx_id <- CharacterList(split(map$TXID, map$GENEID))
	tx2gene <- data.frame(tx=unlist(unname(tx_id)), gene=rep(names(tx_id), elementLengths(tx_id)))
	txidRep <- rep(as.numeric(names(introns)), elementLengths(introns))
	txnameRep <- rep(names(intronsByTranscript(txdb,use.names=TRUE)), elementLengths(introns))
	geneRep <- tx2gene$gene[match(txidRep, tx2gene$tx)]


	#### reduce introns
	intronsRed<- reduce(introns@unlistData, with.revmap=TRUE)
	introns2 <- introns@unlistData
	mcols(introns2) <- DataFrame(txidRep, txnameRep, geneRep)

	mapping <- mcols(intronsRed)$revmap
	addTx <- relist(mcols(introns2)[unlist(mapping), "txidRep"], mapping)
	addTxName <- relist(mcols(introns2)[unlist(mapping), "txnameRep"], mapping)
	addGene <- relist(mcols(introns2)[unlist(mapping), "geneRep"], mapping)
	iRegion <- rep("intronic", length(addTx))
	values(intronsRed) <- DataFrame(iRegion, addTx, addTxName, addGene)
	names(mcols(intronsRed)) <- c("region", "TxID", "TxName", "Gene")

	## promoters next
	promoters <- promoters(txdb)
	mcols(promoters)$Gene <- tx2gene$gene[match(mcols(promoters)$tx_id, tx2gene$tx)]

	promotersRed <- reduce(promoters, with.revmap=TRUE)
	mapping <- mcols(promotersRed)$revmap

	addTx <- relist(mcols(promoters)[unlist(mapping), "tx_id"], mapping)
	addTxName <- relist(mcols(promoters)[unlist(mapping), "tx_name"], mapping)
	addGene <- relist(mcols(promoters)[unlist(mapping), "Gene"], mapping)
	pRegion <- rep("promoter", length(addTx))

	values(promotersRed) <- DataFrame(pRegion, addTx, addTxName, addGene)
	names(mcols(promotersRed)) <- c("region", "TxID", "TxName", "Gene")

	##### UTRs
	# 5'
	utr5 <- fiveUTRsByTranscript(txdb)
	txidRep5 <- rep(as.numeric(names(utr5)), elementLengths(utr5))
	geneRep5 <- tx2gene$gene[match(txidRep5, tx2gene$tx)]
	txnameRep5<- rep(names(fiveUTRsByTranscript(txdb, use.names=TRUE)), elementLengths(utr5))
	utr5a <- utr5@unlistData
	mcols(utr5a) <- DataFrame(txidRep5, txnameRep5, geneRep5)

	utr5red <- reduce(utr5a, with.revmap=TRUE)

	mapping <- mcols(utr5red)$revmap
	addTx <- relist(mcols(utr5a)[unlist(mapping), "txidRep5"], mapping)
	addTxName <- relist(mcols(utr5a)[unlist(mapping), "txnameRep5"], mapping)
	addGene <- relist(mcols(utr5a)[unlist(mapping), "geneRep5"], mapping)
	uRegion <- rep("5UTR", length(addTx))
	values(utr5red) <- DataFrame(uRegion, addTx, addTxName, addGene)
	names(mcols(utr5red)) <- c("region", "TxID", "TxName", "Gene")

	# 3'
	utr3 <- threeUTRsByTranscript(txdb)
	txidRep3 <- rep(as.numeric(names(utr3)), elementLengths(utr3))
	geneRep3 <- tx2gene$gene[match(txidRep3, tx2gene$tx)]
	txnameRep3<- rep(names(threeUTRsByTranscript(txdb, use.names=TRUE)), elementLengths(utr3))
	utr3a <- utr3@unlistData
	mcols(utr3a) <- DataFrame(txidRep3, txnameRep3, geneRep3)

	utr3red <- reduce(utr3a, with.revmap=TRUE)

	mapping <- mcols(utr3red)$revmap
	addTx <- relist(mcols(utr3a)[unlist(mapping), "txidRep3"], mapping)
	addTxName <- relist(mcols(utr3a)[unlist(mapping), "txnameRep3"], mapping)
	addGene <- relist(mcols(utr3a)[unlist(mapping), "geneRep3"], mapping)
	uRegion <- rep("3UTR", length(addTx))
	values(utr3red) <- DataFrame(uRegion, addTx, addTxName, addGene)
	names(mcols(utr3red)) <- c("region", "TxID", "TxName", "Gene")

	### exons
	exons <- exonsBy(txdb)
	txidRep <- rep(as.numeric(names(exons)), elementLengths(exons))
	geneRep <- tx2gene$gene[match(txidRep, tx2gene$tx)]
	txnameRep<- rep(names(exonsBy(txdb, use.names=TRUE)), elementLengths(exons))
	exonsa <- exons@unlistData
	mcols(exonsa) <- DataFrame(txidRep, txnameRep, geneRep)

	exonsRed <- reduce(exonsa, with.revmap=TRUE)

	mapping <- mcols(exonsRed)$revmap
	addTx <- relist(mcols(exonsa)[unlist(mapping), "txidRep"], mapping)
	addTxName <- relist(mcols(exonsa)[unlist(mapping), "txnameRep"], mapping)
	addGene <- relist(mcols(exonsa)[unlist(mapping), "geneRep"], mapping)
	eRegion <- rep("exon", length(addTx))
	values(exonsRed) <- DataFrame(eRegion, addTx, addTxName, addGene)
	names(mcols(exonsRed)) <- c("region","TxID","TxName", "Gene")

	############ merge
	codingGR <- c(intronsRed, exonsRed, promotersRed, utr3red, utr5red)
	codingGR <- sort(codingGR)

	######## disjoin
	dGR <- disjoin(codingGR)

	theTx <- theGene <- IntegerList(vector("list", length(dGR)))
	theTxName <- CharacterList(vector("list", length(dGR)))
	theRegion <- rep(NA, length(dGR))

	## introns
	ooIntrons <- findOverlaps(dGR, intronsRed, type="within")
	theRegion[unique(queryHits(ooIntrons))]<- "intron"
	theTx[unique(queryHits(ooIntrons))] <- intronsRed$TxID[subjectHits(ooIntrons)]
	theTxName[unique(queryHits(ooIntrons))] <- intronsRed$TxName[subjectHits(ooIntrons)]
	theGene[unique(queryHits(ooIntrons))] <- intronsRed$Gene[subjectHits(ooIntrons)]

	## promoters
	ooPromoters <- findOverlaps(dGR, promotersRed, type="within")
	theRegion[unique(queryHits(ooPromoters))]<- "promoter" 
	theTx[unique(queryHits(ooPromoters))] <- promotersRed$TxID[subjectHits(ooPromoters)]
	theTxName[unique(queryHits(ooPromoters))] <- promotersRed$TxName[subjectHits(ooPromoters)]
	theGene[unique(queryHits(ooPromoters))] <- promotersRed$Gene[subjectHits(ooPromoters)]

	## Exons
	ooExons <- findOverlaps(dGR, exonsRed, type="within")
	theRegion[unique(queryHits(ooExons))]<- "exon" 
	theTx[unique(queryHits(ooExons))] <- exonsRed$TxID[subjectHits(ooExons)]
	theTxName[unique(queryHits(ooExons))] <- exonsRed$TxName[subjectHits(ooExons)]
	theGene[unique(queryHits(ooExons))] <- exonsRed$Gene[subjectHits(ooExons)]


	## 5' UTRs
	oo5UTR <- findOverlaps(dGR, utr5red, type="within")
	theRegion[unique(queryHits(oo5UTR))]<- "5UTR" 
	theTx[unique(queryHits(oo5UTR))] <- utr5red$TxID[subjectHits(oo5UTR)]
	theTxName[unique(queryHits(oo5UTR))] <- utr5red$TxName[subjectHits(oo5UTR)]
	theGene[unique(queryHits(oo5UTR))] <- utr5red$Gene[subjectHits(oo5UTR)]

	## 3' UTRs
	oo3UTR <- findOverlaps(dGR, utr3red, type="within")
	theRegion[unique(queryHits(oo3UTR))] <- "3UTR" 
	theTx[unique(queryHits(oo3UTR))] <- utr3red$TxID[subjectHits(oo3UTR)]
	theTxName[unique(queryHits(oo3UTR))] <- utr3red$TxName[subjectHits(oo3UTR)]
	theGene[unique(queryHits(oo3UTR))] <- utr3red$Gene[subjectHits(oo3UTR)]

	## clean up vectors
	theTx <- IntegerList(lapply(theTx, function(x) sort(unique(x))))
	theTxName <- CharacterList(lapply(theTxName, function(x) sort(unique(x))))
	theGene <- IntegerList(lapply(theGene, function(x) sort(unique(x))))


	values(dGR) <- DataFrame(theRegion, theTx, theTxName, theGene)
	names(values(dGR)) <- c("theRegion", "tx_id", "tx_name", "gene")
	codingGR <- dGR

	########
	### introns, exons, promoters
	fullGR <- c(intronsRed, exonsRed)
	######## disjoin
	dGR <- disjoin(fullGR)

	theTx <- theGene <- IntegerList(vector("list", length(dGR)))
	theTxName <- CharacterList(vector("list", length(dGR)))
	theRegion <- rep(NA, length(dGR))

	## introns
	ooIntrons <- findOverlaps(dGR, intronsRed, type="within")
	theRegion[unique(queryHits(ooIntrons))]<- "intron" 
	theTx[unique(queryHits(ooIntrons))] <- intronsRed$TxID[subjectHits(ooIntrons)]
	theTxName[unique(queryHits(ooIntrons))] <- intronsRed$TxName[subjectHits(ooIntrons)]
	theGene[unique(queryHits(ooIntrons))] <- intronsRed$Gene[subjectHits(ooIntrons)]

	## Exons
	ooExons <- findOverlaps(dGR, exonsRed, type="within")
	theRegion[unique(queryHits(ooExons))] <- "exon" 
	theTx[unique(queryHits(ooExons))] <- exonsRed$TxID[subjectHits(ooExons)]
	theTxName[unique(queryHits(ooExons))] <- exonsRed$TxName[subjectHits(ooExons)]
	theGene[unique(queryHits(ooExons))] <- exonsRed$Gene[subjectHits(ooExons)]


	## clean up vectors
	theTx <- IntegerList(lapply(theTx, function(x) sort(unique(x))))
	theTxName <- CharacterList(lapply(theTxName, function(x) sort(unique(x))))
	theGene <- IntegerList(lapply(theGene, function(x) sort(unique(x))))

	values(dGR) <- DataFrame(theRegion, theTx,theTxName,theGene)
	names(values(dGR)) <- c("theRegion", "tx_id", "tx_name", "gene")
	fullGR <- dGR

	######## reduce segments of the same type
	rIndexes <- split(seq_len(length(fullGR)), fullGR$theRegion)
	grList <- lapply(rIndexes, function(i) {
		r <- fullGR[i]
		rr <- reduce(r, with.revmap=TRUE)
	
		mapping <- sapply(mcols(rr)$mapping, "[", 1)
		values(rr) <- mcols(r)[mapping,]
		return(rr)
	})
	fullGR <- unlist(GRangesList(grList))
	fullGR <- sort(fullGR)

	rIndexes2 <- split(seq_len(length(codingGR)), codingGR$theRegion)
	grList2 <- lapply(rIndexes2, function(i) {
		r <- codingGR[i]
		rr <- reduce(r, with.revmap=TRUE)
	
		mapping <- sapply(mcols(rr)$mapping, "[", 1)
		values(rr) <- mcols(r)[mapping,]
		return(rr)
	})
	codingGR <- unlist(GRangesList(grList2))
	codingGR <- sort(codingGR)

	#### gaps to get intragenic candidate regions

	## full
	mcols(fullGR)$Gene_Strand <- strand(fullGR)
	strand(fullGR) <- Rle("*")

	fullIntraGR <- gaps(fullGR)
	strand(fullGR) <- as.character(mcols(fullGR)$Gene_Strand)
	mcols(fullGR)$Gene_Strand <- NULL

	fList <- split(fullIntraGR, seqnames(fullIntraGR))
	for(i in seq(along = fList)) {
		x <- fList[[i]]
		fList[[i]] <- x[width(x) != seqlengths(fullIntraGR)[i]]
	}
	fullIntraGR <- unlist(fList)

	addTxName <- CharacterList(vector("list", length(fullIntraGR)))
	addTx <- addGene <- IntegerList(vector("list", length(fullIntraGR)))
	gRegion <- rep("intragenic", length(fullIntraGR))
	values(fullIntraGR) <- DataFrame(gRegion, addTx, addTxName, addGene)
	names(mcols(fullIntraGR)) <- names(mcols(fullGR))


	fullGenome <- c(fullIntraGR,fullGR)
	fullGenome<- sort(fullGenome) 
	names(fullGenome) <- seq(along = fullGenome)

	## coding
	mcols(codingGR)$Gene_Strand <- strand(codingGR)
	strand(codingGR) <- Rle("*")

	codingIntraGR <- gaps(codingGR)
	strand(codingGR) <- as.character(mcols(codingGR)$Gene_Strand)
	mcols(codingGR)$Gene_Strand <- NULL

	fList <- split(codingIntraGR, seqnames(codingIntraGR))
	for(i in seq(along = fList)) {
		x <- fList[[i]]
		fList[[i]] <- x[width(x) != seqlengths(codingIntraGR)[i]]
	}
	codingIntraGR <- unlist(fList)

	addTxName <- CharacterList(vector("list", length(codingIntraGR)))
	addTx <- addGene <- IntegerList(vector("list", length(codingIntraGR)))
	gRegion <- rep("intragenic", length(codingIntraGR))
	values(codingIntraGR) <- DataFrame(gRegion,	addTx, addTxName, addGene)
	names(mcols(codingIntraGR)) <- names(mcols(codingGR))


	codingGenome <- c(codingIntraGR, codingGR)
	codingGenome <- sort(codingGenome) 
	names(codingGenome) <- seq(along = codingGenome)
	
	## Add chr prefix
	if(addChrPrefix) {
		seqlevels(codingGenome) <- paste0("chr", seqlevels(codingGenome))
		seqlevels(fullGenome) <- paste0("chr", seqlevels(fullGenome))
	}

	GenomicState <- GRangesList(fullGenome = fullGenome, codingGenome = codingGenome)
	
	## Done
	return(GenomicState)	
}
