#' Extract coverage information for exons
#'
#' This function extracts the coverage information calculated by \link{fullCoverage} for a set of exons determined by \link{makeGenomicState}. The underlying code is similar to \link{getRegionCoverage} with additional tweaks for calculating RPKM values.
#' 
#' @param fullCov A list where each element is the result from \link{loadCoverage} used with \code{cutoff=NULL}. The elements of the list should be named according to the chromosome number. Can be generated using \link{fullCoverage}.
#' @param genomicState The output from \link{makeGenomicState}.
#' @param fullOrCoding If \code{full} then the \code{genomicState$fullGenome} genomic state information is used. If \code{coding}, then the \code{genomicState$codingGenome} genomic state information is used.
#' @param L The width of the reads used.
#' @param returnType If \code{raw}, then the raw coverage information per exon is returned. If \code{rpkm}, RPKM values are calculated for each exon.
#' @param verbose If \code{TRUE} basic status updates will be printed along the way.
#'
#' @return A matrix (nrow = number of exons in \code{genomicState} corresponding to the chromosomes in \code{fullCov}, ncol = number of samples) with the number of reads (or RPKM) per exon. The row names correspond to the row indexes of \code{genomicState$fullGenome}  (if \code{fullOrCoding="full"}) or \code{genomicState$codingGenome} (if \code{fullOrCoding="coding"}).
#'
#' @author Andrew Jaffe, Leonardo Collado-Torres
#' @seealso \link{fullCoverage}, \link{getRegionCoverage}
#' @export
#' @importFrom GenomicRanges seqlevels seqnames
#' @importMethodsFrom GenomicRanges names "names<-" length "[" coverage sort width c strand "%in%"
#' @importMethodsFrom IRanges subset as.data.frame as.character runValue "%in%"
#'
#' @examples
#' \dontrun{
#' ## Calculate library size adjustments
#' sampleDepths <- sampleDepth(list(genomeData$coverage), prob=0.5, nonzero=TRUE, center=TRUE, verbose=TRUE)
#' sampleDepths
#' 
#' ## Build the models
#' group <- genomeInfo$pop
#' adjustvars <- data.frame(genomeInfo$gender)
#' models <- makeModels(sampleDepths, testvars=group, adjustvars=adjustvars)
#'
#' ## Preprocess the data
#' ## Automatic chunksize used to then compare 1 vs 4 cores in the 'do not run' section
#' prep <- preprocessCoverage(genomeData, groupInfo=group, cutoff=0, scalefac=32, chunksize=NULL, colsubset=NULL, mc.cores=4)
#' 
#' ## Get the F statistics
#' fstats <- calculateStats(prep, models, mc.cores=1, verbose=TRUE)
#'
#' ## Determine a cutoff from the F-distribution.
#' ## This step is very important and you should consider using quantiles from the observed F statistics
#' n <- dim(prep$coverageSplit[[1]])[2]
#' df1 <- dim(models$mod)[2]
#' df0 <- dim(models$mod0)[2]
#' cutoff <- qf(0.95, df1-df0, n-df1)
#' 
#' ## Low cutoff used for illustrative purposes
#' cutoff <- 1
#'
#' ## Calculate the p-values and define the regions of interest.
#' regsWithP <- calculatePvalues(prep, models, fstats, nPermute=10, seeds=NULL, chr="chr21", cutoff=cutoff, mc.cores=1)
#'
#' ## Obtain fullCov object
#' datadir <- system.file("extdata", "genomeData", package="derfinder")
#' dirs <- makeBamList(datadir=datadir, samplepatt="*accepted_hits.bam$", bamterm=NULL)
#' ## Shorten the column names
#' names(dirs) <- gsub("_accepted_hits.bam", "", names(dirs))
#' 
#' ## Reading the data and filtering it is quite fast.
#' fullCov <- fullCoverage(dirs=dirs, chrnums="21", mc.cores=1)
#' 
#' ## Create GenomicState object:
#' ## Hsapiens.UCSC.hg19.knownGene GenomicState
#' library("TxDb.Hsapiens.UCSC.hg19.knownGene")
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#'
#' ## Creating this GenomicState object takes around 8 min
#' GenomicState.Hsapiens.UCSC.hg19.knownGene <- makeGenomicState(txdb=txdb)
#' 
#' ## Annotate regions
#' annotatedRegions <- annotateRegions(regions=regsWithP$regions, genomicState=GenomicState.Hsapiens.UCSC.hg19.knownGene, minoverlap=1)
#' 
#' ## Finally, get the coverage information for each exon
#' exonCov <- coverageToExon(fullCov=fullCov, genomicState=GenomicState.Hsapiens.UCSC.hg19.knownGene)
#' }


coverageToExon <- function(fullCov, genomicState, fullOrCoding = "full", L = 100, returnType = "raw", verbose=TRUE) {
	stopifnot(length(intersect(names(genomicState), c("fullGenome", "codingGenome"))) == 2)
	stopifnot(length(intersect(fullOrCoding, c("full", "coding"))) == 1)
	stopifnot(length(intersect(returnType, c("raw", "rpkm"))) == 1)
	
	if(fullOrCoding == "full") {
		gs <- genomicState$fullGenome
	} else if (fullOrCoding == "coding") {
		gs <- genomicState$codingGenome
	}

	# just the reduced exons
	etab <- gs[gs$theRegion=="exon"]
	
	## Keep only the exons from the chromosomes in fullCov
	etab <- etab[ seqnames(etab) %in% paste0("chr", names(fullCov)) ]

	# split by strand
	strandIndexes <- split(seq_len(length(etab)), as.character(strand(etab)))
	
	# count reads covering exon on each strand
	exonByStrand <- lapply(strandIndexes, function(ii) {
		e <- etab[ii] # subset
		
		## use logical rle to subset large coverage matrix
		cc <- coverage(e) # first coverage
		for(i in seq(along=cc)) { # then convert to logical
			cc[[i]]@values <- ifelse(cc[[i]]@values > 0, TRUE, FALSE)
		}

		# now count exons
		exonList <- list() # blank list
		for(i in seq(along=fullCov)) { # by chromosome
			chrnum <- names(fullCov)[i]
			chr <- paste0("chr", chrnum)
			if(verbose) message(paste(Sys.time(), "coverageToExon: processing chromosome", chrnum))

			# subset using logical rle (fastest way)
			z <- as.data.frame(subset(fullCov[[chrnum]], cc[[chr]]))
			
			# only exons from this chr
			g <- e[seqnames(e) == chr]
			ind <- rep(names(g), width(g)) # to split
			tmpList <- split(z, ind) # split
			exonList[[i]] <- t(sapply(tmpList, colSums)/L) # get # reads
		}
		out <- do.call("rbind", exonList) # combine
	})

	# combine two strands
	exons <- do.call("rbind", exonByStrand)
	
	# put back in annotation order
	theExons <- exons[names(etab),]

	if(returnType == "rpkm") {
		Mchr <- t(sapply(fullCov, function(z) sapply(z, function(xx) sum(as.numeric(runValue(xx))))))
		M <- colSums(Mchr)/L/1e6
		theExons <- theExons/(width(etab)/1000)/M
	} 
	return(theExons)
}