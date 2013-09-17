#' Generate a HTML report exploring the basic results
#'
#' This function generates a HTML report exploring the basic results. 
#' 
#' @param prefix The main data directory path where \link{mergeResults} was run. It should be the same as \code{mergeResults(prefix)}.
#' @param outdir The name of output directory relative to \code{prefix}.
#' @param output The name output HTML file.
#' @param project The title of the project.
#' @param browse If \code{TRUE} the HTML report is opened in your browser once it's completed.
#' @param nBestRegions The number of region plots to make, ordered by area.
#' @param makeBestClusters If \code{TRUE}, \link{plotCluster} is used on the \code{nBestClusters} regions by area. Note that these plots take some time to make.
#' @param nBestClusters The number of region cluster plots to make, ordered by area.
#' @param fullCov A list where each element is the result from \link{loadCoverage} used with \code{cutoff=NULL}. The elements of the list should be named according to the chromosome number. Can be generated using \link{fullCoverage}.
#' @param hg19 If \code{TRUE} then the reference is assumed to be hg19 and chromosome lengths as well as the default transcription database (TxDb.Hsapiens.UCSC.hg19.knownGene) will be used.
#' @param p.ideos A list where each element is the result of \link[ggbio]{plotIdeogram}. If it's \code{NULL} and \code{hg19=TRUE} then they are created for the hg19 human reference.
#' @param txdb Specify the transcription database to use for making the plots for the top regions by area. If \code{NULL} and \code{hg19=TRUE} then TxDb.Hsapiens.UCSC.hg19.knownGene is used.
#' @param installMissing If \code{TRUE} all missing required packages are installed. Note that some are development versions hosted in GitHub.
#' @param device The graphical device used when knitting. See more at http://yihui.name/knitr/options (dev argument).
#' @param fullRegions Part of the output of \link{mergeResults}. Specify it only if you have already loaded it in memory.
#' @param fullNullSummary Part of the output of \link{mergeResults}. Specify it only if you have already loaded it in memory.
#' @param fullAnnotatedRegions Part of the output of \link{mergeResults}. Specify it only if you have already loaded it in memory.
#' @param optionsStats Part of the output of \link{analyzeChr}. Specify it only if you have already loaded it in memory.
#'
#' @return An HTML report with a basic exploration of the results.
#'
#' @author Leonardo Collado-Torres
#' @seealso \link{mergeResults}, \link{analyzeChr}, \link{fullCoverage}
#' @export
#' @importMethodsFrom IRanges as.numeric
#'
#' @examples
#' \dontrun{
#' generateReport(prefix="run1", makeBestClusters=FALSE)
#' }


generateReport <- function(prefix, outdir="basicExploration", output="basicExploration.html", project=prefix, browse=interactive(), nBestRegions=100, makeBestClusters=TRUE, nBestClusters=2, fullCov=NULL, hg19=TRUE, p.ideos=NULL, txdb=NULL, installMissing=TRUE, device="CairoPNG", fullRegions=NULL, fullNullSummary=NULL, fullAnnotatedRegions=NULL, optionsStats=NULL) {

	## Save start time for getting the total processing time
	startTime <- Sys.time()
	
	## Pleasing R CMD check
	biocLite <- function(x) { NA }
	rm(biocLite)
	
	#### Load/install packages
	## BioC
	if(!suppressMessages(require("IRanges"))) {
		if(installMissing) {
			source("http://bioconductor.org/biocLite.R")
			biocLite("IRanges")
		}
		library("IRanges")
	}
	if(!suppressMessages(require("GenomicRanges"))) {
		if(installMissing) {
			source("http://bioconductor.org/biocLite.R")
			biocLite("GenomicRanges")
		}
		library("GenomicRanges")
	}
	if(hg19) {
		if(!suppressMessages(require("TxDb.Hsapiens.UCSC.hg19.knownGene"))) {
			if(installMissing) {
				source("http://bioconductor.org/biocLite.R")
				biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
			}
			library("TxDb.Hsapiens.UCSC.hg19.knownGene")
		}
	}
	
	## CRAN
	if(!suppressMessages(require("ggplot2"))) {
		if(installMissing) {
			install.packages("ggplot2")
		}
		library("ggplot2")
	}
	if(!suppressMessages(require("gridExtra"))) {
		if(installMissing) {
			install.packages("gridExtra")
		}
		library("gridExtra")
	}
	if(!suppressMessages(require("data.table"))) {
		if(installMissing) {
			install.packages("data.table")
		}
		library("data.table")
	}
	if(!suppressMessages(require("knitcitations"))) {
		if(installMissing) {
			install.packages("knitcitations")
		}
		library("knitcitations")
	}
	if(!suppressMessages(require("xtable"))) {
		if(installMissing) {
			install.packages("xtable")
		}
		library("xtable")
	}
	if(!suppressMessages(require("RColorBrewer"))) {
		if(installMissing) {
			install.packages("RColorBrewer")
		}
		library("RColorBrewer")
	}
	
	## GitHub	
	if(!suppressMessages(require("knitrBootstrap"))) {
		if(installMissing) {
			if(!require("devtools")) {
				if(installMissing) {
					install.packages("devtools")
				}				
				library("devtools")
			}
			install_github(username='jimhester', repo='knitrBootstrap')
		}
		library("knitrBootstrap")
	}
	if(!suppressMessages(require("rCharts"))) {
		if(installMissing) {
			if(!require("devtools")) {
				if(installMissing) {
					install.packages("devtools")
				}				
				library("devtools")
			}
			install_github('rCharts', 'ramnathv', ref='dev')
		}
		library("rCharts")
	}
	
	## Create outdir
	dir.create(file.path(prefix, outdir), showWarnings = FALSE, recursive = TRUE)
	
	## Locate Rmd
	template <- system.file(file.path("basicExploration", "basicExploration.Rmd"), package="derfinder2")
	
	## Load knitcitations with a clean bibliography
	cleanbib()
	cite_options(tooltip=TRUE)

	## I made my own citing function since citep() doesn't work like I want to with
	## urls that are not really pages themselves like part of a GitHub repo.
	mycitep <- function(x, short=NULL, year=substr(date(), 21, 24), tooltip=TRUE) {
		res <- tmp <- citep(x)
		if(!is.null(short)) {
			res <- gsub("></a>", paste0(">", short, "</a>"), tmp)
		}		
		if(tooltip) {
			res <- gsub("\\?\\?\\?\\?", year, res)
		}
		res <- gsub("span> ", "span>", tmp)
		res
	}
	
	## Write bibliography information
	write.bibtex(c("knitcitations" = citation("knitcitations"), "derfinder2" = citation("derfinder2"), "knitrBootstrap" = citation("knitrBootstrap"), "ggbio" = citation("ggbio"), "ggplot2" = citation("ggplot2"), "rCharts" = citation("rCharts"), "knitr" = citation("knitr")[1]), file = file.path(prefix, outdir, "references.bib"))
	bib <- read.bibtex(file.path(prefix, outdir, "references.bib"))
	
	## Load files
	if(is.null(fullRegions)) load(file.path(prefix, "fullRegions.Rdata"))
	if(is.null(fullNullSummary)) load(file.path(prefix, "fullNullSummary.Rdata"))
	if(is.null(fullAnnotatedRegions)) load(file.path(prefix, "fullAnnotatedRegions.Rdata"))
	if(is.null(optionsStats)) load(file.path(prefix, dir(prefix, pattern="chr")[1], "optionsStats.Rdata"))

	## Require fullCov
	if(makeBestClusters) {
		stopifnot(!is.null(fullCov))
	} else {
		nBestClusters <- 0
	}

	##### Setup chunk options
	## Are there any null regions? If not, then there won't be any p-values either.
	nullExist <- length(fullNullSummary) > 0
	## Were permutations used?
	seeds <- optionsStats$seeds
	usedPermutations <- length(optionsStats$nPermute) > 0 & !is.null(seeds)
	## Are there significant (by q-value) regions?
	idx.sig <- which(as.logical(fullRegions$significantQval))
	hasSig <- length(idx.sig) > 0
	## Save the call
	theCall <- match.call()
	
	## Generate report
	tmpdir <- getwd()
	setwd(file.path(prefix, outdir))
	res <- knit_bootstrap(input=template, output=output, code_style='Brown Paper', chooser=c('boot', 'code'), show_code=FALSE)
	
	## Open
	if (browse) browseURL(output)
	setwd(tmpdir)
		
	## Finish
	return(invisible(res))
}