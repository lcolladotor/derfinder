#' Generate a HTML report exploring the basic results
#'
#' This function generates a HTML report exploring the basic results. 
#' 
#' @param prefix The main data directory path where \link{mergeResults} was run. It should be the same as \code{mergeResults(prefix)}.
#' @param outdir The name of output directory relative to \code{prefix}.
#' @param output The name output HTML file.
#' @param project The title of the project.
#' @param browse If \code{TRUE} the HTML report is opened in your browser once it's completed.
#' @param makeBestPlots If \code{TRUE}, \link{plotRegion} is used on the \code{nBest} regions by area. These plots take some time to make
#' @param nBest The number of best by area plots to make.
#' @param fullCov A list where each element is the result from \link{loadCoverage} used with \code{cutoff=NULL}. The elements of the list should be named according to the chromosome number.
#' @param hg19 If \code{TRUE} then the reference is assumed to be hg19 and chromosome lengths as well as the default transcription database (TxDb.Hsapiens.UCSC.hg19.knownGene) will be used.
#' @param p.ideos A list where each element is the result of \link[ggbio]{plotIdeogram}. If it's \code{NULL} and \code{hg19=TRUE} then they are created for the hg19 human reference.
#' @param txdb Specify the transcription database to use for making the plots for the top regions by area. If \code{NULL} and \code{hg19=TRUE} then TxDb.Hsapiens.UCSC.hg19.knownGene is used.
#' @param installMissing If \code{TRUE} all missing required packages are installed. Note that some are development versions hosted in GitHub.
#' @param device The graphical device used when knitting. See more at http://yihui.name/knitr/options (dev argument).
#' @param fullTime Part of the output of \link{mergeResults}. Specify it only if you have already loaded it in memory.
#' @param fullRegions Part of the output of \link{mergeResults}. Specify it only if you have already loaded it in memory.
#' @param fullNullSummary Part of the output of \link{mergeResults}. Specify it only if you have already loaded it in memory.
#' @param optionsStats Part of the output of \link{analyzeChr}. Specify it only if you have already loaded it in memory.
#'
#' @return An HTML report with a basic exploration of the results.
#'
#' @author Leonardo Collado-Torres
#' @seealso \link{mergeResults}, \link{analyzeChr}
#' @export
#' @importMethodsFrom IRanges as.numeric
#'
#' @examples
#' \dontrun{
#' generateReport(prefix="run1", makeBestPlots=FALSE)
#' }


generateReport <- function(prefix, outdir="basicExploration", output="basicExploration.html", project=prefix, browse=interactive(), makeBestPlots=TRUE, nBest=2, fullCov=NULL, hg19=TRUE, p.ideos=NULL, txdb=NULL, installMissing=TRUE, device="CairoPNG", fullTime=NULL, fullRegions=NULL, fullNullSummary=NULL, optionsStats=NULL) {
	## Load all required packages and install them if needed
	pkgs <- c("IRanges", "GenomicRanges", "TxDb.Hsapiens.UCSC.hg19.knownGene", "ggplot2", "gridExtra", "data.table", "knitcitations", "knitrBootstrap", "rCharts")
	pkgs.type <- c(rep("BioC", 3), rep("CRAN", 4), rep("git", 2))
	stopifnot(length(pkgs) == length(pkgs.type))
	
	dir.create(file.path(prefix, outdir), showWarnings = FALSE, recursive = TRUE)
	
	for(i in seq_len(length(pkgs))) {
		if(pkgs[i] != "TxDb.Hsapiens.UCSC.hg19.knownGene" | hg19) {
			## Ignored only if hg19 == FALSE
			tmp <- suppressMessages(require(pkgs[i], character.only=TRUE))
		} else {
			tmp <- TRUE
		}
		
		if(!tmp & installMissing) {
			if(pkgs.type[i] == "BioC") {
				source("http://bioconductor.org/biocLite.R")
				biocLite(pkgs[i])
			} else if (pkgs.type[i] == "CRAN") {
				install.packages(pkgs[i])
				
			} else {
				if(!require("devtools")) {
					install.packages("devtools")
					library("devtools")
				}
				if(pkgs[i] == "rCharts") {
					install_github('rCharts', 'ramnathv', ref='dev')
				} else if (pkgs[i] == "knitrBootstrap") {
					install_github(username='jimhester', repo='knitrBootstrap')
				}
			}
			library(pkgs[i])
		}
	}
	
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
	write.bibtex(c("knitcitations" = citation("knitcitations"), "derfinder2" = citation("derfinder2"), "knitrBootstrap" = citation("knitrBootstrap"), "ggbio" = citation("ggbio"), "ggplot2" = citation("ggplot2"), "rCharts" = citation("rCharts")), file = file.path(prefix, outdir, "references.bib"))
	bib <- read.bibtex(file.path(prefix, outdir, "references.bib"))
	
	## Load files
	if(is.null(fullTime)) load(file.path(prefix, "fullTime.Rdata"))
	if(is.null(fullRegions)) load(file.path(prefix, "fullRegions.Rdata"))
	if(is.null(fullNullSummary)) load(file.path(prefix, "fullNullSummary.Rdata"))
	if(is.null(optionsStats)) load(file.path(prefix, dir(prefix, pattern="chr")[1], "optionsStats.Rdata"))

	## Require fullCov
	if(makeBestPlots) {
		stopifnot(!is.null(fullCov))
	} else {
		nBest <- 0
	}

	##### Setup chunk options
	## Are there any null regions? If not, then there won't be any p-values either.
	nullExist <- length(fullNullSummary) > 0
	## Were permutations used?
	seeds <- optionsStats$seeds
	usedPermutations <- length(optionsStats$nPermute) > 0 & !is.null(seeds)
	
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