# Makefile to use knitrBootstrap for package vignettes
# Modified from https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks/ggbio/vignettes/Makefile
# and
# https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks/BiocCheck/vignettes/Makefile

all:
	cp realVignettes/*.Rmd .
	"$(R_HOME)/bin$(R_ARCH_BIN)/Rscript" -e "library(knitrBootstrap); knitrBoostrapFlag <- packageVersion('knitrBootstrap') < '1.0.0'; if(knitrBoostrapFlag) knit_bootstrap('derfinder.Rmd', chooser=c('boot', 'code'), show_code=TRUE) else {library(rmarkdown); render('derfinder.Rmd', bootstrap_document(theme.chooser=TRUE, highlight.chooser=TRUE))}"
	"$(R_HOME)/bin$(R_ARCH_BIN)/Rscript" -e "library(knitrBootstrap); knitrBoostrapFlag <- packageVersion('knitrBootstrap') < '1.0.0'; if(knitrBoostrapFlag) knit_bootstrap('derfinderAdvanced.Rmd', chooser=c('boot', 'code'), show_code=TRUE) else {library(rmarkdown); render('derfinderAdvanced.Rmd', bootstrap_document(theme.chooser=TRUE, highlight.chooser=TRUE))}"

clean:
	rm -f *.bib *.md *.bw
	rm -fr analysisResults
	cp fakeVignettes/*.Rmd .
