# Makefile to use knitrBoostrap for package vignettes
# Modified https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks/ggbio/vignettes/Makefile

# put all HTML targets here, separated by spaces
HTMLS= derfinder.html

all: $(HTMLS) 

clean:
	cp fakevignettes/*.Rmd .
	rm -f *.bib *.md
	
%.html: %.Rmd
	cp realvignettes/*.Rmd .
	$(R_HOME)/bin/Rscript -e "library(knitrBootstrap); knitrBoostrapFlag <- packageVersion('knitrBootstrap') < '1.0.0'; if(knitrBoostrapFlag) knit_bootstrap('$*.Rmd', chooser=c('boot', 'code'), show_code=TRUE) else library(rmarkdown); render('$*.Rmd', bootstrap_document(theme.chooser=TRUE, highlight.chooser=TRUE))"
