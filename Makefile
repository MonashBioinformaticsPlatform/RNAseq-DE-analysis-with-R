
RMDS=RNAseq_DE_analysis_with_R.Rmd \
	RNAseq_DE_analysis_with_R.Solutions.Rmd \
	index.Rmd

HTMLS=$(patsubst %.Rmd,%.html,$(RMDS))

# Create stripped down versions of .Rmd files
RS=data/RNAseq_DE_analysis_with_R.R \
	data/RNAseq_DE_analysis_with_R.Solutions.R

# Create unevaluated versions (compact teacher's notes)
#UNEVALS=RNAseq_DE_analysis_with_R.html
#topics/programming_uneval.html \
        topics/tidyverse_uneval.html \
        topics/sequences_and_features_uneval.html

#all : $(RS) $(HTMLS) $(UNEVALS) rnaseq-with-r-files.zip
all : $(RS) $(HTMLS) data.zip

%.html : %.Rmd
	Rscript -e 'rmarkdown::render("$<", "all")'
	
#%_uneval.html : %.Rmd Makefile
#	python unevalify.py <$< >topics/temp.Rmd
#	Rscript -e 'rmarkdown::render("topics/temp.Rmd", "all")'
#	mv topics/temp.html $@
#	rm topics/temp.Rmd

data/%.R : %.Rmd purify.py
	python purify.py <$< >$@

data.zip : data/* $(RS)
	zip -FSr data.zip data

clean :
	rm $(HTMLS) $(RS) $(UNEVALS) data.zip
