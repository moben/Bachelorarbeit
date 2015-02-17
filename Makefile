
bib = references.bib
data = $(wildcard data/*.txt)
cleanfiles = 				\
	Bachelorarbeit.aux		\
	Bachelorarbeit.bbl		\
	Bachelorarbeit.bcf		\
	Bachelorarbeit.blg		\
	Bachelorarbeit.log		\
	Bachelorarbeit.out		\
	Bachelorarbeit.run.xml		\
	Bachelorarbeit.synctex.gz	\
	Bachelorarbeit.toc

all: Bachelorarbeit.pdf

%.tex: %.lyx
	lyx -e pdflatex $< 

%.bcf: %.tex $(bib)
	pdflatex $<

%.bbl: %.bcf
	biber $<

%.pdf: %.tex %.bbl $(data)
	pdflatex $<
	pdflatex $<

clean:
	-rm $(cleanfiles)

