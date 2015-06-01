
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

# lyx doesn't update the .tex if it would be the same
# thus `touch`
%.tex: %.lyx
	lyx -e pdflatex $< 
	touch $@

%.bcf: %.tex $(bib)
	pdflatex -interaction=nonstopmode $<

%.bbl: %.bcf
	biber $<

%.pdf: %.tex %.bbl $(data)
	pdflatex -interaction=nonstopmode $<
	pdflatex -interaction=nonstopmode $<

clean:
	-rm $(cleanfiles)

