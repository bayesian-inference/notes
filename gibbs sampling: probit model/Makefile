# author: Ondrej Platek 2012, No warranty
NAME = main

# change all: bib to all: pdf if you want to use it regularly
all: bib    

pdf: ${NAME}.pdf
	# okular --unique "${NAME}.pdf" >/dev/null 2>&1 &
	evince "${NAME}.pdf" >/dev/null 2>&1 &

${NAME}.pdf: ${NAME}.tex $(wildcard *.tex)
	pdflatex -shell-escape $<

.PHONY: clean bib

bib: ${NAME}.tex $(wildcard *.tex) refer.bib
	pdflatex -shell-escape $<
	bibtex ${NAME}
	pdflatex -shell-escape $<
	pdflatex -shell-escape $<

clean:
	rm -f *.log  *.dvi *.aux *.toc *.lof *.out ${NAME}.ps ${NAME}-pics.pdf ${NAME}.pdf *.bbl *.blg
