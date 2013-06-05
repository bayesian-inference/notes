# author: Ondrej Platek 2012, No warranty
# Thanks to Matej Korvas for patches
NAME = main
TEX_CALL = pdflatex -shell-escape
# TEX_CALL = xelatex 

# for fully functional pdf 
# all: bib 
# for regular compilation -> quicker
all: pdf 

pdf: ${NAME}.pdf
	# okular --unique "${NAME}.pdf" >/dev/null 2>&1 &
	evince "${NAME}.pdf" >/dev/null 2>&1 &

${NAME}.pdf: ${NAME}.tex $(wildcard *.tex)
	$(TEX_CALL) $<

.PHONY: clean bib

bib: ${NAME}.tex $(wildcard *.tex) refer.bib
	$(TEX_CALL) $<
	bibtex ${NAME}
	$(TEX_CALL) $<
	$(TEX_CALL) $<

clean:
	rm -f *.log  *.dvi *.aux *.toc *.lof *.out ${NAME}.ps ${NAME}-pics.pdf ${NAME}.pdf *.bbl *.blg
