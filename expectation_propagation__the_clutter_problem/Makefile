.PHONY : all clean

img_files := dividing_gaussians_w_same_var.png interactive_plot.png posteriors.png
img_files_50 := $(foreach file, $(img_files), \
									$(basename $(file))_50$(suffix $(file)))


all : clutter.pdf clean-tex

clutter.pdf : clutter.tex makra.tex clutter.bib $(img_files_50)
	pdflatex $<
	bibtex $(basename $<)
	pdflatex $<
	pdflatex $<

define halving_tpt
$$(basename $(1))_50$$(suffix $(1)) : $(1)
	convert -quality 90 -resize 50% "$$<" "$$@"
endef

$(foreach img, $(img_files), $(eval $(call halving_tpt, $(img))))

clean-tex : | clutter.pdf
	rm -f `ls clutter.* | grep -v '\.\(tex\|bib\|pdf\|py\)$$'`

clean : clean-tex
	rm -f clutter.pdf
	rm -f $(img_files_50)
